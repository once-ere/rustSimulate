//! The `QM` command family: one-dimensional quantum mechanics in the
//! notebook.
//!
//! This is the language front end for the `quantum` crate. The session
//! carries one quantum problem at a time — a grid, a potential sampled
//! on it, and a current wavefunction — built up command by command:
//!
//! ```text
//! def v(x) { 0.5 * x * x }     # the potential, as an ordinary function
//! qm grid -8 8 300             # domain and resolution
//! qm potential v               # sample v(x) onto the grid
//! qm states 5                  # the five lowest bound states
//! qm state 0                   # load the ground state as psi
//! qm packet -3 0.7 4           # ...or launch a Gaussian wavepacket
//! qm run 2 steps 400           # propagate (Crank-Nicolson, unitary)
//! qm energy                    # observables
//! ```
//!
//! # Why the potential is an ordinary user function
//!
//! `DEF` already gives the language first-class functions of one
//! argument, which is exactly what a potential is. `QM POTENTIAL v`
//! evaluates `v(x)` once at every grid point and stores the samples, so
//! there is no need for the expression evaluator to be re-entered during
//! propagation, and no need for a new kind of deferred-expression value.
//! Sampling happens **at command time**, so if you redefine `v` you must
//! re-issue `QM POTENTIAL v` — which is stated in the status output
//! rather than left to be discovered.

use quantum::qm1d::{Grid, Hamiltonian, Propagator, Wavefunction};

use crate::vm::{SimState, Value};

/// A `QM` subcommand. Numeric arguments are compiled as ordinary
/// expressions and popped from the stack, so `qm grid -8 8 2*150` works.
#[derive(Clone, Debug, PartialEq)]
pub enum QmCmd {
    /// Bare `QM`: report what is set up.
    Status,
    /// Pops n, x_max, x_min.
    Grid,
    /// Set the potential. See [`PotentialSpec`].
    Potential(PotentialSpec),
    /// Pops the particle mass.
    Mass,
    /// Pops hbar.
    Hbar,
    /// Pops k: report the k lowest bound-state energies.
    States,
    /// Pops n: load bound state n as the current wavefunction.
    LoadState,
    /// Pops k0, sigma, x0.
    Packet,
    /// Pops dt.
    Step,
    /// Pops the number of steps, then the elapsed time.
    Run,
    Norm,
    Energy,
    Position,
    Momentum,
    /// Pops b, a: probability in [a, b].
    Prob,
    /// The probability density, as a list.
    Density,
    /// Forget the whole quantum problem.
    Reset,
}

/// How a potential is specified.
///
/// A user function is the general case, but this language has **no
/// comparison operators**, so a piecewise potential — a square barrier,
/// a finite well — cannot be written as a `DEF` at all. Those are
/// exactly the canonical 1-D problems, so they are provided as named
/// shapes rather than left unreachable. This is a deliberate
/// workaround for a language limitation, not a preference for built-ins.
#[derive(Clone, Debug, PartialEq)]
pub enum PotentialSpec {
    /// Free particle.
    Zero,
    /// `v0` on `[x1, x2]`, zero elsewhere. Pops x2, x1, v0.
    Barrier,
    /// `-depth` on `[x1, x2]`, zero elsewhere. Pops x2, x1, depth.
    Well,
    /// A `DEF`ined function of one argument, sampled onto the grid.
    Named(String),
}

/// Everything the `QM` family needs to remember between commands.
#[derive(Clone, Debug)]
pub struct QmState {
    pub grid: Option<Grid>,
    pub potential: Option<Vec<f64>>,
    /// The function name the potential came from, for the status line.
    pub potential_name: Option<String>,
    pub mass: f64,
    pub hbar: f64,
    pub psi: Option<Wavefunction>,
    /// Elapsed propagation time, so observables can be reported against
    /// a clock rather than a step count.
    pub time: f64,
    /// Cached bound states, so `QM STATE n` after `QM STATES k` does not
    /// pay for a second diagonalisation.
    pub states: Option<(Vec<f64>, Vec<Vec<f64>>)>,
}

impl Default for QmState {
    fn default() -> Self {
        Self {
            grid: None,
            potential: None,
            potential_name: None,
            mass: 1.0,
            hbar: 1.0,
            psi: None,
            time: 0.0,
            states: None,
        }
    }
}

impl QmState {
    /// The Hamiltonian, if a grid and potential have both been set.
    fn hamiltonian(&self) -> Result<Hamiltonian, String> {
        let grid = self
            .grid
            .clone()
            .ok_or("QM: no grid — use `QM GRID <x_min> <x_max> <n>` first")?;
        let v = self
            .potential
            .clone()
            .ok_or("QM: no potential — use `QM POTENTIAL <function>` (or `QM POTENTIAL zero`)")?;
        Hamiltonian::new(grid, v, self.mass, self.hbar)
    }

    fn wavefunction(&self) -> Result<&Wavefunction, String> {
        self.psi.as_ref().ok_or_else(|| {
            "QM: no wavefunction — use `QM PACKET <x0> <sigma> <k0>` or `QM STATE <n>`".to_string()
        })
    }

    /// Anything that changes the Hamiltonian invalidates cached results.
    fn invalidate(&mut self) {
        self.states = None;
    }
}

fn pop_num(stack: &mut Vec<Value>) -> Result<f64, String> {
    match stack.pop() {
        Some(Value::Num(n)) => Ok(n),
        Some(other) => Err(format!("QM: expected a number, got {other}")),
        None => Err("QM: missing an argument".to_string()),
    }
}

/// A count argument: must be a whole, non-negative number.
fn pop_count(stack: &mut Vec<Value>, what: &str) -> Result<usize, String> {
    let v = pop_num(stack)?;
    if !v.is_finite() || v.fract() != 0.0 || v < 0.0 {
        return Err(format!("QM: {what} must be a whole number >= 0, got {v}"));
    }
    Ok(v as usize)
}

/// Execute one `QM` subcommand, returning the text to display.
pub fn exec_qm(
    cmd: &QmCmd,
    state: &mut SimState,
    stack: &mut Vec<Value>,
) -> Result<String, String> {
    match cmd {
        QmCmd::Status => {
            let q = &state.qm;
            let mut s = String::from("quantum (1-D):\n");
            match &q.grid {
                Some(g) => s.push_str(&format!(
                    "  grid      [{}, {}], {} interior points, h = {:.6}\n",
                    g.x_min,
                    g.x_max,
                    g.n,
                    g.h()
                )),
                None => s.push_str("  grid      (unset — QM GRID <x_min> <x_max> <n>)\n"),
            }
            match &q.potential_name {
                Some(n) => s.push_str(&format!(
                    "  potential {n}  (sampled when the command ran; re-issue \
                     QM POTENTIAL after editing it)\n"
                )),
                None => s.push_str("  potential (unset — QM POTENTIAL <function>)\n"),
            }
            s.push_str(&format!("  mass      {}\n  hbar      {}\n", q.mass, q.hbar));
            match &q.psi {
                Some(w) => s.push_str(&format!(
                    "  psi       set, norm = {:.12}, t = {}\n",
                    w.norm(),
                    q.time
                )),
                None => s.push_str("  psi       (unset — QM PACKET or QM STATE)\n"),
            }
            Ok(s.trim_end().to_string())
        }

        QmCmd::Grid => {
            let n = pop_count(stack, "the point count")?;
            let x_max = pop_num(stack)?;
            let x_min = pop_num(stack)?;
            let g = Grid::new(x_min, x_max, n)?;
            let h = g.h();
            state.qm.grid = Some(g);
            // The grid defines the sampling, so everything downstream is
            // stale: drop it rather than leave a potential of the wrong
            // length to fail confusingly later.
            state.qm.potential = None;
            state.qm.potential_name = None;
            state.qm.psi = None;
            state.qm.time = 0.0;
            state.qm.invalidate();
            Ok(format!(
                "grid [{x_min}, {x_max}] with {n} interior points, h = {h:.6} \
                 (potential and psi cleared)"
            ))
        }

        QmCmd::Potential(spec) => {
            let grid = state
                .qm
                .grid
                .clone()
                .ok_or("QM POTENTIAL: set a grid first (QM GRID <x_min> <x_max> <n>)")?;
            let (v, label): (Vec<f64>, String) = match spec {
                PotentialSpec::Zero => (vec![0.0; grid.n], "zero".to_string()),
                PotentialSpec::Barrier | PotentialSpec::Well => {
                    let x2 = pop_num(stack)?;
                    let x1 = pop_num(stack)?;
                    let amp = pop_num(stack)?;
                    if !amp.is_finite() || !x1.is_finite() || !x2.is_finite() {
                        return Err("QM POTENTIAL: arguments must be finite".to_string());
                    }
                    let (lo, hi) = if x1 <= x2 { (x1, x2) } else { (x2, x1) };
                    let sign = if matches!(spec, PotentialSpec::Well) { -1.0 } else { 1.0 };
                    let vv = (0..grid.n)
                        .map(|i| {
                            let x = grid.x(i);
                            if x >= lo && x <= hi {
                                sign * amp
                            } else {
                                0.0
                            }
                        })
                        .collect();
                    let kind = if sign < 0.0 { "well" } else { "barrier" };
                    // A feature narrower than the grid spacing is not
                    // represented at all, and would silently behave as
                    // if it were absent.
                    if hi - lo < grid.h() {
                        return Err(format!(
                            "QM POTENTIAL {kind}: the region [{lo}, {hi}] is narrower than the \
                             grid spacing h = {:.6}, so it would fall between points — use a \
                             finer grid",
                            grid.h()
                        ));
                    }
                    (vv, format!("{kind} {amp} on [{lo}, {hi}]"))
                }
                PotentialSpec::Named(name) => {
                    if !state.functions.contains_key(name) {
                        return Err(format!(
                            "QM POTENTIAL: no function `{name}` — define one with \
                             `DEF {name}(x) {{ ... }}`, or use one of the built-in shapes \
                             (zero, barrier, well)"
                        ));
                    }
                    let mut out = Vec::with_capacity(grid.n);
                    for i in 0..grid.n {
                        let val = crate::vm::call_user_function_public(
                            name,
                            vec![Value::Num(grid.x(i))],
                            state,
                        )?;
                        match val {
                            Value::Num(y) => out.push(y),
                            other => {
                                return Err(format!(
                                    "QM POTENTIAL: `{name}(x)` must return a number, got {other}"
                                ))
                            }
                        }
                    }
                    (out, name.clone())
                }
            };
            let lo = v.iter().cloned().fold(f64::INFINITY, f64::min);
            let hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let n = grid.n;
            state.qm.potential = Some(v);
            state.qm.potential_name = Some(label.clone());
            state.qm.psi = None;
            state.qm.time = 0.0;
            state.qm.invalidate();
            Ok(format!(
                "potential `{label}` sampled at {n} points, V in [{lo}, {hi}] (psi cleared)"
            ))
        }

        QmCmd::Mass => {
            let m = pop_num(stack)?;
            if !m.is_finite() || m <= 0.0 {
                return Err(format!("QM MASS: must be finite and positive, got {m}"));
            }
            state.qm.mass = m;
            state.qm.invalidate();
            Ok(format!("mass = {m}"))
        }

        QmCmd::Hbar => {
            let h = pop_num(stack)?;
            if !h.is_finite() || h <= 0.0 {
                return Err(format!("QM HBAR: must be finite and positive, got {h}"));
            }
            state.qm.hbar = h;
            state.qm.invalidate();
            Ok(format!("hbar = {h}"))
        }

        QmCmd::States => {
            let k = pop_count(stack, "the state count")?;
            if k == 0 {
                return Err("QM STATES: ask for at least one state".to_string());
            }
            let ham = state.qm.hamiltonian()?;
            let (e, v) = ham.bound_states(k)?;
            let mut s = format!("{k} lowest bound state(s):\n");
            for (i, en) in e.iter().enumerate() {
                s.push_str(&format!("  E[{i}] = {en:.12}\n"));
            }
            state.qm.states = Some((e, v));
            Ok(s.trim_end().to_string())
        }

        QmCmd::LoadState => {
            let n = pop_count(stack, "the state index")?;
            let ham = state.qm.hamiltonian()?;
            // reuse the cached diagonalisation when it is deep enough
            let need = n + 1;
            let have = state.qm.states.as_ref().map(|(e, _)| e.len()).unwrap_or(0);
            if have < need {
                let (e, v) = ham.bound_states(need)?;
                state.qm.states = Some((e, v));
            }
            let (e, v) = state.qm.states.as_ref().expect("just filled");
            let w = Wavefunction::from_real(ham.grid.clone(), &v[n])?;
            let en = e[n];
            state.qm.psi = Some(w);
            state.qm.time = 0.0;
            Ok(format!("psi = bound state {n}, E = {en:.12}, t reset to 0"))
        }

        QmCmd::Packet => {
            let k0 = pop_num(stack)?;
            let sigma = pop_num(stack)?;
            let x0 = pop_num(stack)?;
            let grid = state
                .qm
                .grid
                .clone()
                .ok_or("QM PACKET: set a grid first (QM GRID <x_min> <x_max> <n>)")?;
            let w = Wavefunction::gaussian(grid, x0, sigma, k0)?;
            let edge = w.edge_probability(0.05);
            state.qm.psi = Some(w);
            state.qm.time = 0.0;
            let warn = if edge > 1e-6 {
                format!(
                    "\n  warning: {edge:.3e} of the packet already sits within 5% of a wall; \
                     the walls REFLECT, so widen the domain"
                )
            } else {
                String::new()
            };
            Ok(format!(
                "psi = Gaussian packet at x0 = {x0}, sigma = {sigma}, k0 = {k0}, t = 0{warn}"
            ))
        }

        QmCmd::Step | QmCmd::Run => {
            let (dt, steps) = if matches!(cmd, QmCmd::Step) {
                (pop_num(stack)?, 1usize)
            } else {
                let n = pop_count(stack, "the step count")?;
                let t = pop_num(stack)?;
                if n == 0 {
                    return Err("QM RUN: the step count must be at least 1".to_string());
                }
                (t / n as f64, n)
            };
            if !dt.is_finite() || dt == 0.0 {
                return Err(format!("QM: the time step must be finite and non-zero, got {dt}"));
            }
            let ham = state.qm.hamiltonian()?;
            let mut w = state.qm.wavefunction()?.clone();
            let n0 = w.norm();
            let prop = Propagator::new(ham.clone(), dt)?;
            prop.run(&mut w, steps)?;
            let n1 = w.norm();
            let drift = (n1 / n0 - 1.0).abs();
            let edge = w.edge_probability(0.05);
            state.qm.time += dt * steps as f64;
            let t = state.qm.time;
            let e = w.energy(&ham);
            state.qm.psi = Some(w);
            let warn = if edge > 1e-4 {
                format!(
                    "\n  warning: {edge:.3e} of the probability is within 5% of a wall — \
                     the walls REFLECT, so results past this point are suspect"
                )
            } else {
                String::new()
            };
            Ok(format!(
                "t = {t} ({steps} step(s) of dt = {dt}), <E> = {e:.12}, \
                 norm drift = {drift:.3e}{warn}"
            ))
        }

        QmCmd::Norm => Ok(format!("{:.15}", state.qm.wavefunction()?.norm())),
        QmCmd::Position => Ok(format!("{:.15}", state.qm.wavefunction()?.position())),
        QmCmd::Momentum => {
            let hbar = state.qm.hbar;
            Ok(format!("{:.15}", state.qm.wavefunction()?.momentum(hbar)))
        }
        QmCmd::Energy => {
            let ham = state.qm.hamiltonian()?;
            Ok(format!("{:.15}", state.qm.wavefunction()?.energy(&ham)))
        }
        QmCmd::Prob => {
            let b = pop_num(stack)?;
            let a = pop_num(stack)?;
            Ok(format!("{:.15}", state.qm.wavefunction()?.probability_in(a, b)))
        }
        QmCmd::Density => {
            let d = state.qm.wavefunction()?.density();
            stack.push(Value::List(d.into_iter().map(Value::Num).collect()));
            Ok(String::new())
        }

        QmCmd::Reset => {
            state.qm = QmState::default();
            Ok("quantum state cleared".to_string())
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::vm::{execute_line, SimState};

    fn run(lines: &[&str]) -> (SimState, Vec<String>) {
        let mut st = SimState::default();
        let mut out = Vec::new();
        for l in lines {
            let v = execute_line(l, &mut st).unwrap_or_else(|e| panic!("`{l}` failed: {e}"));
            out.push(v.to_string());
        }
        (st, out)
    }

    /// The harmonic oscillator, entirely through the language: define a
    /// potential as an ordinary function, sample it, diagonalise.
    #[test]
    fn harmonic_oscillator_through_the_language() {
        let (_, out) = run(&[
            "def v(x) { 0.5 * x * x }",
            "qm grid -8 8 250",
            "qm potential v",
            "qm states 4",
        ]);
        let text = out.last().unwrap();
        // E_n = n + 1/2
        for (n, want) in [(0, 0.5), (1, 1.5), (2, 2.5), (3, 3.5)] {
            let line = text
                .lines()
                .find(|l| l.trim_start().starts_with(&format!("E[{n}]")))
                .unwrap_or_else(|| panic!("no E[{n}] in:\n{text}"));
            let got: f64 = line.split('=').nth(1).unwrap().trim().parse().unwrap();
            assert!((got - want).abs() < 5e-3, "E[{n}] = {got}, want {want}");
        }
    }

    /// A bound state loaded into psi must have the eigenvalue as its
    /// energy, and stay put under propagation.
    #[test]
    fn a_loaded_eigenstate_is_stationary() {
        let (_, out) = run(&[
            "def v(x) { 0.5 * x * x }",
            "qm grid -8 8 120",
            "qm potential v",
            "qm state 1",
            "qm energy",
            "qm run 2 steps 200",
            "qm energy",
        ]);
        let e_before: f64 = out[4].trim().parse().unwrap();
        let e_after: f64 = out[6].trim().parse().unwrap();
        assert!((e_before - 1.5).abs() < 5e-3, "E = {e_before}, want ~1.5");
        assert!(
            (e_after - e_before).abs() < 1e-9,
            "energy moved: {e_before} -> {e_after}"
        );
    }

    /// Propagation is unitary, so the reported norm drift must be tiny.
    #[test]
    fn propagation_conserves_the_norm() {
        let (_, out) = run(&[
            "qm grid -40 40 400",
            "qm potential zero",
            "qm packet -10 2 2",
            "qm run 3 steps 300",
            "qm norm",
        ]);
        let norm: f64 = out[4].trim().parse().unwrap();
        assert!((norm - 1.0).abs() < 1e-10, "norm = {norm}");
        assert!(out[3].contains("norm drift"), "run should report drift");
    }

    /// The walls reflect. A packet aimed at one must SAY so rather than
    /// quietly returning a wrong scattering answer.
    #[test]
    fn approaching_a_wall_is_reported() {
        let (_, out) = run(&[
            "qm grid -10 10 200",
            "qm potential zero",
            "qm packet 0 1 8",
            "qm run 2 steps 200",
        ]);
        assert!(
            out[3].contains("warning") && out[3].contains("wall"),
            "expected a wall warning, got: {}",
            out[3]
        );
    }

    /// Changing the grid invalidates the potential: a stale potential of
    /// the wrong length would otherwise fail much later and confusingly.
    #[test]
    fn changing_the_grid_clears_downstream_state() {
        let mut st = SimState::default();
        for l in [
            "def v(x) { 0.5 * x * x }",
            "qm grid -8 8 100",
            "qm potential v",
            "qm state 0",
        ] {
            execute_line(l, &mut st).unwrap();
        }
        assert!(st.qm.psi.is_some());
        execute_line("qm grid -8 8 150", &mut st).unwrap();
        assert!(st.qm.potential.is_none(), "potential should be cleared");
        assert!(st.qm.psi.is_none(), "psi should be cleared");
        let e = execute_line("qm states 2", &mut st).unwrap_err();
        assert!(e.contains("no potential"), "got: {e}");
    }

    /// The built-in shapes exist because the language has no comparison
    /// operators, so a piecewise potential cannot be written as a DEF.
    /// Tunnelling through a square barrier is the canonical 1-D problem,
    /// so it must be reachable.
    #[test]
    fn a_square_barrier_transmits_and_reflects() {
        let (_, out) = run(&[
            "qm grid -60 60 1200",
            "qm potential barrier 2.5 0 1",
            "qm packet -20 2 2",
            "qm run 20 steps 2000",
            "qm prob 1 60",
            "qm prob -60 0",
        ]);
        let t: f64 = out[4].trim().parse().unwrap();
        let r: f64 = out[5].trim().parse().unwrap();
        // E0 = 2 against a 2.5 barrier: substantial tunnelling, and the
        // two channels must exhaust the probability.
        assert!((0.25..0.40).contains(&t), "T = {t}");
        assert!((t + r - 1.0).abs() < 1e-4, "T + R = {}", t + r);
    }

    /// A well binds states below zero; a barrier of the same footprint
    /// binds none.
    #[test]
    fn a_well_binds_states_and_a_barrier_does_not() {
        let (_, out) = run(&[
            "qm grid -20 20 300",
            "qm potential well 5, -2, 2",
            "qm states 3",
        ]);
        let bound = out[2]
            .lines()
            .filter_map(|l| l.split('=').nth(1))
            .filter_map(|v| v.trim().parse::<f64>().ok())
            .filter(|&e| e < 0.0)
            .count();
        assert!(bound >= 2, "expected bound states below 0, got {bound} in:\n{}", out[2]);

        let (_, out2) = run(&[
            "qm grid -20 20 300",
            "qm potential barrier 5, -2, 2",
            "qm states 3",
        ]);
        let neg = out2[2]
            .lines()
            .filter_map(|l| l.split('=').nth(1))
            .filter_map(|v| v.trim().parse::<f64>().ok())
            .filter(|&e| e < 0.0)
            .count();
        assert_eq!(neg, 0, "a barrier must bind nothing below zero");
    }

    /// Negative arguments need commas, because `5 -2` is subtraction.
    /// Both spellings must behave, and the failure must be a parse
    /// error rather than a silently wrong potential.
    #[test]
    fn negative_arguments_need_comma_separation() {
        let mut st = SimState::default();
        execute_line("qm grid -20 20 100", &mut st).unwrap();
        // commas: unambiguous
        assert!(execute_line("qm potential well 5, -2, 2", &mut st).is_ok());
        // spaces with a negative: `5 -2` is subtraction, so an argument
        // goes missing and this MUST fail loudly
        assert!(
            execute_line("qm potential well 5 -2 2", &mut st).is_err(),
            "the ambiguous spelling must be a parse error, not a wrong potential"
        );
        // parentheses work too
        assert!(execute_line("qm potential well 5 (-2) 2", &mut st).is_ok());
        // all-positive space separation still reads fine
        assert!(execute_line("qm potential barrier 2.5 0 1", &mut st).is_ok());
    }

    /// A feature narrower than the grid spacing falls between points and
    /// would silently act as if absent. Refuse instead.
    #[test]
    fn a_subgrid_feature_is_refused() {
        let mut st = SimState::default();
        execute_line("qm grid -10 10 50", &mut st).unwrap();
        let e = execute_line("qm potential barrier 5 0 0.01", &mut st).unwrap_err();
        assert!(e.contains("narrower than the grid spacing"), "got: {e}");
    }

    /// Every ordering mistake gets a message naming the fix.
    #[test]
    fn missing_prerequisites_are_reported_helpfully() {
        let mut st = SimState::default();
        assert!(execute_line("qm states 3", &mut st).unwrap_err().contains("no grid"));
        execute_line("qm grid -5 5 50", &mut st).unwrap();
        assert!(execute_line("qm states 3", &mut st).unwrap_err().contains("no potential"));
        execute_line("qm potential zero", &mut st).unwrap();
        assert!(execute_line("qm norm", &mut st).unwrap_err().contains("no wavefunction"));
        assert!(execute_line("qm potential nosuch", &mut st).unwrap_err().contains("DEF"));
        assert!(execute_line("qm grid 5 -5 10", &mut st).is_err(), "reversed bounds");
        assert!(execute_line("qm grid -5 5 0", &mut st).is_err(), "n = 0");
        assert!(execute_line("qm mass -1", &mut st).is_err());
        assert!(execute_line("qm hbar 0", &mut st).is_err());
    }

    /// Arguments are full expressions, not just literals.
    #[test]
    fn arguments_are_expressions() {
        let (st, _) = run(&[
            "let w = 4",
            "qm grid 0 - w w 2 * 50",
            "qm potential zero",
        ]);
        let g = st.qm.grid.as_ref().unwrap();
        assert_eq!(g.n, 100);
        assert_eq!(g.x_max, 4.0);
    }
}
