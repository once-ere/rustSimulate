//! The `QM3` command family: three-dimensional quantum mechanics.
//!
//! A third family alongside `QM` and `QM2`, for the same reason `QM2` is
//! separate from `QM`: the argument lists differ throughout, and a
//! hidden dimensionality mode that silently reinterprets your commands
//! would be worse than a third word.
//!
//! ```text
//! def v(x, y, z) { 0.5 * (x * x + y * y + z * z) }
//! qm3 grid -6 6 24, -6 6 24, -6 6 24
//! qm3 potential v
//! qm3 states 4            # 1.5, then 2.5 three times
//! qm3 packet -2 0 0, 1 1 1, 2 0 0
//! qm3 run 1 steps 100
//! ```
//!
//! # What is affordable
//!
//! Propagation is `O(nx*ny*nz)` per step and comfortable well past 64³.
//! `QM3 STATES` is not: the Lanczos solver reorthogonalises fully and
//! stores its whole Krylov basis, so it is practical to roughly 40³.
//! Asking for more is refused up front rather than left to exhaust
//! memory.

use quantum::qm3d::{BoundStates3, Grid3, Hamiltonian3, Propagator3, Wavefunction3};

use crate::vm::{SimState, Value};

/// A `QM3` subcommand.
#[derive(Clone, Debug, PartialEq)]
pub enum Qm3Cmd {
    Status,
    /// Pops nz, z_max, z_min, ny, y_max, y_min, nx, x_max, x_min.
    Grid,
    /// `zero`, or a `DEF`ined function of three arguments.
    Potential(String),
    /// Pops kz, ky, kx, sz, sy, sx, z0, y0, x0.
    Packet,
    Step,
    Run,
    Norm,
    Energy,
    Centroid,
    /// Pops zb, za, yb, ya, xb, xa.
    Prob,
    States,
    LoadState,
    Absorb,
    AbsorbOff,
    Reset,
}

/// The 3-D problem carried by a session.
#[derive(Clone, Debug, Default)]
pub struct Qm3State {
    pub grid: Option<Grid3>,
    pub potential: Option<Vec<f64>>,
    pub potential_name: Option<String>,
    pub mass: f64,
    pub hbar: f64,
    pub psi: Option<Wavefunction3>,
    pub time: f64,
    pub absorber: Option<(f64, f64, f64)>,
    pub states: Option<BoundStates3>,
}

impl Qm3State {
    pub fn fresh() -> Self {
        Self { mass: 1.0, hbar: 1.0, ..Default::default() }
    }

    fn hamiltonian(&self) -> Result<Hamiltonian3, String> {
        let grid = self.grid.clone().ok_or(
            "QM3: no grid — use `QM3 GRID <x0> <x1> <nx>, <y0> <y1> <ny>, <z0> <z1> <nz>`",
        )?;
        let v = self
            .potential
            .clone()
            .ok_or("QM3: no potential — use `QM3 POTENTIAL <function of x,y,z>` or `zero`")?;
        let h = Hamiltonian3::new(grid, v, self.mass, self.hbar)?;
        match self.absorber {
            Some((w, s, p)) => h.with_absorber(w, s, p),
            None => Ok(h),
        }
    }

    fn wavefunction(&self) -> Result<&Wavefunction3, String> {
        self.psi
            .as_ref()
            .ok_or_else(|| "QM3: no wavefunction — use `QM3 PACKET ...`".to_string())
    }
}

fn pop_num(stack: &mut Vec<Value>) -> Result<f64, String> {
    match stack.pop() {
        Some(Value::Num(n)) => Ok(n),
        Some(other) => Err(format!("QM3: expected a number, got {other}")),
        None => Err("QM3: missing an argument".to_string()),
    }
}

fn pop_count(stack: &mut Vec<Value>, what: &str) -> Result<usize, String> {
    let v = pop_num(stack)?;
    if !v.is_finite() || v.fract() != 0.0 || v < 0.0 {
        return Err(format!("QM3: {what} must be a whole number >= 0, got {v}"));
    }
    Ok(v as usize)
}

/// Above this the Lanczos basis alone runs to gigabytes, so the command
/// says so instead of trying and failing slowly.
const EIGEN_POINT_LIMIT: usize = 70_000;

/// Execute one `QM3` subcommand.
pub fn exec_qm3(
    cmd: &Qm3Cmd,
    state: &mut SimState,
    stack: &mut Vec<Value>,
) -> Result<String, String> {
    match cmd {
        Qm3Cmd::Status => {
            let q = &state.qm3;
            let mut s = String::from("quantum (3-D, ADI):\n");
            match &q.grid {
                Some(g) => s.push_str(&format!(
                    "  grid      x [{}, {}] x {}, y [{}, {}] x {}, z [{}, {}] x {}\n            \
                     {} points, h = ({:.4}, {:.4}, {:.4})\n",
                    g.x_min, g.x_max, g.nx, g.y_min, g.y_max, g.ny, g.z_min, g.z_max, g.nz,
                    g.len(), g.hx(), g.hy(), g.hz()
                )),
                None => s.push_str("  grid      (unset)\n"),
            }
            s.push_str(&match &q.potential_name {
                Some(n) => format!("  potential {n} (sampled when the command ran)\n"),
                None => "  potential (unset)\n".to_string(),
            });
            s.push_str(&format!("  mass      {}\n  hbar      {}\n", q.mass, q.hbar));
            s.push_str(&match q.absorber {
                Some((w, st, p)) => format!("  absorber  width {w}, strength {st}, power {p}\n"),
                None => "  absorber  off — all six faces REFLECT\n".to_string(),
            });
            s.push_str(&match &q.psi {
                Some(w) => {
                    let (cx, cy, cz) = w.centroid();
                    format!(
                        "  psi       set, norm {:.12}, centroid ({cx:.4}, {cy:.4}, {cz:.4}), \
                         t = {}\n",
                        w.norm(),
                        q.time
                    )
                }
                None => "  psi       (unset — QM3 PACKET)\n".to_string(),
            });
            Ok(s.trim_end().to_string())
        }

        Qm3Cmd::Grid => {
            let nz = pop_count(stack, "nz")?;
            let z_max = pop_num(stack)?;
            let z_min = pop_num(stack)?;
            let ny = pop_count(stack, "ny")?;
            let y_max = pop_num(stack)?;
            let y_min = pop_num(stack)?;
            let nx = pop_count(stack, "nx")?;
            let x_max = pop_num(stack)?;
            let x_min = pop_num(stack)?;
            let g = Grid3::new(x_min, x_max, nx, y_min, y_max, ny, z_min, z_max, nz)?;
            let n = g.len();
            let (hx, hy, hz) = (g.hx(), g.hy(), g.hz());
            state.qm3.grid = Some(g);
            state.qm3.potential = None;
            state.qm3.potential_name = None;
            state.qm3.psi = None;
            state.qm3.time = 0.0;
            state.qm3.states = None;
            Ok(format!(
                "grid {nx} x {ny} x {nz} = {n} points, h = ({hx:.6}, {hy:.6}, {hz:.6}) \
                 (potential and psi cleared)"
            ))
        }

        Qm3Cmd::Potential(name) => {
            let grid = state.qm3.grid.clone().ok_or("QM3 POTENTIAL: set a grid first")?;
            let v: Vec<f64> = if name == "zero" || name == "free" {
                vec![0.0; grid.len()]
            } else {
                if !state.functions.contains_key(name) {
                    return Err(format!(
                        "QM3 POTENTIAL: no function `{name}` — define one with \
                         `DEF {name}(x, y, z) {{ ... }}`, or use `QM3 POTENTIAL zero`"
                    ));
                }
                let mut out = Vec::with_capacity(grid.len());
                for iz in 0..grid.nz {
                    for iy in 0..grid.ny {
                        for ix in 0..grid.nx {
                            let val = crate::vm::call_user_function_public(
                                name,
                                vec![
                                    Value::Num(grid.x(ix)),
                                    Value::Num(grid.y(iy)),
                                    Value::Num(grid.z(iz)),
                                ],
                                state,
                            )?;
                            match val {
                                Value::Num(y) => out.push(y),
                                other => {
                                    return Err(format!(
                                        "QM3 POTENTIAL: `{name}(x, y, z)` must return a number, \
                                         got {other}"
                                    ))
                                }
                            }
                        }
                    }
                }
                out
            };
            let lo = v.iter().cloned().fold(f64::INFINITY, f64::min);
            let hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let n = grid.len();
            state.qm3.potential = Some(v);
            state.qm3.potential_name = Some(name.clone());
            state.qm3.psi = None;
            state.qm3.time = 0.0;
            state.qm3.states = None;
            Ok(format!(
                "potential `{name}` sampled at {n} points, V in [{lo}, {hi}] (psi cleared)"
            ))
        }

        Qm3Cmd::Packet => {
            let kz = pop_num(stack)?;
            let ky = pop_num(stack)?;
            let kx = pop_num(stack)?;
            let sz = pop_num(stack)?;
            let sy = pop_num(stack)?;
            let sx = pop_num(stack)?;
            let z0 = pop_num(stack)?;
            let y0 = pop_num(stack)?;
            let x0 = pop_num(stack)?;
            let grid = state.qm3.grid.clone().ok_or("QM3 PACKET: set a grid first")?;
            let w = Wavefunction3::gaussian(grid, (x0, y0, z0), (sx, sy, sz), (kx, ky, kz))?;
            let edge = w.edge_probability(0.05);
            state.qm3.psi = Some(w);
            state.qm3.time = 0.0;
            let warn = if edge > 1e-6 {
                format!("\n  warning: {edge:.3e} already sits within 5% of a face")
            } else {
                String::new()
            };
            Ok(format!(
                "psi = 3-D Gaussian at ({x0}, {y0}, {z0}), sigma ({sx}, {sy}, {sz}), \
                 k ({kx}, {ky}, {kz}), t = 0{warn}"
            ))
        }

        Qm3Cmd::Step | Qm3Cmd::Run => {
            let (dt, steps) = if matches!(cmd, Qm3Cmd::Step) {
                (pop_num(stack)?, 1usize)
            } else {
                let n = pop_count(stack, "the step count")?;
                let t = pop_num(stack)?;
                if n == 0 {
                    return Err("QM3 RUN: the step count must be at least 1".to_string());
                }
                (t / n as f64, n)
            };
            if !dt.is_finite() || dt == 0.0 {
                return Err(format!("QM3: dt must be finite and non-zero, got {dt}"));
            }
            let ham = state.qm3.hamiltonian()?;
            let mut w = state.qm3.wavefunction()?.clone();
            let n0 = w.norm();
            Propagator3::new(ham.clone(), dt)?.run(&mut w, steps)?;
            let n1 = w.norm();
            let drift = (n1 / n0 - 1.0).abs();
            let edge = w.edge_probability(0.05);
            state.qm3.time += dt * steps as f64;
            let t = state.qm3.time;
            let e = w.energy(&ham);
            let absorbing = ham.is_absorbing();
            state.qm3.psi = Some(w);
            let warn = if edge > 1e-4 && !absorbing {
                format!("\n  warning: {edge:.3e} is within 5% of a face — the faces REFLECT")
            } else {
                String::new()
            };
            Ok(format!(
                "t = {t} ({steps} ADI step(s) of dt = {dt}), <E> = {e:.12}, \
                 norm drift = {drift:.3e}{warn}"
            ))
        }

        Qm3Cmd::Norm => Ok(format!("{:.15}", state.qm3.wavefunction()?.norm())),
        Qm3Cmd::Centroid => {
            let (x, y, z) = state.qm3.wavefunction()?.centroid();
            Ok(format!("[{x}, {y}, {z}]"))
        }
        Qm3Cmd::Energy => {
            let ham = state.qm3.hamiltonian()?;
            Ok(format!("{:.15}", state.qm3.wavefunction()?.energy(&ham)))
        }
        Qm3Cmd::Prob => {
            let zb = pop_num(stack)?;
            let za = pop_num(stack)?;
            let yb = pop_num(stack)?;
            let ya = pop_num(stack)?;
            let xb = pop_num(stack)?;
            let xa = pop_num(stack)?;
            Ok(format!(
                "{:.15}",
                state
                    .qm3
                    .wavefunction()?
                    .probability_in((xa, xb), (ya, yb), (za, zb))
            ))
        }

        Qm3Cmd::States => {
            let k = pop_count(stack, "the state count")?;
            if k == 0 {
                return Err("QM3 STATES: ask for at least one state".to_string());
            }
            let ham = state.qm3.hamiltonian()?;
            let n = ham.grid.len();
            if n > EIGEN_POINT_LIMIT {
                return Err(format!(
                    "QM3 STATES: {n} grid points is beyond what the eigensolver can do. Lanczos \
                     reorthogonalises fully and stores its whole Krylov basis, so cost grows as \
                     O(m^2 n) in time and O(m n) in memory — about {EIGEN_POINT_LIMIT} points \
                     (roughly 40^3) is the practical ceiling. Propagation has no such limit; use \
                     a coarser grid for the spectrum."
                ));
            }
            let b = ham.bound_states(k, 0)?;
            let mut s = format!(
                "{k} lowest bound state(s) — Lanczos, {} iterations{}:\n",
                b.iterations,
                if b.converged { "" } else { " (NOT converged)" }
            );
            for (i, e) in b.energies.iter().enumerate() {
                s.push_str(&format!("  E[{i}] = {e:.10}   residual {:.2e}\n", b.residuals[i]));
            }
            state.qm3.states = Some(b);
            Ok(s.trim_end().to_string())
        }

        Qm3Cmd::LoadState => {
            let n = pop_count(stack, "the state index")?;
            let ham = state.qm3.hamiltonian()?;
            let need = n + 1;
            let have = state.qm3.states.as_ref().map(|b| b.energies.len()).unwrap_or(0);
            if have < need {
                if ham.grid.len() > EIGEN_POINT_LIMIT {
                    return Err(format!(
                        "QM3 STATE: {} grid points is beyond the eigensolver's practical ceiling \
                         of {EIGEN_POINT_LIMIT}",
                        ham.grid.len()
                    ));
                }
                state.qm3.states = Some(ham.bound_states(need, 0)?);
            }
            let b = state.qm3.states.as_ref().expect("just filled");
            let e = b.energies[n];
            let mut w = Wavefunction3::new(
                ham.grid.clone(),
                b.states[n].iter().map(quantum::qm3d::real_to_complex).collect(),
            )?;
            w.normalise()?;
            state.qm3.psi = Some(w);
            state.qm3.time = 0.0;
            Ok(format!("psi = 3-D bound state {n}, E = {e:.10}, t reset to 0"))
        }

        Qm3Cmd::Absorb => {
            let power = pop_num(stack)?;
            let strength = pop_num(stack)?;
            let width = pop_num(stack)?;
            if let Some(g) = state.qm3.grid.clone() {
                let probe = Hamiltonian3::new(g.clone(), vec![0.0; g.len()], 1.0, 1.0)?;
                probe.with_absorber(width, strength, power)?;
            }
            state.qm3.absorber = Some((width, strength, power));
            state.qm3.states = None;
            Ok(format!(
                "absorbing faces on all six sides: width {width}, strength {strength}, \
                 power {power}. Propagation is no longer unitary — the norm decays by design."
            ))
        }
        Qm3Cmd::AbsorbOff => {
            state.qm3.absorber = None;
            state.qm3.states = None;
            Ok("absorbing faces removed — all six faces reflect again".to_string())
        }

        Qm3Cmd::Reset => {
            state.qm3 = Qm3State::fresh();
            Ok("3-D quantum state cleared".to_string())
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

    /// A potential of three arguments, sampled with the right axes. A
    /// NON-CUBIC grid, because an axis mix-up is invisible on a cube.
    #[test]
    fn a_three_argument_potential_is_sampled_with_correct_axes() {
        let (st, _) = run(&[
            "def v(x, y, z) { x + 10 * y + 100 * z }",
            "qm3 grid -3 3 9, -2 2 6, -1 1 4",
            "qm3 potential v",
        ]);
        let g = st.qm3.grid.as_ref().unwrap();
        assert_eq!((g.nx, g.ny, g.nz), (9, 6, 4));
        let p = st.qm3.potential.as_ref().unwrap();
        assert_eq!(p.len(), 9 * 6 * 4);
        for iz in 0..g.nz {
            for iy in 0..g.ny {
                for ix in 0..g.nx {
                    let want = g.x(ix) + 10.0 * g.y(iy) + 100.0 * g.z(iz);
                    let got = p[g.idx(ix, iy, iz)];
                    assert!(
                        (got - want).abs() < 1e-10,
                        "V({ix},{iy},{iz}) = {got}, want {want}"
                    );
                }
            }
        }
    }

    /// ADI in 3-D is exactly unitary.
    #[test]
    fn propagation_conserves_the_norm() {
        let (_, out) = run(&[
            "qm3 grid -6 6 20, -6 6 20, -6 6 20",
            "qm3 potential zero",
            "qm3 packet -2 0 0, 1 1 1, 1 0 0",
            "qm3 run 0.5 steps 50",
            "qm3 norm",
        ]);
        let norm: f64 = out[4].trim().parse().unwrap();
        assert!((norm - 1.0).abs() < 1e-10, "norm = {norm}");
    }

    /// The packet moves along the axis it was given momentum on, and
    /// only that one.
    #[test]
    fn a_packet_drifts_along_its_momentum() {
        let (st, _) = run(&[
            "qm3 grid -12 12 32, -12 12 32, -12 12 32",
            "qm3 potential zero",
            "qm3 packet -4 0 0, 1.5 1.5 1.5, 2 0 0",
            "qm3 run 1.5 steps 150",
        ]);
        let (x, y, z) = st.qm3.psi.as_ref().unwrap().centroid();
        assert!(x > -4.0 + 1.0, "<x> = {x}, should have moved");
        assert!(y.abs() < 1e-9, "<y> = {y}, must not move");
        assert!(z.abs() < 1e-9, "<z> = {z}, must not move");
    }

    /// The 3-D oscillator spectrum, with its three-fold degenerate first
    /// excited level, through the language.
    #[test]
    fn the_oscillator_spectrum_through_the_language() {
        let (_, out) = run(&[
            "def v(x, y, z) { 0.5 * (x * x + y * y + z * z) }",
            "qm3 grid -5.5 5.5 20, -5.5 5.5 20, -5.5 5.5 20",
            "qm3 potential v",
            "qm3 states 4",
        ]);
        let vals: Vec<f64> = out[3]
            .lines()
            .filter_map(|l| l.split('=').nth(1))
            .filter_map(|v| v.split_whitespace().next())
            .filter_map(|v| v.parse().ok())
            .collect();
        assert_eq!(vals.len(), 4, "expected four energies in:\n{}", out[3]);
        // E = 3/2 then 5/2 three times, up to the grid's discretisation
        assert!((vals[0] - 1.5).abs() < 0.06, "E0 = {}", vals[0]);
        for v in &vals[1..] {
            assert!((v - 2.5).abs() < 0.08, "excited level {v}, want ~2.5");
        }
        // the triplet must be degenerate with ITSELF far more tightly
        assert!(
            (vals[3] - vals[1]).abs() < 1e-8,
            "the triplet split by {}",
            vals[3] - vals[1]
        );
    }

    /// Past the eigensolver's ceiling the command must refuse up front
    /// rather than exhaust memory slowly.
    #[test]
    fn a_too_large_grid_is_refused_for_eigenstates() {
        let mut st = SimState::default();
        execute_line("qm3 grid -6 6 60, -6 6 60, -6 6 60", &mut st).unwrap();
        execute_line("qm3 potential zero", &mut st).unwrap();
        let e = execute_line("qm3 states 2", &mut st).unwrap_err();
        assert!(e.contains("ceiling") || e.contains("beyond"), "got: {e}");
        // ...but propagation on the same grid is fine
        execute_line("qm3 packet 0 0 0, 1 1 1, 1 0 0", &mut st).unwrap();
        assert!(execute_line("qm3 run 0.05 steps 2", &mut st).is_ok());
    }

    #[test]
    fn missing_prerequisites_are_reported() {
        let mut st = SimState::default();
        assert!(execute_line("qm3 energy", &mut st).unwrap_err().contains("no grid"));
        execute_line("qm3 grid -4 4 10, -4 4 10, -4 4 10", &mut st).unwrap();
        assert!(execute_line("qm3 energy", &mut st).unwrap_err().contains("no potential"));
        execute_line("qm3 potential zero", &mut st).unwrap();
        assert!(execute_line("qm3 norm", &mut st).unwrap_err().contains("no wavefunction"));
        assert!(execute_line("qm3 potential nosuch", &mut st).unwrap_err().contains("DEF"));
        assert!(execute_line("qm3 grid 4 -4 10, -4 4 10, -4 4 10", &mut st).is_err());
    }
}
