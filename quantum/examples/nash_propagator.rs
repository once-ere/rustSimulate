//! The Nash propagator, measured against things that can contradict it.
//!
//! `EVOLVE_NASH` is the one piece of original numerical work in the
//! SolveIt C++: a split-operator scheme whose kinetic half is applied
//! by a **Bessel stencil** rather than by an FFT or a linear solve. This
//! example does not describe it — it measures it, three ways:
//!
//! 1. against a **closed form**, on a plane wave, where the scheme is
//!    exact and the answer is a single known phase;
//! 2. against **diagonalisation**, with a potential switched on, where
//!    the Lie–Trotter splitting is the only error left;
//! 3. against **its own error bound**, for the Bessel truncation.
//!
//! Run: cargo run -p quantum --release --example nash_propagator

use quantum::nash::{norm, order_for, truncation_bound, NashPropagator, PeriodicGrid};
use special_functions::complex::Complex64 as C;
use special_functions::eigen::jacobi_eigen;

fn packet(grid: &PeriodicGrid, x0: f64, k0: f64, width: f64) -> Vec<C> {
    let mut psi: Vec<C> = (0..grid.n)
        .map(|i| {
            let x = grid.x(i);
            C::from_polar((-((x - x0) / width).powi(2) / 2.0).exp(), k0 * x)
        })
        .collect();
    let s = norm(&psi, grid.h()).sqrt();
    for z in &mut psi {
        *z = *z * (1.0 / s);
    }
    psi
}

/// `exp(-i H t / hbar)` by diagonalising the periodic Hamiltonian.
/// Shares no code with the propagator.
fn exact(grid: &PeriodicGrid, v: &[f64], t: f64, psi0: &[C]) -> Vec<C> {
    let n = grid.n;
    let h = grid.h();
    let kappa = 1.0 / (2.0 * h * h);
    let mut a = vec![vec![0.0; n]; n];
    for i in 0..n {
        a[i][i] = 2.0 * kappa + v[i];
        a[i][(i + 1) % n] += -kappa;
        a[i][(i + n - 1) % n] += -kappa;
    }
    let (e, vecs) = jacobi_eigen(&a).unwrap();
    let mut out = vec![C::ZERO; n];
    for (m, phi) in vecs.iter().enumerate() {
        let mut c = C::ZERO;
        for (i, &p) in phi.iter().enumerate() {
            c = c + psi0[i] * p;
        }
        let ph = C::from_polar(1.0, -e[m] * t);
        for (i, &p) in phi.iter().enumerate() {
            out[i] = out[i] + c * ph * p;
        }
    }
    out
}

fn max_diff(a: &[C], b: &[C]) -> f64 {
    a.iter().zip(b).map(|(p, q)| (*p - *q).abs()).fold(0.0, f64::max)
}

fn main() {
    println!("The Nash propagator: exp(-i H dt) as a Bessel stencil\n");
    println!("    psi_j <- e^(-i L (1+v_j)) [ J_0(L) psi_j");
    println!("                + sum_M i^M J_M(L) (psi_(j-M) + psi_(j+M)) ]\n");
    println!("with L = hbar dt / (m h^2). One pass over the grid, no solver.\n");

    // -----------------------------------------------------------------
    println!("1. Free particle: the scheme is EXACT, so this is a closed form.\n");
    println!("With V = 0 the two factors commute and no splitting error exists.");
    println!("A plane wave exp(i k x) must be multiplied by exactly");
    println!("exp(-i L (1 - cos k h)) — the lattice dispersion, not the continuum one.\n");
    let grid = PeriodicGrid::new(0.0, 1.0, 64).unwrap();
    let free = vec![0.0; grid.n];
    let p = NashPropagator::new(grid.clone(), &free, 1.0, 1.0, 2.0e-4, None).unwrap();
    println!("   lambda = {:.4},  order = {},  truncation <= {:.1e}\n",
             p.lambda(), p.order(), p.truncation_error());
    println!("   {:>6} {:>14} {:>14}", "mode", "|error|", "phase applied");
    for m in [1_i32, 3, 7, 16, 31] {
        let k = 2.0 * std::f64::consts::PI * f64::from(m);
        let mut psi: Vec<C> = (0..grid.n).map(|i| C::from_polar(1.0, k * grid.x(i))).collect();
        let before = psi.clone();
        p.step(&mut psi).unwrap();
        let want = C::from_polar(1.0, -p.lambda() * (1.0 - (k * grid.h()).cos()));
        let got: Vec<C> = before.iter().map(|z| *z * want).collect();
        println!("   {m:>6} {:>14.2e} {:>14.6}", max_diff(&psi, &got), want.arg());
    }

    // -----------------------------------------------------------------
    println!("\n\n2. With a potential: what is left is the SPLITTING, and it is O(dt).\n");
    println!("Harmonic well on [-6, 6], measured against diagonalising H.\n");
    let grid = PeriodicGrid::new(-6.0, 6.0, 48).unwrap();
    let v: Vec<f64> = grid.points().iter().map(|x| 0.5 * x * x).collect();
    let psi0 = packet(&grid, -1.0, 2.0, 0.8);
    let t = 0.05;
    let want = exact(&grid, &v, t, &psi0);
    println!("   {:>8} {:>12} {:>12} {:>10} {:>16}", "steps", "dt", "error", "ratio", "norm drift");
    let mut prev = f64::NAN;
    for steps in [25_usize, 50, 100, 200, 400] {
        let dt = t / steps as f64;
        let p = NashPropagator::new(grid.clone(), &v, 1.0, 1.0, dt, None).unwrap();
        let mut psi = psi0.clone();
        p.run(&mut psi, steps).unwrap();
        let e = max_diff(&psi, &want);
        let drift = (norm(&psi, grid.h()) - 1.0).abs();
        if prev.is_finite() {
            println!("   {steps:>8} {dt:>12.2e} {e:>12.3e} {:>10.2} {drift:>16.2e}", prev / e);
        } else {
            println!("   {steps:>8} {dt:>12.2e} {e:>12.3e} {:>10} {drift:>16.2e}", "-");
        }
        prev = e;
    }
    println!("\n   The ratio is 2: halving dt halves the error, which is Lie-Trotter.");
    println!("   The norm does NOT drift with dt — each factor is unitary whatever");
    println!("   the step size, so accuracy and stability fail independently here.");
    println!("   A Strang arrangement would make this second order for one extra");
    println!("   pointwise multiply; the original is Lie and the port is faithful.");

    // -----------------------------------------------------------------
    println!("\n\n3. The Bessel truncation, which is NOT what limits the scheme.\n");
    println!("Truncating Jacobi-Anger at K leaves 2 sum_(M>K) |J_M(L)|, and");
    println!("|J_M(L)| <= (L/2)^M / M!, so it falls superexponentially past K ~ L.\n");
    println!("   {:>8} {:>12} {:>12} {:>12} {:>12}", "K", "L=0.92", "L=3", "L=8", "L=20");
    for k in [2_usize, 4, 6, 8, 12, 16, 24, 32] {
        print!("   {k:>8}");
        for l in [0.92_f64, 3.0, 8.0, 20.0] {
            print!(" {:>12.1e}", truncation_bound(l, k).unwrap());
        }
        println!();
    }
    println!("\n   SolveIt shipped L = 0.92 with K = 16. The order actually needed");
    println!("   for the stencil to be exact to rounding is {}, so the shipped value",
             order_for(0.92, f64::EPSILON).unwrap());
    println!("   carried a small margin over it — harmless, and worth knowing");
    println!("   rather than guessing at.");
    println!("\n   Read the table the other way and it is a step-size limit: L grows");
    println!("   like dt/h^2, so refining the grid at fixed dt costs stencil width.");

    // -----------------------------------------------------------------
    println!("\n\nPractical guidance:");
    println!("  boundaries    periodic — a packet leaving one edge re-enters at the");
    println!("                other. qm1d's Grid is Dirichlet; they are not the same");
    println!("                domain and results are not interchangeable");
    println!("  accuracy      first order in dt. Halve dt to halve the error");
    println!("  stability     unconditional: norm is conserved to rounding at any dt,");
    println!("                so a wrong answer here stays a normalised wrong answer");
    println!("  cost          O(n K) per step and no solver, against Crank-Nicolson's");
    println!("                O(n) with a tridiagonal solve");
    println!("  choosing K    leave it to order_for; it is linear in cost and");
    println!("                superexponential in accuracy, so there is nothing to save");
}
