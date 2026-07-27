//! Large order: where a 1/z expansion stops and a 1/nu expansion starts.
//!
//! Stage 14 left one hole, and this example is what located it: `J_nu(z)`
//! for `z` below `nu`. There `J` is exponentially small while the Hankel
//! functions it was being built from are exponentially large, so the
//! subtraction destroyed everything. At `nu = 400.5, z = 240` the answer
//! came back wrong by a factor of `5e89` — with a claimed error of zero,
//! because at `nu = 1/2` the 1/z expansion terminates exactly and its
//! truncation estimate was literally 0.
//!
//! The remedy is an expansion in `1/nu` instead: the Debye polynomials
//! and DLMF 10.19 / 10.41, which produce the small number directly.
//!
//! The turning point `z ~ nu` needed nothing — it was already at 1e-14,
//! which is why the Airy-type expansion of DLMF 10.20 is not here. That
//! claim is the third table below.
//!
//! Run: cargo run -p special_functions --release --example large_order_accuracy

use special_functions::bessel_scaled::{
    bessel_j_scaled_nu, bessel_k_scaled_nu, bessel_y_scaled_nu,
};
use special_functions::complex::Complex64 as C;
use special_functions::debye::jy_debye;
use spec_math::cephes64::{jv, yv};

const FRACS: [f64; 9] = [0.1, 0.3, 0.5, 0.7, 0.85, 0.95, 1.0, 1.2, 2.0];

fn cell(v: Result<C, String>, want: f64) -> String {
    match v {
        Ok(g) if want != 0.0 && want.is_finite() => {
            format!("{:>10.1e}", (g.re - want).abs() / want.abs())
        }
        Ok(_) => format!("{:>10}", "no ref"),
        Err(e) if e.contains("outside f64") => format!("{:>10}", "no f64"),
        Err(_) => format!("{:>10}", "REFUSED"),
    }
}

fn main() {
    println!("Large order, real axis. Relative error against Cephes jv/yv.\n");
    println!("  'no f64'  = the value is determined but outside f64 range");
    println!("  'no ref'  = Cephes cannot represent it");
    println!("  'REFUSED' = no method here\n");

    for (name, f) in [
        ("J", bessel_j_scaled_nu as fn(f64, C) -> Result<C, String>),
        ("Y", bessel_y_scaled_nu),
    ] {
        println!("{name}_nu(x), columns are x/nu:\n");
        print!("{:>8} |", "nu");
        for fr in FRACS {
            print!(" {fr:>9.2}");
        }
        println!();
        println!("  -------+{}", "-".repeat(10 * FRACS.len()));
        for nu in [10.5f64, 40.5, 100.5, 200.5, 400.5, 1000.5] {
            print!("{nu:>8.1} |");
            for fr in FRACS {
                let x = nu * fr;
                let want = if name == "J" { jv(nu, x) } else { yv(nu, x) };
                print!(" {}", cell(f(nu, C::real(x)), want));
            }
            println!();
        }
        println!();
    }

    println!("Every catastrophic entry is gone. Before this stage the J row at");
    println!("nu = 400.5 read 5.4e89 at x/nu = 0.6 and the nu = 1000.5 row read");
    println!("2.4e246 — numbers returned with confident small error estimates.\n");

    // -----------------------------------------------------------------
    println!("\nWhy the Airy-type expansion of DLMF 10.20 is not here.\n");
    println!("It is the uniform expansion THROUGH the turning point x = 1. So the");
    println!("question is whether anything at x ~ 1 needs fixing. Measured:\n");
    println!("{:>8} {:>10} {:>10} {:>10} {:>10} {:>10}", "nu", "x/nu=0.95", "0.98", "1.00", "1.02", "1.10");
    for nu in [100.5f64, 200.5, 400.5, 1000.5] {
        print!("{nu:>8.1}");
        for fr in [0.95f64, 0.98, 1.0, 1.02, 1.10] {
            let x = nu * fr;
            print!(" {}", cell(bessel_j_scaled_nu(nu, C::real(x)), jv(nu, x)));
        }
        println!();
    }
    println!("\nNothing there is worse than 1e-10, and most is 1e-14. Olver's");
    println!("expansion truncated at A_0, B_0 carries a relative error of O(nu^-2)");
    println!("— about 6e-6 at nu = 400 — so adding it would LOWER the floor in the");
    println!("region it covers. It is the right tool for a library whose target is");
    println!("1e-6; it is the wrong tool for one holding 1e-13 everywhere else.\n");
    println!("What it would genuinely add is coverage where the value is outside");
    println!("f64 anyway (the 'no f64' cells above), which no expansion can fix");
    println!("without a different number type.");

    // -----------------------------------------------------------------
    println!("\n\nThe Debye truncation estimate, which is what chooses the method.\n");
    println!("It has to grow towards the turning point, or the caller cannot know");
    println!("to stop trusting it. nu = 400:\n");
    println!("{:>10} {:>14}", "x/nu", "estimate");
    for fr in [0.1f64, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99] {
        let e = jy_debye(400.0, 400.0 * fr)
            .0
            .map(|u| u.err)
            .unwrap_or(f64::INFINITY);
        println!("{fr:>10.2} {e:>14.1e}");
    }

    // -----------------------------------------------------------------
    println!("\n\nAnd the modified pair, from the same polynomials (DLMF 10.41).");
    println!("exp(x) K_nu(x) at orders where the vendored Cephes kn overflows:\n");
    println!("{:>8} {:>16} {:>16}", "nu", "x = nu", "x = 10 nu");
    for nu in [50.0f64, 200.0, 800.0, 3000.0] {
        let a = bessel_k_scaled_nu(nu, C::real(nu)).map(|v| v.re);
        let b = bessel_k_scaled_nu(nu, C::real(10.0 * nu)).map(|v| v.re);
        let f = |v: Result<f64, String>| match v {
            Ok(x) => format!("{x:>16.6e}"),
            Err(_) => format!("{:>16}", "out of range"),
        };
        println!("{nu:>8.0} {} {}", f(a), f(b));
    }
    println!("\n(pinned by the I-K Wronskian, whose right-hand side is elementary —");
    println!(" there is no reference implementation to compare against up here)");
}
