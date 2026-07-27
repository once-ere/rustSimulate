//! Bessel functions of **complex argument**, integer order.
//!
//! This closes the item deferred at the original DLMF milestone. It was
//! postponed then because the obvious route — porting AMOS/TOMS 644 —
//! needs `num-complex` and `num-traits`, and this project takes no
//! external dependencies. With [`crate::complex::Complex64`] in place
//! that objection is gone, and the algorithm turns out not to need AMOS
//! at all.
//!
//! # Miller's recurrence works unchanged for complex z
//!
//! The real-argument implementation in [`crate::bessel`] uses downward
//! recurrence normalised by an identity. **Both ingredients hold for
//! complex argument**, which is why no new algorithm is required:
//!
//! * the three-term recurrence
//!   `J_{n-1}(z) + J_{n+1}(z) = (2n/z) J_n(z)`
//!   (DLMF 10.6.1, <https://dlmf.nist.gov/10.6.E1>) is an identity in
//!   `z`, real or not;
//! * `J_n(z) ~ (z/2)^n / n!` for fixed `z`, so the wanted solution still
//!   decays with order while the unwanted one grows — the whole basis of
//!   recurring downward;
//! * the normalisation `J_0(z) + 2[J_2(z) + J_4(z) + ...] = 1` follows
//!   from the generating function `exp((z/2)(t - 1/t)) = sum J_n(z) t^n`
//!   evaluated at `t = 1`, where the left side is `exp(0) = 1`. Also an
//!   identity in `z`.
//!
//! # Accuracy, measured rather than claimed
//!
//! The limit is **cancellation in the normalising sum**. For
//! `z = x + iy` the individual `J_n(z)` grow like `exp(|y|)` while the
//! sum they must reproduce is exactly 1, so about `|y| / ln(10)` decimal
//! digits are lost out of the ~16 an `f64` carries.
//!
//! `examples/bessel_complex_accuracy.rs` measures this against the
//! generating-function identity across the plane, and the measurement is
//! what the numbers below come from:
//!
//! | `|Im z|` | 0 | 5 | 8 | 12 | 18 | 25 |
//! |---|---|---|---|---|---|---|
//! | relative error | 1e-16 | 1e-14 | 1e-13 | 1e-11 | 1e-9 | 1e-6 |
//!
//! The error depends on `|Im z|` almost independently of `Re z`, exactly
//! as the cancellation argument predicts, and the digits retained track
//! `16 - |Im z|/ln(10)` closely.
//!
//! **An earlier version of this note claimed the result was "worthless
//! past `|Im z| ~ 20`". It is not** — at `|Im z| = 25` five or six good
//! digits remain. The claim was a plausible guess written before the
//! measurement existed; the measurement corrected it. Beyond about 30 a
//! scaled or asymptotic method really is needed, and none is implemented
//! here.
//!
//! # What is here and what is not
//!
//! `J_n` and `I_n` for integer `n` and complex `z`. `Y_n` and `K_n` are
//! **not** provided: they need a logarithmic term and a different
//! algorithm, and claiming them on the strength of the same recurrence
//! would be wrong.

use crate::complex::Complex64 as C;

/// `J_0(z) ... J_{n_max}(z)` in one pass, for complex `z`.
///
/// # Errors
/// A non-finite `z`.
///
/// # Examples
/// ```
/// use special_functions::bessel_complex::bessel_j_array_c;
/// use special_functions::complex::Complex64 as C;
/// // At z = 0 only J_0 survives.
/// let j = bessel_j_array_c(4, C::ZERO).unwrap();
/// assert!((j[0].re - 1.0).abs() < 1e-15 && j[0].im.abs() < 1e-15);
/// assert!(j[1..].iter().all(|v| v.abs() < 1e-15));
/// ```
pub fn bessel_j_array_c(n_max: usize, z: C) -> Result<Vec<C>, String> {
    if !z.is_finite() {
        return Err(format!("bessel_j_array_c: z must be finite, got {z:?}"));
    }
    if z.abs() == 0.0 {
        let mut v = vec![C::ZERO; n_max + 1];
        v[0] = C::ONE;
        return Ok(v);
    }

    // Same seeding rule as the real case, driven by |z|: the start must
    // sit well above BOTH the wanted order and |z|, because J_n(z) only
    // begins to decay once n exceeds |z|.
    let az = z.abs();
    let start = n_max + 30 + (1.5 * az + 12.0 * az.sqrt()) as usize;

    let mut jp1 = C::ZERO; // J_{k+1}
    let mut j = C::new(1.0e-290, 0.0); // J_k, arbitrary seed
    let mut out = vec![C::ZERO; n_max + 1];
    // Accumulates J_0 + 2(J_2 + J_4 + ...) in the unnormalised scale.
    let mut sum = C::ZERO;
    let inv_z = z.inv();

    for k in (0..=start).rev() {
        let jm1 = inv_z * j * (2 * (k + 1)) as f64 - jp1;
        jp1 = j;
        j = jm1;

        if j.abs() > 1.0e250 {
            let s = 1.0e-250;
            j = j * s;
            jp1 = jp1 * s;
            sum = sum * s;
            for e in out.iter_mut() {
                *e = *e * s;
            }
        }
        if k <= n_max {
            out[k] = j;
        }
        if k.is_multiple_of(2) {
            sum = sum + if k == 0 { j } else { j * 2.0 };
        }
    }

    let s = sum.abs();
    if s == 0.0 || !sum.is_finite() {
        return Err(format!(
            "bessel_j_array_c: normalisation failed for n_max={n_max}, z={z:?}"
        ));
    }
    let scale = sum.inv();
    for e in out.iter_mut() {
        *e = *e * scale;
    }
    Ok(out)
}

/// A single `J_n(z)` for integer `n >= 0` and complex `z`.
///
/// # Errors
/// Negative order, or a non-finite `z`.
///
/// # Examples
/// ```
/// use special_functions::bessel_complex::bessel_j_c;
/// use special_functions::complex::Complex64 as C;
/// // Real argument must reproduce the real routine.
/// let v = bessel_j_c(0, C::real(2.404_825_557_695_773)).unwrap();
/// assert!(v.abs() < 1e-12);
/// ```
pub fn bessel_j_c(n: i32, z: C) -> Result<C, String> {
    if n < 0 {
        return Err(format!("bessel_j_c: order n must be >= 0, got {n}"));
    }
    Ok(bessel_j_array_c(n as usize, z)?[n as usize])
}

/// The modified Bessel function `I_n(z)` for complex `z`.
///
/// Obtained from `I_n(z) = i^{-n} J_n(i z)` (DLMF 10.27.6), which is an
/// identity rather than a separate algorithm — the whole reason a
/// complex `J` is worth having.
///
/// # Errors
/// As [`bessel_j_c`].
///
/// # Examples
/// ```
/// use special_functions::bessel_complex::bessel_i_c;
/// use special_functions::complex::Complex64 as C;
/// // I_0(0) = 1
/// let v = bessel_i_c(0, C::ZERO).unwrap();
/// assert!((v.re - 1.0).abs() < 1e-15);
/// ```
pub fn bessel_i_c(n: i32, z: C) -> Result<C, String> {
    if n < 0 {
        return Err(format!("bessel_i_c: order n must be >= 0, got {n}"));
    }
    let j = bessel_j_c(n, C::I * z)?;
    Ok(j * i_pow(-n))
}

/// `i^k` for any integer `k`, by cycling rather than by `powi` — exact,
/// with no rounding at all.
fn i_pow(k: i32) -> C {
    match k.rem_euclid(4) {
        0 => C::ONE,
        1 => C::I,
        2 => C::real(-1.0),
        _ => C::new(0.0, -1.0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bessel::bessel_j_array;
    use crate::rel_err;
    use spec_math::cephes64::{i0, i1, jv};

    fn close(a: C, b: C, tol: f64) -> bool {
        (a - b).abs() <= tol * (1.0_f64).max(a.abs().max(b.abs()))
    }

    /// On the real axis the complex routine must reproduce the real one
    /// exactly enough to be interchangeable — and both must agree with
    /// the independently written vendored Cephes.
    #[test]
    fn real_axis_matches_the_real_routine_and_cephes() {
        for &x in &[0.1, 0.7, 2.5, 5.0, 9.0, 20.0, 45.0] {
            let cx = bessel_j_array_c(12, C::real(x)).unwrap();
            let rx = bessel_j_array(12, x).unwrap();
            for n in 0..=12 {
                assert!(
                    cx[n].im.abs() < 1e-14,
                    "J_{n}({x}) should be real, got imaginary part {}",
                    cx[n].im
                );
                assert!(
                    rel_err(cx[n].re, rx[n]) < 1e-10 || (cx[n].re - rx[n]).abs() < 1e-15,
                    "J_{n}({x}): complex {} vs real {}",
                    cx[n].re,
                    rx[n]
                );
                let ceph = jv(n as f64, x);
                let tol = if ceph.abs() < 1e-14 { 1e-6 } else { 1e-10 };
                assert!(
                    rel_err(cx[n].re, ceph) < tol || (cx[n].re - ceph).abs() < 1e-15,
                    "J_{n}({x}): ours {} vs cephes {ceph}",
                    cx[n].re
                );
            }
        }
    }

    /// **The generating function is an absolute test.** For any complex
    /// `z` and any `t != 0`,
    ///
    /// ```text
    ///     exp((z/2)(t - 1/t)) = sum_{n=-inf}^{inf} J_n(z) t^n
    /// ```
    ///
    /// with `J_{-n} = (-1)^n J_n`. The right-hand side is built entirely
    /// from the values under test and the left-hand side from `exp`, so
    /// no reference implementation is involved at all.
    #[test]
    fn the_generating_function_identity_holds_off_the_real_axis() {
        let zs = [
            C::new(1.0, 1.0),
            C::new(3.0, -2.0),
            C::new(-4.0, 1.5),
            C::new(0.5, 5.0),
            C::new(8.0, 3.0),
            C::new(-2.0, -6.0),
        ];
        let ts = [C::new(1.3, 0.0), C::new(0.7, 0.4), C::new(-1.1, 0.9)];
        for z in zs {
            let n_max = 60;
            let j = bessel_j_array_c(n_max, z).unwrap();
            for t in ts {
                let mut sum = j[0];
                let mut tp = C::ONE; // t^n
                let mut tm = C::ONE; // t^-n
                let t_inv = t.inv();
                for (n, jn) in j.iter().enumerate().take(n_max + 1).skip(1) {
                    tp = tp * t;
                    tm = tm * t_inv;
                    let sign = if n.is_multiple_of(2) { 1.0 } else { -1.0 };
                    // J_n t^n + J_{-n} t^{-n}
                    sum = sum + *jn * tp + *jn * tm * sign;
                }
                let want = ((z * 0.5) * (t - t_inv)).exp();
                assert!(
                    close(sum, want, 1e-9),
                    "z = {z:?}, t = {t:?}: sum {sum:?} vs exp {want:?}"
                );
            }
        }
    }

    /// The addition theorem `J_n(z+w) = sum_k J_k(z) J_{n-k}(w)`
    /// (DLMF 10.23.1) — another absolute identity, and one that couples
    /// different arguments so a systematic scale error cannot hide.
    #[test]
    fn the_addition_theorem_holds() {
        let cases = [
            (C::new(1.5, 0.8), C::new(-0.7, 1.1)),
            (C::new(3.0, -1.0), C::new(2.0, 2.0)),
            (C::new(-2.5, 0.0), C::new(1.0, -3.0)),
        ];
        for (z, w) in cases {
            let m = 70;
            let jz = bessel_j_array_c(m, z).unwrap();
            let jw = bessel_j_array_c(m, w).unwrap();
            let jsum = bessel_j_array_c(6, z + w).unwrap();
            let jn = |arr: &[C], k: i32| -> C {
                let a = k.unsigned_abs() as usize;
                if a >= arr.len() {
                    C::ZERO
                } else if k >= 0 || a.is_multiple_of(2) {
                    arr[a]
                } else {
                    arr[a] * -1.0
                }
            };
            for n in 0..=4i32 {
                let mut acc = C::ZERO;
                for k in -(m as i32)..=(m as i32) {
                    acc = acc + jn(&jz, k) * jn(&jw, n - k);
                }
                assert!(
                    close(acc, jsum[n as usize], 1e-9),
                    "n = {n}, z = {z:?}, w = {w:?}: {acc:?} vs {:?}",
                    jsum[n as usize]
                );
            }
        }
    }

    /// On the imaginary axis `J_n(iy) = i^n I_n(y)`. The vendored
    /// Cephes provides `I_0` and `I_1` only — no general `I_v` — so the
    /// cross-check against independent machinery covers those two
    /// orders, and the modified generating function below covers the
    /// rest without needing a reference at all.
    #[test]
    fn the_imaginary_axis_matches_the_vendored_modified_bessel() {
        for &y in &[0.3, 1.0, 2.5, 5.0, 9.0] {
            let j = bessel_j_array_c(4, C::new(0.0, y)).unwrap();
            for (n, want_real) in [(0usize, i0(y)), (1usize, i1(y))] {
                let want = i_pow(n as i32) * want_real;
                assert!(
                    close(j[n], want, 1e-10),
                    "J_{n}(i*{y}) = {:?}, want i^{n} I_{n}({y}) = {want:?}",
                    j[n]
                );
            }
        }
    }

    /// **The modified generating function**, absolute and reference-free:
    /// `exp((z/2)(t + 1/t)) = sum_n I_n(z) t^n` with `I_{-n} = I_n`
    /// (DLMF 10.35.1). This covers `I_n` at every order, which the
    /// vendored library cannot.
    #[test]
    fn the_modified_generating_function_holds() {
        let zs = [
            C::new(1.5, 0.0),
            C::new(2.0, 1.0),
            C::new(-1.0, 2.0),
            C::new(3.0, -1.5),
        ];
        let ts = [C::new(1.4, 0.0), C::new(0.8, 0.5)];
        for z in zs {
            let n_max = 50;
            // I_n(z) = i^-n J_n(i z), taken from one J pass
            let j = bessel_j_array_c(n_max, C::I * z).unwrap();
            let i_of: Vec<C> = (0..=n_max).map(|n| j[n] * i_pow(-(n as i32))).collect();
            for t in ts {
                let t_inv = t.inv();
                let mut sum = i_of[0];
                let (mut tp, mut tm) = (C::ONE, C::ONE);
                for inv in i_of.iter().take(n_max + 1).skip(1) {
                    tp = tp * t;
                    tm = tm * t_inv;
                    // I_{-n} = I_n, so no alternating sign here
                    sum = sum + *inv * tp + *inv * tm;
                }
                let want = ((z * 0.5) * (t + t_inv)).exp();
                assert!(
                    close(sum, want, 1e-9),
                    "z = {z:?}, t = {t:?}: sum {sum:?} vs exp {want:?}"
                );
            }
        }
    }

    /// `I_n(z)` for complex `z` must reduce to the vendored real `I_n`
    /// on the real axis — the identity route is only worth having if it
    /// lands in the right place.
    #[test]
    fn modified_bessel_reduces_correctly_on_the_real_axis() {
        for &x in &[0.2, 1.0, 3.0, 7.0] {
            for (n, want) in [(0i32, i0(x)), (1i32, i1(x))] {
                let got = bessel_i_c(n, C::real(x)).unwrap();
                assert!(got.im.abs() < 1e-10, "I_{n}({x}) should be real, got {got:?}");
                assert!(
                    rel_err(got.re, want) < 1e-9,
                    "I_{n}({x}) = {} vs cephes {want}",
                    got.re
                );
            }
        }
    }

    /// Symmetries: `J_n(conj z) = conj J_n(z)` because the coefficients
    /// are real, and `J_n(-z) = (-1)^n J_n(z)`.
    #[test]
    fn conjugation_and_parity() {
        for z in [C::new(2.0, 3.0), C::new(-1.5, 0.7), C::new(4.0, -2.5)] {
            let a = bessel_j_array_c(6, z).unwrap();
            let b = bessel_j_array_c(6, z.conj()).unwrap();
            let c = bessel_j_array_c(6, z * -1.0).unwrap();
            for n in 0..=6 {
                assert!(close(b[n], a[n].conj(), 1e-11), "conjugation at n = {n}");
                let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
                assert!(close(c[n], a[n] * sign, 1e-11), "parity at n = {n}");
            }
        }
    }

    /// The three-term recurrence, satisfied by the returned values.
    #[test]
    fn the_recurrence_is_satisfied() {
        for z in [C::new(1.0, 2.0), C::new(-3.0, 1.0), C::new(6.0, -4.0)] {
            let j = bessel_j_array_c(25, z).unwrap();
            let inv = z.inv();
            for n in 1..24 {
                let lhs = j[n - 1] + j[n + 1];
                let rhs = inv * j[n] * (2 * n) as f64;
                assert!(close(lhs, rhs, 1e-10), "recurrence n = {n}, z = {z:?}");
            }
        }
    }

    /// The documented accuracy limit is cancellation in the normalising
    /// sum, losing roughly `|Im z| / ln 10` digits. Pinned as a
    /// MEASUREMENT: agreement must hold at `|Im z| = 8` and is allowed
    /// to be much worse at `|Im z| = 25`, which is the honest statement
    /// of where the method stops working.
    #[test]
    fn accuracy_degrades_with_imaginary_part_as_documented() {
        let err_at = |z: C| -> f64 {
            let n_max = 80;
            let j = bessel_j_array_c(n_max, z).unwrap();
            let t = C::new(1.2, 0.3);
            let t_inv = t.inv();
            let mut sum = j[0];
            let (mut tp, mut tm) = (C::ONE, C::ONE);
            for (n, jn) in j.iter().enumerate().take(n_max + 1).skip(1) {
                tp = tp * t;
                tm = tm * t_inv;
                let sign = if n.is_multiple_of(2) { 1.0 } else { -1.0 };
                sum = sum + *jn * tp + *jn * tm * sign;
            }
            let want = ((z * 0.5) * (t - t_inv)).exp();
            (sum - want).abs() / want.abs().max(1.0)
        };
        let good = err_at(C::new(2.0, 8.0));
        assert!(good < 1e-8, "at |Im z| = 8 the error is {good}, expected < 1e-8");
        // and the degradation is real, not imagined
        let bad = err_at(C::new(2.0, 25.0));
        assert!(
            bad > good,
            "error should grow with |Im z|: {good} at 8 vs {bad} at 25"
        );
    }

    #[test]
    fn invalid_input_is_reported() {
        assert!(bessel_j_c(-1, C::ONE).is_err(), "negative order");
        assert!(bessel_i_c(-1, C::ONE).is_err(), "negative order");
        assert!(
            bessel_j_array_c(3, C::new(f64::NAN, 0.0)).is_err(),
            "non-finite z"
        );
        assert!(
            bessel_j_array_c(3, C::new(0.0, f64::INFINITY)).is_err(),
            "non-finite z"
        );
    }
}
