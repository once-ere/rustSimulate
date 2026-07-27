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
//! # `Y_n` and `K_n` need a different method
//!
//! `Y_n` has a **logarithmic branch point** at the origin, so no amount
//! of recurrence produces it from `J_n` alone — the earlier version of
//! this module said as much and left them out. They are here now, by the
//! route that difficulty actually dictates:
//!
//! * `Y_0` and `Y_1` from the ascending series (DLMF 10.8.1), which
//!   carries the `ln(z/2) J_n(z)` term and digamma coefficients;
//! * higher orders by **upward** recurrence, which for `Y` is the stable
//!   direction — precisely the opposite of `J`, because `Y_n` *grows*
//!   with order while `J_n` decays. Using the same downward sweep for
//!   both would destroy one of them;
//! * `K_n` by the identity
//!   `K_n(z) = (pi/2) i^{n+1} [J_n(iz) + i Y_n(iz)]` (DLMF 10.27.8),
//!   so it needs no third algorithm.
//!
//! **The branch cut matters.** `Y_n` and `K_n` inherit the cut of `ln`
//! along the negative real axis and are discontinuous across it. `J_n`
//! and `I_n` are entire and have no such restriction.

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

/// `psi(m)` for a positive integer, exactly: `psi(1) = -gamma` and
/// `psi(m+1) = psi(m) + 1/m`.
fn digamma_int(m: usize) -> f64 {
    // Euler-Mascheroni, to full double precision.
    const GAMMA: f64 = 0.577_215_664_901_532_9;
    let mut v = -GAMMA;
    for k in 1..m {
        v += 1.0 / k as f64;
    }
    v
}

/// `Y_n(z)` for `n = 0` or `n = 1` by the ascending series, DLMF 10.8.1:
///
/// ```text
///  Y_n(z) = -(1/pi) (z/2)^-n sum_{k=0}^{n-1} (n-k-1)!/k! (z^2/4)^k
///         + (2/pi) ln(z/2) J_n(z)
///         - (1/pi) (z/2)^n sum_{k>=0} [psi(k+1)+psi(n+k+1)] (-z^2/4)^k
///                                      / (k! (n+k)!)
/// ```
fn y_series(n: usize, z: C) -> Result<C, String> {
    let half = z * 0.5;
    let q = half * half; // (z/2)^2
    let jn = bessel_j_c(n as i32, z)?;
    let inv_pi = 1.0 / std::f64::consts::PI;

    // finite sum, empty for n = 0
    let mut finite = C::ZERO;
    if n >= 1 {
        let mut qk = C::ONE; // (z^2/4)^k
        let mut fact_k = 1.0f64;
        for k in 0..n {
            let coeff = factorial(n - k - 1) / fact_k;
            finite = finite + qk * coeff;
            qk = qk * q;
            fact_k *= (k + 1) as f64;
        }
        // multiply by (z/2)^-n
        let mut p = C::ONE;
        for _ in 0..n {
            p = p * half;
        }
        finite = finite * p.inv();
    }

    // infinite sum
    let mut infinite = C::ZERO;
    let mut term_pow = C::ONE; // (-z^2/4)^k
    let neg_q = q * -1.0;
    let mut fact_k = 1.0f64;
    let mut fact_nk = factorial(n);
    for k in 0..200 {
        let coeff = (digamma_int(k + 1) + digamma_int(n + k + 1)) / (fact_k * fact_nk);
        let add = term_pow * coeff;
        infinite = infinite + add;
        if k > 4 && add.abs() <= 1e-18 * infinite.abs().max(1e-300) {
            break;
        }
        term_pow = term_pow * neg_q;
        fact_k *= (k + 1) as f64;
        fact_nk *= (n + k + 1) as f64;
    }
    let mut p = C::ONE;
    for _ in 0..n {
        p = p * half;
    }
    infinite = infinite * p;

    Ok(finite * -inv_pi + half.ln() * jn * (2.0 * inv_pi) - infinite * inv_pi)
}

/// `k!` as an `f64`. Only ever called with small `k` here.
fn factorial(k: usize) -> f64 {
    (1..=k).map(|i| i as f64).product::<f64>().max(1.0)
}

/// `Y_0(z) ... Y_{n_max}(z)`, complex argument, integer order.
///
/// # Errors
/// A non-finite `z`, or `z == 0` where `Y` is infinite.
///
/// # Examples
/// ```
/// use special_functions::bessel_complex::bessel_y_array_c;
/// use special_functions::complex::Complex64 as C;
/// // Y_0(1) is about 0.08825696
/// let y = bessel_y_array_c(1, C::real(1.0)).unwrap();
/// assert!((y[0].re - 0.088_256_964_215_676_96).abs() < 1e-10);
/// ```
pub fn bessel_y_array_c(n_max: usize, z: C) -> Result<Vec<C>, String> {
    if !z.is_finite() {
        return Err(format!("bessel_y_array_c: z must be finite, got {z:?}"));
    }
    if z.abs() == 0.0 {
        return Err(
            "bessel_y_array_c: Y_n has a logarithmic singularity at z = 0".to_string()
        );
    }
    let y0 = y_series(0, z)?;
    if n_max == 0 {
        return Ok(vec![y0]);
    }
    let y1 = y_series(1, z)?;
    let mut out = Vec::with_capacity(n_max + 1);
    out.push(y0);
    out.push(y1);
    // Upward recurrence is the STABLE direction for Y, because Y_n grows
    // with order. (For J it is the unstable one — the two functions need
    // opposite sweeps, which is the whole reason they cannot share an
    // implementation.)
    let inv_z = z.inv();
    for n in 1..n_max {
        let next = inv_z * out[n] * (2 * n) as f64 - out[n - 1];
        out.push(next);
    }
    Ok(out)
}

/// A single `Y_n(z)`.
///
/// # Errors
/// Negative order, non-finite `z`, or `z == 0`.
pub fn bessel_y_c(n: i32, z: C) -> Result<C, String> {
    if n < 0 {
        return Err(format!("bessel_y_c: order n must be >= 0, got {n}"));
    }
    Ok(bessel_y_array_c(n as usize, z)?[n as usize])
}

/// The modified Bessel function of the second kind, `K_n(z)`, complex
/// argument.
///
/// From `K_n(z) = (pi/2) i^{n+1} [J_n(iz) + i Y_n(iz)]` (DLMF 10.27.8) —
/// the Hankel function `H^(1)_n` evaluated on the rotated argument, so
/// no third algorithm is needed.
///
/// # Errors
/// Negative order, non-finite `z`, or `z == 0`.
///
/// # Examples
/// ```
/// use special_functions::bessel_complex::bessel_k_c;
/// use special_functions::complex::Complex64 as C;
/// // K_0(1) is about 0.4210244382
/// let k = bessel_k_c(0, C::real(1.0)).unwrap();
/// assert!((k.re - 0.421_024_438_240_708_3).abs() < 1e-9);
/// ```
pub fn bessel_k_c(n: i32, z: C) -> Result<C, String> {
    if n < 0 {
        return Err(format!("bessel_k_c: order n must be >= 0, got {n}"));
    }
    if !z.is_finite() {
        return Err(format!("bessel_k_c: z must be finite, got {z:?}"));
    }
    if z.abs() == 0.0 {
        return Err("bessel_k_c: K_n has a singularity at z = 0".to_string());
    }
    let iz = C::I * z;
    let j = bessel_j_c(n, iz)?;
    let y = bessel_y_c(n, iz)?;
    Ok((j + C::I * y) * i_pow(n + 1) * (std::f64::consts::PI * 0.5))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bessel::bessel_j_array;
    use crate::rel_err;
    use spec_math::cephes64::{i0, i1, jv, k0, k1, yn};

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

    /// **The Wronskian is the absolute test for `Y`.**
    ///
    /// `J_{n+1}(z) Y_n(z) - J_n(z) Y_{n+1}(z) = 2/(pi z)` (DLMF 10.5.2).
    /// The right-hand side is elementary, so this pins `Y` in both scale
    /// and phase against a `J` that is itself independently verified —
    /// no table, no reference library.
    #[test]
    fn the_j_y_wronskian_holds() {
        let zs = [
            C::real(0.7),
            C::real(4.0),
            C::new(1.0, 1.0),
            C::new(3.0, -2.0),
            C::new(-2.5, 1.5),
            C::new(0.4, 3.0),
            C::new(6.0, 2.0),
        ];
        for z in zs {
            let j = bessel_j_array_c(7, z).unwrap();
            let y = bessel_y_array_c(7, z).unwrap();
            let want = z.inv() * (2.0 / std::f64::consts::PI);
            for n in 0..6 {
                let w = j[n + 1] * y[n] - j[n] * y[n + 1];
                assert!(
                    close(w, want, 1e-9),
                    "Wronskian at n = {n}, z = {z:?}: {w:?} vs 2/(pi z) = {want:?}"
                );
            }
        }
    }

    /// `Y` satisfies the same three-term recurrence as `J` — which is
    /// how the higher orders are produced, so this checks the sweep
    /// rather than merely restating it.
    #[test]
    fn the_y_recurrence_is_satisfied() {
        for z in [C::new(1.5, 0.5), C::new(-2.0, 1.0), C::new(4.0, -3.0)] {
            let y = bessel_y_array_c(10, z).unwrap();
            let inv = z.inv();
            for n in 1..9 {
                let lhs = y[n - 1] + y[n + 1];
                let rhs = inv * y[n] * (2 * n) as f64;
                assert!(close(lhs, rhs, 1e-8), "Y recurrence n = {n}, z = {z:?}");
            }
        }
    }

    /// On the real axis `Y_n` must match the independently written
    /// vendored Cephes.
    #[test]
    fn y_matches_cephes_on_the_real_axis() {
        for &x in &[0.3, 1.0, 2.5, 5.0, 9.0] {
            let y = bessel_y_array_c(5, C::real(x)).unwrap();
            for (n, yv) in y.iter().enumerate() {
                let want = yn(n as isize, x);
                assert!(yv.im.abs() < 1e-11, "Y_{n}({x}) should be real, got {yv:?}");
                assert!(
                    rel_err(yv.re, want) < 1e-9,
                    "Y_{n}({x}) = {} vs cephes {want}",
                    yv.re
                );
            }
        }
    }

    /// **The absolute test for `K`**: `I_n K_{n+1} + I_{n+1} K_n = 1/z`
    /// (DLMF 10.28.2). Again elementary on the right, so it fixes `K`'s
    /// normalisation with no reference at all.
    #[test]
    fn the_i_k_wronskian_holds() {
        let zs = [
            C::real(0.6),
            C::real(3.0),
            C::new(1.0, 0.5),
            C::new(2.0, -1.5),
            C::new(0.8, 2.0),
        ];
        for z in zs {
            let want = z.inv();
            for n in 0..4i32 {
                let i_n = bessel_i_c(n, z).unwrap();
                let i_n1 = bessel_i_c(n + 1, z).unwrap();
                let k_n = bessel_k_c(n, z).unwrap();
                let k_n1 = bessel_k_c(n + 1, z).unwrap();
                let w = i_n * k_n1 + i_n1 * k_n;
                assert!(
                    close(w, want, 1e-8),
                    "I-K Wronskian at n = {n}, z = {z:?}: {w:?} vs 1/z = {want:?}"
                );
            }
        }
    }

    /// `K_n` on the real axis against the vendored Cephes `k0`/`k1`.
    #[test]
    fn k_matches_cephes_on_the_real_axis() {
        for &x in &[0.4, 1.0, 2.0, 4.0] {
            for (n, want) in [(0i32, k0(x)), (1i32, k1(x))] {
                let got = bessel_k_c(n, C::real(x)).unwrap();
                assert!(
                    got.im.abs() < 1e-9 * got.re.abs().max(1.0),
                    "K_{n}({x}) should be real, got {got:?}"
                );
                assert!(
                    rel_err(got.re, want) < 1e-8,
                    "K_{n}({x}) = {} vs cephes {want}",
                    got.re
                );
            }
        }
    }

    /// `Y` and `K` inherit the branch cut of `ln` along the negative
    /// real axis, so they are DISCONTINUOUS across it while `J` is not.
    /// Pinned so nobody assumes otherwise.
    #[test]
    fn y_is_discontinuous_across_the_negative_real_axis() {
        let eps = 1e-8;
        let above = C::new(-2.0, eps);
        let below = C::new(-2.0, -eps);
        // J is entire: the two sides agree
        let ja = bessel_j_c(0, above).unwrap();
        let jb = bessel_j_c(0, below).unwrap();
        assert!(close(ja, jb, 1e-6), "J should be continuous: {ja:?} vs {jb:?}");
        // Y is not
        let ya = bessel_y_c(0, above).unwrap();
        let yb = bessel_y_c(0, below).unwrap();
        let jump = (ya - yb).abs();
        assert!(
            jump > 0.1,
            "Y should jump across the cut, but moved only {jump}"
        );
        // The size of the jump follows from the series. Y carries a
        // (2/pi) ln(z/2) J term and nothing else discontinuous; crossing
        // the cut takes arg from +pi to -pi, a change of 2*pi, so
        //
        //     Y(above) - Y(below) = (2/pi) * (2 pi i) * J = 4i J.
        //
        // (A first draft of this test predicted 2i J and failed by
        // exactly a factor of two — the derivation above is the one the
        // measurement agrees with.)
        let predicted = ja * C::I * 4.0;
        assert!(
            close(ya - yb, predicted, 1e-5),
            "jump {:?} vs predicted 4i J_0 = {predicted:?}",
            ya - yb
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
        // Y and K are singular at the origin — an error, not an infinity
        assert!(bessel_y_c(0, C::ZERO).is_err(), "Y at z = 0");
        assert!(bessel_k_c(0, C::ZERO).is_err(), "K at z = 0");
        assert!(bessel_y_c(-1, C::ONE).is_err(), "negative order");
        assert!(bessel_k_c(-1, C::ONE).is_err(), "negative order");
    }
}
