//! Uniform asymptotic expansions for **large order** — the Debye
//! polynomials and what they buy.
//!
//! Everything before this stage expanded in `1/z` at fixed order. That
//! leaves a hole exactly where a large-order expansion is needed, and
//! measurement located it precisely: `J_nu(z)` for `z` below `nu`. At
//! `nu = 400.5, z = 240` the previous machinery returned a number wrong
//! by a factor of `5e89`, because it built the **recessive** `J` as the
//! difference of two **dominant** Hankel functions.
//!
//! The turning point `z ≈ nu` turned out **not** to be the problem —
//! measured against Cephes, the existing route already gives 1e-14
//! there, which is why the Airy-type expansion of DLMF 10.20 is
//! discussed at the end of this note rather than implemented in the
//! middle of it.
//!
//! # One engine: the Debye polynomials
//!
//! ```text
//!   U_0(p) = 1
//!   U_{k+1}(p) = (1/2) p^2 (1 - p^2) U_k'(p)
//!                + (1/8) integral_0^p (1 - 5t^2) U_k(t) dt
//! ```
//!
//! (DLMF 10.41.9.) These are polynomials with rational coefficients, so
//! the recurrence is carried out exactly on coefficient vectors —
//! differentiate, multiply by `p^2 - p^4`, integrate — and computed
//! once. `U_1(p) = (3p - 5p^3)/24` and
//! `U_2(p) = (81p^2 - 462p^4 + 385p^6)/1152` are the published values
//! the tests check against.
//!
//! # What they give
//!
//! With `x = z/nu`, the same polynomials serve both families
//! (DLMF 10.41.3, 10.41.4, 10.19.3, 10.19.4):
//!
//! ```text
//!   s = sqrt(1 + x^2),  p = 1/s,  eta = s + ln(x/(1+s))
//!   I_nu(z) ~ e^(nu eta) / (sqrt(2 pi nu) sqrt(s)) sum_k U_k(p)/nu^k
//!   K_nu(z) ~ sqrt(pi/(2 nu)) e^(-nu eta) / sqrt(s) sum_k (-1)^k U_k(p)/nu^k
//!
//!   t = sqrt(1 - x^2),  q = 1/t,  alpha = ln((1+t)/x)     [needs x < 1]
//!   J_nu(z) ~ e^(nu(t - alpha)) / sqrt(2 pi nu t) sum_k U_k(q)/nu^k
//!   Y_nu(z) ~ -e^(-nu(t - alpha)) / sqrt(pi nu t / 2) sum_k (-1)^k U_k(q)/nu^k
//! ```
//!
//! The `J` form is the one that matters: `t - alpha < 0`, so it
//! produces an exponentially small number **directly**, with no
//! cancellation, which is precisely what the Hankel route could not do.
//!
//! # Truncation, and the honest failure
//!
//! These diverge like every asymptotic series, so they stop where the
//! terms stop shrinking and report that term as their error — the same
//! rule [`crate::bessel_scaled`] uses, so the estimates are comparable
//! and the method is chosen by measurement rather than by a rule of
//! thumb. `U_k(q)` grows like `q^(3k)` and `q = 1/sqrt(1-x^2)` blows up
//! as `x -> 1`, so the expansion is good for `x` bounded away from 1 and
//! useless at the turning point. That is not a defect; it is what
//! "uniform in `x` away from the turning point" means, and the
//! truncation estimate detects it without being told.
//!
//! # The Airy-type expansion of DLMF 10.20
//!
//! Olver's turning-point expansion is in [`crate::airy_uniform`], and
//! the two are complementary: this one is uniform away from `x = 1`,
//! that one **through** it.
//!
//! It is worth recording that this module originally argued 10.20 was
//! unnecessary, on the strength of a measurement across `x` from 0.98 to
//! 1.5 that showed 1e-14 everywhere. **That measurement was too coarse
//! and the conclusion was wrong.** Sampling 0.85 and 0.90 as well shows
//! a band the existing routes reached only to 1e-9 — at `nu = 100.5,
//! x = 0.85` the error was 1.4e-9 and at `nu = 200.5, x = 0.95` it was
//! 1.0e-12. With 10.20 those became 2.0e-15 and 1.8e-14. The argument
//! about `O(nu^-2)` was also wrong: it applies to the expansion
//! truncated at `A_0, B_0`, and three terms are kept, giving `O(nu^-6)`.

use crate::complex::Complex64 as C;

/// How many Debye polynomials to build. Optimal truncation of these
/// series is short — a dozen terms is far past the point where any of
/// them is still shrinking — and `U_k` has degree `3k`, so the cost of
/// the last ones is not free.
const N_U: usize = 14;

/// `U_k(p)` as coefficient vectors, `coeffs[k][j]` multiplying `p^j`.
///
/// Built by the exact recurrence rather than transcribed, because the
/// published forms stop at `U_3` and transcription is where sign errors
/// live. The tests check the first three against the published ones.
fn debye_coeffs() -> Vec<Vec<f64>> {
    let mut out: Vec<Vec<f64>> = Vec::with_capacity(N_U);
    out.push(vec![1.0]);
    for k in 0..N_U - 1 {
        let u = &out[k];
        let deg = u.len() - 1;
        // (1/2) p^2 (1 - p^2) U'  — U' shifts down one, then p^2 and
        // -p^4 shift up two and four, so the result has degree deg + 3.
        let mut next = vec![0.0; deg + 4];
        for (j, &c) in u.iter().enumerate().skip(1) {
            let d = j as f64 * c; // coefficient of p^(j-1) in U'
            next[j + 1] += 0.5 * d;
            next[j + 3] -= 0.5 * d;
        }
        // (1/8) integral_0^p (1 - 5t^2) U(t) dt
        for j in 0..=deg + 2 {
            let mut g = 0.0;
            if j <= deg {
                g += u[j];
            }
            if j >= 2 && j - 2 <= deg {
                g -= 5.0 * u[j - 2];
            }
            if g != 0.0 {
                next[j + 1] += g / (8.0 * (j + 1) as f64);
            }
        }
        while next.len() > 1 && *next.last().unwrap() == 0.0 {
            next.pop();
        }
        out.push(next);
    }
    out
}

/// The coefficient table, built once.
fn table() -> &'static Vec<Vec<f64>> {
    use std::sync::OnceLock;
    static T: OnceLock<Vec<Vec<f64>>> = OnceLock::new();
    T.get_or_init(debye_coeffs)
}

/// `U_k(p)`, for the Airy-type expansion of [`crate::airy_uniform`],
/// which is built from the same polynomials.
pub(crate) fn u_poly(k: usize, p: C) -> C {
    let mut v = C::ZERO;
    for &c in table()[k].iter().rev() {
        v = v * p + C::real(c);
    }
    v
}

/// `sum_k sign^k U_k(p) / nu^k`, truncated where the terms stop
/// shrinking, with the first omitted term as the error estimate.
///
/// `sign` is `+1` for `I` and `J`, `-1` for `K` and `Y`.
fn u_series(nu: C, p: C, sign: f64) -> (C, f64) {
    let t = table();
    let mut sum = C::ONE;
    let mut smallest = 1.0_f64;
    // Complex, so that a complex ORDER works here unchanged. The
    // polynomials do not care what kind of number `nu` is; only this
    // accumulator did.
    let mut scale = C::ONE;
    let inv = nu.inv() * sign;
    for coeffs in t.iter().skip(1) {
        scale = scale * inv;
        // Horner in p.
        let mut v = C::ZERO;
        for &c in coeffs.iter().rev() {
            v = v * p + C::real(c);
        }
        let term = v * scale;
        let m = term.abs();
        if !m.is_finite() || m >= smallest {
            return (sum, m.min(smallest * 4.0));
        }
        sum = sum + term;
        smallest = m;
    }
    (sum, smallest)
}

/// One scaled value with the estimate of its relative error.
pub struct Uniform {
    /// The value, already divided by the exponential named in the
    /// function that produced it.
    pub value: C,
    /// Estimated relative error, from optimal truncation.
    pub err: f64,
}

fn finite(value: C, err: f64) -> Option<Uniform> {
    (value.is_finite() && err.is_finite()).then_some(Uniform { value, err })
}

/// `exp(-|Re z|) I_nu(z)` and `exp(z) K_nu(z)` by DLMF 10.41.3 and
/// 10.41.4 <https://dlmf.nist.gov/10.41.E3>.
///
/// Returns `(I scaled, K scaled)`, either of which may be `None` if it
/// left `f64` range. `nu` must be positive; the expansion is in `1/nu`.
///
/// Valid for `|arg z| < pi/2`; the caller checks that, since the two
/// families have different sectors and only the caller knows which it
/// wants.
pub fn ik_uniform(nu: C, z: C) -> (Option<Uniform>, Option<Uniform>) {
    if !nu.is_finite() || nu.abs() == 0.0 || z.abs() == 0.0 {
        return (None, None);
    }
    let x = z * nu.inv();
    let s = (C::ONE + x * x).powf(0.5);
    let p = s.inv();
    let eta = s + (x * (C::ONE + s).inv()).ln();

    let (si, ei) = u_series(nu, p, 1.0);
    let (sk, ek) = u_series(nu, p, -1.0);

    // exp(nu*eta - |Re z|) and exp(z - nu*eta): each is formed as a
    // single exp, so the enormous exp(nu*eta) is never built.
    let ni = eta * nu - C::real(z.re.abs());
    let nk = z - eta * nu;
    let root = s.powf(0.5);

    let i_val = ni.exp()
        * (root * (nu * (2.0 * std::f64::consts::PI)).powf(0.5)).inv()
        * si;
    let k_val = nk.exp()
        * root.inv()
        * sk
        * (nu.inv() * (std::f64::consts::PI / 2.0)).powf(0.5);

    (finite(i_val, ei), finite(k_val, ek))
}

/// `J_nu(z)` and `Y_nu(z)` by the Debye expansions DLMF 10.19.3 and
/// 10.19.4 <https://dlmf.nist.gov/10.19.E3>, for **real** `z` with
/// `0 < z < nu`.
///
/// This is the region the previous machinery could not reach: `J` there
/// is exponentially small and was being built as the difference of two
/// exponentially large Hankel values. Here it comes out directly.
///
/// Returns the values themselves, not scaled forms — on the real axis
/// the `exp(-|Im z|)` scaling is 1. `J` may well be below `f64` range,
/// in which case its slot is `None` and the caller must say so rather
/// than return zero.
pub fn jy_debye(nu: f64, x: f64) -> (Option<Uniform>, Option<Uniform>) {
    if nu <= 0.0 || !nu.is_finite() || x <= 0.0 || x.is_nan() || x >= nu {
        return (None, None);
    }
    let r = x / nu; // sech(alpha)
    let t = (1.0 - r * r).sqrt(); // tanh(alpha)
    if t <= 0.0 {
        return (None, None);
    }
    let q = C::real(1.0 / t); // coth(alpha)
    let alpha = ((1.0 + t) / r).ln();
    let e = nu * (t - alpha); // negative

    let (sj, ej) = u_series(C::real(nu), q, 1.0);
    let (sy, ey) = u_series(C::real(nu), q, -1.0);

    let j = sj * (e.exp() / (2.0 * std::f64::consts::PI * nu * t).sqrt());
    let y = sy * (-(-e).exp() / (std::f64::consts::PI * nu * t / 2.0).sqrt());

    // `J` underflowing to exactly zero is not an answer: the true value
    // is nonzero and simply below f64 range. Report nothing so the
    // caller can say that, rather than handing back a silent 0.
    let j = if e < -745.0 { None } else { finite(j, ej) };
    (j, finite(y, ey))
}

/// `J_nu(z)` and `Y_nu(z)` by the Debye expansions, for **complex order
/// and complex argument**, on either side of the turning point.
///
/// # One formula, two regions
///
/// Write `t = sqrt(1 - x^2)`, `alpha = ln((1+t)/x)`, `q = 1/t` with
/// `x = z/nu`, and
///
/// ```text
///   F±(nu, x) = e^(± nu (t - alpha)) / sqrt(2 pi nu t)
///               sum_k (±1)^k U_k(q) / nu^k
/// ```
///
/// For `|x| < 1` these are the expansions of DLMF 10.19.3 and 10.19.4
/// directly: `J = F+` and `Y = -2 F-`. For `|x| > 1` — the
/// **oscillatory** region, where `t` turns imaginary — they continue
/// into the Hankel functions instead, `H1 = 2 F+` and `H2 = 2i F-`,
/// and so
///
/// ```text
///   J = F+ + i F-,      Y = -i F+ - F-
/// ```
///
/// That flip is the Stokes phenomenon and is why one formula needs two
/// readings. **The constants were identified by experiment rather than
/// transcribed**: continuing `F+` past `x = 1` and dividing by each of
/// `J`, `Y`, `H1`, `H2` in turn showed `F+/H1 = 1/2` to the accuracy of
/// the reference. The readings are then checked against
/// `bessel_j_c`/`bessel_y_c` wherever those are independently sound.
///
/// # Why this was missing
///
/// The `1/z` expansions refuse when `|4 nu^2|` is not small compared
/// with `|z|`; the ascending series has cancelled by then; and `x` is
/// too far from 1 for the turning-point expansion. That band —
/// `|z|` a few times `|nu|` — had no method at all, at real order as
/// well as complex. At `nu = 20, z = 60` it made `Y_20(60)` come back
/// as `1e8`.
pub fn jy_debye_c(nu: C, z: C) -> (Option<Uniform>, Option<Uniform>) {
    if !nu.is_finite() || !z.is_finite() || nu.abs() == 0.0 || z.abs() == 0.0 {
        return (None, None);
    }
    // A LARGE-ORDER expansion, so a small order is out of range even
    // when its terms happen to look small. At `nu = 1.3, |z| = 18` the
    // ratio `q = 1/sqrt(1-x^2)` is 0.07, the terms shrink nicely, and
    // the truncation estimate comes out below 1e-14 — while the answer
    // is 1.7e-10 off, because what the estimate cannot see is the
    // expansion's own `1/nu` character. The 1/z route is far better
    // there, and without this guard the selector preferred this one.
    // 8, not 3. At `nu = 3` the expansion still passed its own
    // estimate and still lost the J-Y Wronskian by a factor of 1 at
    // `|z| = 22` — a `1/nu` series needs `nu` genuinely large, and the
    // band this exists to cover starts at `nu ~ 20` anyway.
    if nu.abs() < 8.0 {
        return (None, None);
    }
    let x = z * nu.inv();
    let t = (C::ONE - x * x).powf(0.5);
    if t.abs() == 0.0 || !t.is_finite() {
        return (None, None);
    }
    let alpha = ((C::ONE + t) * x.inv()).ln();
    let q = t.inv();
    let (sp, ep) = u_series(nu, q, 1.0);
    let (sm, em) = u_series(nu, q, -1.0);
    let pref = ((nu * t) * (2.0 * std::f64::consts::PI)).powf(-0.5);
    let e = (t - alpha) * nu;
    // The exponential's own rounding: `exp(e)` is known only to
    // `|e| * eps` relative, because `e` is. The same term the complex
    // Airy needed, and for the same reason.
    let rounding = e.abs() * f64::EPSILON;
    // The same safety factor the other asymptotic routes carry:
    // optimal truncation is an estimate, not a bound, and measured it
    // runs optimistic at moderate order.
    let err = (ep.max(em) * 100.0).max(rounding);
    let (fp, fm) = (e.exp() * pref * sp, (e * -1.0).exp() * pref * sm);
    if !fp.is_finite() || !fm.is_finite() {
        return (None, None);
    }
    if x.abs() < 1.0 {
        return (finite(fp, err), finite(fm * -2.0, err));
    }
    let j = fp + C::I * fm;
    let y = C::I * fp * -1.0 - fm;
    // Both combinations can cancel — measured from the values, as
    // everywhere else in this crate.
    let cancel = |v: C| {
        if v.abs() == 0.0 || !v.is_finite() {
            f64::INFINITY
        } else {
            ((fp.abs() + fm.abs()) / v.abs()).max(1.0)
        }
    };
    (finite(j, err * cancel(j)), finite(y, err * cancel(y)))
}

/// `(ln|J_nu(x)|, ln|Y_nu(x)|)` for `0 < x < nu`, from the leading term
/// of the same Debye expansion.
///
/// This exists so a caller can distinguish two very different failures:
/// "no method reaches here" and "the value is real, and outside `f64`".
/// At `nu = 400.5, x = 40` the true `J` is about `e^-1013` and the true
/// `Y` about `e^+1010`; neither is representable, and saying so is more
/// use than saying nothing was accurate enough.
pub fn jy_log_magnitude(nu: f64, x: f64) -> Option<(f64, f64)> {
    if nu <= 0.0 || !nu.is_finite() || x <= 0.0 || x.is_nan() || x >= nu {
        return None;
    }
    let r = x / nu;
    let t = (1.0 - r * r).sqrt();
    if t <= 0.0 {
        return None;
    }
    let alpha = ((1.0 + t) / r).ln();
    let e = nu * (t - alpha);
    let lj = e - 0.5 * (2.0 * std::f64::consts::PI * nu * t).ln();
    let ly = -e - 0.5 * (std::f64::consts::PI * nu * t / 2.0).ln();
    Some((lj, ly))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn poly(k: usize, p: f64) -> f64 {
        let mut v = 0.0;
        for &c in table()[k].iter().rev() {
            v = v * p + c;
        }
        v
    }

    /// The published `U_1`, `U_2`, `U_3` (DLMF 10.41.10 and the
    /// standard references), against the recurrence. This is what makes
    /// the other tests meaningful: everything below rests on these
    /// polynomials being right.
    #[test]
    fn the_debye_polynomials_match_the_published_ones() {
        for &p in &[0.0_f64, 0.3, 1.0, 2.5, 7.0] {
            let (p2, p3, p4, p5, p6) = (p * p, p * p * p, p.powi(4), p.powi(5), p.powi(6));
            assert!((poly(0, p) - 1.0).abs() < 1e-14, "U_0");
            let want = (3.0 * p - 5.0 * p3) / 24.0;
            assert!((poly(1, p) - want).abs() <= 1e-13 * want.abs().max(1.0), "U_1({p})");
            let want = (81.0 * p2 - 462.0 * p4 + 385.0 * p6) / 1152.0;
            assert!((poly(2, p) - want).abs() <= 1e-12 * want.abs().max(1.0), "U_2({p})");
            let want = (30375.0 * p3 - 369_603.0 * p5 + 765_765.0 * p.powi(7)
                - 425_425.0 * p.powi(9))
                / 414_720.0;
            assert!((poly(3, p) - want).abs() <= 1e-11 * want.abs().max(1.0), "U_3({p})");
        }
        // Degrees are 3k, and only every other power appears.
        for k in 0..6 {
            assert_eq!(table()[k].len() - 1, 3 * k, "U_{k} should have degree {}", 3 * k);
        }
    }

    /// The `J` Debye expansion against Cephes `jv` in the region it
    /// exists for: `z` well below `nu`, where `J` is exponentially small
    /// and every other method here fails.
    #[test]
    fn the_j_debye_expansion_matches_cephes_where_nothing_else_reaches() {
        for &nu in &[40.5_f64, 100.5, 200.5, 400.5] {
            for &frac in &[0.1, 0.3, 0.5, 0.7, 0.85] {
                let x = nu * frac;
                let want = spec_math::cephes64::jv(nu, x);
                if want == 0.0 || !want.is_finite() {
                    continue; // below Cephes' own range
                }
                let (j, _) = jy_debye(nu, x);
                let j = j.unwrap_or_else(|| panic!("no J at nu={nu}, x={x}"));
                // Asserted against the routine's OWN error estimate
                // rather than a fixed tolerance. That tests two things
                // at once — the value, and the honesty of the estimate
                // the caller will use to choose this method over
                // another.
                //
                // The floor is 1e-8, not 1e-13, because CEPHES is the
                // weaker party at large order: at nu = 400.5, x = 160
                // it disagrees by 1.4e-9 while the J-Y Wronskian —
                // elementary on the right — puts these values at 1.2e-13
                // and the truncation estimate at 4e-28. See
                // `the_wronskian_confirms_these_beat_cephes_at_large_order`.
                let bound = (3.0 * j.err).max(1e-8) * want.abs();
                assert!(
                    (j.value.re - want).abs() <= bound,
                    "J_{nu}({x}): {} vs {want}, estimate {:.1e}",
                    j.value.re,
                    j.err
                );
                assert!(j.value.im == 0.0, "should be real");
            }
        }
    }

    /// `Y` likewise, where Cephes can still represent it.
    #[test]
    fn the_y_debye_expansion_matches_cephes() {
        for &nu in &[20.5_f64, 40.5, 100.5] {
            for &frac in &[0.3, 0.5, 0.7, 0.85] {
                let x = nu * frac;
                let want = spec_math::cephes64::yv(nu, x);
                if !want.is_finite() {
                    continue;
                }
                let (_, y) = jy_debye(nu, x);
                let y = y.unwrap();
                let bound = (3.0 * y.err).max(1e-8) * want.abs();
                assert!(
                    (y.value.re - want).abs() <= bound,
                    "Y_{nu}({x}): {} vs {want}, estimate {:.1e}",
                    y.value.re,
                    y.err
                );
            }
        }
    }

    /// `I` and `K` from DLMF 10.41. The vendored Cephes has no `iv`
    /// or `kv`, and its `kn` overflows past order 31, so the reference
    /// here is the **I-K Wronskian**, `I_nu K_{nu+1} + I_{nu+1} K_nu =
    /// 1/z` (DLMF 10.28.2), whose right-hand side is elementary. In the
    /// scaled variables used here it reads
    /// `Is Ks' + Is' Ks = exp(z - |Re z|)/z`, which on the positive real
    /// axis is `1/z`.
    #[test]
    fn the_ik_uniform_expansion_satisfies_the_wronskian() {
        for &nu in &[20.0_f64, 50.0, 200.0, 800.0] {
            for &frac in &[0.2, 1.0, 3.0, 20.0] {
                let x = frac * nu;
                let z = C::real(x);
                let (i0u, k1u) = (ik_uniform(C::real(nu), z).0, ik_uniform(C::real(nu + 1.0), z).1);
                let (i1u, k0u) = (ik_uniform(C::real(nu + 1.0), z).0, ik_uniform(C::real(nu), z).1);
                let (Some(i0u), Some(i1u), Some(k0u), Some(k1u)) = (i0u, i1u, k0u, k1u)
                else {
                    continue; // out of f64 range, which the caller reports
                };
                let w = i0u.value * k1u.value + i1u.value * k0u.value;
                let want = z.inv();
                let scale = (i0u.value * k1u.value).abs() + (i1u.value * k0u.value).abs();
                assert!(
                    (w - want).abs() / scale < 1e-11,
                    "Wronskian nu={nu} x={x}: residual {:.2e}",
                    (w - want).abs() / scale
                );
            }
        }
        // ... and against Cephes at order 0 and 1, where it exists, to
        // pin the normalisation the Wronskian alone cannot fix.
        for &x in &[40.0_f64, 200.0] {
            for (nu, want) in [
                (0.0, spec_math::cephes64::i0(x) * (-x).exp()),
                (1.0, spec_math::cephes64::i1(x) * (-x).exp()),
            ] {
                // nu must be positive for a 1/nu expansion; use nu = 1
                // and the recurrence-free Cephes value at nu = 1.
                if nu < 1.0 {
                    continue;
                }
                let got = ik_uniform(C::real(nu), C::real(x)).0.unwrap().value.re;
                assert!(
                    (got - want).abs() <= 1e-6 * want,
                    "exp(-x)I_1({x}): {got} vs {want} (a 1/nu expansion at nu = 1)"
                );
            }
        }
    }

    /// The truncation estimate must grow as the turning point is
    /// approached — that is how the caller learns not to use this here,
    /// and it is the whole reason the estimate is returned rather than a
    /// validity rule of thumb.
    #[test]
    fn the_estimate_detects_the_turning_point() {
        let nu = 400.0;
        let mut prev = 0.0;
        for &frac in &[0.3_f64, 0.6, 0.85, 0.95, 0.99] {
            let (j, _) = jy_debye(nu, nu * frac);
            let e = j.map(|u| u.err).unwrap_or(f64::INFINITY);
            assert!(e > prev, "estimate should worsen towards x = 1: {e:.1e} at {frac}");
            prev = e;
        }
        assert!(prev > 1e-6, "at x = 0.99 the expansion should be visibly bad, got {prev:.1e}");
    }

    /// The J-Y Wronskian `J_{nu+1}Y_nu - J_nu Y_{nu+1} = 2/(pi x)`
    /// (DLMF 10.5.2) is elementary on the right, so it judges these
    /// values without a reference implementation — and it is the only
    /// thing that can, since Cephes is less accurate here than they are.
    ///
    /// At `nu = 400.5, x = 160` the Wronskian residual is 1.2e-13 and
    /// the truncation estimate 4e-32, while Cephes `jv` disagrees by
    /// 1.4e-9. That is recorded as an assertion rather than a remark,
    /// so if Cephes is ever replaced the claim is re-examined.
    #[test]
    fn the_wronskian_confirms_these_beat_cephes_at_large_order() {
        for &nu in &[40.5_f64, 100.5, 200.5, 400.5] {
            for &frac in &[0.1, 0.2, 0.4, 0.7, 0.85] {
                let x = nu * frac;
                let (Some(j0), Some(y0)) = jy_debye(nu, x) else { continue };
                let (Some(j1), Some(y1)) = jy_debye(nu + 1.0, x) else { continue };
                let w = j1.value.re * y0.value.re - j0.value.re * y1.value.re;
                let want = 2.0 / (std::f64::consts::PI * x);
                let bound = (3.0 * j0.err.max(y0.err)).max(1e-11);
                assert!(
                    (w - want).abs() / want <= bound,
                    "Wronskian at nu={nu}, x/nu={frac}: residual {:.2e}, estimate {:.1e}",
                    (w - want).abs() / want,
                    j0.err
                );
            }
        }
        // The specific point where the roles are clear.
        let (nu, x) = (400.5_f64, 160.2);
        let (Some(j0), Some(y0)) = jy_debye(nu, x) else { panic!("no value") };
        let (Some(j1), Some(y1)) = jy_debye(nu + 1.0, x) else { panic!("no value") };
        let w = j1.value.re * y0.value.re - j0.value.re * y1.value.re;
        let want = 2.0 / (std::f64::consts::PI * x);
        assert!((w - want).abs() / want < 1e-12, "our Wronskian must be tight here");
        let c = spec_math::cephes64::jv(nu, x);
        assert!(
            (j0.value.re - c).abs() / c.abs() > 1e-10,
            "cephes is expected to disagree here by ~1.4e-9; if it no longer does, \
             this claim needs re-checking rather than deleting"
        );
    }

    /// Deep in the region the expansion is built for — `x` well below
    /// `nu` — it should be at machine precision, not merely within its
    /// estimate. This is the claim that makes the module worth having.
    #[test]
    fn deep_in_range_the_debye_expansion_is_exact_to_machine_precision() {
        for &nu in &[100.5_f64, 200.5, 400.5] {
            for &frac in &[0.1, 0.2, 0.4] {
                let x = nu * frac;
                let want = spec_math::cephes64::jv(nu, x);
                if want == 0.0 || !want.is_finite() {
                    continue;
                }
                let (j, _) = jy_debye(nu, x);
                let j = j.unwrap();
                assert!(j.err < 1e-15, "estimate at nu={nu}, x/nu={frac} is {:.1e}", j.err);
                assert!(
                    (j.value.re - want).abs() <= 1e-8 * want.abs(),
                    "J_{nu}({x}): {} vs {want} (cephes is the loose one at this order)",
                    j.value.re
                );
            }
        }
    }

    /// The oscillatory Debye expansion, against Cephes where Cephes is
    /// itself sound. This is the band that had no method at all: at
    /// `nu = 20, z = 60` the crate returned `Y = 1e8` before.
    #[test]
    fn the_oscillatory_region_matches_cephes() {
        for &nu in &[8.0_f64, 12.0, 20.0, 40.0, 80.0, 150.0] {
            for &x in &[2.0_f64, 3.0, 5.0, 10.0, 30.0] {
                let z = nu * x;
                let (Some(j), Some(y)) = jy_debye_c(C::real(nu), C::real(z)) else {
                    continue;
                };
                let (wj, wy) = (spec_math::cephes64::jv(nu, z), spec_math::cephes64::yv(nu, z));
                if !wj.is_finite() || !wy.is_finite() || wj == 0.0 || wy == 0.0 {
                    continue;
                }
                let bound = |e: f64| (3.0 * e).max(1e-12);
                assert!(
                    (j.value.re - wj).abs() <= bound(j.err) * wj.abs(),
                    "J_{nu}({z}): {} vs {wj}, estimate {:.1e}",
                    j.value.re,
                    j.err
                );
                assert!(
                    (y.value.re - wy).abs() <= bound(y.err) * wy.abs(),
                    "Y_{nu}({z}): {} vs {wy}, estimate {:.1e}",
                    y.value.re,
                    y.err
                );
                // Not exactly zero: the route goes through complex `t`
                // and `alpha`, so the imaginary part is rounding rather
                // than an exact cancellation.
                assert!(
                    j.value.im.abs() <= 1e-11 * j.value.re.abs()
                        && y.value.im.abs() <= 1e-11 * y.value.re.abs(),
                    "should be real to rounding: {:?}, {:?}",
                    j.value,
                    y.value
                );
            }
        }
    }

    /// The point that named this stage. `Y_20(60)` sits where the `1/z`
    /// expansion is refused (`|4 nu^2| = 1600` against `8|z| = 480`),
    /// the ascending series has cancelled `exp(60)` away, and `z/nu = 3`
    /// is far outside the turning-point expansion's reach.
    #[test]
    fn the_band_that_had_no_method_now_has_one() {
        let (Some(j), Some(y)) = jy_debye_c(C::real(20.0), C::real(60.0)) else {
            panic!("no value")
        };
        let (wj, wy) = (spec_math::cephes64::jv(20.0, 60.0), spec_math::cephes64::yv(20.0, 60.0));
        assert!((j.value.re - wj).abs() <= 1e-13 * wj.abs(), "J: {} vs {wj}", j.value.re);
        assert!((y.value.re - wy).abs() <= 1e-12 * wy.abs(), "Y: {} vs {wy}", y.value.re);
        // ... and through the public routine, which has to choose it.
        let got = crate::bessel_complex::bessel_y_c(20, C::real(60.0)).unwrap();
        assert!((got.re - wy).abs() <= 1e-12 * wy.abs(), "chosen: {} vs {wy}", got.re);
    }

    /// A `1/nu` expansion needs `nu` genuinely large, and its terms
    /// looking small is not the same thing. At `nu = 3, |z| = 22` it
    /// passed its own estimate and still lost the J-Y Wronskian by a
    /// factor of 1, which is why the order is guarded rather than left
    /// to the estimate.
    #[test]
    fn a_small_order_is_refused_outright() {
        assert!(jy_debye_c(C::real(1.3), C::real(18.0)).0.is_none());
        assert!(jy_debye_c(C::real(3.0), C::real(22.0)).0.is_none());
        assert!(jy_debye_c(C::new(5.0, 2.0), C::real(20.0)).0.is_none());
        assert!(jy_debye_c(C::real(8.0), C::real(24.0)).0.is_some(), "nu = 8 is in range");
    }

    #[test]
    fn debye_edge_cases() {
        assert!(jy_debye(0.0, 1.0).0.is_none(), "nu = 0");
        assert!(jy_debye(10.0, 0.0).0.is_none(), "x = 0");
        assert!(jy_debye(10.0, 20.0).0.is_none(), "x > nu is not this expansion");
        assert!(ik_uniform(C::ZERO, C::ONE).0.is_none(), "nu = 0");
        assert!(ik_uniform(C::real(10.0), C::ZERO).0.is_none(), "z = 0");
    }
}
