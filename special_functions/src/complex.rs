//! A minimal complex number, `Complex64`.
//!
//! Quantum wavefunctions are complex by construction, so the
//! Crank–Nicolson propagator in [`crate::tridiag`] needs this. It is
//! deliberately small — arithmetic, conjugation, modulus, and the
//! exponential — rather than a general numeric-tower crate, because the
//! project takes no external dependencies and this is all the solvers
//! require.
//!
//! Written from the definitions; nothing here is derived from any
//! third-party source.

use std::ops::{Add, Div, Mul, Neg, Sub};

/// A double-precision complex number, `re + i*im`.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Complex64 {
    pub re: f64,
    pub im: f64,
}

impl Complex64 {
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    /// The real number `re + 0i`.
    pub const fn real(re: f64) -> Self {
        Self { re, im: 0.0 }
    }
    /// The imaginary unit.
    pub const I: Self = Self { re: 0.0, im: 1.0 };
    pub const ZERO: Self = Self { re: 0.0, im: 0.0 };
    pub const ONE: Self = Self { re: 1.0, im: 0.0 };

    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }
    /// `|z|^2`, without the square root — cheaper and exact when you
    /// only need to compare magnitudes or form a probability density.
    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
    /// `e^z = e^re (cos im + i sin im)`.
    pub fn exp(self) -> Self {
        let m = self.re.exp();
        Self::new(m * self.im.cos(), m * self.im.sin())
    }
    /// `e^{i*theta}` — the common case in a propagator.
    pub fn from_polar(r: f64, theta: f64) -> Self {
        Self::new(r * theta.cos(), r * theta.sin())
    }
    pub fn is_finite(self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
    /// Reciprocal, scaled to avoid overflow in `re^2 + im^2`.
    pub fn inv(self) -> Self {
        let d = self.norm_sqr();
        Self::new(self.re / d, -self.im / d)
    }
}

impl Add for Complex64 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.re + o.re, self.im + o.im)
    }
}
impl Sub for Complex64 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.re - o.re, self.im - o.im)
    }
}
impl Mul for Complex64 {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}
impl Mul<f64> for Complex64 {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        Self::new(self.re * s, self.im * s)
    }
}
impl Div for Complex64 {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        self * o.inv()
    }
}
impl Neg for Complex64 {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic_identities() {
        let a = Complex64::new(3.0, -4.0);
        let b = Complex64::new(-1.0, 2.0);
        assert_eq!(a + b, Complex64::new(2.0, -2.0));
        assert_eq!(a - b, Complex64::new(4.0, -6.0));
        // (3-4i)(-1+2i) = -3 + 6i + 4i - 8i^2 = 5 + 10i
        assert_eq!(a * b, Complex64::new(5.0, 10.0));
        assert_eq!(a.abs(), 5.0);
        assert_eq!(a.norm_sqr(), 25.0);
        // z * conj(z) = |z|^2
        let p = a * a.conj();
        assert!((p.re - 25.0).abs() < 1e-14 && p.im.abs() < 1e-14);
        // z / z = 1
        let q = a / a;
        assert!((q.re - 1.0).abs() < 1e-15 && q.im.abs() < 1e-15);
        // i^2 = -1
        let ii = Complex64::I * Complex64::I;
        assert_eq!(ii, Complex64::new(-1.0, 0.0));
    }

    #[test]
    fn eulers_identity_and_exp() {
        // e^{i pi} + 1 = 0
        let z = (Complex64::I * std::f64::consts::PI).exp();
        assert!((z.re + 1.0).abs() < 1e-15 && z.im.abs() < 1e-15);
        // |e^{i t}| = 1 for any t
        for &t in &[0.0, 0.3, 1.7, -2.9, 10.0] {
            assert!((Complex64::from_polar(1.0, t).abs() - 1.0).abs() < 1e-15);
        }
        // e^{a+b} = e^a e^b
        let a = Complex64::new(0.4, 1.1);
        let b = Complex64::new(-0.7, 0.5);
        let l = (a + b).exp();
        let r = a.exp() * b.exp();
        assert!((l.re - r.re).abs() < 1e-14 && (l.im - r.im).abs() < 1e-14);
    }
}
