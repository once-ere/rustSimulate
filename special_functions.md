# Special Functions — Reference and Examples

*The special-function layer of rustSimulate, cited to the NIST DLMF.
Every function below has at least one **intermediate** example (direct
evaluation you can check by hand) and one **expert** example (the
function doing real work inside a computation).*

This document grows one module at a time, alongside the code. Sections
marked **planned** are not implemented yet and say so rather than
pretending otherwise.

**Citation policy.** Functions cite DLMF equations by number and
permalink. The DLMF is copyright NIST and is *cited, never reproduced*;
formulas displayed here are from Abramowitz & Stegun (1964), a
public-domain US Government work, or written in our own notation. See
[THIRD_PARTY.md](THIRD_PARTY.md).

---

## Status at a glance

| module | source | state |
|---|---|---|
| `sph_bessel` — spherical Bessel jₙ, yₙ, derivatives | **native** | ✅ implemented, 7 tests + cross-validation |
| gamma, erf, Ei/Si/Ci, incomplete γ/B, Airy, cylindrical Bessel, elliptic, Jacobi, ζ | vendored Cephes | ✅ available, 11 identity tests |
| `legendre` — Pₙ, Pₗᵐ, spherical harmonics | native | **planned** |
| `orthopoly` — Hermite, Laguerre, Chebyshev, Gegenbauer, Jacobi | native | **planned** |
| `wigner` — 3j, 6j, Clebsch–Gordan | native | **planned** |

---

## 1. `sph_bessel` — spherical Bessel functions

Solutions of the radial Helmholtz equation in spherical coordinates.
They appear in every partial-wave expansion, in the free-particle
solution of the Schrödinger equation, and in Mie scattering.

Relation to the cylindrical functions — `DLMF 10.47.3`
(<https://dlmf.nist.gov/10.47.E3>):

```
    j_n(x) = sqrt(pi/2x) · J_{n+1/2}(x)
    y_n(x) = sqrt(pi/2x) · Y_{n+1/2}(x)
```

### API

| function | meaning | domain |
|---|---|---|
| `sph_j(n, x)` | jₙ(x), first kind | n ≥ 0, x finite (any sign) |
| `sph_y(n, x)` | yₙ(x), second kind | n ≥ 0, **x > 0** (singular at 0) |
| `sph_j_prime(n, x)` | jₙ′(x) | as `sph_j` |
| `sph_y_prime(n, x)` | yₙ′(x) | as `sph_y` |

All return `Result<f64, String>`; a domain violation is an `Err` with a
message naming the argument, never a silent `NaN`.

### Intermediate example — check against the closed forms

The two lowest orders have elementary closed forms (A&S 10.1.11):
j₀(x) = sin x / x and j₁(x) = sin x / x² − cos x / x.

```rust
use special_functions::sph_bessel::{sph_j, sph_y};

let x = 1.3_f64;
assert!((sph_j(0, x)? - x.sin() / x).abs() < 1e-15);
assert!((sph_j(1, x)? - (x.sin()/(x*x) - x.cos()/x)).abs() < 1e-14);

// y_0(x) = -cos x / x, and it diverges at the origin:
assert!((sph_y(0, x)? + x.cos() / x).abs() < 1e-15);
assert!(sph_y(0, 0.0).is_err());          // reported, not NaN
# Ok::<(), String>(())
```

Note the decay once the order exceeds the argument — this is the
behaviour that makes naive evaluation fail:

```rust
# use special_functions::sph_bessel::sph_j;
assert!(sph_j(20, 1.0)?.abs() < 1e-25);   // j_20(1) ~ 7.6e-26
# Ok::<(), String>(())
```

### Expert example — hard-sphere scattering phase shifts

For a hard sphere of radius *a*, the ℓ-th partial-wave phase shift of a
wave with momentum *k* is

```
    tan(delta_l) = j_l(ka) / y_l(ka)
```

and the total elastic cross-section is
σ = (4π/k²) Σ_ℓ (2ℓ+1) sin²(δ_ℓ). In the long-wavelength limit ka → 0
only the ℓ = 0 term survives and σ → 4πa², four times the geometric
cross-section — the classic result this example reproduces.

```rust
use special_functions::sph_bessel::{sph_j, sph_y};

fn hard_sphere_cross_section(k: f64, a: f64, l_max: i32) -> Result<f64, String> {
    let ka = k * a;
    let mut sigma = 0.0;
    for l in 0..=l_max {
        let delta = (sph_j(l, ka)? / sph_y(l, ka)?).atan();
        sigma += (2 * l + 1) as f64 * delta.sin().powi(2);
    }
    Ok(4.0 * std::f64::consts::PI / (k * k) * sigma)
}

// Low energy: the cross-section approaches 4 pi a^2.
let a = 1.0;
let sigma = hard_sphere_cross_section(1e-3, a, 6)?;
let geometric = std::f64::consts::PI * a * a;
assert!((sigma / geometric - 4.0).abs() < 1e-4);
# Ok::<(), String>(())
```

Two things worth noticing. The sum truncates safely at small `l_max`
because jₗ(ka) vanishes like (ka)^ℓ, so high partial waves contribute
nothing at low energy — physically, a slow particle cannot resolve the
sphere. And `sph_y` is the function that *diverges* at small argument,
which is exactly why the ratio jₗ/yₗ → 0 and the phase shifts vanish.

### Why the implementation looks the way it does

Both functions satisfy the same recurrence (`DLMF 10.51.1`,
<https://dlmf.nist.gov/10.51.E1>):

```
    f_{n+1}(x) = ((2n+1)/x) · f_n(x) - f_{n-1}(x)
```

but its numerical stability differs completely:

* **yₙ grows** with n, so recurring *upward* is stable — that is what
  `sph_y` does, seeded by the closed forms for y₀ and y₁.
* **jₙ decays** once n > x. Recurring upward there amplifies the seed's
  rounding error until the answer is noise. `sph_j` therefore uses
  **Miller's algorithm**: seed an artificial value far above the wanted
  order, recur *downward* (the direction in which the wanted solution
  dominates), rescaling to stay in range, then fix the overall scale
  against the exactly known j₀ = sin x / x.

A regression test, `sph_j_upward_recurrence_is_unstable_for_n_gt_x`,
runs the naive upward version alongside the real one and asserts the
naive result is wrong by more than six orders of magnitude at
n = 20, x = 1. It exists so nobody "simplifies" the implementation.

Small arguments (x < 10⁻⁴) use the leading series jₙ ≈ xⁿ/(2n+1)!!
(A&S 10.1.2) with a second-order correction, because the closed forms
suffer catastrophic cancellation there.

### How it is verified

| check | what it pins |
|---|---|
| closed forms, n ≤ 2 | absolute correctness at low order |
| **Wronskian** jₙyₙ′ − jₙ′yₙ = 1/x² (A&S 10.1.31) | 12 orders × 5 arguments — the strongest single identity |
| three-term recurrence | consistency across n for both functions |
| naive-upward blow-up | the stability property itself |
| small-x series limit | the x → 0 branch |
| parity jₙ(−x) = (−1)ⁿjₙ(x), origin values | sign and special-value handling |
| domain errors | `Err`, never silent `NaN` |
| **cross-validation** vs vendored J₍ₙ₊½₎ | two independent implementations agree |

That last one is the most valuable: our Miller-recurrence code and
Cephes's cylindrical Bessel of half-integer order were written
independently, and `DLMF 10.47.3` says they must agree. They do, to
1e-9 relative across n ∈ [0,10) and x ∈ [0.25, 40].

---

## 2. Vendored classical chapters (Cephes translation)

Available today through `special_functions::cephes`, covering DLMF
chapters 5–10, 19, 22 and 25 on **real arguments**. Provenance and the
licence caveat are in [THIRD_PARTY.md](THIRD_PARTY.md).

These are not our implementations, so we verify them rather than trust
them — `tests/vendored_identities.rs` checks:

| family | identity checked |
|---|---|
| Γ | Γ(½)=√π; Γ(n+1)=n! to n=15; reflection Γ(z)Γ(1−z)=π/sin πz; duplication |
| ψ | ψ(1)=−γ; ψ(z+1)=ψ(z)+1/z |
| erf | erf+erfc=1; odd symmetry; erf(6)→1 |
| Ai, Bi | Wronskian AiBi′−Ai′Bi = 1/π across x ∈ [−6,6]; exact values at 0 |
| Jᵥ, Yᵥ | Wronskian = 2/(πx); J₍±½₎ closed forms; recurrence |
| K, E | Legendre relation E(m)K(1−m)+E(1−m)K(m)−K(m)K(1−m) = π/2 |
| sn, cn, dn | sn²+cn²=1; m·sn²+dn²=1 |
| ζ | ζ(2)=π²/6, ζ(4)=π⁴/90, ζ(6)=π⁶/945 |

**Argument conventions**, pinned empirically before the tests were
written (getting these wrong yields tests that pass for the wrong
reason):

* `ellpe(m)` takes the parameter **m** directly — E(0)=π/2, E(1)=1.
* `ellpk(m₁)` takes the **complement** m₁ = 1 − m — so K(m) is
  `ellpk(1.0 - m)`.

---

*Sections for `legendre`, `orthopoly` and `wigner` will be added as
those modules land, in the same shape: API, intermediate example,
expert example, implementation rationale, verification table.*
