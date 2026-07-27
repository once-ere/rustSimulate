# Special functions — provenance, coverage, and test results

**Date:** 2026-07-26
**Scope:** the `special_functions` crate and its exposure through the
posim language.

This is the report the DLMF task asked for. It leads with the coverage
number because that is the fact most likely to be overstated, and
everything else here is only worth as much as that number is honest.

---

## 1. Coverage, stated first

> **13 of the DLMF's 33 function chapters have some coverage. None of
> the 33 is complete.**

That is roughly 40 % of chapters *touched* and considerably less than
40 % of DLMF content implemented. Twenty chapters have nothing at all.

The original task was "port the complete NIST DLMF". That is not what
this is, and the reason is structural rather than an excuse:

- **No pure-Rust DLMF port exists.** `DLMF language:Rust` on GitHub
  returns zero repositories; a crates.io full-text search for "dlmf"
  returns two unrelated crates. The alternative branch of the task —
  import rather than port — has nothing to import.
- **No complete DLMF implementation exists in any language.**
  Mathematica reaches roughly 90 % of the function chapters (Heun only
  since 2020), Maple ~85 %, mpmath ~70 %, Arb/FLINT ~65 %, SciPy ~55 %,
  GSL ~48 %, Boost.Math ~42 %.
- **The DLMF is a reference handbook of properties, not an algorithm
  specification.** Most chapters state identities and asymptotics with
  no numerically stable evaluation recipe. "Porting" a chapter is
  numerical-analysis research, not transcription. Chapters 29 (Lamé),
  32 (Painlevé), 35 (matrix argument) and 36 (coalescing saddles) have
  no solid general implementation in *any* system.

So this milestone implements what a physics simulator actually consumes,
and the matrix below states plainly what is absent.

### The 33-chapter matrix

| ch. | title | status | what exists here |
|---|---|---|---|
| 4 | Elementary Functions | partial | `sqrt abs sin cos exp log` in the language; Rust `std` underneath |
| 5 | Gamma Function | partial | **vendored** Cephes: gamma, lgamma, rgamma, digamma, poch |
| 6 | Exponential, Log, Sine and Cosine Integrals | partial | **vendored**: Ei, Eₙ, Si, Ci, Shi, Chi |
| 7 | Error Functions, Dawson, Fresnel | partial | **vendored**: erf, erfc + inverses, Dawson, Fresnel |
| 8 | Incomplete Gamma and Related | partial | **vendored**: incomplete gamma/beta + inverses |
| 9 | Airy and Related | partial | **vendored**: Ai, Bi and derivatives, real argument only |
| 10 | Bessel Functions | partial | **vendored** cylindrical Jᵥ Yᵥ Iᵥ Kᵥ (real arg); **native** spherical jₙ yₙ and derivatives; **native** integer-order Jₙ whole-table |
| 11 | Struve and Related | **none** | — |
| 12 | Parabolic Cylinder | **none** | — |
| 13 | Confluent Hypergeometric | **none** | — |
| 14 | Legendre and Related | partial | **native**: Pₙ, P′ₙ, Pₗᵐ, normalised P̄ₗᵐ, Yₗᵐ complex and real |
| 15 | Hypergeometric Function | **none** | — |
| 16 | Generalized Hypergeometric, Meijer G | **none** | — |
| 17 | q-Hypergeometric | **none** | — |
| 18 | Orthogonal Polynomials | partial | **native**: Hermite Hₙ/Heₙ, Laguerre Lₙ/Lₙ^α, Chebyshev Tₙ/Uₙ, Gegenbauer Cₙ^α, Jacobi Pₙ^{α,β}. The discrete Askey-scheme families are absent |
| 19 | Elliptic Integrals | partial | **vendored**: Legendre and Carlson forms |
| 20 | Theta Functions | **none** | — |
| 21 | Multidimensional Theta | **none** | — |
| 22 | Jacobian Elliptic Functions | partial | **vendored**: sn, cn, dn |
| 23 | Weierstrass Elliptic and Modular | **none** | — |
| 24 | Bernoulli and Euler Polynomials | **none** | — |
| 25 | Zeta and Related | partial | **vendored**: ζ, Hurwitz ζ, dilogarithm |
| 26 | Combinatorial Analysis | **none** | — |
| 27 | Functions of Number Theory | **none** | — |
| 28 | Mathieu Functions and Hill's Equation | **none** | — |
| 29 | Lamé Functions | **none** | no general implementation exists anywhere |
| 30 | Spheroidal Wave Functions | **none** | — |
| 31 | Heun Functions | **none** | — |
| 32 | Painlevé Transcendents | **none** | no general implementation exists anywhere |
| 33 | Coulomb Functions | **none** | — |
| 34 | 3j, 6j, 9j Symbols | partial | **native**: 3-j, 6-j, **9-j**, Clebsch–Gordan. The chapter's asymptotics and generating functions are absent |
| 35 | Functions of Matrix Argument | **none** | no general implementation exists anywhere |
| 36 | Integrals with Coalescing Saddles | **none** | no general implementation exists anywhere |

Beyond the DLMF chapters, the crate also carries the numerical
infrastructure a simulator needs and the DLMF does not catalogue:
`eigen` (cyclic Jacobi, dense), `lanczos` (matrix-free symmetric
eigensolver for problems too large to form), `quadrature`
(Gauss–Legendre, adaptive Simpson, Brent roots), `tridiag` (Thomas,
Sherman–Morrison) and `complex`.

---

## 2. Licensing and provenance

### DLMF itself

The DLMF is **not** public domain — its authors assigned copyright to
NIST, which permits limited copying for research and teaching and
forbids commercial reproduction. NIST publishes no reference
implementation.

The policy followed here: **cite DLMF equation numbers and implement the
mathematics independently.** No DLMF prose is reproduced and no equation
displays are bulk-transcribed. Where a formula needed to be *visible* in
documentation it is taken from Abramowitz & Stegun (1964), a genuine US
Government work in the public domain.

### Vendored code

`vendor/spec_math/` — from <https://github.com/matthew-romanowicz/spec_math>,
a Rust translation of Cephes (Stephen L. Moshier). 102 files, **101
byte-identical to upstream**; only `src/lib.rs` differs, carrying a
provenance header plus `#![forbid(unsafe_code)]` and
`#![allow(dead_code)]`. Zero runtime dependencies, no `unsafe`.

**Licensing caveat, recorded rather than papered over:** the upstream
repository ships *no LICENSE file*. `MIT OR Apache-2.0` is declared only
in its `Cargo.toml`. That declaration is reproduced in `THIRD_PARTY.md`
along with the Cephes lineage, and it is flagged for resolution upstream
if certainty is wanted. See [THIRD_PARTY.md](THIRD_PARTY.md).

**Rejected, with reasons:** `puruspe` — algorithms and identifiers derive
from *Numerical Recipes*, whose licence forbids redistribution.
`scilib` — GPL-3.0, would infect the repository, and dead since 2023.
`libm` — 21 `unsafe` sites.

### Clean-room replacements

`bessel`, `tridiag` and `complex` replace licence-encumbered routines
from the SolveIt C++ sources (Numerical Recipes and GSL). Their audit
record is separate and detailed:
[CLEANROOM_PROVENANCE.md](CLEANROOM_PROVENANCE.md).

---

## 3. Per-module provenance

| module | origin | algorithm | citation |
|---|---|---|---|
| `sph_bessel` | native | Miller downward recurrence below `n > x`, upward above | DLMF 10.51, A&S 10.1 |
| `legendre` | native | ascending recurrence in ℓ seeded by Pₘᵐ; normalised form computed *directly* in the normalised basis so it never overflows | DLMF 14.10, A&S 8.5 |
| `orthopoly` | native | three-term recurrences; Clenshaw for series | DLMF 18.9, A&S 22.7 |
| `wigner` | native | Racah single-sum for 3-j and 6-j, all factorials in logarithms; 9-j as a single sum over 6-j | DLMF 34.2.4, 34.4.1, 34.6.1; Edmonds 1957 §3.6 |
| `bessel` | native, **clean-room** | Miller downward recurrence, scale fixed by `J₀+2(J₂+J₄+…)=1` | DLMF 10.6.1, 10.12.4; A&S 9.1.27, 9.1.46 |
| `tridiag` | native, **clean-room** | Thomas algorithm; Sherman–Morrison for the cyclic case | textbook |
| `eigen` | native | cyclic Jacobi, ≤100 sweeps | textbook |
| `lanczos` | native | Lanczos with full reorthogonalisation + deflation, matrix-free | textbook |
| `quadrature` | native | Gauss–Legendre via Newton on Legendre roots; adaptive Simpson; Brent | textbook |
| `complex` | native | from the definitions | — |
| chapters 5–9, 19, 22, 25 | **vendored** Cephes | Moshier's rational approximations and continued fractions | upstream |

---

## 4. Test results, and what they actually establish

**256 tests pass workspace-wide; zero failures; zero build warnings;
`cargo clippy --workspace --all-targets` reports zero errors and zero
warnings.**

| suite | count |
|---|---|
| `special_functions` unit tests | 98 |
| `special_functions` vendored-identity tests | 11 |
| doctests | 31 |
| `posim` (incl. 12 special-function bridge tests) | 51 |
| `physical_object` unit | 40 |
| collision integration | 16 |
| conservation integration | 9 |
| **total** | **256** |

Per module: `sph_bessel` 7, `legendre` 10, `orthopoly` 24, `eigen` 8,
`quadrature` 24, `tridiag` 6, `bessel` 6, `wigner` 11, `complex` 2.

### What the tests are designed to catch

The suite deliberately avoids leaning on fixed tolerances against
remembered constants, because **that is exactly what failed repeatedly
during development**. Four times a hand-supplied "expected" value was
wrong while the code was right:

1. a tridiagonal 3×3 solution (`[0.25, 0.5, 1.25]` asserted; the true
   answer is `[0.5, 0, 1.5]`);
2. `legendre_p(2, 5.0)` asserted to be an error, when Pₙ is a polynomial
   defined on all of ℝ and P₂(5) = 37;
3. the 6-j value `{1 1 2; 1 1 2}` asserted as −1/30; it is **+1/30**;
4. a 6-j orthogonality test whose `p` range ignored the (a,d,p) and
   (c,b,p) triangle conditions, reporting a failure at a legitimate zero.

In every case a *structural* check — a residual, an identity, an
orthogonality sum — was already passing and was correct. The lesson is
recorded here because it shaped the suite:

- **Residuals over expected values.** `‖Ax − b‖ → 0` cannot be
  misremembered.
- **Orthogonality and normalisation sums.** For `wigner`, the
  Clebsch–Gordan rows summing to 1 and the 6-j relation
  `Σₓ(2x+1){a b x; c d p}{a b x; c d q} = δₚq/(2p+1)` fix sign *and*
  normalisation absolutely, with no table lookup. They are also the
  tests most sensitive to cancellation in the alternating Racah sum,
  which is the real accuracy risk in that module.
- **Closed forms swept over a range**, not spot-checked. `(j j 0; m −m 0)
  = (−1)^{j−m}/√(2j+1)` is verified across integer and half-integer j
  and every m; `{a b c; 0 c b}` likewise.
- **Convergence laws instead of tolerances.** The particle-in-a-box
  error is asserted to equal `(kπh)²/12`; the oscillator's
  finite-difference error to track `2n²+2n+1`; halving `h` to drop the
  error 4×. These fail if the *physics* is wrong, and pass at whatever
  absolute magnitude the grid dictates.
- **Cross-validation against independent machinery.** The native
  integer-order `bessel_j_array` is checked against the vendored Cephes
  `jv` — two entirely unrelated algorithms — and against the native
  spherical Bessel via the half-integer identity.
- **Mutation testing.** `orthopoly` and `quadrature` were mutation-tested
  to confirm the suites are non-vacuous.

### Observed accuracy, and where the worst error lives

| module | observed | note |
|---|---|---|
| `sph_bessel` | ≤1e-14 against closed forms | Miller recurrence; the upward direction is *proved unstable* by a test rather than merely asserted |
| `legendre` | ≤1e-13 | `assoc_legendre_p` **overflows f64** and returns `Err`; the driver is the ORDER m, via the (2m−1)!! seed — not ℓ, as an earlier draft of the docs wrongly claimed. `norm_assoc_legendre_p` stays O(1) and is the fix |
| `orthopoly` | ≤1e-13 | worst at high degree with large argument, as the recurrences predict |
| `bessel` | ≤1e-10 vs Cephes | **weakest regime:** the seed order must sit well above x. At `x = 45` an under-sized seed gave only ~9 correct digits — caught by the cross-check, fixed, and the measurement recorded in the source |
| `wigner` | ≤1e-12 on orthogonality sums | **weakest regime:** large j, where the alternating Racah sum cancels catastrophically. This is a property of the formula, not the implementation, and it is stated in the module docs |
| `tridiag` | residual ≤1e-11; CN norm drift **1.47e-12** over 6000 steps | no pivoting — stable for diagonally dominant systems only, documented rather than hidden |
| `eigen` | ≤1e-12 | dense, `O(n^3)` — practical to a few hundred rows |
| `lanczos` | residuals ≤1e-8 on 4900-dim problems | matrix-free; **weakest regime:** clustered-but-not-equal eigenvalues, where deflation needs more passes to separate them. Cross-checked against `eigen` on problems small enough for both |
| `quadrature` | degree-exactness to 2n−1, **and the converse** (not exact at 2n, so the bound is sharp) | |

### End-to-end evidence

Unit tests prove the pieces; three examples prove they do the job:

- `scatter_1d` — Gaussian packet off a rectangular barrier, 6000
  Crank–Nicolson steps. Norm drift **1.47e-12** (the Cayley operator is
  unitary for *any* dt, so this measures the solver, not the timestep);
  transmitted fraction within **0.25 %** of the analytic coefficient
  averaged over the packet's own momentum distribution.
- `harmonic_oscillator` — solved two independent ways (analytic Hermite
  eigenfunctions vs finite-difference diagonalisation) that share no
  code.
- `hard_sphere` — the documented expert example, σ → 4πa².

### What is NOT covered

- Complex arguments for the vendored Bessel family (chapter 10). The
  AMOS/TOMS 644 route needs `num-complex`, so it is deferred.
- 9-j symbols (chapter 34).
- Accuracy at very large j in `wigner`, as above.
- The twenty chapters marked **none** in §1.

---

## 5. Language integration

All 30 registered functions are callable from posim; see
[grammar.md](grammar.md) §4.1 and §4.2 and Example 15.

The genuine parse-time work was **argument-domain checking**, not
syntax: the call production already admitted any builtin. An integer
order must be a whole number — `hermite_h(2.5, 1)` is refused rather
than truncated to `hermite_h(2, 1)`, which would return a confident
wrong answer. Angular momenta are the deliberate exception: they may be
half-integral, so the wigner entry points take plain numbers and
validate in the library.

`Value::Complex` and the imaginary literal `3i` were added to reach the
complex Crank–Nicolson solvers — a real lexer, parser and VM change
rather than a registration.

**The lockstep rule is now mechanical.** A test asserts that every
registered name appears in `HELP_TEXT`, the `parser.rs` EBNF comment,
`grammar.md` and `grammar.tex`; the build fails otherwise. It was
verified to actually fail by registering an undocumented name. This
replaced enforcement-by-discipline, which had already let the entire
integration be skipped once.

---

## 6. Honest summary

**Pros.** Everything present is tested by structural properties rather
than remembered constants, and cross-validated where two independent
routes exist. Zero `unsafe`, zero external dependencies, zero warnings,
zero clippy findings. The clean-room replacements remove real licensing
blockers. The mathematical infrastructure for 1-D quantum mechanics —
complex arithmetic, a unitary propagator, a dense eigensolver,
quadrature, the orthogonal polynomial families — is present and
demonstrated end to end.

**Cons.** Coverage is ~40 % of chapters and less than that of content.
Complex-argument Bessel, 9-j symbols, and every one of the twenty
untouched chapters are absent. `wigner` degrades at large j and
`assoc_legendre_p` overflows at large order — both documented, neither
fixed. The vendored dependency has an unresolved LICENSE-file question
upstream. Two-dimensional quantum problems will need a different
linear-algebra path, since a full 2-D Crank–Nicolson operator is not
tridiagonal.

**The claim being made is narrow and, I believe, defensible:** this is a
tested, documented, licence-clean foundation sufficient for the 1-D
quantum mechanics the SolveIt port needs. It is not the DLMF.
