# dynamic_notebooks/ — one live-scene notebook per documented example

A **dynamic notebook** is a notebook that opens a GUI that executes a
simulation: a `.posim` file that builds a system, prints its analytic
baselines, and ends with `SCENE CREATE` — so launching it opens the
graphical scene window with the simulation loaded, Stopped, and ready.
Press **Start** in the window (or type `SCENE START` at the prompt) to
run and display it; the terminal stays in the interactive notebook, so
you can `GET`/`SET`, `SCENE PAUSE`/`REVERSE`/`RESET`, or extend the
session — the loaded cells keep their `In[n]` numbers and your next
command continues the numbering.

Launch any of them by bare name with `tools/posim_notebook`, or with
cargo directly. The launcher resolves the name against this directory,
works from any working directory, and takes `--list` to enumerate what
is here. To get the bare `posim_notebook` form, put it on your PATH:

```bash
ln -s "$PWD/tools/posim_notebook" /usr/local/bin/posim_notebook
```

Then either form works:

```bash
posim_notebook kepler_orbit
```

```bash
cargo run -p posim --release -- --notebook dynamic_notebooks/kepler_orbit.posim
```

Every notebook was executed and verified: it loads with zero failing
cells, its baseline outputs match the documented analytic values, and
the scene playback genuinely advances the simulation (checked
headlessly via `SCENE START` → `SCENE STATUS`, and live in a browser
for the dumbbell impact — the hud's conserved E, P, L read identically
before and after the collision).

`MANIFEST.md` beside this file is the audit trail: where each notebook
came from, and the result of actually running it.

## The catalogue

| notebook (launch: `posim_notebook <name>`) | what you watch | anchor |
|---|---|---|
| `kepler_orbit` | an e = 0.6 ellipse; fast perihelion, slow aphelion | E = −0.5 exactly; LAPLACE 1 reads [0.6, 0, 0]; period 2π |
| `outer_solar_system` | Sun + the four giants + Pluto (AU/day/solar-mass units) | at 50 days/tick Jupiter laps in ~3 s, Pluto in ~1 min; E = −3.2155e-8 |
| `thrown_ball` | the textbook parabola, camera framed on the full arc | x(t) = 30t; apex 46.9 at t = 3.06; range 183.6 at t = 6.12 |
| `charged_in_b_field` | cyclotron circle of a negative charge | T = 2π·2/(1.5·4) ≈ 2.094; gyroradius mv/(\|q\|B) = 1 |
| `tumbling_body` | the tennis-racket (Dzhanibekov) flip, over and over | moments [5, 4.25, 1.25]; **body-frame** w_y reverses about every 12-13 units (t ≈ 4, 18, 28, 40, …) while world-frame w_y never does; L and E pinned |
| `three_bodies` | a 1–4–256 gravitational trio | total momentum exactly [0, 0, 0] (type `del 2` live to break it) |
| `magnetic_spin_up` | torque pumping spin into a sphere | L(t) = 0.1t exactly; hud E grows 0 → 0.8 by t = 4 — honestly |
| `charged_in_e_field` | uniform-field acceleration | x(t) = 0.375 t²; x(4) = 6 exactly |
| `static_anchor` | a pinned body ignoring gravity beside a falling one | inverse_mass = 0 pins the anchor at y = 10 forever |
| `newtons_cradle` | five spheres; only the far ball exits | v₄ = 1 exactly, the rest stop; watch the contact arrows chain |
| `bouncing_ball_restitution` | an e = 0.8 bounce | impact t = √0.9 ≈ 0.949; rebound apex 3.38 = e²·4.5 + 0.5 |
| `head_on_exchange` | equal masses swap velocities | touch at t = 1.5; E = 1, P = [0, 0, 0] throughout |
| `unequal_masses` | 1-vs-3 head-on | outgoing v₁′ = −1/2, v₂′ = +1/2 (textbook 1-D formulas) |
| `restitution_ladder` | e = 1.0/0.8/0.5/0.2 in four side-by-side lanes | separation speed = e × 2 per lane, all striking at once |
| `billiard_break` | a cue into a two-ball rack | momentum [2, 0, 0] conserved; exits along lines of centers |
| `spin_up` | an off-center hit spins up a box | L = r × J = (0, 0, −2.22) appears from nothing linear |
| `thin_wall_toi` | a 100 u/s bullet vs a 1 cm plate — no tunneling | impact at t = 0.01985 exactly; velocity flips to −100 |
| `colliding_binary` | gravity + collision: a binary bouncing at pericenter | pericenter 0.5 < touch 0.6; e = 0.6 bounce shrinks the orbit |
| `spinning_target` | hitting a tumbling box whose surface is moving | impact near t ≈ 0.97; the K-matrix angular terms do the work |
| `billiard_box` | an elastic ball ping-ponging between static walls | E = 1.445 forever; a wall strike every ≈ 2.94 t-units |
| `box_of_shapes` | the manager's demo: all six body types in BOX 4 | E₀ = 30000 exactly; the point threads the torus hole |
| `box_of_shapes_m32` | the same box, but the only mover is a mass-32 point | E₀ = ½·32·60000 = 960000 exactly; E conserved, P **not** — the infinite walls absorb it |
| `two_dumbbells` | two user-function dumbbells colliding off-center | hud E, P **and** L identical before and after the impact |
| `routh_double_star_period` | three binaries of one total mass but ratios 1:1, 3:1, 19:1 — the lopsided one orbits in the scene | Routh I.Art.400: T = π√2 = 4.442882938158366 closes every ratio to ~1e-11; total mass 2.5 misses by 0.696 |
| `routh_rectangle_diagonal` | a 2×1 plate spun about one diagonal, tumbling onto the other | Routh II.Art.150b Ex.1: T = 2K(sin α)/√cos 2α = 4.285129705594264; ω_body y-component flips sign to 1e-12 |
| `routh_p1_two_trajectories` | two particles leave one point at one speed; one cuts across, one loops the long way, both hit the target | Routh I.Art.339: a = 1.476014760147601; arrivals at t = 2.16112 and t = 8.89752, both within 1.4e-6 of Q |
| `routh_p1_kepler_equation` | the e = 0.6 ellipse: a fast swing through perihelion, a long crawl to aphelion | Routh I.Arts.342-346: u - 0.6 sin u = m predicts r to 1e-11 at t = 1, 2, 3; equation of the centre falls to 0.097 near aphelion |
| `routh_p1_equilateral_three_body` | a triangle turning as though welded — then the same triangle, one corner nudged 0.001, coming apart | Routh I.Arts.407-408: omega = sqrt(3), period 3.6275987284684357, all three home to 3e-10; Art.412 unstable (1/3 vs 1/27) — sides 0.405/1.802/1.409 after three periods |
| `routh_p2_spin_stabilisation` | three identical boxes given identical sideways kicks; the slow one lurches, the fast one barely nods | Routh II.Art.156: same d = 0.3 at n = 1, 3, 9 gives wobble 20.171, 8.326, 2.854 deg — each cosine bit-identical after 10 units |
| `routh_p2_two_quadrics` | a box tumbling irregularly while the hud's E and L do not move at all | Routh II.Arts.140-142: G = 12.838540610209558 bit-identical and E = 19.35625 to 12 digits, while |w| wanders 3.017 -> 3.407 |
| `routh_p2_fixed_couple` | a box winding itself up from dead rest into an ever faster tumble | Routh II.Art.148: L(t) = (0,0,2t) exactly — (0,0,10) at t=5, (0,0,20) at t=10 — while energy climbs 28.3 -> 113.3 |
| `routh_p2_mean_axis_instability` | three identical boxes spun about their three principal axes; only the middle one tumbles end over end | Routh II.Art.155: growth rate k = 1.8; body-frame w_y reverses every ~8 units, while the stable axes stay within 0.005 and 0.018 of theirs |
| `routh_p1_three_projectiles` | three particles on wildly different parabolas; the triangle they span swells but never tilts | Routh I.Art.158 Ex.1: unit normal (0.1788, -0.9298, -0.3218) identical at t = 1, 2, 4 to 1e-15 |
| `routh_p1_expanding_sphere` | a shell of particles blowing outward while the whole shell falls | Routh I.Art.167 Ex.1: radius V*t about the free-fall centre — 5,5,5 at t=1 and 10,10,10 at t=2 |
| `routh_p1_equal_periods` | four orbits of obviously different shape sweeping home together | Routh I.Art.335: one speed 1.1 -> one a = 1.2658227848101269 -> T = 8.948273124536605; all return within 1e-9 |
| `routh_p1_oblique_impact` | an off-centre clip that sends both spheres away at a perfect right angle | Routh I.Art.89: dot(v0,v1) = -8.3e-17; speeds 0.5 and 0.8660254; |dE/E| = 2.2e-16 |
| `routh_p1_escape_velocity` | three particles fired straight out at 0.9, 1.0 and 1.1 times escape speed | Routh I.Arts.312+335: specific energy -0.19 / exactly 0 / +0.21, held to 1e-11; only the first returns |
| `routh_p1_sphere_exchange` | a light sphere hits a heavier one and stops absolutely dead | Routh I.Arts.85+87: equal masses exchange velocities; m = e*m' gives v = 0, v' = 0.5, momentum 1 |
| `routh_p1_geometric_progression` | a pulse walking down four balls of doubling mass | Routh I.Art.88: velocities 1, 2/3, 4/9, 8/27 and rebounds -1/3, -2/9, -4/27; |dE/E| = 1.1e-16 |
| `routh_p1_centre_of_gravity` | four bodies attracting, colliding and scattering while their mean position slides in a straight line | Routh I.Art.92: through 8 impacts the centre of gravity stays 5e-15 from P/M*t |
| `routh_p1_hodograph_circle` | an e = 0.6 ellipse whose velocity vector never leaves a circle | Routh I.Arts.394-398: hodograph radius mu/h = 1.25 at every sample while the speed runs 2 -> 0.5 |

### Not scene notebooks

Two `.posim` files here are quantum notebooks. They have **no**
`SCENE CREATE` — they visualise by writing an HTML file you open
yourself — so they are listed apart from the catalogue above:

| notebook | what it produces | anchor |
|---|---|---|
| `double_slit` | `double_slit.html` — 2-D `QM2` wavepacket through a two-slit wall built from comparison operators | slit separation d = 4, λ = 2π/k = 0.785; fringes come purely from amplitudes adding |
| `tunneling` | `scatter.html` — 1-D `QM` packet against a barrier written as an ordinary user function | k₀ = 2 so E₀ = 2, **below** the 2.5 barrier; classically nothing gets through |

They still launch the same way: `posim_notebook double_slit`.

## Which documented example is which notebook

Every example in the documentation maps here:

- **Rust self-checking examples** (`physical_object/examples/`):
  `kepler_orbit`, `outer_solar_system`, `tumbling_body`,
  `charged_in_b_field`, `newtons_cradle`, `bouncing_ball_restitution` —
  same names.
- **Collision scripts** (`scripts/collisions/01–12`, documented in
  `collision_detection.md` §9): `head_on_exchange` (01),
  `unequal_masses` (02), `restitution_ladder` (03, adapted — the
  script's sequential RESET rungs run here as four simultaneous
  lanes), `newtons_cradle` (04), `billiard_break` (05), `spin_up`
  (06), `thin_wall_toi` (07), `colliding_binary` (08),
  `spinning_target` (09), `billiard_box` (10), `box_of_shapes` (11),
  `two_dumbbells` (12).
- **grammar.md §9 worked examples**: Ex1/Ex11 → `kepler_orbit`,
  Ex2 → `thrown_ball`, Ex3 → `charged_in_b_field`,
  Ex4 → `tumbling_body`, Ex5 → `three_bodies`,
  Ex6 → `magnetic_spin_up`, Ex13 → `box_of_shapes`,
  Ex14 → `two_dumbbells`. Ex7 (%edit), Ex8 (vector algebra),
  Ex9 (hand-built tensors), Ex10 (%save/%load), Ex12 (camera
  driving) are language/notebook-mechanics demos with no simulation
  to display — Ex10's mechanism *is* the `--notebook` loader itself,
  and Ex12's camera commands appear inside these notebooks' framing
  cells.
- **User-guide §8 examples (S1–S14)**: S2 → `kepler_orbit`,
  S4 → `charged_in_e_field`, S5 → `static_anchor` (adapted — a free
  test body added beside the pinned anchor so the contrast is
  visible), S11 → `tumbling_body`, S13 → `box_of_shapes`,
  S14 → `two_dumbbells`. S1/S3/S6/S7/S10 are Rust API
  demonstrations and S8/S9/S12 are protocol/JupyterLab
  demonstrations — no notebook-language simulation to display.

(`outer_solar_system`'s 500,000-day certification — Pluto's position
to eight decimals — is the self-checking Rust example's job; the
notebook is the same system watched live at 50 days per tick.)

## Conventions

Every notebook follows one shape: a header stating the physics, the
numeric anchors and both launch lines; the example's **exact** setup
numbers; baseline observable cells (`energy`, `momentum`, `angmom`,
`laplace`, `get …`) whose printed values you can check against the
header before anything moves; then `scene create` plus, where needed,
`scene set_time_step` (so the interesting span plays out in seconds to
a minute of wall clock) and camera framing (`scene translate` /
`scene zoom`). Any deliberate deviation from a source example (the
ladder's lanes, the anchor's companion) is called out in that
notebook's header. No notebook runs the simulation in batch — that is
what the window's Start button is for.
