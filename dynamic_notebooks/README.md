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

Launch any of them by bare name with the shell alias (see the
repository README for the Windows 11 / macOS / Linux installs), or
with cargo directly:

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

## The catalogue

| notebook (launch: `posim_notebook <name>`) | what you watch | anchor |
|---|---|---|
| `kepler_orbit` | an e = 0.6 ellipse; fast perihelion, slow aphelion | E = −0.5 exactly; LAPLACE 1 reads [0.6, 0, 0]; period 2π |
| `outer_solar_system` | Sun + the four giants + Pluto (AU/day/solar-mass units) | at 50 days/tick Jupiter laps in ~3 s, Pluto in ~1 min; E = −3.2155e-8 |
| `thrown_ball` | the textbook parabola, camera framed on the full arc | x(t) = 30t; apex 46.9 at t = 3.06; range 183.6 at t = 6.12 |
| `charged_in_b_field` | cyclotron circle of a negative charge | T = 2π·2/(1.5·4) ≈ 2.094; gyroradius mv/(\|q\|B) = 1 |
| `tumbling_body` | the tennis-racket (Dzhanibekov) flip | moments [5, 4.25, 1.25]; intermediate-axis spin flips near t ≈ 40 |
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
| `two_dumbbells` | two user-function dumbbells colliding off-center | hud E, P **and** L identical before and after the impact |

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
