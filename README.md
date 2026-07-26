# rustSimulate

Pure-Rust physics simulator. **All** numerical integration runs through
[`sundials_rs/`](sundials_rs) — a pure-Rust translation of SUNDIALS
7.7.0 vendored in this repository (CVODE Adams/BDF, ARKODE symplectic
SPRK). Zero `unsafe`, zero external crate dependencies, zero warnings.

The repository is **self-contained**: an ordinary clone is all you need
— no submodules, no network access during the build, nothing from
crates.io.

```bash
git clone https://github.com/once-ere/rustSimulate.git
cd rustSimulate
cargo run                 # the notebook REPL (type HELP)
```

This project is a standalone export of the simulator; its lineage,
byte-identity manifest and full verification transcript are recorded in
[EXPORT_PROVENANCE.md](EXPORT_PROVENANCE.md).

Latest release: a programmable notebook and a compound rigid body —
user-defined functions (`DEF name(param = default, ...) { body }`,
every body line syntax-checked at definition, `FUNCS`/`SHOW` to list
and edit), named objects (`NEW ... AS name`, plus `LET` variables and
string literals), and the rigid `DUMBBELL` (two solid spheres plus a
rod as ONE rigid body, exact part-wise collisions conserving E, P and
L through real solver events); the scene window gains a permanent
Reset button (with `SCENE RESET` — bit-identical re-initialization,
Start re-runs) and a live labeled conserved-quantities readout (E, P
and L); 104 tests green (40 lib + 16 collision + 9 conservation +
39 posim).

- `physical_object/` — library: `pub struct physical_object`, the
  unique union of the legacy `PointParticle`, `RigidBody` and
  `RigidBody3D`, with get/set for every field; `PhysicalObjectSystem`;
  the sundials integration drivers; validated examples.
- `posim/` — the simulator front end: lexer → grammar compiler → stack
  machine, a notebook REPL (`In[n]`/`Out[n]` cells), script batch mode,
  and a JSON machine mode.
- `jupyter/` — JupyterLab wrapper kernel so notebooks can get/set the
  simulator's data (see `jupyter/README.md`).
- `sundials_rs/` — the numerical engine: a pure-Rust, zero-`unsafe`,
  dependency-free translation of SUNDIALS 7.7.0, vendored here as a
  self-contained workspace (read-only; upstream any changes).
- `PLAN.md` — the integration plan / design record (union mapping,
  grammar, solver mapping, verification results).

## Documentation

- [grammar.md](grammar.md) / [grammar.pdf](grammar.pdf) — the complete
  command-language and notebook specification, with fourteen worked
  examples.
- [physical_object_simulator.md](physical_object_simulator.md) /
  [physical_object_simulator.pdf](physical_object_simulator.pdf) — the
  full solution guide for new users, with fourteen more worked examples.
- [scene_info.md](scene_info.md) / [scene_info.pdf](scene_info.pdf) —
  the graphical scene window: the simulator research survey, the
  protocol, and the UI.
- [collision_detection.md](collision_detection.md) /
  [collision_detection.pdf](collision_detection.pdf) — the collision
  science reference, with documented example scripts in
  `scripts/collisions/` (01–12).
- [ARCHITECTURE.md](ARCHITECTURE.md) — module responsibilities and
  pinned cross-module contracts.
- [CLAUDE.md](CLAUDE.md) — working rules for contributors and agents.

## Quick start

```bash
cargo run                 # notebook REPL (type HELP)
cargo test --workspace    # all tests
cargo run -p physical_object --release --example kepler_orbit
cargo run -p physical_object --release --example outer_solar_system
cargo run -p physical_object --release --example tumbling_body
cargo run -p physical_object --release --example charged_in_b_field
cargo run -p posim -- --script scripts/collisions/12_two_dumbbells.posim
cargo run -p posim -- --script my_session.posim
cargo run -p posim --release -- --notebook dynamic_notebooks/kepler_orbit.posim
cargo run -p posim -- --machine   # JSON protocol for front ends
```

Example session:

```
In[1]:= new sphere { mass = 2, radius = 0.5, position = [0, 10, 0], velocity = [1, 0, 0] }
Out[1]= obj0
In[2]:= set system.gravity = [0, -9.81, 0]
In[3]:= step 1
Out[3]= t = 1 (advanced by 1, 12 solver steps)
In[4]:= get obj0.position
Out[4]= [1, 5.095000000000006, 0]
In[5]:= method sprk leapfrog_2_2 0.001
In[6]:= help
```
