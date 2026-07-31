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

## The Index of Functions

`index_of_entities.html` is a browsable catalog of **every named entity in this
repository** — 4,830 of them — with a definition, the `file:line` where it is
defined, its complete syntax, and examples you can paste into a notebook and
run. Open it directly; it needs no server and fetches nothing.

```bash
open index_of_entities.html          # macOS  (xdg-open on Linux)
```

Keep `catalog-c.js` beside it: the page loads that second payload on demand
when you open a bucket or search, which is what keeps the main file at 1.2 MB
instead of 5.

| tier | what | entries | examples |
|---|---|---|---|
| A | the notebook surface — commands, keywords, field paths, builtins, types, notebooks, worked examples | 452 | 1,141 posim + 24 machine-mode fragments, **all executed** |
| B | the first-party Rust API — `physical_object`, `special_functions`, `quantum`, `posim` | 491 | 177 snippets, **all compiled** |
| C | the vendored `sundials_rs` workspace | 3,887 | linked to the shipped example programs that call them |

**Every example is checked, and the page says which kind of check it got:**
a posim fragment is *executed* (`posim --script`, output captured), a
machine-mode fragment is *executed* through the JSONL protocol, a Rust snippet
is *compiled*, and a shell command reads "run this yourself" because it points
at a program verified elsewhere rather than one this page ran. 1,342 of 1,342
runnable examples pass.

The 68 documented worked examples now carry their real transcripts. Where the
typed half replays, it is **executed** and shown beside the document's own
published output so you can compare; where it cannot — a deliberate refusal
that *is* the lesson, an elided session, or a `SCENE REVERSE` that depends on
wall-clock playback — the transcript is quoted and the reason is stated.

Navigation is A–Z / 0–9 / special-character buckets, with search and kind
filters. The whole thing is operable from the keyboard — press <kbd>?</kbd> for
the map — and has BACK/FORWARD history that agrees with the browser's own
buttons. You can edit, add and hide entries; changes live in a local overlay
you can export and re-import, and the generated catalog is never mutated.

### What it does not claim

The status page inside the app is the authority, and it is blunt:

- **53 Tier-A/B entries are stubs** — down from 390. Each is catalogued with a
  definition and a location, but has no example and no call site outside its
  own file. Three are prose-only sections of `special_functions.md` that carry
  no code block at all. A further 264 carry status `reference`: no generated
  snippet, but a link to the tests, examples or documents that genuinely call
  them.
- **3,887 Tier-C entries carry status `reference`**, not `complete`.
  `sundials_rs` is a faithful translation of a C library whose API is
  `&mut CVodeMem` plus a context, a matrix, a linear solver and callbacks; a
  one-line snippet would misrepresent how any of it is reached. Instead each
  entry links to the shipped example *programs* that actually call it — those
  are diffed byte-for-byte against the upstream C references
  ([sundials_rs/VERIFICATION.md](sundials_rs/VERIFICATION.md)).
- **All 254 Tier-A commands, builtins, field paths and types carry the full
  four-rung ladder.** Closing the last 24 meant giving the machine-mode JSON
  ops real executed examples (`posim --machine` is as runnable as the
  notebook), and giving each shape a rung that *checks* its documented inertia
  formula against the printed tensor rather than restating it.
- Captured output has genuine run-to-run variation normalised — the
  OS-assigned scene port becomes `<port>`, playback counters become
  `<varies>`. Solver step counts are **not** normalised: those are
  deterministic, and they are the anchors the documentation pins.

[index_data/DIVERGENCES.md](index_data/DIVERGENCES.md) records six places where
the prose and the code disagreed, each settled by running a probe rather than
by reading. One was a defect in code-adjacent material and has been fixed with
a test that now gates it in both directions; the rest are documentation drift,
and the index carries the code's behaviour.

### Rebuilding it

```bash
python3 tools/extract_rust_items.py > index_data/rust_items.jsonl
python3 tools/build_commands.py && python3 tools/build_tierb.py && python3 tools/build_tierc.py
python3 tools/build_catalog.py && python3 tools/build_app.py
python3 tools/verify_index_examples.py      # runs every posim fragment
python3 tools/verify_tierb_examples.py      # compiles every Rust snippet
```

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
