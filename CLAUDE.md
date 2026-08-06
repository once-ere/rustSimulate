# CLAUDE.md — physical_object_simulator workspace

Pure-Rust physics simulator. The **only** numerical-integration backend
is the vendored `sundials_rs/` workspace (pure-Rust SUNDIALS 7.7.0) —
no exceptions. Treat `sundials_rs/` as a **read-only vendored library**:
it is byte-identical to `once-ere/sundials_rs@faabb7f` and that identity
is a verifiable property (see `EXPORT_PROVENANCE.md`) — never edit it
in place; upstream it instead.

This repository is a standalone export of the simulator (see
`EXPORT_PROVENANCE.md`). The legacy donor sources it was refactored
from (`solver.rs`'s `PointParticle`, `RigidBody.rs`, `RigidBody3D.rs`)
and the 61 MB SUNDIALS C reference tree were deliberately **not**
exported; the union mapping they produced is recorded in `PLAN.md`.

## Commands

- Build: `cargo build --workspace 2>&1 | tee /tmp/build.log`
- Tests: `cargo test --workspace 2>&1 | tee /tmp/test.log`
- Notebook: `cargo run` (type `HELP`); batch: `cargo run -p posim -- --script <f>`
- Dynamic notebook (loads a file, opens its scene window, stays
  interactive): `cargo run -p posim --release -- --notebook
  dynamic_notebooks/<name>.posim` — one per documented example, see
  `dynamic_notebooks/README.md`
- Scene window: type `SCENE CREATE` in the notebook (opens a browser
  page; `SCENE START/PAUSE/REVERSE/RESET`, arrows/drag/wheel in the
  window).
  Headless runs: set `POSIM_NO_BROWSER=1` to suppress `xdg-open`.
- Self-checking physics examples:
  `cargo run -p physical_object --release --example
  {kepler_orbit|outer_solar_system|tumbling_body|charged_in_b_field|newtons_cradle|bouncing_ball_restitution}`
  (each prints SUCCESS/FAILURE and exits nonzero on failure)
- Collision example scripts: `cargo run -p posim -- --script
  scripts/collisions/NN_name.posim` (01–12; documented with captured
  output in collision_detection.md §9)
- Wire protocol test: `python3 jupyter/test_protocol.py` (needs
  `cargo build --release` first — it prefers `target/release/posim`,
  so a stale release binary silently shadows your debug build; stdlib
  only, covers the SCENE command family and the `events` op)
- Kernel test: `POSIM_NO_BROWSER=1 jupyter/.venv/bin/python jupyter/test_kernel.py`
  (needs `uv venv jupyter/.venv && uv pip install -p jupyter/.venv/bin/python ipykernel jupyter_client`)
- Docs: `pdflatex -interaction=nonstopmode grammar.tex` (run twice);
  same for `physical_object_simulator.tex`, `scene_info.tex` and
  `collision_detection.tex`. Keep `.md` and `.tex` versions in sync —
  the `.md` is the source of truth.

## Layout

- `physical_object/` — library: `linalg` (Vec3/Mat3/Quat + skew/outer),
  `boundary` (enum + `Sdf` trait), `physical_object` (the union
  struct, get/set), `system` (collection + 13N pack/unpack),
  `integrate` (**all** time integration: CVODE Adams/BDF + ARKODE
  SPRK, with collision event rootfinding armed when collidable pairs
  exist), `collide` (contact geometry + impulse response; conventions
  pinned in ARCHITECTURE §3.8 — normal points i→j, separation > 0
  apart and IS the root function).
- `posim/` — binary: `lexer` → `parser` (EBNF in its header) → `vm`
  (stack machine) → `notebook` (REPL) / `machine` (JSONL) / `scene`
  (graphical window: `mod.rs` server+playback, `ws.rs` hand-rolled
  SHA-1/base64/RFC 6455, `scene.html` embedded page).
- `jupyter/` — Python wrapper kernel (outside the Rust constraints);
  a reader thread streams async `{"event":...}` lines to the notebook;
  `.venv/` and `.kernels/` are gitignored scratch.
- `ARCHITECTURE.md` — pinned cross-module contracts (state layout,
  solver driving, setter invariants, wire protocol, scene subsystem
  §3.7). **Read it before touching anything cross-module; update it
  when a contract moves.**
- `PLAN.md` — design record (union mapping table, decisions).
- `grammar.md`/`.pdf` — the command-language spec; user docs in
  `physical_object_simulator.md`/`.pdf`; the scene window + the
  seven-simulator research survey in `scene_info.md`/`.pdf`; the
  collision science reference (research, chart/tree, porting
  evaluation, recommendation, 12 examples) in
  `collision_detection.md`/`.pdf`.
- `.backups/` — pre-modification file backups (gitignored). Back up
  before modifying; this repo's git history is the real undo.

## Hard rules

1. **Sundials-only integration.** All stepping goes through
   `physical_object/src/integrate.rs` calling `cvode_rs`/`arkode_rs`.
   Never add a hand-rolled Euler/Verlet/RK stepper anywhere — including
   in examples, docs, and the scene playback thread (it calls
   `integrate::step`; reverse is snapshot replay from the history ring,
   never negative-dt integration).
2. **Zero `unsafe`, zero external dependencies, zero warnings.** Every
   crate root carries `#![forbid(unsafe_code)]`, `#![deny(warnings)]`,
   and allows `non_snake_case`, `non_camel_case_types`,
   `non_upper_case_globals`. `Cargo.lock` must list only local crates
   (`physical_object`, `posim`, `sundials_core`, `cvode_rs`,
   `arkode_rs`); if it grows, you broke the dependency rule. This
   applies to networking too: the scene server's HTTP/WebSocket/SHA-1/
   base64 are hand-rolled on `std::net` in `posim/src/scene/ws.rs`,
   and `scene.html` uses vanilla JS + canvas (no CDN fetches).
3. **Fidelity to donor physics.** Formulas ported from the legacy types
   (softened gravity, Laplace vector incl. its `r = 0` guard, inertia
   formulas, kinetic energy, SDF central-difference normal) keep their
   arithmetic order — floating point is not associative.
4. **All state access through get/set.** The VM, machine mode, the
   scene subsystem, and any new front end reach `physical_object`
   fields only via the setters (they enforce the coupled invariants:
   mass↔inverse, inertia↔inverse, unit quaternions,
   momentum-canonical velocity).
5. **Missing sundials symbols are reported, not invented.** If an API
   you need is absent from `sundials_rs/`, stop and say exactly which
   symbol is missing, naming the C original's file (`src/` or
   `include/` of the upstream SUNDIALS 7.7.0 release — that reference
   tree is not vendored here). Do not reimplement solver numerics
   locally.
6. **Quaternions are w-first everywhere** (packing, VM literals, JSON,
   scene frame messages, docs). The 13-per-object state layout in
   `system.rs` is a pinned contract (ARCHITECTURE.md §3.2).
7. Public error style: `Result<_, String>` with actionable messages
   naming the field/column/feature. No panics in library code paths;
   scene threads turn bad window input into error events, never crashes.
8. **The struct name is `physical_object`** (lower-case, by
   specification). Import as
   `use physical_object::physical_object::physical_object;`; in files
   that import it, prefix other crate paths with `::`. Do not rename it
   and do not try to re-export it at the crate root (namespace clash).
9. **The scene evolves a copy.** The playback thread owns a synced
   clone of the system (`SCENE CREATE`/`SCENE REFRESH`/`RESET` sync
   it); notebook `STEP`/`RUN` never move the window and window playback
   never moves the notebook. Do not "fix" this by sharing the system —
   it is the isolation that keeps the VM lock-free (ARCHITECTURE §3.7).

## Workflow

- **Make backups before modifying files** (copy into
  `.backups/<date>/`), then edit.
- After EVERY build/test/run command: `2>&1 | tee <log>`, then read the
  log before editing. ≤2 attempts per failing command, then switch
  strategy.
- Commit after every coherent file group; keep
  `cargo build --workspace` warning-free and `cargo test --workspace`
  green at every commit (568 tests: 40 physical_object lib +
  19 collision + 9 conservation + 109 posim + 92 quantum +
  233 special_functions + 11 vendored identities + 55 doctests).
  Phase gates are tagged (`phase-posim-green`).
- New solver features need: a unit or conservation test with an
  analytic expectation, a grammar hook if user-facing (lexer keyword →
  parser production → VM instruction → HELP_TEXT → grammar.md/.tex),
  and an ARCHITECTURE.md update if a contract changed.
- When touching the language, update **all four** in lockstep: the EBNF
  comment in `parser.rs`, `HELP_TEXT` in `vm.rs`, `grammar.md`, and
  `grammar.tex` (then recompile the PDF). Scene-visible behavior also
  updates `scene_info.md`/`.tex`.
- Scene argument parsing is **term-level** (so `scene rotate 15 -5` is
  two arguments); if you change it, `parser::scene_command` and the
  grammar docs move together.
- Verified physical facts to protect (regression anchors): outer solar
  system matches the donor `solar_system` example (Pluto to 8 decimals
  at t = 500,000 days, energy drift < 1e-6); Kepler e = 0.6 conserves
  E/L/Laplace < 1e-6; gyroradius = mv/(qB) to 1e-4; tumbling body
  conserves L exactly; scene reverse replays history to **exactly**
  t = 0 (bit-identical snapshots, `playback_forward_then_reverse_
  restores_state`); collision analytic suite green (exchange TOI
  < 1e-9, thin-wall no-tunnel at 100 m/s vs 5 mm plate, apex = e²h,
  cradle propagation) and zero-collidable-pair systems **bit-identical**
  with COLLIDE ON vs OFF (structural invariance — protect it);
  the tilted-torus fit (an axis-aligned torus of outer radius 2
  exactly inscribes BOX 4; tilted to axis (1,1,1)/√3 its per-axis
  extent is 1.5·√(2/3) + 0.5 ≈ 1.7247, clearing every wall —
  `tilted_torus_fits_a_4x4_box_where_flat_does_not`); a ball rattling
  in a rigid BOX conserves E with the walls **bit-identically at
  rest** (`ball_in_a_rigid_box_conserves_energy_and_walls_never_move`);
  a point threads the torus hole while the fat ball bounces with
  TOI = (3−√1.45)/4 (`point_threads_the_torus_hole_but_a_fat_ball_
  bounces`); the mixed-shape box rattle holds |dE|/E < 1e-6 through
  ~50 events (`mixed_shapes_rattle_in_the_box_conserving_energy`);
  and the BOX grammar family (BOX <size>/OFF/status, `system.box`,
  `[wall: static, inverse_mass=0]` LIST tags, machine
  `box`/`wall`/`inverse_mass`, scene init `box` + wall flags)
  round-trips end-to-end (`box_family_and_infinite_mass_walls`,
  `state_reports_box_walls_and_inverse_mass`,
  `box_shapes_and_wall_flags_reach_the_init_message`);
  parallel round shapes touch side-on at the exact lateral gap —
  the radial rejection axes make cylinder-cylinder side contacts
  exact (`parallel_cylinders_side_contact_is_exact`); a ball released
  FROM REST above a thin plate is caught at the analytic TOI to 1e-6
  — the anti-tunneling cap is acceleration-aware
  (`ball_released_from_rest_does_not_tunnel_the_thin_plate`); two
  disks with PARALLEL planes have separation |dz|, which touches zero
  without a sign change — face-on disk-disk crossings are invisible
  to downward-crossing rootfinding, a KNOWN limitation pinned
  deliberately
  (`parallel_disk_disk_separation_is_the_documented_limitation`; tilt
  one disk or model a thin cylinder); NEW is transactional (a failing
  initializer or final validation leaves no ghost object) and the
  torus inner/outer pair is genuinely order-independent, horn torus
  `inner_radius = 0` valid
  (`torus_pair_is_order_independent_and_new_is_transactional`); and
  `BOX <size>` after a wall deletion removes the surviving tracked
  slabs before building the new box — no orphan leak
  (`box_recreate_after_wall_deletion_leaks_nothing`); two tumbling
  dumbbells colliding off-center conserve E, P **and L (about the
  origin)** to 1e-8 through real CVODE events — the part-wise exact
  narrow phase puts the impulse pair at one shared contact point
  (`colliding_dumbbells_conserve_energy_momentum_and_angular_
  momentum`); scene Reset (toolbar button / window `reset` /
  `SCENE RESET` — one primitive) restores the playback's initial
  state **bit-identically**, clears history and the step counter,
  returns the mode to Stopped, and Start re-runs from the beginning
  (`reset_restores_the_initial_state_and_start_reruns`); the
  `create_dumbell` flow — define/call/members/shorthands/renumber/
  errors/redefine/LET-defaults/ghost-free failing calls —
  round-trips end-to-end
  (`def_call_named_objects_and_dumbbell_members`); and the directed
  supports are exact for BOTH dumbbell ends — asymmetric wall gaps
  exact, the light end pokes farther (|z1| > |z2| when m2 > m1),
  ball-vs-pole and ball-vs-rod contacts exact
  (`dumbbell_wall_gaps_and_ball_contacts_are_exact`,
  `dumbbell_constructor_com_sdf_and_supports`); dumbbell part masses
  whose SUM overflows f64 are refused, and a fat-rod dumbbell
  (rod_radius > both sphere radii) reports honest side-on supports —
  rank 1 with footprint = the rod half-length, not a sphere point
  (`fat_rod_dumbbell_reports_honest_support_ranks`); DEF is hardened:
  parameters cannot be reserved words, builtins, pi/tau or duplicates,
  defaults are expressions only (a command like `reset` is refused at
  definition time), a `}` inside a `#` comment does not close the
  body, `NEW ... AS n` errors when `n` is bound to a number, and a
  function that runs NEW inside another NEW's initializer
  stashes/restores the outer NEW context — both objects exist, a
  failing inner call rolls the outer back
  (`def_hardening_reserved_words_duplicates_defaults_and_nesting`);
  and contact records arriving with a notebook system never reach the
  scene playback copy — cleared at create, refresh-sync and reset
  (`stale_notebook_contacts_never_reach_the_playback_copy`).
- Verified UI facts to protect (re-check with a headless-Chrome CDP
  session after touching `scene.html`): arrow keys translate the view
  right/left (and up/down), left-drag orbits yaw+pitch, mouse wheel
  and `+`/`-` zoom, Space toggles start/pause, toolbar
  Start/Pause/Stop/Reverse work, statusbar shows mode/t/dt/E/bodies/
  contacts/camera/fps, Contacts button + `C` toggle the golden
  contact-normal arrows and the frame protocol carries the exact
  analytic normal/point/impulse for a head-on impact, the dashed
  `#5d84a8` interior box wireframe draws when a BOX exists (the six
  wall slabs are never drawn as bodies), and torus (outer/inner
  equators + tube rings + 4 cross-sections), disk (rim + 2
  diameters) and cylinder (2 rims + 4 side lines) wireframes render
  quaternion-rotated so spin is visible; the permanent toolbar
  Reset button (`bt-reset`) re-initializes the playback (every value
  and the time bit-back to their initial values) and Start re-runs
  the simulation from the beginning; the labeled 'conserved
  quantities' readout (`hud`: E, P and L with components and
  magnitudes) updates live and reads identically before and after a
  dumbbell impact; entity labels show the registered user names
  (`dumbell0`) instead of `objN`; and dumbbells render as one rigid
  body — two shaded spheres at their rotated COM offsets joined by
  the rod's four silhouette lines.
