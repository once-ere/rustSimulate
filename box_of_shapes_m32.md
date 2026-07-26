# The Box of Shapes with a Heavy Point (mass 32)

*A complete, reproducible walkthrough: what the scenario is, every CLI
command and window action needed to run it, the numbers it must
produce, and the one setting that will silently ruin it if you get it
wrong. Written for a reader who has never run this simulator.*

Scenario file: [`dynamic_notebooks/box_of_shapes_m32.posim`](dynamic_notebooks/box_of_shapes_m32.posim)
Parent scenario: [`dynamic_notebooks/box_of_shapes.posim`](dynamic_notebooks/box_of_shapes.posim)
(identical except the moving point has **mass 32** instead of mass 1)

---

## 1. What this scenario is

Six bodies — one of **every** shape the simulator supports — sit inside
a rigid cube of side 4:

| handle | body | mass | notes |
|---|---|---|---|
| `obj6` | torus | 1 | inner radius 1, outer 2, axis tilted to (1,1,1)/√3 |
| `obj7` | **point** | **32** | the only thing moving: **v = (100, 200, 100)** |
| `obj8` | sphere | 2 | radius ½ |
| `obj9` | disk | 2/3 | radius 1, ideal zero thickness |
| `obj10` | cuboid | 5/3 | a ½-cube |
| `obj11` | cylinder | 2 | radius ½, height 3/2, tilted |
| `obj0`–`obj5` | the six box walls | — | `inverse_mass = 0` (infinitely massive) |

Three settings make this a clean physics experiment:

- **`system.g_constant = 0`** — no gravity at all, so nothing accelerates
  between impacts and every change in motion comes from a collision.
- **`restitution = 1` on every body** — this is the default, so the
  notebook never sets it; collisions are perfectly elastic.
- **walls with `inverse_mass = 0`** — infinitely massive, so they bounce
  things without ever moving.

Two consequences you can check while it runs:

- **Energy is conserved exactly.** `E = ½·m·|v|² = ½·32·60000 = 960000`,
  forever, through every impact.
- **Momentum is NOT conserved.** `P₀ = 32·(100, 200, 100) = (3200, 6400,
  3200)`, but every wall bounce hands momentum to an infinitely massive
  object, so `|P|` wanders downward. That is correct, not a bug.

### What changing the mass to 32 does

The point now outweighs every other body (the heaviest is 2). So it
barely deflects when it hits one — instead it *kicks* it. Those light
bodies then ricochet around the box on their own, and the scene becomes
far busier: **9293 collisions** in `0 ≤ t ≤ 1`, versus **1622** for the
mass-1 original — a 5.7× increase.

---

## 2. Reproducing it — every command

### Step 1 — get the code and build

```bash
git clone https://github.com/once-ere/rustSimulate.git
cd rustSimulate
cargo build --release
```

No submodules, no network access during the build, nothing from
crates.io. The release build matters: this scenario runs tens of
thousands of collision events.

### Step 2 — launch the notebook with its scene window

```bash
cargo run -p posim --release -- --notebook dynamic_notebooks/box_of_shapes_m32.posim
```

This executes the file and then leaves you at a live `In[n]:=` prompt.
It prints the setup, the baselines, and the window address:

```
Out[11]= 0                          <- obj0.inverse_mass: the wall is static
Out[12]= 1                          <- obj7.restitution: perfectly elastic
Out[13]= 960000                     <- E0, exactly
Out[14]= [3200, 6400, 3200]         <- P0, exactly
Out[15]= scene window created: http://127.0.0.1:41737/
Out[16]= scene time step dt = 0.0002
Out[17]= camera distance = 6
```

The port is chosen by the OS and will differ each run. Your browser
should open automatically; if it does not, paste that address into it.
(Set `POSIM_NO_BROWSER=1` to suppress the automatic launch — useful on
a headless machine.)

### Step 3 — run it

Either press **Start** in the window, or type at the prompt:

```
scene start
```

Watch the **conserved quantities** panel at the top left. `E` must stay
pinned at `9.60000e+5` while `P` and `L` wander.

### Step 4 — inspect while it runs

The notebook stays responsive because the window animates its own
synchronized *copy* — your notebook state never moves during playback.

```
scene status            # URL, connected windows, mode, t, dt, steps, history
scene events            # drain the window's asynchronous messages
energy                  # notebook-state energy (still 960000 — the copy moved, not this)
get obj7.restitution    # 1
get obj0.inverse_mass   # 0
```

To stop and leave:

```
scene close
%quit
```

---

## 3. The scene window controls

Every control, with the exact label and its tooltip:

| control | tooltip / action |
|---|---|
| **▶ Start** | Start forward evolution (Space) |
| **❚❚ Pause** | Pause evolution (Space) — history is kept |
| **■ Stop** | Stop evolution and clear history |
| **◀ Reverse** | Play backward through recorded history |
| **↺ Reset** | Re-initialize: every value and the time return to their initial values; Start then replays from the beginning |
| **⇤** | One step backward |
| **⇥** | One step forward |
| **dt** box + **Set** | Time step for the playback loop / apply it |
| **🔍+ / 🔍−** | Zoom in / out (also `+` `-` or the mouse wheel) |
| **⌂ View** | Reset the view (`R`) |
| **Grid** | Toggle the ground grid (`G`) |
| **Trails** | Toggle motion trails (`T`) |
| **Labels** | Toggle object labels (`L`) |
| **Contacts** | Toggle contact-normal arrows at collision points (`C`) — the golden arrows |
| **?** | Show the controls cheat-sheet (`H` or `?`) |

Mouse and keyboard inside the canvas: **arrow keys** translate the view,
**left-drag** orbits, **wheel** zooms, **shift-drag** or **right-drag**
pans, **Space** starts/pauses.

The status bar along the bottom reads: connection light, `mode`, `t`,
`dt`, `E`, `bodies`, `hidden`, `contacts`, `history` depth, the camera
readout, and frames per second.

---

## 4. Verified results

### Live scene window (playback `dt = 0.0002`)

Captured from a real browser window attached to the running simulation:

```
conserved quantities
E  = 9.59999981e+5
P  = [1639.58853, -1568.09500, -2446.43369]  |.| = 3336.49671
L  = [3872.52136, 6343.70255, -680.81051]    |.| = 7463.40986

connected to posim   mode running   t = 0.2676   dt = 0.0002
E = 9.60000e+5   bodies 12   hidden 0   contacts 1   history 1338
cam yaw -60° pitch 55° dist 6.00   116 fps
```

`E` holds at 960000 to a relative error of about **2×10⁻⁸**, while `|P|`
has fallen from 7838.4 to 3336.5 — the walls absorbing momentum, exactly
as intended.

### Headless integration over `0 ≤ t ≤ 1`

Same scenario without the window (the same solver and the same collision
handling — only faster, because nothing is throttled to the frame rate):

```
In[15]:= energy
Out[15]= 960000
In[16]:= momentum
Out[16]= [3200, 6400, 3200]
In[17]:= run 1 steps 1000
Out[17]= t = 1 (119334 solver steps, 1000 snapshots, |dE/E| = 6.940e-8,
         9293 collision(s) — CONTACTS lists them)
In[18]:= energy
Out[18]= 959999.9333768479
In[19]:= momentum
Out[19]= [-242.27746642638098, -464.9012023666894, -342.43832783294954]
In[20]:= get system.collisions
Out[20]= 9293
```

**Energy conserved to 6.9×10⁻⁸ through 9293 elastic collisions.**

---

## 5. The setting that will silently ruin this

> **Ask for output more often than collisions happen.**

The point crosses the 4-unit box in about **0.016 time units**. If a
single solver output interval spans several collisions, they are not
resolved correctly and **energy collapses** — with no error message, no
warning, and a perfectly smooth-looking animation.

Measured on this exact scenario, `0 ≤ t ≤ 1`, changing *only* how many
snapshots were requested:

| command | interval | collisions | final E | \|dE/E\| | |
|---|---|---|---|---|---|
| `run 1 steps 1000` | 0.001 | 9293 | 959999.93 | 6.9×10⁻⁸ | correct |
| `run 1 steps 8` | 0.125 | 898 | **741962.41** | 2.3×10⁻¹ | **23% of the energy gone** |

The same effect measured on the mass-1 parent scenario, which makes the
threshold clear:

| snapshots over `0 ≤ t ≤ 1` | interval | \|dE/E\| | verdict |
|---|---|---|---|
| 8 | 0.125 | 9.6×10⁻¹ | collapses |
| 25 | 0.040 | 3.3×10⁻¹ | collapses |
| 100 | 0.010 | 1.7×10⁻⁸ | correct |
| 400 | 0.0025 | 1.0×10⁻⁸ | correct |
| 1000 | 0.0010 | 1.8×10⁻⁸ | correct |
| 4000 | 0.00025 | 2.2×10⁻⁸ | correct |

Conservation holds while the interval stays **below the mean free time
between collisions (~0.016)** and fails above it.

The same rule governs the window's playback step. The notebook ships
`scene set_time_step 0.0002` for this reason. Raising it to `0.1` — 500×
larger, about 24.5 units of travel per step inside a 4-unit box — runs
happily at 37 fps and destroys the physics:

| playback dt | E after the run |
|---|---|
| 0.0002 (shipped) | 29999.99986915 — correct |
| 0.1 | **78.66** — 99.7% of the energy gone |

*(that pair was measured on the mass-1 parent scenario)*

**Practical rule:** keep the playback `dt` at 0.0002, and when
integrating headlessly ask for at least ~100 snapshots per unit of
time. If you need a long span, raise the snapshot count with it —
`run 128 steps 12800`, not `run 128 steps 8`.

### Known limitation

This is a genuine defect, not merely a tuning knob: **the trajectory
should not depend on how often you ask for output.** The simulator
currently accepts a step size or snapshot cadence that makes collision
resolution unsound, reports `mode = running` as though all is well, and
gives a wrong answer silently. A guard — warning when `|v|·dt` exceeds a
fraction of the smallest body dimension, or when energy drifts in a
system where every restitution is 1 — would turn a silent wrong answer
into a visible one.

---

## 6. Why `0 ≤ t ≤ 128` is not an interactive request

At the correct playback step, `t = 128` needs `128 / 0.0002 = 640 000`
frames. The window renders 116–206 frames per second, so watching the
whole span takes **roughly 50–90 minutes** of wall clock. The animation
is meant for inspecting seconds of simulated time, not minutes.

For the long span, integrate headlessly instead — it is the same solver
and the same collision handling, just not throttled to a frame rate:

```
run 128 steps 12800      # interval 0.01: below the ~0.016 collision spacing
energy
get system.collisions
```

---

## 7. One-command summary

```bash
git clone https://github.com/once-ere/rustSimulate.git
cd rustSimulate
cargo build --release
cargo run -p posim --release -- --notebook dynamic_notebooks/box_of_shapes_m32.posim
# then press Start in the browser window, or type:  scene start
```

Expected at launch: `E0 = 960000`, `P0 = [3200, 6400, 3200]`,
`obj0.inverse_mass = 0`, `obj7.restitution = 1`, 12 entities.
Expected while running: `E` pinned at `9.60000e+5`, `|P|` decaying,
golden contact arrows firing as the heavy point scatters the light
bodies.
