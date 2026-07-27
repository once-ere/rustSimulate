# StageNbks — one standalone notebook per stage

Eleven posim notebooks, one for each stage considered. **Each is complete
in itself**: none refers you to another notebook, and every one carries its
own run instructions, its own browser-GUI procedure, its own explanation of
whether the stage can be animated, and the commands that recreate the
scenario.

| notebook | stage | animates? |
|---|---|---|
| `Stage_24.posim` | The sliver at 4 ≤ \|ν\| ≤ 8 | no — pure numerics, reason given inside |
| `Stage_2A.posim` | A faithful port of EVOLVE_NASH, with a Strang variant | **yes** — HTML wavefunction film |
| `Stage_2B.posim` | T(E) by transfer matrix, and a designed absorber | partly — reason given inside |
| `Stage_2C.posim` | The output-granularity defect, repaired | **yes** — full 3-D scene window |
| `Stage_2D.posim` | DLMF 10.20 past the turning point | no — pure numerics, reason given inside |
| `Stage_2E.posim` | Mechanising the staleness rule | no — build tooling, reason given inside |
| `Stage_2F.posim` | Mutation-probing the numerical core | no — build tooling, reason given inside |
| `Stage_2G.posim` | Closing the mutation survivors | partly — corner impact animates |
| `Stage_2H.posim` | Resolving the last survivor | **yes** — Newton's cradle and settling ball |
| `Stage_2I.posim` | The ridge is a real defect | no — pure numerics, reason given inside |
| `Stage_2J.posim` | A measured guard for the reflection route | no — pure numerics, reason given inside |

## Running one

From the repository root:

```bash
cargo run -p posim --release -- --notebook StageNbks/Stage_2C.posim
```

`--script` instead of `--notebook` prints the same transcript with a
different prompt layout. For an interactive session to paste into:

```bash
cargo run -p posim --release
```

Every notebook has been executed end to end; all eleven run with zero
errors.

## Running all of them

```bash
cargo build -p posim --release
for f in StageNbks/*.posim; do
  echo "=== $f ==="
  ./target/release/posim --notebook "$f"
done
```
