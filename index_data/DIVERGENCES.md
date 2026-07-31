# DIVERGENCES — where the prose and the code disagree

*Phase-2 deliverable of the Index of Entities (see `prompt_01.md` §1 defect 2,
§8 phase 2). Extraction date: **2026-07-30**.*

**Policy.** Code is authoritative; prose is corroborative. Every entity in the
index takes its definition from source and is checked against the
documentation. Where the two disagree, the disagreement is recorded here
rather than silently resolved — the index will carry the code's behaviour and
a note pointing at this file.

Each finding states the claim, the code truth, and — where behaviour is at
stake — the verbatim probe that settled it. **Nothing below is asserted from
reading alone.**

---

## D1 — `grammar.md` §2.2 omits three keywords that the lexer defines

**Claim.** `grammar.md:63-80` enumerates the reserved keywords in four
groups (core, scene, collision, names/functions).

**Code.** `posim/src/lexer.rs:36-38` defines `Qm`, `Qm2`, `Qm3` in the
`Keyword` enum, and `lexer.rs:104-106` maps the spellings `"qm"`, `"qm2"`,
`"qm3"` to them. They are reserved exactly as `NEW` is.

**Effect.** 57 keyword spellings exist; §2.2 lists 54. The three quantum
family heads are missing from the list, though §§5.10–5.12 document the
families themselves at length.

**Index treatment.** `QM`, `QM2`, `QM3` get keyword entries, sourced to
`lexer.rs`, cross-referenced to §§5.10–5.12.

---

## D2 — `parser.rs`'s own EBNF comment advertised a `QM2 ISO` the parser rejects — **FIXED**

**Claim.** `posim/src/parser.rs:136`, inside the module-level EBNF:

```
//! qm2cmd   := …
//!           | "ISO" STRING expr [ "FRAMES" expr ] [ "LEVEL" expr ]
```

**Code.** `Parser::qm2_command()` has no `iso` arm. The sub-commands it
accepts are `status, grid, potential, packet, step, run, norm, energy,
centroid|position, prob|probability, absorb, drive, states, state, reset,
animate`. `qm3_command()` *does* have `"iso" | "isosurface"`.

**Verified 2026-07-30** — `posim --script`:

```
In[3]:= qm2 iso "x.html" 0.1
Err[3]: QM2: unknown subcommand `iso` (grid, potential, packet, states, state, step, run, norm, energy, centroid, prob, absorb, animate, status, reset)
```

**Who is right.** `HELP_TEXT` (`vm.rs`) and `grammar.md` §5.11 are both
correct — neither documents `QM2 ISO`. The defect is confined to the EBNF
comment in `parser.rs`, which is the one place a reader would trust most.

**Index treatment.** No `QM2 ISO` entry. The `QM3 ISO` entry notes that the
2-D family has no isosurface command, and why (a 2-D density is already
drawable as a heat map).

### Fixed, 2026-07-30

The phantom production is gone from `parser.rs`, and — more to the point —
a gate now makes it uncatchable-by-reading no longer uncatchable at all.
`every_qm2_subcommand_is_documented_in_lockstep` and its `QM3` twin check
their production in **both** directions:

- *forward*, that every declared subcommand appears in `HELP_TEXT`, the
  EBNF, `grammar.md` and `grammar.tex`;
- *converse*, that every word the EBNF **quotes** is either a declared
  subcommand or a declared argument word.

The converse direction is the one that earns its keep, and it is the one
the pre-existing `QM` gate never had. A forward-only check asks whether
the documents mention what the code does; it can never ask whether the
code does what the documents promise. Re-injecting the deleted line makes
the test fail with:

```
QM2 grammar lockstep is broken:
  the qm2cmd EBNF quotes `ISO`, which is neither a declared subcommand nor a
  declared argument word — either the parser is missing an arm or the comment
  is promising a command that does not exist
  the qm2cmd EBNF quotes `LEVEL`, which is neither a declared subcommand nor a
  declared argument word — …
```

`LEVEL` is the tell: it was only ever an argument of the isosurface
command, so it had been orphaned by the same mistake and nothing noticed.

**Extending the gate found four more gaps**, all of the same shape — the
compressed spelling `QM2 NORM | ENERGY | CENTROID` in `HELP_TEXT` and
`grammar.tex`, which names the family once and then lists bare words. A
reader parses that correctly; a string search does not, and neither does
anyone grepping for `QM2 CENTROID`. Now spelled out in full in both files,
matching how `grammar.md` already wrote them.

The subcommand lists are also no longer duplicated: `parser.rs` builds its
`QM2`/`QM3: unknown subcommand` error messages from `QM2_SUBCOMMANDS` and
`QM3_SUBCOMMANDS`, exactly as it already did for `QM_SUBCOMMANDS`. A
duplicated list is a list that drifts.

---

## D3 — six path aliases work but appear in no documentation

`posim/src/vm.rs` accepts a second spelling for several read paths. None of
these appear in `grammar.md` §5.2/§5.7, `physical_object_simulator.md` or
`collision_detection.md` as paths. (`approach` occurs in those files only as
an English word — `grammar.md:1228`, `:1401`, `:1420`; `impulse_n` occurs
once, at `collision_detection.md:343`, as a *Rust struct field*, not a path.)

| alias | canonical | source |
|---|---|---|
| `contactK.time` | `contactK.t` | `vm.rs:2392` |
| `contactK.approach` | `contactK.rel_vel_n` | `vm.rs:2396` |
| `contactK.impulse_n` | `contactK.impulse` | `vm.rs:2397` |
| `QM2 position` / `QM3 position` | `… centroid` | `parser.rs` qm2/qm3 `centroid \| position` |
| `QM2 probability` / `QM3 probability` | `… prob` | `parser.rs` qm2/qm3 `prob \| probability` |
| `QM3 isosurface` | `QM3 iso` | `parser.rs` qm3 `"iso" \| "isosurface"` |

**Verified 2026-07-30** — head-on exchange of two unit spheres, contact at
TOI 1.5:

```
In[5]:= get contact0.t
Out[5]= 1.5
In[6]:= get contact0.time
Out[6]= 1.5
In[7]:= get contact0.rel_vel_n
Out[7]= -2
In[8]:= get contact0.approach
Out[8]= -2
In[9]:= get contact0.impulse
Out[9]= 2
In[10]:= get contact0.impulse_n
Out[10]= 2
```

and, on a 16×16 2-D grid with a Gaussian packet:

```
In[8]:= qm2 position
Out[8]= [0.0000000000000000878205635495564, 0.00000000000000010353961926829087]
In[9]:= qm2 probability -1 1, -1 1
Out[9]= 0.432954441903126
```

**Index treatment.** Each is recorded in its canonical entry's `aliases`
array and resolves through search, flagged *undocumented alias* so the entry
is honest about where it came from.

---

## D4 — `grammar.md` §3's EBNF predates the comparison operators

**Claim.** `grammar.md:173`:

```ebnf
expr     := term { ("+" | "-") term } ;
```

**Code.** `posim/src/parser.rs:68,79`:

```
expr     := sum { ("<" | "<=" | ">" | ">=" | "==" | "!=") sum } ;
sum      := term { ("+" | "-") term } ;
```

**Effect.** The grammar block in §3 has no comparison level at all, so it
cannot derive `(x > 0) * (x < 1)` — the indicator idiom that §4.2, §5.10 and
Examples 17 and 18 are built on. §4.2 documents comparisons correctly,
including their lowest precedence; only the formal grammar block is stale.

**Index treatment.** The index quotes `parser.rs`'s EBNF as the syntax of
record for every expression-level entry.

---

## D5 — `special_functions.md` calls `wigner` "planned"; it shipped

**Claim.** `special_functions.md:33`:

| `wigner` — 3j, 6j, Clebsch–Gordan | native | **planned** |

and `special_functions.md:568`: *"A `wigner` module (3j, 6j,
Clebsch–Gordan) is planned and will be added in the same shape."*

**Code.** `special_functions/src/wigner.rs` is 871 lines and defines
`wigner_3j` (:133), `clebsch_gordan` (:239), `wigner_6j` (:275) and
`wigner_9j` (:398) — the last of which the doc's own "planned" list does not
even mention. All four are registered notebook builtins in
`posim/src/special.rs`, and `grammar.md` §4.1 documents them as available.

**Index treatment.** All four get full builtin entries. The status page notes
that `special_functions.md`'s status table is behind its own §4.1 sibling.

---

## D6 — ten `special_functions` modules exist but are absent from `special_functions.md`

The document's "Status at a glance" table and section list cover
`sph_bessel`, `legendre`, `orthopoly`, `eigen`, `quadrature`, `complex`,
`tridiag`, `bessel`, `wigner`. The crate also contains:

`airy_complex` · `airy_uniform` · `bessel_cnu` · `bessel_cnu_large` ·
`bessel_complex` · `bessel_scaled` · `debye` · `gamma_complex` · `hankel` ·
`lanczos`

That is 10 of 20 modules undocumented in the file nominally dedicated to
them — 8,006 of the crate's lines. Their *notebook-facing* surface (`airy_z`,
`gamma_z`, `hankel_h1_z`, `bessel_j_nu`, `bessel_k_scaled`, …) **is**
thoroughly documented, in `grammar.md` §4.1 instead. The gap is at module
level, not at API level.

**Index treatment.** Module entries are sourced from the crate; the accuracy
laws and route-selection discussion come from `grammar.md` §4.1, which is the
live document for this material.

---

## D7 — `special_functions.md`'s own header states the policy D5 breaks

Recorded not as a separate defect but as the standard D5/D6 should be read
against — `special_functions.md:8-10`:

> "This document grows one module at a time, alongside the code. Sections
> marked **planned** are not implemented yet and say so rather than
> pretending otherwise."

The intent is right; D5 is the case where the marker was not updated after
the module landed.

---

## Non-divergences checked and cleared

Recorded so a later pass does not re-investigate them.

- **`system.collide` is read-only.** `grammar.md` §5.2 marks it `R`. A write
  arm exists at `vm.rs:2481` but returns `use COLLIDE ON / COLLIDE OFF to
  switch collision detection`. The doc is correct; the arm is a guiding
  refusal, not a setter.
- **`collisions` is deliberately not a keyword**, so `system.collisions`
  keeps its spelling (`grammar.md:85`). Confirmed absent from
  `lexer.rs::Keyword::from_ident`.
- **`DEF` is not a keyword.** Confirmed: no `"def"` arm in `from_ident`. It is
  a line form recognised before the grammar, exactly as documented.
- **All 61 registered special-function names** in `posim/src/special.rs` are
  present in `grammar.md` §4.1's tables. No omissions either way.
- **All seven shapes** (`POINT SPHERE CUBOID TORUS DISK CYLINDER DUMBBELL`)
  and all seven documented aliases (`DELETE CUBE DISC DESTROY SETTIMESTEP
  FUNCTIONS DUMBELL`) match `lexer.rs` exactly.
- **17 `SCENE` sub-commands** in `parser.rs` match `grammar.md` §5.6's table
  and `HELP_TEXT` one for one, `RESET` included.
- **23 `QM` sub-commands** match: `posim/src/qm.rs:102` exposes
  `QM_SUBCOMMANDS`, and the test `every_qm_subcommand_is_documented_in_lockstep`
  already checks each word against the parser, `HELP_TEXT`, the EBNF comment
  and both grammar documents. That gate is why the `QM` family has no
  D1-class defect — and its absence for `QM2`/`QM3` is why D2 survived.

---

## Status

| # | finding | state |
|---|---|---|
| D1 | `grammar.md` §2.2 omits `QM`/`QM2`/`QM3` | open — documentation |
| D2 | phantom `QM2 ISO` in the `parser.rs` EBNF | **fixed 2026-07-30**, with a gate |
| D3 | six undocumented path aliases | open — documentation |
| D4 | `grammar.md` §3's EBNF has no comparison level | open — documentation |
| D5 | `special_functions.md` calls `wigner` "planned" | open — documentation |
| D6 | ten `special_functions` modules undocumented | open — documentation |

D2 was the only finding that was a defect in *code-adjacent* material rather
than in prose, and the only one a reader could act on and be wrong. The five
that remain are documentation drift: the index carries the code's behaviour
and points here, so no reader is misled by them in the meantime.
