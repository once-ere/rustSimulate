#!/usr/bin/env bash
#
# Certify that this repository carries no licence-encumbered source.
#
# Run it against a FRESH PLAIN CLONE, which is the only thing that tells
# you what a user actually receives:
#
#     git clone https://github.com/once-ere/rustSimulate.git /tmp/cert
#     /tmp/cert/scripts/certify_clean.sh
#
# WHY THIS EXISTS
#
# On 2026-07-26, 749 files of the SolveIt 2002 C/C++ sources — including
# Numerical Recipes and GSL-derived code — were pushed to the then-public
# repository, because `.gitignore` said `./obsolete_or_historic` and a
# leading `./` matches nothing. Three separate ad-hoc checks failed to
# catch it: one was scoped to the wrong file set, one used `ls ... ||`
# so it could not fail loudly, and one set a flag inside `$(...)` where
# the assignment died with the subshell.
#
# So this script is written to fail LOUDLY and to be the single place
# that knows what "clean" means. Every check asserts on a count, prints
# what it found when it fails, and the exit status is real.
#
# See CLEANROOM_PROVENANCE.md §7.

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 2

fail=0
pass() { printf '  \033[32mPASS\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; fail=1; }

# Paths that must never be tracked. These are local-only reference trees.
FORBIDDEN_PATHS='^(obsolete_or_historic|SolveIt|SolveIt_2026_MFC)(/|$)'

# The SAME set, matched against `git rev-list --objects` output, whose
# lines are "<sha> <path>" — so the path is NOT at the start of the line.
# Anchoring this one with `^` is precisely the bug that made check 2
# vacuous the first time. Keep the two patterns separate and obviously
# different so the distinction cannot be lost in an edit.
FORBIDDEN_OBJECTS='[[:space:]](obsolete_or_historic|SolveIt|SolveIt_2026_MFC)(/|$)'

# Signatures of encumbered CODE. Deliberately matched only against
# SOURCE files: the documentation legitimately *discusses* Numerical
# Recipes by name (CLEANROOM_PROVENANCE.md, THIRD_PARTY.md), and a check
# that cannot tell "mentions" from "contains" is a check that cries wolf.
CODE_SIGNATURES='bessj0|bessj1|nrutil|gsl_linalg_solve|gsl_vector_|NR_END|FREE_ARG'
SOURCE_EXT='\.(rs|c|h|cpp|hpp|cc|f|f90|py)$'

echo "Certifying $(pwd)"
echo "HEAD: $(git log --oneline -1 2>/dev/null || echo '<not a git repo>')"
echo

# ---- 0. SELF-TEST: prove each detector can actually fire ---------------
#
# The first version of this script shipped a VACUOUS gate. Check 2 piped
# `git rev-list --objects`, whose lines are "<sha> <path>", into a
# pattern anchored with `^` — so it matched against the SHA and never
# fired. 676 forbidden objects were present and it reported PASS.
#
# A check that cannot fail is worse than no check, because it manufactures
# confidence. So every detector below is first run against a synthetic
# line that MUST match. If a detector does not fire on known-bad input,
# this script aborts rather than printing reassuring PASS lines.
selftest() { # <description> <text that must match> <pattern>
  if ! printf '%s\n' "$2" | grep -qE "$3"; then
    printf '  \033[31mABORT\033[0m  self-test failed: %s\n' "$1"
    printf '         the detector below cannot fire, so its PASS would be meaningless\n'
    exit 2
  fi
}
selftest "tracked-path detector"  "obsolete_or_historic/SolveIt/QM/QMEvolve.h" "$FORBIDDEN_PATHS"
selftest "history-object detector" \
         "9da888ae98f3ac8b3e997384da2654a298f5dcd3 obsolete_or_historic/SolveIt/QM/QMEvolve.h" \
         "$FORBIDDEN_OBJECTS"
selftest "code-signature detector" "float bessj0(float x)" "$CODE_SIGNATURES"

# ---- 1. no forbidden path tracked at the tip --------------------------
n=$(git ls-files | grep -cE "$FORBIDDEN_PATHS")
if [ "$n" -eq 0 ]; then
  pass "no reference-tree files tracked"
else
  bad "$n reference-tree files are TRACKED:"
  git ls-files | grep -E "$FORBIDDEN_PATHS" | head -10 | sed 's/^/        /'
fi

# ---- 2. no forbidden path anywhere in history -------------------------
# Removing files at the tip is not enough: old commits stay fetchable.
n=$(git rev-list --objects --all | grep -cE "$FORBIDDEN_OBJECTS" || true)
if [ "$n" -eq 0 ]; then
  pass "no reference-tree objects reachable in history"
else
  bad "$n reference-tree objects are still REACHABLE IN HISTORY"
  git rev-list --objects --all | grep -E "$FORBIDDEN_OBJECTS" \
    | head -5 | sed 's/^/        /'
  echo "        Reachable from these refs:"
  git for-each-ref --format='          %(refname:short)' | sed 's/$//' | head -10
  echo "        (removing files at the tip is NOT enough — the objects stay"
  echo "         fetchable by SHA. Rewrite history, or drop the ref holding them.)"
fi

# ---- 3. no encumbered code signatures in source files -----------------
hits=$(git ls-files | grep -E "$SOURCE_EXT" \
       | xargs -r grep -lniE "$CODE_SIGNATURES" 2>/dev/null || true)
if [ -z "$hits" ]; then
  pass "no encumbered-code signatures in tracked source"
else
  bad "encumbered-code signatures found in tracked SOURCE files:"
  echo "$hits" | sed 's/^/        /'
fi

# ---- 4. the ignore rules must actually match --------------------------
# The whole incident was an ignore pattern that looked right and matched
# nothing. Never assume — interrogate git.
for p in obsolete_or_historic SolveIt SolveIt_2026_MFC; do
  if git check-ignore -q "$p/probe.h"; then
    pass "gitignore matches $p/"
  else
    bad "gitignore does NOT match $p/ — a pattern like './$p' silently ignores nothing"
  fi
done

# ---- 4b. Stages_output must never be published ------------------------
# An internal record of every stage response. It is deliberately not part
# of the repository, and one mechanism is not evidence: the ignore rule is
# interrogated, AND the index is checked for tracked files. Either failing
# is a defect.
if git check-ignore -q "Stages_output/probe.md"; then
  pass "gitignore matches Stages_output/"
else
  bad "gitignore does NOT match Stages_output/ — the transcripts could be published"
fi
tracked_stages=$(git ls-files Stages_output/ | head -20)
if [ -z "$tracked_stages" ]; then
  pass "no Stages_output files tracked"
else
  bad "Stages_output files are TRACKED and would be pushed:"
  echo "$tracked_stages" | sed 's/^/        /'
fi

# ---- 5. the build and test gates --------------------------------------
if cargo build --workspace --release >/tmp/cert_build.$$ 2>&1; then
  w=$(grep -c '^warning' /tmp/cert_build.$$ || true)
  if [ "$w" -eq 0 ]; then pass "release build, 0 warnings"
  else bad "release build emitted $w warnings"; fi
else
  bad "release build FAILED"; tail -20 /tmp/cert_build.$$ | sed 's/^/        /'
fi
rm -f /tmp/cert_build.$$

# clippy across ALL targets, not just the library: examples, tests and
# benches are shipped code too, and lints only surface in the crate the
# run reaches before the first failure — so a partial pass here is
# genuinely uninformative.
if cargo clippy --workspace --all-targets >/tmp/cert_clippy.$$ 2>&1; then
  pass "clippy --workspace --all-targets, 0 errors"
else
  bad "clippy reported problems:"
  grep -E '^(error|warning)' /tmp/cert_clippy.$$ | head -10 | sed 's/^/        /'
fi
rm -f /tmp/cert_clippy.$$

if cargo test --workspace >/tmp/cert_test.$$ 2>&1; then
  p=$(grep -E 'test result' /tmp/cert_test.$$ | awk '{s+=$4} END {print s}')
  pass "all tests pass ($p assertions across the workspace)"
else
  bad "tests FAILED"
  grep -A4 'panicked' /tmp/cert_test.$$ | head -20 | sed 's/^/        /'
fi
rm -f /tmp/cert_test.$$

echo
if [ "$fail" -eq 0 ]; then
  printf '\033[32mCLEAN\033[0m — nothing encumbered is tracked or reachable.\n'
else
  printf '\033[31mDEFECTS REMAIN\033[0m — see the FAIL lines above.\n'
fi
exit "$fail"
