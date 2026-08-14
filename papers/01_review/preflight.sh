#!/usr/bin/env bash
# Pre-submission checks for Paper 01.
#
# Every item here failed at least once during preparation, which is why
# it is a script rather than a checklist: the manuscript, the pinned
# commit and the LaTeX build drift apart silently, and the build keeps
# succeeding while they do.
#
# Usage:  bash papers/01_review/preflight.sh
# Exit:   0 if every check passes, 1 otherwise.

set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

MS=papers/01_review/manuscript.md
TEX=papers/01_review/latex
FIG=papers/01_review/figures/out
fail=0
ok()   { printf '  \033[32mOK\033[0m   %s\n' "$1"; }
bad()  { printf '  \033[31mFAIL\033[0m %s\n' "$1"; fail=1; }
warn() { printf '  \033[33mWARN\033[0m %s\n' "$1"; }

echo "== pinned commit =="
PIN=$(grep -oP '(?<=\(commit `)[a-f0-9]+' "$MS" | head -1)
[ -n "$PIN" ] && ok "manuscript pins $PIN" || bad "no pinned commit found"

if [ -n "$PIN" ]; then
  # Every example the manuscript names must exist at the pinned commit.
  miss=0
  for e in $(grep -oP '(?<=--example )[a-z0-9_]+' "$MS" | sort -u); do
    git cat-file -e "$PIN:solver-2d/examples/$e.rs" 2>/dev/null \
      || { bad "example '$e' is cited but absent from $PIN"; miss=1; }
  done
  [ $miss -eq 0 ] && ok "all cited examples exist at the pin"

  # Flags printed in the paper must exist at the pinned commit too — a
  # command that parses but ignores an unknown flag silently produces
  # the wrong run.
  fmiss=0
  # Extract every flag on each printed cargo-run line, not just the
  # first after ' -- ': the earlier version of this check missed --cfl
  # precisely because it stopped at the first match on the line.
  for line in $(grep -n 'cargo run' "$MS" | cut -d: -f1); do
    cmd=$(sed -n "${line}p" "$MS")
    ex=$(echo "$cmd" | grep -oP '(?<=--example )[a-z0-9_]+')
    [ -z "$ex" ] && continue
    src=$(git show "$PIN:solver-2d/examples/$ex.rs" 2>/dev/null) || continue
    for f in $(echo "${cmd#*-- }" | grep -oP '(?<![a-z-])--[a-z][a-z-]*'); do
      case "$f" in --release|--example) continue;; esac
      echo "$src" | grep -q -- "\"$f\"" \
        || { bad "flag '$f' printed for '$ex' is absent from $PIN"; fmiss=1; }
    done
  done
  [ $fmiss -eq 0 ] && ok "all printed flags exist at the pin"

  # Code must not have moved since the pin.
  n=$(git log --oneline "$PIN"..HEAD -- solver-2d/ autograd/ examples/ 2>/dev/null | wc -l)
  [ "$n" -eq 0 ] && ok "no code commits after the pin" \
                 || bad "$n code commit(s) after the pin — re-pin needed"
fi

echo "== journal limits =="
AW=$(python3 -c "
s=open('$MS').read(); a=s[s.index('# Abstract'):s.index('# Key Points')]
print(len(a.replace('# Abstract','').split()))")
[ "$AW" -le 150 ] && ok "abstract $AW words (limit 150)" || bad "abstract $AW words > 150"

python3 - "$MS" << 'PY'
import sys
s=open(sys.argv[1]).read()
h=[l[2:].strip() for l in s[s.index('# Highlights'):s.index('# Abstract')].split('\n') if l.startswith('- ')]
n=len(h); mx=max(len(x) for x in h)
print(f"  {'OK  ' if 3<=n<=5 else 'FAIL'} {n} highlights (3-5)")
print(f"  {'OK  ' if mx<=85 else 'FAIL'} longest highlight {mx} chars (limit 85)")
PY

echo "== required artifacts =="
for f in "$TEX/paper.pdf" "$TEX/highlights.txt" "$FIG/graphical_abstract.pdf"; do
  [ -f "$f" ] && ok "$(basename "$f") present" || bad "$(basename "$f") missing"
done
grep -qi "data.*availab\|data statement" "$MS" \
  && ok "data availability statement present" \
  || bad "no data availability statement (EMS applies Option C)"

echo "== build integrity =="
if [ -f "$TEX/paper.log" ]; then
  u=$(grep -c "Citation.*undefined" "$TEX/paper.log")
  [ "$u" -eq 0 ] && ok "0 undefined citations" || bad "$u undefined citations"
fi
if [ -f "$TEX/paper.pdf" ]; then
  im=$(pdfimages -list "$TEX/paper.pdf" 2>/dev/null | tail -n +3 | wc -l)
  [ "$im" -gt 0 ] && ok "$im images embedded" || bad "no images embedded"
fi

echo "== staleness =="
# The LaTeX tree carries hand-written sections; if it predates the
# manuscript it is shipping older claims than the ones we corrected.
if [ "$MS" -nt "$TEX/paper.pdf" ]; then
  bad "manuscript is newer than paper.pdf — rebuild"
else
  ok "paper.pdf at least as new as the manuscript"
fi
for term in "one-day peak" "Aluvión" "no ergonomic path" "Major Revision"; do
  if pdftotext "$TEX/paper.pdf" - 2>/dev/null | grep -qi -- "$term"; then
    bad "retracted phrasing still in the PDF: '$term'"
  fi
done
ok "no retracted phrasing found in the PDF"

echo
[ $fail -eq 0 ] && echo "PREFLIGHT PASSED" || echo "PREFLIGHT FAILED"
exit $fail
