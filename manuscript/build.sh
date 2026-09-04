#!/usr/bin/env bash
# Build the arXiv PDF, and leave this directory as a self-contained upload.
#
# arXiv wants a flat-ish source tree with the figures alongside the .tex, so the figures are
# COPIED from results/figures/ rather than referenced out of the repo. That means the upload
# cannot go stale relative to a figure that was regenerated: rerun this script.
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p figures
# DERIVED FROM THE MANUSCRIPT, not maintained by hand. This was a literal list, and adding
# f16 to the text without adding it here would have shipped an upload referencing a figure the
# tree does not contain -- which LaTeX reports as a missing-file warning that is easy to miss in
# a long log. Now the list cannot drift from what the sections actually cite.
FIGS=$(grep -ho 'figures/f[0-9_a-z]*' sections/*.tex | sed 's|figures/||' | sort -u)
for f in $FIGS; do
  [ -f "../results/figures/$f.pdf" ] || { echo "missing ../results/figures/$f.pdf" >&2; exit 1; }
  cp "../results/figures/$f.pdf" "figures/$f.pdf"
done

# Supplementary Table S1 is cited in Data availability, so it ships WITH the manuscript and
# not only in results/. A submission whose supplementary file is a repository path is not a
# submission.
cp ../results/tables/supplementary_table_s1.csv supplementary_table_s1.csv

command -v pdflatex >/dev/null || { echo "pdflatex not found; install MacTeX or TeX Live"; exit 1; }
pdflatex -interaction=nonstopmode -halt-on-error paper.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error paper.tex >/dev/null
echo "wrote paper.pdf ($(wc -c < paper.pdf) bytes, $(pdfinfo paper.pdf 2>/dev/null | awk '/^Pages/{print $2}') pages)"

# Undefined references are silent in a nonstopmode build and fatal in a preprint.
if grep -qE "LaTeX Warning: (Citation|Reference).*undefined" paper.log; then
  echo "UNDEFINED references or citations:" >&2
  grep -E "LaTeX Warning: (Citation|Reference).*undefined" paper.log >&2
  exit 1
fi
echo "no undefined citations or references"
