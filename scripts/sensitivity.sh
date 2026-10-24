#!/usr/bin/env bash
# Is "BUB1B ranks 1" a property of the data, or of the thresholds we chose?
#
# rank_biallelic_genes.py states plainly that three filters do nearly all the
# work. That is a disclosed limit and also an open question: a result
# that only survives at ONE setting of three hand-chosen knobs is a result about
# the knobs. Sweep them and report the rank at every combination.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
AF_TSV="${MVA_GNOMAD_TSV:-data/private/work/gnomad_af_genomewide.tsv}"
# Resolve the interpreter the way the Makefile does, so this works from a bare
# shell with nothing activated. Override with PY=... for a different one.
PY="${PY:-$(test -x .venv/bin/python && echo .venv/bin/python || echo python3)}"
[ -s "$AF_TSV" ] || { echo "missing $AF_TSV - run 'make gnomad-af' first" >&2; exit 1; }
printf "%-6s %-6s %-6s %-8s  %-10s %s\n" MQ minDP maxDP maxAF "BUB1B rank" "of N genes"
for MQ in 40 50 55 60; do
 for MINDP in 10 15 20; do
  for MAXDP in 80 100 150; do
   for AF in 0.01 0.001 0.0001; do
     out=$("$PY" scripts/rank_biallelic_genes.py \
             --gnomad-tsv "$AF_TSV" --min-mq "$MQ" --min-dp "$MINDP" \
             --max-dp "$MAXDP" --max-af "$AF" --top 1 2>/dev/null \
           | grep '^>>> BUB1B')
     if [ -z "$out" ]; then
       printf "%-6s %-6s %-6s %-8s  %-10s %s\n" "$MQ" "$MINDP" "$MAXDP" "$AF" "ABSENT" "-"
     else
       r=$(echo "$out" | sed -E 's/.*rank ([0-9]+) of ([0-9]+).*/\1/')
       n=$(echo "$out" | sed -E 's/.*rank ([0-9]+) of ([0-9]+).*/\2/')
       printf "%-6s %-6s %-6s %-8s  %-10s %s\n" "$MQ" "$MINDP" "$MAXDP" "$AF" "$r" "$n"
     fi
   done
  done
 done
done
