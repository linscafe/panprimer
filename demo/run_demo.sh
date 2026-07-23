#!/usr/bin/env bash
# MVP demo: two-stage universal-primer design over 3 diverse HPRC R2 haplotypes
# (AFR/EUR/EAS) — sized to fit in 15 GB RAM (one genome index in RAM at a time). Input and
# output live under demo/.
#
# Usage:  bash demo/run_demo.sh [TARGET] [NAME]
#   TARGET  CHM13 region 'chr:start-end' (default GAPDH, an easy housekeeping locus)
#   NAME    output subdir under demo/results (default gapdh)
set -euo pipefail
cd "$(dirname "$0")/.."
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pangenome-primer

CHM13=hprc-r2/references/chm13v2.0.fa
TARGET="${1:-chr12:6544868-6548730}"   # GAPDH
NAME="${2:-gapdh}"

# one-time: build the whole-genome projection PAF cache for the demo haplotypes (item 1)
bash scripts/prepare_haplotypes.sh demo/samples.tsv "$CHM13"

# two-stage run (stage A coverage on all candidates; stage B specificity on top-K)
pangenome-primer run \
  --target "$TARGET" --chm13 "$CHM13" \
  --samples demo/samples.tsv \
  --outdir "demo/results/$NAME" \
  --mode thermo --top-k 5

echo "report: demo/results/$NAME/report.html"
