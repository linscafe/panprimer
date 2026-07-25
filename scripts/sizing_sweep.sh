#!/usr/bin/env bash
# Measure the scanner's thread-scaling on THIS host and print the settings to use.
#
# The numbers in docs/sizing.md and ISSUE-002 were measured on a 16-core WSL2 laptop, most
# likely 8 physical cores plus hyperthreading. A server with different NUMA topology and a
# different memory system will have a different curve -- the optimum could reasonably be 8
# rather than 4. Do not copy those numbers onto a new machine; spend the two minutes.
#
#   scripts/sizing_sweep.sh hprc-r2/assemblies/HG00097_hap1.fa.gz
#   scripts/sizing_sweep.sh <hap.fa.gz> "1 2 4 8 16 32"     # custom thread list
#
# Part 1 answers "how well does ONE scan scale?" -- that is the number for a single-haplotype
# run. Part 2 answers "what shape is fastest when many haplotypes run at once?" -- that is the
# number for a real 464-haplotype job, and it is usually NOT the fastest single-scan setting.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null && conda activate pangenome-primer 2>/dev/null

GZ="${1:-}"
if [ -z "$GZ" ] || [ ! -s "$GZ" ]; then
  echo "usage: $0 <haplotype.fa.gz> [\"thread list\"]" >&2
  echo "  e.g. $0 hprc-r2/assemblies/HG00097_hap1.fa.gz" >&2
  exit 1
fi

CORES=$(nproc)
RAM_GB=$(awk '/MemTotal/ {printf "%.0f", $2/1048576}' /proc/meminfo)
# Default sweep: powers of two up to the core count, so the list adapts to the host instead
# of hard-coding the laptop's 1..16.
default_list=""
n=1; while [ "$n" -le "$CORES" ]; do default_list="$default_list $n"; n=$((n*2)); done
[ "$(echo "$default_list" | wc -w)" -gt 0 ] || default_list="1 2 4"
THREADS="${2:-$default_list}"

echo "host: ${CORES} logical cores, ${RAM_GB} GB RAM"
echo "file: $GZ ($(du -h "$GZ" | cut -f1))"
echo

PRIMERS='["ACGTACGTACGTACGTACGT","TTGCACAGTCCAGATTGCAA","GGATCCTTAAGGCCTTAAGC","TTGCACAGTCCAGATTGCAT"]'
scan_once() {  # $1 = threads
  PGP_SCAN_THREADS="$1" /usr/bin/time -f "%e %P %M" python -c "
import json
from pangenome_primer import rust_backend
rust_backend.find_binding_sites_batch(json.loads('$PRIMERS'), '$GZ', 'H#hap1', 3)
" 2>&1 | tail -1
}

echo "=== Part 1: one scan, varying threads ==================================="
printf "%8s %9s %9s %9s %11s %10s\n" threads wall speedup "core-s" efficiency "peak RSS"
base=""
best_t=1; best_wall=""
for t in $THREADS; do
  read -r wall cpu rss <<<"$(scan_once "$t")"
  cpu_n=${cpu%\%}
  cores_s=$(awk -v w="$wall" -v c="$cpu_n" 'BEGIN{printf "%.1f", w*c/100}')
  [ -z "$base" ] && base="$wall"
  sp=$(awk -v b="$base" -v w="$wall" 'BEGIN{printf "%.2f", b/w}')
  eff=$(awk -v s="$sp" -v t="$t" 'BEGIN{printf "%.0f%%", 100*s/t}')
  printf "%8s %8ss %8sx %9s %11s %9sMB\n" "$t" "$wall" "$sp" "$cores_s" "$eff" "$((rss/1024))"
  if [ -z "$best_wall" ] || awk -v a="$wall" -v b="$best_wall" 'BEGIN{exit !(a<b)}'; then
    best_wall="$wall"; best_t="$t"
  fi
done
echo
echo "fastest single scan: ${best_t} threads (${best_wall}s)  <- use for a ONE-haplotype run"

echo
echo "=== Part 2: many haplotypes at once, threads x concurrency ~= cores ====="
# Constant total work per row, so the wall times compare directly. 8 stands in for "a batch
# of haplotypes"; raise it with SWEEP_N for a longer, steadier measurement.
TOTAL="${SWEEP_N:-8}"
echo "(each row runs ${TOTAL} scans total, so the wall times are comparable)"
printf "%9s %12s %12s\n" threads concurrent "wall (total)"
for t in 2 4 8; do
  [ "$t" -gt "$CORES" ] && continue
  conc=$((CORES / t)); [ "$conc" -lt 1 ] && conc=1
  # Never claim more concurrency than there is work: with a small SWEEP_N the row would
  # otherwise be labelled 8-way while actually running 4, and read as a real comparison.
  [ "$conc" -gt "$TOTAL" ] && conc=$TOTAL
  t0=$(date +%s.%N)
  done_n=0
  while [ "$done_n" -lt "$TOTAL" ]; do
    batch=$conc; [ $((done_n + batch)) -gt "$TOTAL" ] && batch=$((TOTAL - done_n))
    for _ in $(seq 1 "$batch"); do
      PGP_SCAN_THREADS="$t" python -c "
import json
from pangenome_primer import rust_backend
rust_backend.find_binding_sites_batch(json.loads('$PRIMERS'), '$GZ', 'H#hap1', 3)
" >/dev/null 2>&1 &
    done
    wait
    done_n=$((done_n + batch))
  done
  t1=$(date +%s.%N)
  printf "%9s %12s %11.1fs\n" "$t" "$conc" "$(awk -v a="$t0" -v b="$t1" 'BEGIN{print b-a}')"
done

echo
echo "=== Apply the winning row from Part 2 ==================================="
echo "  Nextflow : nextflow run main.nf -profile local --scan_threads <threads>"
echo "             (declaring cpus=<threads> is what bounds concurrency; see nextflow.config)"
echo "  CLI      : set  search.threads: <threads>  in config/defaults.yaml"
echo "  ad hoc   : export PGP_SCAN_THREADS=<threads>"
echo
echo "Anchor-grid builds are RAM-bound, not CPU-bound: ~8.3 GB peak per concurrent build."
echo "  this host fits ~$((RAM_GB / 9)) concurrent build(s); see docs/sizing.md."
echo "Record what you measured in docs/sizing.md so the next person does not re-derive it."
