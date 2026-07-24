#!/usr/bin/env bash
# Prepare the HPRC R2 haplotype subset for the pipeline: download -> gunzip -> faidx ->
# bwa index, for every haplotype in config/samples.tsv. Fully resumable (each step is
# skipped if its output already exists), logs to hprc-r2/prepare.log. Designed to run
# unattended (e.g. overnight) since bwa-indexing ~3 Gb assemblies is the slow step.
#
# On completion, rewrites the samples.tsv local_path column from .fa.gz to the indexed .fa.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pangenome-primer

TSV="${1:-config/samples.tsv}"
CHM13="${2:-hprc-r2/references/chm13v2.0.fa}"
LOG=hprc-r2/prepare.log
mkdir -p hprc-r2/assemblies
log(){ echo "[$(date +'%F %T')] $*" | tee -a "$LOG"; }

log "=== prepare_haplotypes start ($(nproc) cpus) TSV=$TSV ==="
# CHM13 minimap2 index, built once, reused for every haplotype's projection PAF
CHM13_MMI="${CHM13%.fa}.asm5.mmi"
if [ -s "$CHM13" ] && [ ! -s "$CHM13_MMI" ]; then
  log "building CHM13 asm5 index ..."; minimap2 -x asm5 -d "$CHM13_MMI" "$CHM13" 2>>"$LOG"
fi
# columns: sample hap superpop local_path md5 release population source_url genbank
grep -vE '^#|^sample\b' "$TSV" | while IFS=$'\t' read -r sample hap superpop local_path md5 release population url genbank; do
  [ -z "${sample:-}" ] && continue
  base="${local_path%.fa.gz}"; base="${base%.fa}"
  gz="$base.fa.gz"; fa="$base.fa"
  id="${sample}#hap${hap} ($superpop)"

  # 1) download
  if [ ! -s "$gz" ] && [ ! -s "$fa" ]; then
    log "$id downloading ..."
    if curl -fsS -o "$gz.part" "$url"; then mv "$gz.part" "$gz"; else log "$id DOWNLOAD FAILED"; continue; fi
  fi
  # 2) md5 check against the pinned manifest value (only if we still have the .gz)
  if [ -s "$gz" ] && [ "${md5:-PENDING}" != "PENDING" ]; then
    got=$(md5sum "$gz" | cut -d' ' -f1)
    [ "$got" = "$md5" ] && log "$id md5 ok" || log "$id MD5 MISMATCH got=$got want=$md5"
  fi
  # 3) gunzip
  if [ ! -s "$fa" ]; then log "$id gunzip ..."; gunzip -kf "$gz" || { log "$id GUNZIP FAILED"; continue; }; fi
  # 4) faidx
  if [ ! -s "$fa.fai" ]; then log "$id faidx ..."; samtools faidx "$fa"; fi
  # 5) bwa index (genome-wide off-target search; the slow step)
  if [ ! -s "$fa.bwt" ]; then
    t0=$(date +%s); log "$id bwa index ..."
    if bwa index "$fa" 2>>"$LOG"; then log "$id bwa index done ($(( ($(date +%s)-t0)/60 )) min)"; else log "$id BWA INDEX FAILED"; continue; fi
  else
    log "$id already indexed (bwa)"
  fi
  # 6) projection cache. Default = a cheap, memory-safe minimap2 index (.mmi) that the
  # per-query projection loads in seconds. The whole-genome PAF (instant per-query lift, but
  # a ~1-2 hr, >15 GB alignment per haplotype) is opt-in via BUILD_PAF=1 and only worth it
  # on high-RAM / cloud. project_target prefers a PAF when present, else the .mmi.
  if [ -s "$fa.chm13.paf" ]; then
    log "$id projection PAF present"
  elif [ "${BUILD_PAF:-0}" = "1" ] && [ -s "${CHM13_MMI:-/nonexistent}" ]; then
    # whole-genome alignment: memory scales with threads. MM_THREADS=4 fits 15 GB; bump it
    # on high-RAM / cloud. This is the instant-per-query-lift path (item 1) at scale.
    t0=$(date +%s); log "$id projection PAF (BUILD_PAF, -c, -t ${MM_THREADS:-4}) ..."
    # -c emits base-level CIGAR (cg:Z:) so projection lifts coordinates precisely
    if minimap2 -c -x asm5 -t "${MM_THREADS:-4}" --secondary=no "$CHM13_MMI" "$fa" > "$fa.chm13.paf.part" 2>>"$LOG"; then
      mv "$fa.chm13.paf.part" "$fa.chm13.paf"; log "$id PAF done ($(( ($(date +%s)-t0)/60 )) min)"
    else log "$id PAF FAILED"; rm -f "$fa.chm13.paf.part"; fi
  elif [ ! -s "$fa.mmi" ]; then
    t0=$(date +%s); log "$id minimap2 -d (projection index) ..."
    # write to .part then rename so an interrupted build never leaves a truncated .mmi that
    # loads as an empty aligner (would turn every projection uncertain)
    if minimap2 -x asm5 -d "$fa.mmi.part" "$fa" 2>>"$LOG"; then
      mv "$fa.mmi.part" "$fa.mmi"; log "$id mmi done ($(( ($(date +%s)-t0)/60 )) min)"
    else log "$id MMI FAILED"; rm -f "$fa.mmi.part"; fi
  fi
  log "$id READY"
done

# point the pipeline at the indexed plain .fa (leave the https source_url column intact;
# local paths have no spaces, so [^[:space:]] safely matches only the local_path column)
sed -i -E 's#(hprc-r2/assemblies/[^[:space:]]+)\.fa\.gz#\1.fa#' "$TSV"
log "samples.tsv local_path -> .fa"
log "=== prepare_haplotypes complete ==="
