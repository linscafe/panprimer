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

TSV=config/samples.tsv
LOG=hprc-r2/prepare.log
mkdir -p hprc-r2/assemblies
log(){ echo "[$(date +'%F %T')] $*" | tee -a "$LOG"; }

log "=== prepare_haplotypes start ($(nproc) cpus) ==="
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
  # 6) minimap2 asm5 index for locus projection (cached so runs don't rebuild it in-process)
  if [ ! -s "$fa.mmi" ]; then
    log "$id minimap2 -d ..."
    minimap2 -x asm5 -d "$fa.mmi" "$fa" 2>>"$LOG" && log "$id mmi done" || log "$id MINIMAP2 INDEX FAILED"
  fi
  log "$id READY"
done

# point the pipeline at the indexed plain .fa (leave the https source_url column intact;
# local paths have no spaces, so [^[:space:]] safely matches only the local_path column)
sed -i -E 's#(hprc-r2/assemblies/[^[:space:]]+)\.fa\.gz#\1.fa#' "$TSV"
log "samples.tsv local_path -> .fa"
log "=== prepare_haplotypes complete ==="
