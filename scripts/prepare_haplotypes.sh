#!/usr/bin/env bash
# Prepare the HPRC R2 haplotype subset: download -> md5 -> BGZF faidx, for every haplotype
# in config/samples.tsv. Fully resumable (each step is skipped if its output already
# exists), logs to hprc-r2/prepare.log.
#
# Storage (Phase 3). The default path keeps what HPRC already ships plus two small indexes
# -- 0.90 GB .fa.gz + .fa.gz.fai + .fa.gz.gzi -- and the 5.80 GB minimap2 projection index,
# so ~6.7 GB per haplotype rather than the previous ~15 GB. The rust search backend streams
# the BGZF directly and pysam random-accesses it through the .fai/.gzi pair, so neither the
# 3.08 GB uncompressed .fa nor the 5.31 GB bwa index is needed.
#
# The .mmi is still the single largest item and is still built here: projection has no
# replacement until the sparse anchor grid lands (Phase 4). Dropping it now would not save
# space -- project.py would just rebuild an equivalent index in-process on every run. Once
# Phase 4 lands this falls to ~0.91 GB/haplotype, i.e. ~0.4 TB rather than ~7 TB across the
# 464-haplotype HPRC R2.
#
# WITH_BWA=1 restores the old behaviour (gunzip + bwa index, ~8.4 GB/haplotype). It is
# needed only to run `search.backend: bwa`, which is kept as the reference implementation
# the rust path was validated against -- not for normal use.
#
# On completion, rewrites the samples.tsv local_path column to whichever file is now the
# canonical one (.fa.gz by default, .fa under WITH_BWA=1).
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
  log "building CHM13 asm5 index ..."
  if minimap2 -x asm5 -d "$CHM13_MMI.part" "$CHM13" 2>>"$LOG"; then
    mv "$CHM13_MMI.part" "$CHM13_MMI"
  else
    log "CHM13 asm5 index FAILED"; rm -f "$CHM13_MMI.part"
  fi
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
    if curl -fsS -o "$gz.part" "$url"; then mv "$gz.part" "$gz"; else log "$id DOWNLOAD FAILED"; rm -f "$gz.part"; continue; fi
  fi
  # 2) md5 check against the pinned manifest value (only if we still have the .gz)
  if [ -s "$gz" ] && [ "${md5:-PENDING}" != "PENDING" ]; then
    got=$(md5sum "$gz" | cut -d' ' -f1)
    [ "$got" = "$md5" ] && log "$id md5 ok" || log "$id MD5 MISMATCH got=$got want=$md5"
  fi
  # 3) BGZF index: .fa.gz.fai + .fa.gz.gzi, ~0.76 MB together. This is what replaces the
  # gunzip+faidx pair — pysam random-accesses the compressed file through them, so the
  # 3.08 GB .fa is never materialised. samtools writes both, and refuses outright if the
  # file is plain gzip rather than BGZF, so a non-BGZF source fails loudly here instead of
  # silently degrading later.
  if [ ! -s "$gz.gzi" ] || [ ! -s "$gz.fai" ]; then
    log "$id faidx (bgzf) ..."
    samtools faidx "$gz" || { log "$id BGZF FAIDX FAILED"; rm -f "$gz.fai" "$gz.gzi"; continue; }
  fi
  # 4) uncompressed .fa + bwa index — ONLY for the reference bwa backend (WITH_BWA=1).
  # Together they cost 8.39 GB/haplotype and ~55 min of bwa indexing; the default rust
  # backend needs neither.
  if [ "${WITH_BWA:-0}" = "1" ]; then
    if [ ! -s "$fa" ]; then
      log "$id gunzip (WITH_BWA) ..."
      gunzip -kf "$gz" || { log "$id GUNZIP FAILED"; rm -f "$fa"; continue; }
    fi
    if [ ! -s "$fa.fai" ]; then
      log "$id faidx ..."
      samtools faidx "$fa" || { log "$id FAIDX FAILED"; rm -f "$fa.fai"; continue; }
    fi
    # bwa writes .bwt/.sa/.pac/.amb/.ann directly (no single renameable output), so on
    # failure remove all of them — otherwise a partial .bwt looks "present" to a resumed run.
    if [ ! -s "$fa.bwt" ]; then
      t0=$(date +%s); log "$id bwa index ..."
      if bwa index "$fa" 2>>"$LOG"; then
        log "$id bwa index done ($(( ($(date +%s)-t0)/60 )) min)"
      else
        log "$id BWA INDEX FAILED"; rm -f "$fa.bwt" "$fa.sa" "$fa.pac" "$fa.amb" "$fa.ann"; continue
      fi
    else
      log "$id already indexed (bwa)"
    fi
  fi
  # 6) projection cache. minimap2 reads the BGZF directly, so this works with or without
  # the uncompressed .fa. Sidecars are written next to whichever file is canonical, and an
  # already-built .fa-named one is reused rather than rebuilt — matching the resolution
  # order in `samples.sidecar_path`, which is what the Python side uses to find them.
  src="$gz"; [ -s "$fa" ] && src="$fa"
  paf="$src.chm13.paf"; [ -s "$fa.chm13.paf" ] && paf="$fa.chm13.paf"
  mmi="$src.mmi";       [ -s "$fa.mmi" ]       && mmi="$fa.mmi"
  # Default = a cheap, memory-safe minimap2 index (.mmi) that the per-query projection
  # loads in seconds. The whole-genome PAF (instant per-query lift, but a ~1-2 hr, >15 GB
  # alignment per haplotype) is opt-in via BUILD_PAF=1 and only worth it on high-RAM /
  # cloud. project_target prefers a PAF when present, else the .mmi.
  if [ -s "$paf" ]; then
    log "$id projection PAF present"
  elif [ "${BUILD_PAF:-0}" = "1" ] && [ -s "${CHM13_MMI:-/nonexistent}" ]; then
    # whole-genome alignment: memory scales with threads. MM_THREADS=4 fits 15 GB; bump it
    # on high-RAM / cloud. This is the instant-per-query-lift path (item 1) at scale.
    t0=$(date +%s); log "$id projection PAF (BUILD_PAF, -c, -t ${MM_THREADS:-4}) ..."
    # -c emits base-level CIGAR (cg:Z:) so projection lifts coordinates precisely; -K 100M
    # caps minimap2's per-batch query buffer so peak RSS stays bounded (unbounded -K hit
    # 14.5 GB / OOM'd on a 15 GB box — see hprc-r2/prepare.log)
    if minimap2 -c -x asm5 -t "${MM_THREADS:-4}" --secondary=no -K 100M "$CHM13_MMI" "$src" > "$paf.part" 2>>"$LOG"; then
      mv "$paf.part" "$paf"; log "$id PAF done ($(( ($(date +%s)-t0)/60 )) min)"
    else log "$id PAF FAILED"; rm -f "$paf.part"; fi
  elif [ ! -s "$mmi" ]; then
    t0=$(date +%s); log "$id minimap2 -d (projection index) ..."
    # write to .part then rename so an interrupted build never leaves a truncated .mmi that
    # loads as an empty aligner (would turn every projection uncertain)
    if minimap2 -x asm5 -d "$mmi.part" "$src" 2>>"$LOG"; then
      mv "$mmi.part" "$mmi"; log "$id mmi done ($(( ($(date +%s)-t0)/60 )) min)"
    else log "$id MMI FAILED"; rm -f "$mmi.part"; fi
  fi
  log "$id READY"
done

# Point the pipeline at whichever file is now canonical. Leave the https source_url column
# intact: local paths have no spaces, so anchoring on the hprc-r2/assemblies/ prefix and
# matching [^[:space:]] keeps the substitution inside the local_path column.
if [ "${WITH_BWA:-0}" = "1" ]; then
  sed -i -E 's#(hprc-r2/assemblies/[^[:space:]]+)\.fa\.gz#\1.fa#' "$TSV"
  log "samples.tsv local_path -> .fa (WITH_BWA)"
else
  sed -i -E 's#(hprc-r2/assemblies/[^[:space:]]+)\.fa(\.gz)?([[:space:]]|$)#\1.fa.gz\3#' "$TSV"
  log "samples.tsv local_path -> .fa.gz"
fi
log "=== prepare_haplotypes complete ==="
