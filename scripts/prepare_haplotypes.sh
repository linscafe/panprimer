#!/usr/bin/env bash
# Prepare the HPRC R2 haplotype subset: download -> md5 -> BGZF faidx -> anchor grid, for
# every haplotype in config/samples.tsv. Fully resumable (each step is skipped if its output already
# exists), logs to hprc-r2/prepare.log.
#
# Storage. Everything kept is what HPRC already ships plus small sidecars: the 0.90 GB
# .fa.gz, its .fa.gz.fai + .fa.gz.gzi (<1 MB), and a ~4 MB projection anchor grid --
# ~0.91 GB per haplotype, against ~15 GB before. The search backend streams the BGZF
# directly and pysam random-accesses it through the .fai/.gzi pair, so no uncompressed .fa
# and no per-haplotype search or alignment index is ever built. Across the 464-haplotype
# HPRC R2 that is ~0.4 TB rather than ~7 TB.
#
# On completion, rewrites the samples.tsv local_path column to the .fa.gz.
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
  # 6) projection cache. minimap2 reads the BGZF directly, so this works with or without
  # the uncompressed .fa. Sidecars are written next to whichever file is canonical, and an
  # already-built .fa-named one is reused rather than rebuilt — matching the resolution
  # order in `samples.sidecar_path`, which is what the Python side uses to find them.
  src="$gz"; [ -s "$fa" ] && src="$fa"
  paf="$src.chm13.paf"; [ -s "$fa.chm13.paf" ] && paf="$fa.chm13.paf"
  mmi="$src.mmi";       [ -s "$fa.mmi" ]       && mmi="$fa.mmi"
  grid="$src.anchors.tsv.gz"; [ -s "$fa.anchors.tsv.gz" ] && grid="$fa.anchors.tsv.gz"
  # Default = the sparse anchor grid (a few MB): probes mapped against the SHARED CHM13
  # index, so one big index serves every haplotype instead of one 5.80 GB .mmi each.
  # project_target prefers a PAF, then the grid, then a whole-haplotype .mmi. The .mmi
  # branch survives only for the case where no CHM13 index is available to build a grid
  # against; it is no longer the normal path. The whole-genome PAF (instant per-query lift,
  # but a ~1-2 hr, >15 GB alignment per haplotype) stays opt-in via BUILD_PAF=1.
  if [ -s "$paf" ]; then
    log "$id projection PAF present"
  elif [ -s "$grid" ]; then
    log "$id anchor grid present ($(( $(stat -c%s "$grid") / 1000000 )) MB)"
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
  elif [ -s "${CHM13_MMI:-/nonexistent}" ]; then
    t0=$(date +%s); log "$id anchor grid (-t ${MM_THREADS:-4}) ..."
    # build_grid writes to <out>.part and renames, so an interrupted build leaves no
    # half-written grid that would load as a sparse one and quietly mis-project.
    if python -c "
import sys
from pangenome_primer import anchor_grid
s = anchor_grid.build_grid(sys.argv[1], sys.argv[2], sys.argv[3], threads=int(sys.argv[4]))
print(f\"  {s['anchors']}/{s['probes']} anchored ({s['anchored_fraction']:.1%}), {s['bytes']/1e6:.1f} MB\")
" "$src" "$CHM13_MMI" "$grid" "${MM_THREADS:-4}" >>"$LOG" 2>&1; then
      log "$id anchor grid done ($(( ($(date +%s)-t0)/60 )) min)"
    else log "$id ANCHOR GRID FAILED"; rm -f "$grid.part"; fi
  elif [ ! -s "$mmi" ]; then
    t0=$(date +%s); log "$id minimap2 -d (projection index; no CHM13 index for a grid) ..."
    # write to .part then rename so an interrupted build never leaves a truncated .mmi that
    # loads as an empty aligner (would turn every projection uncertain)
    if minimap2 -x asm5 -d "$mmi.part" "$src" 2>>"$LOG"; then
      mv "$mmi.part" "$mmi"; log "$id mmi done ($(( ($(date +%s)-t0)/60 )) min)"
    else log "$id MMI FAILED"; rm -f "$mmi.part"; fi
  fi
  log "$id READY"
done

# Point the pipeline at the .fa.gz. Leave the https source_url column intact: local paths
# have no spaces, so anchoring on the hprc-r2/assemblies/ prefix and matching
# [^[:space:]] keeps the substitution inside the local_path column.
sed -i -E 's#(hprc-r2/assemblies/[^[:space:]]+)\.fa(\.gz)?([[:space:]]|$)#\1.fa.gz\3#' "$TSV"
log "samples.tsv local_path -> .fa.gz"
log "=== prepare_haplotypes complete ==="
