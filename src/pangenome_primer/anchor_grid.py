"""Sparse CHM13<->haplotype anchor grid: locus projection without a per-haplotype index.

The problem this replaces. `project_locus` aligns the CHM13 template against a whole
haplotype, which needs a minimap2 index of that haplotype -- 5.80 GB on disk per haplotype,
the single largest artifact left after Phase 3. Building it in-process instead costs minutes
per haplotype per run, so it cannot simply be dropped.

The observation. We do not need to be able to align *anything* to the haplotype. We need to
answer one narrow question -- "where does this CHM13 interval live on this haplotype?" -- for
intervals a few kb wide. A coarse coordinate correspondence answers it to within a few kb,
and a few kb of sequence can then be aligned exactly, in memory, in milliseconds.

So: sample 1 kb probes every 10 kb along the haplotype, map them **against the shared CHM13
index** (one 5.9 GB index for all 464 haplotypes, not one per haplotype), and keep the
resulting (CHM13 position -> haplotype position) pairs. That is ~303k anchors, a few MB
gzipped. At query time, bracket the target between anchors, fetch just that window out of the
BGZF, and align the template to it with an in-memory `mappy.Aligner` -- the idiom
`mask.variability_counts` already uses.

Direction matters: probes go haplotype -> CHM13, so the big index is the *shared* one. The
reverse (CHM13 probes -> haplotype) would need a per-haplotype index, which is the thing
being eliminated.

Accuracy is unchanged where it counts. The final coordinate still comes from a base-level
`asm5` alignment of the real template against real haplotype sequence; the grid only decides
*which* few kb to align. A grid error does not silently shift the locus -- it makes the
window miss, the alignment fail, and the haplotype report `uncertain`, which is the same
failure mode as a failed whole-genome projection.
"""
from __future__ import annotations

import bisect
import gzip
import os
import subprocess
import tempfile
from dataclasses import dataclass

from .model import Locus
from .project import Projection
from .samples import sidecar_path

#: Probe length. Long enough for `asm5` to place uniquely through normal divergence, short
#: enough that 303k of them are only ~303 Mb of query.
PROBE_BP = 1_000
#: Probe spacing. The grid resolves the target to +/- one step before local realignment.
STEP_BP = 10_000
#: Contigs shorter than this are skipped. On HPRC R2 this drops 13 of 75 contigs while
#: keeping 99.97% of bases -- short unplaced contigs cost probes and map ambiguously.
MIN_CONTIG_BP = 100_000
#: Anchors below this mapping quality are discarded as unreliable placements.
MIN_MAPQ = 30
#: Padding applied to the bracketed window before realignment; absorbs local indels between
#: the flanking anchors.
WINDOW_PAD_BP = 20_000
#: Refuse to realign a window larger than this. A window this wide means the anchors
#: disagree (a scaffolding difference or a real structural change), and aligning megabases
#: in-process would blow the time budget silently. Better to report `uncertain`.
MAX_WINDOW_BP = 5_000_000


def grid_path(hap_fasta: str) -> str:
    """Anchor-grid sidecar for a haplotype, resolved like every other sidecar."""
    return sidecar_path(hap_fasta, ".anchors.tsv.gz")


@dataclass(frozen=True)
class Anchor:
    """One probe placement. `t_*` are CHM13 (target) coords, `q_*` haplotype (query)."""

    t_start: int
    t_end: int
    q_contig: str
    q_start: int
    q_end: int
    strand: str


# ----------------------------------------------------------------------------- build


def iter_probes(fa, *, probe_bp=PROBE_BP, step_bp=STEP_BP, min_contig=MIN_CONTIG_BP):
    """Yield (probe_name, sequence) along every contig long enough to be worth probing.

    Probes containing any `N` are skipped: assembly gaps map nowhere useful and would only
    add unplaceable queries to the minimap2 run.
    """
    for contig in fa.references:
        length = fa.get_reference_length(contig)
        if length < min_contig:
            continue
        for start in range(0, length - probe_bp + 1, step_bp):
            seq = fa.fetch(contig, start, start + probe_bp).upper()
            if "N" in seq:
                continue
            yield f"{contig}:{start}", seq


def build_grid(
    hap_fasta: str,
    chm13_mmi: str,
    out_path: str | None = None,
    *,
    threads: int = 4,
    probe_bp: int = PROBE_BP,
    step_bp: int = STEP_BP,
    min_contig: int = MIN_CONTIG_BP,
    min_mapq: int = MIN_MAPQ,
    progress=lambda s: None,
) -> dict:
    """Build the anchor grid for one haplotype. Returns build statistics.

    The probe FASTA is written to a temp file and removed afterwards; only the gzipped TSV
    survives. Writing to `<out>.part` and renaming means an interrupted build never leaves a
    truncated grid that would look usable -- the same discipline `prepare_haplotypes.sh`
    applies to `.mmi` and PAF outputs.
    """
    import pysam

    out_path = out_path or grid_path(hap_fasta)
    fa = pysam.FastaFile(hap_fasta)
    try:
        n_probes = 0
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tmp:
            probe_fa = tmp.name
            for name, seq in iter_probes(
                fa, probe_bp=probe_bp, step_bp=step_bp, min_contig=min_contig
            ):
                tmp.write(f">{name}\n{seq}\n")
                n_probes += 1
        genome_bp = sum(fa.get_reference_length(c) for c in fa.references)
    finally:
        fa.close()

    progress(f"  {n_probes:,} probes -> minimap2 against the shared CHM13 index ...")
    try:
        proc = subprocess.run(
            ["minimap2", "-x", "asm5", "--secondary=no", "-t", str(threads),
             "-K", "100M", chm13_mmi, probe_fa],
            capture_output=True, text=True, check=True,
        )
    finally:
        os.unlink(probe_fa)

    anchors: list[tuple] = []
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qname, strand, tname = f[0], f[4], f[5]
        if int(f[11]) < min_mapq:
            continue
        contig, _, off = qname.rpartition(":")   # contig names contain '#' but never ':'
        probe_off = int(off)
        # PAF q_start/q_end are offsets *within the probe*; lift them to contig coords.
        anchors.append(
            (tname, int(f[7]), int(f[8]), contig,
             probe_off + int(f[2]), probe_off + int(f[3]), strand)
        )
    anchors.sort(key=lambda a: (a[0], a[1]))

    tmp_out = out_path + ".part"
    with gzip.open(tmp_out, "wt") as fh:
        fh.write(
            f"# anchor grid: probe={probe_bp} step={step_bp} min_contig={min_contig} "
            f"min_mapq={min_mapq}\n"
            f"# probes={n_probes} anchored={len(anchors)} genome_bp={genome_bp}\n"
            "# t_name\tt_start\tt_end\tq_contig\tq_start\tq_end\tstrand\n"
        )
        for a in anchors:
            fh.write("\t".join(map(str, a)) + "\n")
    os.replace(tmp_out, out_path)

    anchored_frac = len(anchors) / n_probes if n_probes else 0.0
    return {
        "path": out_path,
        "probes": n_probes,
        "anchors": len(anchors),
        "anchored_fraction": anchored_frac,
        # Phase 5 needs the unanchored fraction: probes that found no confident CHM13
        # placement mark haplotype sequence that CHM13-once discovery cannot reason about.
        "unanchored_fraction": 1.0 - anchored_frac,
        "genome_bp": genome_bp,
        "bytes": os.path.getsize(out_path),
    }


# ------------------------------------------------------------------------------ load


class Grid:
    """Anchors for one haplotype, indexed by CHM13 chromosome and sorted by CHM13 start."""

    def __init__(self, by_chrom: dict[str, list[Anchor]], meta: dict):
        self._by_chrom = by_chrom
        self._starts = {c: [a.t_start for a in v] for c, v in by_chrom.items()}
        self.meta = meta

    def __bool__(self) -> bool:
        return bool(self._by_chrom)

    def near(self, chrom: str, start: int, end: int, flank: int) -> list[Anchor]:
        """Anchors on `chrom` whose CHM13 interval falls within `flank` of [start, end)."""
        arr = self._by_chrom.get(chrom)
        if not arr:
            return []
        starts = self._starts[chrom]
        lo = bisect.bisect_left(starts, start - flank)
        hi = bisect.bisect_right(starts, end + flank)
        return arr[lo:hi]


_CACHE: dict[tuple[str, float], Grid] = {}


def load_grid(path: str) -> Grid:
    """Load and cache a grid. Keyed on (path, mtime) so a rebuilt grid is picked up rather
    than served stale from a long-lived process."""
    key = (path, os.path.getmtime(path))
    hit = _CACHE.get(key)
    if hit is not None:
        return hit
    by_chrom: dict[str, list[Anchor]] = {}
    meta: dict = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                for tok in line[1:].split():
                    k, _, v = tok.partition("=")
                    if v:
                        meta[k] = v
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            by_chrom.setdefault(f[0], []).append(
                Anchor(int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]), f[6])
            )
    for v in by_chrom.values():
        v.sort(key=lambda a: a.t_start)
    grid = Grid(by_chrom, meta)
    _CACHE[key] = grid
    return grid


# --------------------------------------------------------------------------- project


def window_for(
    grid: Grid, chrom: str, tstart: int, tend: int, *, pad: int = WINDOW_PAD_BP
) -> tuple[str, int, int] | None:
    """Bracket a CHM13 interval into a (haplotype_contig, start, end) window, or None.

    The contig is chosen by majority vote over nearby anchors rather than by taking the
    single closest one: a lone probe that mis-placed into a segmental duplication would
    otherwise drag the whole window to the wrong contig, and the realignment would then
    quietly succeed against a paralogue.
    """
    flank = max(pad, 5 * STEP_BP)
    near = grid.near(chrom, tstart, tend, flank)
    if not near:
        return None
    votes: dict[str, int] = {}
    for a in near:
        votes[a.q_contig] = votes.get(a.q_contig, 0) + 1
    contig = max(votes, key=lambda c: votes[c])
    kept = [a for a in near if a.q_contig == contig]

    # Majority vote settles the contig but not the position: a probe that mis-placed
    # *within* the right contig (a tandem repeat, a segmental duplication) would stretch the
    # min/max window across everything between the true locus and the bad anchor. Drop
    # anchors far from the median before taking the extremes. The threshold is deliberately
    # loose -- an order of magnitude beyond the anchors we asked for -- so genuine large
    # indels and local rearrangements survive and only true outliers are cut.
    if len(kept) > 2:
        mids = sorted((a.q_start + a.q_end) // 2 for a in kept)
        median = mids[len(mids) // 2]
        limit = 10 * (flank + pad)
        clustered = [a for a in kept if abs((a.q_start + a.q_end) // 2 - median) <= limit]
        if clustered:
            kept = clustered

    lo = min(min(a.q_start, a.q_end) for a in kept) - pad
    hi = max(max(a.q_start, a.q_end) for a in kept) + pad
    return contig, max(0, lo), hi


def project_from_grid(
    grid_file: str,
    chrom: str,
    tstart: int,
    tend: int,
    template_seq: str,
    hap_fasta: str,
    *,
    fa=None,
    min_frac: float = 0.8,
    min_mapq: int = 20,
) -> Projection:
    """Project a CHM13 interval onto a haplotype via the anchor grid.

    Mirrors `project_locus`'s contract exactly -- same `min_frac`/`min_mapq` gating, same
    unrotated haplotype subsequence in `haplotype_seq` -- so swapping the `.mmi` path for
    this one does not shift any downstream result. The only difference is that the aligner is
    built over a bracketed window instead of the whole assembly.
    """
    import mappy
    import pysam

    grid = load_grid(grid_file)
    win = window_for(grid, chrom, tstart, tend)
    if win is None:
        return Projection(None, reason="no anchor covers the target on this haplotype")
    contig, lo, hi = win

    owns = fa is None
    if owns:
        fa = pysam.FastaFile(hap_fasta)
    try:
        clen = fa.get_reference_length(contig)
        lo, hi = max(0, lo), min(clen, hi)
        if hi - lo > MAX_WINDOW_BP:
            return Projection(
                None,
                reason=(f"anchor window {(hi - lo) / 1e6:.1f} Mb exceeds the "
                        f"{MAX_WINDOW_BP / 1e6:.0f} Mb cap; anchors disagree"),
            )
        window = fa.fetch(contig, lo, hi).upper()
    finally:
        if owns:
            fa.close()

    aligner = mappy.Aligner(seq=window, preset="asm5")
    if not aligner:
        return Projection(None, reason="failed to build the in-memory window aligner")
    hits = [
        h for h in aligner.map(template_seq)
        if h.mapq >= min_mapq and (h.q_en - h.q_st) >= min_frac * len(template_seq)
    ]
    if not hits:
        return Projection(None, reason="target did not map confidently to this haplotype")
    best = max(hits, key=lambda h: h.q_en - h.q_st)
    # Window-local -> contig-global coordinates.
    start, end = lo + best.r_st, lo + best.r_en
    return Projection(
        Locus(hap_fasta, contig, start, end),
        haplotype_seq=window[best.r_st:best.r_en],
        reason="anchor-grid",
    )
