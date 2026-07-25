"""CHM13-once candidate discovery: scan the reference once, not every haplotype.

The observation. Off-target amplification needs *both* primers bound within `max_product`
bp of each other and pointing inwards. Any locus that can produce a product is therefore
homologous to somewhere a primer already binds. So the set of places worth examining on a
haplotype is not the whole assembly -- it is the handful of windows homologous to the places
those primers bind on CHM13.

So: scan CHM13 **once per run**, merge the hits into candidate regions, lift each region onto
each haplotype through the Phase 4 anchor grid, and run the *exact* comparator over only
those windows. For the demo's 8 primers that is ~1 Mb examined per haplotype instead of
3.03 Gb -- roughly three orders of magnitude less sequence.

The many-to-one lift is a feature, not a compromise: a haplotype carrying three CYP2D7-like
copies that all align to one CHM13 locus yields three windows, and all three are scanned.

**Correctness envelope.** Sequence present in a haplotype but absent from CHM13 is
unreachable this way, and Phase 4 measured that at ~8% of each assembly (91.7-92.2% of
probes anchor). A primer binding only in haplotype-private sequence will not be found. That
is why this is a *mode* and not a replacement: `search.scope: exhaustive` remains the
default and the reference this is validated against. See `docs/implementation-plan.md`
Phase 5.

Why the exhaustive path is still the default. Phase 2 made a full scan cost ~6 s/haplotype,
so the saving here is ~6 s x (N-1) -- about 12 s at N=3, but ~46 min at N=464. Trading an
8% blind spot for 12 s is a bad deal; trading it for 46 min may not be. The threshold is a
judgement the caller makes, not one this module should make for them.
"""
from __future__ import annotations

from dataclasses import dataclass

from . import anchor_grid
from .model import BindingSite

#: Extra bp added each side of a CHM13 hit when forming a candidate region. A product needs
#: both primers within `max_product`, so a partner site can sit that far away; padding by it
#: guarantees the region spans any amplicon the hit could participate in.
DEFAULT_REGION_PAD = 2_000

#: Slack added each side when lifting a region onto a haplotype. Absorbs the difference
#: between CHM13 and haplotype coordinates left over after the anchor translation -- indels
#: between the nearest anchor and the region. Generous relative to the ~10 kb anchor spacing
#: it has to cover, because a window that misses costs a site silently.
LIFT_PAD = 5_000

#: Spacer between concatenated windows. `N` is a mismatch everywhere in this codebase, so a
#: run this long cannot be spanned within any realistic mismatch budget; hits that still clip
#: a boundary are rejected explicitly rather than left to the budget.
SEPARATOR = "N" * 64


@dataclass(frozen=True)
class Region:
    """A candidate interval in CHM13 coordinates."""

    chrom: str
    start: int
    end: int


def merge_regions(
    sites: dict[str, list[BindingSite]], *, pad: int = DEFAULT_REGION_PAD
) -> list[Region]:
    """Collapse CHM13 binding sites into padded, non-overlapping candidate regions.

    Merging matters for cost, not correctness: primer pairs bind close together by
    construction, so their padded windows overlap heavily and would otherwise be fetched and
    scanned several times over.
    """
    spans: dict[str, list[tuple[int, int]]] = {}
    for site_list in sites.values():
        for s in site_list:
            spans.setdefault(s.chrom, []).append(
                (max(0, s.start - pad), s.end + pad)
            )
    out: list[Region] = []
    for chrom, iv in spans.items():
        iv.sort()
        cur_lo, cur_hi = iv[0]
        for lo, hi in iv[1:]:
            if lo <= cur_hi:                 # touching or overlapping -> extend
                cur_hi = max(cur_hi, hi)
            else:
                out.append(Region(chrom, cur_lo, cur_hi))
                cur_lo, cur_hi = lo, hi
        out.append(Region(chrom, cur_lo, cur_hi))
    out.sort(key=lambda r: (r.chrom, r.start))
    return out


def _merge_windows(windows: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Collapse overlapping (contig, lo, hi) haplotype windows."""
    by_contig: dict[str, list[tuple[int, int]]] = {}
    for contig, lo, hi in windows:
        by_contig.setdefault(contig, []).append((lo, hi))
    out: list[tuple[str, int, int]] = []
    for contig, iv in by_contig.items():
        iv.sort()
        cur_lo, cur_hi = iv[0]
        for lo, hi in iv[1:]:
            if lo <= cur_hi:
                cur_hi = max(cur_hi, hi)
            else:
                out.append((contig, cur_lo, cur_hi))
                cur_lo, cur_hi = lo, hi
        out.append((contig, cur_lo, cur_hi))
    return out


def haplotype_windows(
    regions: list[Region], grid_file: str
) -> tuple[list[tuple[str, int, int]], int]:
    """Lift CHM13 regions onto one haplotype. Returns (merged windows, n_unanchored).

    Uses `anchor_grid.lift_interval`, not `window_for` or `project_from_grid`. We need a
    window *containing* the homologous sequence, not a precise locus, because the exact
    comparator runs over whatever is fetched -- so no per-region realignment is needed. But
    the window must also be tight: `window_for` spans every anchor within 50 kb, which turned
    1,669 regions into 212 Mb of sequence to rescan and left this mode only ~2x faster than
    just scanning the genome. `lift_interval` returns interval-plus-pad instead.

    `n_unanchored` counts regions the grid could not place. Those are the correctness
    envelope made countable -- the caller should surface it rather than let it pass silently.
    """
    grid = anchor_grid.load_grid(grid_file)
    windows: list[tuple[str, int, int]] = []
    unanchored = 0
    for r in regions:
        win = anchor_grid.lift_interval(grid, r.chrom, r.start, r.end, pad=LIFT_PAD)
        if win is None:
            unanchored += 1
            continue
        windows.append(win)
    return _merge_windows(windows), unanchored


def find_binding_sites_in_windows(
    seqs: list[str],
    hap_fasta: str,
    haplotype_id: str,
    max_mismatches: int,
    windows: list[tuple[str, int, int]],
    *,
    fa=None,
) -> dict[str, list[BindingSite]]:
    """Run the exact comparator over `windows` only, returning contig-global coordinates.

    Same return shape as `search.find_binding_sites_batch`: `{primer_sequence: [BindingSite]}`
    keyed by sequence, `primer_name` set to the sequence, deduplicated on
    `(chrom, start, strand)`. Within a window the search is exhaustive, so a site is missed
    only by not being in any window -- never by an approximation inside one.
    """
    import bisect

    import pysam

    from . import rust_backend

    uniq = sorted(set(seqs))
    out: dict[str, list[BindingSite]] = {s: [] for s in uniq}
    seen: dict[str, set[tuple[str, int, str]]] = {s: set() for s in uniq}
    if not windows:
        return out

    owns = fa is None
    if owns:
        fa = pysam.FastaFile(hap_fasta)
    try:
        # Concatenate every window into ONE sequence and scan it once. The comparator builds
        # a seed index over the primers on each call, so per-window calls pay that fixed cost
        # thousands of times: 1,616 windows totalling only 22.9 Mb took 4.5 s, while the
        # scanning itself is ~0.05 s at ~2 ns/base. One call amortises the index build to
        # one, and 22.9 Mb of concatenated sequence costs a few MB of RAM.
        parts: list[str] = []
        starts: list[int] = []                      # concat offset of each window
        meta: list[tuple[str, int, int]] = []       # (contig, contig_lo, window_len)
        pos = 0
        for contig, lo, hi in windows:
            clen = fa.get_reference_length(contig)
            lo, hi = max(0, lo), min(clen, hi)
            if hi <= lo:
                continue
            ref = fa.fetch(contig, lo, hi).upper()
            parts.append(ref)
            starts.append(pos)
            meta.append((contig, lo, len(ref)))
            pos += len(ref)
            parts.append(SEPARATOR)
            pos += len(SEPARATOR)
        if not starts:
            return out
        joined = "".join(parts)

        if rust_backend.available():
            found = rust_backend.find_binding_sites_in_seq_batch(
                uniq, joined, haplotype_id, "concat", max_mismatches
            )
        else:
            from .binding import find_binding_sites_naive

            found = {
                s: find_binding_sites_naive(
                    s, s, joined, haplotype_id, "concat", max_mismatches
                )
                for s in uniq
            }

        for s, sites in found.items():
            for site in sites:
                i = bisect.bisect_right(starts, site.start) - 1
                if i < 0:
                    continue
                contig, contig_lo, wlen = meta[i]
                off = site.start - starts[i]
                # Reject anything overlapping a separator. `N` counts as a mismatch, so a
                # window straddling a boundary is usually rejected by the mismatch budget --
                # but one that clips only 1-3 Ns would not be, and would report a site at a
                # junction that does not exist in the genome.
                # Length must be captured BEFORE `start` is reassigned -- computing it from
                # the mutated field leaves `end` in concatenated coordinates while `start` is
                # in contig coordinates, which surfaces far downstream as a pysam
                # "start > stop" on an unrelated fetch.
                length = site.end - site.start
                if off < 0 or off + length > wlen:
                    continue
                site.chrom = contig
                site.start = contig_lo + off
                site.end = site.start + length
                key = (site.chrom, site.start, site.strand.value)
                if key in seen[s]:
                    continue                        # windows can overlap after clamping
                seen[s].add(key)
                out[s].append(site)
    finally:
        if owns:
            fa.close()
    return out


class Chm13OnceSearcher:
    """Holds the single CHM13 scan for a run and serves per-haplotype site lists.

    Built once by the caller, then queried per haplotype -- which is what makes the CHM13
    scan a fixed cost rather than one paid N times.
    """

    def __init__(
        self,
        seqs: list[str],
        chm13_fasta: str,
        max_mismatches: int,
        *,
        backend: str | None = None,
        region_pad: int = DEFAULT_REGION_PAD,
        progress=lambda s: None,
    ):
        from .search import find_binding_sites_batch

        progress(f"genome-wide search on CHM13 (once for all haplotypes) ...")
        chm13_sites = find_binding_sites_batch(
            seqs, chm13_fasta, "CHM13", max_mismatches, backend=backend
        )
        self.seqs = list(seqs)
        self.max_mismatches = max_mismatches
        self.regions = merge_regions(chm13_sites, pad=region_pad)
        self.n_chm13_sites = sum(len(v) for v in chm13_sites.values())
        progress(
            f"  {self.n_chm13_sites} CHM13 sites -> {len(self.regions)} candidate region(s)"
        )

    def sites_for(
        self, hap_fasta: str, haplotype_id: str, *, fa=None
    ) -> tuple[dict[str, list[BindingSite]], dict]:
        """Binding sites on one haplotype, plus a stats dict.

        Raises `FileNotFoundError` when the haplotype has no anchor grid: this mode cannot
        work without one, and silently falling back to a whole-genome scan would hide the
        very cost the caller opted in to avoid.
        """
        from pathlib import Path

        grid_file = anchor_grid.grid_path(hap_fasta)
        if not Path(grid_file).exists():
            raise FileNotFoundError(
                f"search.scope='chm13-once' needs an anchor grid for {haplotype_id}, but "
                f"{grid_file} is missing. Build it with `pangenome-primer build-anchor-grid` "
                f"or use search.scope='exhaustive'."
            )
        windows, unanchored = haplotype_windows(self.regions, grid_file)
        scanned_bp = sum(hi - lo for _, lo, hi in windows)
        sites = find_binding_sites_in_windows(
            self.seqs, hap_fasta, haplotype_id, self.max_mismatches, windows, fa=fa
        )
        return sites, {
            "regions": len(self.regions),
            "windows": len(windows),
            "unanchored_regions": unanchored,
            "scanned_bp": scanned_bp,
        }
