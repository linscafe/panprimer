"""Primer binding-site search across a haplotype sequence.

Two backends behind one interface (`find_binding_sites`):

* **naive** — pure Python, no deps. O(N*L) sliding window with a mismatch budget. Correct
  and exact; used for tiny fixtures, small projected windows, and the correctness tests.
* **bwa** — for genome-wide search on multi-Gb assemblies (added when the conda env is
  built). `bwa aln -n` reports short-seed hits with mismatches; we re-derive 3' offsets
  from the CIGAR/MD. Not needed for the subset correctness path.

Strand convention (reference given as top strand T, 5'->3'):
  * PLUS  site: primer sequence matches T[i:i+L]; primer 3' base at i+L-1.
  * MINUS site: primer matches revcomp(T[i:i+L]); primer 3' base at i.
Mismatch offsets are measured from the primer 3' end (0 == 3' terminal base), because the
dropout model weights 3'-proximal mismatches most heavily.
"""
from __future__ import annotations

from .model import BindingSite, Strand

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _scan_strand(
    primer_name: str,
    primer: str,
    ref: str,
    haplotype_id: str,
    chrom: str,
    strand: Strand,
    max_mismatches: int,
) -> list[BindingSite]:
    """Slide `primer` (already oriented to match the top strand) across `ref`."""
    L = len(primer)
    sites: list[BindingSite] = []
    ref = ref.upper()
    primer = primer.upper()
    for i in range(0, len(ref) - L + 1):
        window = ref[i : i + L]
        mm = 0
        offsets: list[int] = []
        for k in range(L):
            if window[k] != primer[k]:
                mm += 1
                if mm > max_mismatches:
                    break
                # offset from the primer's 3' end
                offsets.append((L - 1) - k)
        if mm <= max_mismatches:
            sites.append(
                BindingSite(
                    primer_name=primer_name,
                    haplotype_id=haplotype_id,
                    chrom=chrom,
                    start=i,
                    end=i + L,
                    strand=strand,
                    mismatches=mm,
                    mismatch_offsets_3p=sorted(offsets),
                )
            )
    return sites


def find_binding_sites_naive(
    primer_name: str,
    primer_seq: str,
    ref: str,
    haplotype_id: str,
    chrom: str,
    max_mismatches: int,
) -> list[BindingSite]:
    """Find all binding sites of `primer_seq` in `ref` on both strands (naive backend)."""
    plus = _scan_strand(
        primer_name, primer_seq, ref, haplotype_id, chrom, Strand.PLUS, max_mismatches
    )
    # On the minus strand the primer matches revcomp(window); equivalently the primer's
    # revcomp matches the window on the top strand. We scan with the revcomp'd primer so
    # index math stays identical, then relabel the 3' offsets.
    rc_primer = revcomp(primer_seq)
    L = len(primer_seq)
    minus_raw = _scan_strand(
        primer_name, rc_primer, ref, haplotype_id, chrom, Strand.MINUS, max_mismatches
    )
    minus: list[BindingSite] = []
    for s in minus_raw:
        # rc_primer[k] corresponds to primer[L-1-k]; its offset from primer 3' end is k.
        offsets = sorted((L - 1) - o for o in s.mismatch_offsets_3p)
        minus.append(
            BindingSite(
                primer_name=s.primer_name,
                haplotype_id=s.haplotype_id,
                chrom=s.chrom,
                start=s.start,
                end=s.end,
                strand=Strand.MINUS,
                mismatches=s.mismatches,
                mismatch_offsets_3p=offsets,
            )
        )
    return plus + minus


def find_binding_sites(
    primer_name: str,
    primer_seq: str,
    ref: str,
    haplotype_id: str,
    chrom: str,
    max_mismatches: int,
    backend: str = "naive",
) -> list[BindingSite]:
    if backend == "naive":
        return find_binding_sites_naive(
            primer_name, primer_seq, ref, haplotype_id, chrom, max_mismatches
        )
    if backend == "bwa":  # pragma: no cover - requires the conda env
        raise NotImplementedError(
            "bwa backend is wired in modules/binding_search.nf once the conda env exists"
        )
    raise ValueError(f"unknown binding backend: {backend!r}")
