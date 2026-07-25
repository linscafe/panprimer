"""Synthetic BGZF mini-genome for Phase 0 backend-conformance tests.

Builds a few contigs of deterministic random sequence with primer binding sites planted
at KNOWN coordinates -- including sites at contig edges, sites within `slop` distance of
each other, multiple contigs, and both strands -- then bgzip-compresses and indexes it
exactly the way real HPRC haplotypes are indexed: `samtools faidx` on the `.gz` produces
both the `.fai` and the BGZF `.gzi` sidecar. Nothing here is committed; `build_mini_genome`
writes fresh files under whatever `tmp_path` it is given, at test time.

Sizing note: `binding.find_binding_sites_naive` is a pure-Python O(N*L) sliding-window
scan -- correct by construction, and explicitly documented as meant "for tiny fixtures,
small projected windows" rather than genome-scale search (that's exactly the gap the
compiled rust backend fills). Measured on this machine, scanning a genuine multi-Mb
genome with it took >100s for the conformance suite's battery of assertions -- far past
CI's budget, and itself a small illustration of the runtime problem the optimization plan
exists to fix. So CONTIG_LEN below is sized in the hundreds-of-kb, not multi-Mb: still
multiple contigs with edge and close-pair coverage, still exercises real BGZF/`.gzi`
random access, but keeps the naive-backend conformance pass in the low single-digit
seconds. `test_backend_conformance.py` caches one scan per backend rather than
re-scanning per assertion, for the same reason.

Requires `bgzip` + `samtools` on PATH (both ship in the `pangenome-primer` conda env).
"""
from __future__ import annotations

import random
import subprocess
from dataclasses import dataclass
from pathlib import Path

from pangenome_primer.binding import revcomp
from pangenome_primer.model import Strand

BASES = "ACGT"
PRIMER_LEN = 20
DEFAULT_SEED = 20260724

CONTIG_LEN = 150_000
CONTIGS = ("ctgA", "ctgB", "ctgC")


@dataclass(frozen=True)
class PlantedSite:
    """A binding site planted at a known coordinate, using `binding.py`'s convention:
    PLUS matches T[start:end] exactly (once `mismatches` engineered flips are applied);
    MINUS matches revcomp(T[start:end])."""

    label: str
    primer_seq: str  # 5'->3' primer sequence
    chrom: str
    start: int  # 0-based, haplotype coordinates
    end: int  # start + len(primer_seq)
    strand: Strand
    mismatches: int = 0  # engineered mismatch count vs. the top-strand window

    @property
    def expected_offsets_3p(self) -> list[int]:
        """Mismatches are planted at the deep (5'-most-from-3') end of the site by
        construction (see `_plant_plus`/`_plant_minus`), so they land at the *largest*
        3'-offsets: range(L - mismatches, L)."""
        L = len(self.primer_seq)
        return sorted(range(L - self.mismatches, L))


@dataclass
class MiniGenome:
    fasta_gz: str  # bgzip-compressed FASTA path (.fai/.gzi are siblings)
    contigs: dict[str, str]  # chrom -> full top-strand sequence, for cross-checks
    planted: list[PlantedSite]
    haplotype_id: str = "MINI#hap1"


def _rand_seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choices(BASES, k=n))


def _flip(base: str) -> str:
    """Return a base guaranteed to differ from `base` (never special-cases N)."""
    return "C" if base != "C" else "G"


def _plant_plus(
    contigs: dict[str, list[str]],
    planted: list[PlantedSite],
    chrom: str,
    start: int,
    primer_seq: str,
    mismatches: int = 0,
) -> None:
    """Write `primer_seq` (optionally with `mismatches` flips at its 5' end, i.e. the
    largest 3'-offsets) into the top strand at [start, start+len) -- a PLUS site."""
    L = len(primer_seq)
    top = list(primer_seq)
    for k in range(mismatches):
        top[k] = _flip(top[k])
    contigs[chrom][start : start + L] = top
    planted.append(
        PlantedSite(
            label=f"{chrom}:{start}_plus_mm{mismatches}",
            primer_seq=primer_seq,
            chrom=chrom,
            start=start,
            end=start + L,
            strand=Strand.PLUS,
            mismatches=mismatches,
        )
    )


def _plant_minus(
    contigs: dict[str, list[str]],
    planted: list[PlantedSite],
    chrom: str,
    start: int,
    primer_seq: str,
    mismatches: int = 0,
) -> None:
    """Write revcomp(primer_seq) (optionally with `mismatches` flips at the *last*
    indices, which -- per the minus-strand offset relabeling in binding.py -- also land
    at the largest 3'-offsets) into the top strand at [start, start+len) -- a MINUS site."""
    L = len(primer_seq)
    top = list(revcomp(primer_seq))
    for k in range(L - mismatches, L):
        top[k] = _flip(top[k])
    contigs[chrom][start : start + L] = top
    planted.append(
        PlantedSite(
            label=f"{chrom}:{start}_minus_mm{mismatches}",
            primer_seq=primer_seq,
            chrom=chrom,
            start=start,
            end=start + L,
            strand=Strand.MINUS,
            mismatches=mismatches,
        )
    )


def build_mini_genome(tmp_path: Path, seed: int = DEFAULT_SEED) -> MiniGenome:
    """Write a bgzip-compressed, faidx+gzi-indexed synthetic FASTA under `tmp_path`, with
    primer binding sites planted at known coordinates. Returns a `MiniGenome` describing
    where everything landed."""
    rng = random.Random(seed)
    contigs: dict[str, list[str]] = {
        c: list(_rand_seq(CONTIG_LEN, rng)) for c in CONTIGS
    }
    planted: list[PlantedSite] = []

    # --- contig-edge sites -----------------------------------------------------------
    # ctgA, position 0: PLUS, exact match, right at the contig start.
    _plant_plus(contigs, planted, "ctgA", 0, "TTGCACAGTCCAGATTGCAA")
    # ctgA, last 20 bases: MINUS, exact match, right at the contig end.
    _plant_minus(contigs, planted, "ctgA", CONTIG_LEN - PRIMER_LEN, "CATTGCGATTGACACTTGCG")
    # ctgC, position 0: MINUS at the contig start.
    _plant_minus(contigs, planted, "ctgC", 0, "GGCTTAACGGTCCAATGCTA")
    # ctgC, last 20 bases: PLUS with 1 engineered mismatch, right at the contig end.
    _plant_plus(
        contigs, planted, "ctgC", CONTIG_LEN - PRIMER_LEN, "AACGGTTGCAGGTTACCGGA",
        mismatches=1,
    )

    # --- ordinary interior sites, both strands, both remaining contigs ---------------
    _plant_plus(contigs, planted, "ctgB", 30_000, "GCGTGGGGGAAGATGGTTTA")
    _plant_minus(contigs, planted, "ctgB", 90_000, "CAGGGACCCATCCTCAAACC")
    _plant_plus(contigs, planted, "ctgA", 80_000, "ACTCCTGGGCTCAAGCAATC", mismatches=2)

    # --- a close pair: two independent PLUS sites separated by a 2 bp gap (< slop=3) --
    _plant_plus(contigs, planted, "ctgB", 120_000, "AGGCTGAACCTGATGCAGTT")
    _plant_plus(contigs, planted, "ctgB", 120_000 + PRIMER_LEN + 2, "TCCAGGTTAACGGCATTGAC")

    joined = {c: "".join(seq) for c, seq in contigs.items()}

    fa_path = Path(tmp_path) / "mini.fa"
    with open(fa_path, "w") as fh:
        for c in CONTIGS:
            fh.write(f">{c}\n{joined[c]}\n")

    subprocess.run(["bgzip", "-f", str(fa_path)], check=True)
    gz_path = fa_path.parent / (fa_path.name + ".gz")  # bgzip appends .gz to "mini.fa"
    subprocess.run(["samtools", "faidx", str(gz_path)], check=True)

    return MiniGenome(fasta_gz=str(gz_path), contigs=joined, planted=planted)
