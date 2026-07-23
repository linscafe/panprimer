"""Core domain types. See CONTEXT.md for the ubiquitous language.

These dataclasses are the vocabulary every pipeline stage speaks. They carry no I/O and
no external-tool dependencies so they stay trivially testable.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class Strand(str, Enum):
    PLUS = "+"
    MINUS = "-"


class HaplotypeStatus(str, Enum):
    """Exactly one of these per (primer pair, haplotype). See CONTEXT.md."""

    PASS = "pass"
    DROPOUT = "dropout"
    OFF_TARGET = "off_target"
    MULTI_PRODUCT = "multi_product"
    UNCERTAIN = "uncertain"


@dataclass(frozen=True)
class Locus:
    """A genomic interval, half-open [start, end), 0-based. Coordinates are on whatever
    assembly `assembly` names — CHM13 v2.0 for the anchor, a haplotype id otherwise."""

    assembly: str
    chrom: str
    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 0 or self.end < self.start:
            raise ValueError(f"invalid interval {self.chrom}:{self.start}-{self.end}")

    @property
    def length(self) -> int:
        return self.end - self.start

    def overlaps(self, other: "Locus") -> bool:
        return (
            self.chrom == other.chrom
            and self.start < other.end
            and other.start < self.end
        )


@dataclass(frozen=True)
class Haplotype:
    """One phased HPRC R2 assembly. `sample` is the individual, `hap` is 1 or 2."""

    sample: str
    hap: int
    superpop: str = ""
    fasta: str = ""  # local path once fetched
    sha256: str = ""
    release: str = "R2"

    @property
    def id(self) -> str:
        return f"{self.sample}#hap{self.hap}"


@dataclass(frozen=True)
class Primer:
    """A single oligo. `sequence` is 5'->3'."""

    name: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class PrimerPair:
    name: str
    forward: Primer
    reverse: Primer
    product_size_chm13: int = 0  # designed product size on the CHM13 template


@dataclass
class BindingSite:
    """Where a primer anneals in a haplotype.

    `mismatch_offsets_3p` lists mismatch distances from the primer 3' end (0 = the very
    3' base). This is the field the dropout model cares about most: a mismatch at a small
    offset is far more likely to abolish extension than one deep in the 5' half.
    """

    primer_name: str
    haplotype_id: str
    chrom: str
    start: int  # 0-based, on the haplotype
    end: int
    strand: Strand
    mismatches: int
    mismatch_offsets_3p: list[int] = field(default_factory=list)
    has_indel: bool = False


@dataclass
class Amplicon:
    """A predicted product from a convergent forward/reverse binding-site pair."""

    haplotype_id: str
    chrom: str
    start: int
    end: int
    size: int
    forward: BindingSite
    reverse: BindingSite
    on_target: bool


@dataclass
class HaplotypeResult:
    """Outcome of one primer pair against one haplotype."""

    haplotype_id: str
    status: HaplotypeStatus
    amplicons: list[Amplicon] = field(default_factory=list)
    reason: str = ""  # human-readable "why", esp. for dropout / uncertain


@dataclass
class PairResult:
    """Aggregated outcome of one primer pair across all evaluated haplotypes."""

    pair: PrimerPair
    per_haplotype: list[HaplotypeResult] = field(default_factory=list)
    primer3_penalty: float = 0.0

    @property
    def evaluable(self) -> list[HaplotypeResult]:
        return [
            r for r in self.per_haplotype if r.status is not HaplotypeStatus.UNCERTAIN
        ]

    @property
    def n_uncertain(self) -> int:
        return sum(
            1 for r in self.per_haplotype if r.status is HaplotypeStatus.UNCERTAIN
        )

    def _count(self, status: HaplotypeStatus) -> int:
        return sum(1 for r in self.per_haplotype if r.status is status)

    @property
    def on_target_coverage(self) -> float:
        """Fraction of evaluable haplotypes that produce the intended target amplicon *at
        all* — i.e. have >=1 on-target amplicon. This credits `multi_product` haplotypes
        (target amplifies, but with extra bands), so it answers "does the target amplify
        across diversity?" and is distinct from `unique_product_rate`. It must NOT be
        `pass`-only, or a locus that always throws an extra band reads as 0% coverage."""
        ev = self.evaluable
        if not ev:
            return 0.0
        covered = sum(1 for r in ev if any(a.on_target for a in r.amplicons))
        return covered / len(ev)

    @property
    def unique_product_rate(self) -> float:
        """Fraction of evaluable haplotypes with exactly one product that is the on-target
        one and no off-target — precisely the `pass` definition. The specificity metric."""
        n = len(self.evaluable)
        return self._count(HaplotypeStatus.PASS) / n if n else 0.0

    @property
    def n_dropout(self) -> int:
        return self._count(HaplotypeStatus.DROPOUT)

    @property
    def n_off_target(self) -> int:
        return self._count(HaplotypeStatus.OFF_TARGET)

    @property
    def n_multi_product(self) -> int:
        return self._count(HaplotypeStatus.MULTI_PRODUCT)
