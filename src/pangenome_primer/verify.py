"""Verify mode: screen user-supplied primer pairs across the pangenome (Decision 1D).

Given a CSV of (primer_id, target region, forward, reverse), report per haplotype whether the
pair makes the intended product, its size, and any off-target products — reusing the design
engine's binding/competence/pairing, but skipping Primer3 and masking entirely.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field

from . import config as cfgmod
from .classify import pair_tm
from .engine import EvalConfig, evaluate_with_sites
from .model import Primer, PrimerPair
from .project import (
    make_aligner,
    needs_whole_genome_aligner,
    project_target,
    resolve_target,
)
from .samples import load_haplotypes

# target column aliases accepted in the CSV header (case-insensitive)
_ALIASES = {
    "primer_id": {"primer_id", "id", "name", "pair_id"},
    "target": {"target", "coords", "coordinates", "grch38", "region"},
    "forward": {"forward", "fwd", "forward_primer", "f", "left"},
    "reverse": {"reverse", "rev", "reverse_primer", "r", "right"},
}


@dataclass
class PrimerSpec:
    primer_id: str
    target: str      # chr:start-end (GRCh38 by default)
    forward: str
    reverse: str


@dataclass
class VerifyCell:
    haplotype_id: str
    status: str                                   # pass/dropout/off_target/multi_product/uncertain
    on_target: list[int] = field(default_factory=list)   # correct amplicon sizes
    off_target: list[int] = field(default_factory=list)  # off-target amplicon sizes
    size_flag: bool = False                       # on-target size deviates from expected
    reason: str = ""
    site_cap: int | None = None                   # set when the binding-site cap tripped;
                                                  # cell renders ">N binding sites"


@dataclass
class VerifySummary:
    """Per-row counts across haplotypes — the row's verdict, independent of column count.

    At 3 haplotypes the matrix is read by eye; at 30 it is not, and the number a reader
    actually wants ("does this pair work across diversity?") is not recoverable by scanning
    30 cells. The definitions mirror `model.PairResult` on purpose, so verify and design
    report the same quantity under the same name:

    - `coverage` is `on_target_coverage`: of the **evaluable** haplotypes, the fraction with
      at least one on-target amplicon. It credits `multi_product` deliberately — a pass-only
      coverage makes a locus that always throws an extra band read as 0%, the mislabel
      `model.py` warns about at `on_target_coverage`.
    - `unique_rate` is `unique_product_rate`: pass-only, the specificity counterpart.

    `evaluable` excludes `uncertain` (locus not projectable), matching the rule that every
    rank gate is applied over total − uncertain.
    """

    total: int
    evaluable: int
    n_covered: int         # evaluable cells with >=1 on-target amplicon (coverage numerator)
    n_pass: int
    n_dropout: int
    n_off_target: int
    n_multi_product: int
    n_uncertain: int
    n_capped: int          # cells whose binding-site cap tripped; see note below
    coverage: float
    unique_rate: float


def summarize(cells: list[VerifyCell]) -> VerifySummary:
    """Derive a row summary from its cells. Always computed, never stored alongside the
    cells it describes — a summary that can drift from its matrix is worse than none."""
    n = len(cells)
    n_unc = sum(1 for c in cells if c.status == "uncertain")
    ev = [c for c in cells if c.status != "uncertain"]
    n_cov = sum(1 for c in ev if c.on_target)
    # A capped cell enumerates no products by design, so it counts as not-covered here.
    # Reported separately because that is an *unknown*, not a measured failure to amplify:
    # without it a promiscuous primer's row is indistinguishable from a genuinely dead one.
    return VerifySummary(
        total=n,
        evaluable=len(ev),
        n_covered=n_cov,
        n_pass=sum(1 for c in cells if c.status == "pass"),
        n_dropout=sum(1 for c in cells if c.status == "dropout"),
        n_off_target=sum(1 for c in cells if c.status == "off_target"),
        n_multi_product=sum(1 for c in cells if c.status == "multi_product"),
        n_uncertain=n_unc,
        n_capped=sum(1 for c in cells if c.site_cap),
        coverage=(n_cov / len(ev)) if ev else 0.0,
        unique_rate=(sum(1 for c in ev if c.status == "pass") / len(ev)) if ev else 0.0,
    )


@dataclass
class VerifyRow:
    primer_id: str
    forward: str
    reverse: str
    target_input: str
    target_chm13: str
    expected_size: int
    cells: list[VerifyCell] = field(default_factory=list)

    @property
    def summary(self) -> VerifySummary:
        return summarize(self.cells)


def _resolve_headers(fieldnames: list[str]) -> dict[str, str]:
    lut = {}
    lowered = {fn.strip().lower(): fn for fn in fieldnames}
    for canon, aliases in _ALIASES.items():
        match = next((lowered[a] for a in aliases if a in lowered), None)
        if match is None:
            raise ValueError(
                f"primers CSV missing a '{canon}' column (accepted: {sorted(aliases)})"
            )
        lut[canon] = match
    return lut


def parse_primer_csv(path: str) -> list[PrimerSpec]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        cols = _resolve_headers(reader.fieldnames or [])
        specs = []
        for rec in reader:
            if not rec.get(cols["primer_id"], "").strip():
                continue
            specs.append(PrimerSpec(
                rec[cols["primer_id"]].strip(),
                rec[cols["target"]].strip(),
                rec[cols["forward"]].strip().upper(),
                rec[cols["reverse"]].strip().upper(),
            ))
    return specs


def run_verify(
    specs: list[PrimerSpec],
    chm13_fasta: str,
    samples_tsv: str,
    *,
    grch38_fasta: str | None = None,
    assembly: str = "grch38",
    max_amplicon: int = 2000,
    size_tolerance: int = 20,
    mode: str | None = None,
    config_path: str | None = None,
    pad: int = 250,
    progress=lambda s: None,
) -> list[VerifyRow]:
    import pysam

    from .search import backend_from_config, find_binding_sites_batch, threads_from_config

    raw = cfgmod.load_raw(config_path)
    mode = mode or raw["dropout"]["mode"]
    max_mm = raw["search"]["max_mismatches"]
    # Absent from an older config file => fall back to the EvalConfig default rather than
    # KeyError, so a user's pinned defaults.yaml keeps working.
    max_sites = raw["search"].get("max_binding_sites", EvalConfig.max_binding_sites)
    backend = backend_from_config(raw, warn=progress)
    threads = threads_from_config(raw)
    tm_opt = cfgmod.design_config(raw).tm_opt

    haplos = load_haplotypes(samples_tsv)
    all_seqs = [s.forward for s in specs] + [s.reverse for s in specs]

    # Resolve each pair's target + build its eval config ONCE (independent of haplotype).
    chm = pysam.FastaFile(chm13_fasta)
    pctx = []  # (spec, chrom, start, end, tstart, tend, expected_size, template, pair, ecfg)
    for spec in specs:
        chrom, start, end = resolve_target(
            spec.target, chm13_fasta, assembly=assembly, grch38_fasta=grch38_fasta)
        expected_size = end - start
        clen = chm.get_reference_length(chrom)
        tstart, tend = max(0, start - pad), min(clen, end + pad)
        template = chm.fetch(chrom, tstart, tend)
        max_prod = max(max_amplicon, expected_size + pad)  # never filter out a larger correct one
        pair = PrimerPair(spec.primer_id, Primer(f"{spec.primer_id}_F", spec.forward),
                          Primer(f"{spec.primer_id}_R", spec.reverse),
                          product_size_chm13=expected_size)
        ecfg = EvalConfig(
            mode=mode, max_mismatches=max_mm, min_product=40, max_product=max_prod,
            max_binding_sites=max_sites,
            rule_cfg=cfgmod.rule_config(raw),
            thermo_cfg=cfgmod.thermo_config(raw, pair_tm(spec.forward, spec.reverse, tm_opt)))
        pctx.append((spec, chrom, start, end, tstart, tend, expected_size, template, pair, ecfg))
    chm.close()

    # Haplotype-outer: one genome-wide search AND one projection-index load per haplotype,
    # reused across all pairs; the index is released before the next haplotype (bounds RAM).
    # One FASTA handle per haplotype is opened here and threaded into both the batch search
    # and the per-pair PAF-lift projection, instead of each reopening its own. (The rust
    # backend reads the BGZF file itself and ignores the handle; projection uses it.)
    cells: dict[tuple[str, str], VerifyCell] = {}
    for hid, fasta in haplos:
        fa = pysam.FastaFile(fasta)
        progress(f"genome-wide search ({backend}) on {hid} ...")
        sites = find_binding_sites_batch(all_seqs, fasta, hid, max_mm, fa=fa,
                                         backend=backend, threads=threads)
        aligner = make_aligner(fasta) if needs_whole_genome_aligner(fasta) else None
        for (spec, chrom, start, end, tstart, tend, expected_size, template, pair, ecfg) in pctx:
            proj = project_target(chrom, tstart, tend, template, fasta, aligner=aligner, fa=fa)
            if proj.locus is None:
                cells[(spec.primer_id, hid)] = VerifyCell(hid, "uncertain", reason=proj.reason)
                continue
            res = evaluate_with_sites(
                pair, hid, sites.get(spec.forward, []), sites.get(spec.reverse, []),
                proj.locus, lambda c, s, e, _fa=fa: _fa.fetch(c, s, e), ecfg)
            on = sorted(a.size for a in res.amplicons if a.on_target)
            off = sorted(a.size for a in res.amplicons if not a.on_target)
            flag = any(abs(sz - expected_size) > size_tolerance for sz in on)
            # The engine reports a tripped binding-site cap in `reason`; surface it as a
            # first-class cell field so the matrix can render ">N binding sites" rather
            # than an empty multi_product cell.
            capped = ecfg.max_binding_sites if "binding sites" in res.reason else None
            cells[(spec.primer_id, hid)] = VerifyCell(
                hid, res.status.value, on, off, flag, res.reason, site_cap=capped
            )
        fa.close()
        aligner = None  # release the projection index before the next haplotype

    rows = []
    for (spec, chrom, start, end, tstart, tend, expected_size, template, pair, ecfg) in pctx:
        rows.append(VerifyRow(
            spec.primer_id, spec.forward, spec.reverse, spec.target,
            f"{chrom}:{start}-{end}", expected_size,
            cells=[cells[(spec.primer_id, hid)] for hid, _ in haplos]))
    return rows
