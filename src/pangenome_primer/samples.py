"""Shared haplotype-manifest loading. Columns 0-3 of samples.tsv are
sample / hap / superpop / local_path; comment (#) and header rows are skipped."""
from __future__ import annotations

from pathlib import Path


def sidecar_path(seq_path: str, suffix: str) -> str:
    """Locate a sidecar file (`.mmi`, `.chm13.paf`, ...) for a haplotype sequence path.

    Phase 3 repoints `samples.tsv` `local_path` from `X.fa` to the BGZF `X.fa.gz`, but the
    projection sidecars are named after the *uncompressed* stem (`X.fa.mmi`, built by
    `scripts/prepare_haplotypes.sh`). Deriving them by plain string append -- as this code
    used to -- would look for `X.fa.gz.mmi`, miss, and fall through to building a full
    minimap2 index in-process from the gzipped FASTA on every haplotype. That does not fail;
    it just costs minutes per haplotype, which is exactly the kind of regression that reads
    as "Phase 3 works, only slower".

    Resolution order:

    1. `<seq_path><suffix>` when it exists -- an exact-named sidecar always wins, so a
       genuinely `.fa.gz`-derived index is never shadowed;
    2. the `.gz`-stripped stem + `suffix` when that exists -- the `.fa`-named sidecar sitting
       beside the `.fa.gz`;
    3. otherwise `<seq_path><suffix>` unchanged, so *writers* still create the sidecar next
       to the path they were handed.
    """
    exact = seq_path + suffix
    if Path(exact).exists():
        return exact
    if seq_path.endswith(".gz"):
        stem = seq_path[:-3] + suffix
        if Path(stem).exists():
            return stem
    return exact


def load_haplotypes(samples_tsv: str) -> list[tuple[str, str]]:
    """Return (haplotype_id, fasta_path) for each row; the FASTA must exist locally."""
    rows: list[tuple[str, str]] = []
    for line in Path(samples_tsv).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("sample\t"):
            continue
        f = line.split("\t")  # columns: sample, hap, superpop, local_path, ...
        sample, hap, path = f[0], f[1], f[3]
        hid = f"{sample}#hap{hap}"
        if not Path(path).exists():
            raise FileNotFoundError(
                f"haplotype FASTA not found for {hid}: {path}\n"
                "Fetch/prepare the subset first (see config/samples.tsv)."
            )
        rows.append((hid, path))
    return rows
