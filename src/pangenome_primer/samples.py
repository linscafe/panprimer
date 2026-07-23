"""Shared haplotype-manifest loading. Columns 0-3 of samples.tsv are
sample / hap / superpop / local_path; comment (#) and header rows are skipped."""
from __future__ import annotations

from pathlib import Path


def load_haplotypes(samples_tsv: str) -> list[tuple[str, str]]:
    """Return (haplotype_id, fasta_path) for each row; the FASTA must exist locally."""
    rows: list[tuple[str, str]] = []
    for line in Path(samples_tsv).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("sample\t"):
            continue
        f = line.split("\t")
        sample, hap, _superpop, path = f[0], f[1], f[2], f[3]
        hid = f"{sample}#hap{hap}"
        if not Path(path).exists():
            raise FileNotFoundError(
                f"haplotype FASTA not found for {hid}: {path}\n"
                "Fetch/prepare the subset first (see config/samples.tsv)."
            )
        rows.append((hid, path))
    return rows
