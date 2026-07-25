"""Fetch an HPRC Release 2 haplotype subset and write the pipeline's samples.tsv.

Authoritative machine-readable sources (human-pangenomics/hprc_intermediate_assembly):
  * assemblies_release2_v1.0.index.csv  — 465 haplotypes, S3 URLs + md5 + genbank accession
  * hprc_release2_sample_metadata.csv   — sample_id -> population_abbreviation

Assemblies live in the public, no-egress-fee `human-pangenomics` S3 bucket and are served
over plain HTTPS, so no aws CLI or credentials are needed. The multi-GB assembly download is
gated behind an explicit flag — by default we only build the manifest (subset-first).
"""
from __future__ import annotations

import csv
import hashlib
import io
import os
import shutil
import urllib.request
from dataclasses import dataclass
from pathlib import Path

# the default subset: hap1 of these three samples (AFR/EUR/EAS) — identical to the demo
# (demo-design-pipeline/samples.tsv)
DEMO_SAMPLES = ("HG01884", "HG00097", "HG00408")

# rough per-haplotype footprints (see the README data-setup warning) used for the resource
# check before a real download
DOWNLOAD_GB_PER_HAP = 1.0     # compressed .fa.gz pulled now (0.90 GB measured)
# After `prepare_haplotypes.sh`: the BGZF stays (0.90 GB), its .fai/.gzi cost <1 MB, and the
# sparse projection anchor grid adds ~4.3 MB. The 3.08 GB uncompressed .fa and the 5.31 GB
# bwa index are opt-in via WITH_BWA=1 (which returns this to ~15 GB); the 5.80 GB minimap2
# .mmi the grid replaces is no longer built at all.
INDEXED_GB_PER_HAP = 0.91
CHM13_GB = 8.5                # one-time reference + its indexes
PIPELINE_PEAK_GB = 12.0       # peak RAM: one genome index in memory at a time
RAM_HEADROOM_GB = 3.0

_RAW = "https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables"
ASSEMBLY_INDEX_URL = f"{_RAW}/assemblies_release2_v1.0.index.csv"
SAMPLE_META_URL = f"{_RAW}/sample/hprc_release2_sample_metadata.csv"
RELEASE = "R2"

# 1000 Genomes population code -> superpopulation. HPRC-specific codes (e.g. ASL) mapped by
# the population they represent; anything unknown falls back to "UNK".
POP_TO_SUPERPOP = {
    # AFR
    "YRI": "AFR", "LWK": "AFR", "GWD": "AFR", "MSL": "AFR", "ESN": "AFR",
    "ASW": "AFR", "ACB": "AFR", "ASL": "AFR", "MKK": "AFR",
    # EUR
    "CEU": "EUR", "TSI": "EUR", "FIN": "EUR", "GBR": "EUR", "IBS": "EUR",
    # EAS
    "CHB": "EAS", "JPT": "EAS", "CHS": "EAS", "CDX": "EAS", "KHV": "EAS",
    # SAS
    "GIH": "SAS", "PJL": "SAS", "BEB": "SAS", "STU": "SAS", "ITU": "SAS",
    # AMR
    "MXL": "AMR", "PUR": "AMR", "CLM": "AMR", "PEL": "AMR",
}


@dataclass
class HaplotypeRecord:
    sample_id: str
    hap: int
    superpop: str
    population: str
    assembly_url: str  # https
    md5_url: str       # https
    fai_url: str       # https
    genbank: str

    @property
    def hap_id(self) -> str:
        return f"{self.sample_id}#hap{self.hap}"

    def local_path(self, data_dir: str) -> str:
        return str(Path(data_dir) / f"{self.sample_id}_hap{self.hap}.fa.gz")


def s3_to_https(url: str) -> str:
    if url.startswith("s3://"):
        bucket, _, key = url[5:].partition("/")
        return f"https://{bucket}.s3.amazonaws.com/{key}"
    return url


def _read_csv(url_or_path: str) -> list[dict]:
    if str(url_or_path).startswith("http"):
        with urllib.request.urlopen(url_or_path) as r:  # noqa: S310 (trusted host)
            text = r.read().decode()
    else:
        text = Path(url_or_path).read_text()
    return list(csv.DictReader(io.StringIO(text)))


def load_records(
    index_url: str = ASSEMBLY_INDEX_URL, meta_url: str = SAMPLE_META_URL
) -> list[HaplotypeRecord]:
    """Join the assembly index with sample metadata into annotated haplotype records."""
    meta = {r["sample_id"]: r for r in _read_csv(meta_url)}
    records: list[HaplotypeRecord] = []
    for row in _read_csv(index_url):
        sid = row["sample_id"]
        pop = (meta.get(sid, {}).get("population_abbreviation") or "").strip()
        records.append(
            HaplotypeRecord(
                sample_id=sid,
                hap=int(row["haplotype"]),
                superpop=POP_TO_SUPERPOP.get(pop, "UNK"),
                population=pop or "NA",
                assembly_url=s3_to_https(row["assembly"]),
                md5_url=s3_to_https(row["assembly_md5"]),
                fai_url=s3_to_https(row["assembly_fai"]),
                genbank=row.get("genbank_accession", ""),
            )
        )
    return records


def select_diverse(
    records: list[HaplotypeRecord],
    per_superpop: int = 1,
    superpops: tuple[str, ...] = ("AFR", "EUR", "EAS", "SAS", "AMR"),
    haps: tuple[int, ...] | None = None,
) -> list[HaplotypeRecord]:
    """Pick `per_superpop` samples from each of `superpops` for a subset that tests
    cross-population universality. Both haplotypes by default; pass `haps=(1,)` to restrict
    to hap1 (halves the footprint)."""
    by_sample: dict[str, list[HaplotypeRecord]] = {}
    for r in records:
        by_sample.setdefault(r.sample_id, []).append(r)
    chosen: list[HaplotypeRecord] = []
    for superpop in superpops:
        samples = sorted(
            {s for s, rs in by_sample.items() if rs[0].superpop == superpop}
        )
        for sid in samples[:per_superpop]:
            recs = sorted(by_sample[sid], key=lambda r: r.hap)
            if haps is not None:
                recs = [r for r in recs if r.hap in haps]
            chosen.extend(recs)
    return chosen


def select_default(records: list[HaplotypeRecord]) -> list[HaplotypeRecord]:
    """The default subset: hap1 of the three demo samples — HG01884 (AFR), HG00097 (EUR),
    HG00408 (EAS). Identical to demo-design-pipeline/samples.tsv."""
    return [r for r in select_samples(records, list(DEMO_SAMPLES)) if r.hap == 1]


def select_all(records: list[HaplotypeRecord]) -> list[HaplotypeRecord]:
    """Every R2 haplotype (~460), sorted. Huge — gate behind the resource check."""
    return sorted(records, key=lambda r: (r.sample_id, r.hap))


def select_samples(
    records: list[HaplotypeRecord], sample_ids: list[str]
) -> list[HaplotypeRecord]:
    wanted = set(sample_ids)
    return sorted(
        [r for r in records if r.sample_id in wanted], key=lambda r: (r.sample_id, r.hap)
    )


# --- resource check (before a real download) --------------------------------

@dataclass
class ResourceCheck:
    """Estimated footprint of a subset vs the host's free disk and RAM."""

    n_haps: int
    download_gb: float   # compressed pull now
    indexed_gb: float    # after gunzip + bwa + minimap2 indexing
    free_disk_gb: float
    total_ram_gb: float  # 0.0 => could not determine

    @property
    def disk_ok(self) -> bool:
        return self.free_disk_gb >= self.indexed_gb + CHM13_GB

    @property
    def ram_ok(self) -> bool:
        # can't determine => don't block on it
        return self.total_ram_gb == 0.0 or self.total_ram_gb >= PIPELINE_PEAK_GB + RAM_HEADROOM_GB

    @property
    def risky(self) -> bool:
        return not self.disk_ok or not self.ram_ok


def free_disk_gb(path: str) -> float:
    """Free space (GB) on the filesystem holding `path`, walking up to the nearest existing
    ancestor (the data dir may not exist yet)."""
    p = Path(path)
    while not p.exists() and p != p.parent:
        p = p.parent
    try:
        return shutil.disk_usage(p).free / 1e9
    except OSError:
        return 0.0


def total_ram_gb() -> float:
    """Total system RAM (GB), or 0.0 if it can't be determined on this platform."""
    try:
        return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 1e9
    except (ValueError, OSError, AttributeError):
        return 0.0


def assess_resources(n_haps: int, data_dir: str) -> ResourceCheck:
    return ResourceCheck(
        n_haps=n_haps,
        download_gb=n_haps * DOWNLOAD_GB_PER_HAP,
        indexed_gb=n_haps * INDEXED_GB_PER_HAP,
        free_disk_gb=free_disk_gb(data_dir),
        total_ram_gb=total_ram_gb(),
    )


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as r, open(tmp, "wb") as fh:  # noqa: S310
        while chunk := r.read(1 << 20):
            fh.write(chunk)
    tmp.replace(dest)


def _md5(path: Path) -> str:
    h = hashlib.md5()  # noqa: S324 (matching HPRC-published md5, not for security)
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def fetch_md5(rec: HaplotypeRecord) -> str:
    """Read the published .md5 (format: 'hash  filename')."""
    with urllib.request.urlopen(rec.md5_url) as r:  # noqa: S310
        return r.read().decode().split()[0]


def download_haplotype(rec: HaplotypeRecord, data_dir: str, verify: bool = True) -> str:
    """Download the assembly (+ .fai), verify md5, return the local path."""
    dest = Path(rec.local_path(data_dir))
    expected = fetch_md5(rec)
    if not (dest.exists() and _md5(dest) == expected):
        _download(rec.assembly_url, dest)
        got = _md5(dest)
        if verify and got != expected:
            dest.unlink(missing_ok=True)
            raise RuntimeError(f"md5 mismatch for {rec.hap_id}: got {got} != {expected}")
    _download(rec.fai_url, Path(str(dest) + ".fai"))
    return str(dest)


SAMPLES_HEADER = [
    "sample", "hap", "superpop", "local_path",
    "md5", "release", "population", "source_url", "genbank",
]


def write_samples_tsv(
    records: list[HaplotypeRecord], out_path: str, data_dir: str, md5s: dict[str, str]
) -> None:
    """Write the pipeline manifest. Columns 0-3 (sample/hap/superpop/local_path) are what
    the pipeline reads; the rest pin provenance (md5, release, population, URL, accession)."""
    lines = [
        "# HPRC R2 subset — provenance pinned. local_path is populated by "
        "`pangenome-primer fetch-subset --download`.",
        "\t".join(SAMPLES_HEADER),
    ]
    for r in records:
        lines.append("\t".join([
            r.sample_id, str(r.hap), r.superpop, r.local_path(data_dir),
            md5s.get(r.hap_id, "PENDING"), RELEASE, r.population,
            r.assembly_url, r.genbank,
        ]))
    Path(out_path).write_text("\n".join(lines) + "\n")
