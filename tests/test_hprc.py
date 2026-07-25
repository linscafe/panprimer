"""Unit tests for the HPRC subset selection/join logic. No network: the index and metadata
are monkeypatched with tiny inline CSVs so the join, superpop mapping, and diverse-pick are
deterministic."""
from __future__ import annotations

import io

from pangenome_primer import hprc

_INDEX = """sample_id,haplotype,assembly,assembly_md5,assembly_fai,genbank_accession
S_AFR,1,s3://human-pangenomics/x/S_AFR_pat.fa.gz,s3://b/S_AFR.md5,s3://b/S_AFR.fai,GCA_1
S_AFR,2,s3://human-pangenomics/x/S_AFR_mat.fa.gz,s3://b/S_AFR2.md5,s3://b/S_AFR2.fai,GCA_2
S_EUR,1,s3://human-pangenomics/x/S_EUR_h1.fa.gz,s3://b/S_EUR.md5,s3://b/S_EUR.fai,GCA_3
S_EUR,2,s3://human-pangenomics/x/S_EUR_h2.fa.gz,s3://b/S_EUR2.md5,s3://b/S_EUR2.fai,GCA_4
S_UNK,1,s3://human-pangenomics/x/S_UNK.fa.gz,s3://b/S_UNK.md5,s3://b/S_UNK.fai,GCA_5
"""
_META = """sample_id,population_abbreviation
S_AFR,YRI
S_EUR,GBR
S_UNK,ZZZ
"""


def _patch(monkeypatch):
    def fake(url):
        import csv

        text = _INDEX if "index" in url else _META
        return list(csv.DictReader(io.StringIO(text)))

    monkeypatch.setattr(hprc, "_read_csv", fake)


def test_s3_to_https():
    assert hprc.s3_to_https("s3://human-pangenomics/a/b.fa.gz") == \
        "https://human-pangenomics.s3.amazonaws.com/a/b.fa.gz"


def test_join_assigns_superpop(monkeypatch):
    _patch(monkeypatch)
    recs = hprc.load_records("index", "meta")
    by = {r.hap_id: r for r in recs}
    assert by["S_AFR#hap1"].superpop == "AFR"
    assert by["S_EUR#hap1"].superpop == "EUR"
    assert by["S_UNK#hap1"].superpop == "UNK"  # unmapped population
    assert by["S_AFR#hap1"].assembly_url.startswith("https://human-pangenomics.s3")


def test_select_diverse_picks_per_superpop(monkeypatch):
    _patch(monkeypatch)
    recs = hprc.load_records("index", "meta")
    chosen = hprc.select_diverse(recs, per_superpop=1)
    ids = {r.sample_id for r in chosen}
    assert ids == {"S_AFR", "S_EUR"}          # UNK excluded; AMR/EAS/SAS absent here
    assert len(chosen) == 4                    # both haplotypes of each


def test_select_samples_explicit(monkeypatch):
    _patch(monkeypatch)
    recs = hprc.load_records("index", "meta")
    chosen = hprc.select_samples(recs, ["S_UNK"])
    assert [r.hap for r in chosen] == [1]
    assert chosen[0].sample_id == "S_UNK"


def test_select_default_is_the_three_demo_hap1():
    # default is exactly hap1 of the three demo samples (a decoy sample is ignored)
    recs = [
        hprc.HaplotypeRecord(sid, hap, "NA", "NA", f"u/{sid}_{hap}.fa.gz", "m", "f", "g")
        for sid in ("HG01884", "HG00097", "HG00408", "HG99999")
        for hap in (1, 2)
    ]
    chosen = hprc.select_default(recs)
    assert [(r.sample_id, r.hap) for r in chosen] == [
        ("HG00097", 1), ("HG00408", 1), ("HG01884", 1)]  # sorted; hap1 only; decoy excluded


def test_select_all_returns_everything(monkeypatch):
    _patch(monkeypatch)
    recs = hprc.load_records("index", "meta")
    chosen = hprc.select_all(recs)
    assert len(chosen) == len(recs)                       # includes UNK + both haplotypes
    assert [(r.sample_id, r.hap) for r in chosen] == sorted((r.sample_id, r.hap) for r in recs)


def test_assess_resources_ok(monkeypatch, tmp_path):
    monkeypatch.setattr(hprc, "free_disk_gb", lambda p: 100.0)
    monkeypatch.setattr(hprc, "total_ram_gb", lambda: 16.0)
    rc = hprc.assess_resources(3, str(tmp_path))
    # Derived from the constants rather than hardcoded: the per-haplotype footprint is
    # expected to fall as storage phases land (15.0 -> 6.7 with BGZF-only -> ~0.91 once the
    # anchor grid replaces the .mmi), and this test is about the per-hap arithmetic, not the
    # current figure.
    assert (rc.download_gb, rc.indexed_gb) == (
        3 * hprc.DOWNLOAD_GB_PER_HAP,
        3 * hprc.INDEXED_GB_PER_HAP,
    )
    assert rc.disk_ok and rc.ram_ok and not rc.risky


def test_assess_resources_flags_shortfalls(monkeypatch, tmp_path):
    # Half of what 3 haplotypes + CHM13 need. Derived, not hardcoded: a fixed 20.0 GB stops
    # being a shortfall once the per-haplotype footprint drops far enough, and the test would
    # then pass for the wrong reason -- silently no longer exercising the failure path.
    short = (3 * hprc.INDEXED_GB_PER_HAP + hprc.CHM13_GB) / 2
    monkeypatch.setattr(hprc, "free_disk_gb", lambda p: short)
    monkeypatch.setattr(hprc, "total_ram_gb", lambda: 8.0)     # < 12 + 3 headroom
    rc = hprc.assess_resources(3, str(tmp_path))
    assert not rc.disk_ok and not rc.ram_ok and rc.risky


def test_assess_resources_unknown_ram_not_blocking(monkeypatch, tmp_path):
    monkeypatch.setattr(hprc, "free_disk_gb", lambda p: 100.0)
    monkeypatch.setattr(hprc, "total_ram_gb", lambda: 0.0)     # undeterminable
    rc = hprc.assess_resources(1, str(tmp_path))
    assert rc.ram_ok and not rc.risky
