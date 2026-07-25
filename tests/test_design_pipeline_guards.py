"""Guards on the design pipeline that mirror ones the verify pipeline already had.

The two pipelines share `search.find_binding_sites_batch` and `EvalConfig`, so most of the
storage/search rework reached design for free. Two things did not, because they live in the
caller rather than the seam:

* `search.max_binding_sites` was read from config by `verify.run_verify` but not by `cli`,
  so the design pipeline silently used the dataclass default and ignored the config file.
* a truncated bwa hit list was refused by verify but scored by design. That one is worse
  than it sounds: truncation leaves ONE site, which is far *below* the binding-site cap, so
  the pair looks cleanly unique and ranks near the top. The cap cannot catch it.

Both are exercised through the `evaluate` subcommand, which is the per-haplotype unit the
Nextflow path calls and needs no CHM13 anchoring to drive.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner

from pangenome_primer import cli

FWD = "TTGCACAGTCCAGATTGCAA"
REV = "CATTGCGATTGACACTTGCG"


def _write_bgzf(tmp_path: Path, seq: str) -> str:
    """A one-contig BGZF FASTA with .fai/.gzi, the layout the pipeline reads."""
    plain = tmp_path / "hap.fa"
    with open(plain, "w") as fh:
        fh.write(">ctg1\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    gz = tmp_path / "hap.fa.gz"
    subprocess.run(["bgzip", "-f", "-c", str(plain)], stdout=open(gz, "wb"), check=True)
    subprocess.run(["samtools", "faidx", str(gz)], check=True)
    return str(gz)


@pytest.fixture
def bundle(tmp_path):
    """candidates.json + projection.json + a real BGZF haplotype."""
    import random

    rng = random.Random(5)
    body = "".join(rng.choice("ACGT") for _ in range(4000))
    seq = body[:1000] + FWD + body[1000:1200] + REV + body[1200:]
    fasta = _write_bgzf(tmp_path, seq)

    cand = {
        "pair": {
            "name": "p1",
            "forward": {"name": "p1_F", "sequence": FWD},
            "reverse": {"name": "p1_R", "sequence": REV},
            "product_size_chm13": 240,
        },
        "left_start": 0, "left_len": 20, "right_start": 220, "right_len": 20,
        "penalty": 0.5,
    }
    cj = tmp_path / "cands.json"
    cj.write_text(json.dumps({"candidates": [cand]}))
    pj = tmp_path / "proj.json"
    pj.write_text(json.dumps({
        "hap_id": "H#hap1",
        "fasta": fasta,
        "locus": {"assembly": fasta, "chrom": "ctg1", "start": 900, "end": 1400},
    }))
    return {"cands": str(cj), "proj": str(pj), "fasta": fasta, "tmp": tmp_path}


def _run(bundle, extra=None):
    out = Path(bundle["tmp"]) / "out.json"
    args = ["evaluate", "--candidates", bundle["cands"], "--projection", bundle["proj"],
            "--out", str(out)]
    return CliRunner().invoke(cli.cli, args + (extra or [])), out


class TestTruncationGuard:
    """A truncated hit list must stop the run, not quietly rank the primer."""

    def test_truncation_aborts_with_a_clear_message(self, bundle, monkeypatch):
        from pangenome_primer import search

        def fake(seqs, fasta, hid, max_mm, *, slop=3, fa=None, backend=None, truncated=None):
            if truncated is not None:
                truncated[seqs[0]] = 329_968
            return {s: [] for s in seqs}

        monkeypatch.setattr(search, "find_binding_sites_batch", fake)
        res, _ = _run(bundle)
        assert res.exit_code != 0, "a truncated hit list was scored instead of refused"
        assert "329968" in res.output
        assert "incomplete" in res.output.lower()

    def test_clean_search_is_not_flagged(self, bundle, monkeypatch):
        """The guard must fire on truncation only -- not on every bwa run."""
        from pangenome_primer import search

        def fake(seqs, fasta, hid, max_mm, *, slop=3, fa=None, backend=None, truncated=None):
            return {s: [] for s in seqs}

        monkeypatch.setattr(search, "find_binding_sites_batch", fake)
        res, out = _run(bundle)
        assert res.exit_code == 0, res.output
        assert out.exists()


class TestMaxBindingSitesIsReadFromConfig:
    def test_config_value_reaches_the_design_evaluation(self, bundle, tmp_path):
        """Previously `cli` relied on EvalConfig's default, so this key was honoured in
        verify and silently ignored in design -- config that works in one pipeline and not
        the other."""
        captured = {}
        # `load_raw` does not merge a partial file over the defaults, so start from the real
        # defaults and override the single key under test.
        import re

        defaults = Path("config/defaults.yaml").read_text()
        cfg = tmp_path / "cap.yaml"
        cfg.write_text(re.sub(r"^  max_binding_sites: \d+", "  max_binding_sites: 7",
                              defaults, count=1, flags=re.M))

        from pangenome_primer import engine

        real = engine.evaluate_with_sites

        def spy(pair, hid, f_sites, r_sites, locus, window_fn, ecfg):
            captured["cap"] = ecfg.max_binding_sites
            return real(pair, hid, f_sites, r_sites, locus, window_fn, ecfg)

        import pangenome_primer.cli as clim
        orig = clim.evaluate_with_sites
        clim.evaluate_with_sites = spy
        try:
            res, _ = _run(bundle, ["--config", str(cfg)])
        finally:
            clim.evaluate_with_sites = orig
        assert res.exit_code == 0, res.output
        assert captured.get("cap") == 7, (
            f"design used cap={captured.get('cap')!r}; config/defaults.yaml said 7"
        )
