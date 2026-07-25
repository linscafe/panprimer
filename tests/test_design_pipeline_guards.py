"""The design pipeline must honour `search:` config the same way verify does.

Both pipelines share `search.find_binding_sites_batch` and `EvalConfig`, so most of the
storage/search rework reached design for free. `search.max_binding_sites` did not: it was
read from config by `verify.run_verify` but not by `cli`, so design silently used the
dataclass default and ignored the config file -- a key that worked in one pipeline and not
the other.

Exercised through the `evaluate` subcommand, the per-haplotype unit the Nextflow path calls,
which needs no CHM13 anchoring to drive.
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


def _config_with(tmp_path, pattern, replacement, name="cfg.yaml"):
    """A full copy of config/defaults.yaml with one key rewritten.

    `load_raw` does not merge a partial file over the defaults, so a two-line test config
    would KeyError on everything it omits. Start from the real defaults and override.
    """
    import re

    cfg = tmp_path / name
    cfg.write_text(re.sub(pattern, replacement, Path("config/defaults.yaml").read_text(),
                          count=1, flags=re.M))
    return cfg


class TestSearchThreadsIsReadFromConfig:
    """`search.threads` caps the compiled scanner's rayon pool (ISSUE-002). It is the knob
    that keeps a concurrent fan-out from spawning cores x tasks threads, so a version that
    parses cleanly and is then dropped on the floor is worse than no knob at all -- exactly
    how `max_binding_sites` failed below."""

    # The scanner's pool is per-process and built by whichever test scanned first, so asking
    # for 3 here legitimately warns that it will not be honoured. That warning is the subject
    # of its own test in tests/test_scan_stack_depth.py; what *this* test checks is the
    # config -> call-site wiring, which is upstream of the pool entirely.
    @pytest.mark.filterwarnings("ignore:search.threads=:RuntimeWarning")
    def test_config_value_reaches_the_scanner(self, bundle, tmp_path):
        cfg = _config_with(tmp_path, r"^  threads: .*$", "  threads: 3")

        from pangenome_primer import search

        captured = {}
        real = search.find_binding_sites_batch

        def spy(*a, **kw):
            captured["threads"] = kw.get("threads")
            return real(*a, **kw)

        search.find_binding_sites_batch = spy
        try:
            res, _ = _run(bundle, ["--config", str(cfg)])
        finally:
            search.find_binding_sites_batch = real
        assert res.exit_code == 0, res.output
        assert captured.get("threads") == 3, (
            f"design passed threads={captured.get('threads')!r}; config said 3"
        )

    def test_default_is_none_not_a_fabricated_number(self, bundle, tmp_path):
        """Shipped default is null = "let the scanner choose", which on the sequential CLI
        loop means all cores. A test that pinned a number here would lock in the wrong
        behaviour for a solo scan."""
        from pangenome_primer import search

        assert search.threads_from_config({"search": {}}) is None
        assert search.threads_from_config({"search": {"threads": None}}) is None


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
