"""Command-line interface.

* `selftest` — run the synthetic mini-pangenome end to end (no data/network). Verifies the
  engine + reporting path anywhere.
* `run` — the real subset pipeline: anchor the target on CHM13, project onto each
  haplotype, mask variable sites, design with Primer3, evaluate genome-wide per haplotype,
  rank, and report. Needs local CHM13 + haplotype FASTAs (see config/samples.tsv).

The Nextflow workflow shells out to these same subcommands, one per stage, for scale-out.
"""
from __future__ import annotations

import sys
from pathlib import Path

import click

from . import __version__, config as cfgmod, report
from .engine import EvalConfig, evaluate_pair, evaluate_with_sites
from .model import HaplotypeStatus, PairResult
from .rank import rank_pairs


@click.group()
@click.version_option(__version__)
def cli() -> None:
    """Design universal PCR primers against the HPRC R2 pangenome."""


@cli.command()
@click.option("--outdir", default="selftest_out", help="where to write results")
def selftest(outdir: str) -> None:
    """Run the synthetic mini-pangenome (rule mode) and write a report."""
    from . import fixtures as fx
    from .classify import RuleConfig

    cfg = EvalConfig(mode="rule", min_product=80, max_product=200, rule_cfg=RuleConfig())
    result = evaluate_pair(fx.demo_pair(), fx.demo_contexts(), cfg)
    ranked = rank_pairs([result])
    paths = report.write_all(ranked, outdir, provenance={"release": "R2-selftest-demo"})
    counts = {s.value: 0 for s in HaplotypeStatus}
    for hr in result.per_haplotype:
        counts[hr.status.value] += 1
    click.echo(f"selftest statuses: {counts}")
    click.echo(f"coverage={result.on_target_coverage:.0%} over {len(result.evaluable)} evaluable")
    for k, v in paths.items():
        click.echo(f"  {k}: {v}")
    expected = {"pass": 1, "dropout": 1, "off_target": 1, "multi_product": 1, "uncertain": 1}
    if counts != expected:
        click.echo(f"UNEXPECTED: wanted {expected}", err=True)
        sys.exit(1)
    click.echo("selftest OK")


def _load_haplotypes(samples_tsv: str) -> list[tuple[str, str]]:
    """Return (haplotype_id, fasta_path) from samples.tsv; requires local FASTAs present."""
    rows: list[tuple[str, str]] = []
    for line in Path(samples_tsv).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        sample, hap, _superpop, url = f[0], f[1], f[2], f[3]
        # url is a local path once fetched (prepare-haplotypes stage); accept local paths.
        fasta = url
        hid = f"{sample}#hap{hap}"
        if not Path(fasta).exists():
            raise click.ClickException(
                f"haplotype FASTA not found for {hid}: {fasta}\n"
                "Fetch the subset first (see config/samples.tsv provenance columns)."
            )
        rows.append((hid, fasta))
    return rows


def _extract_template(chm13_fasta: str, chrom: str, start: int, end: int) -> str:
    import pysam

    fa = pysam.FastaFile(chm13_fasta)
    seq = fa.fetch(chrom, start, end)
    fa.close()
    return seq


@cli.command()
@click.option("--target", required=True,
              help="CHM13 region 'chr:start-end' or a FASTA path of the locus")
@click.option("--chm13", "chm13_fasta", required=True, help="CHM13 v2.0 FASTA (indexed)")
@click.option("--samples", "samples_tsv", required=True, help="haplotype subset TSV")
@click.option("--outdir", default="results", help="output directory")
@click.option("--mode", type=click.Choice(["thermo", "rule"]), default=None,
              help="dropout model (default: from config)")
@click.option("--config", "config_path", default=None, help="params YAML (default: config/defaults.yaml)")
def run(target, chm13_fasta, samples_tsv, outdir, mode, config_path) -> None:
    """Design + evaluate + rank primers for a target across the haplotype subset."""
    import pysam

    from .bwa_backend import ensure_index, find_binding_sites_bwa
    from .design import design_candidates
    from .mask import build_excluded_regions
    from .project import AmbiguousAnchor, anchor_sequence, project_locus

    raw = cfgmod.load_raw(config_path)
    mode = mode or raw["dropout"]["mode"]
    dcfg = cfgmod.design_config(raw)
    flank = raw["design"]["flank"]
    max_mm = raw["search"]["max_mismatches"]

    # 1) anchor target on CHM13 -> (chrom, start, end) + template sequence with flank
    if ":" in target and Path(target).suffix == "":
        chrom, span = target.split(":")
        start, end = (int(x.replace(",", "")) for x in span.split("-"))
    else:
        try:
            anchor = anchor_sequence(Path(target).read_text().split("\n", 1)[1].replace("\n", ""), chm13_fasta)
        except AmbiguousAnchor as e:
            raise click.ClickException(str(e))
        chrom, start, end = anchor.chrom, anchor.start, anchor.end
    fa = pysam.FastaFile(chm13_fasta)
    clen = fa.get_reference_length(chrom)
    fa.close()
    tstart, tend = max(0, start - flank), min(clen, end + flank)
    template = _extract_template(chm13_fasta, chrom, tstart, tend)
    click.echo(f"target anchored on CHM13 {chrom}:{start}-{end} (template {tstart}-{tend})")

    # 2) project onto each haplotype
    haplos = _load_haplotypes(samples_tsv)
    projections = {}
    hap_seqs = []
    for hid, fasta in haplos:
        proj = project_locus(template, fasta)
        projections[hid] = (fasta, proj)
        if proj.haplotype_seq:
            hap_seqs.append(proj.haplotype_seq)
        click.echo(f"  {hid}: {'projected' if proj.locus else 'UNCERTAIN — ' + proj.reason}")

    # 3) variability mask -> excluded regions (coords relative to the template)
    excluded = build_excluded_regions(
        template, hap_seqs,
        min_allele_freq=raw["mask"]["min_allele_freq"],
    )
    click.echo(f"masked {len(excluded)} variable region(s)")

    # 4) design candidates
    candidates = design_candidates(template, excluded, dcfg, seq_id=f"{chrom}_{start}")
    click.echo(f"designed {len(candidates)} candidate pair(s)")

    # 5) evaluate each candidate across haplotypes (genome-wide bwa search)
    thermo_default_tm = dcfg.tm_opt
    results: list[PairResult] = []
    hap_fastas = {hid: fa for hid, fa in haplos}
    for cand in candidates:
        design_tm = _pair_tm(cand.pair.forward.sequence, cand.pair.reverse.sequence, thermo_default_tm)
        ecfg = EvalConfig(
            mode=mode, max_mismatches=max_mm,
            min_product=dcfg.product_size_min, max_product=dcfg.product_size_max,
            rule_cfg=cfgmod.rule_config(raw),
            thermo_cfg=cfgmod.thermo_config(raw, design_tm),
        )
        per_hap = []
        for hid, (fasta, proj) in projections.items():
            if proj.locus is None:
                per_hap.append(evaluate_with_sites(cand.pair, hid, [], [], None, lambda *a: "", ecfg))
                continue
            ensure_index(fasta)
            f_sites = find_binding_sites_bwa(cand.pair.forward.name, cand.pair.forward.sequence, fasta, hid, max_mm)
            r_sites = find_binding_sites_bwa(cand.pair.reverse.name, cand.pair.reverse.sequence, fasta, hid, max_mm)
            fa = pysam.FastaFile(fasta)
            res = evaluate_with_sites(
                cand.pair, hid, f_sites, r_sites, proj.locus,
                lambda c, s, e, _fa=fa: _fa.fetch(c, s, e), ecfg,
            )
            fa.close()
            per_hap.append(res)
        results.append(PairResult(cand.pair, per_hap, primer3_penalty=cand.penalty))

    # 6) rank + report
    ranked = rank_pairs(results, cfgmod.rank_config(raw))
    provenance = {
        "reference_build": "CHM13v2.0",
        "target": f"{chrom}:{start}-{end}",
        "dropout_mode": mode,
        "n_haplotypes": len(haplos),
    }
    paths = report.write_all(ranked, outdir, provenance)
    n_pass = sum(1 for rp in ranked if rp.passed)
    click.echo(f"\n{n_pass}/{len(ranked)} pairs pass filters. outputs:")
    for k, v in paths.items():
        click.echo(f"  {k}: {v}")


def _pair_tm(fwd: str, rev: str, fallback: float) -> float:
    try:
        import primer3

        return min(primer3.calc_tm(fwd), primer3.calc_tm(rev))
    except Exception:
        return fallback


# --- granular stage subcommands (one per Nextflow process) -------------------
import json  # noqa: E402


def _read_target(target: str, chm13_fasta: str):
    from .project import AmbiguousAnchor, anchor_sequence

    if ":" in target and Path(target).suffix == "":
        chrom, span = target.split(":")
        start, end = (int(x.replace(",", "")) for x in span.split("-"))
        return chrom, start, end
    seq = "".join(l for l in Path(target).read_text().splitlines() if not l.startswith(">"))
    anchor = anchor_sequence(seq, chm13_fasta)  # may raise AmbiguousAnchor
    return anchor.chrom, anchor.start, anchor.end


@cli.command(name="anchor")
@click.option("--target", required=True)
@click.option("--chm13", "chm13_fasta", required=True)
@click.option("--config", "config_path", default=None)
@click.option("--out", required=True)
def anchor_cmd(target, chm13_fasta, config_path, out) -> None:
    """Stage 1: normalize the target to a CHM13 anchor + template sequence."""
    import pysam

    from .project import AmbiguousAnchor

    raw = cfgmod.load_raw(config_path)
    flank = raw["design"]["flank"]
    try:
        chrom, start, end = _read_target(target, chm13_fasta)
    except AmbiguousAnchor as e:
        raise click.ClickException(str(e))
    fa = pysam.FastaFile(chm13_fasta)
    clen = fa.get_reference_length(chrom)
    tstart, tend = max(0, start - flank), min(clen, end + flank)
    template = fa.fetch(chrom, tstart, tend)
    fa.close()
    Path(out).write_text(json.dumps({
        "chrom": chrom, "start": start, "end": end,
        "tstart": tstart, "tend": tend, "template": template,
    }))
    click.echo(f"anchored {chrom}:{start}-{end}")


@cli.command(name="project")
@click.option("--anchor", "anchor_json", required=True)
@click.option("--hap-fasta", required=True)
@click.option("--hap-id", required=True)
@click.option("--out", required=True)
def project_cmd(anchor_json, hap_fasta, hap_id, out) -> None:
    """Stage 2: project the target onto one haplotype."""
    from .project import project_locus
    from .serialize import locus_to_dict

    anchor = json.loads(Path(anchor_json).read_text())
    proj = project_locus(anchor["template"], hap_fasta)
    Path(out).write_text(json.dumps({
        "hap_id": hap_id, "fasta": hap_fasta,
        "locus": locus_to_dict(proj.locus),
        "hap_seq": proj.haplotype_seq, "reason": proj.reason,
    }))
    click.echo(f"{hap_id}: {'ok' if proj.locus else 'uncertain'}")


@cli.command(name="design")
@click.option("--anchor", "anchor_json", required=True)
@click.option("--projection", "proj_files", multiple=True, required=True)
@click.option("--config", "config_path", default=None)
@click.option("--out", required=True)
def design_cmd(anchor_json, proj_files, config_path, out) -> None:
    """Stage 3+4: variability mask + Primer3 candidate generation."""
    from .design import design_candidates
    from .mask import build_excluded_regions
    from .serialize import candidate_to_dict

    raw = cfgmod.load_raw(config_path)
    anchor = json.loads(Path(anchor_json).read_text())
    template = anchor["template"]
    hap_seqs = [
        p["hap_seq"] for p in (json.loads(Path(f).read_text()) for f in proj_files)
        if p.get("hap_seq")
    ]
    excluded = build_excluded_regions(template, hap_seqs,
                                      min_allele_freq=raw["mask"]["min_allele_freq"])
    cands = design_candidates(template, excluded, cfgmod.design_config(raw),
                              seq_id=f"{anchor['chrom']}_{anchor['start']}")
    Path(out).write_text(json.dumps({
        "excluded": excluded, "candidates": [candidate_to_dict(c) for c in cands],
    }))
    click.echo(f"masked {len(excluded)} region(s), designed {len(cands)} candidate(s)")


@cli.command(name="evaluate")
@click.option("--candidates", "candidates_json", required=True)
@click.option("--projection", "proj_json", required=True)
@click.option("--hap-fasta", default=None, help="override the FASTA path in the projection")
@click.option("--config", "config_path", default=None)
@click.option("--mode", type=click.Choice(["thermo", "rule"]), default=None)
@click.option("--out", required=True)
def evaluate_cmd(candidates_json, proj_json, hap_fasta, config_path, mode, out) -> None:
    """Stage 5: evaluate all candidate pairs against one haplotype (genome-wide bwa)."""
    import pysam

    from .bwa_backend import ensure_index, find_binding_sites_bwa
    from .serialize import candidate_from_dict, locus_from_dict, result_to_dict

    raw = cfgmod.load_raw(config_path)
    mode = mode or raw["dropout"]["mode"]
    dcfg = cfgmod.design_config(raw)
    max_mm = raw["search"]["max_mismatches"]
    cands = [candidate_from_dict(d) for d in json.loads(Path(candidates_json).read_text())["candidates"]]
    proj = json.loads(Path(proj_json).read_text())
    hid = proj["hap_id"]
    locus = locus_from_dict(proj["locus"])
    fasta = hap_fasta or proj["fasta"]

    per_pair = {}
    if locus is None:
        for c in cands:
            r = evaluate_with_sites(c.pair, hid, [], [], None, lambda *a: "",
                                    EvalConfig(mode=mode))
            per_pair[c.pair.name] = result_to_dict(r)
    else:
        ensure_index(fasta)
        fa = pysam.FastaFile(fasta)
        for c in cands:
            design_tm = _pair_tm(c.pair.forward.sequence, c.pair.reverse.sequence, dcfg.tm_opt)
            ecfg = EvalConfig(
                mode=mode, max_mismatches=max_mm,
                min_product=dcfg.product_size_min, max_product=dcfg.product_size_max,
                rule_cfg=cfgmod.rule_config(raw),
                thermo_cfg=cfgmod.thermo_config(raw, design_tm),
            )
            f_sites = find_binding_sites_bwa(c.pair.forward.name, c.pair.forward.sequence, fasta, hid, max_mm)
            r_sites = find_binding_sites_bwa(c.pair.reverse.name, c.pair.reverse.sequence, fasta, hid, max_mm)
            r = evaluate_with_sites(c.pair, hid, f_sites, r_sites, locus,
                                    lambda c2, s, e, _fa=fa: _fa.fetch(c2, s, e), ecfg)
            per_pair[c.pair.name] = result_to_dict(r)
        fa.close()
    Path(out).write_text(json.dumps({"hap_id": hid, "per_pair": per_pair}))
    click.echo(f"evaluated {len(cands)} pair(s) on {hid}")


@cli.command(name="aggregate")
@click.option("--anchor", "anchor_json", required=True)
@click.option("--candidates", "candidates_json", required=True)
@click.option("--result", "result_files", multiple=True, required=True)
@click.option("--config", "config_path", default=None)
@click.option("--outdir", required=True)
def aggregate_cmd(anchor_json, candidates_json, result_files, config_path, outdir) -> None:
    """Stage 6: gather per-haplotype results, rank, and report."""
    from .serialize import candidate_from_dict, result_from_dict

    raw = cfgmod.load_raw(config_path)
    anchor = json.loads(Path(anchor_json).read_text())
    cands = [candidate_from_dict(d) for d in json.loads(Path(candidates_json).read_text())["candidates"]]
    per_hap_files = [json.loads(Path(f).read_text()) for f in result_files]

    results: list[PairResult] = []
    for c in cands:
        hap_results = [
            result_from_dict(h["per_pair"][c.pair.name])
            for h in per_hap_files if c.pair.name in h["per_pair"]
        ]
        results.append(PairResult(c.pair, hap_results, primer3_penalty=c.penalty))
    ranked = rank_pairs(results, cfgmod.rank_config(raw))
    provenance = {
        "reference_build": "CHM13v2.0",
        "target": f"{anchor['chrom']}:{anchor['start']}-{anchor['end']}",
        "dropout_mode": raw["dropout"]["mode"],
        "n_haplotypes": len(per_hap_files),
    }
    paths = report.write_all(ranked, outdir, provenance)
    n_pass = sum(1 for rp in ranked if rp.passed)
    click.echo(f"{n_pass}/{len(ranked)} pairs pass. report: {paths['html']}")


if __name__ == "__main__":
    cli()
