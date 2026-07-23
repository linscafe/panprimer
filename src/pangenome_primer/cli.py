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
        if not line or line.startswith("#") or line.startswith("sample\t"):
            continue  # skip comments and the header row
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
@click.option("--top-k", default=5, show_default=True,
              help="genome-wide off-target search (stage B) runs only on the top-K by "
                   "on-target coverage; 0 = all candidates")
@click.option("--target-assembly", type=click.Choice(["chm13", "grch38"]), default="chm13",
              show_default=True, help="coordinate system of a chr:start-end --target")
@click.option("--grch38", "grch38_fasta", default=None,
              help="GRCh38 FASTA (required when --target-assembly grch38)")
def run(target, chm13_fasta, samples_tsv, outdir, mode, config_path, top_k,
        target_assembly, grch38_fasta) -> None:
    """Design + evaluate + rank primers for a target across the haplotype subset.

    --target accepts a CHM13 chr:start-end (default), a GRCh38 chr:start-end
    (--target-assembly grch38 --grch38 <fa>), or a FASTA of the locus.

    Two-stage: stage A evaluates on-target coverage/dropout cheaply against each projected
    homologous window (no genome-wide index); stage B runs the expensive genome-wide bwa
    off-target search only on the top-K coverage candidates."""
    import pysam

    from .bwa_backend import find_binding_sites_batch
    from .design import design_candidates
    from .engine import HaplotypeContext, evaluate_pair
    from .mask import build_excluded_regions
    from .model import Locus
    from .project import AmbiguousAnchor, project_target, resolve_target

    raw = cfgmod.load_raw(config_path)
    mode = mode or raw["dropout"]["mode"]
    dcfg = cfgmod.design_config(raw)
    flank = raw["design"]["flank"]
    max_mm = raw["search"]["max_mismatches"]

    # 1) resolve target to a CHM13 (chrom, start, end) + template sequence with flank
    try:
        chrom, start, end = resolve_target(
            target, chm13_fasta, assembly=target_assembly, grch38_fasta=grch38_fasta)
    except AmbiguousAnchor as e:
        raise click.ClickException(str(e))
    except (ValueError, FileNotFoundError) as e:
        raise click.ClickException(str(e))
    if target_assembly == "grch38":
        click.echo(f"GRCh38 {target} -> CHM13 {chrom}:{start}-{end}")
    fa = pysam.FastaFile(chm13_fasta)
    clen = fa.get_reference_length(chrom)
    fa.close()
    tstart, tend = max(0, start - flank), min(clen, end + flank)
    template = _extract_template(chm13_fasta, chrom, tstart, tend)
    click.echo(f"target anchored on CHM13 {chrom}:{start}-{end} (template {tstart}-{tend})")

    # 2) project onto each haplotype (PAF cache lift, or on-the-fly fallback)
    haplos = _load_haplotypes(samples_tsv)
    projections = {}
    hap_seqs = []
    for hid, fasta in haplos:
        proj = project_target(chrom, tstart, tend, template, fasta)
        projections[hid] = (fasta, proj)
        if proj.haplotype_seq:
            hap_seqs.append(proj.haplotype_seq)
        click.echo(f"  {hid}: {'projected (' + proj.reason + ')' if proj.locus else 'UNCERTAIN — ' + proj.reason}")

    # 3) variability mask -> excluded regions (coords relative to the template)
    excluded = build_excluded_regions(
        template, hap_seqs,
        min_allele_freq=raw["mask"]["min_allele_freq"],
        merge_gap=raw["mask"].get("merge_gap", 25),
        max_regions=raw["mask"].get("max_regions", 200),
    )
    click.echo(f"masked {len(excluded)} variable region(s)")

    # 4) design candidates
    candidates = design_candidates(template, excluded, dcfg, seq_id=f"{chrom}_{start}")
    click.echo(f"designed {len(candidates)} candidate pair(s)")

    def _ecfg(cand):
        design_tm = _pair_tm(cand.pair.forward.sequence, cand.pair.reverse.sequence, dcfg.tm_opt)
        return EvalConfig(
            mode=mode, max_mismatches=max_mm,
            min_product=dcfg.product_size_min, max_product=dcfg.product_size_max,
            rule_cfg=cfgmod.rule_config(raw), thermo_cfg=cfgmod.thermo_config(raw, design_tm),
        )

    # 5a) STAGE A — cheap on-target coverage against each projected window (no genome-wide
    # index). Distinguishes pass/dropout; off-target elsewhere is a stage-B concern.
    windows = {
        hid: HaplotypeContext(
            hid, (proj.locus.chrom if proj.locus else "?"),
            proj.haplotype_seq,
            Locus(hid, proj.locus.chrom, 0, len(proj.haplotype_seq)) if proj.locus else None,
            projection_ok=proj.locus is not None,
        )
        for hid, (_f, proj) in projections.items()
    }
    stageA = [
        evaluate_pair(c.pair, list(windows.values()), _ecfg(c), primer3_penalty=c.penalty)
        for c in candidates
    ]
    order = sorted(range(len(candidates)),
                   key=lambda i: (-stageA[i].on_target_coverage, stageA[i].primer3_penalty))
    k = len(candidates) if not top_k else min(top_k, len(candidates))
    shortlist = order[:k]
    click.echo(f"stage A done; genome-wide specificity (stage B) on top {k} by coverage")

    # 5b) STAGE B — genome-wide bwa off-target search for the shortlist only. One batched
    # search per haplotype covers all shortlisted primers (index loaded once).
    short_cands = [candidates[i] for i in shortlist]
    short_seqs = ([c.pair.forward.sequence for c in short_cands]
                  + [c.pair.reverse.sequence for c in short_cands])
    site_cache: dict[str, dict] = {}

    def _sites(fasta, hid):
        if hid not in site_cache:
            click.echo(f"    stage B: genome-wide search on {hid} ...")
            site_cache[hid] = find_binding_sites_batch(short_seqs, fasta, hid, max_mm)
        return site_cache[hid]

    results: list[PairResult] = []
    for cand in short_cands:
        ecfg = _ecfg(cand)
        per_hap = []
        for hid, (fasta, proj) in projections.items():
            if proj.locus is None:
                per_hap.append(evaluate_with_sites(cand.pair, hid, [], [], None, lambda *a: "", ecfg))
                continue
            sbs = _sites(fasta, hid)
            fa = pysam.FastaFile(fasta)
            per_hap.append(evaluate_with_sites(
                cand.pair, hid,
                sbs.get(cand.pair.forward.sequence, []), sbs.get(cand.pair.reverse.sequence, []),
                proj.locus, lambda c, s, e, _fa=fa: _fa.fetch(c, s, e), ecfg,
            ))
            fa.close()
        results.append(PairResult(cand.pair, per_hap, primer3_penalty=cand.penalty))

    # 6) rank + report
    ranked = rank_pairs(results, cfgmod.rank_config(raw))
    provenance = {
        "reference_build": "CHM13v2.0", "target": f"{chrom}:{start}-{end}",
        "dropout_mode": mode, "n_haplotypes": len(haplos),
        "candidates": len(candidates), "specificity_checked": k,
    }
    paths = report.write_all(ranked, outdir, provenance)
    n_pass = sum(1 for rp in ranked if rp.passed)
    click.echo(f"\n{n_pass}/{len(ranked)} specificity-checked pairs pass filters. outputs:")
    for kk, v in paths.items():
        click.echo(f"  {kk}: {v}")


def _pair_tm(fwd: str, rev: str, fallback: float) -> float:
    try:
        import primer3

        return min(primer3.calc_tm(fwd), primer3.calc_tm(rev))
    except Exception:
        return fallback


# --- granular stage subcommands (one per Nextflow process) -------------------
import json  # noqa: E402


@cli.command(name="anchor")
@click.option("--target", required=True)
@click.option("--chm13", "chm13_fasta", required=True)
@click.option("--target-assembly", type=click.Choice(["chm13", "grch38"]), default="chm13")
@click.option("--grch38", "grch38_fasta", default=None)
@click.option("--config", "config_path", default=None)
@click.option("--out", required=True)
def anchor_cmd(target, chm13_fasta, target_assembly, grch38_fasta, config_path, out) -> None:
    """Stage 1: normalize the target (CHM13/GRCh38 coords or FASTA) to a CHM13 anchor."""
    import pysam

    from .project import AmbiguousAnchor, resolve_target

    raw = cfgmod.load_raw(config_path)
    flank = raw["design"]["flank"]
    try:
        chrom, start, end = resolve_target(
            target, chm13_fasta, assembly=target_assembly, grch38_fasta=grch38_fasta)
    except AmbiguousAnchor as e:
        raise click.ClickException(str(e))
    except (ValueError, FileNotFoundError) as e:
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
    """Stage 2: project the target onto one haplotype (cache-aware: PAF lift or .mmi)."""
    from .project import project_target
    from .serialize import locus_to_dict

    a = json.loads(Path(anchor_json).read_text())
    proj = project_target(a["chrom"], a["tstart"], a["tend"], a["template"], hap_fasta)
    Path(out).write_text(json.dumps({
        "hap_id": hap_id, "fasta": hap_fasta,
        "locus": locus_to_dict(proj.locus),
        "hap_seq": proj.haplotype_seq, "reason": proj.reason,
    }))
    click.echo(f"{hap_id}: {proj.reason if proj.locus else 'uncertain'}")


@cli.command(name="stage-a")
@click.option("--candidates", "candidates_json", required=True)
@click.option("--projection", "proj_files", multiple=True, required=True)
@click.option("--config", "config_path", default=None)
@click.option("--mode", type=click.Choice(["thermo", "rule"]), default=None)
@click.option("--top-k", default=5, help="keep this many candidates by on-target coverage; 0 = all")
@click.option("--out", required=True)
def stage_a_cmd(candidates_json, proj_files, config_path, mode, top_k, out) -> None:
    """Stage A: rank candidates by on-target coverage against each projected window (no
    genome-wide index), and write a shortlist for the expensive stage-B specificity search."""
    from .engine import HaplotypeContext, evaluate_pair
    from .model import Locus
    from .serialize import candidate_from_dict, candidate_to_dict, locus_from_dict

    raw = cfgmod.load_raw(config_path)
    mode = mode or raw["dropout"]["mode"]
    dcfg = cfgmod.design_config(raw)
    cands = [candidate_from_dict(d) for d in json.loads(Path(candidates_json).read_text())["candidates"]]
    projs = [json.loads(Path(f).read_text()) for f in proj_files]
    windows = []
    for p in projs:
        loc = locus_from_dict(p["locus"])
        windows.append(HaplotypeContext(
            p["hap_id"], loc.chrom if loc else "?", p.get("hap_seq", ""),
            Locus(p["hap_id"], loc.chrom, 0, len(p["hap_seq"])) if loc else None,
            projection_ok=loc is not None))

    def ecfg(c):
        design_tm = _pair_tm(c.pair.forward.sequence, c.pair.reverse.sequence, dcfg.tm_opt)
        return EvalConfig(mode=mode, min_product=dcfg.product_size_min,
                          max_product=dcfg.product_size_max, rule_cfg=cfgmod.rule_config(raw),
                          thermo_cfg=cfgmod.thermo_config(raw, design_tm))

    cov = {c.pair.name: evaluate_pair(c.pair, windows, ecfg(c)).on_target_coverage for c in cands}
    order = sorted(cands, key=lambda c: (-cov[c.pair.name], c.penalty))
    k = len(cands) if not top_k else min(top_k, len(cands))
    shortlist = order[:k]
    Path(out).write_text(json.dumps({
        "candidates": [candidate_to_dict(c) for c in shortlist],
        "stage_a_coverage": {c.pair.name: round(cov[c.pair.name], 4) for c in cands},
    }))
    click.echo(f"stage A: {len(cands)} candidates -> shortlist {k} by coverage")


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

    from .bwa_backend import find_binding_sites_batch
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
        # one genome-wide bwa search for ALL primers on this haplotype (index loaded once)
        all_seqs = [c.pair.forward.sequence for c in cands] + [c.pair.reverse.sequence for c in cands]
        sites_by_seq = find_binding_sites_batch(all_seqs, fasta, hid, max_mm)
        fa = pysam.FastaFile(fasta)
        for c in cands:
            design_tm = _pair_tm(c.pair.forward.sequence, c.pair.reverse.sequence, dcfg.tm_opt)
            ecfg = EvalConfig(
                mode=mode, max_mismatches=max_mm,
                min_product=dcfg.product_size_min, max_product=dcfg.product_size_max,
                rule_cfg=cfgmod.rule_config(raw),
                thermo_cfg=cfgmod.thermo_config(raw, design_tm),
            )
            f_sites = sites_by_seq.get(c.pair.forward.sequence, [])
            r_sites = sites_by_seq.get(c.pair.reverse.sequence, [])
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
        "specificity_checked": len(cands),
    }
    paths = report.write_all(ranked, outdir, provenance)
    n_pass = sum(1 for rp in ranked if rp.passed)
    click.echo(f"{n_pass}/{len(ranked)} pairs pass. report: {paths['html']}")


@cli.command(name="rerender")
@click.option("--results-json", required=True, help="an existing results.json")
@click.option("--outdir", required=True)
@click.option("--config", "config_path", default=None)
def rerender_cmd(results_json, outdir, config_path) -> None:
    """Regenerate TSV/JSON/HTML from a saved results.json with the current code — applies
    metric/template fixes without re-running the pipeline."""
    raw = cfgmod.load_raw(config_path)
    paths = report.rerender_from_json(results_json, outdir, cfgmod.rank_config(raw))
    click.echo("re-rendered:")
    for k, v in paths.items():
        click.echo(f"  {k}: {v}")


@cli.command(name="fetch-subset")
@click.option("--sample", "sample_ids", multiple=True,
              help="explicit HPRC R2 sample id(s); repeatable. Overrides --per-superpop.")
@click.option("--per-superpop", default=1, show_default=True,
              help="samples to pick from each of AFR/EUR/EAS/SAS/AMR (diverse default)")
@click.option("--data-dir", default="hprc-r2/assemblies", show_default=True,
              help="where assembly FASTAs are/should be stored")
@click.option("--out", "out_path", default="config/samples.tsv", show_default=True)
@click.option("--download/--manifest-only", default=False,
              help="actually download (multi-GB each) + verify md5, or just write the manifest")
@click.option("--limit", default=0, help="cap number of haplotypes (0 = no cap)")
def fetch_subset_cmd(sample_ids, per_superpop, data_dir, out_path, download, limit) -> None:
    """Select a diverse HPRC R2 haplotype subset and write config/samples.tsv.

    Default is manifest-only: fetch the official index + sample metadata, pick the subset,
    pin provenance (md5 + S3 URL + accession), and write the manifest. Pass --download to
    pull the assemblies (≈1 GB each) into --data-dir and verify their md5s."""
    from . import hprc

    click.echo("fetching HPRC R2 assembly index + sample metadata ...")
    records = hprc.load_records()
    if sample_ids:
        chosen = hprc.select_samples(records, list(sample_ids))
        if not chosen:
            raise click.ClickException(f"no R2 haplotypes for samples: {list(sample_ids)}")
    else:
        chosen = hprc.select_diverse(records, per_superpop=per_superpop)
    if limit:
        chosen = chosen[:limit]

    click.echo(f"selected {len(chosen)} haplotype(s):")
    for r in chosen:
        click.echo(f"  {r.hap_id:18s} {r.superpop:4s} {r.population:6s} {r.assembly_url}")

    md5s: dict[str, str] = {}
    if download:
        for r in chosen:
            click.echo(f"downloading {r.hap_id} ...")
            hprc.download_haplotype(r, data_dir)
            md5s[r.hap_id] = hprc.fetch_md5(r)
            click.echo(f"  ok -> {r.local_path(data_dir)}")
    else:
        click.echo("pinning md5s (manifest-only; no assembly download) ...")
        for r in chosen:
            try:
                md5s[r.hap_id] = hprc.fetch_md5(r)
            except Exception:
                md5s[r.hap_id] = "PENDING"

    hprc.write_samples_tsv(chosen, out_path, data_dir, md5s)
    click.echo(f"\nwrote {out_path}")
    if not download:
        click.echo("next: `pangenome-primer fetch-subset --download` (or re-run with --download) "
                   "to pull the assemblies before `run`/Nextflow.")


if __name__ == "__main__":
    cli()
