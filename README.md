# Pangenome PCR Primer Design

Design **universal PCR primers against the HPRC R2 human pangenome** instead of a single
reference. Primers are evaluated across many phased human haplotypes, so allele **dropout**
(a binding-site SNP/indel under a primer 3′ end) and **off-target** amplification (paralogs
a single reference collapses) are caught before you order oligos.

See [`CONTEXT.md`](CONTEXT.md) for the ubiquitous language and
[the plan](../.claude/plans) for the architecture and the decisions behind it.

## How it works

Target (CHM13 region, GRCh38 region, or FASTA) → anchor on **CHM13 v2.0** → project onto each haplotype
(the *expected homologous locus*) → mask variable sites → **Primer3** candidates → genome-wide
binding search per haplotype (**bwa**) → **thermodynamic** dropout classification (3′-end
aware) → pair into amplicons → per-haplotype status (`pass/dropout/off_target/multi_product/
uncertain`) → rank by coverage + specificity → TSV/JSON + HTML report.

## Setup

```bash
mamba env create -f env/environment.yml
mamba activate pangenome-primer
pip install -e .
```

## Quick check (no data needed)

```bash
pangenome-primer selftest --outdir selftest_out   # synthetic mini-pangenome, all 5 statuses
pytest -q                                          # engine correctness anchor
```

## Get the data

Select a diverse HPRC R2 subset (one sample per superpopulation, both haplotypes) and pin
provenance from the official index — manifest first, then the multi-GB pull:

```bash
pangenome-primer fetch-subset --per-superpop 1            # writes config/samples.tsv (manifest)
pangenome-primer fetch-subset --per-superpop 1 --download # pulls ~1 GB/haplotype + verifies md5
# or explicit samples:  pangenome-primer fetch-subset --sample HG00408 --sample HG02602 --download
```

CHM13 v2.0 reference (for `--chm13`), from the same no-egress bucket:

```bash
curl -O https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gunzip chm13v2.0.fa.gz && samtools faidx chm13v2.0.fa
```

**Prepare the haplotypes** (download → gunzip → faidx → `bwa index`). Indexing a ~3 Gb
assembly is slow (tens of minutes, ~5 GB RAM each), so for the full subset run this
unattended — it is resumable and logs to `hprc-r2/prepare.log`:

```bash
bash scripts/prepare_haplotypes.sh
```

This is the local bottleneck; the containerized Nextflow `container` profile is the
scale-out path when you outgrow one machine.

## Run

Single-process (subset, local dev):

```bash
pangenome-primer run --target chr1:1200-1400 --chm13 CHM13v2.0.fa \
    --samples config/samples.tsv --outdir results
```

The `--target` may be a **CHM13 region** (default), a **GRCh38 region** (the slice is
extracted and aligned to CHM13), or a **FASTA** of the locus:

```bash
# GRCh38 coordinates (needs a GRCh38 FASTA)
pangenome-primer run --target chr22:42126499-42130810 --target-assembly grch38 \
    --grch38 GRCh38.fa --chm13 CHM13v2.0.fa --samples config/samples.tsv --outdir results
```

Nextflow (parallel per-haplotype; the cloud-scale path):

```bash
nextflow run main.nf -profile local \
    --target chr1:1200-1400 --chm13 CHM13v2.0.fa --samples config/samples.tsv --outdir results
```

Fill `config/samples.tsv` with the HPRC R2 subset (diverse superpopulations) and pin
provenance (URL + sha256 + release) before a real run. Tune thresholds in
`config/defaults.yaml`.

## Status

**Done**
- **Engine + reporting**: implemented and tested (18 unit/integration tests; `selftest`).
- **Real HPRC R2 data**: `fetch-subset` downloader (official index, md5-pinned); 10-haplotype
  subset downloaded + `bwa`-indexed; CHM13 v2.0 reference in place.
- **Real-locus validation**: **GAPDH** (easy) → clean universal primer at 100% coverage;
  **CYP2D6** (hard, full-specificity `--top-k 20`) → 14/20 clean, and the tool *catches both*
  failure modes: CYP2D7 pseudogene co-amplification (multi-product) and reference-bias
  dropout (CHM13-designed primers that fail on real haplotypes). Reports in `demo/results/`.
- **Runtime**: two-stage search (`--top-k`: cheap coverage → genome-wide specificity on the
  shortlist) + `.mmi`-cached projection. 3-haplotype demo ≈ 9 min, ~11.8 GB peak
  (`demo/run_demo.sh`). See `docs/runtime_plan.md`.
- **Nextflow pipeline**: **two-stage + cache-aware projection** (STAGE_A coverage shortlist →
  genome-wide EVALUATE on the top-K), with prebuilt per-haplotype caches staged so no index is
  rebuilt in a work dir. Validated end-to-end on a synthetic mini-pangenome (`--top_k`).
- **PAF projection at scale**: opt-in `BUILD_PAF=1` in `scripts/prepare_haplotypes.sh`
  (threads via `MM_THREADS`) builds the whole-genome projection cache; the Nextflow/`run`
  projection uses it automatically when present. Too heavy for 15 GB locally — a cloud lever.

- **Input**: CHM13 region, **GRCh38 region** (`--target-assembly grch38 --grch38 <fa>`;
  extracted + aligned to CHM13), or FASTA of the locus. Ambiguous mappings fail loud.

**Not yet done**

**Cloud run (not yet done)**
- Full ~460-haplotype / cloud run; the graph annotation layer; automated primer rescue.
  See the plan's "out of scope (v2+)".
