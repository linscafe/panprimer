# Pangenome PCR Primer Design

Design **universal PCR primers against the HPRC R2 human pangenome** instead of a single
reference. Primers are evaluated across many phased human haplotypes, so allele **dropout**
(a binding-site SNP/indel under a primer 3′ end) and **off-target** amplification (paralogs
a single reference collapses) are caught before you order oligos.

**New here?** Start with the plain-language introductions:
[design pipeline](docs/intro_design_pipeline.md) · [verify pipeline](docs/intro_verify_pipeline.md)
(they explain what a pangenome is and why the pipelines work the way they do).

See also [`CONTEXT.md`](CONTEXT.md) for the ubiquitous language and the design docs in
[`docs/`](docs/) — [execution plan (design)](docs/execution_plan_design_pipeline.md), [execution plan (verify)](docs/execution_plan_verify_pipeline.md),
[config_reference](docs/config_reference.md).

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
extracted from the GRCh38 reference and aligned to CHM13), or a **FASTA** of the locus:

```bash
# GRCh38 coordinates: --grch38 is the GRCh38 *reference genome* (indexed .fai) — coordinates
# are absolute, so the slice is fetched from it (random-access; the whole genome is not read).
pangenome-primer run --target chr22:42126499-42130810 --target-assembly grch38 \
    --grch38 GRCh38.fa --chm13 CHM13v2.0.fa --samples config/samples.tsv --outdir results

# Just a locus sequence, no coordinates? Use FASTA input (no reference genome needed):
pangenome-primer run --target my_locus.fa --chm13 CHM13v2.0.fa --samples config/samples.tsv
```

Nextflow (parallel per-haplotype; the cloud-scale path):

```bash
nextflow run main.nf -profile local \
    --target chr1:1200-1400 --chm13 CHM13v2.0.fa --samples config/samples.tsv --outdir results
```

Fill `config/samples.tsv` with the HPRC R2 subset (diverse superpopulations) and pin
provenance (URL + sha256 + release) before a real run. Tune thresholds in
`config/defaults.yaml`.

## Verify existing primers

Screen primer pairs you already have (no design). Input CSV — `primer_id,target,forward,reverse`
where `target` is the **GRCh38** intended-amplicon region:

```csv
primer_id,target,forward,reverse
GAPDH_ex,chr12:6534000-6534240,CGCTTCATGCTGCACATCTC,TTCAGTAATGGCTGCCTGGG
```

```bash
pangenome-primer verify --primers primers.csv --chm13 CHM13v2.0.fa --grch38 GRCh38.fa \
    --samples config/samples.tsv --outdir verify_out       # --target-assembly chm13 to use CHM13 coords
```

Produces `verify_matrix.html` (+ `verify.json`/`.tsv`): a matrix, one row per pair, one column
per haplotype. Cells show predicted product sizes — **green** correct (on-target), **red**
off-target, grey **dropout** (a binding-site variant kills the product), `?` not projectable.
See `docs/execution_plan_verify_pipeline.md`.

## Status

**Done**
- **Engine + reporting**: implemented and tested (27 unit/integration tests; `selftest`).
- **Input**: CHM13 region, **GRCh38 region** (`--target-assembly grch38 --grch38 <fa>`;
  extracted + aligned to CHM13), or FASTA of the locus. Ambiguous mappings fail loud.
- **Real HPRC R2 data**: `fetch-subset` downloader (official index, md5-pinned); 10-haplotype
  subset downloaded + `bwa`-indexed; CHM13 v2.0 reference in place.
- **Real-locus validation**: **GAPDH** (easy) → clean universal primer at 100% coverage;
  **CYP2D6** (hard, full-specificity `--top-k 20`) → 14/20 clean, and the tool *catches both*
  failure modes: CYP2D7 pseudogene co-amplification (multi-product) and reference-bias
  dropout (CHM13-designed primers that fail on real haplotypes). Reports in `demo-design-pipeline/`.
- **Runtime**: two-stage search (`--top-k`: cheap coverage → genome-wide specificity on the
  shortlist) + `.mmi`-cached projection. 3-haplotype demo ≈ 9 min, ~11.8 GB peak
  (`demo/run_demo.sh`). See `docs/execution_plan_design_pipeline.md`.
- **Nextflow pipeline**: **two-stage + cache-aware projection** (STAGE_A coverage shortlist →
  genome-wide EVALUATE on the top-K), with prebuilt per-haplotype caches staged so no index is
  rebuilt in a work dir. Validated end-to-end on a synthetic mini-pangenome (`--top_k`).
- **PAF projection at scale**: opt-in `BUILD_PAF=1` in `scripts/prepare_haplotypes.sh`
  (threads via `MM_THREADS`) builds the whole-genome projection cache; the Nextflow/`run`
  projection uses it automatically when present. Too heavy for 15 GB locally — a cloud lever.
- **Verify mode** (`pangenome-primer verify`): screen user-supplied primer pairs from a CSV;
  amplicon-size matrix (green on-target / red off-target / grey dropout). See `docs/execution_plan_verify_pipeline.md`.

**Not yet done**
- Full ~460-haplotype / cloud run; the graph annotation layer; automated primer rescue.
  See the plan's "out of scope (v2+)".
