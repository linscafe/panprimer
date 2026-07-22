# Pangenome PCR Primer Design

Design **universal PCR primers against the HPRC R2 human pangenome** instead of a single
reference. Primers are evaluated across many phased human haplotypes, so allele **dropout**
(a binding-site SNP/indel under a primer 3′ end) and **off-target** amplification (paralogs
a single reference collapses) are caught before you order oligos.

See [`CONTEXT.md`](CONTEXT.md) for the ubiquitous language and
[the plan](../.claude/plans) for the architecture and the decisions behind it.

## How it works

Target (CHM13 region or FASTA) → anchor on **CHM13 v2.0** → project onto each haplotype
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

Nextflow (parallel per-haplotype; the cloud-scale path):

```bash
nextflow run main.nf -profile local \
    --target chr1:1200-1400 --chm13 CHM13v2.0.fa --samples config/samples.tsv --outdir results
```

Fill `config/samples.tsv` with the HPRC R2 subset (diverse superpopulations) and pin
provenance (URL + sha256 + release) before a real run. Tune thresholds in
`config/defaults.yaml`.

## Status

- **Engine + reporting**: implemented and tested (8 unit/integration tests; `selftest`).
- **Nextflow pipeline**: validated end-to-end on a synthetic mini-pangenome.
- **Not yet done**: fetching the real HPRC R2 subset; validation on real loci
  (housekeeping gene, then a hard case like HLA/CYP2D6); GRCh38-coord input; the graph
  annotation layer; automated primer rescue. See the plan's "out of scope (v2+)".
