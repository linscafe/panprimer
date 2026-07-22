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
