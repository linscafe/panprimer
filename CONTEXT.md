# CONTEXT — Ubiquitous Language

Glossary for the pangenome PCR primer design tool. Terms here are canonical: use them in code (class/function names), CLI, config, and reports. This file is a glossary only — no implementation details, no decisions (those live in the plan / ADRs).

## Core terms

- **Reference genome** — a single consensus human sequence (GRCh38, CHM13). Convenient but lossy: it hides variation and collapses repeated regions into one copy.

- **Pangenome** — many complete genomes together, replacing that single sequence. Here HPRC Release 2: near-gapless assemblies from ~232 individuals across global populations, ~460 haplotypes.

- **CHM13 v2.0** — a complete, gap-free telomere-to-telomere reference, used here **only as the coordinate backbone** (a common ruler for naming positions). It never decides whether a primer works; that verdict comes from the haplotypes.

- **In-silico PCR** — predicting products by pairing forward and reverse binding sites that point at each other within a plausible distance. A single binding site is not a product.

- **Target locus** — the genomic region we want to amplify. Anchored on **CHM13 v2.0**, the primary coordinate system. All other inputs (GRCh38 coords, raw FASTA) are normalized to a CHM13 anchor before anything else happens.

- **Haplotype** — one phased assembly (hap1 or hap2) of one HPRC R2 individual. The atomic unit of evaluation: coverage and specificity are counted over haplotypes, not individuals.

- **Expected homologous locus** — the region in a given haplotype that the target locus projects onto (by lifting the CHM13 target ± flank through a cached alignment, or aligning it on the fly). An amplicon is **on-target** only if it falls within this window; anything else is **off-target**.

- **Binding site** — a location in a haplotype where a primer anneals, with its strand, coordinates, per-position mismatches, and the mismatch offset from the primer's 3′ end.

- **Amplicon** — a predicted PCR product: a forward and reverse binding site on one haplotype in convergent orientation within the allowed product-size window.

- **Dropout** — a haplotype where the target *should* amplify but a SNP/indel under a primer (typically near the 3′ end) prevents it.

- **Off-target product** — an amplicon outside the expected homologous locus.

## Per-haplotype status

Every evaluated haplotype gets exactly one status for a given primer pair:

- **pass** — one valid on-target amplicon, no acceptable off-target product.
- **dropout** — no valid on-target amplicon (binding-site variant prevents amplification).
- **off_target** — an acceptable amplicon forms outside the expected homologous locus.
- **multi_product** — more than one acceptable amplicon (on- and/or off-target).
- **uncertain** — the locus is broken/unmappable in this assembly, or its projection is unreliable. Never counted as failure: excluded from the denominator and reported separately.

## Metrics

- **On-target coverage** — haplotypes that produce the intended target amplicon *at all* (≥1 on-target amplicon: `pass` **or** `multi_product` that includes an on-target band) ÷ **evaluable** haplotypes (evaluable = total − uncertain). The anti-dropout metric: "does the target amplify across diversity?" Always report the denominator, not just the percentage.

- **Unique product rate** — haplotypes with exactly one target-specific product and no acceptable off-target (i.e. `pass`) ÷ evaluable haplotypes. The specificity metric: "is the reaction clean?" Coverage can be high while this is low (target amplifies everywhere but always with extra bands).
