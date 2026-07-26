# An introduction to the primer **verification** pipeline

> New to pangenomes? [`intro_design_pipeline.md`](intro_design_pipeline.md) §2 has the fuller explanation. Both pipelines share an engine but answer opposite questions.

## 1. The problem

Sometimes you don't need to *design* primers — you already **have** a pair, from a paper, an inherited protocol, or a kit. The question is:

> *Will these primers behave across the diversity of real human genomes — or drop out in some people, or throw extra bands from a look-alike locus?*

You hand over a list of pairs; the pipeline predicts, for each pair in each of many real genomes, what you would see on a gel: the size of the correct product and of any off-target products. **It designs nothing.**

Two things a single-reference check cannot tell you:

- **Does it drop out?** A common variant under a primer's 3′ end means no band in carriers.
- **Does it amplify a look-alike?** Paralogs and pseudogenes (CYP2D6 vs CYP2D7) give extra bands, invisible on a reference that collapses them.

## 2. Pangenome in one minute

The usual reference is a *single* consensus sequence — an average matching no real person, merging repeats into one copy. A **pangenome** is a collection of many complete real genomes: here **HPRC Release 2**, ~232 people, each contributing two phased **haplotypes**, so up to **~460 genomes** to test against.

**CHM13 v2.0** serves only as a **coordinate ruler**; verdicts come from the haplotypes.

## 3. How it works

**Input:** a CSV, one row per pair.

| column | meaning |
|---|---|
| `primer_id` | a label you choose |
| `target` | intended amplicon region, **GRCh38** coordinates (how published primers are usually annotated) |
| `forward` / `reverse` | primer sequences, 5′→3′ |

```
primers.csv
    │  (1) locate each target on CHM13 ── from GRCh38 coordinates
    ▼
    │  (2) project onto every haplotype ── the expected "correct" window per person
    ▼
    │  (3) search all primers genome-wide, per haplotype ── one pass covers every pair
    ▼
    │  (4) keep sites that would truly amplify (thermodynamic, 3′-aware)
    ▼
    │  (5) pair into products; record size + on/off-target
    ▼
  matrix report
```

**(1) Locate.** The GRCh38 slice is fetched and aligned to CHM13, recording where the intended amplicon lives on the common ruler. CHM13 coordinates work directly too.

**(2) Project.** For every genome, find the homologous window. Products **there** are on-target; anywhere else is off-target. If the locus can't be confidently found, the cell is **uncertain** (`?`).

**(3) Search genome-wide.** The scanner reads each haplotype's compressed sequence straight through and finds everywhere each primer could bind. The search is **exhaustive** — every position is checked, so no site is missed to a heuristic. One pass covers all your pairs at once, so screening a panel costs little more than screening one pair.

**(4) Keep what amplifies.** A thermodynamic, 3′-aware model filters the many harmless near-matches.

**(5) Pair into products.** Surviving forward/reverse sites pointing at each other within range (default ≤ 2 kb) become predicted amplicons, each with a size and an on/off-target label.

## 4. Reading the output

The matrix has **one column per primer pair and one row per haplotype**. Haplotypes go on the vertical axis because that is the count that grows — 3, 30, or all ~460 — while pairs stay few. The first rows carry each pair's expected size and a **summary**: on-target coverage plus a count of each failure mode.

| cell | meaning |
|---|---|
| **green** number | the correct on-target amplicon. A few bases' difference between people is real biology (a small indel inside the amplicon); large deviations are underlined |
| **red** number | an off-target product |
| grey **dropout** | no correct product — a binding-site variant killed it |
| **`?`** | locus not projectable in that assembly: unknown, not a failure |
| **`>100 binding sites`** | a repeat-derived primer (e.g. on an *Alu*). Sizes are not reported, because a primer annealing in ~100,000 places amplifies indiscriminately and listing bands would imply precision that isn't there |

Read a column top to bottom and you see the gel you'd get across many people. A machine-readable **TSV/JSON** is written alongside.

> [!CAUTION]
> **A cell is an individual genotype, not a population frequency.** With few haplotypes a dropout will inevitably land on some superpopulation label, and it is tempting to call the variant specific to it. Usually it is not — it is a common variant your sample was too small to place. The bundled 3-haplotype demo shows the trap exactly: its CYP2D6 dropout falls on the European haplotype, yet at 30 haplotypes the same pair fails in 11 of 24 evaluable, in every superpopulation. For a claim about a population, size the panel for it.

## 5. Design decisions

**Why separate from design?** A different question deserves a different shape. Verification skips candidate generation and masking entirely, and its output is a **size matrix**, because the practical question about existing primers is "what bands, what sizes, in whom?" — not "which of 20 candidates is best?"

**Why sizes rather than pass/fail?** A biologist reads a gel by band size. A size matrix mirrors the experiment, and shows *how* a pair misbehaves — a 282 bp paralog band sitting just above the correct 280 bp product is a different problem from a missing band.

**Why accept GRCh38 coordinates?** That is how existing primers are almost always documented. Translation to the CHM13 ruler happens internally.

**Why one genome-wide search for all pairs?** The expensive step is reading a whole genome, not comparing primers against it — measured, one haplotype takes ~5.8 s for a single primer and ~6.1 s for 128. The cost is per *genome*, not per *pair*.

**Why keep "uncertain" distinct from "dropout"?** "No band because of a variant" and "we couldn't examine this locus" are different statements. Collapsing them makes a pair look better or worse than the evidence supports.

Terms used here are defined in [`../CONTEXT.md`](../CONTEXT.md).
