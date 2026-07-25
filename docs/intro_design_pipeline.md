# An introduction to the primer **design** pipeline

---

## 1. What problem does this solve?

You want to design a PCR primer pair that amplifies some target region — a gene, an exon, a variant of interest. The classic way is to open the human **reference genome** (GRCh38, or the newer CHM13), pick primers in your region, run a quick specificity check, and order the oligos.

The catch: **the reference genome is one sequence, but humans are not.** Any two people differ at millions of positions (SNPs, insertions/deletions, larger structural changes). The reference is essentially a genetic "average" that matches no single real person exactly. So a primer that looks perfect on the reference can quietly fail in the lab, in two classic ways:

- **Dropout.** A common SNP or indel sits **under a primer — especially near its 3′ end.** The polymerase can't extend off a mismatched 3′ end, so in people carrying that variant the band simply doesn't appear. Your assay works on some samples and mysteriously fails on others.
- **Off-target amplification.** The reference collapses **paralogs, pseudogenes, and duplicated regions** into one copy. Your primer may also bind a near-identical gene elsewhere (e.g. CYP2D6's pseudogene CYP2D7), giving extra bands or amplifying the wrong locus.

This pipeline designs **"universal" primers** — primers chosen and vetted against the genetic diversity of *hundreds of real human genomes at once*, so dropout and off-target problems are caught **before** you spend money and bench time.

---

## 2. First — what is a pangenome?

A **reference genome** is a single consensus sequence (~3 billion bases) used as *the* human genome. It's convenient but lossy: it hides variation and merges repeated regions.

A **pangenome** replaces that one sequence with a **collection of many complete genomes** from different people. This project uses the **Human Pangenome Reference Consortium, Release 2 (HPRC R2)**: high-quality, near-gapless genome assemblies from **~232 individuals sampled across global populations** (African, European, East Asian, South Asian, American ancestries).

Two terms you'll see:

- **Haplotype.** You inherit one copy of each chromosome from each parent. Modern long-read assembly can separate ("phase") those two copies into two independent sequences. So each person contributes **two haplotypes**, and ~232 people give **~460 haplotype genomes** — 460 real, complete versions of the human genome instead of one average.
- **CHM13 v2.0.** A single, **complete telomere-to-telomere** reference (no gaps, including hard regions GRCh38 leaves blank). We use it purely as a **coordinate backbone** — a common ruler to name positions — while the *biology* is judged on the 460 haplotypes.

> **Why this matters for PCR:** a primer's real question isn't "does it match the reference?" but "does it amplify the intended locus, cleanly and uniquely, in as many *actual people* as possible?" A pangenome lets us ask exactly that.

For a laptop-scale demo we use a **diverse subset** (e.g. 3–10 haplotypes spanning populations) rather than all 460 — enough to expose diversity problems, small enough to run.

---

## 3. How the pipeline works, step by step

```
target region
    │  (1) anchor on CHM13  ── accepts CHM13 coords, GRCh38 coords, or a FASTA sequence
    ▼
    │  (2) project onto each haplotype ── find the same locus in every person's genome
    ▼
    │  (3) mask variable sites ── mark positions that differ across people
    ▼
    │  (4) design candidates (Primer3) ── on the reference template, avoiding variable sites
    ▼
    │  (5) search each primer genome-wide in every haplotype
    ▼
    │  (6) decide which binding sites actually work (thermodynamics, 3′-aware)
    ▼
    │  (7) pair forward+reverse into predicted PCR products (in-silico PCR)
    ▼
    │  (8) score each haplotype: pass / dropout / off-target / multi-product / uncertain
    ▼
    │  (9) rank primer pairs by coverage + specificity
    ▼
  report (table + colour-coded HTML)
```

**(1) Anchor the target.** You give a region — as CHM13 coordinates, GRCh38 coordinates (we fetch that slice from the GRCh38 reference and align it to CHM13), or just a FASTA sequence of the locus. Everything is normalized to one CHM13 anchor so later steps share a common ruler.

**(2) Project onto each haplotype.** The "same" gene sits at slightly different coordinates in each person (because of upstream insertions/deletions). We locate the homologous window in every haplotype — this is the **expected on-target region** for that person. If a locus is broken or unmappable in someone's assembly, that haplotype is marked **uncertain** (not a failure — see §4).

**(3) Mask variable sites.** Comparing the projected windows across people reveals which positions are polymorphic. We tell the primer designer to **avoid placing primers over those positions**, weighting the 3′-end zone most heavily — because that's where a mismatch is most damaging. In effect we steer primers toward the **conserved** parts of the locus.

**(4) Generate candidates.** We run **Primer3** (the standard primer-design engine — it handles Tm, GC content, hairpins, dimers, product size) on the reference template, with the variable sites excluded. Out come ~20 candidate forward/reverse pairs.

**(5) Search genome-wide.** For each candidate we ask: *where in each haplotype could these primers bind?* A built-in scanner reads the whole genome and finds every binding site, allowing a few mismatches — not just at the target, but everywhere a paralog or repeat might lurk.

**(6) Decide which sites actually amplify.** A binding site in a table is not the same as amplification. Using a **thermodynamic (nearest-neighbor) model** — the same physics that predicts melting temperature from sequence, salt, and Mg²⁺ — we ask whether each primer would *really* stay annealed and, critically, whether its **3′ end** is matched enough to extend. A mismatch in the last couple of bases is a hard "no."

**(7) Pair into products (in-silico PCR).** A PCR product needs a **forward and a reverse primer pointing at each other, within a plausible distance.** We combine the surviving binding sites into predicted amplicons, record each product's **size**, and label it **on-target** (inside the expected window) or **off-target** (anywhere else).

**(8) Score each haplotype.** For a given primer pair, each person's genome gets exactly one verdict: **pass** (one clean on-target product), **dropout** (target should amplify but a binding-site variant kills it), **off-target**, **multi-product** (extra bands), or **uncertain**.

**(9) Rank the pairs.** We summarize across people with two independent numbers:
- **On-target coverage** — in what fraction of people does the target amplify at all? (the anti-dropout score)
- **Unique-product rate** — in what fraction is it a single, clean, on-target band? (the specificity score)

Pairs must clear hard filters (e.g. ≥95 % coverage, zero off-target) and are then ordered so the best universal primer floats to the top.

---

## 4. Why it's built this way (and not another)

**Why a pangenome instead of one reference?** Because the failures we care about — dropout and paralog off-target — are **invisible on a single reference by definition.** You only see a SNP-under-the-primer or a collapsed paralog when you look at many real genomes. On a real hard locus (CYP2D6) this pipeline flagged both a pseudogene co-amplification *and* primer pairs that fail on everyone because the reference happened to carry a rare allele — neither is detectable from the reference alone.

**Why use CHM13 only as a coordinate backbone, not as truth?** CHM13 is complete and gap-free, which makes it a clean ruler. But it's still *one* genome, so we never let it decide whether a primer works — that verdict always comes from the 460 haplotypes.

**Why *both* mask at design time *and* test at evaluation time?** Masking (step 3) steers primers to conserved regions *by construction* — cheap insurance. But masking only knows about variation *among the haplotypes*; it can't catch a position where the **reference itself is the odd one out** (everyone agrees with each other but differs from CHM13). Only the full evaluation (steps 5–8) catches that. Design reduces the problem; evaluation proves it.

**Why a thermodynamic, 3′-aware test — instead of a simple exact-match / BLAST count?** Because raw sequence matches badly overstate risk. A good primer can have *hundreds* of approximate matches genome-wide, yet be perfectly specific — because almost all of them have a 3′ mismatch and won't extend, and almost none have a partner primer pointing back at them in range. A naïve "many hits → reject" rule throws away excellent primers. Modeling the actual chemistry (3′ extension + Tm) is what separates a real product from a harmless near-match.

**Why in-silico PCR (pair the primers) instead of scoring each primer alone?** A single primer binding somewhere is not a product. Amplification requires **two** primers, **convergent**, **within range.** Scoring primers individually would flag risks that can never actually produce a band.

**Why keep "uncertain" separate from "failure"?** If a locus can't be confidently located in someone's assembly, calling that a dropout would be dishonest. We exclude uncertain haplotypes from the denominator and report them separately, so coverage percentages mean what they say.

**Why the two-stage (cheap-then-expensive) design?** The genome-wide search (steps 5–6) is the costly part. Checking *coverage* only needs the small projected window, so we do that cheap check on **all** candidates first, then run the expensive genome-wide specificity search only on the **top-ranked shortlist.** Same answer for the primers that matter, a fraction of the compute.

---

## 5. What you get

- A ranked table (**TSV/JSON**) of primer pairs with coverage, unique-product rate, dropout and off-target counts, and product sizes.
- A **colour-coded HTML report** with a per-pair × per-haplotype status matrix, so you can *see* which populations a pair covers and where it stumbles.

The top pair is one you can order with real confidence that it will behave across human diversity — not just against an average that no patient actually is.

---

## Mini-glossary

- **Reference genome** — a single consensus human sequence (GRCh38 / CHM13).
- **Pangenome** — many complete genomes together (here HPRC R2, ~460 haplotypes).
- **Haplotype** — one phased copy of a person's genome (each person has two).
- **CHM13 v2.0** — a complete, gap-free reference used here as the coordinate backbone.
- **Dropout** — no product because a variant sits under a primer (usually near the 3′ end).
- **Off-target** — an extra/incorrect product from a paralog, pseudogene, or repeat.
- **On-target coverage** — fraction of people in whom the target amplifies at all.
- **Unique-product rate** — fraction in whom it's a single clean band (specificity).
- **Primer3** — the standard tool for primer candidate design (Tm, GC, hairpins, dimers).
- **In-silico PCR** — predicting products by pairing forward+reverse binding sites.
