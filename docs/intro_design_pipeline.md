# **DESIGN** — an introduction to the primer-design pipeline

> Target region in, ranked primer pairs out. Run it with `pangenome-primer design`.
> Already have primers? That is **VERIFY** — see [`intro_verify_pipeline.md`](intro_verify_pipeline.md).

## 1. The problem

You want primers for some target region. The classic route: open the reference genome, pick primers, run a specificity check, order oligos.

The catch is that **the reference is one sequence and humans are not.** Any two people differ at millions of positions, so the reference is an average matching no real person. A primer that looks perfect on it can fail in the lab two ways:

- **Dropout.** A common SNP or indel sits under a primer, especially near its **3′ end**. The polymerase cannot extend off a mismatched 3′ end, so carriers show no band. The assay works on some samples and mysteriously fails on others.
- **Off-target amplification.** The reference collapses paralogs and pseudogenes into one copy. Your primer may also bind a near-identical gene elsewhere (CYP2D6's pseudogene CYP2D7 is the classic case), giving extra bands.

**DESIGN** vets primers against hundreds of real genomes at once, so both are caught before you spend bench time.

## 2. What is a pangenome?

A **pangenome** replaces the single consensus sequence with a **collection of complete genomes** from many people. This project uses **HPRC Release 2**: near-gapless assemblies from ~232 individuals across global populations.

Because long-read assembly separates ("phases") your two inherited copies of each chromosome, each person contributes **two haplotypes** — so ~232 people give **~460 real genomes** instead of one average.

**CHM13 v2.0**, a complete telomere-to-telomere reference, is used purely as a **coordinate backbone** — a common ruler for naming positions. The biology is judged on the haplotypes, never on CHM13.

> The real question about a primer is not "does it match the reference?" but "does it amplify the intended locus, cleanly, in as many actual people as possible?" A pangenome lets you ask exactly that.

A laptop-scale run uses a diverse **subset** (3–30 haplotypes spanning populations) — enough to expose diversity problems, small enough to finish.

## 3. How it works

```
target region
    │  (1) anchor on CHM13  ── CHM13 coords, GRCh38 coords, or a FASTA
    ▼
    │  (2) project onto each haplotype ── find the same locus in every genome
    ▼
    │  (3) mask variable sites ── positions that differ across people
    ▼
    │  (4) design candidates (Primer3) ── avoiding those sites
    ▼
    │  (5) search each primer genome-wide, in every haplotype
    ▼
    │  (6) decide which binding sites actually amplify (thermodynamic, 3′-aware)
    ▼
    │  (7) pair forward+reverse into predicted products (in-silico PCR)
    ▼
    │  (8) score each haplotype: pass / dropout / off-target / multi-product / uncertain
    ▼
    │  (9) rank pairs by coverage + specificity
    ▼
  report (TSV/JSON + colour-coded HTML)
```

**(1) Anchor.** Everything is normalised to one CHM13 anchor so later steps share a ruler.

**(2) Project.** The same gene sits at different coordinates in each person, because of upstream insertions and deletions. We locate the homologous window in every haplotype — the **expected on-target region**. If a locus is broken or unmappable in an assembly, that haplotype is **uncertain**, not a failure.

**(3) Mask.** Comparing projected windows across people reveals which positions are polymorphic. Primer3 is told to avoid them, weighting the 3′ zone most heavily. This steers primers toward conserved sequence by construction.

**(4) Generate candidates.** Primer3 handles Tm, GC, hairpins, dimers and product size, and returns ~20 pairs.

**(5) Search genome-wide.** A compiled scanner finds every place each primer could bind in each haplotype, allowing mismatches — not just at the target, but wherever a paralog or repeat lurks.

**(6) Decide what actually amplifies.** A binding site is not a product. A **nearest-neighbor thermodynamic model** asks whether each primer would stay annealed and, critically, whether its **3′ end** is matched enough to extend. A mismatch in the last couple of bases is a hard no.

**(7) Pair into products.** A product needs forward and reverse primers pointing at each other within a plausible distance. Surviving sites are combined into amplicons, each recorded with its **size** and labelled on- or off-target.

**(8) Score.** Each haplotype gets exactly one verdict per pair.

**(9) Rank.** Two independent numbers summarise across people — **on-target coverage** (does it amplify at all? the anti-dropout score) and **unique-product rate** (is it a single clean band? the specificity score). Pairs must clear hard filters, then are ordered so the best universal pair floats to the top.

## 4. Why it works this way

**Why a pangenome, not one reference?** The failures that matter are invisible on a single reference *by definition*. You only see a SNP-under-the-primer or a collapsed paralog by looking at many genomes. On CYP2D6 this pipeline flagged both a pseudogene co-amplification and pairs that fail on everyone because the reference happens to carry a rare allele.

**Why is CHM13 only a ruler?** It is complete and gap-free, which makes it a clean coordinate system — but it is still one genome, so it never decides whether a primer works.

**Why mask at design time *and* test at evaluation time?** Masking is cheap insurance, but it only knows about variation *among the haplotypes*. It cannot catch a position where the **reference itself is the odd one out** — everyone agrees with each other and differs from CHM13. Masking reduces the problem; evaluation proves it.

**Why thermodynamics instead of counting sequence matches?** Raw matches wildly overstate risk. A good primer can have hundreds of approximate genome-wide matches and still be perfectly specific, because almost all have a 3′ mismatch and won't extend, and almost none have a partner primer pointing back in range. A "many hits → reject" rule discards excellent primers.

**Why pair the primers?** One primer binding somewhere is not a product. Amplification needs **two**, **convergent**, **in range**. Scoring primers alone would flag risks that can never produce a band.

**Why keep "uncertain" separate from failure?** If a locus cannot be confidently located in an assembly, calling that a dropout is dishonest. Uncertain haplotypes leave the denominator and are reported separately, so coverage percentages mean what they say.

**Why two-stage, cheap-then-expensive?** The genome-wide search is the costly part, and checking *coverage* only needs the small projected window. So the cheap check runs on all candidates, and the expensive specificity search only on the top-ranked shortlist — same answer for the primers that matter, a fraction of the compute.

## 5. What you get

- A ranked **TSV/JSON** table of pairs with coverage, unique-product rate, dropout and off-target counts, and product sizes.
- A colour-coded **HTML report** with a per-haplotype status matrix, so you can see where a pair stumbles.

Terms used here are defined in [`../CONTEXT.md`](../CONTEXT.md).
