# An introduction to the primer **verification** pipeline

> New to pangenomes? Read §2 below for the one-minute version, or [`intro_design_pipeline.md`](intro_design_pipeline.md) for the fuller explanation. The two pipelines share the same engine but answer opposite questions.

---

## 1. What problem does this solve?

Sometimes you don't need to *design* primers — you already **have** a pair. Maybe from a published paper, an inherited protocol, a commercial kit, or a colleague's notebook. The question is different:

> *"Will these particular primers behave across the diversity of real human genomes — or will they drop out in some people, or throw extra bands from a look-alike locus?"*

The verification pipeline answers exactly that. You hand it a list of primer pairs; it predicts, **for each pair in each of many real human genomes**, what you would see on a gel: the size of the correct product, and the sizes of any off-target products. **It designs nothing** — it only evaluates what you give it.

Two things a single-reference check can't tell you, but this can:

- **Does the pair drop out in some people?** A common variant under a primer's 3′ end means no band in carriers — an assay that silently works on some samples and fails on others.
- **Does it amplify a look-alike locus?** Paralogs and pseudogenes (e.g. CYP2D6 vs its pseudogene CYP2D7) can give extra bands or the wrong product — invisible on a reference that collapses them into one copy.

---

## 2. Pangenome in one minute

The usual **reference genome** is a *single* consensus sequence — a genetic "average" that matches no real person exactly, and that merges repeated regions into one copy.

A **pangenome** is instead a **collection of many complete real genomes.** We use **HPRC Release 2**: near-gapless assemblies from **~232 people across global populations.** Because each person's two chromosome copies are separated ("phased") into two **haplotypes**, that's up to **~460 real human genomes** to test against, instead of one.

**CHM13 v2.0** — a single, complete, gap-free genome — is used only as a **coordinate ruler** so every genome can be lined up; the actual verdicts come from the many haplotypes. (A laptop-scale run uses a diverse **subset** of haplotypes spanning populations.)

Everything below reduces to one idea: **test the primers you have against many real people, not against one average.**

---

## 3. How the pipeline works, step by step

**Input:** a CSV, one row per primer pair, four columns:

| column | meaning |
|---|---|
| `primer_id` | a label you choose |
| `target` | the **intended amplicon region** in **GRCh38** coordinates (where the correct product should land) — most published primers are annotated in GRCh38 |
| `forward` | forward primer sequence, 5′→3′ |
| `reverse` | reverse primer sequence, 5′→3′ |

```
primers.csv
    │  (1) locate each target on CHM13 ── from GRCh38 coordinates
    ▼
    │  (2) project each target onto every haplotype ── the expected "correct" window per person
    ▼
    │  (3) search all primers genome-wide in each haplotype (bwa) ── one pass covers every pair
    ▼
    │  (4) keep binding sites that would truly amplify (thermodynamics, 3′-aware)
    ▼
    │  (5) pair forward+reverse into predicted products; record size + on/off-target
    ▼
  matrix report: rows = primer pairs, columns = haplotypes
```

**(1) Locate the target.** The `target` region is given in GRCh38 coordinates. We fetch that slice from the GRCh38 reference, align it to CHM13, and record where the **intended amplicon** lives on the common ruler. (You can also give CHM13 coordinates directly.)

**(2) Project onto each haplotype.** For every person's genome we find the homologous window — the place where the *correct* product is expected. Anything the primers produce **there** is **on-target**; anything **elsewhere** is **off-target.** If the locus can't be confidently found in someone's assembly, that cell is **uncertain** (`?`).

**(3) Search genome-wide.** Using **bwa**, we find everywhere each primer could bind in each haplotype — target and look-alikes alike. One genome-wide pass per haplotype covers **all** your primer pairs at once, so checking many pairs costs little more than checking one.

**(4) Keep sites that really amplify.** As in the design pipeline, a binding site in a table isn't a product. A **thermodynamic, 3′-end-aware** model decides whether each primer would actually stay annealed and extend, filtering out the many harmless near-matches.

**(5) Pair into products.** We combine surviving forward/reverse sites that point at each other within a plausible distance (default ≤ 2 kb) into predicted amplicons, and record each product's **size** and whether it's **on-** or **off-target.**

---

## 4. Reading the output

The primary deliverable is a **matrix** — one **row per primer pair**, one **column per haplotype** — where each cell shows the predicted PCR product size(s):

- **Green** number = the **correct (on-target) amplicon** size. Sizes may differ by a few bases between people — that's a real biological signal (a small insertion/deletion inside the amplicon), not an error; large deviations are underlined.
- **Red** number = an **off-target** product size (a spurious band).
- Grey **dropout** = the primers do **not** make the correct product in that person (a binding-site variant killed it).
- **`?`** = the locus couldn't be projected in that assembly (unknown, not a failure).

Read a row left-to-right and you effectively see the gel you'd get across many people: all-green means a robust pair; a red here or a "dropout" there tells you precisely where — and in which population — the pair would let you down. A machine-readable **TSV/JSON** is written alongside.

---

## 5. Why it's built this way (and not another)

**Why a separate pipeline from design?** It's a different *question*, so it deserves a different *shape.* Verification skips primer generation and variable-site masking entirely — there's nothing to design. And its output is a **size matrix**, because the practical question a bench scientist asks about existing primers is "what bands, what sizes, in whom?" — not "which of 20 candidates is best?"

**Why report product sizes, not just pass/fail?** Because a biologist reads a gel by **band size.** A size matrix mirrors the actual experiment: a clean single green band, an extra red band, or a missing band map directly onto what you'd load and image. Encoding on/off-target as sizes (not just a verdict) also reveals *how* a pair misbehaves — e.g. a ~282 bp paralog band sitting just above the correct 280 bp product.

**Why accept GRCh38 coordinates for the target?** Because that's how existing primers are almost always documented. Meeting people where their data already lives lowers the barrier; we translate to the CHM13 ruler internally.

**Why one genome-wide search for all pairs at once?** The expensive step is loading a whole genome's index into memory. Doing it **once per haplotype** and screening every primer against it makes the cost depend on the number of *genomes*, not the number of *primer pairs* — so screening a whole panel is nearly as cheap as screening one pair.

**Why the thermodynamic, in-silico-PCR core (shared with the design pipeline)?** Same reasoning as design: a primer with many approximate genome-wide matches can still be perfectly specific, because a mismatched 3′ end won't extend and most "hits" have no partner primer pointing back. Predicting *products* from the actual chemistry — rather than counting raw sequence matches — is what makes the off-target calls trustworthy.

**Why keep "uncertain" distinct from "dropout"?** "No correct band because of a variant" and "we couldn't confidently examine this locus in this assembly" are different statements. Collapsing them would make a pair look worse (or better) than the evidence supports, so we keep them apart.

---

## Mini-glossary

- **Pangenome / haplotype / CHM13** — see §2 (many real genomes; one phased genome copy; the coordinate ruler).
- **On-target product** — the intended amplicon, at the expected locus (shown green).
- **Off-target product** — a spurious amplicon from a paralog/pseudogene/repeat (shown red).
- **Dropout** — no correct product because a variant sits under a primer (usually near the 3′ end).
- **In-silico PCR** — predicting products by pairing forward+reverse binding sites in silico.
- **bwa** — a fast genome-wide sequence search tool, used to find primer binding sites.
