"""Correctness anchor (plan verification step 2): a synthetic mini-pangenome where each
haplotype is engineered to land in exactly one status. Runs in pure Python (rule mode),
no primer3/network needed.
"""
from __future__ import annotations

from pangenome_primer.binding import revcomp
from pangenome_primer.classify import RuleConfig, ThermoConfig
from pangenome_primer.engine import EvalConfig, HaplotypeContext, evaluate_pair
from pangenome_primer.model import HaplotypeStatus, Locus, Primer, PrimerPair
from pangenome_primer.rank import RankConfig, rank_pairs

# --- fixture building blocks -------------------------------------------------
F = "TTGCACAGTCCAGATTGCAA"        # forward primer (binds + strand)
R = "CATTGCGATTGACACTTGCG"        # reverse primer (binds - strand)
MID = "GATACCATGCTGACGTTAAC" * 5  # 100 bp spacer, low chance of spurious hits
LEFT = "AACCGGTTAACCGGTTAACC"
RIGHT = "TTGGCCAATTGGCCAATTGG"
SPACER = "ACGACTACGT" * 50        # 500 bp between paralog copies


def region() -> str:
    """One intact amplifiable copy: LEFT F MID revcomp(R) RIGHT."""
    return LEFT + F + MID + revcomp(R) + RIGHT


REGION = region()
F_START = len(LEFT)
PRODUCT = 20 + len(MID) + 20  # F..end of revcomp(R)


def mutate_forward_3p(seq: str, at: int) -> str:
    """Flip the 3' terminal base of the forward primer site starting at `at`."""
    pos = at + len(F) - 1
    orig = seq[pos]
    new = "C" if orig != "C" else "G"
    return seq[:pos] + new + seq[pos + 1 :]


def ctx(hid: str, seq: str, locus: Locus | None, ok: bool = True) -> HaplotypeContext:
    return HaplotypeContext(hid, "chr1", seq, locus, projection_ok=ok)


def full_locus(hid: str, seq: str) -> Locus:
    return Locus(hid, "chr1", 0, len(seq))


PAIR = PrimerPair("p1", Primer("p1_F", F), Primer("p1_R", R), product_size_chm13=PRODUCT)
CFG = EvalConfig(mode="rule", max_mismatches=3, min_product=80, max_product=200,
                 rule_cfg=RuleConfig(max_total_mm=2, max_3prime_mm=0, suffix=5))


def _status(contexts) -> dict[str, HaplotypeStatus]:
    res = evaluate_pair(PAIR, contexts, CFG)
    return {r.haplotype_id: r.status for r in res.per_haplotype}


def test_pass_perfect_match():
    seq = REGION
    st = _status([ctx("A", seq, full_locus("A", seq))])
    assert st["A"] is HaplotypeStatus.PASS


def test_dropout_on_3prime_snp():
    seq = mutate_forward_3p(REGION, F_START)
    st = _status([ctx("B", seq, full_locus("B", seq))])
    assert st["B"] is HaplotypeStatus.DROPOUT


def test_off_target_when_target_drops_but_paralog_intact():
    # region1 target has a 3' SNP (drops out); an intact paralog sits far downstream.
    target = mutate_forward_3p(REGION, F_START)
    seq = target + SPACER + REGION
    expected = Locus("C", "chr1", 0, len(target))  # only the (dead) target region
    st = _status([ctx("C", seq, expected)])
    assert st["C"] is HaplotypeStatus.OFF_TARGET


def test_multi_product_two_intact_copies():
    seq = REGION + SPACER + REGION
    expected = Locus("D", "chr1", 0, len(REGION))  # first copy is on-target
    st = _status([ctx("D", seq, expected)])
    assert st["D"] is HaplotypeStatus.MULTI_PRODUCT


def test_uncertain_when_projection_fails():
    st = _status([ctx("E", REGION, None, ok=False)])
    assert st["E"] is HaplotypeStatus.UNCERTAIN


def test_coverage_excludes_uncertain_from_denominator():
    contexts = [
        ctx("A", REGION, full_locus("A", REGION)),
        ctx("B", mutate_forward_3p(REGION, F_START), full_locus("B", REGION)),
        ctx("E", REGION, None, ok=False),
    ]
    res = evaluate_pair(PAIR, contexts, CFG)
    assert res.n_uncertain == 1
    assert len(res.evaluable) == 2
    assert res.on_target_coverage == 0.5  # 1 pass / 2 evaluable, not / 3


def test_thermo_mode_hard_gates_3prime_mismatch():
    """Thermo mode must still call a 3'-terminal SNP a dropout even though whole-duplex Tm
    barely moves — the non-negotiable 3' gate. Skips if primer3 is unavailable."""
    import pytest

    primer3 = pytest.importorskip("primer3")
    tm = primer3.calc_tm(F)
    thermo = ThermoConfig(design_tm=tm, tm_drop_max=5.0, three_prime_hard_nt=2)
    cfg = EvalConfig(mode="thermo", max_mismatches=3, min_product=80, max_product=200,
                     thermo_cfg=thermo)
    perfect = evaluate_pair(PAIR, [ctx("A", REGION, full_locus("A", REGION))], cfg)
    dead = mutate_forward_3p(REGION, F_START)
    dropped = evaluate_pair(PAIR, [ctx("B", dead, full_locus("B", REGION))], cfg)
    assert perfect.per_haplotype[0].status is HaplotypeStatus.PASS
    assert dropped.per_haplotype[0].status is HaplotypeStatus.DROPOUT


def test_coverage_and_uniqueness_diverge():
    """On-target coverage credits multi_product haplotypes (target amplifies, extra bands);
    unique_product_rate counts only clean single-product (pass). They must differ — the
    real GAPDH pair7 case: target amplifies everywhere but always with off-target bands."""
    from pangenome_primer import fixtures as fx

    res = evaluate_pair(fx.demo_pair(), fx.demo_contexts(), CFG)
    by = {r.haplotype_id: r.status for r in res.per_haplotype}
    assert by["HG_multi#hap1"] is HaplotypeStatus.MULTI_PRODUCT
    assert by["HG_off#hap1"] is HaplotypeStatus.OFF_TARGET
    # evaluable = pass, drop, off, multi (uncertain excluded) = 4
    assert len(res.evaluable) == 4
    # coverage: pass + multi (both have an on-target amplicon) = 2/4
    assert res.on_target_coverage == 0.5
    # uniqueness: only pass = 1/4
    assert res.unique_product_rate == 0.25
    assert res.on_target_coverage != res.unique_product_rate


def test_ranking_filters_and_orders():
    good = [ctx(f"h{i}", REGION, full_locus(f"h{i}", REGION)) for i in range(4)]
    bad = [ctx("h0", REGION, full_locus("h0", REGION))] + [
        ctx(f"h{i}", mutate_forward_3p(REGION, F_START), full_locus(f"h{i}", REGION))
        for i in range(1, 4)
    ]
    r_good = evaluate_pair(PrimerPair("good", PAIR.forward, PAIR.reverse), good, CFG)
    r_bad = evaluate_pair(PrimerPair("bad", PAIR.forward, PAIR.reverse), bad, CFG)
    ranked = rank_pairs([r_bad, r_good], RankConfig(min_coverage=0.95))
    assert ranked[0].result.pair.name == "good" and ranked[0].passed
    assert not ranked[1].passed  # bad pair rejected with reasons
    assert ranked[1].reject_reasons
