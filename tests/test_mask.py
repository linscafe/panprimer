"""Tests for the variability mask's region coalescing and Primer3 cap (the CYP2D6 case:
hundreds of variant sites must not overflow SEQUENCE_EXCLUDED_REGION)."""
from __future__ import annotations

from pangenome_primer import mask


def test_merge_coalesces_within_gap():
    # positions 0,1 (adjacent) and 5 (gap 3 from 1) -> merged at gap>=3; 100 stays separate
    assert mask._merge([0, 1, 5, 100], gap=0) == [(0, 2), (5, 1), (100, 1)]
    assert mask._merge([0, 1, 5, 100], gap=3) == [(0, 6), (100, 1)]


def test_build_excluded_caps_and_merges(monkeypatch):
    # 500 isolated variant sites, spaced 10 bp apart -> would be 500 regions.
    positions = list(range(0, 5000, 10))
    counts = [0] * 5000
    for p in positions:
        counts[p] = 3
    monkeypatch.setattr(mask, "variability_counts", lambda *a, **k: (counts, 3))

    # gap 25 coalesces the 10-bp-spaced sites into one long zone
    regions = mask.build_excluded_regions("x" * 5000, ["y"], min_allele_freq=0.5,
                                          merge_gap=25, max_regions=200)
    assert len(regions) == 1

    # with no coalescing, 500 regions must be capped to Primer3's limit
    capped = mask.build_excluded_regions("x" * 5000, ["y"], min_allele_freq=0.5,
                                         merge_gap=0, max_regions=200)
    assert len(capped) <= 200
