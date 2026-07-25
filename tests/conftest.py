"""Shared pytest configuration for this suite (Phase 0 safety net).

Registers the `slow` marker and skips `slow`-marked tests by default: they need real HPRC
genome data under `hprc-r2/` and take ~9 minutes, so they must not run in the default
`pytest -q` pass. Run them explicitly with `pytest -m slow --runslow` (or select the file
directly). Lives in tests/ rather than pyproject.toml so it doesn't touch config a parallel
agent may also be editing.
"""
from __future__ import annotations

import pytest


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers",
        "slow: needs real HPRC genome data and takes minutes; skipped by default, "
        "run with `pytest -m slow --runslow`.",
    )


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run tests marked 'slow'"
    )


def pytest_collection_modifyitems(config: pytest.Config, items: list) -> None:
    if config.getoption("--runslow"):
        return
    markexpr = config.getoption("-m") or ""
    if "slow" in markexpr:
        return  # user explicitly selected slow tests via -m; let normal filtering apply
    skip_slow = pytest.mark.skip(reason="needs --runslow (real HPRC genome data, ~9 min)")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
