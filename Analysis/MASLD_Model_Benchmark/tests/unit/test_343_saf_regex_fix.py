"""Unit tests for the SAF regex fix in 343_harvest_documented_fstage.py.

The fix makes the deposited S (steatosis) and A (activity) grades survive
parsing. Previously the SAF pattern captured only the F group.

Every test here is written so that it FAILS against the pre-fix behaviour.
test_prepatch_pattern_would_fail_these makes that explicit rather than assumed.
"""
from __future__ import annotations
import gzip
import importlib.util
import os
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
H343 = ROOT/"Analysis/SingleCell/scripts/343_harvest_documented_fstage.py"
SERIES = (ROOT/"Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
              /"geo_cache/GSE202379_series_matrix.txt.gz")


@pytest.fixture(scope="module")
def h343():
    spec = importlib.util.spec_from_file_location("h343_under_test", H343)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _row(*vals):
    return pd.Series({f"char_{i}": v for i, v in enumerate(vals)} | {"sample": "GSMX"})


# ---------------------------------------------------------------- the fix
def test_saf_pattern_captures_three_groups(h343):
    assert h343.SAF_FULL_PATTERN.groups == 3
    m = h343.SAF_FULL_PATTERN.search("saf score: S2A3F1")
    assert m.groups() == ("2", "3", "1")


def test_extract_saf_returns_all_three_aspects(h343):
    assert h343.extract_saf(_row("tissue: Liver", "saf score: S2A3F1")) == (2, 3, 1)


def test_extract_saf_handles_S0_the_fourth_steatosis_level(h343):
    """Steatosis is S0-S3. A 1-3 assumption drops the S0 donor (P30)."""
    s, a, f = h343.extract_saf(_row("saf score: S0A1F0"))
    assert (s, a, f) == (0, 1, 0)
    assert s not in (1, 2, 3)


def test_extract_saf_returns_none_triple_when_absent(h343):
    assert h343.extract_saf(_row("disease status: Healthy control")) == (None, None, None)


# -------------------------------------------------- F value must not move
def test_extract_fstage_reads_the_LAST_group_not_the_first(h343):
    """With a 3-group pattern, group(1) is S. Reading it as F is the bug this
    guards against: S2A3F1 must give F=1, not 2."""
    assert h343.extract_fstage(_row("saf score: S2A3F1")) == 1


def test_extract_fstage_other_patterns_still_single_group(h343):
    assert h343.extract_fstage(_row("fibrosis stage: F3")) == 3
    assert h343.extract_fstage(_row("Kleiner fibrosis stage: 2")) == 2


# ------------------------------------------- the 2026-08-30 precedence fix
# A deposited SAF grade is a MEASUREMENT; "healthy control" / "end stage" are
# free-text CATEGORY LABELS. The measurement must win. These tests exist so a
# revert to label-first fails loudly: every assertion below returns the OLD
# label-derived value under the pre-fix source.
def test_saf_measurement_beats_cooccurring_healthy_control_label(h343):
    """THE guard. Pre-fix this returned 0 (the label). It must return 1."""
    row = _row("disease status: Healthy control", "saf score: S1A1F1")
    assert h343.extract_fstage(row) == 1
    # ... and independently of which column the label sits in
    assert h343.extract_fstage(_row("saf score: S1A1F1",
                                    "disease status: Healthy control")) == 1
    # the measurement itself is unchanged
    assert h343.extract_saf(row) == (1, 1, 1)


def test_saf_measurement_beats_cooccurring_end_stage_label(h343):
    """Pre-fix this returned 4 (the label). The deposited grade is F2."""
    assert h343.extract_fstage(_row("disease status: end stage",
                                    "saf score: S2A3F2")) == 2


def test_P98_and_P30_as_deposited(h343):
    """The two GSE202379 records that carry BOTH a label and a real SAF.
    P98 moves 0 -> 1. P30 hit the identical short-circuit and does not move,
    because its measured grade is F0 either way -- a disagreement-only audit
    would have missed it, so it is pinned explicitly."""
    p98 = _row("disease status: Healthy control", "saf score: S1A1F1")
    p30 = _row("disease status: Healthy control", "saf score: S0A1F0")
    assert h343.extract_fstage(p98) == 1        # pre-fix: 0
    assert h343.extract_fstage(p30) == 0        # pre-fix: 0 as well
    assert h343.extract_saf(p98) == (1, 1, 1)
    assert h343.extract_saf(p30) == (0, 1, 0)


def test_label_map_is_still_the_last_resort_when_nothing_was_measured(h343):
    """Demoting the label map must not delete it: the 7 GSE202379 donors with
    a label and NO SAF string still get their only available value."""
    assert h343.extract_fstage(_row("disease status: Healthy control",
                                    "saf score: healthy control")) == 0
    assert h343.extract_fstage(_row("disease status: end stage",
                                    "saf score: end stage")) == 4
    assert h343.extract_saf(_row("saf score: end stage")) == (None, None, None)


def test_f_stage_source_names_the_pass_that_produced_the_value(h343):
    """Provenance must stay readable, per pass."""
    assert h343.extract_fstage_with_source(
        _row("disease status: Healthy control",
             "saf score: S1A1F1")) == (1, "saf_regex")
    assert h343.extract_fstage_with_source(
        _row("fibrosis stage: F3")) == (3, "other_fstage_regex")
    assert h343.extract_fstage_with_source(
        _row("saf score: end stage")) == (4, "label_map:end stage")
    assert h343.extract_fstage_with_source(_row("tissue: liver")) == (None, "none")


def test_label_map_hit_still_reports_the_label(h343):
    assert h343.label_map_hit(_row("disease status: Healthy control")) == "healthy control"
    assert h343.label_map_hit(_row("disease status: NAFLD")) is None


# ------------------------------------------- the pre-fix behaviour fails these
def test_prepatch_pattern_would_fail_these(h343):
    """Construct the violation on purpose: the old single-group pattern
    cannot produce S or A at all."""
    import re
    old = re.compile(r"saf\s*score[:\s]*S\d+A\d+F(\d)", re.I)
    m = old.search("saf score: S2A3F1")
    assert m.groups() == ("1",), "the old pattern captured F only"
    assert old.groups == 1
    with pytest.raises(IndexError):
        m.group(2)


# ---------------------------------------------------- end-to-end, offline
@pytest.mark.skipif(not SERIES.exists(), reason="cached series matrix absent")
def test_harvest_emits_saf_columns_offline(h343, monkeypatch):
    """Exercise the harvest emit branch without touching the network, and
    without writing the frozen output file."""
    monkeypatch.setattr(h343, "GEO_URLS", {"GSE202379": "unused"})
    monkeypatch.setattr(h343, "download_geo", lambda acc, url: [SERIES])
    monkeypatch.setattr(h343, "parse_geo_bioproject", lambda p: None)
    monkeypatch.setattr(h343, "resolve_gsm_to_srr", lambda prj: {})

    df = h343.harvest_geo()
    for c in ("saf_S", "saf_A", "saf_F", "f_stage_source"):
        assert c in df.columns, f"{c} missing from the harvest output"

    graded = df[df["saf_F"].notna()]
    assert len(graded) == 40, f"expected 40 SAF-graded GSMs, got {len(graded)}"
    assert sorted(graded["saf_S"].astype(int).unique()) == [0, 1, 2, 3]

    # provenance vocabulary after the precedence fix: the label map can no
    # longer override a measurement, and the two records where it used to
    # (P30, P98) are now marked as measurement-over-demoted-label.
    srcs = set(df["f_stage_source"])
    assert not any("OVERRODE" in s for s in srcs), \
        "a label-map override reappeared; precedence has been reverted"
    assert any(s.startswith("saf_regex_over_label_map:") for s in srcs), srcs
    assert any(s.startswith("label_map:") and "no_saf_string" in s for s in srcs), srcs
    assert "saf_regex" in srcs

    # exactly the two GSMs that carry a label AND a measurement
    over = df[df["f_stage_source"].str.startswith("saf_regex_over_label_map:")]
    assert sorted(over["sample"]) == ["GSM6112243", "GSM6112244"], sorted(over["sample"])
    assert sorted(over["F_stage_documented"]) == [0, 1]      # P30 F0, P98 F1

    # a zero recovery must never look like a legitimate result
    assert graded["saf_A"].notna().all()
    assert graded["saf_S"].notna().all()
