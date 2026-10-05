"""Known-unit checks for Resource-specific identity and source-label handling."""
import importlib.util
from pathlib import Path

import pandas as pd

spec = importlib.util.spec_from_file_location("resource_identity", Path(__file__).parents[1] / "11_resource_identity.py")
identity = importlib.util.module_from_spec(spec)
spec.loader.exec_module(identity)


def base_roster():
    return pd.DataFrame({
        "sample_id": ["a", "b", "c"], "dataset": ["test"] * 3,
        "group_binary": ["Disease"] * 3, "source_fibrosis_category": ["F0–F1", "F0–F1", "F2"],
        "nas_score": [3, 3, 3], "reported_sex": ["F"] * 3, "diagnosis_harmonized": ["MASH"] * 3,
        "source_person_id": [None] * 3, "first_biopsy": ["False", "True", "unknown"],
        "full_called": [2000, 1100, 1800], "filtered_component": ["a", "a", "c"],
        "identity_qualified": [True] * 3, "exclusion_reason": [""] * 3})


def test_full_histology_conflict_is_not_hidden_by_coarse_stage_code():
    assert identity.source_conflicts(base_roster()) == ["source_fibrosis_category"]


def test_source_known_repeat_biopsy_may_change_histology():
    roster = base_roster()
    roster["dataset"] = "GSE193066"
    roster["source_person_id"] = "source_person"
    assert identity.source_conflicts(roster) == []
    roster.loc[2, "reported_sex"] = "M"
    assert identity.source_conflicts(roster) == ["reported_sex"]


def test_representative_is_first_biopsy_and_is_inside_actual_context():
    roster = base_roster()
    selected, ids = identity.select_context(roster, ["a", "b"], "both")
    assert ids.sample_id.tolist() == ["b"]
    assert selected.loc[selected.primary_selected, "unit_choice"].tolist() == ["documented_first_biopsy"]
    _, ids = identity.select_context(roster, ["a"], "one")
    assert ids.sample_id.tolist() == ["a"]


def test_unresolved_identity_is_not_rescued_by_selecting_another_library():
    roster = base_roster()
    roster.loc[:1, "identity_qualified"] = False
    _, ids = identity.select_context(roster, ["a", "b", "c"], "test")
    assert ids.sample_id.tolist() == ["c"]


def test_pair_abstains_on_low_coverage_or_discordant_independent_partition():
    original = pd.DataFrame({"#IID1": ["a"], "IID2": ["b"], "relation": ["same_individual"]})
    table = pd.DataFrame({"#IID1": ["a"], "IID2": ["b"], "NSNP": [600], "KINSHIP": [.49], "IBS0": [0]})
    filtered = {label: table.copy() for label in ["full", "half0", "half1"]}
    counts = pd.DataFrame({"sample_id": ["a", "b"], "full_called": [1200, 1200]})
    _, stable, unresolved = identity.qualify_pairs(original, filtered, counts)
    assert stable == {("a", "b")} and not unresolved
    filtered["half1"].loc[0, "KINSHIP"] = .25
    _, stable, unresolved = identity.qualify_pairs(original, filtered, counts)
    assert not stable and unresolved == {("a", "b")}
    filtered["half1"].loc[0, "KINSHIP"] = .49
    filtered["half1"].loc[0, "NSNP"] = 249
    _, stable, unresolved = identity.qualify_pairs(original, filtered, counts)
    assert not stable and unresolved == {("a", "b")}


def test_new_same_identity_edge_with_tiny_overlap_cannot_join_a_graph():
    original = pd.DataFrame(columns=["#IID1", "IID2", "relation"])
    table = pd.DataFrame({"#IID1": ["a"], "IID2": ["b"], "NSNP": [5], "KINSHIP": [.5], "IBS0": [0]})
    filtered = {label: table.copy() for label in ["full", "half0", "half1"]}
    counts = pd.DataFrame({"sample_id": ["a", "b"], "full_called": [1200, 1200]})
    checks, stable, unresolved = identity.qualify_pairs(original, filtered, counts)
    assert checks.empty and not stable and not unresolved


def test_context_rejects_missing_sample_ids():
    try:
        identity.select_context(base_roster(), ["missing"], "test")
    except ValueError as error:
        assert "Unjoined" in str(error)
    else:
        raise AssertionError("Missing sample was silently accepted")
