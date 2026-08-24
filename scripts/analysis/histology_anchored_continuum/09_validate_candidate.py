#!/usr/bin/env python3
"""Fail-closed validator for a histology-anchored continuum candidate."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import pathlib
import sys
from collections import defaultdict
from typing import Callable, Iterable


PROJECT = pathlib.Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
).resolve()
OUT_VALUE = os.environ.get("HAC_OUT_ROOT", "")
if not OUT_VALUE:
    raise SystemExit("HAC_OUT_ROOT is unset")
OUT = pathlib.Path(OUT_VALUE).resolve()
EXPECTED_PARENT = (PROJECT / "RNA-seq/results/histology_anchored_continuum/candidates").resolve()
if EXPECTED_PARENT not in OUT.parents:
    raise SystemExit(f"Candidate is outside the isolated result namespace: {OUT}")

VALIDATION = OUT / "validation"
if VALIDATION.exists():
    raise SystemExit(f"Refusing to overwrite validation namespace: {VALIDATION}")
VALIDATION.mkdir(parents=True)

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
with (SCRIPT_DIR / "00_prespecification.json").open() as handle:
    PRESPEC = json.load(handle)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists():
        raise AssertionError(f"missing file: {path}")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value: str) -> float:
    if value in ("", "NA", "NaN", "nan", None):
        return math.nan
    return float(value)


def as_bool(value: str) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"true", "t", "1"}:
        return True
    if normalized in {"false", "f", "0"}:
        return False
    raise AssertionError(f"not a logical value: {value!r}")


def unique(rows: Iterable[dict[str, str]], columns: tuple[str, ...]) -> set[tuple[str, ...]]:
    return {tuple(row[column] for column in columns) for row in rows}


def assert_columns(rows: list[dict[str, str]], columns: Iterable[str], label: str) -> None:
    if not rows:
        raise AssertionError(f"{label} has no rows")
    absent = set(columns) - set(rows[0])
    if absent:
        raise AssertionError(f"{label} lacks columns: {sorted(absent)}")


CHECKS: list[dict[str, str]] = []


def check(check_id: str, fn: Callable[[], str]) -> None:
    try:
        detail = fn()
        CHECKS.append({"check_id": check_id, "status": "PASS", "detail": str(detail)})
    except Exception as error:  # retain every failure in one audit table
        CHECKS.append(
            {
                "check_id": check_id,
                "status": "FAIL",
                "detail": f"{type(error).__name__}: {error}",
            }
        )


def check_reproduction() -> str:
    summary = read_tsv(OUT / "reproduction/reproduction_summary.tsv")
    by_metric = {row["metric"]: row for row in summary}
    rho = abs(as_float(by_metric["released_pc1_spearman"]["value"]))
    assert rho >= PRESPEC["promotion_gates"]["reproduction_abs_spearman_minimum"]
    assert by_metric["vendor_commit"]["value"] == PRESPEC["competitor_release"]["commit"]
    assert as_float(by_metric["n_discovery_participants"]["value"]) == 135
    multiplier = as_float(by_metric["operational_pc1_multiplier"]["value"])
    assert multiplier == PRESPEC["orientation"]["operational_continuum_multiplier"] == -1
    assert as_float(by_metric["operational_pc1_negation_lines"]["value"]) >= 1
    assert as_float(by_metric["operational_rf_inverse_pc1_lines"]["value"]) >= 1
    assert as_bool(by_metric["released_pc1_file_strictly_descending"]["value"])
    assert as_float(
        by_metric["operational_percentile_file_order_max_abs_difference"]["value"]
    ) == 0
    scores = read_tsv(OUT / "reproduction/participant_scores.tsv")
    assert len(scores) == 135
    for row in scores:
        raw_reproduced = as_float(row["pc1_reproduced_raw"])
        raw_released = as_float(row["pc1_released_raw"])
        assert math.isclose(
            as_float(row["pc1_reproduced_operational"]), -raw_reproduced,
            rel_tol=0, abs_tol=1e-12,
        )
        assert math.isclose(
            as_float(row["pc1_released_operational"]), -raw_released,
            rel_tol=0, abs_tol=1e-12,
        )
    ordered = sorted(scores, key=lambda row: int(row["released_file_position"]))
    raw_in_file_order = [as_float(row["pc1_released_raw"]) for row in ordered]
    assert all(left > right for left, right in zip(raw_in_file_order, raw_in_file_order[1:]))
    denominator = len(ordered) - 1
    for index, row in enumerate(ordered):
        assert math.isclose(
            as_float(row["released_operational_rank_percentile"]),
            index / denominator, rel_tol=0, abs_tol=1e-12,
        )
    return (
        f"abs Spearman rho={rho:.8f}; operational multiplier={multiplier:.0f}; "
        f"commit={by_metric['vendor_commit']['value']}"
    )


def check_signature_mapping() -> str:
    signature = read_tsv(OUT / "reproduction/signature_genes.tsv")
    missing = read_tsv(OUT / "reproduction/missing_signature_genes.tsv")
    mapping = read_tsv(OUT / "programs/signature_gencode_v49_mapping.tsv")
    assert len(signature) == 145 and len(unique(signature, ("gene_id_base",))) == 145
    assert sum(as_bool(row["in_resource"]) for row in signature) == 139
    assert len(missing) == 6 and all(not as_bool(row["in_resource"]) for row in missing)
    assert len(mapping) == 145 and all(as_bool(row["symbol_concordant"]) for row in mapping)
    return "145 mapped by Ensembl base ID; 139 present; six explicitly missing"


def check_score_registry() -> str:
    unsup = read_tsv(OUT / "unsupervised/participant_scores.tsv")
    projection = read_tsv(OUT / "projection/participant_scores.tsv")
    native = read_tsv(OUT / "programs/native_nmf_participant_scores.tsv")
    assert len(unsup) == 2 * 844
    assert len(projection) == 3 * 844
    assert len(native) == 2 * 844
    assert {row["axis_id"] for row in native} == {"nmf_k4", "nmf_k6"}
    combined = unsup + projection
    assert len(unique(combined, ("sample_id", "axis_id"))) == 5 * 844
    expected_axes = set(
        PRESPEC["axes"]["co_primary"]
        + PRESPEC["axes"]["secondary"]
        + PRESPEC["axes"]["display_only"]
    )
    assert {row["axis_id"] for row in combined} == expected_axes
    for row in combined:
        percentile = as_float(row["axis_percentile"])
        assert math.isfinite(percentile) and 0 <= percentile <= 1
    return "844 donors have five bounded, nonduplicated scores"


def check_coverage() -> str:
    rows = read_tsv(OUT / "unsupervised/coverage.tsv")
    assert len(rows) == 5
    minimum = PRESPEC["signature"]["minimum_coverage_fraction"]
    observed = [as_float(row["coverage_fraction"]) for row in rows]
    assert all(value >= minimum for value in observed)
    return f"five cohorts; minimum signature coverage={min(observed):.4f}"


def check_rf_isolation() -> str:
    manifest = read_tsv(OUT / "projection/rf_training_manifest.tsv")
    assert len(manifest) == 1
    row = manifest[0]
    expected = PRESPEC["random_forest"]
    assert row["training_scope"] == "published_discovery_135_only"
    assert int(row["n_training_donors"]) == 135
    assert int(row["num_trees"]) == expected["num_trees"]
    assert int(row["mtry"]) == expected["mtry"]
    assert int(row["min_node_size"]) == expected["min_node_size"]
    assert int(row["seed"]) == PRESPEC["seeds"]["rf"]
    assert int(row["n_predictors"]) == 139
    assert row["predictor_freeze_scope"] == "published_discovery_signature_only"
    assert row["missing_target_predictor_rule"] == "centered_value_zero_equals_discovery_center"
    assert not as_bool(row["target_histology_used_for_training"])
    overlap = read_tsv(OUT / "projection/rf_training_overlap_audit.tsv")
    independent = [r for r in overlap if r["evaluation_role"] == "independent_validation"]
    assert len(independent) == 2
    assert all(int(r["n_exact_sample_id_overlap_with_discovery"]) == 0 for r in independent)
    return "ranger frozen at 1000/11/5; 135 discovery donors; independent overlap=0"


def check_program_exclusion_and_testability() -> str:
    rows = read_tsv(OUT / "programs/program_testability.tsv")
    family = PRESPEC["programs"]["family_size"]
    assert len(rows) == family * 5
    assert len(unique(rows, ("dataset", "program_uid"))) == family * 5
    assert all("testable" in row for row in rows)
    assert all("observed_l1_fraction" in row for row in rows)
    focal = set(PRESPEC["programs"]["focal_programs"].values())
    focal_rows = [row for row in rows if row["program_uid"] in focal]
    assert len(focal_rows) == 2 * 5
    assert all(as_bool(row["testable"]) for row in focal_rows)
    assert all(
        as_float(row["retained_l1_fraction"])
        >= PRESPEC["programs"]["minimum_retained_l1_weight"]
        for row in focal_rows
    )
    expected_retention = {
        "hotspot_hepatocytes_48f39dd4d817a10e": 0.957397465312792,
        "hotspot_hepatocytes_f05c535ae5bbc0b9": 0.804964557414518,
    }
    assert all(
        math.isclose(
            as_float(row["retained_l1_fraction"]),
            expected_retention[row["program_uid"]],
            rel_tol=0,
            abs_tol=1e-12,
        )
        for row in focal_rows
    )
    assert any(int(row["n_excluded_signature_genes"]) > 0 for row in rows)
    return "117 programs explicit in five cohorts; focal programs remain testable"


def bh_with_complete_family(values: list[float], family_size: int) -> list[float]:
    finite = [(index, value) for index, value in enumerate(values) if math.isfinite(value)]
    result = [math.nan] * len(values)
    ordered = sorted(finite, key=lambda pair: pair[1], reverse=True)
    running = 1.0
    m = len(ordered)
    for rank_from_largest, (index, value) in enumerate(ordered):
        denominator = m - rank_from_largest
        adjusted = min(1.0, family_size * value / denominator)
        running = min(running, adjusted)
        result[index] = running
    return result


def check_program_family() -> str:
    rows = read_tsv(OUT / "programs/program_meta_analysis.tsv")
    family = PRESPEC["programs"]["family_size"]
    expected_axes = set(
        PRESPEC["axes"]["co_primary"]
        + PRESPEC["axes"]["secondary"]
        + PRESPEC["axes"]["display_only"]
    )
    assert len(rows) == family * len(expected_axes)
    assert len(unique(rows, ("axis_id", "program_uid"))) == len(rows)
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["axis_id"]].append(row)
    assert set(grouped) == expected_axes
    for axis, block in grouped.items():
        observed_p = [as_float(row["p_value"]) for row in block]
        expected_q = bh_with_complete_family(observed_p, family)
        observed_q = [as_float(row["q_value"]) for row in block]
        for expected, observed in zip(expected_q, observed_q):
            if math.isnan(expected):
                assert math.isnan(observed)
            else:
                assert math.isclose(expected, observed, rel_tol=1e-10, abs_tol=1e-12), axis
    return "all five scorer families contain 117 rows; BH recomputed with n=117"


def check_windows_and_permutations() -> str:
    windows = read_tsv(OUT / "programs/fixed_program_windows.tsv.gz")
    expected_centers = {f"{value:.1f}" for value in PRESPEC["windows"]["centers"]}
    focal = set(PRESPEC["programs"]["focal_programs"].values())
    expected_axes = set(PRESPEC["axes"]["co_primary"] + PRESPEC["axes"]["display_only"])
    independent = set(PRESPEC["cohorts"]["independent_primary"])
    for axis in expected_axes:
        for dataset in independent:
            for program in focal:
                centers = {
                    f"{as_float(row['window_center']):.1f}"
                    for row in windows
                    if row["axis_id"] == axis
                    and row["dataset"] == dataset
                    and row["program_uid"] == program
                }
                assert centers == expected_centers, (axis, dataset, program, centers)
    assert all(as_bool(row["lower_inclusive"]) for row in windows)
    for row in windows:
        expected_upper = math.isclose(as_float(row["window_center"]), 0.9)
        assert as_bool(row["upper_inclusive"]) == expected_upper
    permutations = read_tsv(OUT / "programs/focal_within_stage_permutations.tsv")
    assert len(permutations) == 2 * 2 * 2
    assert all(
        int(row["permutation_replicates"])
        == PRESPEC["resampling"]["permutation_replicates"]
        for row in permutations
    )
    assert all(row["permutation_strata"] == "fibrosis_stage_by_sex" for row in permutations)
    return "nine fixed windows; eight position-preserving focal permutations"


def check_paired() -> str:
    scores = read_tsv(OUT / "paired/participant_scores.tsv")
    deltas = read_tsv(OUT / "paired/donor_deltas.tsv")
    tests = read_tsv(OUT / "paired/paired_tests.tsv")
    offset = read_tsv(OUT / "paired/stable_visit_offset_sensitivity.tsv")
    assert len(scores) == 54 * 2 * 2
    assert len(deltas) == 54 * 2
    assert len(tests) == 4
    assert len(offset) == 2 and all(as_bool(row["algebraically_rank_invariant"]) for row in offset)
    fibrosis = [row for row in tests if row["endpoint"] == "delta_fibrosis"]
    assert len(fibrosis) == 2
    assert all(row["alternative"] == "greater" for row in fibrosis)
    return "54 paired donors retained, including stable participants; fibrosis and NAS reported"


def check_gate_outputs() -> str:
    decisions = read_tsv(OUT / "decision/gate_decisions.tsv")
    promotion = read_tsv(OUT / "decision/promotion_decision.tsv")
    assert len(promotion) == 1
    row = promotion[0]
    overall = as_bool(row["overall_pass"])
    assert as_bool(row["explicit_approval_required"])
    expected_route = "figure3_proposal" if overall else "figure_s3_only"
    assert row["figure_route"] == expected_route
    required_gate_fragments = {
        "reproduction",
        "coverage",
        "agreement",
        "anchor",
        "program",
        "paired",
    }
    gate_ids = " ".join(row["gate_id"].lower() for row in decisions)
    assert all(fragment in gate_ids for fragment in required_gate_fragments)
    required = [row for row in decisions if as_bool(row["required_for_overall"])]
    assert required
    assert overall == all(as_bool(row["pass"]) for row in required)
    return f"overall_pass={overall}; route={expected_route}; no automatic promotion"


def check_figures() -> str:
    manifest = read_tsv(OUT / "figures/figure_manifest.tsv")
    assert len(manifest) == 4
    for row in manifest:
        path = pathlib.Path(row["path"])
        assert path.is_file() and OUT in path.resolve().parents
        with path.open("rb") as handle:
            assert handle.read(4) == b"%PDF"
        assert row["sha256"] == sha256(path)
        assert not as_bool(row["automatic_promotion"])
        assert not as_bool(row["useDingbats"])
    return "four candidate-only Cairo PDFs; all scorer and all-program controls present"


def check_protected_scopes() -> str:
    snapshot = read_tsv(OUT / "inputs/protected_scope_snapshot.tsv")
    assert snapshot
    for row in snapshot:
        path = pathlib.Path(row["protected_root"]) / row["relative_path"]
        assert path.is_file(), f"protected file disappeared: {path}"
        assert int(row["bytes"]) == path.stat().st_size, f"protected file size changed: {path}"
        assert row["sha256"] == sha256(path), f"protected file changed: {path}"
    scopes = {row["protected_scope"] for row in snapshot}
    assert scopes == {"current_figure3", "prior_continuum_benchmark"}
    return f"{len(snapshot)} protected Figure 3/benchmark files unchanged"


def check_input_hashes() -> str:
    rows = read_tsv(OUT / "inputs/input_manifest.tsv")
    checked = 0
    for row in rows:
        path = pathlib.Path(row["path"])
        if path.is_file():
            assert row["sha256"] == sha256(path), f"input changed: {path}"
            checked += 1
    assert checked >= 9
    return f"recomputed {checked} immutable input hashes"


def check_session_and_tests() -> str:
    sessions = [
        OUT / "inputs/sessionInfo.txt",
        OUT / "reproduction/sessionInfo.txt",
        OUT / "unsupervised/sessionInfo.txt",
        OUT / "projection/sessionInfo.txt",
        OUT / "programs/sessionInfo.txt",
        OUT / "paired/sessionInfo.txt",
        OUT / "decision/sessionInfo.txt",
        OUT / "figures/sessionInfo.txt",
        OUT / "report/sessionInfo.txt",
    ]
    assert all(path.is_file() and path.stat().st_size > 0 for path in sessions)
    test_log = OUT / "tests/contract_tests.log"
    text = test_log.read_text()
    assert "ALL_CONTRACT_TESTS_PASS" in text
    return "nine sessionInfo files and complete compute-node contract-test log"


def check_report() -> str:
    report = OUT / "report/candidate_report.md"
    text = report.read_text()
    required_phrases = [
        "candidate",
        "Slingshot",
        "alternative sensitivity",
        "cross-sectional",
        "explicit approval",
    ]
    assert all(phrase.lower() in text.lower() for phrase in required_phrases)
    assert PRESPEC["scientific_claim"] in text
    return "report preserves method roles, claim language, and approval boundary"


for check_id, function in [
    ("released_pc1_reproduction", check_reproduction),
    ("gencode_v49_signature_mapping", check_signature_mapping),
    ("participant_score_registry", check_score_registry),
    ("signature_coverage", check_coverage),
    ("rf_training_isolation", check_rf_isolation),
    ("program_exclusion_testability", check_program_exclusion_and_testability),
    ("bh_family_117", check_program_family),
    ("fixed_windows_and_permutations", check_windows_and_permutations),
    ("paired_54_donors", check_paired),
    ("promotion_gate_decision", check_gate_outputs),
    ("candidate_figures", check_figures),
    ("protected_scopes_unchanged", check_protected_scopes),
    ("input_hashes", check_input_hashes),
    ("session_info_and_contract_tests", check_session_and_tests),
    ("candidate_report_language", check_report),
]:
    check(check_id, function)


# Freeze every candidate artifact except the vendor checkout's internal .git
# database and this validation directory. Inputs from the released checkout are
# already commit-pinned and separately hashed by the reproduction stage.
manifest_rows: list[dict[str, str]] = []
for path in sorted(OUT.rglob("*")):
    if not path.is_file() or VALIDATION in path.parents:
        continue
    relative = path.relative_to(OUT)
    if relative.parts[:3] == ("external", "kamzolas_v1", ".git"):
        continue
    manifest_rows.append(
        {
            "relative_path": str(relative),
            "bytes": str(path.stat().st_size),
            "sha256": sha256(path),
        }
    )

with (VALIDATION / "candidate_file_manifest.tsv").open("x", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["relative_path", "bytes", "sha256"],
                            delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest_rows)

with (VALIDATION / "validation_checks.tsv").open("x", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["check_id", "status", "detail"],
                            delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(CHECKS)

passed = sum(row["status"] == "PASS" for row in CHECKS)
failed = len(CHECKS) - passed
summary = {
    "validator": "09_validate_candidate.py",
    "candidate_root": str(OUT),
    "checks_passed": passed,
    "checks_failed": failed,
    "candidate_files_hashed": len(manifest_rows),
    "status": "PASS" if failed == 0 else "FAIL",
}
with (VALIDATION / "validation_summary.json").open("x") as handle:
    json.dump(summary, handle, indent=2, sort_keys=True)
    handle.write("\n")

print(f"CANDIDATE_VALIDATION_{summary['status']}: {passed}/{len(CHECKS)} checks passed")
if failed:
    for row in CHECKS:
        if row["status"] == "FAIL":
            print(f"FAIL\t{row['check_id']}\t{row['detail']}", file=sys.stderr)
    raise SystemExit(1)
