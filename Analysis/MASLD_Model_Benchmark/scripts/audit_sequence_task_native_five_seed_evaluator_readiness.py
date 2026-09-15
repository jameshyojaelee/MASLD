#!/usr/bin/env python3
"""Outcome-blind readiness check for the full task-native five-seed rectangle."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import verify_frozen_tree


SCHEMA = "masld-bench-sequence-task-native-five-seed-evaluator-contract-v1"
FIT_FAMILIES = (
    "bpnet",
    "chrombpnet",
    "sequence_cnn_control",
    "sequence_transformer_control",
)
CANDIDATE_VIEWS = (
    "bpnet",
    "chrombpnet_full",
    "chrombpnet_nobias",
    "sequence_cnn_control",
    "sequence_transformer_control",
)
FIT_DISPOSITIONS = {
    "complete_compatible",
    "not_dispatched_or_pending",
    "terminal_fit_failure",
    "frozen_incompatible",
}
PASS_RECEIPT_STATES = {"passed", "recovered_passed"}
FAIL_RECEIPT_STATES = {
    "failed",
    "existing_failed",
    "failed_invalid_artifact",
    "failed_missing_artifact",
}
FIT_FIELDS = (
    "model_id",
    "lineage_id",
    "outer_fold",
    "split_id",
    "seed",
    "logical_work_id",
    "disposition",
    "canonical_artifact_path",
    "canonical_artifacts_sha256",
    "receipt_count",
    "receipt_statuses",
)
VIEW_FIELDS = (
    "candidate_view",
    "root_model_id",
    "lineage_id",
    "outer_fold",
    "split_id",
    "seed",
    "disposition",
    "canonical_artifact_path",
    "canonical_artifacts_sha256",
)


class ReadinessError(RuntimeError):
    """Raised when the frozen rectangle or evaluator separation differs."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--correction-root", type=Path, required=True)
    parser.add_argument("--queue-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
    path: Path,
    rows: Iterable[Mapping[str, object]],
    fields: Sequence[str],
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def load_contract(path: Path) -> dict[str, Any]:
    contract = json.loads(path.read_text(encoding="utf-8"))
    if contract.get("schema_version") != SCHEMA:
        raise ReadinessError("evaluator contract schema differs")
    rectangle = contract.get("rectangle", {})
    if (
        rectangle.get("expected_fits") != 500
        or rectangle.get("expected_prediction_views") != 625
        or rectangle.get("expected_units_per_candidate_view") != 125
        or rectangle.get("diagonal_outer_folds") != list(range(5))
        or rectangle.get("fixed_seeds")
        != [20260824, 20260825, 20260826, 20260827, 20260828]
        or len(rectangle.get("lineages", ())) != 5
        or tuple(rectangle.get("fit_families", ())) != FIT_FAMILIES
        or tuple(rectangle.get("candidate_views", ())) != CANDIDATE_VIEWS
    ):
        raise ReadinessError("evaluator rectangle cardinality differs")
    gate = contract.get("readiness_gate", {})
    if (
        set(gate.get("allowed_fit_dispositions", ())) != FIT_DISPOSITIONS
        or gate.get("required_complete_fits") != 500
        or gate.get("required_complete_prediction_views") != 625
        or gate.get("model_specific_partial_ranking_allowed") is not False
        or gate.get("outcome_access_before_full_lock_allowed") is not False
    ):
        raise ReadinessError("readiness or partial-ranking contract differs")
    firewall = contract.get("firewall", {})
    if (
        firewall.get("prediction_or_metric_values_may_be_read_by_readiness_audit") is not False
        or firewall.get("observed_outcome_files_may_be_read_by_readiness_audit") is not False
        or firewall.get("test_role_authorized") is not False
        or firewall.get("sealed_data_authorized") is not False
    ):
        raise ReadinessError("outcome firewall differs")
    uncertainty = contract.get("uncertainty", {})
    if (
        uncertainty.get("bootstrap_replicates") != 10_000
        or uncertainty.get("confidence_level") != 0.95
        or uncertainty.get("method") != "paired_two_way_donor_genomic_block_bootstrap"
        or uncertainty.get("cells_windows_and_seeds_are_independent_units") is not False
        or uncertainty.get("cross_fold_zero_cells_may_be_imputed") is not False
    ):
        raise ReadinessError("uncertainty contract differs")
    return contract


def load_admission_module(root: Path, contract: Mapping[str, Any]) -> Any:
    source = root / "scripts/build_sequence_task_native_five_seed_rectangle.py"
    expected = contract["frozen_sources"]["artifact_admission_source_sha256"]
    if digest(source) != expected:
        raise ReadinessError("frozen artifact-admission source differs")
    spec = importlib.util.spec_from_file_location("five_seed_rectangle_admission", source)
    if spec is None or spec.loader is None:
        raise ReadinessError("artifact-admission source cannot be loaded")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_frozen_sources(
    *,
    root: Path,
    contract: Mapping[str, Any],
    plan_root: Path,
    correction_root: Path,
    queue_root: Path,
) -> None:
    sources = contract["frozen_sources"]
    bindings = (
        (plan_root, sources["rectangle_plan_artifacts_sha256"]),
        (correction_root, sources["execution_correction_artifacts_sha256"]),
        (queue_root, sources["corrected_queue_artifacts_sha256"]),
    )
    for path, expected in bindings:
        if root not in path.parents:
            raise ReadinessError("frozen source escapes benchmark root")
        if digest(path / "ARTIFACTS.json") != expected:
            raise ReadinessError(f"frozen source ARTIFACTS differs: {path.name}")
        verify_frozen_tree(path)
    if digest(root / "config/evaluation/inference.toml") != sources["inference_contract_sha256"]:
        raise ReadinessError("project inference contract differs")
    if digest(root / "src/masld_bench/evaluators/stats.py") != sources["two_way_bootstrap_source_sha256"]:
        raise ReadinessError("two-way bootstrap implementation differs")

    plan_summary = json.loads((plan_root / "plan/summary.json").read_text(encoding="utf-8"))
    correction = json.loads((correction_root / "ARTIFACTS.json").read_text(encoding="utf-8"))["metadata"]
    queue = json.loads((queue_root / "ARTIFACTS.json").read_text(encoding="utf-8"))["metadata"]
    if (
        plan_summary.get("expected_fits") != 500
        or plan_summary.get("expected_prediction_views") != 625
        or plan_summary.get("prediction_or_metric_values_read") is not False
        or plan_summary.get("cross_model_ranking_performed") is not False
        or correction.get("scientific_task_contract_changed") is not False
        or correction.get("v1_queue_items_claimed") != 0
        or correction.get("v1_queue_items_submitted") != 0
        or queue.get("logical_fits") != 438
        or queue.get("submission_authority") != "dispatcher_only_no_manual_gpu_sbatch"
        or queue.get("qos") != "nslab"
        or queue.get("manual_gpu_sbatch_allowed") is not False
    ):
        raise ReadinessError("plan, correction, or dispatcher contract differs")


def fit_key(row: Mapping[str, Any]) -> tuple[str, str, int, int]:
    return (
        str(row["model_id"]),
        str(row["lineage_id"]),
        int(row["outer_fold"]),
        int(row["seed"]),
    )


def logical_work_id(key: tuple[str, str, int, int]) -> str:
    model, lineage, fold, seed = key
    return f"rect__{model}__{lineage}__d{fold}g{fold}__s{seed}"


def expected_keys(contract: Mapping[str, Any]) -> list[tuple[str, str, int, int]]:
    rectangle = contract["rectangle"]
    result = [
        (model, lineage, fold, seed)
        for model in rectangle["fit_families"]
        for lineage in rectangle["lineages"]
        for fold in rectangle["diagonal_outer_folds"]
        for seed in rectangle["fixed_seeds"]
    ]
    if len(result) != 500 or len(set(result)) != 500:
        raise ReadinessError("expected fit rectangle is not exactly 500 unique units")
    return result


def receipt_index(root: Path, keys: Sequence[tuple[str, str, int, int]]) -> tuple[dict[tuple[str, str, int, int], list[dict[str, Any]]], list[dict[str, object]]]:
    by_logical = {logical_work_id(key): key for key in keys}
    indexed: dict[tuple[str, str, int, int], list[dict[str, Any]]] = defaultdict(list)
    rows: list[dict[str, object]] = []
    for path in sorted((root / "executions").glob("*/task_receipts/*.json")):
        try:
            receipt = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        logical_id = receipt.get("logical_work_id")
        if logical_id not in by_logical:
            continue
        key = by_logical[logical_id]
        if (
            receipt.get("model_id") != key[0]
            or receipt.get("lineage_id") != key[1]
            or receipt.get("outer_fold") != key[2]
            or receipt.get("seed") != key[3]
            or receipt.get("status") not in PASS_RECEIPT_STATES | FAIL_RECEIPT_STATES
        ):
            raise ReadinessError(f"rectangle task receipt identity or status differs: {path}")
        indexed[key].append(receipt)
        rows.append({
            "logical_work_id": logical_id,
            "model_id": key[0],
            "lineage_id": key[1],
            "outer_fold": key[2],
            "seed": key[3],
            "status": receipt["status"],
            "expected_artifact": receipt.get("expected_artifact", "not_available"),
            "artifacts_sha256": receipt.get("artifacts_sha256", "not_available"),
            "receipt_path": str(path.resolve()),
        })
    return indexed, rows


def relevant_rejections(
    root: Path,
    plan_contract: Mapping[str, Any],
    admission: Any,
    compatible_paths: set[str],
) -> dict[tuple[str, str, int, int], list[dict[str, str]]]:
    result: dict[tuple[str, str, int, int], list[dict[str, str]]] = defaultdict(list)
    for path in sorted((root / "executions").glob("*/ARTIFACTS.json")):
        artifact_root = str(path.parent.resolve())
        if artifact_root in compatible_paths:
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        metadata = payload.get("metadata", {})
        split = re.fullmatch(r"donor([0-4])_genomic([0-4])", str(metadata.get("split_id")))
        if (
            metadata.get("dataset_id") != plan_contract["dataset_id"]
            or metadata.get("model_id") not in plan_contract["fit_families"]
            or metadata.get("lineage_id") not in plan_contract["lineages"]
            or metadata.get("seed") not in plan_contract["fixed_seeds"]
            or split is None
            or split.group(1) != split.group(2)
        ):
            continue
        key = (
            metadata["model_id"],
            metadata["lineage_id"],
            int(split.group(1)),
            int(metadata["seed"]),
        )
        candidate, reason = admission.artifact_candidate(path, plan_contract)
        if candidate is None:
            result[key].append({
                "artifact_path": str(path.parent.resolve()),
                "reason": reason,
            })
    return result


def classify_rectangle(
    keys: Sequence[tuple[str, str, int, int]],
    compatible: Mapping[tuple[str, str, int, int], Mapping[str, Any]],
    receipts: Mapping[tuple[str, str, int, int], Sequence[Mapping[str, Any]]],
    incompatible: Mapping[tuple[str, str, int, int], Sequence[Mapping[str, str]]],
) -> tuple[list[dict[str, object]], str]:
    rows: list[dict[str, object]] = []
    for key in keys:
        fit_receipts = list(receipts.get(key, ()))
        statuses = sorted({str(row["status"]) for row in fit_receipts})
        artifact = compatible.get(key)
        if artifact is not None:
            disposition = "complete_compatible"
        elif incompatible.get(key) or any(status in PASS_RECEIPT_STATES for status in statuses):
            disposition = "frozen_incompatible"
        elif any(status in FAIL_RECEIPT_STATES for status in statuses):
            disposition = "terminal_fit_failure"
        else:
            disposition = "not_dispatched_or_pending"
        rows.append({
            "model_id": key[0],
            "lineage_id": key[1],
            "outer_fold": key[2],
            "split_id": f"donor{key[2]}_genomic{key[2]}",
            "seed": key[3],
            "logical_work_id": logical_work_id(key),
            "disposition": disposition,
            "canonical_artifact_path": artifact["artifact_path"] if artifact else "not_available",
            "canonical_artifacts_sha256": artifact["artifacts_sha256"] if artifact else "not_available",
            "receipt_count": len(fit_receipts),
            "receipt_statuses": ",".join(statuses) if statuses else "not_available",
        })
    counts = Counter(str(row["disposition"]) for row in rows)
    if counts["frozen_incompatible"]:
        status = "blocked_incompatible_artifact"
    elif counts["terminal_fit_failure"]:
        status = "blocked_terminal_failure_requires_new_immutable_retry"
    elif counts["not_dispatched_or_pending"]:
        status = "waiting_for_full_rectangle"
    elif counts["complete_compatible"] == len(keys) == 500:
        status = "ready_full_rectangle_locked"
    else:
        raise ReadinessError("fit dispositions do not resolve to a gate status")
    return rows, status


def expand_views(
    fit_rows: Sequence[Mapping[str, object]],
    prediction_views: Mapping[str, Sequence[str]],
) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for row in fit_rows:
        model = str(row["model_id"])
        for view in prediction_views[model]:
            result.append({
                "candidate_view": view,
                "root_model_id": model,
                "lineage_id": row["lineage_id"],
                "outer_fold": row["outer_fold"],
                "split_id": row["split_id"],
                "seed": row["seed"],
                "disposition": row["disposition"],
                "canonical_artifact_path": row["canonical_artifact_path"],
                "canonical_artifacts_sha256": row["canonical_artifacts_sha256"],
            })
    if len(result) != 625:
        raise ReadinessError("prediction-view rectangle is not exactly 625 units")
    return result


def readiness_summary(
    *,
    fit_rows: Sequence[Mapping[str, object]],
    view_rows: Sequence[Mapping[str, object]],
    gate_status: str,
) -> dict[str, Any]:
    fit_counts = Counter(str(row["disposition"]) for row in fit_rows)
    view_counts = Counter(str(row["disposition"]) for row in view_rows)
    complete_by_view = Counter(
        str(row["candidate_view"])
        for row in view_rows
        if row["disposition"] == "complete_compatible"
    )
    full = (
        gate_status == "ready_full_rectangle_locked"
        and fit_counts["complete_compatible"] == 500
        and view_counts["complete_compatible"] == 625
        and set(complete_by_view.values()) == {125}
        and len(complete_by_view) == 5
    )
    return {
        "schema_version": "masld-bench-sequence-task-native-five-seed-readiness-v1",
        "status": gate_status,
        "expected_fits": 500,
        "expected_prediction_views": 625,
        "fit_disposition_counts": dict(sorted(fit_counts.items())),
        "prediction_view_disposition_counts": dict(sorted(view_counts.items())),
        "complete_units_by_candidate_view": dict(sorted(complete_by_view.items())),
        "full_compatible_rectangle": full,
        "outcome_access_authorized": full,
        "development_cross_model_ranking_authorized": full,
        "partial_cross_model_ranking_authorized": False,
        "cross_model_ranking_performed": False,
        "prediction_or_metric_values_read": False,
        "observed_outcome_files_read": False,
        "test_role_authorized": False,
        "sealed_data_authorized": False,
        "seeds_are_biological_replicates": False,
        "terminal_failures_are_zero_filled_or_dropped": False,
    }


def main() -> None:
    args = parse_args()
    root = args.root.resolve(strict=True)
    contract_path = args.contract.resolve(strict=True)
    plan_root = args.plan_root.resolve(strict=True)
    correction_root = args.correction_root.resolve(strict=True)
    queue_root = args.queue_root.resolve(strict=True)
    output = args.output.resolve(strict=False)
    if output.exists():
        raise ReadinessError("output already exists")
    for path in (contract_path, plan_root, correction_root, queue_root, output):
        if root not in path.parents:
            raise ReadinessError("input or output escapes benchmark root")

    contract = load_contract(contract_path)
    validate_frozen_sources(
        root=root,
        contract=contract,
        plan_root=plan_root,
        correction_root=correction_root,
        queue_root=queue_root,
    )
    admission = load_admission_module(root, contract)
    plan_contract = admission.load_contract(plan_root / "plan/contract_snapshot.json")
    if (
        plan_contract["lineages"] != contract["rectangle"]["lineages"]
        or plan_contract["diagonal_outer_folds"] != contract["rectangle"]["diagonal_outer_folds"]
        or plan_contract["fixed_seeds"] != contract["rectangle"]["fixed_seeds"]
        or plan_contract["fit_families"] != contract["rectangle"]["fit_families"]
    ):
        raise ReadinessError("evaluator rectangle differs from frozen fit rectangle")

    keys = expected_keys(contract)
    plan_rows = read_tsv(plan_root / "plan/expected_fit_matrix.tsv")
    if len(plan_rows) != 500 or {fit_key(row) for row in plan_rows} != set(keys):
        raise ReadinessError("frozen expected-fit matrix differs")
    canonical, duplicates, _rejected = admission.scan_existing(root, plan_contract)
    compatible = {fit_key(row): row for row in canonical}
    receipts, receipt_rows = receipt_index(root, keys)
    compatible_paths = {
        str(row["artifact_path"]) for row in canonical
    } | {
        str(row["artifact_path"]) for row in duplicates
    }
    incompatible = relevant_rejections(
        root,
        plan_contract,
        admission,
        compatible_paths,
    )
    fit_rows, gate_status = classify_rectangle(keys, compatible, receipts, incompatible)
    view_rows = expand_views(fit_rows, plan_contract["prediction_views"])
    summary = readiness_summary(
        fit_rows=fit_rows,
        view_rows=view_rows,
        gate_status=gate_status,
    )

    output.mkdir(mode=0o750)
    (output / "contract_snapshot.json").write_bytes(contract_path.read_bytes())
    write_tsv(output / "fit_dispositions.tsv", fit_rows, FIT_FIELDS)
    write_tsv(output / "prediction_view_dispositions.tsv", view_rows, VIEW_FIELDS)
    write_tsv(output / "terminal_receipts.tsv", receipt_rows, (
        "logical_work_id", "model_id", "lineage_id", "outer_fold", "seed", "status",
        "expected_artifact", "artifacts_sha256", "receipt_path",
    ))
    write_tsv(output / "duplicate_compatible_aliases.tsv", duplicates, (
        "model_id", "lineage_id", "outer_fold", "split_id", "seed", "artifact_path",
        "artifacts_sha256", "terminal_slurm_job_id", "canonical_artifact_path", "canonical_rule",
    ))
    rejected_rows = [
        {
            "model_id": key[0],
            "lineage_id": key[1],
            "outer_fold": key[2],
            "seed": key[3],
            **row,
        }
        for key in sorted(incompatible)
        for row in incompatible[key]
    ]
    write_tsv(output / "incompatible_relevant_artifacts.tsv", rejected_rows, (
        "model_id", "lineage_id", "outer_fold", "seed", "artifact_path", "reason",
    ))
    summary.update({
        "contract_sha256": digest(contract_path),
        "fit_dispositions_sha256": digest(output / "fit_dispositions.tsv"),
        "prediction_view_dispositions_sha256": digest(output / "prediction_view_dispositions.tsv"),
        "compatible_duplicate_aliases": len(duplicates),
        "relevant_incompatible_artifacts": len(rejected_rows),
        "terminal_receipts": len(receipt_rows),
        "readiness_audit_scope": "artifact_manifests_frozen_tree_hashes_and_task_receipts_only",
    })
    (output / "readiness_lock.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
