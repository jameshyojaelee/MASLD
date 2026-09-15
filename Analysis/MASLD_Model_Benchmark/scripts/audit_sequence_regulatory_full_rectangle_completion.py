#!/usr/bin/env python3
"""Outcome-blind completion census for sequence-regulatory rectangles."""

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


SCHEMA = "masld-bench-sequence-regulatory-full-rectangle-completion-v1"
SEEDS = (20260824, 20260825, 20260826, 20260827, 20260828)
FOLDS = tuple(range(5))
SCBASSET_FIELDS = (
    "outer_fold",
    "split_id",
    "seed",
    "disposition",
    "model_root",
    "model_artifacts_sha256",
    "prediction_root",
    "prediction_artifacts_sha256",
)
SEQUENCE_SOURCE_FIELDS = (
    "artifact_role",
    "lineage_id",
    "split_id",
    "root",
    "artifacts_sha256",
)


class CompletionError(RuntimeError):
    """Raised when a frozen source, rectangle, or separation differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _write_tsv(
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


def _resolve(root: Path, value: str) -> Path:
    path = (root / value).resolve(strict=True)
    if root not in path.parents:
        raise CompletionError(f"bound path escapes benchmark root: {value}")
    return path


def _verify_bound_root(root: Path, value: str, expected: str) -> Path:
    path = _resolve(root, value)
    if digest(path / "ARTIFACTS.json") != expected:
        raise CompletionError(f"bound ARTIFACTS differs: {value}")
    verify_frozen_tree(path)
    return path


def load_contract(path: Path, root: Path) -> dict[str, Any]:
    contract = json.loads(path.read_text(encoding="utf-8"))
    if contract.get("schema_version") != SCHEMA:
        raise CompletionError("completion contract schema differs")
    task = contract.get("task_families", {}).get("task_native_profile", {})
    scbasset = contract.get("task_families", {}).get("scbasset_sequence_only", {})
    gate = contract.get("gate", {})
    claims = contract.get("claims", {})
    if (
        task.get("expected_fits") != 500
        or task.get("expected_prediction_views") != 625
        or task.get("required_units_per_view") != 125
        or scbasset.get("expected_fits") != 25
        or scbasset.get("expected_prediction_views") != 25
        or scbasset.get("outer_folds") != list(FOLDS)
        or scbasset.get("fixed_seeds") != list(SEEDS)
        or gate.get("require_both_full_rectangles") is not True
        or gate.get("partial_ranking_allowed") is not False
        or gate.get("zero_fill_allowed") is not False
        or gate.get("drop_failed_unit_allowed") is not False
        or gate.get("prediction_values_may_be_read") is not False
        or gate.get("metric_values_may_be_read") is not False
        or gate.get("observed_outcomes_may_be_read") is not False
        or gate.get("sealed_data_may_be_read") is not False
        or claims.get("cross_task_family_ranking_allowed") is not False
        or claims.get("universal_claim_allowed") is not False
    ):
        raise CompletionError("completion cardinality, firewall, or claim contract differs")
    for section, key, expected_schema in (
        (
            task,
            "evaluator_contract",
            "masld-bench-sequence-task-native-five-seed-evaluator-contract-v1",
        ),
        (scbasset, "contract", "masld-bench-scbasset-five-seed-rectangle-v1"),
        (
            scbasset,
            "evaluator_contract",
            "masld-bench-scbasset-five-seed-evaluator-contract-v1",
        ),
    ):
        bound = _resolve(root, section[key])
        if digest(bound) != section[f"{key}_sha256"]:
            raise CompletionError(f"bound contract differs: {bound}")
        payload = json.loads(bound.read_text(encoding="utf-8"))
        if payload.get("schema_version") != expected_schema:
            raise CompletionError(f"bound contract schema differs: {bound}")
    return contract


def _load_readiness_module(root: Path) -> Any:
    source = root / "scripts/audit_sequence_task_native_five_seed_evaluator_readiness.py"
    spec = importlib.util.spec_from_file_location("sequence_five_seed_readiness", source)
    if spec is None or spec.loader is None:
        raise CompletionError("sequence readiness source cannot be loaded")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sequence_census(
    *, root: Path, contract: Mapping[str, Any]
) -> tuple[dict[str, Any], list[dict[str, object]], list[dict[str, object]]]:
    task = contract["task_families"]["task_native_profile"]
    readiness = _load_readiness_module(root)
    evaluator_contract_path = _resolve(root, task["evaluator_contract"])
    evaluator_contract = readiness.load_contract(evaluator_contract_path)
    plan_root = _verify_bound_root(root, task["plan_root"], task["plan_artifacts_sha256"])
    correction_root = _verify_bound_root(
        root, task["correction_root"], task["correction_artifacts_sha256"]
    )
    queue_root = _verify_bound_root(root, task["queue_root"], task["queue_artifacts_sha256"])
    readiness.validate_frozen_sources(
        root=root,
        contract=evaluator_contract,
        plan_root=plan_root,
        correction_root=correction_root,
        queue_root=queue_root,
    )
    admission = readiness.load_admission_module(root, evaluator_contract)
    plan_contract = admission.load_contract(plan_root / "plan/contract_snapshot.json")
    keys = readiness.expected_keys(evaluator_contract)
    canonical, duplicates, _rejected = admission.scan_existing(root, plan_contract)
    compatible = {readiness.fit_key(row): row for row in canonical}
    receipts, _receipt_rows = readiness.receipt_index(root, keys)
    compatible_paths = {str(row["artifact_path"]) for row in canonical} | {
        str(row["artifact_path"]) for row in duplicates
    }
    incompatible = readiness.relevant_rejections(
        root, plan_contract, admission, compatible_paths
    )
    fit_rows, status = readiness.classify_rectangle(
        keys, compatible, receipts, incompatible
    )
    view_rows = readiness.expand_views(fit_rows, plan_contract["prediction_views"])
    summary = readiness.readiness_summary(
        fit_rows=fit_rows,
        view_rows=view_rows,
        gate_status=status,
    )
    if summary["cross_model_ranking_performed"] is not False:
        raise CompletionError("sequence readiness performed ranking")
    return summary, fit_rows, view_rows


def _terminal_job_id(path: Path) -> int:
    match = re.search(r"-(\d+)$", path.name)
    return int(match.group(1)) if match else 2**63 - 1


def _scbasset_artifacts(root: Path) -> tuple[
    dict[tuple[int, int], list[dict[str, Any]]],
    dict[tuple[int, int], list[dict[str, Any]]],
]:
    models: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    predictions: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for manifest_path in sorted((root / "executions").glob("*/ARTIFACTS.json")):
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        metadata = payload.get("metadata", {})
        artifact_class = metadata.get("artifact_class")
        if artifact_class not in {"scbasset_fold_model", "scbasset_fold_predictions"}:
            continue
        split = re.fullmatch(r"donor([0-4])_genomic([0-4])", str(metadata.get("split_id")))
        seed = metadata.get("seed")
        if (
            metadata.get("dataset_id") != "gse296875"
            or metadata.get("model_id") != "scbasset"
            or metadata.get("status") != "passed"
            or split is None
            or split.group(1) != split.group(2)
            or seed not in SEEDS
        ):
            continue
        fold = int(split.group(1))
        record = {
            "root": manifest_path.parent.resolve(),
            "artifacts_sha256": digest(manifest_path),
            "metadata": metadata,
            "terminal_job_id": _terminal_job_id(manifest_path.parent),
        }
        target = models if artifact_class == "scbasset_fold_model" else predictions
        target[(fold, int(seed))].append(record)
    return models, predictions


def _scbasset_terminal_failures(root: Path) -> set[tuple[int, int]]:
    failures: set[tuple[int, int]] = set()
    for path in sorted((root / "executions").glob("scbasset-*/fit_receipts/*.json")):
        try:
            row = json.loads(path.read_text(encoding="utf-8"))
            key = (int(row["outer_fold"]), int(row["seed"]))
        except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
            continue
        if key[0] not in FOLDS or key[1] not in SEEDS:
            continue
        train_rc = row.get("train_returncode")
        predict_rc = row.get("predict_returncode")
        if train_rc != 0 or predict_rc != 0:
            failures.add(key)
    return failures


def scbasset_census(
    *, root: Path, contract: Mapping[str, Any]
) -> tuple[dict[str, Any], list[dict[str, object]]]:
    section = contract["task_families"]["scbasset_sequence_only"]
    _verify_bound_root(
        root, section["admission_root"], section["admission_artifacts_sha256"]
    )
    models, predictions = _scbasset_artifacts(root)
    terminal_failures = _scbasset_terminal_failures(root)
    rows: list[dict[str, object]] = []
    for fold in FOLDS:
        split_id = f"donor{fold}_genomic{fold}"
        for seed in SEEDS:
            key = (fold, seed)
            fit_candidates = sorted(
                models.get(key, ()),
                key=lambda row: (row["terminal_job_id"], str(row["root"])),
            )
            prediction_candidates = sorted(
                predictions.get(key, ()),
                key=lambda row: (row["terminal_job_id"], str(row["root"])),
            )
            pairs = []
            for candidate in prediction_candidates:
                metadata = candidate["metadata"]
                for fit_candidate in fit_candidates:
                    if (
                        metadata.get("model_artifacts_sha256")
                        == fit_candidate["artifacts_sha256"]
                        and metadata.get("input_artifacts_sha256")
                        == fit_candidate["metadata"].get("input_artifacts_sha256")
                        and metadata.get("evaluation_role") == "valid"
                        and metadata.get("held_donor_atac_used") is False
                        and metadata.get("held_cell_embedding_available") is False
                    ):
                        pairs.append((fit_candidate, candidate))
            pairs.sort(
                key=lambda pair: (
                    pair[1]["terminal_job_id"],
                    pair[0]["terminal_job_id"],
                    str(pair[1]["root"]),
                    str(pair[0]["root"]),
                )
            )
            model, prediction = pairs[0] if pairs else (None, None)
            disposition = (
                "complete_compatible"
                if model is not None and prediction is not None
                else (
                    "terminal_fit_failure"
                    if key in terminal_failures
                    else "not_dispatched_or_pending"
                )
            )
            if disposition == "complete_compatible":
                verify_frozen_tree(model["root"])
                verify_frozen_tree(prediction["root"])
            rows.append(
                {
                    "outer_fold": fold,
                    "split_id": split_id,
                    "seed": seed,
                    "disposition": disposition,
                    "model_root": str(model["root"]) if model else "not_available",
                    "model_artifacts_sha256": (
                        model["artifacts_sha256"] if model else "not_available"
                    ),
                    "prediction_root": (
                        str(prediction["root"]) if prediction else "not_available"
                    ),
                    "prediction_artifacts_sha256": (
                        prediction["artifacts_sha256"]
                        if prediction
                        else "not_available"
                    ),
                }
            )
    counts = Counter(str(row["disposition"]) for row in rows)
    full = len(rows) == 25 and counts.get("complete_compatible", 0) == 25
    status = (
        "ready_full_rectangle_locked"
        if full
        else (
            "blocked_terminal_failure_requires_new_immutable_retry"
            if counts.get("terminal_fit_failure", 0)
            else "waiting_for_full_rectangle"
        )
    )
    return {
        "schema_version": "masld-bench-scbasset-five-seed-completion-v1",
        "status": status,
        "expected_fits": 25,
        "expected_prediction_views": 25,
        "fit_disposition_counts": dict(sorted(counts.items())),
        "full_compatible_rectangle": full,
        "outcome_access_authorized": full,
        "partial_ranking_authorized": False,
        "prediction_values_read": False,
        "metric_values_read": False,
        "observed_outcomes_read": False,
        "sealed_data_read": False,
        "rna_conditioned_claim_allowed": False,
    }, rows


def combined_status(sequence_status: str, scbasset_status: str) -> str:
    if (
        sequence_status == "ready_full_rectangle_locked"
        and scbasset_status == "ready_full_rectangle_locked"
    ):
        return "ready_both_rectangles_locked"
    if sequence_status.startswith("blocked_") or scbasset_status.startswith("blocked_"):
        return "blocked_rectangle_requires_new_immutable_retry"
    return "waiting_for_full_rectangles"


def verify_separate_family(
    *, root: Path, contract: Mapping[str, Any]
) -> list[dict[str, object]]:
    rows = []
    for binding in contract["separate_gse281364_development_family"]:
        artifact = _verify_bound_root(
            root, binding["root"], binding["artifacts_sha256"]
        )
        metadata = json.loads((artifact / "ARTIFACTS.json").read_text(encoding="utf-8"))[
            "metadata"
        ]
        rows.append(
            {
                "family": binding["family"],
                "root": str(artifact),
                "artifacts_sha256": binding["artifacts_sha256"],
                "status": metadata.get("status", "not_available"),
                "task_family": "gse281364_variant_regulatory_development",
                "comparable_to_gse296875_rectangles": False,
            }
        )
    return rows


def write_sequence_source_manifest(
    *,
    root: Path,
    contract: Mapping[str, Any],
    fit_rows: Sequence[Mapping[str, object]],
    output: Path,
) -> None:
    rows: list[dict[str, object]] = []
    for binding in contract["evaluator_inputs"]:
        _verify_bound_root(root, binding["root"], binding["artifacts_sha256"])
        rows.append(
            {
                "artifact_role": "evaluator_input_bundle",
                "lineage_id": binding["lineage_id"],
                "split_id": "all",
                "root": str(_resolve(root, binding["root"])),
                "artifacts_sha256": binding["artifacts_sha256"],
            }
        )
    for row in fit_rows:
        if row["disposition"] != "complete_compatible":
            continue
        rows.append(
            {
                "artifact_role": "prediction_artifact",
                "lineage_id": row["lineage_id"],
                "split_id": row["split_id"],
                "root": row["canonical_artifact_path"],
                "artifacts_sha256": row["canonical_artifacts_sha256"],
            }
        )
    _write_tsv(output, rows, SEQUENCE_SOURCE_FIELDS)


def run(*, root: Path, contract_path: Path, output: Path) -> dict[str, Any]:
    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    output = output.resolve(strict=False)
    if output.exists() or root not in contract_path.parents or root not in output.parents:
        raise CompletionError("completion input/output path differs")
    contract = load_contract(contract_path, root)
    sequence_summary, fit_rows, view_rows = sequence_census(root=root, contract=contract)
    scbasset_summary, scbasset_rows = scbasset_census(root=root, contract=contract)
    separate_family = verify_separate_family(root=root, contract=contract)
    status = combined_status(sequence_summary["status"], scbasset_summary["status"])
    both_full = status == "ready_both_rectangles_locked"
    output.mkdir(mode=0o750)
    _write_tsv(output / "sequence_fit_dispositions.tsv", fit_rows, (
        "model_id", "lineage_id", "outer_fold", "split_id", "seed",
        "logical_work_id", "disposition", "canonical_artifact_path",
        "canonical_artifacts_sha256", "receipt_count", "receipt_statuses",
    ))
    _write_tsv(output / "sequence_prediction_view_dispositions.tsv", view_rows, (
        "candidate_view", "root_model_id", "lineage_id", "outer_fold",
        "split_id", "seed", "disposition", "canonical_artifact_path",
        "canonical_artifacts_sha256",
    ))
    _write_tsv(output / "scbasset_fit_dispositions.tsv", scbasset_rows, SCBASSET_FIELDS)
    _write_tsv(output / "separate_gse281364_family.tsv", separate_family, (
        "family", "root", "artifacts_sha256", "status", "task_family",
        "comparable_to_gse296875_rectangles",
    ))
    write_sequence_source_manifest(
        root=root,
        contract=contract,
        fit_rows=fit_rows,
        output=output / "sequence_evaluation_sources.tsv",
    )
    summary = {
        "schema_version": "masld-bench-sequence-regulatory-full-rectangle-completion-receipt-v1",
        "status": status,
        "contract_sha256": digest(contract_path),
        "task_native_sequence": sequence_summary,
        "scbasset": scbasset_summary,
        "separate_gse281364_development_families": len(separate_family),
        "both_full_rectangles": both_full,
        "evaluation_launch_authorized": both_full,
        "partial_ranking_authorized": False,
        "cross_task_family_ranking_authorized": False,
        "prediction_values_read": False,
        "metric_values_read": False,
        "observed_outcomes_read": False,
        "sealed_data_read": False,
        "champion_claim_allowed": False,
        "universal_claim_allowed": False,
    }
    (output / "completion_lock.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = run(root=args.root, contract_path=args.contract, output=args.output)
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
