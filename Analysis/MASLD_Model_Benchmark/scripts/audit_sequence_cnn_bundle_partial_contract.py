#!/usr/bin/env python3
"""Audit the exact eight-unit sequence-CNN partial-development coverage."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
from typing import Any

from masld_bench.artifacts import verify_frozen_tree


SCHEMA = "masld-bench-sequence-cnn-bundle-partial-evaluation-v1"
SOURCE_FIELDS = (
    "artifact_role",
    "lineage_id",
    "split_id",
    "root",
    "artifacts_sha256",
)


class SequenceCNNContractError(RuntimeError):
    """Raised when bundle coverage or safety differs from its selection record."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SequenceCNNContractError(f"JSON object differs: {path}")
    return value


def resolve_project_path(root: Path, value: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise SequenceCNNContractError("unsafe project-relative path")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    return path


def validate_expected_units(units: object) -> list[dict[str, Any]]:
    if not isinstance(units, list) or len(units) != 8:
        raise SequenceCNNContractError("expected unit count differs")
    required = {"logical_work_id", "lineage_id", "outer_fold", "split_id", "seed"}
    result: list[dict[str, Any]] = []
    seen: set[str] = set()
    for raw in units:
        if not isinstance(raw, dict) or set(raw) != required:
            raise SequenceCNNContractError("expected unit schema differs")
        fold, seed = raw["outer_fold"], raw["seed"]
        if (
            isinstance(fold, bool)
            or not isinstance(fold, int)
            or fold not in range(5)
            or isinstance(seed, bool)
            or not isinstance(seed, int)
            or seed < 0
            or raw["split_id"] != f"donor{fold}_genomic{fold}"
            or re.fullmatch(
                rf"seq__sequence_cnn_control__{re.escape(str(raw['lineage_id']))}__d{fold}g{fold}__s{seed}",
                str(raw["logical_work_id"]),
            )
            is None
        ):
            raise SequenceCNNContractError("unit donor/genomic split or seed differs")
        logical_id = str(raw["logical_work_id"])
        if logical_id in seen:
            raise SequenceCNNContractError("logical unit is duplicated")
        seen.add(logical_id)
        result.append(dict(raw))
    return result


def read_source_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != SOURCE_FIELDS:
            raise SequenceCNNContractError("source manifest schema differs")
        rows = [dict(row) for row in reader]
    if len(rows) != 5:
        raise SequenceCNNContractError("source manifest row count differs")
    return rows


def audit(*, config_path: Path, project_root: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise SequenceCNNContractError("refusing to overwrite coverage contract")
    config = load_json(config_path)
    safety = config.get("safety", {})
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != "prespecified_partial_development_evaluation"
        or config.get("dataset_id") != "gse296875"
        or config.get("evaluation_role") != "valid"
        or safety
        != {
            "partial_development_only": True,
            "expected_units": 8,
            "expected_independent_evaluations": 8,
            "donor_safe_diagonal_splits_only": True,
            "actual_receipt_seed_required": True,
            "one_candidate_per_evaluator_call": True,
            "coverage_may_be_filled_post_hoc": False,
            "cross_model_ranking_authorized": False,
            "cross_family_ranking_authorized": False,
            "cross_unit_metric_aggregation_authorized": False,
            "promotion_authorized": False,
            "test_outcome_access_authorized": False,
            "sealed_data_access_authorized": False,
            "external_claim_authorized": False,
            "champion_claim_authorized": False,
        }
    ):
        raise SequenceCNNContractError("partial evaluation safety contract differs")
    expected_units = validate_expected_units(config.get("expected_units"))
    source_binding = config["source_manifest"]
    source_path = resolve_project_path(project_root, source_binding["path"])
    if sha256_file(source_path) != source_binding["sha256"]:
        raise SequenceCNNContractError("source manifest hash differs")
    source_rows = read_source_manifest(source_path)

    bundle_binding = config["bundle"]
    bundle_root = resolve_project_path(project_root, bundle_binding["path"])
    if sha256_file(bundle_root / "ARTIFACTS.json") != bundle_binding["artifacts_sha256"]:
        raise SequenceCNNContractError("bundle ARTIFACTS hash differs")
    verify_frozen_tree(bundle_root)
    bundle_meta = load_json(bundle_root / "ARTIFACTS.json")["metadata"]
    bundle_summary = load_json(bundle_root / "bundle_summary.json")
    if (
        bundle_meta.get("artifact_class") != "sequence_gpu_bundle_allocation_receipt"
        or bundle_meta.get("bundle_id") != bundle_binding["bundle_id"]
        or bundle_meta.get("slurm_job_id") != bundle_binding["slurm_job_id"]
        or bundle_meta.get("logical_tasks") != 8
        or bundle_meta.get("terminal_status") != "passed"
        or bundle_meta.get("bundle_plan_artifacts_sha256")
        != bundle_binding["bundle_plan_artifacts_sha256"]
        or bundle_summary.get("logical_tasks") != 8
        or bundle_summary.get("status_counts") != {"passed": 8}
        or bundle_summary.get("terminal_status") != "passed"
    ):
        raise SequenceCNNContractError("completed bundle contract differs")
    prediction_rows = [row for row in source_rows if row["artifact_role"] == "prediction_bundle"]
    input_rows = [row for row in source_rows if row["artifact_role"] == "evaluator_input_bundle"]
    expected_lineages = {row["lineage_id"] for row in expected_units}
    if (
        len(prediction_rows) != 1
        or prediction_rows[0]["lineage_id"] != "all"
        or prediction_rows[0]["split_id"] != "all"
        or resolve_project_path(project_root, prediction_rows[0]["root"]) != bundle_root
        or prediction_rows[0]["artifacts_sha256"] != bundle_binding["artifacts_sha256"]
        or {row["lineage_id"] for row in input_rows} != expected_lineages
        or any(row["split_id"] != "all" for row in input_rows)
    ):
        raise SequenceCNNContractError("source coverage binding differs")

    receipt_paths = sorted((bundle_root / "task_receipts").glob("*.json"))
    expected_ids = {row["logical_work_id"] for row in expected_units}
    if {path.stem for path in receipt_paths} != expected_ids:
        raise SequenceCNNContractError("bundle receipt coverage differs")
    coverage: list[dict[str, Any]] = []
    for expected in expected_units:
        receipt = load_json(
            bundle_root / "task_receipts" / f"{expected['logical_work_id']}.json"
        )
        artifact = Path(str(receipt.get("expected_artifact", ""))).resolve(strict=True)
        artifact_hash = str(receipt.get("artifacts_sha256", ""))
        if sha256_file(artifact / "ARTIFACTS.json") != artifact_hash:
            raise SequenceCNNContractError("prediction receipt hash differs")
        verify_frozen_tree(artifact)
        metadata = load_json(artifact / "ARTIFACTS.json")["metadata"]
        summary = load_json(artifact / "predictions/fixed_ccre/summary.json")
        expected_root = (
            project_root
            / "executions"
            / (
                f"sequence_cnn_control-{expected['lineage_id']}-{expected['split_id']}"
                f"-seed{expected['seed']}-{bundle_binding['slurm_job_id']}"
            )
        ).resolve(strict=True)
        if (
            receipt.get("status") != "passed"
            or receipt.get("return_code") != 0
            or receipt.get("model_id") != "sequence_cnn_control"
            or receipt.get("stage") != "train_predict"
            or receipt.get("logical_work_id") != expected["logical_work_id"]
            or receipt.get("lineage_id") != expected["lineage_id"]
            or receipt.get("outer_fold") != expected["outer_fold"]
            or receipt.get("seed") != expected["seed"]
            or artifact != expected_root
            or metadata.get("artifact_class") != "sequence_control_production"
            or metadata.get("model_id") != "sequence_cnn_control"
            or metadata.get("dataset_id") != "gse296875"
            or metadata.get("lineage_id") != expected["lineage_id"]
            or metadata.get("split_id") != expected["split_id"]
            or metadata.get("seed") != expected["seed"]
            or metadata.get("status") != "passed"
            or metadata.get("held_donor_atac_exposed") is not False
            or metadata.get("test_outcomes_used") is not False
            or summary.get("status") != "pass"
            or summary.get("model_id") != "sequence_cnn_control"
            or summary.get("seed") != expected["seed"]
            or summary.get("observed_atac_input_exposed") is not False
            or summary.get("benchmark_metrics_calculated") is not False
        ):
            raise SequenceCNNContractError("prediction unit identity or firewall differs")
        coverage.append(
            {
                **expected,
                "prediction_root": str(artifact),
                "prediction_artifacts_sha256": artifact_hash,
            }
        )
    result = {
        "schema_version": "masld-bench-sequence-cnn-eight-unit-coverage-v1",
        "status": "pass_complete_eligible_eight_unit_partial_coverage",
        "dataset_id": "gse296875",
        "evaluation_role": "valid",
        "bundle_artifacts_sha256": bundle_binding["artifacts_sha256"],
        "units": coverage,
        "unit_count": 8,
        "independent_evaluations_planned": 8,
        "actual_seeds": sorted({row["seed"] for row in coverage}),
        "partial_development_only": True,
        "coverage_filled_post_hoc": False,
        "one_candidate_per_evaluator_call": True,
        "cross_model_ranking_authorized": False,
        "cross_family_ranking_authorized": False,
        "cross_unit_metric_aggregation_authorized": False,
        "promotion_authorized": False,
        "test_outcomes_authorized": False,
        "sealed_data_authorized": False,
        "external_claim_authorized": False,
        "champion_claim_authorized": False,
        "config_sha256": sha256_file(config_path),
        "source_manifest_sha256": sha256_file(source_path),
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = audit(
        config_path=arguments.config.resolve(strict=True),
        project_root=arguments.project_root.resolve(strict=True),
        output=arguments.output,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
