#!/usr/bin/env python3
"""Bind an eight-unit ChromBPNet bundle to separate development evaluators."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import verify_frozen_tree


CONFIG_SCHEMA = "masld-bench-chrombpnet-bundle-partial-evaluation-v1"
MANIFEST_FIELDS = (
    "candidate_id",
    "root_model_id",
    "root",
    "artifacts_sha256",
    "prediction_subdir",
    "seed",
)
PLAN_FIELDS = (
    "dataset_id",
    "lineage_id",
    "split_id",
    "outer_fold",
    "seed",
    "logical_work_id",
    "view_id",
    "prediction_manifest",
    "prediction_manifest_sha256",
    "prediction_root",
    "prediction_artifacts_sha256",
    "observed_root",
    "observed_artifacts_sha256",
    "baseline_root",
    "baseline_artifacts_sha256",
)
VIEWS = {
    "chrombpnet_full": "predictions/full_model",
    "chrombpnet_nobias": "predictions/nobias_model",
}


class PartialEvaluationError(RuntimeError):
    """Raised when bundle, split, seed, or outcome separations differ."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise PartialEvaluationError(f"JSON object differs: {path}")
    return value


def write_tsv(
    path: Path,
    fields: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def resolve_project_path(root: Path, value: str) -> Path:
    configured = Path(value)
    if configured.is_absolute() or ".." in configured.parts:
        raise PartialEvaluationError("unsafe project-relative path")
    path = (root / configured).resolve(strict=True)
    path.relative_to(root)
    return path


def validate_expected_units(units: object) -> list[dict[str, Any]]:
    if not isinstance(units, list) or len(units) != 8:
        raise PartialEvaluationError("expected unit count differs")
    required = {"logical_work_id", "lineage_id", "outer_fold", "split_id", "seed"}
    result: list[dict[str, Any]] = []
    identities: set[str] = set()
    for raw in units:
        if not isinstance(raw, dict) or set(raw) != required:
            raise PartialEvaluationError("expected unit schema differs")
        fold = raw["outer_fold"]
        seed = raw["seed"]
        split_id = raw["split_id"]
        if (
            isinstance(fold, bool)
            or not isinstance(fold, int)
            or fold not in range(5)
            or isinstance(seed, bool)
            or not isinstance(seed, int)
            or seed < 0
            or split_id != f"donor{fold}_genomic{fold}"
            or re.fullmatch(
                rf"seq__chrombpnet__{re.escape(str(raw['lineage_id']))}__d{fold}g{fold}__s{seed}",
                str(raw["logical_work_id"]),
            )
            is None
        ):
            raise PartialEvaluationError("donor/genomic split or seed differs")
        identity = str(raw["logical_work_id"])
        if identity in identities:
            raise PartialEvaluationError("expected logical work ID is duplicated")
        identities.add(identity)
        result.append(dict(raw))
    return result


def candidate_row(
    *,
    view_id: str,
    prediction_subdir: str,
    root: Path,
    artifacts_sha256: str,
    seed: int,
) -> dict[str, str]:
    if VIEWS.get(view_id) != prediction_subdir:
        raise PartialEvaluationError("native view binding differs")
    return {
        "candidate_id": f"{view_id}__s{seed}",
        "root_model_id": "chrombpnet",
        "root": str(root),
        "artifacts_sha256": artifacts_sha256,
        "prediction_subdir": prediction_subdir,
        "seed": str(seed),
    }


def verify_tree(root: Path, expected_hash: str) -> dict[str, Any]:
    if sha256_file(root / "ARTIFACTS.json") != expected_hash:
        raise PartialEvaluationError(f"ARTIFACTS hash differs: {root}")
    verify_frozen_tree(root)
    return load_json(root / "ARTIFACTS.json").get("metadata", {})


def load_lineage_inputs(
    *,
    project_root: Path,
    lineage_id: str,
    binding: Mapping[str, str],
) -> dict[str, dict[str, str]]:
    root = resolve_project_path(project_root, binding["path"])
    metadata = verify_tree(root, binding["artifacts_sha256"])
    if (
        metadata.get("artifact_class") != "sequence_evaluator_inputs_lineage_bundle"
        or metadata.get("dataset_id") != "gse296875"
        or metadata.get("lineage_id") != lineage_id
        or metadata.get("logical_tasks") != 5
        or metadata.get("failures") != 0
        or metadata.get("status") != "passed"
        or metadata.get("biological_unit") != "donor"
        or metadata.get("model_inference_input_eligible") is not False
    ):
        raise PartialEvaluationError("evaluator input bundle differs")
    result: dict[str, dict[str, str]] = {}
    for fold in range(5):
        split_id = f"donor{fold}_genomic{fold}"
        receipt = load_json(root / "receipts" / f"{split_id}.json")
        if (
            receipt.get("status") not in {"passed", "reused_frozen_artifacts"}
            or receipt.get("split_id") != split_id
        ):
            raise PartialEvaluationError("evaluator input receipt differs")
        observed = Path(str(receipt["profile_root"])).resolve(strict=True)
        baseline = Path(str(receipt["baseline_root"])).resolve(strict=True)
        observed_meta = verify_tree(observed, str(receipt["profile_artifacts_sha256"]))
        baseline_meta = verify_tree(baseline, str(receipt["baseline_artifacts_sha256"]))
        identity = ("gse296875", split_id, lineage_id)
        fields = ("dataset_id", "split_id", "lineage_id")
        if (
            tuple(observed_meta.get(field) for field in fields) != identity
            or tuple(baseline_meta.get(field) for field in fields) != identity
            or observed_meta.get("artifact_class")
            != "chrombpnet_development_base_profiles"
            or observed_meta.get("model_inference_input_eligible") is not False
            or baseline_meta.get("artifact_class")
            != "sequence_training_mean_baseline"
            or baseline_meta.get("held_donor_outcomes_exposed") is not False
        ):
            raise PartialEvaluationError("evaluator input identity differs")
        result[split_id] = {
            "observed_root": str(observed),
            "observed_artifacts_sha256": str(receipt["profile_artifacts_sha256"]),
            "baseline_root": str(baseline),
            "baseline_artifacts_sha256": str(receipt["baseline_artifacts_sha256"]),
        }
    return result


def prepare(*, config_path: Path, project_root: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise PartialEvaluationError(f"output exists: {output}")
    config = load_json(config_path)
    safety = config.get("safety", {})
    if (
        config.get("schema_version") != CONFIG_SCHEMA
        or config.get("status") != "prespecified_partial_development_evaluation"
        or config.get("dataset_id") != "gse296875"
        or config.get("evaluation_role") != "valid"
        or safety
        != {
            "partial_development_only": True,
            "expected_units": 8,
            "expected_independent_evaluations": 16,
            "donor_safe_diagonal_splits_only": True,
            "actual_artifact_seed_required": True,
            "one_candidate_view_per_evaluator_call": True,
            "cross_view_ranking_authorized": False,
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
        raise PartialEvaluationError("partial evaluation safety contract differs")
    expected_units = validate_expected_units(config.get("expected_units"))
    views = config.get("views")
    if (
        not isinstance(views, list)
        or {row.get("view_id"): row.get("prediction_subdir") for row in views if isinstance(row, dict)}
        != VIEWS
    ):
        raise PartialEvaluationError("two-view roster differs")

    bundle_binding = config["bundle"]
    plan_root = resolve_project_path(project_root, bundle_binding["bundle_plan_path"])
    verify_tree(plan_root, bundle_binding["bundle_plan_artifacts_sha256"])
    bundle_root = resolve_project_path(project_root, bundle_binding["path"])
    bundle_hash = sha256_file(bundle_root / "ARTIFACTS.json")
    bundle_metadata = verify_tree(bundle_root, bundle_hash)
    bundle_summary = load_json(bundle_root / "bundle_summary.json")
    if (
        bundle_metadata.get("artifact_class") != "sequence_gpu_bundle_allocation_receipt"
        or bundle_metadata.get("bundle_id") != bundle_binding["bundle_id"]
        or bundle_metadata.get("slurm_job_id") != bundle_binding["slurm_job_id"]
        or bundle_metadata.get("logical_tasks") != 8
        or bundle_metadata.get("terminal_status") != "passed"
        or bundle_metadata.get("bundle_plan_artifacts_sha256")
        != bundle_binding["bundle_plan_artifacts_sha256"]
        or bundle_summary.get("status_counts") != {"passed": 8}
        or bundle_summary.get("terminal_status") != "passed"
    ):
        raise PartialEvaluationError("bundle completion contract differs")
    receipt_paths = sorted((bundle_root / "task_receipts").glob("*.json"))
    expected_ids = {row["logical_work_id"] for row in expected_units}
    if {path.stem for path in receipt_paths} != expected_ids:
        raise PartialEvaluationError("bundle task receipt roster differs")

    input_bundles = config.get("evaluator_input_bundles")
    expected_lineages = {row["lineage_id"] for row in expected_units}
    if not isinstance(input_bundles, dict) or set(input_bundles) != expected_lineages:
        raise PartialEvaluationError("evaluator input lineage roster differs")
    inputs = {
        lineage_id: load_lineage_inputs(
            project_root=project_root,
            lineage_id=lineage_id,
            binding=input_bundles[lineage_id],
        )
        for lineage_id in sorted(expected_lineages)
    }

    plan_rows: list[dict[str, object]] = []
    artifact_hashes: dict[str, str] = {}
    for expected in expected_units:
        receipt = load_json(
            bundle_root / "task_receipts" / f"{expected['logical_work_id']}.json"
        )
        prediction_root = Path(str(receipt.get("expected_artifact", ""))).resolve(strict=True)
        prediction_hash = str(receipt.get("artifacts_sha256", ""))
        metadata = verify_tree(prediction_root, prediction_hash)
        expected_root = (
            project_root
            / "executions"
            / (
                f"chrombpnet-{expected['lineage_id']}-{expected['split_id']}"
                f"-seed{expected['seed']}-{bundle_binding['slurm_job_id']}"
            )
        ).resolve(strict=True)
        identity = ("gse296875", expected["split_id"], expected["lineage_id"])
        fields = ("dataset_id", "split_id", "lineage_id")
        if (
            receipt.get("status") != "passed"
            or receipt.get("return_code") != 0
            or receipt.get("model_id") != "chrombpnet"
            or receipt.get("stage") != "train_predict"
            or receipt.get("logical_work_id") != expected["logical_work_id"]
            or receipt.get("lineage_id") != expected["lineage_id"]
            or receipt.get("outer_fold") != expected["outer_fold"]
            or receipt.get("seed") != expected["seed"]
            or prediction_root != expected_root
            or metadata.get("artifact_class")
            != "chrombpnet_full_depth_training_and_predictions"
            or metadata.get("model_id") != "chrombpnet"
            or tuple(metadata.get(field) for field in fields) != identity
            or metadata.get("seed") != expected["seed"]
            or metadata.get("status") != "passed"
            or metadata.get("benchmark_metrics_calculated") is not False
            or metadata.get("test_outcomes_used") is not False
            or metadata.get("champion_eligible") is not False
        ):
            raise PartialEvaluationError("prediction unit identity or firewall differs")
        artifact_hashes[expected["logical_work_id"]] = prediction_hash
        binding = inputs[expected["lineage_id"]][expected["split_id"]]
        for view in sorted(views, key=lambda row: row["view_id"]):
            candidate = candidate_row(
                view_id=view["view_id"],
                prediction_subdir=view["prediction_subdir"],
                root=prediction_root,
                artifacts_sha256=prediction_hash,
                seed=expected["seed"],
            )
            relative_manifest = (
                Path("manifests")
                / expected["lineage_id"]
                / expected["split_id"]
                / f"{view['view_id']}.tsv"
            )
            manifest_path = output / relative_manifest
            write_tsv(manifest_path, MANIFEST_FIELDS, [candidate])
            plan_rows.append(
                {
                    "dataset_id": "gse296875",
                    "lineage_id": expected["lineage_id"],
                    "split_id": expected["split_id"],
                    "outer_fold": expected["outer_fold"],
                    "seed": expected["seed"],
                    "logical_work_id": expected["logical_work_id"],
                    "view_id": view["view_id"],
                    "prediction_manifest": relative_manifest.as_posix(),
                    "prediction_manifest_sha256": sha256_file(manifest_path),
                    "prediction_root": str(prediction_root),
                    "prediction_artifacts_sha256": prediction_hash,
                    **binding,
                }
            )
    if len(plan_rows) != 16:
        raise PartialEvaluationError("independent evaluation count differs")
    output.mkdir(parents=True, exist_ok=True)
    write_tsv(output / "evaluation_plan.tsv", PLAN_FIELDS, plan_rows)
    summary = {
        "schema_version": "masld-bench-chrombpnet-bundle-partial-plan-v1",
        "status": "pass",
        "dataset_id": "gse296875",
        "evaluation_role": "valid",
        "bundle_artifacts_sha256": bundle_hash,
        "unit_artifacts_sha256": dict(sorted(artifact_hashes.items())),
        "units": 8,
        "independent_evaluations": 16,
        "views": sorted(VIEWS),
        "actual_seeds": sorted({row["seed"] for row in expected_units}),
        "one_candidate_view_per_evaluator_call": True,
        "cross_view_ranking_authorized": False,
        "cross_model_ranking_authorized": False,
        "cross_family_ranking_authorized": False,
        "cross_unit_metric_aggregation_authorized": False,
        "promotion_authorized": False,
        "test_outcomes_authorized": False,
        "sealed_data_authorized": False,
        "champion_claim_authorized": False,
        "config_sha256": sha256_file(config_path),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = prepare(
        config_path=arguments.config.resolve(strict=True),
        project_root=arguments.project_root.resolve(strict=True),
        output=arguments.output,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
