#!/usr/bin/env python3
"""Freeze an outcome-safe evaluation plan from exact sequence-model receipts."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Iterable

from masld_bench.artifacts import verify_frozen_tree


SOURCE_FIELDS = (
    "artifact_role",
    "lineage_id",
    "split_id",
    "root",
    "artifacts_sha256",
)
PREDICTION_FIELDS = (
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
    "prediction_manifest",
    "prediction_manifest_sha256",
    "observed_root",
    "observed_artifacts_sha256",
    "baseline_root",
    "baseline_artifacts_sha256",
    "candidate_count",
    "candidate_ids",
)
ALLOWED_ROLES = {
    "evaluator_input_bundle",
    "prediction_bundle",
    "prediction_artifact",
}
MODEL_VIEWS = {
    "bpnet": (("bpnet", "predictions"),),
    "chrombpnet": (
        ("chrombpnet_full", "predictions/full_model"),
        ("chrombpnet_nobias", "predictions/nobias_model"),
    ),
    "sequence_cnn_control": (
        ("sequence_cnn_control", "predictions/fixed_ccre"),
    ),
    "sequence_transformer_control": (
        ("sequence_transformer_control", "predictions/fixed_ccre"),
    ),
}


class SequenceMatrixError(RuntimeError):
    """Raised when a receipt or donor-safe evaluation binding differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SequenceMatrixError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def write_tsv(
    path: Path, fields: Iterable[str], rows: Iterable[dict[str, object]]
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


def _resolve_root(project_root: Path, value: str) -> Path:
    configured = Path(value)
    root = configured if configured.is_absolute() else project_root / configured
    root = root.resolve(strict=True)
    execution_root = (project_root / "executions").resolve(strict=True)
    if execution_root not in root.parents:
        raise SequenceMatrixError(f"artifact is outside executions: {root}")
    return root


def _verify_artifact(root: Path, expected: str) -> dict[str, Any]:
    if sha256_file(root / "ARTIFACTS.json") != expected:
        raise SequenceMatrixError(f"ARTIFACTS hash differs: {root}")
    verify_frozen_tree(root)
    payload = json.loads((root / "ARTIFACTS.json").read_text(encoding="utf-8"))
    return payload.get("metadata", {})


def _prediction_receipts(bundle: Path) -> list[tuple[Path, str]]:
    receipt_dir = bundle / "task_receipts"
    if not receipt_dir.is_dir():
        raise SequenceMatrixError(f"prediction bundle lacks task receipts: {bundle}")
    artifacts: list[tuple[Path, str]] = []
    for receipt in sorted(receipt_dir.iterdir()):
        if receipt.suffix == ".json":
            row = json.loads(receipt.read_text(encoding="utf-8"))
        elif receipt.suffix == ".tsv":
            fields, rows = read_tsv(receipt)
            if len(rows) != 1 or "status" not in fields:
                raise SequenceMatrixError(f"bundle receipt schema differs: {receipt}")
            row = rows[0]
        else:
            continue
        if row.get("status") not in {"passed", "recovered_passed"}:
            continue
        path_value = row.get("expected_artifact")
        digest = row.get("artifacts_sha256")
        if not isinstance(path_value, str) or not isinstance(digest, str):
            raise SequenceMatrixError(f"bundle receipt binding differs: {receipt}")
        artifacts.append((Path(path_value).resolve(strict=True), digest))
    if not artifacts:
        raise SequenceMatrixError(f"prediction bundle has no passed tasks: {bundle}")
    return artifacts


def _load_evaluator_inputs(
    *,
    root: Path,
    expected_hash: str,
    expected_lineage: str,
) -> dict[str, dict[str, str]]:
    metadata = _verify_artifact(root, expected_hash)
    if (
        metadata.get("artifact_class") != "sequence_evaluator_inputs_lineage_bundle"
        or metadata.get("dataset_id") != "gse296875"
        or metadata.get("lineage_id") != expected_lineage
        or metadata.get("logical_tasks") != 5
        or metadata.get("failures") != 0
        or metadata.get("status") != "passed"
        or metadata.get("biological_unit") != "donor"
        or metadata.get("model_inference_input_eligible") is not False
    ):
        raise SequenceMatrixError(f"evaluator-input bundle contract differs: {root}")
    result: dict[str, dict[str, str]] = {}
    for fold in range(5):
        split_id = f"donor{fold}_genomic{fold}"
        receipt_path = root / "receipts" / f"{split_id}.json"
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if (
            receipt.get("status") not in {"passed", "reused_frozen_artifacts"}
            or receipt.get("split_id") != split_id
        ):
            raise SequenceMatrixError(f"evaluator-input receipt differs: {receipt_path}")
        profile = Path(receipt["profile_root"]).resolve(strict=True)
        baseline = Path(receipt["baseline_root"]).resolve(strict=True)
        result[split_id] = {
            "observed_root": str(profile),
            "observed_artifacts_sha256": receipt["profile_artifacts_sha256"],
            "baseline_root": str(baseline),
            "baseline_artifacts_sha256": receipt["baseline_artifacts_sha256"],
        }
    return result


def _prediction_views(root: Path, expected_hash: str) -> list[dict[str, str]]:
    metadata = _verify_artifact(root, expected_hash)
    model_id = metadata.get("model_id")
    if model_id not in MODEL_VIEWS:
        raise SequenceMatrixError(f"unsupported sequence model artifact: {root}")
    if (
        metadata.get("dataset_id") != "gse296875"
        or metadata.get("status") != "passed"
        or metadata.get("test_outcomes_used") is True
        or metadata.get("held_donor_atac_exposed") is True
        or metadata.get("evaluator_outcomes_exposed") is True
        or metadata.get("benchmark_metrics_calculated") is True
    ):
        raise SequenceMatrixError(f"prediction outcome firewall differs: {root}")
    lineage_id = metadata.get("lineage_id")
    split_id = metadata.get("split_id")
    seed = metadata.get("seed")
    if not isinstance(lineage_id, str) or not isinstance(split_id, str):
        raise SequenceMatrixError(f"prediction identity is incomplete: {root}")
    if not isinstance(seed, int) or seed < 0:
        raise SequenceMatrixError(f"prediction seed is invalid: {root}")
    result = []
    for view_id, subdir in MODEL_VIEWS[model_id]:
        if not (root / subdir / "summary.json").is_file():
            raise SequenceMatrixError(f"prediction view is absent: {root / subdir}")
        result.append(
            {
                "candidate_id": f"{view_id}__s{seed}",
                "root_model_id": model_id,
                "root": str(root),
                "artifacts_sha256": expected_hash,
                "prediction_subdir": subdir,
                "seed": str(seed),
                "lineage_id": lineage_id,
                "split_id": split_id,
            }
        )
    return result


def prepare(*, source_manifest: Path, project_root: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise SequenceMatrixError(f"output exists: {output}")
    fields, source_rows = read_tsv(source_manifest)
    if fields != SOURCE_FIELDS or not source_rows:
        raise SequenceMatrixError("source manifest schema/count differs")
    if any(row["artifact_role"] not in ALLOWED_ROLES for row in source_rows):
        raise SequenceMatrixError("source manifest has an unsupported artifact role")

    inputs: dict[tuple[str, str], dict[str, str]] = {}
    prediction_bindings: list[tuple[Path, str, str, str]] = []
    seen_roots: set[Path] = set()
    for row in source_rows:
        root = _resolve_root(project_root, row["root"])
        if root in seen_roots:
            raise SequenceMatrixError(f"source artifact root is duplicated: {root}")
        seen_roots.add(root)
        role = row["artifact_role"]
        if role == "evaluator_input_bundle":
            if row["split_id"] != "all":
                raise SequenceMatrixError("evaluator-input bundle split must be all")
            lineage_inputs = _load_evaluator_inputs(
                root=root,
                expected_hash=row["artifacts_sha256"],
                expected_lineage=row["lineage_id"],
            )
            for split_id, binding in lineage_inputs.items():
                key = (row["lineage_id"], split_id)
                if key in inputs:
                    raise SequenceMatrixError(f"evaluator inputs are duplicated: {key}")
                inputs[key] = binding
        elif role == "prediction_bundle":
            if row["lineage_id"] != "all" or row["split_id"] != "all":
                raise SequenceMatrixError("prediction bundle identity must be all/all")
            metadata = _verify_artifact(root, row["artifacts_sha256"])
            if metadata.get("artifact_class") != "sequence_gpu_bundle_allocation_receipt":
                raise SequenceMatrixError("prediction bundle artifact class differs")
            prediction_bindings.extend(
                (artifact, digest, "all", "all")
                for artifact, digest in _prediction_receipts(root)
            )
        else:
            if row["lineage_id"] == "all" or row["split_id"] == "all":
                raise SequenceMatrixError("prediction artifact identity must be explicit")
            prediction_bindings.append(
                (
                    root,
                    row["artifacts_sha256"],
                    row["lineage_id"],
                    row["split_id"],
                )
            )

    views_by_key: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    prediction_root_hashes: dict[Path, str] = {}
    for root, expected_hash, declared_lineage, declared_split in prediction_bindings:
        prior = prediction_root_hashes.get(root)
        if prior is not None:
            if prior != expected_hash:
                raise SequenceMatrixError(f"prediction root hash is inconsistent: {root}")
            continue
        prediction_root_hashes[root] = expected_hash
        for view in _prediction_views(root, expected_hash):
            actual_lineage = view.pop("lineage_id")
            actual_split = view.pop("split_id")
            if (
                declared_lineage not in {"all", actual_lineage}
                or declared_split not in {"all", actual_split}
            ):
                raise SequenceMatrixError(f"prediction source identity differs: {root}")
            key = (actual_lineage, actual_split)
            views_by_key[key].append(view)

    output.mkdir(mode=0o750)
    plan_rows: list[dict[str, object]] = []
    candidate_coverage: dict[str, list[str]] = defaultdict(list)
    for key in sorted(views_by_key):
        lineage_id, split_id = key
        if key not in inputs:
            raise SequenceMatrixError(f"prediction lacks evaluator inputs: {key}")
        binding = inputs[key]
        profile_metadata = _verify_artifact(
            Path(binding["observed_root"]), binding["observed_artifacts_sha256"]
        )
        baseline_metadata = _verify_artifact(
            Path(binding["baseline_root"]), binding["baseline_artifacts_sha256"]
        )
        identity = ("gse296875", split_id, lineage_id)
        identity_fields = ("dataset_id", "split_id", "lineage_id")
        if (
            tuple(profile_metadata.get(field) for field in identity_fields) != identity
            or tuple(baseline_metadata.get(field) for field in identity_fields) != identity
            or profile_metadata.get("model_inference_input_eligible") is not False
            or baseline_metadata.get("held_donor_outcomes_exposed") is not False
        ):
            raise SequenceMatrixError(f"evaluator-input identity differs: {key}")
        views = sorted(views_by_key[key], key=lambda row: row["candidate_id"])
        candidate_ids = [row["candidate_id"] for row in views]
        if len(set(candidate_ids)) != len(candidate_ids):
            raise SequenceMatrixError(f"candidate identity is duplicated: {key}")
        relative_manifest = Path("manifests") / lineage_id / f"{split_id}.tsv"
        manifest_path = output / relative_manifest
        write_tsv(manifest_path, PREDICTION_FIELDS, views)
        for candidate_id in candidate_ids:
            candidate_coverage[candidate_id].append(f"{lineage_id}:{split_id}")
        plan_rows.append(
            {
                "dataset_id": "gse296875",
                "lineage_id": lineage_id,
                "split_id": split_id,
                "prediction_manifest": relative_manifest.as_posix(),
                "prediction_manifest_sha256": sha256_file(manifest_path),
                **binding,
                "candidate_count": len(views),
                "candidate_ids": ",".join(candidate_ids),
            }
        )
    if not plan_rows:
        raise SequenceMatrixError("no evaluable prediction groups were found")
    write_tsv(output / "matrix.tsv", PLAN_FIELDS, plan_rows)
    summary = {
        "schema_version": "masld-bench-sequence-evaluation-matrix-plan-v1",
        "status": "pass",
        "dataset_id": "gse296875",
        "evaluation_role": "valid",
        "evaluation_groups": len(plan_rows),
        "prediction_artifacts": len(prediction_root_hashes),
        "prediction_views": sum(len(rows) for rows in views_by_key.values()),
        "candidate_coverage": dict(sorted(candidate_coverage.items())),
        "lineages_with_inputs": sorted({lineage for lineage, _ in inputs}),
        "input_splits": len(inputs),
        "test_outcomes_authorized": False,
        "sealed_data_authorized": False,
        "cross_fold_ranking_authorized": False,
        "folds_or_seeds_are_biological_replicates": False,
        "incomplete_candidate_coverage_is_rankable": False,
        "source_manifest_sha256": sha256_file(source_manifest),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(
        source_manifest=args.source_manifest,
        project_root=args.project_root.resolve(strict=True),
        output=args.output,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
