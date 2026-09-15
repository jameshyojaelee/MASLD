#!/usr/bin/env python3
"""Independently rederive and verify blind GSE244832 baseline predictions."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import (
    freeze_tree,
    verify_frozen_tree,
    write_json_exclusive,
    write_text_exclusive,
)
from masld_bench.gse244832_task_native_baselines import MODEL_IDS, fit, load_state, predict
from scripts.commit_gse244832_atac_transport_predictions import validate_bundle
from scripts.run_gse244832_task_native_blind_baselines import (
    ROLE_INDEX,
    ROTATIONS,
    digest,
    read_json,
    read_query_axis,
    validate_config,
    validate_source,
)


class BlindBaselineVerificationError(ValueError):
    """Raised when a frozen prediction cannot be independently reproduced."""


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_runtime_lock(output: Path) -> None:
    completed = subprocess.run(
        [sys.executable, "-m", "pip", "freeze"],
        check=True,
        capture_output=True,
        text=True,
    )
    write_text_exclusive(output / "pip_freeze.txt", completed.stdout)
    write_text_exclusive(
        output / "python_version.txt",
        f"{platform.python_version()}\n{sys.executable}\n",
    )


def _index(records: Sequence[Mapping[str, Any]]) -> dict[tuple[str, str], Mapping[str, Any]]:
    output = {
        (str(record["model_id"]), str(record["rotation_id"])): record
        for record in records
    }
    if len(output) != len(records):
        raise BlindBaselineVerificationError("run child records are duplicated")
    return output


def verify(
    root: Path,
    config_path: Path,
    source_input: Path,
    source_artifacts_sha256: str,
    run_input: Path,
    run_artifacts_sha256: str,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise BlindBaselineVerificationError("refusing to overwrite verification output")
    root = root.resolve(strict=True)
    config, paths = validate_config(root, config_path)
    source = source_input.resolve(strict=True)
    run_root = run_input.resolve(strict=True)
    source.relative_to(root / "executions")
    run_root.relative_to(root / "executions")
    source_counts, source_lineages = validate_source(
        source, source_artifacts_sha256, config
    )
    run_manifest = verify_frozen_tree(run_root)
    receipt = read_json(run_root / "receipt.json")
    if (
        digest(run_root / "ARTIFACTS.json") != run_artifacts_sha256
        or run_manifest.get("metadata", {}).get("artifact_class")
        != "gse244832_task_native_blind_baseline_run"
        or receipt.get("status") != "passed"
        or receipt.get("predictions") != 6
        or receipt.get("commits") != 6
        or receipt.get("all_predictions_committed_before_outcome_open") is not True
        or receipt.get("development_outcomes_read") is not False
        or receipt.get("metrics_calculated") is not False
    ):
        raise BlindBaselineVerificationError("blind baseline run authority differs")
    prediction_records = _index(receipt["prediction_records"])
    commit_records = _index(receipt["commit_records"])
    expected = {(model, rotation) for model in MODEL_IDS for rotation in ROTATIONS}
    if set(prediction_records) != expected or set(commit_records) != expected:
        raise BlindBaselineVerificationError("blind baseline child roster differs")
    registration = read_json(paths["registration"] / "execution_contract.json")
    unit_states = np.load(
        paths["registration"] / registration["query_lineage_missing_state_path"],
        allow_pickle=False,
    )
    eligible = unit_states.reshape(-1) == 0
    output.mkdir(parents=True, mode=0o750)
    verification_rows: list[dict[str, Any]] = []
    for rotation_id, (context_role, target_role) in ROTATIONS.items():
        source_context = np.asarray(
            source_counts[:, :, ROLE_INDEX[context_role], :], dtype=np.float64
        ).reshape(156, 16_000)
        source_target = np.asarray(
            source_counts[:, :, ROLE_INDEX[target_role], :], dtype=np.float64
        ).reshape(156, 16_000)
        query_root = paths[f"query_{context_role}"]
        query_lineages = read_query_axis(query_root / "query_axis.tsv")
        query_memmap = np.load(
            query_root / "query_unit_fragment_total.uint32.npy",
            mmap_mode="r",
            allow_pickle=False,
        )
        query_context = np.asarray(
            query_memmap.reshape(72, 16_000)[eligible], dtype=np.float64
        )
        eligible_lineages = query_lineages[eligible]
        for model_id in MODEL_IDS:
            key = (model_id, rotation_id)
            prediction_record = prediction_records[key]
            prediction_root = (
                run_root / str(prediction_record["path"])
            ).resolve(strict=True)
            prediction_root.relative_to(run_root)
            prediction_manifest = verify_frozen_tree(prediction_root)
            if (
                digest(prediction_root / "ARTIFACTS.json")
                != prediction_record["artifacts_sha256"]
                or prediction_manifest.get("metadata", {}).get("artifact_class")
                != "gse244832_task_native_blind_prediction"
            ):
                raise BlindBaselineVerificationError("prediction child authority differs")
            frozen = np.load(
                prediction_root / "predictions.float32.npy",
                mmap_mode="r",
                allow_pickle=False,
            )
            missing = np.load(
                prediction_root / "missing_state.uint8.npy",
                mmap_mode="r",
                allow_pickle=False,
            )
            if (
                frozen.shape != (18, 4, 16_000)
                or missing.shape != frozen.shape
                or not np.array_equal(
                    missing,
                    np.broadcast_to(unit_states[:, :, None], missing.shape),
                )
                or not np.isnan(frozen.reshape(72, 16_000)[~eligible]).all()
            ):
                raise BlindBaselineVerificationError("prediction missingness differs")
            restored = load_state(prediction_root / "fit_state")
            restored_prediction = predict(
                restored,
                context_counts=query_context,
                lineage_indices=eligible_lineages,
            )
            refitted = fit(
                model_id,
                context_counts=source_context,
                target_counts=source_target,
                lineage_indices=source_lineages,
                lsi_components=int(config["fit"]["lsi_components"]),
                ridge_alpha=float(config["fit"]["ridge_alpha"]),
                standardized_clip=float(config["fit"]["standardized_clip"]),
            )
            refitted_prediction = predict(
                refitted,
                context_counts=query_context,
                lineage_indices=eligible_lineages,
            )
            observed = np.asarray(
                frozen.reshape(72, 16_000)[eligible], dtype=np.float64
            )
            state_max = float(np.max(np.abs(observed - restored_prediction)))
            refit_max = float(np.max(np.abs(observed - refitted_prediction)))
            tolerance = 1.0e-3 + 1.0e-5 * float(np.max(np.abs(observed)))
            if state_max > tolerance or refit_max > tolerance:
                raise BlindBaselineVerificationError("blind prediction reproduction differs")
            structural_receipt = validate_bundle(
                prediction_root=prediction_root,
                prediction_artifacts_sha256=str(prediction_record["artifacts_sha256"]),
                axis_root=paths["axis"],
                axis_artifacts_sha256=config["exchange_axis"]["artifacts_sha256"],
                registration_root=paths["registration"],
                registration_artifacts_sha256=config["execution_registration"]["artifacts_sha256"],
            )
            commit_record = commit_records[key]
            commit_root = (run_root / str(commit_record["path"])).resolve(strict=True)
            commit_root.relative_to(run_root)
            commit_manifest = verify_frozen_tree(commit_root)
            commit_receipt = read_json(commit_root / "commit.json")
            if (
                digest(commit_root / "ARTIFACTS.json")
                != commit_record["artifacts_sha256"]
                or commit_manifest.get("metadata", {}).get("artifact_class")
                != "gse244832_atac_transport_prediction_commit"
                or commit_receipt != structural_receipt
                or commit_receipt.get("development_outcomes_read") is not False
                or commit_receipt.get("metrics_calculated") is not False
            ):
                raise BlindBaselineVerificationError("prediction commit differs")
            verification_rows.append(
                {
                    "model_id": model_id,
                    "rotation_id": rotation_id,
                    "eligible_units": int(eligible.sum()),
                    "ineligible_units": int((~eligible).sum()),
                    "state_reproduction_max_abs": f"{state_max:.12g}",
                    "refit_reproduction_max_abs": f"{refit_max:.12g}",
                    "tolerance": f"{tolerance:.12g}",
                    "prediction_artifacts_sha256": prediction_record["artifacts_sha256"],
                    "commit_artifacts_sha256": commit_record["artifacts_sha256"],
                    "status": "passed",
                }
            )
    write_tsv(
        output / "verification.tsv",
        (
            "model_id",
            "rotation_id",
            "eligible_units",
            "ineligible_units",
            "state_reproduction_max_abs",
            "refit_reproduction_max_abs",
            "tolerance",
            "prediction_artifacts_sha256",
            "commit_artifacts_sha256",
            "status",
        ),
        verification_rows,
    )
    write_json_exclusive(
        output / "receipt.json",
        {
            "schema_version": "masld-bench-gse244832-task-native-blind-baseline-verification-v1",
            "status": "passed",
            "verified_predictions": len(verification_rows),
            "verified_commits": len(verification_rows),
            "source_refit_reproduced": True,
            "frozen_state_reproduced": True,
            "mask_applied_before_query_transform": True,
            "eligible_units": 48,
            "structurally_missing_units": 4,
            "below_qc_units": 20,
            "evaluator_outcomes_opened": False,
            "condition_or_phenotype_read": False,
            "development_outcomes_read": False,
            "sealed_data_read": False,
            "metrics_calculated": False,
            "scoring_performed": False,
            "champion_claim_allowed": False,
        },
    )
    write_runtime_lock(output)
    artifact_sha = freeze_tree(
        output,
        {
            "artifact_class": "gse244832_task_native_blind_baseline_verification",
            "verified_predictions": 6,
            "verified_commits": 6,
            "evaluator_outcomes_opened": False,
            "metrics_calculated": False,
            "status": "passed",
        },
    )
    return {"output": str(output), "artifacts_sha256": artifact_sha, "verified": 6}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", dest="config_path", required=True, type=Path)
    parser.add_argument("--source-input", required=True, type=Path)
    parser.add_argument("--source-artifacts-sha256", required=True)
    parser.add_argument("--run-input", required=True, type=Path)
    parser.add_argument("--run-artifacts-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    print(json.dumps(verify(**vars(arguments)), sort_keys=True))


if __name__ == "__main__":
    main()
