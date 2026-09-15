#!/usr/bin/env python3
"""Launch only a fully included native Borzoi prediction entrypoint."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
from typing import Any

from masld_bench.artifacts import write_json_exclusive
from scripts.preflight_gse281364_borzoi_native_execution import (
    BorzoiExecutionError,
    load_json,
    preflight,
    project_path,
    sha256_file,
)


def validate_prediction_bundle(output: Path) -> dict[str, Any]:
    required = {
        "manifest.tsv",
        "member_track_deltas.npz",
        "receipt.json",
    }
    for relative in required:
        path = output / relative
        if path.is_symlink() or not path.is_file():
            raise BorzoiExecutionError(f"native prediction output is absent: {relative}")
    receipt = load_json(output / "receipt.json", label="native prediction receipt")
    if (
        receipt.get("schema_version")
        != "masld-bench-gse281364-borzoi-native-predictions-v1"
        or receipt.get("status") != "pass_native_predictions"
        or receipt.get("dataset_id") != "gse281364"
        or receipt.get("model_id") != "borzoi_ensemble"
        or receipt.get("fixture_elements") != 1_033
        or receipt.get("ensemble_members") != [0, 1, 2, 3]
        or receipt.get("orientations") != ["forward", "reverse_complement"]
        or receipt.get("orientation_predictions_per_allele") != 8
        or receipt.get("shift_bp") != [0]
        or receipt.get("human_tracks") != 7_611
        or receipt.get("released_prediction_bins") != 16_352
        or receipt.get("allele_effect_sign") != "ALT_minus_REF"
    ):
        raise BorzoiExecutionError("native prediction receipt identity differs")
    for field in (
        "outcomes_read",
        "reporter_counts_read",
        "sealed_assets_read",
        "metrics_calculated",
        "head_fit",
    ):
        if receipt.get(field) is not False:
            raise BorzoiExecutionError(f"native prediction firewall differs: {field}")
    with (output / "manifest.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1_033 or not rows or "fixture_id" not in rows[0]:
        raise BorzoiExecutionError("native prediction manifest census differs")
    fixture_ids = [row["fixture_id"] for row in rows]
    if len(set(fixture_ids)) != 1_033 or any(not value for value in fixture_ids):
        raise BorzoiExecutionError("native prediction fixture IDs differ")
    return {
        "receipt_sha256": sha256_file(output / "receipt.json"),
        "manifest_sha256": sha256_file(output / "manifest.tsv"),
        "member_track_deltas_sha256": sha256_file(
            output / "member_track_deltas.npz"
        ),
        "fixture_elements": len(rows),
    }


def run(*, root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    if output.exists() or output.is_symlink():
        raise BorzoiExecutionError(f"refusing to overwrite {output}")
    audit = preflight(root=root, config_path=config_path)
    if not audit["executable"] or audit["status"] != "pass_executable":
        raise BorzoiExecutionError(
            "native Borzoi execution remains blocked: " + ",".join(audit["blockers"])
        )
    root = root.resolve(strict=True)
    config = load_json(config_path.resolve(strict=True), label="Borzoi execution config")
    gates = config["required_execution_artifacts"]
    fixed = config["frozen_authorities"]
    runtime_root = project_path(root, gates["native_runtime"]["path"], must_exist=True)
    checkpoint_root = project_path(
        root, gates["official_checkpoints"]["path"], must_exist=True
    )
    runtime = audit["execution_gates"]["native_runtime"]["details"]
    python = runtime_root / runtime["python_executable"]["path"]
    entrypoint = runtime_root / runtime["prediction_entrypoint"]["path"]
    fixture = project_path(root, fixed["native_fixture"]["path"], must_exist=True)
    reference = (
        project_path(root, fixed["project_reference"]["path"], must_exist=True)
        / "GRCh38.p14.sequence_model.fa"
    )
    targets = (
        project_path(root, fixed["admission_sources"]["path"], must_exist=True)
        / "sources/borzoi_targets_human.txt.gz"
    )
    command = [
        str(python),
        str(entrypoint),
        "--contract-version",
        "masld-bench-borzoi-native-predictor-cli-v1",
        "--fixture-root",
        str(fixture / "fixture"),
        "--checkpoint-root",
        str(checkpoint_root),
        "--reference-fasta",
        str(reference),
        "--target-manifest",
        str(targets),
        "--output",
        str(output),
        "--member-order",
        "0,1,2,3",
        "--orientations",
        "forward,reverse_complement",
        "--shift-bp",
        "0",
        "--offline",
    ]
    environment = os.environ.copy()
    offline_guard = root / "scripts/offline_guard"
    environment.update(
        {
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONNOUSERSITE": "1",
            "PYTHONPATH": str(offline_guard),
            "HF_HUB_OFFLINE": "1",
            "TRANSFORMERS_OFFLINE": "1",
            "WANDB_MODE": "offline",
        }
    )
    subprocess.run(command, check=True, env=environment)
    bundle = validate_prediction_bundle(output)
    write_json_exclusive(output / "execution_preflight.json", audit, mode=0o640)
    launcher = {
        "schema_version": "masld-bench-gse281364-borzoi-native-launcher-v1",
        "status": "pass_prediction_bundle",
        "dataset_id": "gse281364",
        "model_id": "borzoi_ensemble",
        "config_sha256": audit["config_sha256"],
        "runtime_artifacts_sha256": gates["native_runtime"]["artifacts_sha256"],
        "checkpoint_artifacts_sha256": gates["official_checkpoints"][
            "artifacts_sha256"
        ],
        "python_executable_sha256": runtime["python_executable"]["sha256"],
        "prediction_entrypoint_sha256": runtime["prediction_entrypoint"]["sha256"],
        "prediction_bundle": bundle,
        "checkpoint_bytes_opened_by_launcher": False,
        "model_modules_imported_by_launcher": False,
        "model_forward_executed_by_native_entrypoint": True,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "metrics_calculated": False,
    }
    write_json_exclusive(output / "launcher_receipt.json", launcher, mode=0o640)
    return launcher


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    result = run(root=arguments.root, config_path=arguments.config, output=arguments.output)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
