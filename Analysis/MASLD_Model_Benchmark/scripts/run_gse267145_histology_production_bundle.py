#!/usr/bin/env python3
"""Run and freeze the resumable 25-unit GSE267145 prediction campaign."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import (
    freeze_tree,
    publish_directory_noreplace,
    verify_frozen_tree,
    write_json_exclusive,
    write_text_exclusive,
)


OUTER_FOLDS = tuple(range(5))
SEEDS = (1701, 1709, 1721, 1723, 1733)
MODEL_COUNT = 11
FITTER_SHA256 = "b9832e9ddb1fe794103167db23467099288e3d9c351971716e043d6405d9dd27"
DIAGNOSTIC_ARTIFACTS_SHA256 = "da2e5cb0cd6da1d7f7ecb15cf9e38cd0d11a7ea5a8a0e8c989cf4f1174b6f6ef"
SURFACE_SHA256 = "23367a28bd16266c79adde854534b124cfad7d2062788625839e78f3267a8c0b"
MOLECULAR_ARTIFACTS_SHA256 = "2a60ef5d5d904f5f833833ebef75e7946e7db3b6b9aaf502cb84ddedf97eb456"
FOLDS_ARTIFACTS_SHA256 = "9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3"
FIT_VIEWS_ARTIFACTS_SHA256 = "8254274f1b11089656f4b6d5ec446a6ba491e2f28ae8abcc088fc6e78ac75d6a"
PREFLIGHT_ARTIFACTS_SHA256 = "2b95408bc92719706b728fd795ae918bb03280ff99ff8c549965318fee5d5bf3"
RUNTIME_ARTIFACTS_SHA256 = "f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb"


class ProductionBundleError(RuntimeError):
    """Raised when the campaign cannot proceed without violating being read-only."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ProductionBundleError(f"invalid JSON document: {path}: {error}") from error
    if not isinstance(value, dict):
        raise ProductionBundleError(f"JSON document is not an object: {path}")
    return value


def _check_hash(path: Path, expected: str, *, label: str) -> None:
    if not path.is_file() or path.is_symlink() or sha256_file(path) != expected:
        raise ProductionBundleError(f"{label} SHA-256 differs: {path}")


def _input_contract(arguments: argparse.Namespace) -> None:
    checks = (
        (arguments.runtime / "ARTIFACTS.json", RUNTIME_ARTIFACTS_SHA256, "runtime"),
        (arguments.molecular / "ARTIFACTS.json", MOLECULAR_ARTIFACTS_SHA256, "molecular fixture"),
        (arguments.folds / "ARTIFACTS.json", FOLDS_ARTIFACTS_SHA256, "fold fixture"),
        (arguments.fit_views / "ARTIFACTS.json", FIT_VIEWS_ARTIFACTS_SHA256, "fit views"),
        (arguments.preflight / "ARTIFACTS.json", PREFLIGHT_ARTIFACTS_SHA256, "validated preflight"),
        (arguments.diagnostic / "ARTIFACTS.json", DIAGNOSTIC_ARTIFACTS_SHA256, "corrected convergence diagnostic"),
        (arguments.surface, SURFACE_SHA256, "benchmark surface"),
        (arguments.fitter, FITTER_SHA256, "validated one-unit fitter"),
    )
    for path, expected, label in checks:
        _check_hash(path, expected, label=label)
    if arguments.fit_views.parent != arguments.preflight:
        raise ProductionBundleError("fit views are not inside the validated preflight tree")
    verify_frozen_tree(arguments.diagnostic)
    diagnostic = _read_json_object(
        arguments.diagnostic / "units/outer_1--seed_1701/stdout.json"
    )
    solver = diagnostic.get("elastic_logistic_solver")
    if (
        not isinstance(solver, dict)
        or solver.get("solver") != "saga"
        or solver.get("tolerance") != 0.0001
        or solver.get("previous_max_iter") != 20_000
        or solver.get("intermediate_diagnostic_max_iter") != 100_000
        or solver.get("max_iter") != 500_000
        or diagnostic.get("status") != "passed_timing_predictions_unscored"
        or diagnostic.get("outer_test_outcomes_read") is not False
        or diagnostic.get("outer_test_metrics_calculated") is not False
    ):
        raise ProductionBundleError("corrected convergence diagnostic differs")
    try:
        arguments.python.resolve(strict=True).relative_to(arguments.runtime.resolve(strict=True))
    except (OSError, ValueError) as error:
        raise ProductionBundleError("runtime Python executable is outside the frozen runtime") from error
    if not arguments.python.resolve(strict=True).is_file():
        raise ProductionBundleError("runtime Python executable is not a regular file")
    if arguments.fold_workers != 5 or arguments.blas_threads != 3:
        raise ProductionBundleError("production resources must remain five workers by three threads")


def _source_manifest(paths: Sequence[Path]) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for path in paths:
        if not path.is_file() or path.is_symlink():
            raise ProductionBundleError(f"campaign source is missing or symlinked: {path}")
        records.append({"path": path.as_posix(), "sha256": sha256_file(path)})
    return records


def _campaign_spec(arguments: argparse.Namespace) -> dict[str, Any]:
    return {
        "schema_version": "masld-bench-gse267145-production-campaign-spec-v1",
        "campaign_id": arguments.campaign.name,
        "logical_units": len(OUTER_FOLDS) * len(SEEDS),
        "outer_folds": list(OUTER_FOLDS),
        "model_seeds": list(SEEDS),
        "model_family_fits_per_unit": MODEL_COUNT,
        "fold_workers": arguments.fold_workers,
        "blas_threads_per_worker": arguments.blas_threads,
        "source_files": _source_manifest(arguments.source),
        "inputs": {
            "runtime": arguments.runtime.as_posix(),
            "runtime_artifacts_sha256": RUNTIME_ARTIFACTS_SHA256,
            "molecular": arguments.molecular.as_posix(),
            "molecular_artifacts_sha256": MOLECULAR_ARTIFACTS_SHA256,
            "folds": arguments.folds.as_posix(),
            "folds_artifacts_sha256": FOLDS_ARTIFACTS_SHA256,
            "fit_views": arguments.fit_views.as_posix(),
            "fit_views_artifacts_sha256": FIT_VIEWS_ARTIFACTS_SHA256,
            "validated_preflight": arguments.preflight.as_posix(),
            "validated_preflight_artifacts_sha256": PREFLIGHT_ARTIFACTS_SHA256,
            "corrected_convergence_diagnostic": arguments.diagnostic.as_posix(),
            "corrected_convergence_diagnostic_artifacts_sha256": DIAGNOSTIC_ARTIFACTS_SHA256,
            "surface": arguments.surface.as_posix(),
            "surface_sha256": SURFACE_SHA256,
            "fitter": arguments.fitter.as_posix(),
            "fitter_sha256": FITTER_SHA256,
        },
        "outcomes_read": False,
        "metrics_calculated": False,
        "scorer_called": False,
        "seeds_are_biological_replicates": False,
    }


def _prepare_spec(arguments: argparse.Namespace, job_id: str) -> str:
    target = arguments.campaign / "campaign_spec"
    expected = _campaign_spec(arguments)
    if target.exists():
        verify_frozen_tree(target)
        observed = _read_json_object(target / "campaign_spec.json")
        if observed != expected:
            raise ProductionBundleError("frozen campaign specification differs")
        return sha256_file(target / "ARTIFACTS.json")
    stage = arguments.campaign / f"campaign_spec.staging-{job_id}"
    if stage.exists():
        raise ProductionBundleError(f"refusing to reuse campaign-spec staging path: {stage}")
    stage.mkdir()
    write_json_exclusive(stage / "campaign_spec.json", expected)
    freeze_tree(
        stage,
        {
            "artifact_class": "gse267145_histology_production_campaign_spec",
            "campaign_id": arguments.campaign.name,
            "status": "locked_before_fit",
        },
    )
    publish_directory_noreplace(stage, target)
    return sha256_file(target / "ARTIFACTS.json")


def _verify_unit(unit: Path, outer_fold: int, seed: int) -> None:
    manifest = verify_frozen_tree(unit)
    metadata = manifest["metadata"]
    expected = {
        "artifact_class": "gse267145_histology_production_outer_seed_unit",
        "outer_fold": outer_fold,
        "seed": seed,
        "outcomes_read": False,
        "metrics_calculated": False,
        "status": "passed_unscored",
    }
    for field, value in expected.items():
        if metadata.get(field) != value:
            raise ProductionBundleError(f"frozen production unit differs: {unit}: {field}")
    receipt = _read_json_object(unit / "unit_receipt.json")
    if (
        receipt.get("outer_fold") != outer_fold
        or receipt.get("seed") != seed
        or receipt.get("production_unit_complete") is not True
        or receipt.get("prediction_files") != MODEL_COUNT
        or receipt.get("outer_test_outcomes_read") is not False
        or receipt.get("outer_test_metrics_calculated") is not False
        or receipt.get("elastic_logistic_solver", {}).get("max_iter") != 500_000
    ):
        raise ProductionBundleError(f"frozen production unit receipt differs: {unit}")
    verify_frozen_tree(unit / "fit/predictions")


def _preserve_failed_unit(
    stage: Path, failed: Path, *, outer_fold: int, seed: int, returncode: int
) -> None:
    if not stage.exists():
        return
    try:
        if (stage / "COMPLETE").is_file():
            verify_frozen_tree(stage)
            publish_directory_noreplace(stage, failed)
            return
        write_text_exclusive(stage / "EXIT_STATUS", f"{returncode}\n")
        freeze_tree(
            stage,
            {
                "artifact_class": "gse267145_histology_production_failed_attempt",
                "outer_fold": outer_fold,
                "seed": seed,
                "returncode": returncode,
                "status": "failed_preserved",
            },
        )
        publish_directory_noreplace(stage, failed)
    except Exception as error:
        raise ProductionBundleError(
            f"fit failed and failure preservation also failed: {stage}: {error}"
        ) from error


def _run_unit(
    arguments: argparse.Namespace,
    *,
    outer_fold: int,
    seed: int,
    job_id: str,
    campaign_spec_sha256: str,
) -> str:
    fold_root = arguments.campaign / "units" / f"outer_{outer_fold}"
    fold_root.mkdir(parents=True, exist_ok=True)
    target = fold_root / f"seed_{seed}"
    if target.exists():
        _verify_unit(target, outer_fold, seed)
        return "reused"
    stage = fold_root / f"seed_{seed}.staging-{job_id}"
    failed = fold_root / f"seed_{seed}.failed-{job_id}"
    if stage.exists() or failed.exists():
        raise ProductionBundleError(f"refusing to reuse unit attempt path: {stage}")
    stage.mkdir()
    fit_output = stage / "fit"
    command = [
        arguments.python.as_posix(),
        arguments.fitter.as_posix(),
        "--molecular",
        arguments.molecular.as_posix(),
        "--molecular-artifacts-sha256",
        MOLECULAR_ARTIFACTS_SHA256,
        "--folds",
        arguments.folds.as_posix(),
        "--folds-artifacts-sha256",
        FOLDS_ARTIFACTS_SHA256,
        "--fit-views",
        arguments.fit_views.as_posix(),
        "--fit-views-artifacts-sha256",
        FIT_VIEWS_ARTIFACTS_SHA256,
        "--surface",
        arguments.surface.as_posix(),
        "--surface-sha256",
        SURFACE_SHA256,
        "--preflight-outer-fold",
        str(outer_fold),
        "--preflight-seed",
        str(seed),
        "--output",
        fit_output.as_posix(),
    ]
    timed_command = [
        "/usr/bin/time",
        "-f",
        '{"elapsed_seconds": %e, "maximum_resident_kib": %M, "user_seconds": %U, "system_seconds": %S}',
        "-o",
        (stage / "fit_resource.json").as_posix(),
        *command,
    ]
    environment = dict(os.environ)
    thread_count = str(arguments.blas_threads)
    environment.update(
        {
            "OMP_NUM_THREADS": thread_count,
            "MKL_NUM_THREADS": thread_count,
            "OPENBLAS_NUM_THREADS": thread_count,
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONNOUSERSITE": "1",
            "PYTHONPATH": (arguments.root / "src").as_posix(),
        }
    )
    with (stage / "fit_stdout.json").open("x", encoding="utf-8") as stdout:
        completed = subprocess.run(
            timed_command,
            cwd=arguments.root,
            env=environment,
            stdout=stdout,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        _preserve_failed_unit(
            stage,
            failed,
            outer_fold=outer_fold,
            seed=seed,
            returncode=completed.returncode,
        )
        raise ProductionBundleError(
            f"outer-fold/seed fit failed: {outer_fold}/{seed}: exit {completed.returncode}"
        )

    fit_receipt = _read_json_object(fit_output / "receipt.json")
    if (
        fit_receipt.get("outer_fold") != outer_fold
        or fit_receipt.get("seed") != seed
        or fit_receipt.get("exact_production_grids_used") is not True
        or fit_receipt.get("model_family_fits") != MODEL_COUNT
        or fit_receipt.get("prediction_files") != MODEL_COUNT
        or fit_receipt.get("outer_test_outcomes_read") is not False
        or fit_receipt.get("outer_test_metrics_calculated") is not False
    ):
        raise ProductionBundleError("one-unit fitter receipt violates production contract")
    predictions = fit_output / "predictions"
    observed_predictions = sorted(path.name for path in predictions.glob("*.tsv"))
    if len(observed_predictions) != MODEL_COUNT:
        raise ProductionBundleError("one-unit fitter did not emit eleven prediction files")
    prediction_artifacts_sha256 = freeze_tree(
        predictions,
        {
            "artifact_class": "gse267145_histology_production_unit_predictions",
            "outer_fold": outer_fold,
            "seed": seed,
            "prediction_files": MODEL_COUNT,
            "outcomes_read": False,
            "metrics_calculated": False,
            "status": "passed_unscored",
        },
    )
    unit_receipt = {
        "schema_version": "masld-bench-gse267145-production-unit-v1",
        "status": "passed_unscored",
        "outer_fold": outer_fold,
        "seed": seed,
        "logical_outer_seed_tasks": 1,
        "model_family_fits": MODEL_COUNT,
        "prediction_files": MODEL_COUNT,
        "prediction_artifacts_sha256": prediction_artifacts_sha256,
        "campaign_spec_artifacts_sha256": campaign_spec_sha256,
        "validated_preflight_artifacts_sha256": PREFLIGHT_ARTIFACTS_SHA256,
        "corrected_convergence_diagnostic_artifacts_sha256": DIAGNOSTIC_ARTIFACTS_SHA256,
        "elastic_logistic_solver": fit_receipt.get("elastic_logistic_solver"),
        "exact_production_grids_used": True,
        "production_unit_complete": True,
        "outer_test_outcomes_read": False,
        "outer_test_metrics_calculated": False,
        "scorer_called": False,
        "seeds_are_biological_replicates": False,
        "fit_counts": fit_receipt.get("fit_counts"),
        "invalid_candidate_count": fit_receipt.get("invalid_candidate_count"),
        "invalid_candidate_reasons": fit_receipt.get("invalid_candidate_reasons"),
    }
    write_json_exclusive(stage / "unit_receipt.json", unit_receipt)
    freeze_tree(
        stage,
        {
            "artifact_class": "gse267145_histology_production_outer_seed_unit",
            "outer_fold": outer_fold,
            "seed": seed,
            "campaign_spec_artifacts_sha256": campaign_spec_sha256,
            "outcomes_read": False,
            "metrics_calculated": False,
            "status": "passed_unscored",
        },
    )
    publish_directory_noreplace(stage, target)
    _verify_unit(target, outer_fold, seed)
    return "created"


def _run_fold(
    arguments: argparse.Namespace,
    *,
    outer_fold: int,
    job_id: str,
    campaign_spec_sha256: str,
) -> dict[str, int]:
    statuses = {"created": 0, "reused": 0}
    for seed in SEEDS:
        try:
            status = _run_unit(
                arguments,
                outer_fold=outer_fold,
                seed=seed,
                job_id=job_id,
                campaign_spec_sha256=campaign_spec_sha256,
            )
        except Exception:
            fold_root = arguments.campaign / "units" / f"outer_{outer_fold}"
            _preserve_failed_unit(
                fold_root / f"seed_{seed}.staging-{job_id}",
                fold_root / f"seed_{seed}.failed-{job_id}",
                outer_fold=outer_fold,
                seed=seed,
                returncode=70,
            )
            raise
        statuses[status] += 1
    return statuses


def _aggregate(arguments: argparse.Namespace, job_id: str) -> str:
    target = arguments.campaign / "aggregate"
    if target.exists():
        verify_frozen_tree(target)
        return "reused"
    stage = arguments.campaign / f"aggregate.staging-{job_id}"
    failed = arguments.campaign / f"aggregate.failed-{job_id}"
    if stage.exists() or failed.exists():
        raise ProductionBundleError(f"refusing to reuse aggregation attempt path: {stage}")
    command = [
        arguments.python.as_posix(),
        arguments.aggregator.as_posix(),
        "--units",
        (arguments.campaign / "units").as_posix(),
        "--folds",
        arguments.folds.as_posix(),
        "--folds-artifacts-sha256",
        FOLDS_ARTIFACTS_SHA256,
        "--output",
        stage.as_posix(),
    ]
    environment = dict(os.environ)
    environment.update(
        {
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONNOUSERSITE": "1",
            "PYTHONPATH": (arguments.root / "src").as_posix(),
        }
    )
    completed = subprocess.run(
        command,
        cwd=arguments.root,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        if stage.exists():
            write_text_exclusive(stage / "aggregate_stdout.txt", completed.stdout)
            write_text_exclusive(stage / "EXIT_STATUS", f"{completed.returncode}\n")
            freeze_tree(
                stage,
                {
                    "artifact_class": "gse267145_histology_production_aggregation_failure",
                    "returncode": completed.returncode,
                    "status": "failed_preserved",
                },
            )
            publish_directory_noreplace(stage, failed)
        raise ProductionBundleError(
            f"outcome-blind aggregation failed with exit {completed.returncode}: {completed.stdout}"
        )
    write_text_exclusive(stage / "aggregate_stdout.txt", completed.stdout)
    freeze_tree(
        stage,
        {
            "artifact_class": "gse267145_histology_production_aggregation",
            "logical_outer_seed_units": 25,
            "prediction_files": 66,
            "outcomes_read": False,
            "metrics_calculated": False,
            "status": "passed_unscored",
        },
    )
    publish_directory_noreplace(stage, target)
    verify_frozen_tree(target)
    return "created"


def _preserve_failed_aggregation(stage: Path, failed: Path) -> None:
    if not stage.exists():
        return
    if (stage / "COMPLETE").is_file():
        verify_frozen_tree(stage)
        publish_directory_noreplace(stage, failed)
        return
    write_text_exclusive(stage / "AGGREGATION_EXIT_STATUS", "70\n")
    freeze_tree(
        stage,
        {
            "artifact_class": "gse267145_histology_production_aggregation_failure",
            "returncode": 70,
            "status": "failed_preserved",
        },
    )
    publish_directory_noreplace(stage, failed)


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    _input_contract(arguments)
    job_id = os.environ.get("SLURM_JOB_ID")
    if not job_id or not job_id.isdigit():
        raise ProductionBundleError("production bundle must run inside a SLURM job")
    arguments.campaign.parent.mkdir(parents=True, exist_ok=True)
    lock_path = arguments.campaign.parent / f".{arguments.campaign.name}.lock"
    with lock_path.open("a+", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ProductionBundleError("another production bundle process holds the campaign lock") from error
        if arguments.campaign.exists() and (arguments.campaign / "COMPLETE").exists():
            manifest = verify_frozen_tree(arguments.campaign)
            return {
                "status": "reused_complete_campaign",
                "campaign_artifacts_sha256": sha256_file(arguments.campaign / "ARTIFACTS.json"),
                "metadata": manifest["metadata"],
            }
        arguments.campaign.mkdir(exist_ok=True)
        campaign_spec_sha256 = _prepare_spec(arguments, job_id)
        created = 0
        reused = 0
        errors: list[str] = []
        with ThreadPoolExecutor(max_workers=arguments.fold_workers) as executor:
            futures = {
                executor.submit(
                    _run_fold,
                    arguments,
                    outer_fold=outer_fold,
                    job_id=job_id,
                    campaign_spec_sha256=campaign_spec_sha256,
                ): outer_fold
                for outer_fold in OUTER_FOLDS
            }
            for future in as_completed(futures):
                outer_fold = futures[future]
                try:
                    status = future.result()
                    created += status["created"]
                    reused += status["reused"]
                except Exception as error:
                    errors.append(f"outer_{outer_fold}: {type(error).__name__}: {error}")
        if errors:
            raise ProductionBundleError("one or more fold workers failed: " + " | ".join(errors))
        for outer_fold in OUTER_FOLDS:
            for seed in SEEDS:
                _verify_unit(
                    arguments.campaign / "units" / f"outer_{outer_fold}" / f"seed_{seed}",
                    outer_fold,
                    seed,
                )
        try:
            aggregation_status = _aggregate(arguments, job_id)
        except Exception:
            _preserve_failed_aggregation(
                arguments.campaign / f"aggregate.staging-{job_id}",
                arguments.campaign / f"aggregate.failed-{job_id}",
            )
            raise
        aggregate_sha256 = sha256_file(arguments.campaign / "aggregate/ARTIFACTS.json")
        receipt = {
            "schema_version": "masld-bench-gse267145-production-bundle-v1",
            "status": "passed_predictions_unscored",
            "campaign_id": arguments.campaign.name,
            "logical_outer_seed_units": 25,
            "model_family_unit_fits": 275,
            "prediction_files": 66,
            "campaign_spec_artifacts_sha256": campaign_spec_sha256,
            "aggregate_artifacts_sha256": aggregate_sha256,
            "outcomes_read": False,
            "metrics_calculated": False,
            "scorer_called": False,
            "seeds_are_biological_replicates": False,
        }
        receipt_path = arguments.campaign / "campaign_receipt.json"
        if receipt_path.exists():
            if _read_json_object(receipt_path) != receipt:
                raise ProductionBundleError("existing campaign receipt differs")
        else:
            write_json_exclusive(receipt_path, receipt)
        campaign_artifacts_sha256 = freeze_tree(
            arguments.campaign,
            {
                "artifact_class": "gse267145_histology_production_prediction_bundle",
                "logical_outer_seed_units": 25,
                "model_family_unit_fits": 275,
                "prediction_files": 66,
                "outcomes_read": False,
                "metrics_calculated": False,
                "status": "passed_unscored",
            },
        )
        receipt["campaign_artifacts_sha256"] = campaign_artifacts_sha256
        receipt["attempt_units_created"] = created
        receipt["attempt_units_reused"] = reused
        receipt["attempt_aggregation_status"] = aggregation_status
        return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--runtime", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--campaign", required=True, type=Path)
    parser.add_argument("--molecular", required=True, type=Path)
    parser.add_argument("--folds", required=True, type=Path)
    parser.add_argument("--fit-views", required=True, type=Path)
    parser.add_argument("--preflight", required=True, type=Path)
    parser.add_argument("--diagnostic", required=True, type=Path)
    parser.add_argument("--surface", required=True, type=Path)
    parser.add_argument("--fitter", required=True, type=Path)
    parser.add_argument("--aggregator", required=True, type=Path)
    parser.add_argument("--source", action="append", type=Path, default=[])
    parser.add_argument("--fold-workers", required=True, type=int)
    parser.add_argument("--blas-threads", required=True, type=int)
    arguments = parser.parse_args()
    if not arguments.source:
        raise ProductionBundleError("at least one source file must be frozen in the campaign spec")
    result = run(arguments)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
