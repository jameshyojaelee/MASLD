"""Fail-closed audits for GPU repeatability and descriptive-only projection."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .embedding import load_embedding
from .execution import require_execution_ownership, verify_execution_record
from .firewall import load_selection_lock


def _difference_stats(left: np.ndarray, right: np.ndarray) -> tuple[float, float]:
    if left.shape != right.shape:
        raise ContractError("determinism comparison requires matched embedding dimensions")
    maximum, total, count = 0.0, 0.0, 0
    for start in range(0, len(left), 50_000):
        difference = np.abs(
            np.asarray(left[start:start + 50_000], dtype=np.float64)
            - np.asarray(right[start:start + 50_000], dtype=np.float64)
        )
        maximum = max(maximum, float(difference.max(initial=0)))
        total += float(difference.sum())
        count += int(difference.size)
    return maximum, total / count if count else 0.0


def _repeatability_source(
    config: dict[str, Any], selection: dict[str, Any], embedding_value: str | Path,
    record_value: str | Path, seed: int,
) -> tuple[dict[str, Any], np.ndarray, Any, dict[str, Any]]:
    embedding_path = Path(embedding_value).resolve()
    info, latent, cells = load_embedding(embedding_path)
    run_path = embedding_path.parent / "update_manifest.json"
    environment_path = embedding_path.parent / "environment.json"
    with run_path.open() as handle:
        run = json.load(handle)
    selected = selection["selected"]
    if (
        info.get("model_kind") != "all_lineage"
        or run.get("schema_version") != "masld-cl-update-v1"
        or run.get("config_sha256") != config["_config_sha256"]
        or run.get("selection_lock_sha256") != selection["lock_sha256"]
        or run.get("embedding") != info
        or run.get("seed") != seed
        or run.get("method") != "continual_learning"
        or run.get("replay_mode") != "random"
        or run.get("sensitivity_only") is True
        or float(run.get("ewc_lambda")) != float(selected["ewc_lambda"])
        or float(run.get("replay_fraction")) != float(selected["replay_fraction"])
    ):
        raise ContractError("GPU repeatability source is not the selected all-lineage CL fit")
    with environment_path.open() as handle:
        environment = json.load(handle)
    record_path = Path(record_value).resolve()
    pipeline_root = Path(config["_config_path"]).resolve().parent
    record = verify_execution_record(record_path, pipeline_root, config["_config_sha256"])
    if not record.get("slurm_gpus") or environment.get("cuda_available") is not True:
        raise ContractError("GPU repeatability source was not produced on a CUDA GPU")
    require_execution_ownership(
        record, [embedding_path, run_path, environment_path], role="GPU repeatability"
    )
    source = {
        "embedding_manifest": str(embedding_path),
        "embedding_manifest_sha256": sha256_path(embedding_path),
        "run_manifest": str(run_path.resolve()),
        "run_manifest_sha256": sha256_path(run_path),
        "environment": str(environment_path.resolve()),
        "environment_sha256": sha256_path(environment_path),
        "execution_record": str(record_path),
        "execution_record_sha256": sha256_path(record_path),
    }
    return source, latent, cells, run


def measure_gpu_tolerance(
    config: dict[str, Any], selection_lock: str | Path,
    left_embedding: str | Path, right_embedding: str | Path,
    left_record: str | Path, right_record: str | Path,
    seed: int, output: str | Path,
) -> dict[str, Any]:
    selection = load_selection_lock(selection_lock, config)
    if Path(left_record).resolve() == Path(right_record).resolve():
        raise ContractError("GPU tolerance requires two independent execution records")
    left, left_latent, left_cells, left_run = _repeatability_source(
        config, selection, left_embedding, left_record, seed
    )
    right, right_latent, right_cells, right_run = _repeatability_source(
        config, selection, right_embedding, right_record, seed
    )
    context_keys = (
        "seed", "production", "query_datasets", "held_out_datasets",
        "ewc_lambda", "replay_fraction", "method", "replay_mode",
    )
    if any(left_run.get(key) != right_run.get(key) for key in context_keys):
        raise ContractError("GPU tolerance fits have different training contexts")
    if not np.array_equal(
        left_cells["cell_id"].astype(str).to_numpy(),
        right_cells["cell_id"].astype(str).to_numpy(),
    ):
        raise ContractError("GPU tolerance fits have different cell rosters or order")
    maximum, mean = _difference_stats(left_latent, right_latent)
    result = {
        "schema_version": "masld-cl-seed-tolerance-v2",
        "config_sha256": config["_config_sha256"],
        "selection_lock_sha256": selection["lock_sha256"],
        "device": "gpu",
        "seed": seed,
        "shape": list(left_latent.shape),
        "maximum_absolute_difference": maximum,
        "mean_absolute_difference": mean,
        "sources": [left, right],
        "gpu_tolerance_is_measured_not_assumed": True,
    }
    write_json_exclusive(output, result)
    return result


def verify_gpu_tolerance(
    config: dict[str, Any], selection: dict[str, Any], path_value: str | Path,
) -> dict[str, Any]:
    path = Path(path_value).resolve()
    with path.open() as handle:
        result = json.load(handle)
    if (
        result.get("schema_version") != "masld-cl-seed-tolerance-v2"
        or result.get("config_sha256") != config["_config_sha256"]
        or result.get("selection_lock_sha256") != selection["lock_sha256"]
        or result.get("device") != "gpu"
        or result.get("gpu_tolerance_is_measured_not_assumed") is not True
        or len(result.get("sources", [])) != 2
    ):
        raise ContractError("GPU determinism tolerance was not measured from locked selected fits")
    loaded = []
    for source in result["sources"]:
        for path_key, hash_key in (
            ("embedding_manifest", "embedding_manifest_sha256"),
            ("run_manifest", "run_manifest_sha256"),
            ("environment", "environment_sha256"),
            ("execution_record", "execution_record_sha256"),
        ):
            if sha256_path(source[path_key]) != source[hash_key]:
                raise ContractError("GPU tolerance source changed after measurement")
        loaded.append(_repeatability_source(
            config, selection, source["embedding_manifest"],
            source["execution_record"], int(result["seed"]),
        ))
    if not np.array_equal(
        loaded[0][2]["cell_id"].astype(str).to_numpy(),
        loaded[1][2]["cell_id"].astype(str).to_numpy(),
    ):
        raise ContractError("GPU tolerance cell roster changed")
    maximum, mean = _difference_stats(loaded[0][1], loaded[1][1])
    if maximum != float(result.get("maximum_absolute_difference")) or not np.isclose(
        mean, float(result.get("mean_absolute_difference")), rtol=0, atol=1e-15
    ):
        raise ContractError("GPU tolerance values do not reproduce")
    return result


def audit_descriptive_projection(
    config: dict[str, Any], selection_lock: str | Path,
    update_manifest: str | Path, embedding_manifest: str | Path,
    execution_record: str | Path, output: str | Path,
) -> dict[str, Any]:
    selection = load_selection_lock(selection_lock, config)
    update_path = Path(update_manifest).resolve()
    embedding_path = Path(embedding_manifest).resolve()
    with update_path.open() as handle:
        update = json.load(handle)
    embedding, _, _ = load_embedding(embedding_path)
    selected = selection["selected"]
    if (
        update.get("schema_version") != "masld-cl-update-v1"
        or update.get("config_sha256") != config["_config_sha256"]
        or update.get("selection_lock_sha256") != selection["lock_sha256"]
        or update.get("production") is not True
        or update.get("model_kind") != "all_lineage"
        or update.get("method") != "continual_learning"
        or update.get("replay_mode") != "random"
        or update.get("sensitivity_only") is True
        or float(update.get("ewc_lambda")) != float(selected["ewc_lambda"])
        or float(update.get("replay_fraction")) != float(selected["replay_fraction"])
        or update.get("embedding") != embedding
    ):
        raise ContractError("descriptive projection is not owned by the selected production fit")
    record_path = Path(execution_record).resolve()
    record = verify_execution_record(
        record_path, Path(config["_config_path"]).resolve().parent,
        config["_config_sha256"],
    )
    require_execution_ownership(
        record, [update_path, embedding_path], role="descriptive projection"
    )
    expected = 33
    result = {
        "schema_version": "masld-cl-descriptive-projection-v2",
        "config_sha256": config["_config_sha256"],
        "selection_lock_sha256": selection["lock_sha256"],
        "update_manifest": str(update_path),
        "update_manifest_sha256": sha256_path(update_path),
        "embedding_manifest": str(embedding_path),
        "embedding_manifest_sha256": sha256_path(embedding_path),
        "execution_record": str(record_path),
        "execution_record_sha256": sha256_path(record_path),
        "descriptive_only_cells": expected,
        "trained_descriptive_only_cells": update.get("n_descriptive_only_training_cells"),
        "projected_descriptive_only_cells": embedding.get("n_analysis_ineligible"),
        "passed": (
            update.get("n_descriptive_only_training_cells") == 0
            and embedding.get("n_analysis_ineligible") == expected
        ),
    }
    if not result["passed"]:
        raise ContractError(f"descriptive projection firewall failed: {result}")
    write_json_exclusive(output, result)
    return result


def verify_descriptive_projection(
    config: dict[str, Any], selection: dict[str, Any], path_value: str | Path,
) -> dict[str, Any]:
    path = Path(path_value).resolve()
    with path.open() as handle:
        result = json.load(handle)
    if (
        result.get("schema_version") != "masld-cl-descriptive-projection-v2"
        or result.get("config_sha256") != config["_config_sha256"]
        or result.get("selection_lock_sha256") != selection["lock_sha256"]
        or result.get("passed") is not True
    ):
        raise ContractError("descriptive projection audit is invalid")
    for path_key, hash_key in (
        ("update_manifest", "update_manifest_sha256"),
        ("embedding_manifest", "embedding_manifest_sha256"),
        ("execution_record", "execution_record_sha256"),
    ):
        if sha256_path(result[path_key]) != result[hash_key]:
            raise ContractError("descriptive projection source changed")
    update_path = Path(result["update_manifest"])
    embedding_path = Path(result["embedding_manifest"])
    with update_path.open() as handle:
        update = json.load(handle)
    embedding, _, _ = load_embedding(embedding_path)
    record = verify_execution_record(
        result["execution_record"], Path(config["_config_path"]).resolve().parent,
        config["_config_sha256"],
    )
    require_execution_ownership(record, [update_path, embedding_path], role="descriptive projection")
    if (
        update.get("selection_lock_sha256") != selection["lock_sha256"]
        or update.get("embedding") != embedding
        or update.get("production") is not True
        or update.get("model_kind") != "all_lineage"
        or update.get("method") != "continual_learning"
        or update.get("replay_mode") != "random"
        or update.get("sensitivity_only") is True
        or float(update.get("ewc_lambda")) != float(selection["selected"]["ewc_lambda"])
        or float(update.get("replay_fraction")) != float(selection["selected"]["replay_fraction"])
        or update.get("n_descriptive_only_training_cells") != 0
        or embedding.get("n_analysis_ineligible") != 33
    ):
        raise ContractError("descriptive projection audit no longer reproduces")
    return result
