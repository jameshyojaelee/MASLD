"""Reference-retention and query-control metrics allowed before selection."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .embedding import load_embedding, matched_rows
from .execution import require_execution_ownership, verify_execution_record
from .firewall import assert_control_only_metric_names
from .metrics import donor_centroids, macro_f1, normalized_shift_control
from .sampling import capped_group_indices


def _load_training_identity(
    config: dict[str, Any], embedding_manifest: str | Path, expected_schema: str,
) -> tuple[dict[str, Any], dict[str, Any], Path]:
    path = Path(embedding_manifest).resolve()
    manifest_name = (
        "reference_manifest.json"
        if expected_schema == "masld-cl-reference-v1" else "update_manifest.json"
    )
    training_path = path.parent / manifest_name
    if not training_path.is_file():
        raise ContractError(f"embedding lacks sibling {manifest_name}")
    with path.open() as handle:
        embedding = json.load(handle)
    with training_path.open() as handle:
        training = json.load(handle)
    if training.get("schema_version") != expected_schema:
        raise ContractError(f"unexpected training manifest schema: {training_path}")
    if training.get("config_sha256") != config["_config_sha256"]:
        raise ContractError("training manifest config differs from control evaluation")
    if training.get("embedding") != embedding:
        raise ContractError("training manifest does not own the supplied embedding")
    return training, embedding, training_path


def _stratified_bootstrap_groups(values: np.ndarray, strata: np.ndarray, rng) -> np.ndarray:
    selected = []
    for stratum in sorted(set(strata.astype(str))):
        indices = np.flatnonzero(strata.astype(str) == stratum)
        selected.extend(rng.choice(indices, len(indices), replace=True))
    return np.asarray(selected, dtype=np.int64)


def _donor_f1(cells) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    donors, values, studies = [], [], []
    for donor, group in cells.groupby("donor_id", sort=True):
        donors.append(str(donor))
        studies.append(str(group["dataset"].iloc[0]))
        values.append(macro_f1(group["audit_cell_type"], group["predicted_cell_type"]))
    return np.asarray(donors, object), np.asarray(values, float), np.asarray(studies, object)


def _reference_f1_change_ci(reference_cells, candidate_cells, replicates: int, seed: int):
    left = reference_cells.set_index("cell_id", drop=False)
    right = candidate_cells.set_index("cell_id", drop=False)
    shared = sorted(set(left.index) & set(right.index))
    if len(shared) != len(left) or len(shared) != len(right):
        raise ContractError("reference F1 comparison requires identical reference cell rosters")
    left = left.loc[shared].reset_index(drop=True)
    right = right.loc[shared].reset_index(drop=True)
    if not np.array_equal(left["audit_cell_type"].to_numpy(), right["audit_cell_type"].to_numpy()):
        raise ContractError("frozen reference audit labels changed between embeddings")
    left_donor, left_f1, studies = _donor_f1(left)
    right_donor, right_f1, right_studies = _donor_f1(right)
    if not np.array_equal(left_donor, right_donor) or not np.array_equal(studies, right_studies):
        raise ContractError("reference donor roster changed between embeddings")
    rng = np.random.default_rng(seed)
    draws = np.empty(replicates, dtype=float)
    for replicate in range(replicates):
        indices = _stratified_bootstrap_groups(left_donor, studies, rng)
        draws[replicate] = np.mean(right_f1[indices] - left_f1[indices])
    return {
        "estimate": float(np.mean(right_f1 - left_f1)),
        "ci_low": float(np.quantile(draws, 0.025)),
        "ci_high": float(np.quantile(draws, 0.975)),
    }


def _neighbor_jaccard_loss(
    reference_latent: np.ndarray, candidate_latent: np.ndarray, cells,
    *, k: int, maximum_cells: int, seed: int,
) -> float:
    groups = list(zip(cells["donor_id"].astype(str), cells["audit_cell_type"].astype(str)))
    n_groups = len(set(groups))
    cap = max(1, int(np.ceil(maximum_cells / n_groups)))
    indices = capped_group_indices(groups, cap, seed)
    if len(indices) > maximum_cells:
        indices = np.random.default_rng(seed + 1).choice(
            indices, maximum_cells, replace=False
        )
    if len(indices) <= k:
        raise ContractError("reference neighborhood sample does not contain more cells than k")
    before = np.asarray(reference_latent[indices], dtype=np.float32)
    after = np.asarray(candidate_latent[indices], dtype=np.float32)
    before_n = cKDTree(before).query(before, k=k + 1, workers=-1)[1][:, 1:]
    after_n = cKDTree(after).query(after, k=k + 1, workers=-1)[1][:, 1:]
    values = []
    for left, right in zip(before_n, after_n):
        intersection = len(set(map(int, left)) & set(map(int, right)))
        values.append(intersection / (2 * k - intersection))
    return 1.0 - float(np.mean(values))


def _shift_and_standard_error(cells, latent, replicates: int, seed: int):
    selected = cells["strict_reference"] | (cells["primary_query"] & cells["query_control"])
    cells = cells.loc[selected].reset_index(drop=True)
    latent = np.asarray(latent[selected.to_numpy()], dtype=np.float64)
    donor_ids, centroids = donor_centroids(latent, cells["donor_id"].astype(str))
    donor_meta = cells.drop_duplicates("donor_id").set_index("donor_id").loc[donor_ids]
    is_reference = donor_meta["strict_reference"].to_numpy(dtype=bool)
    is_control = (
        donor_meta["primary_query"].to_numpy(dtype=bool)
        & donor_meta["query_control"].to_numpy(dtype=bool)
    )
    if is_reference.sum() < 3 or is_control.sum() < 3:
        raise ContractError("pooled shift-control requires at least three donors per group")
    estimate = normalized_shift_control(centroids[is_reference], centroids[is_control])
    rng = np.random.default_rng(seed)
    values = []
    reference_studies = donor_meta.loc[is_reference, "dataset"].astype(str).to_numpy()
    control_studies = donor_meta.loc[is_control, "dataset"].astype(str).to_numpy()
    reference_values = centroids[is_reference]
    control_values = centroids[is_control]
    for _ in range(replicates):
        r = _stratified_bootstrap_groups(
            np.arange(len(reference_values)), reference_studies, rng
        )
        c = _stratified_bootstrap_groups(
            np.arange(len(control_values)), control_studies, rng
        )
        try:
            values.append(normalized_shift_control(reference_values[r], control_values[c]))
        except ValueError:
            continue
    if len(values) < max(100, int(0.9 * replicates)):
        raise ContractError("too many degenerate shift-control bootstrap replicates")
    return estimate, float(np.std(values, ddof=1))


def evaluate_control_metrics(
    config: dict[str, Any], reference_manifest: str | Path,
    candidate_manifest: str | Path, setting_id: str, ewc_lambda: float,
    replay_fraction: float, reference_execution_record: str | Path,
    candidate_execution_record: str | Path, output: str | Path,
) -> list[dict[str, Any]]:
    reference_run, _, reference_run_path = _load_training_identity(
        config, reference_manifest, "masld-cl-reference-v1"
    )
    candidate_run, _, candidate_run_path = _load_training_identity(
        config, candidate_manifest, "masld-cl-update-v1"
    )
    pipeline_root = Path(config["_config_path"]).resolve().parent
    reference_record_path = Path(reference_execution_record).resolve()
    candidate_record_path = Path(candidate_execution_record).resolve()
    if reference_record_path == candidate_record_path:
        raise ContractError("reference and candidate require distinct GPU records")
    reference_record = verify_execution_record(
        reference_record_path, pipeline_root, config["_config_sha256"]
    )
    candidate_record = verify_execution_record(
        candidate_record_path, pipeline_root, config["_config_sha256"]
    )
    require_execution_ownership(
        reference_record, [reference_manifest, reference_run_path], role="reference"
    )
    require_execution_ownership(
        candidate_record, [candidate_manifest, candidate_run_path], role="candidate"
    )
    expected_setting = f"lambda_{ewc_lambda:g}__replay_{replay_fraction:g}"
    if setting_id != expected_setting:
        raise ContractError("setting ID does not encode the supplied hyperparameters")
    expected_method = (
        "fine_tune" if ewc_lambda == 0 and replay_fraction == 0
        else "replay_only" if ewc_lambda == 0
        else "ewc_only" if replay_fraction == 0
        else "continual_learning"
    )
    if (
        candidate_run.get("production") is not False
        or candidate_run.get("seed") != config["screen"]["seed"]
        or float(candidate_run.get("ewc_lambda")) != ewc_lambda
        or float(candidate_run.get("replay_fraction")) != replay_fraction
        or candidate_run.get("method") != expected_method
        or set(candidate_run.get("query_datasets", []))
        != set(config["roles"]["primary_query_datasets"])
        or candidate_run.get("held_out_datasets", [])
    ):
        raise ContractError("candidate training manifest does not match the control-only setting")
    reference_info, reference_latent, reference_cells = load_embedding(reference_manifest)
    candidate_info, candidate_latent, candidate_cells = load_embedding(candidate_manifest)
    if reference_info["model_kind"] != candidate_info["model_kind"]:
        raise ContractError("reference and candidate model kinds differ")
    reference_cells = reference_cells.loc[reference_cells["strict_reference"]].reset_index(drop=True)
    candidate_reference = candidate_cells.loc[candidate_cells["strict_reference"]].reset_index(drop=True)
    left, right = matched_rows(reference_cells, candidate_reference)
    reference_cells = reference_cells.iloc[left].reset_index(drop=True)
    candidate_reference = candidate_reference.iloc[right].reset_index(drop=True)
    f1 = _reference_f1_change_ci(
        reference_cells, candidate_reference,
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"],
    )
    jaccard_loss = _neighbor_jaccard_loss(
        np.asarray(reference_latent)[left], np.asarray(candidate_latent)[
            np.flatnonzero(candidate_cells["strict_reference"].to_numpy())[right]
        ],
        reference_cells,
        k=config["evaluation"]["reference_neighborhood_k"],
        maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
        seed=config["screen"]["seed"],
    )
    shift, shift_se = _shift_and_standard_error(
        candidate_cells, candidate_latent,
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 1,
    )
    values = {
        "reference_macro_f1_change_ci_low": f1["ci_low"],
        "reference_neighborhood_jaccard_loss": jaccard_loss,
        "shift_control": shift,
        "shift_control_standard_error": shift_se,
    }
    assert_control_only_metric_names(list(values))
    rows = [{
        "setting_id": setting_id,
        "ewc_lambda": ewc_lambda,
        "replay_fraction": replay_fraction,
        "model_kind": candidate_info["model_kind"],
        "metric": metric,
        "scope": "pooled_primary_query",
        "value": value,
    } for metric, value in values.items()]
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    metric_sha = sha256_path(output)
    write_json_exclusive(output.with_suffix(output.suffix + ".details.json"), {
        "schema_version": "masld-cl-control-evaluation-v1",
        "setting_id": setting_id,
        "model_kind": candidate_info["model_kind"],
        "reference_f1_change": f1,
        "reference_embedding": str(Path(reference_manifest).resolve()),
        "candidate_embedding": str(Path(candidate_manifest).resolve()),
        "metrics_realpath": str(output.resolve()),
        "metrics_sha256": metric_sha,
        "reference_run_manifest": str(reference_run_path.resolve()),
        "reference_run_manifest_sha256": sha256_path(reference_run_path),
        "candidate_run_manifest": str(candidate_run_path.resolve()),
        "candidate_run_manifest_sha256": sha256_path(candidate_run_path),
        "reference_execution_record": str(reference_record_path),
        "reference_execution_record_sha256": sha256_path(reference_record_path),
        "candidate_execution_record": str(candidate_record_path),
        "candidate_execution_record_sha256": sha256_path(candidate_record_path),
        "candidate_seed": candidate_run["seed"],
        "candidate_method": candidate_run["method"],
        "config_sha256": config["_config_sha256"],
    })
    return rows
