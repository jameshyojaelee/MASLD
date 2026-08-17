"""Derive acquisition-order gates from six locked sequence outputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .benchmark import _bundle_centroids
from .config import write_json_exclusive
from .confirmation_evaluation import (
    _derive_passes, _required_run_manifests, _verify_evaluation,
)
from .contracts import ContractError, sha256_path
from .embedding import load_embedding
from .execution import verify_execution_record
from .firewall import load_selection_lock
from .metrics import donor_centroids, donor_distance_spearman, normalized_shift_control


def _execution_files(
    config: dict[str, Any], values: list[str], cache: dict[str, dict[str, Any]],
) -> tuple[set[str], list[dict[str, Any]]]:
    pipeline_root = Path(config["_config_path"]).resolve().parent
    files: set[str] = set()
    records = []
    for value in values:
        path = str(Path(value).resolve())
        if path not in cache:
            cache[path] = verify_execution_record(
                path, pipeline_root, config["_config_sha256"]
            )
        record = cache[path]
        identity = record["result_identity"]
        root = Path(identity["output_realpath"])
        files.add(identity["result_manifest_realpath"])
        files.update(str((root / item["path"]).resolve()) for item in identity["output_files"])
        records.append(record)
    return files, records


def _metric_record_paths(*details_values: dict[str, Any]) -> set[str]:
    result = set()
    for details in details_values:
        for source in details.get("execution_records", []):
            if not isinstance(source, dict) or set(source) != {"path", "sha256"}:
                raise ContractError("order metric provenance has malformed execution records")
            path = str(Path(source["path"]).resolve())
            if sha256_path(path) != source["sha256"]:
                raise ContractError("order metric execution record changed")
            result.add(path)
    return result


def _candidate_donor_representation(details: dict[str, Any]):
    path = details["embedding_sources"]["candidate|all_lineage"]["embedding_manifest"]
    _, latent, cells = load_embedding(path)
    mask = cells["analysis_eligible"].to_numpy(dtype=bool)
    selected = cells.loc[mask].reset_index(drop=True)
    donors, values = donor_centroids(
        np.asarray(latent[mask], dtype=float), selected["donor_id"].astype(str)
    )
    return donors.astype(str), values


def _matched_distance_spearman(left, right) -> float:
    left_donors, left_values = left
    right_donors, right_values = right
    if set(left_donors) != set(right_donors):
        raise ContractError("order and joint donor rosters differ")
    index = {donor: i for i, donor in enumerate(right_donors)}
    order = np.asarray([index[donor] for donor in left_donors], dtype=np.int64)
    return donor_distance_spearman(left_values, right_values[order])


def _control_shift_from_embedding(path: str | Path) -> float:
    bundle = load_embedding(path)
    _, values, _, reference, controls = _bundle_centroids(bundle)
    if reference.sum() < 3 or controls.sum() < 3:
        raise ContractError("de novo seed control shift lacks three donors per group")
    return normalized_shift_control(values[reference], values[controls])


def evaluate_order_bundle(
    config: dict[str, Any], selection_lock: str | Path,
    bundle_path: str | Path, output: str | Path,
) -> list[dict[str, Any]]:
    selection = load_selection_lock(selection_lock, config)
    with Path(bundle_path).open() as handle:
        bundle = json.load(handle)
    if (
        bundle.get("schema_version") != "masld-cl-order-bundle-v2"
        or bundle.get("selection_lock_sha256") != selection["lock_sha256"]
    ):
        raise ContractError("order bundle schema or selection lock is invalid")
    expected_orders = set(config["acquisition_orders"])
    entries = bundle.get("orders", [])
    observed_orders = [entry.get("order_name") for entry in entries]
    if set(observed_orders) != expected_orders or len(observed_orders) != len(expected_orders):
        raise ContractError("order bundle does not contain the six exact acquisition orders")
    cache: dict[str, dict[str, Any]] = {}
    joint = bundle.get("joint", {})
    joint_benchmark, joint_details, _ = _verify_evaluation(
        joint["benchmark_metrics"], "masld-cl-benchmark-v1", selection["lock_sha256"], config
    )
    joint_outcome, joint_outcome_details, _ = _verify_evaluation(
        joint["outcome_metrics"], "masld-cl-outcome-evaluation-v1", selection["lock_sha256"], config
    )
    joint_files, _ = _execution_files(config, joint.get("execution_records", []), cache)
    if {str(Path(value).resolve()) for value in joint.get("execution_records", [])} != _metric_record_paths(
        joint_details, joint_outcome_details
    ):
        raise ContractError("joint order reference supplies different GPU records than its metrics")
    if not _required_run_manifests(joint_details, joint_outcome_details).issubset(joint_files):
        raise ContractError("joint metrics lack their GPU execution records")
    joint_representation = _candidate_donor_representation(joint_details)
    joint_shift = float(
        joint_details["alignment"]["all_lineage"]["harmony"]["candidate_shift"]
    )
    if joint_shift <= 0:
        raise ContractError("joint control shift must be positive")
    order_rows = []
    for entry in entries:
        name = entry["order_name"]
        benchmark, details, _ = _verify_evaluation(
            entry["benchmark_metrics"], "masld-cl-benchmark-v1", selection["lock_sha256"], config
        )
        outcome, outcome_details, _ = _verify_evaluation(
            entry["outcome_metrics"], "masld-cl-outcome-evaluation-v1", selection["lock_sha256"], config
        )
        files, records = _execution_files(config, entry.get("execution_records", []), cache)
        if {str(Path(value).resolve()) for value in entry.get("execution_records", [])} != _metric_record_paths(
            details, outcome_details
        ):
            raise ContractError(f"order supplies different GPU records than its metrics: {name}")
        if not _required_run_manifests(details, outcome_details).issubset(files):
            raise ContractError(f"order metrics lack GPU execution records: {name}")
        sequence_manifests = []
        for record in records:
            result = Path(record["result_identity"]["result_manifest_realpath"])
            with result.open() as handle:
                manifest = json.load(handle)
            if manifest.get("schema_version") == "masld-cl-sequence-v1":
                sequence_manifests.append(manifest)
        if len(sequence_manifests) != len(config["lineages"]) + 1 or any(
            manifest.get("order_name") != name
            or manifest.get("order") != config["acquisition_orders"][name]
            for manifest in sequence_manifests
        ):
            raise ContractError(f"order is not backed by six matching sequence runs: {name}")
        passes = _derive_passes(config, benchmark, outcome)
        shift = float(details["alignment"]["all_lineage"]["harmony"]["candidate_shift"])
        order_rows.append({
            "order_name": name,
            "control_shift": shift,
            "degradation": (shift - joint_shift) / joint_shift,
            "distance_spearman": _matched_distance_spearman(
                _candidate_donor_representation(details), joint_representation
            ),
            **passes,
        })
    de_novo = bundle.get("de_novo_seed_reference", {})
    embedding_paths = [str(Path(value).resolve()) for value in de_novo.get("embeddings", [])]
    if len(embedding_paths) != len(config["screen"]["confirmation_seeds"]):
        raise ContractError("order bundle lacks five de novo seed embeddings")
    de_novo_files, _ = _execution_files(
        config, de_novo.get("execution_records", []), cache
    )
    observed_seeds = set()
    shifts = []
    for embedding_path in embedding_paths:
        run_path = Path(embedding_path).parent / "denovo_manifest.json"
        if str(run_path.resolve()) not in de_novo_files:
            raise ContractError("de novo seed embedding lacks its GPU execution record")
        with run_path.open() as handle:
            run = json.load(handle)
        if (
            run.get("selection_lock_sha256") != selection["lock_sha256"]
            or run.get("model_kind") != "all_lineage"
        ):
            raise ContractError("de novo seed run provenance is invalid")
        observed_seeds.add(run.get("seed"))
        shifts.append(_control_shift_from_embedding(embedding_path))
    if observed_seeds != set(config["screen"]["confirmation_seeds"]):
        raise ContractError("de novo seed embeddings do not use the five exact seeds")
    de_novo_variability = float(np.std(shifts, ddof=1))
    if de_novo_variability <= 0:
        raise ContractError("de novo seed variability must be positive")
    order_variability = float(np.std([row["control_shift"] for row in order_rows], ddof=1))
    metrics = [
        {"metric": "all_orders_reference_and_disease_pass", "scope": "six_orders", "value": int(all(
            row["reference_pass"] and row["disease_pass"] for row in order_rows
        ))},
        {"metric": "passing_order_count", "scope": "six_orders", "value": sum(
            row["alignment_pass"] for row in order_rows
        )},
        {"metric": "worst_order_degradation", "scope": "six_orders", "value": max(
            row["degradation"] for row in order_rows
        )},
        {"metric": "order_distance_spearman", "scope": "worst_order", "value": min(
            row["distance_spearman"] for row in order_rows
        )},
        {"metric": "order_to_seed_variability_ratio", "scope": "six_orders", "value": order_variability / de_novo_variability},
    ]
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "scope", "value"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(metrics)
    write_json_exclusive(output.with_suffix(output.suffix + ".details.json"), {
        "schema_version": "masld-cl-order-evaluation-v2",
        "config_sha256": config["_config_sha256"],
        "selection_lock_sha256": selection["lock_sha256"],
        "bundle_realpath": str(Path(bundle_path).resolve()),
        "bundle_sha256": sha256_path(bundle_path),
        "joint_control_shift": joint_shift,
        "order_variability": order_variability,
        "de_novo_seed_variability": de_novo_variability,
        "orders": order_rows,
        "execution_records": [
            {"path": path, "sha256": sha256_path(path)} for path in sorted(cache)
        ],
        "metrics_realpath": str(output.resolve()),
        "metrics_sha256": sha256_path(output),
    })
    return metrics
