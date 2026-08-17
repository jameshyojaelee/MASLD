"""Reference-only kNN audit of latent cell-identity preservation."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from .embedding import load_embedding
from .firewall import validate_program_firewall
from .metrics import macro_f1, positive_class_f1


class LatentKNNLabelAuditError(ContractError):
    """Raised when the fixed V25 latent-label audit changes identity."""


POLICY_SCHEMA = "masld-cl-latent-knn-label-audit-policy-v25"


def load_latent_knn_label_audit_policy(
    config: dict[str, Any], value: str | Path,
) -> tuple[Path, dict[str, Any]]:
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "latent_knn_label_audit_policy_v25.json"
    )
    if path != expected:
        raise LatentKNNLabelAuditError("V25 requires its source-controlled policy")
    with path.open() as handle:
        policy = json.load(handle)
    classifier = policy.get("classifier", {})
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("status") != "single_reference_only_latent_identity_audit"
        or policy.get("case_stage_program_hero_gene_umap_cas13_used") is not False
        or classifier.get("query_labels_available_to_fit") is not False
        or classifier.get("neighbors") != config["evaluation"]["reference_neighborhood_k"]
        or classifier.get("reference_cap_seed") != config["screen"]["seed"]
        or classifier.get("pynndescent_version")
        != importlib.metadata.version("pynndescent")
    ):
        raise LatentKNNLabelAuditError("V25 policy identity differs")
    root = path.parent
    parent_spec = policy["parent_v14_outcome"]
    parent_path = (root / parent_spec["path"]).resolve()
    with parent_path.open() as handle:
        parent = json.load(handle)
    label = parent.get("label_concordance", {})
    changes = label.get("major_lineage_f1_changes", {})
    worst = min(changes, key=changes.get) if changes else None
    if (
        sha256_path(parent_path) != parent_spec["sha256"]
        or parent.get("pilot_pass") is not parent_spec["required_pilot_pass"]
        or label.get("gates", {}).get(parent_spec["required_failed_gate"]) is not False
        or worst != parent_spec["required_worst_lineage"]
    ):
        raise LatentKNNLabelAuditError("V25 parent V14 outcome differs")
    candidate = policy["candidate"]
    for key, hash_key in (("embedding", "embedding_sha256"), ("selection_lock", "selection_lock_sha256")):
        source = (root / candidate[key]).resolve()
        if sha256_path(source) != candidate[hash_key]:
            raise LatentKNNLabelAuditError(f"V25 candidate {key} changed")
    if set(policy["non_cl_embeddings"]) != {
        "architecture_surgery", "de_novo", "ewc_only", "fine_tune", "replay_only"
    }:
        raise LatentKNNLabelAuditError("V25 non-CL method roster differs")
    for method, source in policy["non_cl_embeddings"].items():
        if sha256_path((root / source["path"]).resolve()) != source["sha256"]:
            raise LatentKNNLabelAuditError(f"V25 non-CL embedding changed: {method}")
    return path, policy


def _majority_vote(neighbor_codes: np.ndarray, classes: np.ndarray) -> np.ndarray:
    neighbor_codes = np.asarray(neighbor_codes)
    if neighbor_codes.ndim != 2 or not len(neighbor_codes):
        raise LatentKNNLabelAuditError("V25 neighbor-code matrix is invalid")
    counts = np.stack(
        [(neighbor_codes == code).sum(axis=1) for code in range(len(classes))],
        axis=1,
    )
    # classes are lexicographically sorted, so argmax supplies the locked tie break.
    return classes[np.argmax(counts, axis=1)]


def _fit_predict_knn(
    reference_latent: np.ndarray, reference_labels: np.ndarray,
    query_latent: np.ndarray, classifier: dict[str, Any],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Fit from reference labels only; query labels are not an accepted argument."""
    from pynndescent import NNDescent

    reference_latent = np.asarray(reference_latent, dtype=np.float32)
    query_latent = np.asarray(query_latent, dtype=np.float32)
    classes = np.asarray(sorted(set(map(str, reference_labels))), dtype=object)
    lookup = {value: index for index, value in enumerate(classes)}
    label_codes = np.asarray([lookup[str(value)] for value in reference_labels], dtype=np.int16)
    index = NNDescent(
        reference_latent,
        n_neighbors=classifier["neighbors"],
        metric=classifier["metric"],
        random_state=classifier["random_state"],
        n_jobs=classifier["n_jobs"],
        low_memory=classifier["low_memory"],
        verbose=False,
    )
    index.prepare()
    neighbors, distances = index.query(
        query_latent, k=classifier["neighbors"],
        epsilon=classifier["query_epsilon"],
    )
    if (
        neighbors.shape != (len(query_latent), classifier["neighbors"])
        or np.any(neighbors < 0)
        or np.any(neighbors >= len(reference_latent))
        or not np.isfinite(distances).all()
    ):
        raise LatentKNNLabelAuditError("V25 neighbor query returned invalid results")
    predictions = _majority_vote(label_codes[neighbors], classes)
    identity = {
        "classes": list(map(str, classes)),
        "neighbor_indices_sha256": hashlib.sha256(
            np.asarray(neighbors, dtype="<i8").tobytes()
        ).hexdigest(),
        "neighbor_distances_sha256_float32_le": hashlib.sha256(
            np.asarray(distances, dtype="<f4").tobytes()
        ).hexdigest(),
        "predictions_sha256": hashlib.sha256(
            "\0".join(map(str, predictions)).encode("utf-8")
        ).hexdigest(),
    }
    return predictions, identity


def _prediction_scores(cells, predictions: np.ndarray, lineages: list[str]):
    truth = cells["audit_cell_type"].astype(str).to_numpy()
    donors = cells["donor_id"].astype(str).to_numpy()
    donor_macro = {}
    donor_lineage = {lineage: {} for lineage in lineages}
    for donor in sorted(set(donors)):
        keep = donors == donor
        donor_macro[donor] = macro_f1(truth[keep], predictions[keep])
        for lineage in lineages:
            positive = truth[keep] == lineage
            if positive.any():
                donor_lineage[lineage][donor] = positive_class_f1(
                    positive, predictions[keep] == lineage
                )
    if any(not values for values in donor_lineage.values()):
        raise LatentKNNLabelAuditError("V25 major lineage is not evaluable")
    return {
        "mean_donor_macro_f1": float(np.mean(list(donor_macro.values()))),
        "mean_donor_lineage_f1": {
            lineage: float(np.mean(list(values.values())))
            for lineage, values in donor_lineage.items()
        },
        "donor_macro_f1": donor_macro,
        "donor_lineage_f1": donor_lineage,
    }


def run_latent_knn_label_audit(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    from .training import _capped_indices

    validate_program_firewall(config)
    policy_path, policy = load_latent_knn_label_audit_policy(config, policy_value)
    root = policy_path.parent
    sources = {
        "candidate": (root / policy["candidate"]["embedding"]).resolve(),
        **{
            method: (root / source["path"]).resolve()
            for method, source in policy["non_cl_embeddings"].items()
        },
    }
    output = Path(output_value)
    output.mkdir(parents=True, exist_ok=False)
    method_scores = {}
    roster_identity = None
    for method, source_path in sources.items():
        info, latent, cells = load_embedding(source_path)
        reference = cells["strict_reference"].to_numpy(dtype=bool)
        query = (
            cells["analysis_eligible"].to_numpy(dtype=bool)
            & ~reference
        )
        if reference.sum() != 216957 or query.sum() != 687559:
            raise LatentKNNLabelAuditError(f"V25 cell roster differs: {method}")
        current_roster = hashlib.sha256(
            "\0".join(cells.loc[query, "cell_id"].astype(str)).encode("utf-8")
        ).hexdigest()
        if roster_identity is None:
            roster_identity = current_roster
        elif current_roster != roster_identity:
            raise LatentKNNLabelAuditError("V25 query rosters differ across methods")
        reference_pool = _capped_indices(
            type("ADataView", (), {"obs": cells, "n_obs": len(cells)})(),
            np.flatnonzero(reference), "all_lineage", config,
            policy["classifier"]["reference_cap_seed"],
        )
        if len(reference_pool) != 61960:
            raise LatentKNNLabelAuditError("V25 capped reference roster differs")
        predictions, knn_identity = _fit_predict_knn(
            np.asarray(latent)[reference_pool],
            cells.iloc[reference_pool]["audit_cell_type"].astype(str).to_numpy(),
            np.asarray(latent)[query], policy["classifier"],
        )
        query_cells = cells.loc[query].reset_index(drop=True)
        scores = _prediction_scores(query_cells, predictions, config["lineages"])
        prediction_path = output / f"{method}_predictions.tsv.gz"
        with DeterministicGzipTextWriter(prediction_path) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["cell_id", "donor_id", "audit_cell_type", "knn_predicted_cell_type"])
            writer.writerows(zip(
                query_cells["cell_id"].astype(str),
                query_cells["donor_id"].astype(str),
                query_cells["audit_cell_type"].astype(str),
                map(str, predictions),
            ))
        method_scores[method] = {
            **scores,
            "embedding": str(source_path),
            "embedding_sha256": sha256_path(source_path),
            "embedding_model_kind": info["model_kind"],
            "n_reference_training_cells": int(len(reference_pool)),
            "n_query_cells": int(query.sum()),
            "knn_identity": knn_identity,
            "predictions_file": str(prediction_path.resolve()),
            "predictions_file_sha256": sha256_path(prediction_path),
        }
        write_json_exclusive(output / f"progress_{method}.json", {
            "schema_version": "masld-cl-latent-knn-label-progress-v25",
            "config_sha256": config["_config_sha256"],
            "method": method,
            "score": method_scores[method],
        })
        del latent, cells, query_cells, predictions

    non_cl = {key: value for key, value in method_scores.items() if key != "candidate"}
    best_method = sorted(
        non_cl, key=lambda method: (-non_cl[method]["mean_donor_macro_f1"], method)
    )[0]
    candidate = method_scores["candidate"]
    baseline = non_cl[best_method]
    macro_change = candidate["mean_donor_macro_f1"] - baseline["mean_donor_macro_f1"]
    lineage_changes = {
        lineage: (
            candidate["mean_donor_lineage_f1"][lineage]
            - baseline["mean_donor_lineage_f1"][lineage]
        )
        for lineage in config["lineages"]
    }
    thresholds = policy["selection"]
    gates = [
        {
            "gate": "query_macro_f1_noninferior",
            "observed": macro_change,
            "threshold": thresholds["candidate_query_macro_f1_change_at_least"],
            "pass": macro_change >= thresholds["candidate_query_macro_f1_change_at_least"],
        },
        {
            "gate": "all_major_lineages_noninferior",
            "observed": min(lineage_changes.values()),
            "threshold": thresholds["every_major_lineage_f1_change_at_least"],
            "pass": min(lineage_changes.values())
            >= thresholds["every_major_lineage_f1_change_at_least"],
        },
    ]
    result = {
        "schema_version": "masld-cl-latent-knn-label-audit-v25",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "reference_only_fit": True,
        "query_labels_available_to_fit": False,
        "query_roster_sha256": roster_identity,
        "best_non_cl_method": best_method,
        "method_scores": method_scores,
        "candidate_query_macro_f1_change": macro_change,
        "candidate_major_lineage_f1_changes": lineage_changes,
        "gates": gates,
        "pass": all(item["pass"] for item in gates),
        "labels_are_audit_targets_and_were_not_overwritten": True,
        "case_stage_program_hero_gene_umap_cas13_read": False,
        "sensitivity_only": True,
    }
    write_json_exclusive(output / "label_audit.json", result)
    return result
