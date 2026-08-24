"""Post-mapping diagnostic and independent evaluation for Router V2."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .external_gse212837 import _read_frame_column
from .metrics import macro_f1
from .router_v2 import collapse_router_labels, load_router_v2_policy


class RouterV2EvaluationError(ContractError):
    """Raised when a Router V2 evaluation input differs."""


def paired_donor_bootstrap(
    candidate: dict[str, float], baseline: dict[str, float],
    replicates: int, seed: int,
) -> dict[str, float]:
    donors = sorted(set(candidate) & set(baseline))
    if set(candidate) != set(baseline) or len(donors) < 3:
        raise RouterV2EvaluationError("paired router bootstrap donor roster differs")
    candidate_values = np.asarray([candidate[donor] for donor in donors], dtype=float)
    baseline_values = np.asarray([baseline[donor] for donor in donors], dtype=float)
    rng = np.random.default_rng(int(seed))
    differences = np.empty(int(replicates), dtype=float)
    for index in range(int(replicates)):
        sampled = rng.integers(0, len(donors), size=len(donors))
        differences[index] = np.mean(candidate_values[sampled] - baseline_values[sampled])
    return {
        "estimate": float(np.mean(candidate_values - baseline_values)),
        "ci_low": float(np.quantile(differences, 0.025)),
        "ci_high": float(np.quantile(differences, 0.975)),
        "valid_replicates": int(replicates),
    }


def _per_donor_f1(frame, prediction: str, minimum_classes: int):
    output = {}
    for donor, group in frame.groupby("donor_id", sort=True):
        if group["truth"].nunique() < minimum_classes:
            continue
        output[str(donor)] = macro_f1(group["truth"], group[prediction])
    return output


def _per_class(frame, prediction: str, policy: dict[str, Any]):
    from sklearn.metrics import f1_score

    rows = []
    validation = policy["validation"]
    for label, group in frame.groupby("truth", sort=True):
        donor_count = int(group["donor_id"].nunique())
        cells = int(len(group))
        evaluable = bool(
            donor_count >= validation["minimum_donors_for_class_evaluation"]
            and cells >= validation["minimum_cells_for_class_evaluation"]
        )
        truth = frame["truth"].astype(str).to_numpy() == str(label)
        predicted = frame[prediction].astype(str).to_numpy() == str(label)
        value = float(f1_score(truth, predicted, zero_division=0))
        rows.append({
            "label": str(label), "cells": cells, "donors": donor_count,
            "f1": value, "evaluable": evaluable,
            "pass": bool(
                not evaluable
                or value >= validation["minimum_evaluable_class_f1"]
            ),
        })
    return rows


def evaluate_gse212837_router_v2_diagnostic(
    config: dict[str, Any], policy_value: str | Path,
    router_manifest_value: str | Path, v40_policy_value: str | Path,
    v40_mapping_manifest_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    import h5py
    import pandas as pd

    policy_path, policy, _, source = load_router_v2_policy(config, policy_value)
    manifest_path = Path(router_manifest_value).resolve()
    with manifest_path.open() as handle:
        manifest = json.load(handle)
    if (
        manifest.get("schema_version") != "masld-cl-router-v2-gse212837-diagnostic-v43"
        or manifest.get("policy", {}).get("sha256") != sha256_path(policy_path)
        or manifest.get("author_labels_or_conditions_opened") is not False
        or manifest.get("diagnostic_only_no_method_or_threshold_revision") is not True
    ):
        raise RouterV2EvaluationError("Router V2 diagnostic mapping identity differs")
    cells_path = manifest_path.parent / manifest["cells_file"]
    if sha256_path(cells_path) != manifest["cells_sha256"]:
        raise RouterV2EvaluationError("Router V2 diagnostic cell file changed")
    routed = pd.read_csv(cells_path, sep="\t", compression="gzip", dtype=str)
    if routed["cell_id"].duplicated().any() or len(routed) != policy["diagnostic_gse212837"]["cells"]:
        raise RouterV2EvaluationError("Router V2 diagnostic cell roster differs")

    v40_policy_path = Path(v40_policy_value).resolve()
    with v40_policy_path.open() as handle:
        v40_policy = json.load(handle)
    author_map = v40_policy.get("evaluation_only", {}).get("author_label_map", {})
    if (
        v40_policy.get("schema_version") != "masld-cl-external-gse212837-policy-v40"
        or set(author_map) != {
            "CentralHep", "InterHep", "PortalHep", "Cholangiocyte",
            "Stellate", "NKTcell", "PortalEndo", "cvLSECs",
        }
    ):
        raise RouterV2EvaluationError("V40 author-label audit map differs")
    collapsed_author_map = {
        "CentralHep": "Hepatocytes", "InterHep": "Hepatocytes",
        "PortalHep": "Hepatocytes", "Cholangiocyte": "Cholangiocytes",
        "Stellate": "Fibroblasts", "NKTcell": "NK-T",
        "PortalEndo": "Endothelial cells", "cvLSECs": "Endothelial cells",
    }
    with h5py.File(source, "r") as handle:
        author = pd.DataFrame({
            "cell_id": _read_frame_column(handle["obs"], "_index"),
            "author_label": _read_frame_column(handle["obs"], "celltype_pred"),
        })
    routed = routed.merge(author, on="cell_id", validate="one_to_one")
    routed["truth"] = routed["author_label"].map(collapsed_author_map)
    if routed["truth"].isna().any():
        raise RouterV2EvaluationError("GSE212837 author label is unmapped")

    baseline_manifest_path = Path(v40_mapping_manifest_value).resolve()
    with baseline_manifest_path.open() as handle:
        baseline_manifest = json.load(handle)
    if baseline_manifest.get("schema_version") != "masld-cl-external-gse212837-model-v40":
        raise RouterV2EvaluationError("V40 baseline mapping identity differs")
    baseline_path = baseline_manifest_path.parent / baseline_manifest["cells_file"]
    if sha256_path(baseline_path) != baseline_manifest["cells_sha256"]:
        raise RouterV2EvaluationError("V40 baseline routing cells changed")
    baseline = pd.read_csv(baseline_path, sep="\t", compression="gzip", dtype=str)
    baseline = baseline[["cell_id", "routing_label"]]
    baseline["baseline_collapsed"] = collapse_router_labels(
        baseline["routing_label"], policy["collapsed_labels"],
    )
    routed = routed.merge(
        baseline[["cell_id", "baseline_collapsed"]], on="cell_id", validate="one_to_one"
    )

    minimum_classes = int(policy["validation"]["minimum_classes_per_donor"])
    donor_scores = {
        prediction: _per_donor_f1(routed, prediction, minimum_classes)
        for prediction in ("majority_collapsed", "raw_collapsed", "baseline_collapsed")
    }
    if len(set(map(tuple, (sorted(value) for value in donor_scores.values())))) != 1:
        raise RouterV2EvaluationError("router diagnostic donor metric rosters differ")
    bootstrap = paired_donor_bootstrap(
        donor_scores["majority_collapsed"], donor_scores["baseline_collapsed"],
        policy["validation"]["bootstrap_replicates"],
        policy["validation"]["bootstrap_seed"],
    )
    majority_mean = float(np.mean(list(donor_scores["majority_collapsed"].values())))
    raw_mean = float(np.mean(list(donor_scores["raw_collapsed"].values())))
    baseline_mean = float(np.mean(list(donor_scores["baseline_collapsed"].values())))
    classes = _per_class(routed, "majority_collapsed", policy)
    result = {
        "schema_version": "masld-cl-router-v2-gse212837-evaluation-v43",
        "config_sha256": config["_config_sha256"],
        "scientific_role": "diagnostic_only_no_selection_or_method_revision",
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "mapping": {"path": str(manifest_path), "sha256": sha256_path(manifest_path)},
        "baseline_mapping": {
            "path": str(baseline_manifest_path), "sha256": sha256_path(baseline_manifest_path),
        },
        "donors_evaluated": len(donor_scores["majority_collapsed"]),
        "donor_balanced_macro_f1": {
            "majority_voting": majority_mean,
            "raw_prediction": raw_mean,
            "v40_latent_lda": baseline_mean,
            "locked_minimum": policy["validation"]["minimum_donor_balanced_macro_f1"],
        },
        "paired_improvement_majority_vs_v40": bootstrap,
        "per_class": classes,
        "diagnostic_gate_pass": bool(
            majority_mean >= policy["validation"]["minimum_donor_balanced_macro_f1"]
            and bootstrap["ci_low"] > 0
            and all(row["pass"] for row in classes)
        ),
        "author_labels_opened_only_after_mapping_complete": True,
        "result_may_change_router_or_thresholds": False,
    }
    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(output / "router_diagnostic.json", result)
    with (output / "per_donor_metrics.tsv").open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["donor_id", "majority_macro_f1", "raw_macro_f1", "v40_macro_f1"])
        for donor in sorted(donor_scores["majority_collapsed"]):
            writer.writerow([
                donor, donor_scores["majority_collapsed"][donor],
                donor_scores["raw_collapsed"][donor], donor_scores["baseline_collapsed"][donor],
            ])
    return result
