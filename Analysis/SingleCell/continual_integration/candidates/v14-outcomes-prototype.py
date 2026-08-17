#!/usr/bin/env python3
"""Post-lock label and matched disease-geometry gates for the V14 pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from masld_cl.benchmark import _label_scores
from masld_cl.config import canonical_json_bytes, load_config, write_json_exclusive
from masld_cl.contracts import ContractError, sha256_path
from masld_cl.embedding import load_embedding
from masld_cl.orthogonal_bridge_outcome import (
    _aligned_candidate_and_raw,
    _scope_centroids,
    evaluate_matched_scope,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--best-non-cl-label-benchmark", required=True)
    parser.add_argument("--raw-pca-embedding", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    selection_path = Path(args.selection_lock).resolve()
    with selection_path.open() as handle:
        selection = json.load(handle)
    payload = {key: value for key, value in selection.items() if key != "lock_sha256"}
    if (
        selection.get("schema_version")
        != "masld-cl-geometry-preserving-continual-adapter-selection-v14"
        or selection.get("config_sha256") != config["_config_sha256"]
        or selection.get("selection_frozen") is not True
        or selection.get("outcomes_unlocked") is not True
        or selection.get("selected") is None
        or hashlib.sha256(canonical_json_bytes(payload)).hexdigest()
        != selection.get("lock_sha256")
    ):
        raise ContractError("V14 selection lock is invalid")

    candidate_path = Path(selection["selected"]["embedding"]).resolve()
    if sha256_path(candidate_path) != selection["selected"]["embedding_sha256"]:
        raise ContractError("selected V14 embedding changed")
    candidate_bundle = load_embedding(candidate_path)

    label_path = Path(args.best_non_cl_label_benchmark).resolve()
    with label_path.open() as handle:
        label_baseline = json.load(handle)
    if (
        label_baseline.get("schema_version")
        != "masld-cl-best-non-cl-label-benchmark-v10"
        or label_baseline.get("config_sha256") != config["_config_sha256"]
        or label_baseline.get("best_non_cl_method") != "replay_only"
    ):
        raise ContractError("best non-CL label authority differs")
    best = label_baseline["best_non_cl_method"]
    baseline = label_baseline["method_scores"][best]
    candidate_overall, candidate_lineages = _label_scores(
        candidate_bundle, config["lineages"]
    )
    macro_change = candidate_overall - baseline["mean_donor_macro_f1"]
    lineage_changes = {
        lineage: candidate_lineages[lineage]
        - baseline["mean_donor_lineage_f1"][lineage]
        for lineage in config["lineages"]
    }
    label_gates = {
        "query_macro_f1_noninferior": macro_change
        >= -config["gates"]["query_macro_f1_margin"],
        "all_major_lineages_noninferior": min(lineage_changes.values())
        >= -config["gates"]["major_lineage_f1_margin"],
    }

    raw_path = Path(args.raw_pca_embedding).resolve()
    candidate_info, candidate, raw, cells = _aligned_candidate_and_raw(
        candidate_path, raw_path
    )
    if candidate_info.get("method") != "geometry_preserving_continual_adapter":
        raise ContractError("selected V14 embedding method differs")
    scopes = []
    disease_gates = []
    for lineage in config["lineages"]:
        for study in [None, *config["evaluation"]["powered_query_studies"]]:
            scope = (
                f"POOLED_PRIMARY|{lineage}" if study is None else f"{study}|{lineage}"
            )
            candidate_centroids, raw_centroids, _, controls, n_cells = _scope_centroids(
                candidate, raw, cells, lineage, study
            )
            result = evaluate_matched_scope(
                candidate_centroids, raw_centroids, controls
            )
            detail = {
                "scope": scope,
                "n_cells": n_cells,
                "n_donors": int(len(controls)),
                "n_controls": int(controls.sum()),
                "n_cases": int((~controls).sum()),
            }
            if result is None:
                scopes.append(
                    {
                        **detail,
                        "testable": False,
                        "reason": "fewer_than_three_donors_per_group",
                    }
                )
                continue
            scopes.append({**detail, "testable": True, **result})
            if study is None:
                disease_gates.append(
                    {
                        "metric": "pooled_disease_retention",
                        "scope": lineage,
                        "observed": result["retention"],
                        "threshold": 0.9,
                        "pass": result["retention"] >= 0.9,
                    }
                )
            else:
                disease_gates.extend(
                    [
                        {
                            "metric": "per_study_disease_retention",
                            "scope": scope,
                            "observed": result["retention"],
                            "threshold": 0.8,
                            "pass": result["retention"] >= 0.8,
                        },
                        {
                            "metric": "within_study_distance_spearman",
                            "scope": scope,
                            "observed": result["distance_spearman"],
                            "threshold": 0.9,
                            "pass": result["distance_spearman"] >= 0.9,
                        },
                    ]
                )

    result = {
        "schema_version": "masld-cl-v14-pilot-outcomes-v1",
        "config_sha256": config["_config_sha256"],
        "selection_lock": {
            "path": str(selection_path),
            "sha256": sha256_path(selection_path),
            "lock_sha256": selection["lock_sha256"],
        },
        "selected_setting": selection["selected"]["setting"],
        "candidate_embedding": {
            "path": str(candidate_path),
            "sha256": sha256_path(candidate_path),
        },
        "label_concordance": {
            "best_non_cl_method": best,
            "baseline_authority": {
                "path": str(label_path),
                "sha256": sha256_path(label_path),
            },
            "candidate_mean_donor_macro_f1": candidate_overall,
            "baseline_mean_donor_macro_f1": baseline["mean_donor_macro_f1"],
            "query_macro_f1_change": macro_change,
            "candidate_mean_donor_lineage_f1": candidate_lineages,
            "major_lineage_f1_changes": lineage_changes,
            "gates": label_gates,
            "pass": all(label_gates.values()),
            "labels_are_audit_targets_and_were_not_overwritten": True,
        },
        "disease_preservation": {
            "raw_pca_embedding": {
                "path": str(raw_path),
                "sha256": sha256_path(raw_path),
            },
            "scopes": scopes,
            "gates": disease_gates,
            "pass": bool(disease_gates) and all(row["pass"] for row in disease_gates),
            "language_constraint": "cross-sectional stage-associated remodeling",
        },
    }
    result["pilot_pass"] = bool(
        result["label_concordance"]["pass"]
        and result["disease_preservation"]["pass"]
    )
    write_json_exclusive(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
