"""Donor-level evaluation of the locked GSE212837 external mapping."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .embedding import load_embedding
from .external_gse212837 import (
    ExternalGSE212837Error,
    _read_frame_column,
    load_external_gse212837_policy,
)
from .firewall import validate_program_firewall
from .metrics import (
    donor_centroids,
    donor_distance_spearman,
    macro_f1,
    normalized_shift_control,
    standardized_case_control_separation,
)
from .secondary_stress import _bootstrap_alignment_improvement


def _load_model_cells(path: Path):
    import pandas as pd

    cells = pd.read_csv(path, sep="\t", compression="gzip", dtype=str)
    if cells["cell_id"].duplicated().any() or not np.array_equal(
        cells["row_index"].astype(int).to_numpy(), np.arange(len(cells))
    ):
        raise ExternalGSE212837Error("external model cell roster differs")
    return cells


def _aligned_donor_centroids(values, cells, expected_donors=None):
    donors, centroids = donor_centroids(values, cells["donor_id"].astype(str))
    donors = donors.astype(str)
    if expected_donors is None:
        return donors, centroids
    lookup = {donor: index for index, donor in enumerate(donors)}
    if set(donors) != set(expected_donors):
        raise ExternalGSE212837Error("external representations have different donor rosters")
    return donors, centroids[np.asarray([lookup[donor] for donor in expected_donors])]


def evaluate_external_gse212837(
    config: dict[str, Any], policy_value: str | Path,
    bundle_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    import h5py
    import pandas as pd

    validate_program_firewall(config)
    policy_path, policy, sources = load_external_gse212837_policy(config, policy_value)
    bundle_path = Path(bundle_value).resolve()
    with bundle_path.open() as handle:
        bundle = json.load(handle)
    if (
        bundle.get("schema_version") != "masld-cl-external-gse212837-bundle-v40"
        or bundle.get("policy", {}).get("sha256") != sha256_path(policy_path)
        or bundle.get("query_labels_or_conditions_read") is not False
        or bundle.get("evaluation_locked_until_bundle_complete") is not True
        or set(bundle.get("models", {})) != set(policy["models"])
    ):
        raise ExternalGSE212837Error("external mapping bundle identity differs")

    evaluation = policy["evaluation_only"]
    control_donors = set(evaluation["control_donors"])
    case_donors = set(evaluation["case_donors"])
    expected_donors = control_donors | case_donors
    if len(control_donors) != 3 or len(case_donors) != 9 or control_donors & case_donors:
        raise ExternalGSE212837Error("external evaluation donor split differs")

    # Evaluation-only author labels are opened only after a complete mapping bundle exists.
    with h5py.File(sources["external_h5ad"], "r") as handle:
        source_cell_ids = _read_frame_column(handle["obs"], "_index")
        source_donors = _read_frame_column(handle["obs"], "Donor")
        author_labels = _read_frame_column(handle["obs"], "celltype_pred")
    if set(source_donors) != expected_donors:
        raise ExternalGSE212837Error("external donor identities differ at evaluation")
    author_map = evaluation["author_label_map"]
    if set(author_labels) != set(author_map):
        raise ExternalGSE212837Error("external author-label vocabulary differs")
    source = pd.DataFrame({
        "cell_id": source_cell_ids,
        "source_donor": source_donors,
        "audit_label": [author_map[value] for value in author_labels],
    }).set_index("cell_id")

    alignment_rows, disease_rows, untestable = [], [], []
    model_payloads: dict[str, tuple[Any, Any, Any]] = {}
    for model_index, (model_kind, spec) in enumerate(policy["models"].items()):
        model_source = bundle["models"][model_kind]
        manifest_path = Path(model_source["manifest"])
        if sha256_path(manifest_path) != model_source["manifest_sha256"]:
            raise ExternalGSE212837Error(f"external model manifest changed: {model_kind}")
        with manifest_path.open() as handle:
            manifest = json.load(handle)
        model_dir = manifest_path.parent
        coordinates_path = model_dir / manifest["coordinates_file"]
        cells_path = model_dir / manifest["cells_file"]
        if (
            manifest.get("schema_version") != "masld-cl-external-gse212837-model-v40"
            or manifest.get("model_kind") != model_kind
            or manifest.get("query_labels_or_conditions_read") is not False
            or manifest.get("author_labels_read") is not False
            or manifest.get("reference_coordinates_written_or_changed") is not False
            or sha256_path(coordinates_path) != manifest["coordinates_sha256"]
            or sha256_path(cells_path) != manifest["cells_sha256"]
        ):
            raise ExternalGSE212837Error(f"external model provenance differs: {model_kind}")
        cells = _load_model_cells(cells_path)
        joined = cells.join(source, on="cell_id", how="left", validate="one_to_one")
        if joined[["source_donor", "audit_label"]].isna().any().any() or not np.array_equal(
            joined["donor_id"].astype(str).to_numpy(),
            joined["source_donor"].astype(str).to_numpy(),
        ):
            raise ExternalGSE212837Error(f"external cell metadata mismatch: {model_kind}")
        arrays = np.load(coordinates_path)
        unadapted = arrays["unadapted"]
        adapted = arrays["adapted"]
        raw = arrays["raw_study_pca"]
        if not (
            unadapted.shape == adapted.shape == raw.shape == (len(cells), 30)
            and np.isfinite(unadapted).all() and np.isfinite(adapted).all()
            and np.isfinite(raw).all()
        ):
            raise ExternalGSE212837Error(f"external coordinate arrays differ: {model_kind}")
        model_payloads[model_kind] = (cells, joined, arrays)
        if not spec["independently_evaluable"]:
            untestable.append({
                "model_kind": model_kind,
                "metric_families": ["control_alignment", "disease_preservation"],
                "reason": spec["untestable_reason"],
                "mapped_cells": len(cells),
            })
            continue

        reference_info, reference, reference_cells = load_embedding(
            sources[f"reference.{model_kind}"]
        )
        expected_reference_kind = (
            "all_lineage" if spec["lineage_label"] is None else spec["lineage_label"]
        )
        if reference_info.get("model_kind") != expected_reference_kind:
            raise ExternalGSE212837Error(f"external reference model differs: {model_kind}")
        reference_ids, reference_centroids = donor_centroids(
            reference, reference_cells["donor_id"].astype(str)
        )
        query_ids, adapted_centroids = donor_centroids(adapted, cells["donor_id"].astype(str))
        baseline_ids, unadapted_centroids = donor_centroids(
            unadapted, cells["donor_id"].astype(str)
        )
        raw_ids, raw_centroids = donor_centroids(raw, cells["donor_id"].astype(str))
        query_ids = query_ids.astype(str)
        if not (
            np.array_equal(query_ids, baseline_ids.astype(str))
            and np.array_equal(query_ids, raw_ids.astype(str))
        ):
            raise ExternalGSE212837Error(f"external donor orders differ: {model_kind}")
        control_mask = np.asarray([donor in control_donors for donor in query_ids])
        case_mask = np.asarray([donor in case_donors for donor in query_ids])
        minimum = int(evaluation["minimum_donors_per_group"])
        if min(int(control_mask.sum()), int(case_mask.sum())) < minimum:
            raise ExternalGSE212837Error(f"evaluable external model is underpowered: {model_kind}")

        candidate_shift = normalized_shift_control(
            reference_centroids, adapted_centroids[control_mask]
        )
        baseline_shift = normalized_shift_control(
            reference_centroids, unadapted_centroids[control_mask]
        )
        improvement = (baseline_shift - candidate_shift) / baseline_shift
        reference_meta = reference_cells.drop_duplicates("donor_id").set_index("donor_id")
        reference_strata = reference_meta.loc[reference_ids, "dataset"].astype(str).to_numpy()
        bootstrap = _bootstrap_alignment_improvement(
            reference_centroids, adapted_centroids[control_mask],
            reference_centroids, unadapted_centroids[control_mask],
            reference_strata, int(evaluation["bootstrap_replicates"]),
            int(evaluation["bootstrap_seed"] + 100 * model_index),
        )
        alignment_pass = bool(
            improvement >= evaluation["gates"]["minimum_alignment_improvement"]
            and bootstrap["ci_low"] > evaluation["gates"]["alignment_ci_low_must_exceed"]
        )
        alignment_rows.append({
            "model_kind": model_kind,
            "n_reference_donors": int(len(reference_ids)),
            "n_control_donors": int(control_mask.sum()),
            "adapted_shift_control": float(candidate_shift),
            "unadapted_shift_control": float(baseline_shift),
            "improvement_vs_condition_blind_unadapted_mapping": float(improvement),
            **bootstrap,
            "pass": alignment_pass,
        })

        candidate_sep = standardized_case_control_separation(
            adapted_centroids[case_mask], adapted_centroids[control_mask]
        )
        raw_sep = standardized_case_control_separation(
            raw_centroids[case_mask], raw_centroids[control_mask]
        )
        retention = candidate_sep / raw_sep
        distance = donor_distance_spearman(adapted_centroids, raw_centroids)
        disease_pass = bool(
            retention >= evaluation["gates"]["minimum_disease_retention"]
            and distance >= evaluation["gates"]["minimum_distance_spearman"]
        )
        disease_rows.append({
            "model_kind": model_kind,
            "n_controls": int(control_mask.sum()),
            "n_cases": int(case_mask.sum()),
            "candidate_separation": float(candidate_sep),
            "external_study_only_raw_pca_separation": float(raw_sep),
            "retention": float(retention),
            "distance_spearman": float(distance),
            "pass": disease_pass,
        })

    # Routing is audited on all cells and was never available to the mapper.
    _, all_joined, _ = model_payloads["all_lineage"]
    donor_f1 = []
    for donor, group in all_joined.groupby("donor_id", sort=True):
        donor_f1.append({
            "donor_id": str(donor),
            "macro_f1": macro_f1(group["audit_label"], group["routing_label"]),
            "cells": int(len(group)),
        })
    routing_mean = float(np.mean([row["macro_f1"] for row in donor_f1]))
    routing_pass = bool(
        routing_mean >= evaluation["gates"]["minimum_donor_balanced_routing_macro_f1"]
    )

    required = set(evaluation["gates"]["required_evaluable_models"])
    observed_alignment = {row["model_kind"] for row in alignment_rows}
    observed_disease = {row["model_kind"] for row in disease_rows}
    required_models_present = observed_alignment == required and observed_disease == required
    gates = {
        "required_evaluable_models_present": required_models_present,
        "all_required_control_alignment_gates": required_models_present and all(
            row["pass"] for row in alignment_rows
        ),
        "all_required_disease_preservation_gates": required_models_present and all(
            row["pass"] for row in disease_rows
        ),
        "reference_only_routing_gate": routing_pass,
        "macrophages_prespecified_untestable": any(
            row["model_kind"] == "macrophages" for row in untestable
        ),
        "t_cells_prespecified_untestable": any(
            row["model_kind"] == "t_cells" for row in untestable
        ),
    }
    result = {
        "schema_version": "masld-cl-external-gse212837-evaluation-v40",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "mapping_bundle": {"path": str(bundle_path), "sha256": sha256_path(bundle_path)},
        "cohort": {
            "cells": len(source), "biological_donors": len(expected_donors),
            "controls": len(control_donors), "cases": len(case_donors),
            "inferential_unit": "biological donor",
        },
        "control_alignment": alignment_rows,
        "disease_preservation": disease_rows,
        "routing_audit": {
            "mean_donor_macro_f1": routing_mean,
            "minimum": evaluation["gates"]["minimum_donor_balanced_routing_macro_f1"],
            "pass": routing_pass,
            "by_donor": donor_f1,
            "author_labels_used_for_mapping": False,
        },
        "untestable": untestable,
        "gates": gates,
        "independent_external_confirmation_pass": all(gates.values()),
        "harmony_comparison": {
            "applicable": False,
            "reason": policy["status"]["reason_harmony_not_applicable"],
        },
        "external_results_may_trigger_method_revision": False,
    }
    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(output / "external_evaluation.json", result)
    rows = []
    for row in alignment_rows:
        rows.extend([
            {"metric": "alignment_improvement", "scope": row["model_kind"], "value": row["improvement_vs_condition_blind_unadapted_mapping"]},
            {"metric": "alignment_improvement_ci_low", "scope": row["model_kind"], "value": row["ci_low"]},
        ])
    for row in disease_rows:
        rows.extend([
            {"metric": "disease_retention", "scope": row["model_kind"], "value": row["retention"]},
            {"metric": "donor_distance_spearman", "scope": row["model_kind"], "value": row["distance_spearman"]},
        ])
    rows.extend([
        {"metric": "routing_macro_f1", "scope": "all_lineage", "value": routing_mean},
        {"metric": "independent_external_confirmation_pass", "scope": "GSE212837", "value": int(result["independent_external_confirmation_pass"])},
    ])
    with (output / "metrics_long.tsv").open("x", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["metric", "scope", "value"],
            delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return result
