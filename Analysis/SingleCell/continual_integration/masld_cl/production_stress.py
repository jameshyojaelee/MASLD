"""Evaluate the six-model V35 production expansion without hiding leakage."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .embedding import load_embedding
from .firewall import validate_program_firewall
from .metrics import (
    donor_centroids,
    donor_distance_spearman,
    normalized_shift_control,
    standardized_case_control_separation,
)
from .secondary_stress import (
    _aligned_centroids,
    _bootstrap_alignment_improvement,
    _centroids,
    _donor_stage_map,
)


POLICY_SCHEMA = "masld-cl-production-stress-policy-v36"


def _load_policy(config: dict[str, Any], value: str | Path):
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "production_stress_policy_v36.json"
    )
    if path != expected:
        raise ContractError("production stress requires the source-controlled V36 policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("studies") != config["evaluation"]["secondary_stress_studies"]
        or policy.get("minimum_donors_per_group") != config["evaluation"]["minimum_donors_per_group"]
        or policy.get("status", {}).get("GSE136103_used_during_method_development") is not True
        or policy.get("status", {}).get("results_cannot_establish_independent_secondary_confirmation") is not True
        or set(policy.get("model_to_scope", {}).values()) != {"all_lineage", *config["lineages"]}
    ):
        raise ContractError("production-stress policy identity differs")
    sources = {}
    for name, source in policy["sources"].items():
        source_path = (path.parent / source["path"]).resolve()
        if sha256_path(source_path) != source["sha256"]:
            raise ContractError(f"production-stress source changed: {name}")
        sources[name] = source_path
    return path, policy, sources


def _full_positions(full_cells, model_cells):
    lookup = {
        value: index for index, value in enumerate(full_cells["cell_id"].astype(str))
    }
    positions = np.asarray(
        [lookup.get(value, -1) for value in model_cells["cell_id"].astype(str)],
        dtype=np.int64,
    )
    if np.any(positions < 0) or len(set(positions.tolist())) != len(positions):
        raise ContractError("production-stress model roster is not a unique full-atlas subset")
    return positions


def evaluate_production_stress(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path, policy, sources = _load_policy(config, policy_value)
    with sources["candidate_bundle"].open() as handle:
        bundle = json.load(handle)
    if (
        bundle.get("schema_version") != "masld-cl-production-expansion-bundle-v35"
        or set(bundle.get("models", {})) != set(policy["model_to_scope"])
        or bundle.get("secondary_stress_status")
        != "method_development_not_independent_confirmation"
    ):
        raise ContractError("production-stress candidate bundle identity differs")
    _, harmony, full_cells = load_embedding(sources["harmony_embedding"])
    _, raw, raw_cells = load_embedding(sources["raw_pca_embedding"])
    if not np.array_equal(
        full_cells["cell_id"].astype(str).to_numpy(),
        raw_cells["cell_id"].astype(str).to_numpy(),
    ):
        raise ContractError("production-stress comparator rosters differ")
    stages = _donor_stage_map(sources["donor_manifest"])
    minimum = int(policy["minimum_donors_per_group"])
    alignments, disease, untestable = [], [], []
    for model_index, (model_kind, scope) in enumerate(policy["model_to_scope"].items()):
        source = bundle["models"][model_kind]
        for key in ("embedding_manifest", "production_manifest"):
            if sha256_path(source[key]) != source[f"{key}_sha256"]:
                raise ContractError(f"production-stress bundle source changed: {model_kind}.{key}")
        with Path(source["production_manifest"]).open() as handle:
            production = json.load(handle)
        candidate_info, candidate, cells = load_embedding(source["embedding_manifest"])
        if (
            production.get("model_kind") != model_kind
            or production.get("embedding") != candidate_info
        ):
            raise ContractError(f"production-stress model provenance differs: {model_kind}")
        positions = _full_positions(full_cells, cells)
        model_harmony = harmony[positions]
        model_raw = raw[positions]
        datasets = cells["dataset"].astype(str).to_numpy()
        preparations = cells["preparation"].astype(str).to_numpy()
        eligible = cells["analysis_eligible"].to_numpy(dtype=bool)
        reference = cells["strict_reference"].to_numpy(dtype=bool)
        healthy = cells["donor_id"].astype(str).map(stages).astype(str).to_numpy() == "Healthy"
        for study_index, study in enumerate(policy["studies"]):
            study_mask = eligible & (datasets == study)
            for protocol_index, protocol in enumerate(["all", *sorted(set(preparations[study_mask]))]):
                protocol_mask = (
                    np.ones(len(cells), dtype=bool)
                    if protocol == "all" else preparations == protocol
                )
                controls = study_mask & healthy & protocol_mask
                compatible_reference = reference & protocol_mask
                n_controls = cells.loc[controls, "donor_id"].astype(str).nunique()
                n_reference = cells.loc[compatible_reference, "donor_id"].astype(str).nunique()
                alignment_scope = f"{study}|{scope}|{protocol}"
                if n_controls < minimum or n_reference < minimum:
                    untestable.append({
                        "metric_family": "control_alignment", "scope": alignment_scope,
                        "n_controls": int(n_controls), "n_reference": int(n_reference),
                        "reason": "fewer_than_three_donors_in_a_required_group",
                    })
                    continue
                reference_ids, candidate_reference = _centroids(
                    candidate, cells, compatible_reference
                )
                control_ids, candidate_controls = _centroids(candidate, cells, controls)
                harmony_reference = _aligned_centroids(
                    model_harmony, cells, compatible_reference, reference_ids
                )
                harmony_controls = _aligned_centroids(
                    model_harmony, cells, controls, control_ids
                )
                candidate_shift = normalized_shift_control(
                    candidate_reference, candidate_controls
                )
                harmony_shift = normalized_shift_control(
                    harmony_reference, harmony_controls
                )
                improvement = (harmony_shift - candidate_shift) / harmony_shift
                reference_meta = cells.loc[compatible_reference].drop_duplicates("donor_id")
                strata = reference_meta.set_index("donor_id").loc[reference_ids]["dataset"]
                bootstrap = _bootstrap_alignment_improvement(
                    candidate_reference, candidate_controls,
                    harmony_reference, harmony_controls,
                    strata.astype(str).to_numpy(),
                    int(policy["bootstrap_replicates"]),
                    int(policy["bootstrap_seed"] + 1000 * model_index
                        + 100 * study_index + protocol_index),
                )
                alignments.append({
                    "scope": alignment_scope, "model_kind": model_kind,
                    "study": study, "lineage": scope, "protocol": protocol,
                    "n_controls": int(n_controls), "n_reference": int(n_reference),
                    "candidate_shift_control": float(candidate_shift),
                    "harmony_shift_control": float(harmony_shift),
                    "improvement_vs_harmony": float(improvement),
                    **bootstrap,
                    "no_more_than_10_percent_worsening_pass": bool(
                        improvement >= -policy["gates"]["maximum_control_alignment_worsening"]
                    ),
                })

            donor_ids, candidate_values = _centroids(candidate, cells, study_mask)
            raw_values = _aligned_centroids(model_raw, cells, study_mask, donor_ids)
            controls = np.asarray([stages[donor] == "Healthy" for donor in donor_ids])
            n_controls, n_cases = int(controls.sum()), int((~controls).sum())
            disease_scope = f"{study}|{scope}"
            supporting_spearman = donor_distance_spearman(candidate_values, raw_values)
            if n_controls < minimum or n_cases < minimum:
                untestable.append({
                    "metric_family": "disease_preservation", "scope": disease_scope,
                    "n_controls": n_controls, "n_cases": n_cases,
                    "supporting_distance_spearman": float(supporting_spearman),
                    "reason": "fewer_than_three_donors_in_a_required_group",
                })
                continue
            candidate_separation = standardized_case_control_separation(
                candidate_values[~controls], candidate_values[controls]
            )
            raw_separation = standardized_case_control_separation(
                raw_values[~controls], raw_values[controls]
            )
            retention = candidate_separation / raw_separation
            disease.append({
                "scope": disease_scope, "model_kind": model_kind,
                "study": study, "lineage": scope,
                "n_controls": n_controls, "n_cases": n_cases,
                "retention": float(retention),
                "distance_spearman": float(supporting_spearman),
                "retention_pass": bool(retention >= policy["gates"]["minimum_disease_retention"]),
                "distance_pass": bool(supporting_spearman >= policy["gates"]["minimum_distance_spearman"]),
            })

    technical_gates = {
        "all_evaluable_alignment_scopes_pass": bool(alignments) and all(
            row["no_more_than_10_percent_worsening_pass"] for row in alignments
        ),
        "all_evaluable_disease_scopes_pass": bool(disease) and all(
            row["retention_pass"] and row["distance_pass"] for row in disease
        ),
        "GSE174748_prespecified_underpowered": all(
            row["study"] != "GSE174748" for row in alignments + disease
        ),
    }
    result = {
        "schema_version": "masld-cl-production-stress-v36",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "alignment": alignments,
        "disease_preservation": disease,
        "untestable": untestable,
        "technical_gates": technical_gates,
        "technical_stress_pass": all(technical_gates.values()),
        "independent_secondary_confirmation_pass": False,
        "independent_confirmation_reason": "GSE136103 informed V35 method development and GSE174748 has only two donors per condition",
    }
    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(output / "production_stress.json", result)
    rows = []
    for row in alignments:
        rows.extend([
            {"metric": "alignment_improvement_vs_harmony", "scope": row["scope"], "value": row["improvement_vs_harmony"]},
            {"metric": "alignment_improvement_vs_harmony_ci_low", "scope": row["scope"], "value": row["ci_low"]},
        ])
    for row in disease:
        rows.extend([
            {"metric": "per_study_disease_retention", "scope": row["scope"], "value": row["retention"]},
            {"metric": "within_study_distance_spearman", "scope": row["scope"], "value": row["distance_spearman"]},
        ])
    rows.extend([
        {"metric": "technical_secondary_stress_pass", "scope": "evaluable_secondary_studies", "value": int(result["technical_stress_pass"])},
        {"metric": "independent_secondary_confirmation_pass", "scope": "independent_secondary_studies", "value": 0},
    ])
    metrics_path = output / "metrics_long.tsv"
    with metrics_path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "scope", "value"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return result
