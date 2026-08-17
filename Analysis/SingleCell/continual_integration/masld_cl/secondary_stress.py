"""Donor-level secondary-cohort stress tests for the frozen V33 projection."""

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
    macro_f1,
    normalized_shift_control,
    standardized_case_control_separation,
)


POLICY_SCHEMA = "masld-cl-secondary-stress-policy-v34"


def _load_policy(config: dict[str, Any], value: str | Path):
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "secondary_stress_policy_v34.json"
    )
    if path != expected:
        raise ContractError("secondary stress requires the source-controlled V34 policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("studies") != config["evaluation"]["secondary_stress_studies"]
        or policy.get("minimum_donors_per_group") != config["evaluation"]["minimum_donors_per_group"]
        or policy.get("firewall", {}).get("no_selection_or_refitting") is not True
        or policy.get("firewall", {}).get("donors_are_the_inferential_unit") is not True
    ):
        raise ContractError("secondary-stress policy identity differs")
    sources = {}
    for name, source in policy["sources"].items():
        source_path = (path.parent / source["path"]).resolve()
        if sha256_path(source_path) != source["sha256"]:
            raise ContractError(f"secondary-stress source changed: {name}")
        sources[name] = source_path
    return path, policy, sources


def _donor_stage_map(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {row["donor_id"]: row["harmonized_stage"] for row in rows}
    if len(result) != 104 or len(rows) != 104:
        raise ContractError("secondary stress donor manifest does not contain 104 unique donors")
    return result


def _centroids(latent, cells, mask: np.ndarray):
    selected = cells.loc[mask].reset_index(drop=True)
    donors, values = donor_centroids(
        np.asarray(latent[mask], dtype=np.float64), selected["donor_id"].astype(str)
    )
    return donors.astype(str), values


def _aligned_centroids(latent, cells, mask: np.ndarray, donors: np.ndarray):
    found, values = _centroids(latent, cells, mask)
    lookup = {donor: index for index, donor in enumerate(found)}
    if set(found) != set(donors):
        raise ContractError("secondary-stress representations have different donor rosters")
    return values[np.asarray([lookup[donor] for donor in donors], dtype=np.int64)]


def _bootstrap_alignment_improvement(
    candidate_reference: np.ndarray, candidate_controls: np.ndarray,
    harmony_reference: np.ndarray, harmony_controls: np.ndarray,
    reference_strata: np.ndarray, replicates: int, seed: int,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(replicates):
        reference_draw = []
        for stratum in sorted(set(reference_strata)):
            positions = np.flatnonzero(reference_strata == stratum)
            reference_draw.extend(rng.choice(positions, len(positions), replace=True))
        reference_draw = np.asarray(reference_draw, dtype=np.int64)
        control_draw = rng.choice(
            np.arange(len(candidate_controls)), len(candidate_controls), replace=True
        )
        try:
            candidate_shift = normalized_shift_control(
                candidate_reference[reference_draw], candidate_controls[control_draw]
            )
            harmony_shift = normalized_shift_control(
                harmony_reference[reference_draw], harmony_controls[control_draw]
            )
        except ValueError:
            continue
        if harmony_shift > 0:
            values.append((harmony_shift - candidate_shift) / harmony_shift)
    if len(values) < max(100, int(0.9 * replicates)):
        raise ContractError("secondary alignment bootstrap produced too many degenerate draws")
    values = np.asarray(values, dtype=float)
    return {
        "ci_low": float(np.quantile(values, 0.025)),
        "ci_high": float(np.quantile(values, 0.975)),
        "bootstrap_standard_error": float(np.std(values, ddof=1)),
        "successful_replicates": int(len(values)),
    }


def evaluate_secondary_stress(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path, policy, sources = _load_policy(config, policy_value)
    _, candidate, cells = load_embedding(sources["candidate_embedding"])
    _, harmony, harmony_cells = load_embedding(sources["harmony_embedding"])
    _, raw, raw_cells = load_embedding(sources["raw_pca_embedding"])
    expected_ids = cells["cell_id"].astype(str).to_numpy()
    if not (
        np.array_equal(expected_ids, harmony_cells["cell_id"].astype(str).to_numpy())
        and np.array_equal(expected_ids, raw_cells["cell_id"].astype(str).to_numpy())
    ):
        raise ContractError("secondary-stress full-atlas cell orders differ")
    stages = _donor_stage_map(sources["donor_manifest"])
    donor_stage = cells["donor_id"].astype(str).map(stages)
    if donor_stage.isna().any():
        raise ContractError("secondary-stress cell roster contains an unknown donor")
    datasets = cells["dataset"].astype(str).to_numpy()
    preparations = cells["preparation"].astype(str).to_numpy()
    labels = cells["audit_cell_type"].astype(str).to_numpy()
    eligible = cells["analysis_eligible"].to_numpy(dtype=bool)
    reference = cells["strict_reference"].to_numpy(dtype=bool)
    healthy = donor_stage.astype(str).to_numpy() == "Healthy"

    alignments, disease, untestable = [], [], []
    minimum = int(policy["minimum_donors_per_group"])
    for study_index, study in enumerate(policy["studies"]):
        study_mask = eligible & (datasets == study)
        study_preparations = sorted(set(preparations[study_mask]))
        for scope_index, scope in enumerate(policy["scopes"]):
            lineage = np.ones(len(cells), dtype=bool) if scope == "all_lineage" else labels == scope
            for protocol in ["all", *study_preparations]:
                protocol_mask = (
                    np.ones(len(cells), dtype=bool)
                    if protocol == "all" else preparations == protocol
                )
                control_mask = study_mask & lineage & protocol_mask & healthy
                compatible_reference = reference & lineage & protocol_mask
                control_donors = cells.loc[control_mask, "donor_id"].astype(str).nunique()
                reference_donors = cells.loc[compatible_reference, "donor_id"].astype(str).nunique()
                alignment_scope = f"{study}|{scope}|{protocol}"
                if control_donors < minimum or reference_donors < minimum:
                    untestable.append({
                        "metric_family": "control_alignment", "scope": alignment_scope,
                        "n_controls": int(control_donors),
                        "n_reference": int(reference_donors),
                        "reason": "fewer_than_three_donors_in_a_required_group",
                    })
                    continue
                reference_ids, candidate_reference = _centroids(
                    candidate, cells, compatible_reference
                )
                control_ids, candidate_controls = _centroids(candidate, cells, control_mask)
                harmony_reference = _aligned_centroids(
                    harmony, cells, compatible_reference, reference_ids
                )
                harmony_controls = _aligned_centroids(
                    harmony, cells, control_mask, control_ids
                )
                candidate_shift = normalized_shift_control(
                    candidate_reference, candidate_controls
                )
                harmony_shift = normalized_shift_control(
                    harmony_reference, harmony_controls
                )
                improvement = (harmony_shift - candidate_shift) / harmony_shift
                reference_meta = cells.loc[compatible_reference].drop_duplicates("donor_id")
                reference_dataset = reference_meta.set_index("donor_id").loc[reference_ids]["dataset"]
                bootstrap = _bootstrap_alignment_improvement(
                    candidate_reference, candidate_controls,
                    harmony_reference, harmony_controls,
                    reference_dataset.astype(str).to_numpy(),
                    int(policy["bootstrap_replicates"]),
                    int(policy["bootstrap_seed"] + 100 * study_index + 10 * scope_index
                        + study_preparations.index(protocol) + 1 if protocol != "all"
                        else policy["bootstrap_seed"] + 100 * study_index + 10 * scope_index),
                )
                alignments.append({
                    "scope": alignment_scope,
                    "study": study, "lineage": scope, "protocol": protocol,
                    "n_controls": int(control_donors),
                    "n_reference": int(reference_donors),
                    "candidate_shift_control": float(candidate_shift),
                    "harmony_shift_control": float(harmony_shift),
                    "improvement_vs_harmony": float(improvement),
                    **bootstrap,
                    "no_more_than_10_percent_worsening_pass": bool(
                        improvement >= -policy["gates"]["maximum_control_alignment_worsening"]
                    ),
                })

            full_scope = study_mask & lineage
            donor_ids = sorted(set(cells.loc[full_scope, "donor_id"].astype(str)))
            healthy_ids = [donor for donor in donor_ids if stages[donor] == "Healthy"]
            case_ids = [donor for donor in donor_ids if stages[donor] != "Healthy"]
            disease_scope = f"{study}|{scope}"
            candidate_ids, candidate_values = _centroids(candidate, cells, full_scope)
            raw_values = _aligned_centroids(raw, cells, full_scope, candidate_ids)
            candidate_control = np.asarray([stages[x] == "Healthy" for x in candidate_ids])
            supporting_spearman = donor_distance_spearman(candidate_values, raw_values)
            if len(healthy_ids) < minimum or len(case_ids) < minimum:
                untestable.append({
                    "metric_family": "disease_preservation", "scope": disease_scope,
                    "n_controls": len(healthy_ids), "n_cases": len(case_ids),
                    "supporting_distance_spearman": float(supporting_spearman),
                    "reason": "fewer_than_three_donors_in_a_required_group",
                })
                continue
            candidate_separation = standardized_case_control_separation(
                candidate_values[~candidate_control], candidate_values[candidate_control]
            )
            raw_separation = standardized_case_control_separation(
                raw_values[~candidate_control], raw_values[candidate_control]
            )
            retention = candidate_separation / raw_separation
            disease.append({
                "scope": disease_scope, "study": study, "lineage": scope,
                "n_controls": len(healthy_ids), "n_cases": len(case_ids),
                "candidate_separation": float(candidate_separation),
                "raw_pca_separation": float(raw_separation),
                "retention": float(retention),
                "distance_spearman": float(supporting_spearman),
                "retention_pass": bool(retention >= policy["gates"]["minimum_disease_retention"]),
                "distance_pass": bool(supporting_spearman >= policy["gates"]["minimum_distance_spearman"]),
            })

    route_scores = {}
    for study in policy["studies"]:
        donor_values = []
        for _, group in cells.loc[eligible & (datasets == study)].groupby("donor_id", sort=True):
            donor_values.append(macro_f1(group["audit_cell_type"], group["routing_label"]))
        route_scores[study] = {
            "n_donors": len(donor_values),
            "mean_donor_macro_f1": float(np.mean(donor_values)),
            "audit_only_not_a_selection_gate": True,
        }

    evaluable_studies = sorted({row["study"] for row in alignments + disease})
    result = {
        "schema_version": "masld-cl-secondary-stress-v34",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "alignment": alignments,
        "disease_preservation": disease,
        "untestable": untestable,
        "routing_audit": route_scores,
        "evaluable_studies": evaluable_studies,
        "mapping_only": {
            "study": policy["mapping_only_study"],
            "used_for_fit_or_gate": False,
            "analyzed_cells": int((eligible & (datasets == policy["mapping_only_study"])).sum()),
            "descriptive_only_cells": int((~eligible & (datasets == policy["mapping_only_study"])).sum()),
        },
        "gates": {
            "all_evaluable_alignment_scopes_pass": bool(alignments) and all(
                row["no_more_than_10_percent_worsening_pass"] for row in alignments
            ),
            "all_evaluable_disease_scopes_pass": bool(disease) and all(
                row["retention_pass"] and row["distance_pass"] for row in disease
            ),
            "GSE174748_prespecified_underpowered": all(
                row["study"] != "GSE174748" for row in alignments + disease
            ),
        },
    }
    result["secondary_stress_pass"] = all(result["gates"].values())
    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(output / "secondary_stress.json", result)
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
    rows.append({"metric": "secondary_stress_pass", "scope": "evaluable_secondary_studies", "value": int(result["secondary_stress_pass"])})
    metrics_path = output / "metrics_long.tsv"
    with metrics_path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "scope", "value"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return result
