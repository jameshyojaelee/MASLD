#!/usr/bin/env python3
"""Compute Plan 43's mechanical scientific and promotion verdicts.

This script does not discover alternative hypotheses.  It consumes only the
sealed analyses and records pass, fail, or source-gated skip states.
"""

from __future__ import annotations

import datetime as dt
import json
import math
import platform
import sys
from pathlib import Path

from bridge_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    SCRIPT_ROOT,
    read_tsv,
    require_validated_seal,
    sha256_file,
    write_tsv,
)


def number(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def truth(value: str) -> bool:
    return str(value).strip().lower() == "true"


def main() -> None:
    seal = require_validated_seal()

    human = read_tsv(CANDIDATE_ROOT / "human_reversal_effects.tsv")
    human_primary = {row["dataset_id"]: row for row in human if row["score_scheme"] == "weighted__primary"}
    human_meta = read_tsv(CANDIDATE_ROOT / "human_reversal/human_reversal_meta_analysis.tsv")[0]
    human_sensitivity = read_tsv(CANDIDATE_ROOT / "human_reversal/human_reversal_sensitivity.tsv")
    rygb = [row for row in human_sensitivity if row["sensitivity"] == "RYGB_direction_check"]
    unseen_human = [human_primary[key] for key in ("GSE83452", "GSE48452")]
    all_human_direction = all(truth(row["direction_agrees"]) for row in human_primary.values())
    unseen_nominal = sum(number(row["p"]) < 0.05 and truth(row["direction_agrees"]) for row in unseen_human)
    human_gate = (
        all_human_direction
        and unseen_nominal >= 1
        and int(human_meta["n_effective_cohorts"]) == 3
        and number(human_meta["p"]) < 0.01
        and truth(human_meta["direction_agrees"])
        and any(number(row.get("p", "")) < 0.05 and truth(row["direction_agrees"]) for row in rygb)
        and all(truth(row["direction_agrees"]) for row in rygb)
    )

    crop = read_tsv(CANDIDATE_ROOT / "source_gates/cropseq_gate.tsv")[0]
    crop_gate = crop["status"] == "pass" and truth(crop["inference_authorized"])
    crop_output = [{
        "dataset_id": "GSE281160",
        "analysis_id": "risk_oriented_state_axis_bridge",
        "status": crop["status"],
        "inference_authorized": crop["inference_authorized"],
        "n_source_variants": crop["n_source_variants"],
        "n_eligible_loci": crop["n_independent_corrected_tier12_loci"],
        "estimate": "", "p": "", "fraction_positive": "", "guide_agreement_fraction": "",
        "claim_scope": "none",
        "reason": crop["detail"],
    }]
    write_tsv(
        CANDIDATE_ROOT / "cropseq_risk_bridge.tsv",
        crop_output,
        list(crop_output[0]),
    )

    protein = read_tsv(CANDIDATE_ROOT / "protein_transportability.tsv")
    protein_map = {row["dataset_id"]: row for row in protein}
    protein_sensitivity = read_tsv(CANDIDATE_ROOT / "protein/protein_sensitivity.tsv")
    residual = {
        row["dataset_id"]: row
        for row in protein_sensitivity
        if row["sensitivity"] == "continuous_bulk_t_residualized_abs_z"
    }
    protein_gate = (
        truth(protein_map["PXD052787"]["direction_agrees"])
        and truth(protein_map["PXD051911"]["direction_agrees"])
        and truth(protein_map["PXD052787"]["matching_gate_pass"])
        and truth(protein_map["PXD051911"]["matching_gate_pass"])
        and number(protein_map["PXD052787"]["q"]) < 0.05
        and number(protein_map["equal_cohort_meta"]["p"]) < 0.01
        and all(number(row["estimate"]) > 0 and number(row["p"]) < 0.05 for row in residual.values())
    )

    context = read_tsv(CANDIDATE_ROOT / "context_topology.tsv")
    context_primary = [row for row in context if row["score_scheme"] == "weighted__primary"]
    pcls = {row["contrast_id"]: row for row in context_primary if row["dataset_id"] == "GSE200418"}
    tgfb = [
        row for row in context_primary
        if row["dataset_id"] == "GSE207889" and row["contrast_id"] == "TGFB_minus_matched_control"
    ]
    culture_exceeds_lipid = (
        number(pcls["culture_CTR_48h_minus_24h"]["estimate"])
        > number(pcls["GFIPO_48h_minus_GFI_48h"]["estimate"])
        and number(pcls["culture_CTR_48h_minus_24h"]["estimate"]) > 0
    )
    tgfb_replication = any(number(row["estimate"]) > 0 and int(row["n_positive"]) == 2 for row in tgfb)
    context_gate = culture_exceeds_lipid and tgfb_replication and human_gate

    spatial = read_tsv(CANDIDATE_ROOT / "spatial_lineage_support.tsv")
    class_spatial = next(row for row in spatial if row["analysis_id"] == "evidence_class_spatial_comparison")
    identity = read_tsv(CANDIDATE_ROOT / "program_identity_audit.tsv")
    spatial_gate = class_spatial["status"] == "pass" and all(truth(row["identity_support_pass"]) for row in identity)

    clcc = read_tsv(CANDIDATE_ROOT / "clcc1_integrity_correction/crispr_verdict.tsv")[0]
    full_gate = human_gate and crop_gate and protein_gate and spatial_gate
    chronic_reversible_partial = human_gate and context_gate and protein_gate and spatial_gate

    verdict_rows = [
        {
            "component_id": "human_histologic_reversal", "status": "fail", "gate_pass": human_gate,
            "evidence": (
                f"expected_direction={all_human_direction}; unseen_nominal={unseen_nominal}; "
                f"effective_cohorts={human_meta['n_effective_cohorts']}; aligned_meta_p={human_meta['p']}"
            ),
            "permitted_claim": "discovery-cohort reversal signal only; external replication failed",
        },
        {
            "component_id": "chronic_context_topology", "status": "fail", "gate_pass": context_gate,
            "evidence": (
                f"culture_effect={pcls['culture_CTR_48h_minus_24h']['estimate']}; "
                f"combined_lipid_effect={pcls['GFIPO_48h_minus_GFI_48h']['estimate']}; "
                f"TGFb_two_of_two_lineage={tgfb_replication}; human_gate={human_gate}"
            ),
            "permitted_claim": "TGF-beta directional context support; no chronic-state topology claim",
        },
        {
            "component_id": "regulatory_risk_cropseq_bridge", "status": crop["status"], "gate_pass": crop_gate,
            "evidence": (
                f"eligible_loci={crop['n_independent_corrected_tier12_loci']}; "
                f"orientation_increase={crop['n_risk_activity_increase']}; orientation_decrease={crop['n_risk_activity_decrease']}"
            ),
            "permitted_claim": "none; state-axis outcome was not inspected",
        },
        {
            "component_id": "protein_translation", "status": "fail", "gate_pass": protein_gate,
            "evidence": (
                f"PXD052787_q={protein_map['PXD052787']['q']}; "
                f"PXD052787_residual_p={residual['PXD052787']['p']}; "
                f"matching_pass={protein_map['PXD052787']['matching_gate_pass']}"
            ),
            "permitted_claim": "ordinary RNA-protein concordance in the frozen tissue cohort only",
        },
        {
            "component_id": "spatial_lineage_remodeling", "status": "fail_source_covariates", "gate_pass": spatial_gate,
            "evidence": class_spatial["status"],
            "permitted_claim": "accepted frozen M8/M20 spatial organization and identity audit only",
        },
        {
            "component_id": "clcc1_integrity_correction", "status": clcc["status"], "gate_pass": False,
            "evidence": "corrected 326-locus retrospective rerun; supplementary and nonconfirmatory",
            "permitted_claim": "valid null/error-correction boundary",
        },
        {
            "component_id": "complete_risk_to_state_architecture", "status": "complete_no_promotion", "gate_pass": full_gate,
            "evidence": (
                f"human={human_gate}; cropseq={crop_gate}; protein={protein_gate}; spatial={spatial_gate}; "
                f"chronic_partial={chronic_reversible_partial}"
            ),
            "permitted_claim": "none beyond pre-existing Resource thesis",
        },
    ]
    write_tsv(
        CANDIDATE_ROOT / "promotion_verdict.tsv",
        verdict_rows,
        ["component_id", "status", "gate_pass", "evidence", "permitted_claim"],
    )

    source_rows = []
    for row in read_tsv(CANDIDATE_ROOT / "source_gates/human_pair_gates_v3.tsv"):
        source_rows.append({
            "source_id": row["dataset_id"], "gate_id": row["gate_id"], "status": row["status"],
            "inference_authorized": row["inference_authorized"], "biological_units": row["n_primary_complete_units"],
            "detail": row["detail"],
        })
    source_rows.append({
        "source_id": "GSE106737+GSE83452", "gate_id": "PARTICIPANT_INDEPENDENCE",
        "status": "fail_confirmed_overlap", "inference_authorized": "false", "biological_units": "1_effective_cohort",
        "detail": "78/111 reciprocal expression fingerprints; all 20 lifestyle participants overlap",
    })
    source_rows.append({
        "source_id": "GSE281160", "gate_id": "CROPSEQ_ORIENTATION", "status": crop["status"],
        "inference_authorized": crop["inference_authorized"], "biological_units": crop["n_independent_corrected_tier12_loci"],
        "detail": crop["detail"],
    })
    for row in read_tsv(CANDIDATE_ROOT / "source_gates/protein_source_gates.tsv"):
        source_rows.append({
            "source_id": row["dataset_id"], "gate_id": row["gate_id"], "status": row["status"],
            "inference_authorized": row["inference_authorized"], "biological_units": row["n_biological_samples"],
            "detail": row["detail"],
        })
    source_rows.extend([
        {
            "source_id": "accepted_spatial_evidence", "gate_id": "PROGRAM_LEVEL_IMPORT", "status": "pass_frozen_M8_M20_only",
            "inference_authorized": "true", "biological_units": "donor_collapsed_source_specific",
            "detail": "four accepted residual spatial-organization rows imported without reinterpretation",
        },
        {
            "source_id": "accepted_spatial_evidence", "gate_id": "EVIDENCE_CLASS_COMPARISON",
            "status": class_spatial["status"], "inference_authorized": "false", "biological_units": "",
            "detail": "prespecified gene-level spatial variability and coverage covariates unavailable",
        },
    ])
    write_tsv(
        CANDIDATE_ROOT / "source_gate_status.tsv",
        source_rows,
        ["source_id", "gate_id", "status", "inference_authorized", "biological_units", "detail"],
    )

    source_manifest_rows = []
    for lane in ("human", "crop", "protein"):
        lane_rows = read_tsv(CANDIDATE_ROOT / f"acquisition/{lane}_source_manifest.tsv")
        source_manifest_rows.extend(lane_rows)
    write_tsv(
        CANDIDATE_ROOT / "source_manifest.tsv",
        source_manifest_rows,
        list(source_manifest_rows[0]),
    )

    correction = read_tsv(CANDIDATE_ROOT / "integrity_correction_v2/integrity_correction_manifest.tsv")
    write_tsv(CANDIDATE_ROOT / "integrity_correction_manifest.tsv", correction, list(correction[0]))

    scripts = sorted(path for path in SCRIPT_ROOT.iterdir() if path.is_file() and path.suffix in {".py", ".R", ".cpp", ".sbatch"})
    execution_rows = [{
        "artifact": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(path), "status": "executed_or_release_code",
    } for path in scripts]
    write_tsv(CANDIDATE_ROOT / "execution_manifest.tsv", execution_rows, ["artifact", "sha256", "status"])
    environment_rows = [
        {"component": "python", "version": platform.python_version(), "detail": sys.executable},
        {"component": "platform", "version": platform.platform(), "detail": ""},
        {"component": "seal", "version": seal["specification_sha256"], "detail": seal["sealed_at_utc"]},
    ]
    write_tsv(CANDIDATE_ROOT / "environment_manifest.tsv", environment_rows, ["component", "version", "detail"])

    summary = {
        "candidate_id": seal["candidate_id"],
        "adjudicated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "complete_no_promotion",
        "human_reversal_gate_pass": human_gate,
        "context_topology_gate_pass": context_gate,
        "cropseq_bridge_gate_pass": crop_gate,
        "protein_translation_gate_pass": protein_gate,
        "spatial_lineage_gate_pass": spatial_gate,
        "complete_architecture_gate_pass": full_gate,
        "figure5_promotion": False,
        "canonical_promotion_authorized": False,
        "journal_posture": "Cell_Genomics_Resource",
    }
    (CANDIDATE_ROOT / "ADJUDICATED.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
