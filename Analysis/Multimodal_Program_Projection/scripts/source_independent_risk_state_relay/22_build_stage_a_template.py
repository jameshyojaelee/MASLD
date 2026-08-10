#!/usr/bin/env python3
"""Seal the target-independent Plan 45 Stage A analysis template."""

from __future__ import annotations

import csv
import json
import os
import shutil
import tempfile
from collections import Counter
from pathlib import Path

from relay_common import PROJECT_ROOT, sha256_file


CANDIDATE_ID = os.environ.get(
    "PLAN45_STAGE_A_TEMPLATE_ID",
    "source-independent-risk-state-relay-stage-a-template-v2-2026-08-10",
).strip()
CANDIDATE_PARENT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates"
ROOT = CANDIDATE_PARENT / CANDIDATE_ID
REFERENCE = (
    CANDIDATE_PARENT
    / "source-independent-risk-state-relay-lineage-disease-reference-2026-08-10"
)
CV_CORRECTION = (
    CANDIDATE_PARENT
    / "source-independent-risk-state-relay-lineage-disease-reference-cv-fdr-correction-2026-08-10"
)
REFERENCE_TABLE = REFERENCE / "frozen_lineage_disease_reference.tsv"
REFERENCE_SEAL = REFERENCE / "LINEAGE_DISEASE_REFERENCE_SEALED.json"
CV_GATE = CV_CORRECTION / "cv_gate_status_multiplicity_corrected.tsv"
CV_SEAL = CV_CORRECTION / "LINEAGE_DISEASE_CV_FDR_CORRECTION_SEALED.json"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, data: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def main() -> None:
    if ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite Stage A template: {ROOT}")
    for path in (REFERENCE_TABLE, REFERENCE_SEAL, CV_GATE, CV_SEAL):
        if not path.is_file():
            raise RuntimeError(f"Stage A template prerequisite absent: {path}")
    reference_seal = json.loads(REFERENCE_SEAL.read_text(encoding="utf-8"))
    cv_seal = json.loads(CV_SEAL.read_text(encoding="utf-8"))
    if reference_seal.get("experimental_targets_frozen") is not False:
        raise RuntimeError("Reference unexpectedly freezes a target")
    if cv_seal.get("experimental_targets_frozen") is not False:
        raise RuntimeError("CV correction unexpectedly freezes a target")

    counts = Counter()
    within_counts = Counter()
    l1 = Counter()
    with REFERENCE_TABLE.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            lineage = row["cell_type"]
            if row["primary_reference_eligible"].lower() == "true":
                counts[lineage] += 1
                l1[lineage] += abs(float(row["primary_raw_loading"]))
            if row["within_dataset_direction_sensitivity_pass"].lower() == "true":
                within_counts[lineage] += 1

    axis_registry = []
    for row in read_rows(CV_GATE):
        lineage = row["cell_type"]
        eligible = row["multiplicity_corrected_cv_gate_pass"].lower() == "true"
        axis_registry.append(
            {
                "recipient_lineage": lineage,
                "primary_recipient_eligible": str(eligible).upper(),
                "heldout_gse244832_effect": row["gse244832_effect"],
                "heldout_gse244832_lineage_family_q": row["gse244832_lineage_family_q"],
                "heldout_gse174748_direction_effect": row["gse174748_effect"],
                "n_primary_reference_genes": counts[lineage],
                "n_within_dataset_direction_genes": within_counts[lineage],
                "primary_reference_l1": f"{l1[lineage]:.17g}",
                "primary_loading_column": "primary_raw_loading",
                "sign_sensitivity_column": "sign_only_raw_loading",
                "primary_use": (
                    "prespecified_recipient_endpoint"
                    if eligible
                    else "secondary_nonconfirmatory_only"
                ),
                "cell_autonomy_interpretation": "prohibited_without_substate_controls",
            }
        )
    axis_registry.sort(key=lambda row: row["recipient_lineage"])

    metadata_schema = [
        ("sample_id", "string", "always", "unique assay-level sample identifier"),
        ("biological_unit_id", "string", "always", "background x independent differentiation; inferential unit"),
        ("background_id", "string", "always", "at least three independent iPSC backgrounds"),
        ("clone_id", "string", "when exact edit is used", "clone nested within background"),
        ("differentiation_id", "string", "always", "at least two independent differentiations/background/condition"),
        ("batch_id", "string", "always", "randomization and batch audit"),
        ("well_id", "string", "always", "technical unit; never biological replicate"),
        ("assay", "enum", "always", "RNA; secretome; imaging; target_expression; viability"),
        ("target_locus_uid", "string", "after target freeze", "must match future immutable target registry"),
        ("target_gene", "string", "after target freeze", "dynamically removed from every state score"),
        ("guide_id", "string", "for CRISPRa/i arms", "two nonoverlapping guides plus controls"),
        ("guide_role", "enum", "for CRISPRa/i arms", "targeting_1; targeting_2; non_targeting; rescue"),
        ("guide_sequence_sha256", "sha256", "for every guide", "sequence provenance without exposing sequence in analysis table"),
        ("risk_direction", "enum", "always", "risk_mimic; control; protective_mimic"),
        ("source_lineage", "enum", "always", "hepatocyte primary; alternatives require independent accessibility replication"),
        ("edited_lineage", "enum", "always", "source; reciprocal_recipient; none"),
        ("recipient_lineage", "enum", "recipient assays", "Cholangiocytes; Endothelial_cells; Fibroblasts; Macrophages"),
        ("culture_context", "enum", "always", "source_only; direct_mosaic; transwell_recipient; conditioned_medium_recipient; reciprocal_recipient_edit"),
        ("time_role", "enum", "always", "early_cis; intermediate_mediator; late_relay"),
        ("time_value", "numeric", "always", "platform-native time with unit"),
        ("time_unit", "enum", "always", "hours; days"),
        ("challenge_id", "string", "always", "basal or one pilot-frozen chronic challenge"),
        ("challenge_role", "enum", "always", "basal; primary_chronic"),
        ("n_cells", "integer", "RNA recipient pseudobulk", "QC and fixed-cell-count sensitivity"),
        ("recipient_substate_fractions_path", "path", "recipient RNA", "composition-adjusted sensitivity"),
        ("viability_pass", "boolean", "always", "release gate, not rescue covariate"),
        ("maturation_pass", "boolean", "always", "release gate, not rescue covariate"),
        ("target_engagement_pass", "boolean", "targeting arms", "cis gate prerequisite"),
        ("randomization_block", "string", "always", "pre-unblinding allocation block"),
        ("blinded_label", "string", "always", "analysis-facing condition identifier"),
        ("unblinding_status", "enum", "always", "sealed; released_after_qc"),
        ("exclusion_reason", "string", "if excluded", "must be assigned before outcome unblinding"),
    ]
    metadata_rows = [
        {
            "field": field,
            "type": field_type,
            "required_when": required_when,
            "purpose": purpose,
        }
        for field, field_type, required_when, purpose in metadata_schema
    ]

    conditions = [
        ("SRC_CTRL", "source_only", "source", "control", "early_cis;intermediate_mediator;late_relay", "cis/metabolic/secretome control; no recipient score"),
        ("SRC_RISK", "source_only", "source", "risk_mimic", "early_cis;intermediate_mediator;late_relay", "cis/metabolic/secretome effect; no recipient score"),
        ("MOS_CTRL", "direct_mosaic", "source", "control", "early_cis;intermediate_mediator;late_relay", "matched recipient-containing control"),
        ("MOS_RISK", "direct_mosaic", "source", "risk_mimic", "early_cis;intermediate_mediator;late_relay", "primary source-to-recipient arm"),
        ("TRW_CTRL", "transwell_recipient", "source", "control", "intermediate_mediator;late_relay", "soluble-transfer control"),
        ("TRW_RISK", "transwell_recipient", "source", "risk_mimic", "intermediate_mediator;late_relay", "soluble-transfer test"),
        ("CM_CTRL", "conditioned_medium_recipient", "source", "control", "intermediate_mediator;late_relay", "alternative if transwell is infeasible"),
        ("CM_RISK", "conditioned_medium_recipient", "source", "risk_mimic", "intermediate_mediator;late_relay", "alternative if transwell is infeasible"),
        ("RECIP_CTRL", "reciprocal_recipient_edit", "reciprocal_recipient", "control", "early_cis;late_relay", "cell-of-origin control"),
        ("RECIP_RISK", "reciprocal_recipient_edit", "reciprocal_recipient", "risk_mimic", "early_cis;late_relay", "direct-recipient alternative"),
    ]
    condition_rows = [
        {
            "arm_id": arm_id,
            "culture_context": context,
            "edited_lineage_role": edited,
            "risk_direction": direction,
            "required_time_roles": times,
            "inferential_role": role,
            "recipient_score_permitted": str(context != "source_only").upper(),
        }
        for arm_id, context, edited, direction, times, role in conditions
    ]

    estimands = [
        {
            "estimand_id": "CIS01",
            "hierarchy": 1,
            "outcome": "source-lineage target expression",
            "contrast": "risk_mimic - control at early_cis in source_only and direct_mosaic",
            "family": "all Stage A loci x two guides",
            "pass_rule": "expected direction for both guides; family q<0.05; agreement in at least two backgrounds",
        },
        {
            "estimand_id": "RELAY01",
            "hierarchy": 2,
            "outcome": "validated recipient-lineage disease alignment",
            "contrast": "[(risk_mimic-control)_late_relay - (risk_mimic-control)_early_cis] within direct_mosaic",
            "family": "all Stage A loci x four eligible recipient axes",
            "pass_rule": "positive; family q<0.05; both guides; at least two backgrounds; cis gate already passed",
        },
        {
            "estimand_id": "ORIGIN01",
            "hierarchy": 3,
            "outcome": "recipient-lineage disease alignment",
            "contrast": "source-perturbed direct_mosaic relay effect - reciprocal_recipient_edit effect",
            "family": "only loci passing CIS01 and RELAY01",
            "pass_rule": "source effect larger in expected direction; family q<0.05",
        },
        {
            "estimand_id": "ROUTE01",
            "hierarchy": 4,
            "outcome": "recipient-lineage disease alignment",
            "contrast": "direct_mosaic versus transwell or conditioned_medium risk-control effect at late_relay",
            "family": "complete passing locus x recipient family",
            "pass_rule": "classify soluble versus contact-dependent route; no forced positive direction",
        },
        {
            "estimand_id": "COMP01",
            "hierarchy": 5,
            "outcome": "RELAY01 after substate covariates and fixed-cell-count pseudobulk",
            "contrast": "same as RELAY01",
            "family": "same as RELAY01",
            "pass_rule": "direction retained in both sensitivities; otherwise composition-mediated only",
        },
        {
            "estimand_id": "PHENO01",
            "hierarchy": 6,
            "outcome": "two pilot-frozen non-RNA phenotypes",
            "contrast": "risk_mimic - control at late_relay",
            "family": "all selected phenotypes x advancing loci",
            "pass_rule": "expected direction in at least two phenotypes; family corrected",
        },
    ]

    score_contract = [
        {"rule": "biological_unit", "value": "background x independent differentiation; cells/wells/lanes are technical"},
        {"rule": "normalization", "value": "RNA-count-native TMM/logCPM or negative-binomial pseudobulk; no cell-level inference"},
        {"rule": "standardization_reference", "value": "control-source early direct-mosaic recipient pseudobulks, frozen before unblinding"},
        {"rule": "primary_formula", "value": "sum(primary_raw_loading * standardized_expression) / sum(abs(retained_primary_raw_loading))"},
        {"rule": "recipient_axes", "value": "Cholangiocytes;Endothelial_cells;Fibroblasts;Macrophages"},
        {"rule": "dynamic_target_exclusion", "value": "remove perturbed target from every axis and every secondary program containing it"},
        {"rule": "guide_response_exclusion", "value": "freeze genes changed by non-targeting guide versus untransduced early QC before target-label unblinding; remove from all axes"},
        {"rule": "coverage_gate", "value": "at least 500 retained genes and at least 80% original axis L1 loading"},
        {"rule": "weight_sensitivities", "value": "sign-only; top-50%-absolute; top-25%-absolute; within-dataset-direction subset"},
        {"rule": "composition_sensitivities", "value": "recipient-substate covariate adjustment; fixed-cell-count pseudobulk subsampling"},
        {"rule": "secondary_program_family", "value": "all 117 frozen programs with complete BH family; none may replace RELAY01"},
    ]

    randomization = [
        {"item": "allocation", "contract": "randomize perturbation/context/time across plate and batch within background/differentiation"},
        {"item": "blinding", "contract": "analysis receives blinded labels until viability, maturation, count, and target-independent QC gates pass"},
        {"item": "exclusions", "contract": "record before unblinding; no outcome-based sample removal"},
        {"item": "guide_handling", "contract": "analyze both guides separately and jointly; one guide cannot carry a gate"},
        {"item": "missingness", "contract": "no imputation of missing biological units; report complete design cells"},
        {"item": "power", "contract": "blinded pilot variance locks Stage B sample size at family-wise alpha; Stage A is nomination only"},
    ]

    template_status = [
        {
            "template_status": "sealed_target_independent_analysis_template",
            "target_freeze_status": "prohibited_pending_registry_orientation_accessibility_guideability",
            "experimental_outcomes_inspected": "FALSE",
            "wetlab_samples_exist": "FALSE",
            "authoritative_recipient_axes": 4,
            "source_lineage_priority": "hepatocyte",
            "stage_a_finalization_required": "TRUE",
        }
    ]

    CANDIDATE_PARENT.mkdir(parents=True, exist_ok=True)
    temp_root = Path(tempfile.mkdtemp(prefix=f".{CANDIDATE_ID}.", dir=CANDIDATE_PARENT))
    try:
        objects = {
            "stage_a_recipient_axis_registry": axis_registry,
            "stage_a_required_metadata_schema": metadata_rows,
            "stage_a_condition_template": condition_rows,
            "stage_a_estimand_registry": estimands,
            "stage_a_score_contract": score_contract,
            "stage_a_randomization_blinding_contract": randomization,
            "stage_a_template_status": template_status,
        }
        output_paths = {}
        for name, data in objects.items():
            path = temp_root / f"{name}.tsv"
            write_rows(path, data, list(data[0]))
            output_paths[name] = path
        manifest = [
            {
                "source_id": source_id,
                "source_path": str(path.relative_to(PROJECT_ROOT)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for source_id, path in (
                ("lineage_disease_reference", REFERENCE_TABLE),
                ("lineage_disease_reference_seal", REFERENCE_SEAL),
                ("lineage_cv_corrected_gate", CV_GATE),
                ("lineage_cv_correction_seal", CV_SEAL),
            )
        ]
        manifest_path = temp_root / "stage_a_template_input_manifest.tsv"
        write_rows(manifest_path, manifest, list(manifest[0]))
        output_paths["stage_a_template_input_manifest"] = manifest_path
        seal = {
            "status": "sealed_target_independent_stage_a_template",
            "candidate_id": CANDIDATE_ID,
            "target_freeze_status": "prohibited",
            "experimental_outcomes_inspected": False,
            "wetlab_samples_exist": False,
            "primary_source_lineage": "hepatocyte",
            "primary_recipient_lineages": [
                "Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Macrophages"
            ],
            "primary_recipient_estimand": "RELAY01",
            "invalid_source_monoculture_recipient_comparison_prohibited": True,
            "finalization_dependencies": [
                "registry_v3", "shared_signal_orientation", "lineage_accessibility",
                "guideability", "target_registry", "blinded_pilot"
            ],
            "output_sha256": {
                name: sha256_file(path) for name, path in output_paths.items()
            },
        }
        (temp_root / "STAGE_A_TEMPLATE_SEALED.json").write_text(
            json.dumps(seal, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        temp_root.rename(ROOT)
    except Exception:
        if temp_root.exists():
            shutil.rmtree(temp_root)
        raise
    print(
        "Plan 45 Stage A template sealed: recipient_axes=4; "
        "estimands=6; target freeze prohibited"
    )


if __name__ == "__main__":
    main()
