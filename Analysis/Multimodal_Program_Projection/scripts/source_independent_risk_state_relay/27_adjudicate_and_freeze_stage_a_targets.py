#!/usr/bin/env python3
"""Apply the sealed outcome-blind policy and freeze the actual Stage-A design.

The only inputs are sealed genetics/routing worklists, external guide-design QC,
a target-independent platform pilot, and the previously sealed recipient-axis
template. No cis, recipient-state, mediator, or phenotype outcome is accepted.
"""

from __future__ import annotations

import datetime as dt
import json
import os
from collections import defaultdict
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_json,
    read_tsv,
    sha256_file,
    write_tsv,
)


DIRECT_STRATA = {
    "direct_masld_mash_diagnosis",
    "mri_pdff_or_histologic_steatosis",
}
BOOLEAN_PROTOCOL_GATES = [
    "pilot_blinded", "viability_gate_pass", "maturation_gate_pass",
    "source_lineage_gate_pass", "lineage_mix_gate_pass",
    "chronic_challenge_gate_pass", "early_intermediate_late_timepoints_locked",
    "randomization_locked", "blinding_locked", "missingness_locked",
    "power_inputs_locked",
]


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def number(value: object) -> float:
    result = float(value)
    if result != result or abs(result) == float("inf"):
        raise RuntimeError(f"Non-finite target-ranking value: {value!r}")
    return result


def candidate_source(variable: str) -> Path:
    raw = os.environ.get(variable, "").strip()
    if not raw:
        raise RuntimeError(f"{variable} is required")
    path = Path(raw)
    path = path if path.is_absolute() else PROJECT_ROOT / path
    path = path.resolve()
    allowed = (
        PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates"
    ).resolve()
    if allowed not in path.parents:
        raise RuntimeError(f"{variable} escapes the candidate root: {path}")
    return path


def stage_a_guide_pass(row: dict[str, str]) -> bool:
    return (
        row["audit_row_type"] == "guide"
        and row["design_status"] == "design_available"
        and yes(row["source_recommended"])
        and yes(row["sequence_qc_pass"])
        and yes(row["off_target_review_pass"])
        and yes(row["review_concordant"])
    )


def stage_b_design_pass(row: dict[str, str]) -> bool:
    return (
        row["audit_row_type"] == "design"
        and row["design_status"] == "design_available"
        and yes(row["source_recommended"])
        and yes(row["complete_design"])
        and yes(row["intended_edit_matches_worklist"])
        and yes(row["off_target_review_pass"])
        and yes(row["review_concordant"])
    )


def select_distinct(
    rows: list[dict[str, str]], id_field: str, sequence_fields: list[str], count: int
) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    identifiers: set[str] = set()
    sequences: set[tuple[str, ...]] = set()
    for row in sorted(
        rows, key=lambda value: (int(value["source_design_rank"]), value[id_field])
    ):
        signature = tuple(row[field].upper() for field in sequence_fields)
        if row[id_field] in identifiers or signature in sequences:
            continue
        selected.append(row)
        identifiers.add(row[id_field])
        sequences.add(signature)
        if len(selected) == count:
            break
    return selected


def rank_key(candidate: dict[str, object]) -> tuple[object, ...]:
    return (
        -int(candidate["n_concordant_gwas_evidence_families"]),
        -float(candidate["best_orientation_consensus"]),
        -float(candidate["best_susie_pp4"]),
        -float(candidate["selected_exact_shared_posterior"]),
        str(candidate["coarse_locus_uid"]),
        str(candidate["ensembl_id"]),
    )


def choose_targets(
    candidates: list[dict[str, object]], maximum: int = 6
) -> tuple[str, list[dict[str, object]]]:
    """Choose a balanced screen when possible, otherwise one locus only."""

    eligible = [row for row in candidates if bool(row["candidate_gate_pass"])]
    # One physical locus is one inferential unit even if a registry defect emits
    # two representative genes. Keep the deterministic top candidate only.
    by_locus: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in eligible:
        by_locus[str(row["coarse_locus_uid"])].append(row)
    pool = [sorted(rows, key=rank_key)[0] for rows in by_locus.values()]
    pool = sorted(pool, key=rank_key)
    if not pool:
        return "no_target_freeze", []

    directions = {str(row["oriented_risk_effect"]) for row in pool}
    strata = {
        stratum
        for row in pool
        for stratum in row["phenotype_strata_set"]  # type: ignore[index]
    }
    screen_eligible = len(pool) >= 4 and len(directions) >= 2 and DIRECT_STRATA <= strata
    if not screen_eligible:
        return "single_locus_mechanism_only", [pool[0]]

    selected: list[dict[str, object]] = []
    uncovered_directions = set(directions)
    uncovered_strata = set(DIRECT_STRATA)
    remaining = list(pool)
    while remaining and len(selected) < min(maximum, len(pool)):
        ordered = sorted(
            remaining,
            key=lambda row: (
                -(
                    int(str(row["oriented_risk_effect"]) in uncovered_directions)
                    + len(set(row["phenotype_strata_set"]) & uncovered_strata)  # type: ignore[arg-type]
                ),
                *rank_key(row),
            ),
        )
        chosen = ordered[0]
        selected.append(chosen)
        remaining.remove(chosen)
        uncovered_directions.discard(str(chosen["oriented_risk_effect"]))
        uncovered_strata -= set(chosen["phenotype_strata_set"])  # type: ignore[arg-type]
    return "balanced_four_to_six_locus_screen", selected


def validate_seal_outputs(root: Path, seal: dict[str, object]) -> None:
    for name, expected in seal["output_sha256"].items():  # type: ignore[union-attr]
        path = root / str(name)
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Sealed input hash mismatch: {path}")


def validate_stage_a_template_outputs(root: Path, seal: dict[str, object]) -> None:
    for stem, expected in seal["output_sha256"].items():  # type: ignore[union-attr]
        path = root / f"{stem}.tsv"
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Stage-A template output hash mismatch: {path}")


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite Stage-A target freeze: {CANDIDATE_ROOT}")
    guide_root = candidate_source("PLAN45_GUIDE_ROOT")
    audit_root = candidate_source("PLAN45_GUIDE_AUDIT_ROOT")
    contract_root = candidate_source("PLAN45_ADJUDICATION_CONTRACT_ROOT")
    template_root = candidate_source("PLAN45_STAGE_A_TEMPLATE_ROOT")

    guide_seal_path = guide_root / "GUIDE_DESIGN_WORKLIST_SEALED.json"
    audit_seal_path = audit_root / "EXTERNAL_GUIDE_AUDIT_SEALED.json"
    contract_seal_path = contract_root / "TARGET_ADJUDICATION_CONTRACT_SEALED.json"
    template_seal_path = template_root / "STAGE_A_TEMPLATE_SEALED.json"
    seals = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in [guide_seal_path, audit_seal_path, contract_seal_path, template_seal_path]
    ]
    guide_seal, audit_seal, contract_seal, template_seal = seals
    expected_status = [
        "separate_stage_a_promoter_and_stage_b_exact_designs_pending",
        "external_guide_and_target_independent_protocol_audits_sealed",
        "sealed_outcome_blind_target_adjudication_contract",
        "sealed_target_independent_stage_a_template",
    ]
    if [seal.get("status") for seal in seals] != expected_status:
        raise RuntimeError("One or more target-freeze dependencies have an invalid status")
    if any(seal.get("scientific_outcomes_inspected") is not False for seal in seals[:3]):
        raise RuntimeError("A target-freeze dependency reports scientific outcome access")
    if any(seal.get("experimental_targets_frozen") is not False for seal in seals[:3]):
        raise RuntimeError("An upstream dependency already froze targets")
    for root, seal in [
        (guide_root, guide_seal), (audit_root, audit_seal),
        (contract_root, contract_seal),
    ]:
        validate_seal_outputs(root, seal)
    validate_stage_a_template_outputs(template_root, template_seal)

    stage_a_work = read_tsv(guide_root / "stage_a_promoter_design_worklist.tsv")
    stage_b_work = read_tsv(guide_root / "stage_b_exact_edit_worklist.tsv")
    stage_a_audit = read_tsv(audit_root / "stage_a_crispick_design_audit.tsv")
    stage_b_audit = read_tsv(audit_root / "stage_b_primedesign_design_audit.tsv")
    protocol = read_tsv(audit_root / "target_independent_protocol_gate.tsv")

    passing_platforms = [
        row for row in protocol
        if not yes(row["scientific_outcomes_opened"])
        and int(row["n_independent_backgrounds"]) >= 3
        and int(row["n_independent_differentiations_per_background_condition"]) >= 2
        and all(yes(row[field]) for field in BOOLEAN_PROTOCOL_GATES)
    ]
    selected_platform = (
        sorted(passing_platforms, key=lambda row: row["platform_id"])[0]
        if passing_platforms else None
    )

    audits_a: dict[str, list[dict[str, str]]] = defaultdict(list)
    audits_b: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in stage_a_audit:
        audits_a[row["promoter_request_uid"]].append(row)
    for row in stage_b_audit:
        audits_b[row["design_uid"]].append(row)
    exact_by_orientation: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in stage_b_work:
        exact_by_orientation[row["orientation_uid"]].append(row)

    by_target: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in stage_a_work:
        by_target[(row["coarse_locus_uid"], row["ensembl_id"])].append(row)

    candidate_rows: list[dict[str, object]] = []
    technical: dict[str, dict[str, object]] = {}
    for (locus_uid, ensembl_id), routes in sorted(by_target.items()):
        target_uid = f"target::{locus_uid}::{ensembl_id}"
        directions = {row["oriented_risk_effect"] for row in routes}
        family_concordant = len(directions) == 1
        phenotype_strata = {row["phenotype_stratum"] for row in routes}
        options: list[dict[str, object]] = []
        for route in routes:
            valid_guides = select_distinct(
                [
                    row for row in audits_a[route["promoter_request_uid"]]
                    if stage_a_guide_pass(row)
                ],
                "guide_id", ["guide_sequence"], 2,
            )
            if len(valid_guides) < 2:
                continue
            for exact in exact_by_orientation[route["orientation_uid"]]:
                if not yes(exact["exact_edit_exportable"]):
                    continue
                valid_pegrnas = select_distinct(
                    [
                        row for row in audits_b[exact["design_uid"]]
                        if stage_b_design_pass(row)
                    ],
                    "pegrna_id",
                    ["spacer_sequence", "pbs_sequence", "rtt_sequence", "ngrna_sequence"],
                    2,
                )
                if len(valid_pegrnas) < 2:
                    continue
                options.append(
                    {
                        "route": route,
                        "exact": exact,
                        "guides": valid_guides,
                        "pegrnas": valid_pegrnas,
                    }
                )
        options.sort(
            key=lambda option: (
                -number(option["route"]["orientation_consensus"]),  # type: ignore[index]
                -number(option["route"]["susie_pp4"]),  # type: ignore[index]
                -number(option["exact"]["shared_posterior"]),  # type: ignore[index]
                option["route"]["gwas_name"],  # type: ignore[index]
                option["exact"]["design_uid"],  # type: ignore[index]
            )
        )
        best = options[0] if options else None
        architecture_eligible = all(
            row["dominant_celltype"] == "Hepatocytes"
            and row["accessibility_evidence_class"]
            == "two_cohort_hepatocyte_accessibility"
            for row in routes
        )
        gate_pass = bool(
            family_concordant and architecture_eligible and best and selected_platform
        )
        if not family_concordant:
            gate_reason = "failed_cross_gwas_orientation_concordance"
        elif not architecture_eligible:
            gate_reason = "failed_independently_replicated_hepatocyte_source_gate"
        elif best is None:
            gate_reason = "failed_stage_a_guide_or_stage_b_exact_editability"
        elif selected_platform is None:
            gate_reason = "failed_target_independent_platform_protocol_gate"
        else:
            gate_reason = "pass_pending_deterministic_screen_selection"

        best_route = best["route"] if best else sorted(routes, key=lambda row: row["orientation_uid"])[0]
        best_exact = best["exact"] if best else None
        candidate = {
            "target_uid": target_uid,
            "coarse_locus_uid": locus_uid,
            "gene_symbol": best_route["gene_symbol"],
            "ensembl_id": ensembl_id,
            "dominant_celltype": best_route["dominant_celltype"],
            "accessibility_evidence_class": best_route["accessibility_evidence_class"],
            "oriented_risk_effect": (
                next(iter(directions)) if family_concordant else "conflicting"
            ),
            "stage_a_perturbation_mode": best_route["stage_a_perturbation_mode"],
            "n_concordant_gwas": len({row["gwas_name"] for row in routes}),
            "n_concordant_gwas_evidence_families": len(
                {
                    row["gwas_evidence_family"]
                    for row in routes
                    if yes(row["counts_for_replication_breadth"])
                }
            ),
            "gwas_names": ";".join(sorted({row["gwas_name"] for row in routes})),
            "gwas_evidence_families": ";".join(
                sorted({row["gwas_evidence_family"] for row in routes})
            ),
            "gwas_independence_classes": ";".join(
                sorted({row["gwas_independence_class"] for row in routes})
            ),
            "phenotype_strata": ";".join(sorted(phenotype_strata)),
            "phenotype_strata_set": phenotype_strata,
            "cross_gwas_orientation_concordant": family_concordant,
            "best_orientation_uid": best_route["orientation_uid"],
            "best_orientation_consensus": number(best_route["orientation_consensus"]),
            "best_susie_pp4": number(best_route["susie_pp4"]),
            "selected_promoter_request_uid": (
                best_route["promoter_request_uid"] if best else ""
            ),
            "selected_exact_design_uid": best_exact["design_uid"] if best_exact else "",
            "selected_exact_snp_hg19": best_exact["snp_hg19"] if best_exact else "",
            "selected_exact_position_hg38": (
                best_exact["position_hg38"] if best_exact else ""
            ),
            "selected_exact_reference_allele": (
                best_exact["reference_allele_forward"] if best_exact else ""
            ),
            "selected_exact_alternate_allele": (
                best_exact["alternate_allele_forward"] if best_exact else ""
            ),
            "selected_exact_reference_to_alternate_role": (
                best_exact["reference_to_alternate_role"] if best_exact else ""
            ),
            "selected_exact_shared_posterior": (
                number(best_exact["shared_posterior"]) if best_exact else -1.0
            ),
            "candidate_gate_pass": gate_pass,
            "candidate_gate_reason": gate_reason,
            "selection_status": "not_selected",
            "selection_mode": "pending",
            "target_freeze_status": "pending_deterministic_adjudication",
        }
        candidate_rows.append(candidate)
        if best:
            technical[target_uid] = best

    selection_mode, selected = choose_targets(candidate_rows, maximum=6)
    selected_uids = {str(row["target_uid"]) for row in selected}
    for rank, row in enumerate(sorted(selected, key=rank_key), start=1):
        row["selection_rank"] = rank
    for row in candidate_rows:
        if row["target_uid"] in selected_uids:
            row["selection_status"] = "selected_for_stage_a"
            row["selection_mode"] = selection_mode
            row["target_freeze_status"] = "frozen_for_stage_a_only"
        else:
            row["selection_status"] = (
                "eligible_not_selected_by_frozen_priority"
                if row["candidate_gate_pass"] else "ineligible"
            )
            row["selection_mode"] = selection_mode
            row["target_freeze_status"] = "not_frozen"
        row["selection_rank"] = row.get("selection_rank", "")
        row["phenotype_strata_set"] = ""

    selected_guides: list[dict[str, object]] = []
    editability: list[dict[str, object]] = []
    for row in selected:
        target_uid = str(row["target_uid"])
        best = technical[target_uid]
        for guide in best["guides"]:  # type: ignore[index]
            selected_guides.append(
                {
                    "target_uid": target_uid,
                    "gene_symbol": row["gene_symbol"],
                    "orientation_uid": row["best_orientation_uid"],
                    "promoter_request_uid": row["selected_promoter_request_uid"],
                    "perturbation_mode": row["stage_a_perturbation_mode"],
                    "guide_id": guide["guide_id"],
                    "guide_sequence": guide["guide_sequence"],
                    "source_design_rank": guide["source_design_rank"],
                    "stage_a_status": "frozen_primary_guide",
                }
            )
        for pegrna in best["pegrnas"]:  # type: ignore[index]
            editability.append(
                {
                    "target_uid": target_uid,
                    "gene_symbol": row["gene_symbol"],
                    "orientation_uid": row["best_orientation_uid"],
                    "design_uid": row["selected_exact_design_uid"],
                    "snp_hg19": row["selected_exact_snp_hg19"],
                    "pegrna_id": pegrna["pegrna_id"],
                    "spacer_sequence": pegrna["spacer_sequence"],
                    "pbs_sequence": pegrna["pbs_sequence"],
                    "rtt_sequence": pegrna["rtt_sequence"],
                    "ngrna_sequence": pegrna["ngrna_sequence"],
                    "source_design_rank": pegrna["source_design_rank"],
                    "stage_b_status": "editability_passed_not_stage_b_frozen",
                }
            )

    conditions = read_tsv(template_root / "stage_a_condition_template.tsv")
    recipient_axes = [
        row for row in read_tsv(template_root / "stage_a_recipient_axis_registry.tsv")
        if yes(row["primary_recipient_eligible"])
    ]
    design_rows: list[dict[str, object]] = []
    guides_by_target: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in selected_guides:
        guides_by_target[str(row["target_uid"])].append(row)
    for target in selected:
        for guide in guides_by_target[str(target["target_uid"])]:
            for condition in conditions:
                axes = recipient_axes if yes(condition["recipient_score_permitted"]) else [None]
                for axis in axes:
                    design_rows.append(
                        {
                            "target_uid": target["target_uid"],
                            "gene_symbol": target["gene_symbol"],
                            "source_lineage": target["dominant_celltype"],
                            "oriented_risk_effect": target["oriented_risk_effect"],
                            "perturbation_mode": target["stage_a_perturbation_mode"],
                            "guide_id": guide["guide_id"],
                            "arm_id": condition["arm_id"],
                            "culture_context": condition["culture_context"],
                            "risk_direction": condition["risk_direction"],
                            "required_time_roles": condition["required_time_roles"],
                            "recipient_score_permitted": condition[
                                "recipient_score_permitted"
                            ],
                            "recipient_lineage": (
                                axis["recipient_lineage"] if axis else "not_applicable"
                            ),
                            "design_status": "frozen_before_scientific_outcome_access",
                        }
                    )

    CANDIDATE_ROOT.mkdir(parents=True)
    registry_path = CANDIDATE_ROOT / "frozen_target_registry.tsv"
    guide_path = CANDIDATE_ROOT / "frozen_stage_a_guides.tsv"
    edit_path = CANDIDATE_ROOT / "stage_b_exact_editability_registry.tsv"
    design_path = CANDIDATE_ROOT / "stage_a_design.tsv"
    candidate_fields = [key for key in candidate_rows[0] if key != "phenotype_strata_set"]
    write_tsv(registry_path, candidate_rows, candidate_fields)
    write_tsv(
        guide_path, selected_guides,
        list(selected_guides[0]) if selected_guides else [
            "target_uid", "gene_symbol", "orientation_uid", "promoter_request_uid",
            "perturbation_mode", "guide_id", "guide_sequence", "source_design_rank",
            "stage_a_status",
        ],
    )
    write_tsv(
        edit_path, editability,
        list(editability[0]) if editability else [
            "target_uid", "gene_symbol", "orientation_uid", "design_uid", "snp_hg19",
            "pegrna_id", "spacer_sequence", "pbs_sequence", "rtt_sequence",
            "ngrna_sequence", "source_design_rank", "stage_b_status",
        ],
    )
    write_tsv(
        design_path, design_rows,
        list(design_rows[0]) if design_rows else [
            "target_uid", "gene_symbol", "source_lineage", "oriented_risk_effect",
            "perturbation_mode", "guide_id", "arm_id", "culture_context",
            "risk_direction", "required_time_roles", "recipient_score_permitted",
            "recipient_lineage", "design_status",
        ],
    )
    status_rows = [{
        "selection_mode": selection_mode,
        "n_candidate_targets": len(candidate_rows),
        "n_gate_pass": sum(bool(row["candidate_gate_pass"]) for row in candidate_rows),
        "n_stage_a_targets_frozen": len(selected),
        "selected_platform_id": selected_platform["platform_id"] if selected_platform else "",
        "architecture_generalization_permitted": str(
            selection_mode == "balanced_four_to_six_locus_screen"
        ).lower(),
        "stage_b_frozen": "false",
        "scientific_outcomes_inspected": "false",
    }]
    status_path = CANDIDATE_ROOT / "target_freeze_status.tsv"
    write_tsv(status_path, status_rows, list(status_rows[0]))

    manifest_rows = []
    for role, path in [
        ("guide_release", guide_seal_path),
        ("external_guide_audit", audit_seal_path),
        ("adjudication_contract", contract_seal_path),
        ("stage_a_template", template_seal_path),
    ]:
        manifest_rows.append({
            "role": role,
            "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    manifest_path = CANDIDATE_ROOT / "target_freeze_input_manifest.tsv"
    write_tsv(manifest_path, manifest_rows, ["role", "source_path", "size_bytes", "sha256"])

    release_status = {
        "balanced_four_to_six_locus_screen": "stage_a_screen_targets_frozen",
        "single_locus_mechanism_only": (
            "single_locus_stage_a_target_frozen_architecture_generalization_prohibited"
        ),
        "no_target_freeze": "no_stage_a_target_frozen",
    }[selection_mode]
    payload = {
        "status": release_status,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "selection_mode": selection_mode,
        "n_candidate_targets": len(candidate_rows),
        "n_gate_pass": sum(bool(row["candidate_gate_pass"]) for row in candidate_rows),
        "n_stage_a_targets_frozen": len(selected),
        "frozen_target_uids": sorted(selected_uids),
        "scientific_outcomes_inspected": False,
        "experimental_targets_frozen": bool(selected),
        "stage_b_exact_editability_preverified": bool(editability),
        "stage_b_design_frozen": False,
        "architecture_generalization_permitted": (
            selection_mode == "balanced_four_to_six_locus_screen"
        ),
        "next_gate": (
            "execute the complete blinded Stage-A screen and adjudicate CIS01/RELAY01"
            if selected else "stop: no target passed the frozen source/guide/protocol gate"
        ),
        "output_sha256": {
            path.name: sha256_file(path) for path in [
                registry_path, guide_path, edit_path, design_path, status_path, manifest_path
            ]
        },
    }
    atomic_write_json(CANDIDATE_ROOT / "FROZEN_STAGE_A_TARGETS.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
