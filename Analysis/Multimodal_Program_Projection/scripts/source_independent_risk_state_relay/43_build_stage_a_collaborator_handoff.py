#!/usr/bin/env python3
"""Build a portable, target-independent Plan 45 experimental handoff."""

from __future__ import annotations

import datetime as dt
import json
import os
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_json,
    atomic_write_text, sha256_file, write_tsv,
)


def candidate_source(variable: str) -> Path:
    raw = os.environ.get(variable, "").strip()
    if not raw:
        raise RuntimeError(f"{variable} is required")
    path = Path(raw)
    path = path if path.is_absolute() else PROJECT_ROOT / path
    path = path.resolve()
    allowed = (PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates").resolve()
    if allowed not in path.parents:
        raise RuntimeError(f"{variable} escapes the candidate root: {path}")
    return path


CONTRACTS = {
    "stage_a_template": (
        "PLAN45_STAGE_A_TEMPLATE_ROOT", "STAGE_A_TEMPLATE_SEALED.json",
        "sealed_target_independent_stage_a_template",
    ),
    "target_adjudication": (
        "PLAN45_TARGET_ADJUDICATION_CONTRACT_ROOT",
        "TARGET_ADJUDICATION_CONTRACT_SEALED.json",
        "sealed_outcome_blind_target_adjudication_contract",
    ),
    "execution": (
        "PLAN45_STAGE_A_EXECUTION_CONTRACT_ROOT",
        "STAGE_A_EXECUTION_CONTRACT_SEALED.json",
        "sealed_outcome_blind_stage_a_execution_contract",
    ),
    "qc_unblinding": (
        "PLAN45_STAGE_A_QC_CONTRACT_ROOT",
        "STAGE_A_QC_UNBLINDING_CONTRACT_SEALED.json",
        "sealed_outcome_blind_stage_a_qc_unblinding_contract",
    ),
    "terminal": (
        "PLAN45_STAGE_A_TERMINAL_CONTRACT_ROOT",
        "STAGE_A_TERMINAL_CONTRACT_SEALED.json",
        "sealed_outcome_blind_stage_a_terminal_and_stage_b_freeze_contract",
    ),
}


QUESTIONS = [
    ("FQ01", "platform", "Can the platform generate a multicellular human liver organoid or microphysiological system with separately introduced source and recipient lineages?", "blocking"),
    ("FQ02", "replication", "Can at least three genetically independent iPSC backgrounds be studied?", "blocking"),
    ("FQ03", "replication", "Can at least two independent differentiations per background and condition be produced?", "blocking"),
    ("FQ04", "lineages", "Can hepatocyte source cells and cholangiocyte, endothelial, fibroblast or stellate, and macrophage recipient compartments be measured separately?", "blocking"),
    ("FQ05", "mosaic", "Can source-restricted and reciprocal-lineage perturbation mosaics be assembled while recipient lineages remain isogenic controls?", "blocking"),
    ("FQ06", "exposure", "Can one continuous or repeatedly dosed chronic metabolic challenge be maintained through a 10-to-21-day late readout?", "blocking"),
    ("FQ07", "route", "Can a direct mosaic and either transwell or conditioned-medium transfer be run from the same source differentiation?", "blocking"),
    ("FQ08", "assays", "Can target expression, lineage-resolved RNA, secretome, composition, viability or maturation, and two biologically distinct non-RNA phenotypes be measured?", "blocking"),
    ("FQ09", "units", "Can cells, wells, organoids, and lanes remain technical units while background by independent differentiation is preserved as the inferential unit?", "blocking"),
    ("FQ10", "blinding", "Can a data manager who is not the analysis lead hold the randomization and condition key through blinded QC?", "blocking"),
    ("FQ11", "stage_b", "If Stage A passes, can exact endogenous prime editing use at least two independent edited clones per allele and background?", "stage_b_blocking"),
    ("FQ12", "stage_b", "If Stage A passes, can target-dosage rescue and a mediator-blockade arm be implemented without changing edit identity?", "stage_b_blocking"),
]


DELIVERABLES = [
    ("D01", "platform feasibility answers", "platform_lead", "before target disclosure", "platform_feasibility_questions.tsv"),
    ("D02", "target-independent platform protocol gate", "wetlab_lead", "before target freeze", "contracts/target_adjudication/target_independent_protocol_gate_template.tsv"),
    ("D03", "blinded pilot and power audit", "wetlab_lead", "before Stage-A randomization", "contracts/execution/pilot_power_audit_template.tsv"),
    ("D04", "biological-unit manifest", "wetlab_lead", "before Stage-A randomization", "contracts/execution/experimental_unit_manifest_template.tsv"),
    ("D05", "assay registry", "wetlab_and_analysis_leads", "before Stage-A randomization", "contracts/execution/assay_registry_template.tsv"),
    ("D06", "two-phenotype registry", "wetlab_and_analysis_leads", "before Stage-A randomization", "contracts/execution/phenotype_registry_template.tsv"),
    ("D07", "complete randomized sample manifest", "data_manager", "before outcome generation", "contracts/execution/randomization_manifest_template.tsv"),
    ("D08", "source-file and dual-review manifests", "data_manager", "before outcome generation", "contracts/execution/execution_source_manifest_template.tsv;contracts/execution/dual_review_signoff_template.tsv"),
    ("D09", "pre-QC assay-threshold freeze", "qc_reviewers", "before raw QC inspection", "contracts/qc_unblinding/blinded_assay_qc_thresholds_template.tsv;contracts/qc_unblinding/threshold_freeze_signoff_template.tsv"),
    ("D10", "blinded files, QC metrics, and guide-response family", "qc_reviewers", "before condition-key release", "contracts/qc_unblinding/blinded_sample_file_manifest_template.tsv;contracts/qc_unblinding/blinded_sample_qc_metrics_template.tsv;contracts/qc_unblinding/target_independent_guide_response_template.tsv"),
    ("D11", "condition key and chronological signoff", "data_manager", "after blinded-QC seal", "contracts/qc_unblinding/unblinding_key_template.tsv;contracts/qc_unblinding/unblinding_signoff_template.tsv"),
    ("D12", "complete Stage-A estimand and sensitivity families", "analysis_lead", "after authorized unblinding", "contracts/terminal/stage_a_primary_estimand_results_template.tsv;contracts/terminal/stage_a_sensitivity_results_template.tsv"),
    ("D13", "complete mediator-discovery family and outcome signoff", "analysis_lead", "after authorized unblinding", "contracts/terminal/stage_a_mediator_candidates_template.tsv;contracts/terminal/stage_a_outcome_signoff_template.tsv"),
    ("D14", "prospective exact-edit Stage-B freeze if eligible", "wetlab_and_analysis_leads", "after Stage-A adjudication and before Stage-B outcomes", "contracts/terminal/stage_b_exact_design_template.tsv;contracts/terminal/stage_b_power_units_template.tsv;contracts/terminal/stage_b_randomization_template.tsv;contracts/terminal/stage_b_freeze_signoff_template.tsv"),
]


ROLES = [
    ("platform_lead", "platform capability and feasibility answers", "no", "no", "cannot select a target or endpoint"),
    ("wetlab_lead", "pilot, culture, perturbation, assay execution, source files", "after immutable target freeze", "only assay QC before handoff", "cannot exclude samples after outcome access"),
    ("data_manager", "randomization, opaque labels, condition key, file checksums", "yes after target freeze", "no scientific outcome interpretation", "must be distinct from analysis lead"),
    ("qc_reviewer_1", "threshold and blinded technical-QC review", "no condition key", "technical QC only", "must be distinct from reviewer 2"),
    ("qc_reviewer_2", "independent threshold and blinded-QC review", "no condition key", "technical QC only", "must be distinct from reviewer 1"),
    ("analysis_lead", "prespecified assay-native inference and terminal adjudication", "only after QC seal", "yes after authorized unblinding", "cannot hold condition key during blinded QC"),
    ("independent_validator", "hash, family, biological-unit, and gate rederivation", "only released manifests", "yes after immutable release", "cannot change targets, contrasts, or exclusions"),
]


HANDOFF = """# Stage-A source-independent MASLD risk-to-state experiment: collaborator handoff

## Ultimate goal

The goal is to determine whether an exactly oriented MASLD regulatory-risk operation first changes its nominated target in a defined source lineage and then produces a delayed, rescuable response in unedited recipient lineages, with concordant non-RNA tissue phenotypes. A passing experiment would provide the source-independent perturbational bridge needed to explain why inherited regulatory susceptibility and the established multicellular MASLD state nominate different genes. A failed cis, relay, composition, rescue, blockade, or phenotype gate is a valid biological boundary and must be reported without substituting another target.

## What this package is

This is a target-independent feasibility and data-governance package. It contains no target, guide, condition allocation, biological outcome, or authorization to start Stage B. Its purpose is to let an experimental collaborator decide whether the platform can execute the frozen design before target disclosure.

## Experimental spine

1. A blinded target-independent pilot freezes one chronic challenge, early/intermediate/late timing, viability, maturation, composition, and paired variance.
2. After the genetic, accessibility, guideability, and exact-edit gates select a target, Stage A uses two promoter guides in the genetically predicted target-expression direction across at least three backgrounds and two independent differentiations per background and condition.
3. Source-only cultures establish cis engagement. Direct mosaics test delayed responses in unedited cholangiocyte, endothelial, fibroblast or stellate, and macrophage recipients. Reciprocal-lineage perturbation tests cell of origin; transwell or conditioned medium classifies transfer route.
4. Blinded QC freezes exclusions before the condition key opens. Cells, wells, organoids, and lanes never become biological replicates.
5. Stage A can nominate at most two loci. Stage B remains prohibited unless a prespecified mediator, blocker, exact allele design, target rescue, two phenotypes, biological units, and randomization are frozen.

## First collaborator action

Treat this sealed directory as read-only and follow `RESPONSE_INSTRUCTIONS.md`. Complete copies of `platform_feasibility_questions.tsv` and `role_and_blinding_firewall.tsv` in a separately identified response package without requesting a target identity. Any blocking `no` or unresolved answer stops experimental planning and leaves the paper in its current Resource posture. If the platform passes, return the target-independent protocol-gate and pilot plan using copies of the immutable templates.

## Non-negotiable interpretation boundaries

- Stage A models risk-oriented target dosage; it does not validate the endogenous regulatory allele.
- Stage B exact endogenous editing is required for a locus mechanism.
- A recipient response without time ordering, source-lineage restriction, composition controls, target rescue, mediator blockade, and two non-RNA phenotypes is supportive, not mechanistic.
- No secondary program, post-outcome mediator, easier replacement locus, or favorable sample exclusion can rescue a failed primary gate.
"""


RESPONSE_INSTRUCTIONS = """# Collaborator response instructions

This handoff candidate is immutable and read-only. Do not edit, rename, replace,
or delete any file inside it.

1. Create a new response directory outside this candidate root with a unique
   identifier, date, and responsible organization.
2. Copy `platform_feasibility_questions.tsv`,
   `role_and_blinding_firewall.tsv`, `collaborator_deliverable_register.tsv`,
   and only the contract templates needed for the response into that directory.
3. Record the parent candidate ID and SHA256 of
   `STAGE_A_COLLABORATOR_HANDOFF_SEALED.json` in a response manifest before
   editing the copies.
4. Complete every feasibility answer, evidence path, owner, role assignment,
   and acceptance timestamp in the response copies. Do not request or infer a
   target identity.
5. Hash every completed response file and return the response manifest, the
   completed copies, and supporting evidence as a new package.
6. Never write the response back. Keep the sealed candidate unchanged.

A response package is feasibility evidence only. It cannot freeze a target,
authorize outcome generation, open Stage B, or promote the paper.
"""


def source_outputs(root: Path, seal: dict[str, object]) -> list[tuple[Path, str]]:
    outputs = []
    for name, expected in seal["output_sha256"].items():  # type: ignore[union-attr]
        candidates = [root / str(name), root / f"{name}.tsv"]
        paths = [path for path in candidates if path.is_file()]
        if len(paths) != 1 or sha256_file(paths[0]) != expected:
            raise RuntimeError(f"Contract output drift or ambiguity: {root} {name}")
        outputs.append((paths[0], str(expected)))
    return outputs


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite collaborator handoff: {CANDIDATE_ROOT}")
    resolved = {}
    for role, (variable, seal_name, status) in CONTRACTS.items():
        root = candidate_source(variable)
        seal_path = root / seal_name
        seal = json.loads(seal_path.read_text(encoding="utf-8"))
        if seal.get("status") != status:
            raise RuntimeError(f"Invalid {role} contract status")
        if seal.get("scientific_outcomes_inspected", seal.get("experimental_outcomes_inspected")) is not False:
            raise RuntimeError(f"{role} contract reports outcome access")
        resolved[role] = (root, seal_path, seal, source_outputs(root, seal))

    CANDIDATE_ROOT.mkdir(parents=True)
    outputs: list[Path] = []
    handoff_path = CANDIDATE_ROOT / "COLLABORATOR_HANDOFF.md"
    atomic_write_text(handoff_path, HANDOFF)
    outputs.append(handoff_path)
    response_path = CANDIDATE_ROOT / "RESPONSE_INSTRUCTIONS.md"
    atomic_write_text(response_path, RESPONSE_INSTRUCTIONS)
    outputs.append(response_path)
    questions_path = CANDIDATE_ROOT / "platform_feasibility_questions.tsv"
    write_tsv(
        questions_path,
        [{"question_id": row[0], "domain": row[1], "frozen_question": row[2], "impact": row[3], "answer": "", "evidence_path": "", "owner": "", "signed_utc": ""} for row in QUESTIONS],
        ["question_id", "domain", "frozen_question", "impact", "answer", "evidence_path", "owner", "signed_utc"],
    )
    outputs.append(questions_path)
    deliverable_path = CANDIDATE_ROOT / "collaborator_deliverable_register.tsv"
    write_tsv(
        deliverable_path,
        [{"deliverable_id": row[0], "deliverable": row[1], "owner_role": row[2], "required_before": row[3], "immutable_template_paths": row[4], "status": "not_started"} for row in DELIVERABLES],
        ["deliverable_id", "deliverable", "owner_role", "required_before", "immutable_template_paths", "status"],
    )
    outputs.append(deliverable_path)
    roles_path = CANDIDATE_ROOT / "role_and_blinding_firewall.tsv"
    write_tsv(
        roles_path,
        [{"role": row[0], "responsibility": row[1], "target_access": row[2], "outcome_access": row[3], "prohibition": row[4], "assigned_person": "", "accepted_utc": ""} for row in ROLES],
        ["role", "responsibility", "target_access", "outcome_access", "prohibition", "assigned_person", "accepted_utc"],
    )
    outputs.append(roles_path)

    bundle_rows = []
    for role, (root, seal_path, _, source_files) in resolved.items():
        for source, expected in [(seal_path, sha256_file(seal_path)), *source_files]:
            relative = Path("contracts") / role / source.name
            destination = CANDIDATE_ROOT / relative
            atomic_write_text(destination, source.read_text(encoding="utf-8"))
            if sha256_file(destination) != expected:
                raise RuntimeError(f"Copied collaborator contract drift: {destination}")
            outputs.append(destination)
            bundle_rows.append({
                "contract_role": role,
                "source_path": str(source.relative_to(PROJECT_ROOT)),
                "bundle_path": str(relative),
                "size_bytes": source.stat().st_size,
                "sha256": expected,
            })
    bundle_path = CANDIDATE_ROOT / "contract_bundle_manifest.tsv"
    write_tsv(bundle_path, bundle_rows, ["contract_role", "source_path", "bundle_path", "size_bytes", "sha256"])
    outputs.append(bundle_path)

    source_manifest = []
    for role, path in [
        ("handoff_producer", SCRIPT_ROOT / "43_build_stage_a_collaborator_handoff.py"),
        ("handoff_validator", SCRIPT_ROOT / "44_validate_stage_a_collaborator_handoff.py"),
        *[(f"{role}_seal", values[1]) for role, values in resolved.items()],
    ]:
        source_manifest.append({
            "role": role, "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size, "sha256": sha256_file(path),
        })
    source_path = CANDIDATE_ROOT / "handoff_source_manifest.tsv"
    write_tsv(source_path, source_manifest, ["role", "source_path", "size_bytes", "sha256"])
    outputs.append(source_path)
    payload = {
        "status": "sealed_target_independent_stage_a_collaborator_handoff",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "n_input_contracts": len(CONTRACTS),
        "n_feasibility_questions": len(QUESTIONS),
        "n_deliverables": len(DELIVERABLES),
        "n_roles": len(ROLES),
        "primary_source_lineage": "Hepatocytes",
        "primary_recipient_lineages": ["Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Macrophages"],
        "experimental_targets_frozen": False,
        "scientific_conditions_opened": False,
        "scientific_outcomes_inspected": False,
        "wetlab_samples_exist": False,
        "stage_b_design_frozen": False,
        "paper_promotion_authorized": False,
        "bundle_read_only": True,
        "next_gate": "named collaborator feasibility answers and role acceptance before target disclosure",
        "output_sha256": {str(path.relative_to(CANDIDATE_ROOT)): sha256_file(path) for path in outputs},
    }
    atomic_write_json(CANDIDATE_ROOT / "STAGE_A_COLLABORATOR_HANDOFF_SEALED.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
