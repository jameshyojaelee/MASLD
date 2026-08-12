#!/usr/bin/env python3
"""Build a machine-readable audit of the full paper-upgrade objective."""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_json, sha256_file,
    write_tsv,
)


PLAN43 = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09"
PLAN44 = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/multicellular-assembly-response-2026-08-09"
PLAN45_HANDOFF = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-stage-a-collaborator-handoff-v3-2026-08-10"
PLAN45_RESPONSE = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-collaborator-response-draft-v1-2026-08-10"


REQUIREMENTS = [
    ("UG01", "A reproducible established-state transcriptomic geometry is frozen independently of validation outcomes.", "partial", "Plan 20/43 freeze 117 stage-associated programs, but this proves cross-sectional established state rather than chronicity.", "External longitudinal behavior must validate the unchanged geometry.", "The same frozen geometry is testable in both new longitudinal cohorts without recalibration.", "Plan20/Plan43"),
    ("UG02", "The established state is demonstrably chronic rather than a generic acute lipid or culture response.", "failed", "Plan 43 global PCLS culture effect was near zero and Plan 44 controlled cue topology failed, including an opposite TGF-beta prediction.", "A source-independent acute-versus-chronic time-course interaction in a multicellular system.", "Frozen state/lineage effects are stronger at prespecified chronic than acute time points with composition and viability controls.", "Plan45 Stage A"),
    ("UG03", "Transcriptomic reversal reproducibly accompanies histologic improvement across independent paired human cohorts.", "failed", "GSE106737 discovery did not replicate; GSE83452 was overlapping/null, GSE48452 was direction-discordant, and the effective-cohort aligned meta p was 0.912.", "Two genuinely independent, outcome-blind paired-biopsy cohorts with improvers and non-improvers.", "Both locked cohorts show expected paired change, each passes its prespecified test, and equal-cohort meta-analysis passes with no participant overlap.", "New longitudinal cohorts A/B"),
    ("UG04", "The state translates into an orthogonal liver-tissue protein phenotype beyond ordinary replication of large RNA effects.", "partial", "PXD051911 showed tissue concordance, but matching balance failed; this cannot establish evidence-role architecture.", "Paired tissue proteomics in a new longitudinal cohort with abundance/coverage controls.", "A frozen protein projection reverses with histology and survives assay observability plus continuous RNA-effect adjustment.", "New nested tissue proteomics"),
    ("UG05", "The state translates into a liver-linked secreted-protein phenotype.", "failed", "Outcome-unseen PXD052787 was null after continuous bulk-effect control and failed variability balance.", "Longitudinal plasma/serum secretome or liver-effluent/PCLS secretome paired to the new human tissue cohort.", "A frozen liver-secreted protein projection reverses in both cohorts or one locked cohort plus a source-independent liver secretion model.", "New nested secretome"),
    ("UG06", "The state has a donor-resolved physical niche and its organization changes with histologic improvement.", "partial", "M8 spatial organization is robust in one independent cohort and source-dependent in another, but the evidence-class comparison lacked prespecified gene-level covariates and no paired reversal was tested.", "Paired spatial transcriptomic/proteomic tissue from a prespecified subset with authoritative participant and histology joins.", "Donor-level niche coupling and paired niche reversal survive zonation, composition, coverage, and spatially blocked nulls.", "New nested spatial arm"),
    ("UG07", "Genetically anchored regulatory operations shift the independently frozen disease-state geometry in the predicted direction.", "missing", "GSE281160 failed the source gate at five oriented loci and its state-axis outcome remained unopened; Plan 45 has no target or outcome.", "Corrected pair-specific allele orientation followed by source-independent direction-matched Stage A and exact endogenous Stage B perturbation.", "At least one locus passes orientation, cis engagement, delayed recipient response, and exact-allele replication; architecture generalization requires a second locus.", "Plan45 RSR-00--06"),
    ("UG08", "The regulatory effect is time ordered, source-lineage restricted, rescuable, blockable, and reproduced in non-RNA phenotypes.", "missing", "Only outcome-blind Plan 45 contracts exist; no wet-lab sample has been generated.", "Mosaic/transfer experiment with target rescue, mediator blockade, reciprocal lineage control, and two orthogonal phenotypes.", "Every frozen cis, timing, origin, route, rescue, blockade, composition, and phenotype gate passes without target substitution.", "Plan45 Stage A/B"),
    ("UG09", "The genetic-transcriptomic non-overlap is explained as an assay-conditioned risk-to-state architecture.", "not_achieved", "The non-overlap is descriptive; Plans 42--44 failed functional partition, chronic reversal, and topology gates.", "Joint satisfaction of the human reversal/translation package and the exact regulatory relay.", "The same frozen state connects independent human resolution to an oriented, rescuable regulatory perturbation without a universal synthetic score.", "Integrated Plan46/Plan45"),
    ("UG10", "The Resource provides explicit follow-up rules for genetic and disease-state nominations.", "achieved_design_only", "The MASLD Gene Catalog and Plan 45 contracts encode observability, lineage, assay, and falsifying next experiments.", "Prospective use of those rules must predict which validation assay succeeds or fails.", "At least one routed genetic nomination and one routed disease-state program are prospectively validated under their assigned assays.", "Plan50/Plan45"),
    ("UG11", "The paper supports a full biological-discovery and metabolism-journal escalation claim.", "not_achieved", "Plans 41--44 are complete_no_promotion and Plan 45 remains pre-outcome.", "All claim-bearing requirements UG02--UG09 must pass under immutable releases.", "Mechanical full-discovery gate passes; Plan 60 independently promotes figures and prose.", "Plan60"),
]


NEW_DATA = [
    ("NDP01", "Locked longitudinal human cohort A", "participant", "Paired baseline/follow-up liver biopsy, centrally read histology, treatment/weight covariates, bulk RNA-seq.", "Participant IDs, responder definition, frozen state geometry, primary paired contrast, exclusions, covariates, multiplicity.", "Expected paired state reversal in improvers versus non-improvers with all scoring sensitivities.", "Valid null closes the chronic-reversal claim; no alternate responder definition."),
    ("NDP02", "Independent locked longitudinal human cohort B", "participant", "Same minimum fields as cohort A from a non-overlapping institution, trial, or biobank.", "Complete protocol and analysis copied from cohort A before cohort-B outcome access.", "Directionally concordant, independently significant replication and equal-cohort meta-analysis.", "Any overlap or protocol retuning collapses this to discovery evidence, not replication."),
    ("NDP03", "Nested paired liver-tissue proteomics", "participant", "Baseline/follow-up DIA-MS or equivalent from cohort A and preferably cohort B.", "Observable protein universe, normalization, frozen projections, abundance/coverage matching, RNA-effect sensitivity.", "Protein-state reversal tracks histology after observability and continuous RNA-effect control.", "Failure after RNA adjustment is ordinary cross-assay replication only."),
    ("NDP04", "Nested longitudinal secreted-protein arm", "participant", "Paired plasma/serum proteomics with liver-secretion annotation or paired liver-effluent/PCLS secretome.", "Source attribution, assay universe, sample timing, hemolysis/QC, frozen secreted-state projection.", "Secreted projection reverses with histology and agrees with tissue direction.", "Systemic-only signal without liver attribution is supportive, not liver-secretome evidence."),
    ("NDP05", "Paired physical-niche arm", "participant", "Paired spatial transcriptomics or spatial proteomics on a prespecified participant subset.", "Section selection, tissue pairing, histology registration, lineage/zonation/composition covariates, spatial null.", "Frozen niche axes disassemble with improvement at donor level and replicate across cohorts or modalities.", "Spot/section pseudoreplication or outcome-selected regions invalidate the niche claim."),
    ("NDP06", "Direction-matched regulatory Stage A", "background_by_independent_differentiation", "At least three genetic backgrounds, two differentiations, chronic challenge, source/recipient mosaics, transfer arms, secretome and two phenotypes.", "Corrected locus orientation, target/guide freeze, pilot, units, assays, timing, randomization, blinding, QC thresholds.", "Both guides pass cis then delayed recipient relay with origin, route, composition, and phenotype support.", "No locus substitution; a failed gate becomes the terminal biological boundary."),
    ("NDP07", "Exact endogenous regulatory Stage B", "background_by_clone_by_independent_differentiation", "Exact risk/protective alleles, independent clones, parental/sham controls, target rescue, mediator blockade.", "Variant, alleles, lineage, challenge, timing, pegRNAs, mediator, blocker, phenotypes, and power units frozen before outcomes.", "Exact allele reproduces the relay; target rescue and mediator blockade interrupt the predicted steps.", "Direction-matched dosage alone cannot support an endogenous locus mechanism."),
    ("NDP08", "Cross-package integration and release", "cohort_or_locus", "Immutable human, protein, secretome, spatial, Stage-A, and Stage-B releases.", "No universal score; assay-native effects and complete testability matrix only.", "Human resolution and regulatory perturbation connect to the same frozen state under the mechanical full-discovery gate.", "Any failed required arm routes to a partial claim or Resource posture."),
]


CLAIMS = [
    ("CL01", "chronic_reversible_remodeling_state", "UG02;UG03", "Both chronic-context and two-cohort human reversal pass.", "Chronic, histologically reversible tissue-remodeling state."),
    ("CL02", "orthogonal_human_translation", "UG04;UG05;UG06", "Tissue protein, secreted protein, and paired physical-niche gates pass.", "Orthogonally translated and physically organized human state."),
    ("CL03", "regulatory_risk_to_state_relay", "UG07;UG08", "Exact oriented allele, cis target, delayed relay, rescue, blockade, and phenotypes pass.", "Validated locus-specific regulatory relay; two loci required for architecture generalization."),
    ("CL04", "assay_conditioned_risk_to_state_architecture", "UG02;UG03;UG04;UG05;UG06;UG07;UG08;UG09", "Every human and perturbational component passes without outcome-selected substitution.", "Full risk-to-state architecture and high-level biological discovery."),
    ("CL05", "experiment_routing_framework", "UG10", "Prospective routed validations succeed under their assigned assay classes.", "Experiment-routing framework with prospective support."),
    ("CL06", "metabolism_journal_escalation", "CL01;CL02;CL03;CL04", "Plan 60 independently verifies and promotes the complete immutable evidence package.", "Presubmission inquiry justified; acceptance remains editorial."),
]


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite goal-completion audit: {CANDIDATE_ROOT}")
    sources = [
        ("plan43_release", PLAN43 / "release_manifest.tsv"),
        ("plan43_promotion", PLAN43 / "promotion_verdict.tsv"),
        ("plan43_human", PLAN43 / "human_reversal_effects.tsv"),
        ("plan43_protein", PLAN43 / "protein_transportability.tsv"),
        ("plan43_spatial", PLAN43 / "spatial_lineage_support.tsv"),
        ("plan43_cropseq", PLAN43 / "cropseq_risk_bridge.tsv"),
        ("plan44_terminal_release", PLAN44 / "release_manifest_terminal.tsv"),
        ("plan44_promotion", PLAN44 / "promotion_verdict.tsv"),
        ("plan45_handoff_seal", PLAN45_HANDOFF / "STAGE_A_COLLABORATOR_HANDOFF_SEALED.json"),
        ("plan45_response_marker", PLAN45_RESPONSE / "DRAFT_EDITABLE.json"),
        ("producer", SCRIPT_ROOT / "49_build_goal_completion_audit.py"),
        ("validator", SCRIPT_ROOT / "50_validate_goal_completion_audit.py"),
    ]
    for _, path in sources:
        if not path.is_file():
            raise RuntimeError(f"Missing goal-audit source: {path}")
    CANDIDATE_ROOT.mkdir(parents=True)
    matrix_path = CANDIDATE_ROOT / "ultimate_goal_requirement_matrix.tsv"
    write_tsv(
        matrix_path,
        [{
            "requirement_id": row[0], "ultimate_goal_requirement": row[1],
            "current_verdict": row[2], "authoritative_evidence": row[3],
            "decisive_missing_evidence": row[4], "completion_test": row[5],
            "owner_workstream": row[6],
        } for row in REQUIREMENTS],
        ["requirement_id", "ultimate_goal_requirement", "current_verdict", "authoritative_evidence", "decisive_missing_evidence", "completion_test", "owner_workstream"],
    )
    new_data_path = CANDIDATE_ROOT / "minimal_new_data_package.tsv"
    write_tsv(
        new_data_path,
        [{
            "package_id": row[0], "component": row[1], "biological_unit": row[2],
            "minimum_inputs": row[3], "freeze_before_outcome": row[4],
            "pass_gate": row[5], "failure_interpretation": row[6],
        } for row in NEW_DATA],
        ["package_id", "component", "biological_unit", "minimum_inputs", "freeze_before_outcome", "pass_gate", "failure_interpretation"],
    )
    claims_path = CANDIDATE_ROOT / "claim_promotion_logic.tsv"
    write_tsv(
        claims_path,
        [{
            "claim_id": row[0], "claim": row[1], "required_components": row[2],
            "mechanical_gate": row[3], "permitted_language": row[4],
        } for row in CLAIMS],
        ["claim_id", "claim", "required_components", "mechanical_gate", "permitted_language"],
    )
    source_path = CANDIDATE_ROOT / "goal_audit_source_manifest.tsv"
    write_tsv(
        source_path,
        [{
            "role": role, "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size, "sha256": sha256_file(path),
        } for role, path in sources],
        ["role", "source_path", "size_bytes", "sha256"],
    )
    outputs = [matrix_path, new_data_path, claims_path, source_path]
    seal = {
        "status": "goal_not_achieved_new_data_required",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "n_requirements": len(REQUIREMENTS),
        "n_achieved_biological_requirements": 0,
        "n_design_only_requirements": 1,
        "n_new_data_components": len(NEW_DATA),
        "n_claim_rules": len(CLAIMS),
        "public_data_escalation_closed": True,
        "paired_human_reversal_achieved": False,
        "orthogonal_translation_achieved": False,
        "regulatory_perturbation_achieved": False,
        "full_risk_to_state_architecture_achieved": False,
        "paper_promotion_authorized": False,
        "new_scientific_outcomes_inspected": False,
        "retrospective_terminal_evidence_audited": True,
        "next_decisive_inputs": [
            "two independent locked paired-human histology/RNA cohorts",
            "nested tissue-protein, secreted-protein, and physical-niche data",
            "Plan 45 oriented and rescuable exact regulatory perturbation",
        ],
        "output_sha256": {
            str(path.relative_to(CANDIDATE_ROOT)): sha256_file(path) for path in outputs
        },
    }
    atomic_write_json(CANDIDATE_ROOT / "GOAL_COMPLETION_AUDIT.json", seal)
    print(json.dumps(seal, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
