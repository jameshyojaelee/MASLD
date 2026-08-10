#!/usr/bin/env python3
"""Build an outcome-free paired-human clinical data collaboration packet."""

from __future__ import annotations

import datetime as dt
import json
import os
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_json,
    atomic_write_text, sha256_file, write_tsv,
)


QUESTIONS = [
    ("CF01", "paired_tissue", "Does the cohort contain authoritative baseline and follow-up liver-biopsy pairs from the same participants?", "cohort_blocking"),
    ("CF02", "participant_key", "Can a deidentified participant key link biopsy, histology, intervention, covariates, and molecular assays without inference from sample order?", "cohort_blocking"),
    ("CF03", "histology", "Are baseline and follow-up NAS, steatohepatitis status, fibrosis, inflammation, and ballooning available from central or equivalently harmonized reads?", "cohort_blocking"),
    ("CF04", "response_groups", "Does the cohort contain both histologic improvers and non-improvers under a definition that can be frozen before molecular outcome access?", "cohort_blocking"),
    ("CF05", "rna", "Is RNA-capable paired tissue or an existing genome-wide expression product available for both time points?", "cohort_blocking"),
    ("CF06", "covariates", "Are intervention, time interval, weight or BMI change, sex, age, and major treatment covariates available?", "cohort_blocking"),
    ("CF07", "independence", "Can publications, institutions, trial IDs, participant IDs, and expression fingerprints be audited for overlap with GSE106737, GSE83452, GSE48452, and the other proposed cohort?", "replication_blocking"),
    ("CF08", "locked_replication", "Can this cohort be kept outcome-locked until the other cohort's immutable analysis release exists?", "replication_blocking"),
    ("CF09", "tissue_protein", "Is paired liver tissue available for DIA-MS or another quantitative proteomic assay?", "orthogonal_blocking"),
    ("CF10", "secreted_protein", "Are paired plasma or serum samples, liver effluent, or tissue suitable for a liver-linked secretome experiment available?", "orthogonal_blocking"),
    ("CF11", "spatial", "Is paired fixed or frozen liver tissue available for spatial transcriptomics or spatial proteomics with section-level histology registration?", "orthogonal_blocking"),
    ("CF12", "biospecimen_timing", "Are collection timing, fasting state where relevant, processing delay, storage duration, and freeze-thaw history available for molecular biospecimens?", "orthogonal_blocking"),
    ("CF13", "metadata_first", "Can deidentified metadata and assay availability be reviewed before expression, protein, spatial, or responder-linked outcomes are released?", "governance_blocking"),
    ("CF14", "permissions", "Do consent, ethics, data-use, and publication terms permit the proposed paired multi-omic analysis and reporting of valid null results?", "governance_blocking"),
    ("CF15", "role_separation", "Can a data manager hold participant and response keys while the analysis team freezes QC, exclusions, and models?", "governance_blocking"),
    ("CF16", "source_files", "Can all transferred files be versioned, checksum-bound, and accompanied by a complete data dictionary?", "governance_blocking"),
]


DATA_DICTIONARY = [
    ("participant_id", "string", "deidentified stable participant identifier", "required", "biological unit; never inferred from sample order"),
    ("cohort_id", "string", "institution/trial/biobank cohort identifier", "required", "used for independence and overlap audit"),
    ("sample_id", "string", "unique biospecimen or assay identifier", "required", "must map one-to-one to assay files"),
    ("timepoint", "enum", "baseline or followup", "required", "paired timing must be authoritative"),
    ("time_interval_days", "numeric", "days from baseline to follow-up collection", "required", "no date or PHI required"),
    ("intervention", "string", "treatment, lifestyle, surgery, or observational exposure", "required", "freeze intervention classes before outcome access"),
    ("weight_kg", "numeric", "weight at each time point", "preferred", "or supply BMI if weight unavailable"),
    ("bmi", "numeric", "BMI at each time point", "preferred", "used only if predeclared"),
    ("nas", "numeric", "NAS at each biopsy", "required", "include source scoring framework"),
    ("steatohepatitis_status", "enum", "histologic steatohepatitis status", "required", "baseline and follow-up"),
    ("fibrosis_stage", "ordinal", "documented fibrosis stage", "required", "no inferred stage"),
    ("steatosis_grade", "ordinal", "documented steatosis grade", "preferred", "baseline and follow-up"),
    ("inflammation_grade", "ordinal", "documented lobular inflammation", "preferred", "baseline and follow-up"),
    ("ballooning_grade", "ordinal", "documented ballooning", "preferred", "baseline and follow-up"),
    ("responder_definition_source", "string", "protocol/source for improvement classification", "required", "must predate molecular outcome access"),
    ("sex", "enum", "deposited sex variable", "preferred", "no expression inference for clinical primary"),
    ("age_at_baseline", "numeric", "age in years or approved bin", "preferred", "follow governance requirements"),
    ("rna_assay_id", "string", "paired RNA assay/file identifier", "required", "genome-wide expression required"),
    ("tissue_protein_assay_id", "string", "paired tissue-proteomics identifier", "conditional", "needed for orthogonal tissue gate"),
    ("secreted_protein_assay_id", "string", "paired plasma/serum/effluent identifier", "conditional", "record liver-attribution basis"),
    ("spatial_assay_id", "string", "paired spatial assay and section identifier", "conditional", "must map section to participant/timepoint"),
    ("processing_qc", "string", "collection, processing, storage, and assay QC fields", "required", "freeze exclusions before outcomes"),
]


ASSAYS = [
    ("bulk_rna", "paired liver biopsy", "required", "state reversal"),
    ("histology", "paired liver biopsy", "required", "clinical improvement anchor"),
    ("tissue_proteomics", "paired liver tissue", "required_for_full_goal", "orthogonal tissue translation"),
    ("secreted_proteomics", "paired plasma/serum, effluent, or PCLS", "required_for_full_goal", "liver-linked secreted phenotype"),
    ("spatial_transcriptomics_or_proteomics", "paired liver sections", "required_for_full_goal", "physical-niche reversal"),
]


ROLES = [
    ("clinical_custodian", "cohort provenance, permissions, participant pairing, and data dictionary", "cannot change responder definition after molecular outcome access"),
    ("histopathology_lead", "histology framework, scoring provenance, and blinded adjudication", "cannot inspect molecular state scores before histology freeze"),
    ("data_manager", "deidentified linkage keys, transfers, checksums, and response lock", "must differ from analysis lead"),
    ("assay_lead", "RNA/protein/spatial specimen inventory and assay QC", "cannot exclude participants using claim-bearing outcomes"),
    ("analysis_lead", "prespecified participant-level inference", "cannot hold response key before model freeze"),
    ("independent_validator", "pairing, overlap, hashes, models, multiplicity, and release rederivation", "must differ from analysis lead"),
]


BRIEF = """# Paired-human MASLD resolution and multi-omic translation collaboration

## Ultimate goal

We seek two genuinely independent longitudinal liver-biopsy cohorts to test a
single frozen hypothesis: whether the established MASLD transcriptomic state
reverses with histologic improvement and whether the same change appears in
liver-tissue protein, liver-linked secreted protein, and physical tissue-niche
measurements. This human-resolution lane will be integrated only after separate
validation with an exactly oriented, rescuable regulatory perturbation. The
collaboration is designed to report positive, null, and untestable results
without redefining responders, selecting programs, or replacing cohorts after
outcome access.

## First request: metadata and availability only

Please do not send participant-level molecular outcomes initially. Complete
`clinical_cohort_feasibility_questions.tsv`, `assay_availability_matrix.tsv`,
and the role assignments. Return a deidentified sample/pair census and a data
dictionary sufficient to verify participant independence, paired timing,
histology availability, response-group resolution, assay availability, and
permissions. Cohort A and Cohort B protocols will be frozen together before
either locked replication outcome is inspected.

## Minimum scientific design

- Authoritative baseline/follow-up participant pairs.
- Both histologic improvers and non-improvers.
- Genome-wide paired liver expression.
- Independently powered discovery and locked replication cohorts.
- Nested paired tissue protein, liver-linked secreted protein, and spatial
  tissue whenever the complete biological-discovery claim is intended.
- Participants, never biopsies, arrays, sections, spots, or proteins, are the
  inferential units.
- Responder definition, exclusions, covariates, scoring, and multiplicity are
  frozen before molecular outcomes.

This packet requests no protected health information, target identity,
experimental guide, responder-linked expression result, or paper-promotion
decision.
"""


def source_root() -> Path:
    raw = os.environ.get("PLAN46_GOAL_AUDIT_ROOT", "").strip()
    if not raw:
        raise RuntimeError("PLAN46_GOAL_AUDIT_ROOT is required")
    path = Path(raw)
    path = path if path.is_absolute() else PROJECT_ROOT / path
    return path.resolve()


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite clinical handoff: {CANDIDATE_ROOT}")
    audit_root = source_root()
    audit_path = audit_root / "GOAL_COMPLETION_AUDIT.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if (
        audit.get("status") != "goal_not_achieved_new_data_required"
        or audit.get("paper_promotion_authorized") is not False
        or audit.get("new_scientific_outcomes_inspected") is not False
    ):
        raise RuntimeError("Invalid goal-completion audit")
    for relative, expected in audit["output_sha256"].items():
        path = audit_root / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Goal-completion audit drift: {relative}")
    CANDIDATE_ROOT.mkdir(parents=True)
    outputs: list[Path] = []
    brief_path = CANDIDATE_ROOT / "CLINICAL_DATA_COLLABORATOR_BRIEF.md"
    atomic_write_text(brief_path, BRIEF)
    outputs.append(brief_path)
    questions_path = CANDIDATE_ROOT / "clinical_cohort_feasibility_questions.tsv"
    write_tsv(
        questions_path,
        [{
            "question_id": row[0], "domain": row[1], "frozen_question": row[2],
            "impact": row[3], "cohort_role_candidate": "", "answer": "",
            "evidence_path": "", "owner": "", "signed_utc": "",
        } for row in QUESTIONS],
        ["question_id", "domain", "frozen_question", "impact", "cohort_role_candidate", "answer", "evidence_path", "owner", "signed_utc"],
    )
    outputs.append(questions_path)
    dictionary_path = CANDIDATE_ROOT / "minimum_data_dictionary.tsv"
    write_tsv(
        dictionary_path,
        [{"field": row[0], "type": row[1], "definition": row[2], "requirement": row[3], "constraint": row[4]} for row in DATA_DICTIONARY],
        ["field", "type", "definition", "requirement", "constraint"],
    )
    outputs.append(dictionary_path)
    assay_path = CANDIDATE_ROOT / "assay_availability_matrix.tsv"
    write_tsv(
        assay_path,
        [{
            "assay_id": row[0], "biospecimen": row[1], "requirement": row[2],
            "scientific_role": row[3], "available": "", "n_baseline": "",
            "n_followup": "", "n_authoritative_pairs": "", "evidence_path": "",
        } for row in ASSAYS],
        ["assay_id", "biospecimen", "requirement", "scientific_role", "available", "n_baseline", "n_followup", "n_authoritative_pairs", "evidence_path"],
    )
    outputs.append(assay_path)
    roles_path = CANDIDATE_ROOT / "clinical_role_firewall.tsv"
    write_tsv(
        roles_path,
        [{
            "role": row[0], "responsibility": row[1], "prohibition": row[2],
            "assigned_person": "", "accepted_utc": "",
        } for row in ROLES],
        ["role", "responsibility", "prohibition", "assigned_person", "accepted_utc"],
    )
    outputs.append(roles_path)
    source_path = CANDIDATE_ROOT / "clinical_handoff_source_manifest.tsv"
    sources = [
        ("goal_audit", audit_path),
        ("producer", SCRIPT_ROOT / "51_build_clinical_data_handoff.py"),
        ("validator", SCRIPT_ROOT / "52_validate_clinical_data_handoff.py"),
    ]
    write_tsv(
        source_path,
        [{
            "role": role, "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size, "sha256": sha256_file(path),
        } for role, path in sources],
        ["role", "source_path", "size_bytes", "sha256"],
    )
    outputs.append(source_path)
    seal = {
        "status": "sealed_outcome_free_clinical_data_collaborator_handoff",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "n_feasibility_questions": len(QUESTIONS),
        "n_data_dictionary_fields": len(DATA_DICTIONARY),
        "n_assay_roles": len(ASSAYS),
        "n_firewall_roles": len(ROLES),
        "required_independent_human_cohorts": 2,
        "participant_data_accessed": False,
        "clinical_outcomes_accessed": False,
        "molecular_outcomes_accessed": False,
        "target_identity_present": False,
        "permission_granted": False,
        "paper_promotion_authorized": False,
        "next_gate": "metadata-only responses from two independent clinical collaborators",
        "output_sha256": {
            str(path.relative_to(CANDIDATE_ROOT)): sha256_file(path) for path in outputs
        },
    }
    atomic_write_json(CANDIDATE_ROOT / "CLINICAL_DATA_HANDOFF_SEALED.json", seal)
    print(json.dumps(seal, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
