#!/usr/bin/env python3
"""Build immutable, metadata-only Plan 46A clinical outreach packets."""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
import shutil
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[4]
SCRIPT_ROOT = Path(__file__).resolve().parent
PLAN = PROJECT_ROOT / "docs/plans/2026-08-07_paper_program/46A_CLINICAL_COLLABORATOR_ROUTING_AND_OUTREACH.md"
PARENT_ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-clinical-data-handoff-v1-2026-08-10"
PARENT_SEAL = PARENT_ROOT / "CLINICAL_DATA_HANDOFF_SEALED.json"
ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-clinical-outreach-v1-2026-08-10"

PARENT_FILES = (
    "CLINICAL_DATA_COLLABORATOR_BRIEF.md",
    "clinical_cohort_feasibility_questions.tsv",
    "assay_availability_matrix.tsv",
    "clinical_role_firewall.tsv",
    "minimum_data_dictionary.tsv",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


ROUTES = {
    "nash_crn_as116": {
        "rank": 1,
        "cohort_role": "candidate_cohort_A",
        "independence_group": "NASH_CRN",
        "organization": "NASH CRN / AS116 with DCC and NIDDK repository",
        "contact_name": "NASH CRN Steering Committee liaison and AS116/DCC custodians",
        "contact_route": "https://jhuccs1.us/nash/open/ancillary/ancstudies.htm",
        "contact_status": "liaison_and_form_required",
        "subject": "Outcome-blind feasibility inquiry: paired NASH CRN liver RNA-seq, central histology, and matched biospecimens",
        "email": """Dear Drs. Chalasani, Diehl, and NASH CRN DCC colleagues,

We have built a donor- and cohort-aware MASLD transcriptomic atlas and frozen a 117-program disease-state geometry before approaching any new longitudinal cohort. We would like to test one prespecified question: does that geometry reverse within participants whose centrally scored liver histology improves, relative to non-improvers, and does the same change appear in matched protein and tissue-organization measurements?

The public NASH CRN materials describe AS116 RNA-seq of 904 liver-tissue specimens from 754 adults and children. Before requesting participant-level data or viewing any molecular outcome, could you advise whether AS116 contains an authoritative adult baseline-follow-up subset linked to central histology, and whether visit-matched serum, plasma, or residual tissue exists for that subset? We are initially requesting only a deidentified census and assay/data dictionary: numbers of distinct biopsies, unique participants, authoritative pairs by parent study, improver and non-improver resolution, RNA product and QC, matched material, overlap with published accessions, and the appropriate NASH CRN liaison and ancillary-study route.

Our protocol precommits the score, responder rule, exclusions, biological unit, multiplicity, locked replication, and reporting of positive and null results. We can provide the sealed collaboration brief and complete the NASH CRN Study Proposal with an appropriate Steering Committee liaison.

Would a 30-minute metadata-only feasibility discussion with the AS116/DCC custodians be possible?

Sincerely,

[PI and team]""",
        "questions": [
            ("NC01", "CF01;CF02", "How many of the reported RNA-seq records are distinct biopsy events, unique participants, and authoritative adult baseline-follow-up participant pairs?"),
            ("NC02", "CF07", "What are the authoritative paired counts by PIVENS, FLINT, Database strata, and other parent study, and what overlap exists with published accessions?"),
            ("NC03", "CF03;CF04", "How many adult pairs have complete central histology at both visits and resolve into both a source-defined improvement group and a non-improvement group?"),
            ("NC04", "CF05;CF16", "Are raw counts or another genome-wide expression product, library and extraction QC, batch fields, and frozen sample identifiers available through AS116 or the DCC?"),
            ("NC05", "CF09;CF10;CF11", "Which authoritative pairs have visit-matched serum or plasma, residual RNA-capable tissue, FFPE material, unstained slides, or digitized slides?"),
            ("NC06", "CF08", "Can one NASH CRN trial stratum remain outcome-locked while another is analyzed under a frozen protocol?"),
            ("NC07", "CF13;CF14", "What ancillary approval, consent, data-use, authorship, publication-review, and null-reporting rules apply?"),
            ("NC08", "CF12;CF15", "Can processing histories and a deidentified pairing key be held by a data manager until QC and models are frozen?"),
        ],
        "inventory": [
            ("hepatic_tissue_RNAseq_records", "904", "public_aggregate", "https://jhuccs1.us/nash/eep/Slides/CHALASANI_Accomplishments.pdf"),
            ("people_represented", "754", "public_aggregate", "https://jhuccs1.us/nash/eep/Slides/CHALASANI_Accomplishments.pdf"),
            ("AS116_liver_biopsy_RNAseq", "reported_and_shipped", "public_aggregate", "https://jhuccs1.us/nash/open/Anc_Proposals_22Aug22_EXT.pdf"),
            ("authoritative_adult_pair_count", "", "unverified_route_gate", "https://jhuccs1.us/nash/open/ancillary/ancstudies.htm"),
            ("both_response_directions", "", "unverified_route_gate", "https://jhuccs1.us/nash/open/ancillary/SP7.pdf"),
            ("paired_multiomic_material", "", "unverified_route_gate", "https://repository.niddk.nih.gov/pages/external_specimen_access_program"),
        ],
    },
    "maestro_nash": {
        "rank": 2,
        "cohort_role": "candidate_cohort_B",
        "independence_group": "MADRIGAL_MAESTRO",
        "organization": "Madrigal / MAESTRO-NASH investigators",
        "contact_name": "Madrigal Medical Affairs and MAESTRO biospecimen/data custodian",
        "contact_route": "https://madrigalmedical.com/investigator-initiated-studies/",
        "contact_status": "IIS_route_public_custodian_unverified",
        "subject": "Frozen paired-biopsy transcriptomic and physical-remodeling analysis in MAESTRO-NASH",
        "email": """Dear Madrigal Medical Affairs and MAESTRO-NASH investigators,

We are seeking one institutionally independent paired-biopsy cohort for an outcome-locked replication of a frozen MASLD tissue-state geometry. The MAESTRO-NASH resource is uniquely suited because published sources report 782 participants with baseline and week-52 biopsies and 693 paired SHG/TPEF slide sets, allowing one participant framework to connect transcriptomic reversal with central histology and regional matrix remodeling.

Before requesting any participant-level or molecular outcome, we would like a metadata-only feasibility determination: whether residual paired tissue can support genome-wide RNA profiling; whether paired serum or plasma, digital histology or qFibrosis, imaging, weight, and treatment fields can be joined by an authoritative deidentified key; and whether a prespecified secondary analysis that reports valid nulls is permitted through Madrigal's investigator-initiated-study or collaboration process.

The program definitions, scoring, responder contrasts, exclusions, and replication gate are already sealed. If feasible, we would submit the complete protocol and collaboration packet for scientific and governance review.

Could you direct us to the MAESTRO biospecimen and data custodian for a short metadata-only feasibility call?

Sincerely,

[PI and team]""",
        "questions": [
            ("MA01", "CF01;CF02;CF05", "Among the reported paired biopsies, how many have authoritative joins and residual tissue at both visits suitable for genome-wide RNA profiling?"),
            ("MA02", "CF03;CF06;CF09;CF10;CF11", "How many joined pairs have central histology, paired SHG/TPEF or qFibrosis, serum or plasma, imaging, weight change, treatment assignment, and residual tissue?"),
            ("MA03", "CF04;CF08", "Can one responder definition be frozen before molecular access and can the cohort remain locked until Cohort A is released?"),
            ("MA04", "CF12;CF16", "Can all biospecimen processing, storage, assay, and file-QC metadata be transferred with versioned checksums?"),
            ("MA05", "CF13;CF14", "Do consent, sponsor review, publication terms, and the IIS route permit source-independent analysis, valid nulls, and derived non-identifying tables?"),
            ("MA06", "CF15", "Can an independent data manager hold participant, treatment, and response keys until exclusions and models are frozen?"),
        ],
        "inventory": [
            ("baseline_week52_biopsy_pairs", "782", "public_aggregate", "https://www.journal-of-hepatology.eu/article/S0168-8278%2826%2900147-9/fulltext"),
            ("paired_SHG_TPEF_slide_sets", "693", "public_aggregate", "https://www.journal-of-hepatology.eu/article/S0168-8278%2826%2900147-9/fulltext"),
            ("primary_trial_population", "966", "public_aggregate", "https://www.nejm.org/doi/10.1056/NEJMoa2309000"),
            ("paired_RNA_capable_tissue", "", "unverified_route_gate", "https://madrigalmedical.com/investigator-initiated-studies/"),
            ("paired_blood_and_spatial_material", "", "unverified_route_gate", "https://madrigalmedical.com/investigator-initiated-studies/"),
            ("external_analysis_timing_and_permission", "", "unverified_route_gate", "https://madrigalmedical.com/investigator-initiated-studies/"),
        ],
    },
    "essence": {
        "rank": 3,
        "cohort_role": "candidate_cohort_B",
        "independence_group": "NOVO_ESSENCE",
        "organization": "Novo Nordisk / ESSENCE investigators",
        "contact_name": "Novo Nordisk clinical transparency and ESSENCE data/biospecimen custodian",
        "contact_route": "https://clinicaltrials.gov/study/NCT04822181",
        "contact_status": "public_sharing_signal_custodian_unverified",
        "subject": "Metadata-only feasibility inquiry for frozen paired-biopsy MASLD state analysis in ESSENCE",
        "email": """Dear Novo Nordisk Clinical Transparency and ESSENCE investigators,

We have frozen a MASLD tissue-state score and seek an independent longitudinal cohort in which to test whether within-participant transcriptomic change tracks centrally adjudicated histologic improvement. ESSENCE offers the necessary clinical design: paired baseline and week-72 liver biopsies, rigorous central pathology, active treatment and placebo, and longitudinal metabolic and noninvasive measurements.

At this stage we request no participant-level data and no outcome. We ask only whether the first-800 cohort has residual paired tissue suitable for genome-wide RNA profiling, visit-matched serum or plasma, authoritative assay joins, and digital pathology or other physical-remodeling measurements; and whether these resources are eligible under Novo Nordisk's external data or lab-sample sharing or research-partnership process.

Our score, responder definition, exclusions, analysis, locked-replication rule, and commitment to report valid nulls are prespecified. We can provide a sealed collaboration brief and complete protocol after a positive feasibility determination.

Could you identify the ESSENCE data and biospecimen custodian or the correct external-research intake route?

Sincerely,

[PI and team]""",
        "questions": [
            ("ES01", "CF01;CF02;CF03", "How many first-800 participants have authoritative baseline-week72 pairs with complete central histology and residual tissue at both visits suitable for genome-wide RNA profiling?"),
            ("ES02", "CF06;CF10", "Which authoritative pairs have visit-matched serum or plasma and complete weight, metabolic, noninvasive-test, and intervention data?"),
            ("ES03", "CF09;CF11", "Can digital pathology or section-level physical-remodeling measurements be shared or generated without consuming RNA-capable material?"),
            ("ES04", "CF04;CF08", "Can one responder definition and one outcome-locked replication protocol be approved before any molecular projection is shown?"),
            ("ES05", "CF13;CF14", "Does the external data or lab-sample sharing program cover ESSENCE tissue and biospecimens, and what scientific, data-use, publication, and turnaround rules apply?"),
            ("ES06", "CF12;CF15;CF16", "Can an independent data manager provide a checksum-bound pair and assay census with complete processing metadata before outcomes are released?"),
        ],
        "inventory": [
            ("interim_population", "800", "public_aggregate", "https://www.nejm.org/doi/10.1056/NEJMoa2413258"),
            ("central_pathology_design", "six_pathologists_three_pairs", "public_aggregate", "https://pmc.ncbi.nlm.nih.gov/articles/PMC11599791/"),
            ("clinicaltrials_IPD_sharing", "yes", "public_aggregate", "https://clinicaltrials.gov/study/NCT04822181"),
            ("paired_RNA_capable_tissue", "", "unverified_route_gate", "https://clinicaltrials.gov/study/NCT04822181"),
            ("paired_blood_and_spatial_material", "", "unverified_route_gate", "https://www.novonordisk.com/science-and-technology/research-technologies/processing-of-personal-data.html"),
            ("ESSENCE_specific_sample_sharing_and_timing", "", "unverified_route_gate", "https://www.novonordisk.com/science-and-technology/research-technologies/processing-of-personal-data.html"),
        ],
    },
}


README = """# Metadata-only Plan 46A clinical outreach

This immutable candidate contains three route-specific clinical feasibility
packets. It contains public aggregate inventory statements and blank response
tables only. It contains no participant-level data, responder-linked molecular
outcome, permission claim, target identity, or paper-promotion authorization.
No message has been sent.

Send all three routes in parallel after PI review. NASH CRN is a candidate
Cohort A contribution. MAESTRO-NASH and ESSENCE are alternative independent
Cohort B routes; they are not counted as two cohorts unless both independently
pass and are actually used under separate locked protocols.
"""


RESPONSE = """# Metadata-only response instructions

This route directory is read-only. Create a separate response directory and
complete copies of the clinical feasibility, assay-availability, role, and
route-specific question tables. Initially return only an authoritative
deidentified pair and assay census, data dictionary, evidence paths,
permissions/governance status, and named custodians.

Do not send participant-level molecular values, responder-linked expression,
protein or spatial outcomes, direct identifiers, or a favorable-result
summary. A positive feasibility response permits only a blinded census and
overlap/power audit. It does not grant permission, open outcomes, authorize an
analysis, or promote the paper.

Return an evidence-backed response within 14 calendar days of outreach. A
meeting, expression of interest, or unsigned estimate is not a passed gate.
"""


def payload_files(root: Path) -> list[Path]:
    return sorted(
        path for path in root.rglob("*")
        if path.is_file() and path.name not in {"CLINICAL_OUTREACH_SEALED.json", "CLINICAL_OUTREACH_SEAL_SHA256.txt"}
    )


def main() -> None:
    if ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite clinical outreach candidate: {ROOT}")
    for source in (PLAN, PARENT_SEAL, *[PARENT_ROOT / name for name in PARENT_FILES], SCRIPT_ROOT / "56_validate_clinical_outreach_packets.py"):
        if not source.is_file():
            raise RuntimeError(f"Missing clinical outreach source: {source}")
    parent = json.loads(PARENT_SEAL.read_text(encoding="utf-8"))
    if parent.get("status") != "sealed_outcome_free_clinical_data_collaborator_handoff":
        raise RuntimeError("Invalid parent clinical handoff")
    for field in ("participant_data_accessed", "clinical_outcomes_accessed", "molecular_outcomes_accessed", "permission_granted", "paper_promotion_authorized"):
        if parent.get(field) is not False:
            raise RuntimeError(f"Parent clinical firewall crossed: {field}")

    ROOT.mkdir(parents=True)
    write_text(ROOT / "README.md", README)
    parent_dir = ROOT / "parent_handoff"
    parent_dir.mkdir()
    for name in PARENT_FILES:
        shutil.copyfile(PARENT_ROOT / name, parent_dir / name)
    route_rows = []
    for route_id, route in ROUTES.items():
        route_root = ROOT / "routes" / route_id
        route_root.mkdir(parents=True)
        write_text(route_root / "EMAIL_DRAFT.md", f"# Email draft\n\n**To:** {route['contact_name']}\n\n**Subject:** {route['subject']}\n\n{route['email']}")
        write_text(route_root / "RESPONSE_INSTRUCTIONS.md", RESPONSE)
        for name in ("clinical_cohort_feasibility_questions.tsv", "assay_availability_matrix.tsv", "clinical_role_firewall.tsv", "minimum_data_dictionary.tsv"):
            shutil.copyfile(PARENT_ROOT / name, route_root / name)
        write_tsv(
            route_root / "route_specific_questions.tsv",
            [{"route_question_id": qid, "frozen_cf_links": links, "metadata_only_question": question, "answer": "", "evidence_path": "", "owner": "", "signed_utc": ""} for qid, links, question in route["questions"]],
            ["route_question_id", "frozen_cf_links", "metadata_only_question", "answer", "evidence_path", "owner", "signed_utc"],
        )
        write_tsv(
            route_root / "public_inventory_audit.tsv",
            [{"inventory_item": item, "public_value": value, "evidence_state": state, "source_url": url, "custodian_confirmation": ""} for item, value, state, url in route["inventory"]],
            ["inventory_item", "public_value", "evidence_state", "source_url", "custodian_confirmation"],
        )
        route_rows.append({
            "route_id": route_id,
            "rank": route["rank"],
            "candidate_cohort_role": route["cohort_role"],
            "independence_group": route["independence_group"],
            "organization": route["organization"],
            "contact_name": route["contact_name"],
            "contact_route": route["contact_route"],
            "contact_status": route["contact_status"],
            "message_sent": "false",
            "participant_data_received": "false",
            "outcome_received": "false",
            "permission_granted": "false",
            "response_received": "false",
        })
    write_tsv(
        ROOT / "route_registry.tsv",
        route_rows,
        ["route_id", "rank", "candidate_cohort_role", "independence_group", "organization", "contact_name", "contact_route", "contact_status", "message_sent", "participant_data_received", "outcome_received", "permission_granted", "response_received"],
    )
    write_tsv(
        ROOT / "send_log_template.tsv",
        [{"route_id": route_id, "recipient_used": "", "sender": "", "sent_utc": "", "message_id": "", "attachment_manifest_sha256": "", "status": "not_sent"} for route_id in ROUTES],
        ["route_id", "recipient_used", "sender", "sent_utc", "message_id", "attachment_manifest_sha256", "status"],
    )
    write_tsv(
        ROOT / "parent_handoff_identity.tsv",
        [{"parent_candidate": PARENT_ROOT.relative_to(PROJECT_ROOT), "parent_seal_path": PARENT_SEAL.relative_to(PROJECT_ROOT), "parent_seal_sha256": sha256_file(PARENT_SEAL), "participant_data_accessed": "false", "outcomes_accessed": "false", "permission_granted": "false"}],
        ["parent_candidate", "parent_seal_path", "parent_seal_sha256", "participant_data_accessed", "outcomes_accessed", "permission_granted"],
    )
    write_tsv(
        ROOT / "source_manifest.tsv",
        [
            {"role": "plan46a", "path": PLAN.relative_to(PROJECT_ROOT), "sha256": sha256_file(PLAN), "size_bytes": PLAN.stat().st_size},
            {"role": "parent_clinical_handoff_seal", "path": PARENT_SEAL.relative_to(PROJECT_ROOT), "sha256": sha256_file(PARENT_SEAL), "size_bytes": PARENT_SEAL.stat().st_size},
            {"role": "producer", "path": Path(__file__).resolve().relative_to(PROJECT_ROOT), "sha256": sha256_file(Path(__file__).resolve()), "size_bytes": Path(__file__).resolve().stat().st_size},
            {"role": "validator", "path": (SCRIPT_ROOT / "56_validate_clinical_outreach_packets.py").relative_to(PROJECT_ROOT), "sha256": sha256_file(SCRIPT_ROOT / "56_validate_clinical_outreach_packets.py"), "size_bytes": (SCRIPT_ROOT / "56_validate_clinical_outreach_packets.py").stat().st_size},
        ],
        ["role", "path", "sha256", "size_bytes"],
    )
    write_tsv(
        ROOT / "outcome_firewall.tsv",
        [{"public_aggregate_inventory_only": "true", "participant_rows": 0, "direct_identifiers": 0, "responder_linked_molecular_outcomes": 0, "permissions_assumed": 0, "messages_sent": 0, "verdict": "pass"}],
        ["public_aggregate_inventory_only", "participant_rows", "direct_identifiers", "responder_linked_molecular_outcomes", "permissions_assumed", "messages_sent", "verdict"],
    )
    files = [path for path in payload_files(ROOT) if path.name != "file_manifest.tsv"]
    write_tsv(
        ROOT / "file_manifest.tsv",
        [{"relative_path": path.relative_to(ROOT).as_posix(), "sha256": sha256_file(path), "size_bytes": path.stat().st_size} for path in files],
        ["relative_path", "sha256", "size_bytes"],
    )
    outputs = payload_files(ROOT)
    seal = {
        "candidate_id": ROOT.name,
        "status": "sealed_metadata_only_clinical_outreach",
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "n_routes": 3,
        "required_independent_human_cohorts": 2,
        "n_candidate_cohort_A_routes": 1,
        "n_candidate_cohort_B_routes": 2,
        "n_route_questions": sum(len(route["questions"]) for route in ROUTES.values()),
        "bundle_read_only": True,
        "message_sent": False,
        "participant_data_accessed": False,
        "clinical_outcomes_accessed": False,
        "molecular_outcomes_accessed": False,
        "permission_granted": False,
        "analysis_authorized": False,
        "paper_promotion_authorized": False,
        "parent_handoff_seal_sha256": sha256_file(PARENT_SEAL),
        "output_sha256": {path.relative_to(ROOT).as_posix(): sha256_file(path) for path in outputs},
    }
    seal_path = ROOT / "CLINICAL_OUTREACH_SEALED.json"
    write_text(seal_path, json.dumps(seal, indent=2, sort_keys=True))
    write_text(ROOT / "CLINICAL_OUTREACH_SEAL_SHA256.txt", sha256_file(seal_path))
    for path in ROOT.rglob("*"):
        if path.is_file():
            os.chmod(path, 0o444)
    for path in sorted((path for path in ROOT.rglob("*") if path.is_dir()), reverse=True):
        os.chmod(path, 0o555)
    os.chmod(ROOT, 0o555)
    print(
        "CLINICAL_OUTREACH_BUILD_PASS "
        f"routes=3 cohort_A=1 cohort_B=2 participant_data=false outcomes=false seal_sha256={sha256_file(seal_path)}"
    )


if __name__ == "__main__":
    main()
