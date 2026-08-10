#!/usr/bin/env python3
"""Build immutable, target-blind Plan 45A platform-outreach packets."""

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
PLAN = PROJECT_ROOT / "docs/plans/2026-08-07_paper_program/45A_EXPERIMENTAL_PLATFORM_ROUTING_AND_OUTREACH.md"
PARENT_ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-stage-a-collaborator-handoff-v3-2026-08-10"
PARENT_SEAL = PARENT_ROOT / "STAGE_A_COLLABORATOR_HANDOFF_SEALED.json"
OUTREACH_ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-platform-outreach-v3-2026-08-10"

TOP_LEVEL_PARENT_FILES = (
    "COLLABORATOR_HANDOFF.md",
    "RESPONSE_INSTRUCTIONS.md",
    "platform_feasibility_questions.tsv",
    "role_and_blinding_firewall.tsv",
    "collaborator_deliverable_register.tsv",
)

# These legacy audit names must never enter an initial platform packet. They are
# checked without being emitted into the candidate.
FORBIDDEN_TARGET_TOKENS = (
    "ACSL5", "ERCC2", "ZBTB41", "IL18R1", "SHMT1", "CENPQ", "NECAB2",
    "PNPLA6", "RANBP17",
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
    "du_tsinghua": {
        "rank": 1,
        "organization": "Yanan Du laboratory, Tsinghua University",
        "role": "core_relay_candidate",
        "send_priority": "primary_parallel",
        "contact_name": "Professor Yanan Du",
        "contact_value": "duyanan@tsinghua.edu.cn",
        "contact_status": "verified_in_2026_pubmed_record",
        "contact_url": "https://pubmed.ncbi.nlm.nih.gov/41714631/",
        "subject": "Target-blind feasibility inquiry: a lineage-restricted MASLD regulatory relay in five-lineage human liver organoids",
        "email": """Dear Professor Du and colleagues,

We have built a donor-aware MASLD atlas that separates inherited regulatory susceptibility from the established multicellular disease state. We are now planning one prospective, source-independent experiment: test whether an exactly oriented regulatory operation in one liver lineage produces a delayed response in unedited recipient lineages, and whether target rescue or relay blockade breaks the phenotype.

Your five-lineage organoid is the closest published match to the required relay architecture because it includes hepatocytes, cholangiocytes, endothelial cells, stellate cells, and Kupffer cells and has resolved a secreted HSC-to-hepatocyte mechanism. Before disclosing any locus or guide, could your team complete a target-independent feasibility review covering three independent iPSC backgrounds, separately assembled lineage mosaics, chronic metabolic challenge, lineage-resolved RNA and secretome, two non-RNA phenotypes, blinded biological replication, and a path to exact endogenous editing?

We have frozen the experimental gates and will report a valid null. At this stage we request no experiment and disclose no target; we ask only for an evidence-backed platform response and a short feasibility discussion.

Sincerely,

[PI and team]""",
        "questions": [
            ("DU01", "FQ01;FQ02;FQ03", "Can all five lineages be generated from the same iPSC background and repeated in at least three genetically independent iPSC backgrounds with two independent differentiations per condition?"),
            ("DU02", "FQ01;FQ04;FQ05", "Can source and recipient lineages be differentiated separately, tracked, and assembled in source-only, source-restricted, reciprocal, and all-lineage configurations?"),
            ("DU03", "FQ05;FQ08", "Can promoter CRISPRa or CRISPRi be restricted to the source lineage using at least two guides without detectable recipient carryover?"),
            ("DU04", "FQ06;FQ08", "Can the platform tolerate a 10-to-21-day repeated metabolic challenge while preserving all five lineages, maturation, viability, and composition?"),
            ("DU05", "FQ07", "Can direct mosaics and conditioned-medium, transwell, or purified-secretome transfer be generated from the same source differentiation?"),
            ("DU06", "FQ08", "Can lineage RNA, secretome, target engagement, composition, lipid or metabolic phenotype, and matrix or injury phenotype be measured from one frozen design?"),
            ("DU07", "FQ11;FQ12", "Can exact-edited iPSC clones generated by a partner be accepted, blinded, differentiated, rescued, and tested under the same core platform?"),
        ],
        "evidence": [
            ("five_required_lineages", "publicly_supported", "https://www.nature.com/articles/s41467-026-69548-0"),
            ("culture_through_35_days", "publicly_supported", "https://pubmed.ncbi.nlm.nih.gov/41714631/"),
            ("non_cell_autonomous_secreted_relay", "publicly_supported", "https://www.nature.com/articles/s41467-026-69548-0"),
            ("three_independent_iPSC_backgrounds", "not_publicly_demonstrated", "https://pubmed.ncbi.nlm.nih.gov/41714631/"),
            ("chronic_MASLD_context", "not_publicly_demonstrated", "https://pubmed.ncbi.nlm.nih.gov/41714631/"),
            ("exact_endogenous_editing_in_platform", "not_publicly_demonstrated", "https://pubmed.ncbi.nlm.nih.gov/41714631/"),
        ],
    },
    "takebe_cincinnati": {
        "rank": 2,
        "organization": "Takanori Takebe laboratory, Cincinnati Children's",
        "role": "integrated_MASLD_genetics_candidate",
        "send_priority": "primary_parallel",
        "contact_name": "Professor Takanori Takebe",
        "contact_value": "takanori.takebe@cchmc.org",
        "contact_status": "verified_on_official_lab_contact_page",
        "contact_url": "https://www.cincinnatichildrens.org/research/divisions/g/gastroenterology/labs/takebe/contact",
        "subject": "Outcome-blind collaboration inquiry: multi-background MASLD organoids for an exact regulatory risk-to-state relay",
        "email": """Dear Professor Takebe and colleagues,

We have frozen a MASLD risk-to-state experiment before selecting its final locus. The required test combines a genetically oriented source-lineage perturbation, reciprocal isogenic mosaics, delayed recipient-lineage responses, route transfer, target rescue, mediator blockade, and orthogonal lipid and tissue-remodeling phenotypes.

Your multi-donor steatohepatitis organoid panel, genetic editing work, and newer multi-lineage immune FLO system uniquely cover the multi-background genetics, MASLD context, and lineage complexity. Could your team advise—without target disclosure—whether separately differentiated protective and edited lineages can be recombined into reciprocal mosaics across at least three backgrounds and whether the same batches can support direct and soluble-transfer arms through a chronic late readout?

We can provide a sealed target-independent handoff and response table. Only after feasibility, corrected-colocalization orientation, accessibility, guideability, and external protocol audits pass would a target be disclosed.

Sincerely,

[PI and team]""",
        "questions": [
            ("TA01", "FQ01;FQ04", "Which current workflow can preserve separately identifiable hepatocyte, cholangiocyte, endothelial, mesenchymal, and macrophage recipients through an adult MASLD challenge?"),
            ("TA02", "FQ02;FQ03;FQ05", "Can lineages from protective-isogenic and edited iPSC clones be differentiated separately and recombined into reciprocal source-recipient mosaics across three backgrounds?"),
            ("TA03", "FQ11", "Can the published editing workflow be generalized to a target-blind noncoding variant with at least two independent clones per allele and three backgrounds?"),
            ("TA04", "FQ07;FQ09", "Can source-only, direct-mosaic, and conditioned-medium or transwell arms share the same differentiation batches and biological-unit structure?"),
            ("TA05", "FQ06;FQ08", "Can the chronic late readout include lineage RNA, secretome, composition, viability or maturation, lipid or metabolic phenotype, and a distinct matrix or injury phenotype?"),
            ("TA06", "FQ12", "Can target-dosage rescue and recipient mediator blockade be performed without changing the endogenous edit?"),
        ],
        "evidence": [
            ("multi_donor_steatohepatitis_panel", "publicly_supported", "https://pmc.ncbi.nlm.nih.gov/articles/PMC9617783/"),
            ("common_risk_variant_editing", "publicly_supported", "https://pmc.ncbi.nlm.nih.gov/articles/PMC9617783/"),
            ("ten_PSC_line_immune_lineage_platform", "publicly_supported", "https://www.sciencedirect.com/science/article/pii/S0168827825026571"),
            ("hepatobiliary_endothelial_mesenchymal_immune_outputs", "publicly_supported", "https://www.sciencedirect.com/science/article/pii/S0168827825026571"),
            ("reciprocal_separately_assembled_mosaics", "not_publicly_demonstrated", "https://www.sciencedirect.com/science/article/pii/S0168827825026571"),
            ("same_batch_direct_and_soluble_transfer", "not_publicly_demonstrated", "https://pmc.ncbi.nlm.nih.gov/articles/PMC9617783/"),
        ],
    },
    "ebrahimkhani_pittsburgh": {
        "rank": 3,
        "organization": "Mo Ebrahimkhani laboratory, University of Pittsburgh",
        "role": "engineering_or_exact_editing_candidate",
        "send_priority": "primary_parallel",
        "contact_name": "Professor Mo Ebrahimkhani",
        "contact_value": "mo.ebr@pitt.edu",
        "contact_status": "verified_on_official_institutional_page",
        "contact_url": "https://engineering.pitt.edu/people/faculty/mo-reza-ebrahimkhani/",
        "subject": "Target-blind engineering inquiry: exact-edited, multilineage liver organoids for a regulatory relay experiment",
        "email": """Dear Professor Ebrahimkhani and colleagues,

We are seeking an engineering partner for a prospective MASLD mechanism test whose target is not yet disclosed. The experiment requires an exact noncoding allele in at least three independent iPSC backgrounds, source-restricted and reciprocal multicellular mosaics, delayed recipient responses, rescue and blockade, and quantitative lipid plus physical matrix or injury phenotypes.

Your work combines CRISPR-based programming of vascularized multilineage liver organoids with a three-background isogenic fibrosis model and physical collagen readouts. Before target selection, could you assess whether your current platform can add a separately introduced macrophage recipient and run the full relay experiment, or instead serve as the blinded exact-editing and clone-QC partner for a five-lineage core platform?

The full design and null-reporting rules are frozen. We ask first for a signed capability response with evidence paths and a brief technical discussion; no biological outcome or target is requested.

Sincerely,

[PI and team]""",
        "questions": [
            ("EB01", "FQ02;FQ03;FQ11", "Can the exact-edit, three-background workflow be applied to a target-blind noncoding allele with at least two independent clones per allele and background?"),
            ("EB02", "FQ01;FQ05", "Can edited source-lineage progenitors be combined with protective-isogenic recipient progenitors without compromising multilineage differentiation?"),
            ("EB03", "FQ04", "Can macrophages be added as a separately introduced and measurable recipient while retaining hepatocyte, biliary, endothelial, and stellate compartments?"),
            ("EB04", "FQ05;FQ06;FQ07", "Can the platform support chronic metabolic injury, source-only cultures, reciprocal mosaics, and direct-versus-soluble route discrimination?"),
            ("EB05", "FQ08", "Can physical matrix imaging be paired with lipid or metabolic, secretome, lineage RNA, viability, maturation, and composition measurements?"),
            ("EB06", "FQ11;FQ12", "Is the preferred role a complete core platform or a blinded clone-generation and genetic-engineering complement to another core organoid team?"),
        ],
        "evidence": [
            ("hepatocyte_biliary_endothelial_stellate_outputs", "publicly_supported", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8164844/"),
            ("CRISPRa_and_synthetic_gene_circuits", "publicly_supported", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8164844/"),
            ("three_background_isogenic_mutation", "publicly_supported", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8536785/"),
            ("physical_collagen_readouts", "publicly_supported", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8536785/"),
            ("separately_introduced_macrophage_recipient", "not_publicly_demonstrated", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8164844/"),
            ("chronic_MASLD_relay", "not_publicly_demonstrated", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8536785/"),
        ],
    },
    "ge_iorgantech": {
        "rank": 4,
        "organization": "Lungen Lu and Qichao Ge teams with iORGANtech",
        "role": "MASLD_condition_calibration_or_support",
        "send_priority": "supporting_optional",
        "contact_name": "Current corresponding author or iORGANtech scientific contact",
        "contact_value": "",
        "contact_status": "requires_current_confirmation",
        "contact_url": "https://onlinelibrary.wiley.com/doi/abs/10.1111/liv.70716",
        "subject": "Target-independent feasibility inquiry: extending a multilineage MASLD organoid to a regulatory-relay design",
        "email": """Dear Professors Lu and Ge and iORGANtech colleagues,

Your multilineage iPSC liver organoid has demonstrated fatty-acid MASLD-like injury, fibrosis, pharmacologic response, and lineage-restricted genetic perturbation. We are planning a frozen regulatory-relay experiment and would like to determine whether the platform can be extended—before any target is disclosed—to include an endothelial recipient, three independent iPSC backgrounds, reciprocal source-recipient mosaics, direct-versus-soluble transfer, and exact-edited clones.

Could your team complete a target-independent capability table and identify which chronic challenge, lineage, secretome, lipid, injury, and matrix readouts can be performed from the same biological differentiations? A partial answer may still make the platform valuable for outcome-blind MASLD condition calibration, but it would not be treated as a complete mechanistic platform.

Sincerely,

[PI and team]""",
        "questions": [
            ("GE01", "FQ04", "Can an endothelial lineage be incorporated and measured without replacing any of the four reported lineages?"),
            ("GE02", "FQ02;FQ03", "How many independent iPSC backgrounds and independent differentiations have passed the full MASLD and fibrosis protocol?"),
            ("GE03", "FQ05", "Can the lineage-specific perturbation workflow support reciprocal source-lineage and recipient-isogenic mosaics rather than only whole-organoid exposure?"),
            ("GE04", "FQ06;FQ08", "Can the fatty-acid challenge be maintained for 10-to-21 days with stable maturation, composition, and two quantitative non-RNA phenotypes?"),
            ("GE05", "FQ07;FQ09", "Can direct mosaic and conditioned-medium or transwell transfer be generated from the same source differentiation?"),
            ("GE06", "FQ11;FQ12", "Can externally exact-edited iPSC clones enter the platform, or can an editing partner be named before target disclosure?"),
        ],
        "evidence": [
            ("four_lineage_iPSC_organoid", "publicly_supported", "https://pubmed.ncbi.nlm.nih.gov/42220250/"),
            ("fatty_acid_MASLD_like_injury", "publicly_supported", "https://pubmed.ncbi.nlm.nih.gov/42220250/"),
            ("lineage_restricted_genetic_perturbation", "publicly_supported", "https://pubmed.ncbi.nlm.nih.gov/42220250/"),
            ("lipid_injury_and_matrix_phenotypes", "publicly_supported", "https://onlinelibrary.wiley.com/doi/abs/10.1111/liv.70716"),
            ("endothelial_recipient", "not_publicly_demonstrated", "https://pubmed.ncbi.nlm.nih.gov/42220250/"),
            ("three_background_exact_editing_and_route_transfer", "not_publicly_demonstrated", "https://pubmed.ncbi.nlm.nih.gov/42220250/"),
        ],
    },
}


README = """# Target-blind Plan 45A platform outreach

This immutable candidate contains four route-specific feasibility packets. It
does not contain a candidate locus, target identity, guide, allocation,
biological sample, or scientific outcome. No message has been sent.

## Use

1. The PI reviews `route_registry.tsv` and each route's `EMAIL_DRAFT.md`.
2. Send the three `primary_parallel` routes on the same day. The supporting
   route may be sent concurrently or retained for condition calibration.
3. Attach the route directory and the copied parent handoff files only. Do not
   attach a genetics registry or experimental worklist.
4. Record the actual recipient and send metadata in a copy of
   `send_log_template.tsv` outside this immutable root.
5. The collaborator creates a separate editable response and completes both
   the frozen FQ01-FQ12 table and the route-specific addendum.
6. A response permits target disclosure only after independent adjudication.
   It never freezes a target or authorizes Stage A or Stage B.

## One-platform rule

One core biological platform must run all claim-bearing cis, mosaic, transfer,
rescue, blockade, and phenotype arms. A second laboratory may generate and QC
edited clones under a frozen transfer and blinding contract, but results from
different organoid systems cannot be assembled into one mechanism.
"""


RESPONSE_ADDENDUM = """# Route-specific response instructions

This directory is a read-only outreach packet. Do not edit it.

1. Create a new response directory outside this candidate.
2. Copy and complete `platform_feasibility_questions.tsv`,
   `role_and_blinding_firewall.tsv`, `collaborator_deliverable_register.tsv`,
   and `route_specific_questions.tsv`.
3. Cite a protocol, publication, pilot manifest, core-facility statement, or
   named responsible owner for every answer.
4. Preserve backgrounds by independent differentiations as biological units;
   cells, wells, organoids, and lanes remain technical units.
5. Do not request or infer a target, guide, allocation, or outcome.
6. Return the completed response with supporting evidence and hashes within
   14 calendar days of outreach.

Feasibility evidence can permit target disclosure only after independent
adjudication. It cannot freeze a target, authorize outcome generation, open
Stage B, or promote the paper.
"""


def all_payload_files(root: Path) -> list[Path]:
    return sorted(
        path for path in root.rglob("*")
        if path.is_file() and path.name not in {"OUTREACH_PACKET_SEALED.json", "PACKET_SEAL_SHA256.txt", "VALIDATION.txt"}
    )


def main() -> None:
    if OUTREACH_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite platform outreach candidate: {OUTREACH_ROOT}")
    for source in (PLAN, PARENT_SEAL, *[PARENT_ROOT / name for name in TOP_LEVEL_PARENT_FILES]):
        if not source.is_file():
            raise RuntimeError(f"Missing source: {source}")
    parent_seal = json.loads(PARENT_SEAL.read_text(encoding="utf-8"))
    if parent_seal.get("status") != "sealed_target_independent_stage_a_collaborator_handoff":
        raise RuntimeError("Parent handoff is not the authoritative sealed target-independent candidate")
    if parent_seal.get("experimental_targets_frozen") is not False or parent_seal.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("Parent handoff crosses the target/outcome firewall")

    OUTREACH_ROOT.mkdir(parents=True)
    write_text(OUTREACH_ROOT / "README.md", README)
    parent_dir = OUTREACH_ROOT / "parent_handoff"
    parent_dir.mkdir()
    for name in TOP_LEVEL_PARENT_FILES:
        shutil.copyfile(PARENT_ROOT / name, parent_dir / name)

    route_registry = []
    for route_id, route in ROUTES.items():
        route_root = OUTREACH_ROOT / "routes" / route_id
        route_root.mkdir(parents=True)
        write_text(route_root / "EMAIL_DRAFT.md", f"# Email draft\n\n**To:** {route['contact_name']}\n\n**Subject:** {route['subject']}\n\n{route['email']}")
        write_text(route_root / "RESPONSE_INSTRUCTIONS.md", RESPONSE_ADDENDUM)
        for name in ("platform_feasibility_questions.tsv", "role_and_blinding_firewall.tsv", "collaborator_deliverable_register.tsv"):
            shutil.copyfile(PARENT_ROOT / name, route_root / name)
        write_tsv(
            route_root / "route_specific_questions.tsv",
            [
                {
                    "route_question_id": qid,
                    "frozen_fq_links": fq,
                    "target_blind_question": question,
                    "answer": "",
                    "evidence_path": "",
                    "owner": "",
                    "signed_utc": "",
                }
                for qid, fq, question in route["questions"]
            ],
            ["route_question_id", "frozen_fq_links", "target_blind_question", "answer", "evidence_path", "owner", "signed_utc"],
        )
        write_tsv(
            route_root / "public_capability_audit.tsv",
            [
                {"capability": capability, "public_evidence_state": state, "source_url": url, "collaborator_confirmation": ""}
                for capability, state, url in route["evidence"]
            ],
            ["capability", "public_evidence_state", "source_url", "collaborator_confirmation"],
        )
        route_registry.append(
            {
                "route_id": route_id,
                "rank": route["rank"],
                "organization": route["organization"],
                "proposed_role": route["role"],
                "send_priority": route["send_priority"],
                "contact_name": route["contact_name"],
                "contact_value": route["contact_value"],
                "contact_status": route["contact_status"],
                "contact_source_url": route["contact_url"],
                "message_sent": "false",
                "target_disclosed": "false",
                "response_received": "false",
            }
        )
    write_tsv(
        OUTREACH_ROOT / "route_registry.tsv",
        route_registry,
        ["route_id", "rank", "organization", "proposed_role", "send_priority", "contact_name", "contact_value", "contact_status", "contact_source_url", "message_sent", "target_disclosed", "response_received"],
    )
    write_tsv(
        OUTREACH_ROOT / "send_log_template.tsv",
        [
            {
                "route_id": route_id,
                "recipient_used": "",
                "sender": "",
                "sent_utc": "",
                "message_id": "",
                "attachment_manifest_sha256": "",
                "status": "not_sent",
            }
            for route_id in ROUTES
        ],
        ["route_id", "recipient_used", "sender", "sent_utc", "message_id", "attachment_manifest_sha256", "status"],
    )
    write_tsv(
        OUTREACH_ROOT / "parent_handoff_identity.tsv",
        [{
            "parent_candidate": PARENT_ROOT.relative_to(PROJECT_ROOT),
            "parent_seal_path": PARENT_SEAL.relative_to(PROJECT_ROOT),
            "parent_seal_sha256": sha256_file(PARENT_SEAL),
            "parent_status": parent_seal["status"],
            "targets_frozen": "false",
            "outcomes_inspected": "false",
        }],
        ["parent_candidate", "parent_seal_path", "parent_seal_sha256", "parent_status", "targets_frozen", "outcomes_inspected"],
    )
    write_tsv(
        OUTREACH_ROOT / "source_manifest.tsv",
        [
            {"role": "plan45a", "path": PLAN.relative_to(PROJECT_ROOT), "sha256": sha256_file(PLAN), "size_bytes": PLAN.stat().st_size},
            {"role": "parent_handoff_seal", "path": PARENT_SEAL.relative_to(PROJECT_ROOT), "sha256": sha256_file(PARENT_SEAL), "size_bytes": PARENT_SEAL.stat().st_size},
            {"role": "producer", "path": Path(__file__).resolve().relative_to(PROJECT_ROOT), "sha256": sha256_file(Path(__file__).resolve()), "size_bytes": Path(__file__).resolve().stat().st_size},
            {"role": "validator", "path": (SCRIPT_ROOT / "54_validate_platform_outreach_packets.py").relative_to(PROJECT_ROOT), "sha256": sha256_file(SCRIPT_ROOT / "54_validate_platform_outreach_packets.py"), "size_bytes": (SCRIPT_ROOT / "54_validate_platform_outreach_packets.py").stat().st_size},
        ],
        ["role", "path", "sha256", "size_bytes"],
    )

    # Scan before writing the seal so that no forbidden identifier is copied
    # into a route packet or into an email draft.
    scanned = all_payload_files(OUTREACH_ROOT)
    hits = []
    for path in scanned:
        content = path.read_text(encoding="utf-8", errors="ignore")
        for token in FORBIDDEN_TARGET_TOKENS:
            if token in content:
                hits.append((path.relative_to(OUTREACH_ROOT).as_posix(), token))
    if hits:
        raise RuntimeError(f"Target firewall failed with {len(hits)} hit(s)")
    write_tsv(
        OUTREACH_ROOT / "target_firewall_scan.tsv",
        [{
            "files_scanned": len(scanned),
            "forbidden_identifier_hits": 0,
            "target_disclosed": "false",
            "guide_disclosed": "false",
            "allocation_disclosed": "false",
            "outcome_disclosed": "false",
            "verdict": "pass",
        }],
        ["files_scanned", "forbidden_identifier_hits", "target_disclosed", "guide_disclosed", "allocation_disclosed", "outcome_disclosed", "verdict"],
    )

    payload = all_payload_files(OUTREACH_ROOT)
    write_tsv(
        OUTREACH_ROOT / "file_manifest.tsv",
        [
            {"relative_path": path.relative_to(OUTREACH_ROOT).as_posix(), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
            for path in payload
        ],
        ["relative_path", "sha256", "size_bytes"],
    )
    # Rebuild once so the manifest also lists and authenticates the firewall
    # report but deliberately does not list itself.
    payload = [path for path in all_payload_files(OUTREACH_ROOT) if path.name != "file_manifest.tsv"]
    write_tsv(
        OUTREACH_ROOT / "file_manifest.tsv",
        [
            {"relative_path": path.relative_to(OUTREACH_ROOT).as_posix(), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
            for path in payload
        ],
        ["relative_path", "sha256", "size_bytes"],
    )
    sealed_outputs = all_payload_files(OUTREACH_ROOT)
    seal = {
        "candidate_id": OUTREACH_ROOT.name,
        "status": "sealed_target_blind_platform_outreach",
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "n_routes": 4,
        "n_primary_parallel_routes": 3,
        "n_supporting_routes": 1,
        "parent_handoff_seal_sha256": sha256_file(PARENT_SEAL),
        "bundle_read_only": True,
        "message_sent": False,
        "target_disclosed": False,
        "guide_disclosed": False,
        "condition_allocation_disclosed": False,
        "scientific_outcome_disclosed": False,
        "target_freeze_authorized": False,
        "stage_a_authorized": False,
        "stage_b_authorized": False,
        "paper_promotion_authorized": False,
        "output_sha256": {
            path.relative_to(OUTREACH_ROOT).as_posix(): sha256_file(path) for path in sealed_outputs
        },
    }
    seal_path = OUTREACH_ROOT / "OUTREACH_PACKET_SEALED.json"
    write_text(seal_path, json.dumps(seal, indent=2, sort_keys=True))
    write_text(OUTREACH_ROOT / "PACKET_SEAL_SHA256.txt", sha256_file(seal_path))

    for path in OUTREACH_ROOT.rglob("*"):
        if path.is_file():
            os.chmod(path, 0o444)
    for path in sorted((p for p in OUTREACH_ROOT.rglob("*") if p.is_dir()), reverse=True):
        os.chmod(path, 0o555)
    os.chmod(OUTREACH_ROOT, 0o555)
    print(
        "PLATFORM_OUTREACH_BUILD_PASS "
        f"routes=4 primary=3 supporting=1 targets_disclosed=false seal_sha256={sha256_file(seal_path)}"
    )


if __name__ == "__main__":
    main()
