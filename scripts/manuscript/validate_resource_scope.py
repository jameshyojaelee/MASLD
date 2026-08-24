#!/usr/bin/env python3
"""Fail-closed checks for the consolidated standalone MASLD Resource scope.

The validator checks document identity and execution boundaries. It does not
claim that scientific results are correct and performs no writes.
"""

from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = ROOT / "config/resource_paper_scope.config.json"
EXEMPTION_PATH = ROOT / "config/method_exemptions.json"

# The five human authorities. A method exemption licenses retired ordering
# vocabulary inside its own workstream; it must never license that vocabulary in
# a Resource-level document, and no Resource document may cite a path produced
# under an exemption. Without this firewall an exemption granted for a benchmark
# would silently become a Resource claim, which is how these terms leaked back
# into the project before.
RESOURCE_AUTHORITIES = [
    "docs/README.md",
    "docs/PAPER.md",
    "docs/STATUS.md",
    "docs/RESULTS.md",
    "docs/ROADMAP.md",
]


def load_scope_policy() -> dict:
    """Return the method-exemption policy, or an empty fail-closed policy.

    A missing config means no exemptions exist, so nothing is firewalled and
    nothing is licensed. Fail closed in both directions.
    """
    if not EXEMPTION_PATH.is_file():
        return {}
    return json.loads(EXEMPTION_PATH.read_text())


def active_exemption_terms(config: dict) -> tuple[list[str], list[str]]:
    """Return (exempt_terms, scope_paths) across active method exemptions."""
    terms: list[str] = []
    scopes: list[str] = []
    for entry in config.get("exemptions", []):
        if entry.get("status") != "active":
            continue
        terms.extend(entry.get("terms", []))
        scopes.extend(entry.get("scope_paths", []))
    return sorted(set(terms)), sorted(set(scopes))


def check_claim_firewall(errors: list[str]) -> None:
    config = load_scope_policy()
    exempt_terms, scope_paths = active_exemption_terms(config)
    if not exempt_terms and not scope_paths:
        return
    never_exemptible = {
        term.lower() for term in config.get("never_exemptible_terms", [])
    }
    approved_by_document: dict[str, set[str]] = {}
    for approval in config.get("resource_language_approvals", []):
        if approval.get("status") != "active":
            continue
        approval_id = approval.get("approval_id", "missing_id")
        decision_record = approval.get("decision_record", "")
        if not decision_record or not (ROOT / decision_record).is_file():
            errors.append(f"missing_language_approval_decision:{approval_id}")
        approved_terms = {
            term.lower() for term in approval.get("terms", [])
        }
        unknown_terms = sorted(
            approved_terms - {term.lower() for term in exempt_terms}
        )
        for term in unknown_terms:
            errors.append(f"language_approval_without_method_exemption:{approval_id}:{term}")
        forbidden = sorted(approved_terms & never_exemptible)
        for term in forbidden:
            errors.append(f"never_exemptible_language_approval:{approval_id}:{term}")
        for relative in approval.get("scope_documents", []):
            if relative not in RESOURCE_AUTHORITIES:
                errors.append(
                    f"language_approval_non_authority_scope:{approval_id}:{relative}"
                )
                continue
            approved_by_document.setdefault(relative, set()).update(approved_terms)
        citation_marker = approval.get("required_citation_marker", "")
        if not citation_marker:
            errors.append(f"missing_language_approval_citation_marker:{approval_id}")
        for relative in approval.get("citation_documents", []):
            path = ROOT / relative
            if not path.is_file():
                errors.append(
                    f"missing_language_approval_citation_document:{approval_id}:{relative}"
                )
            elif citation_marker and citation_marker not in path.read_text(
                errors="replace"
            ):
                errors.append(
                    f"missing_language_approval_citation:{approval_id}:{relative}"
                )
        for requirement in approval.get("required_routing_markers", []):
            relative = requirement.get("document", "")
            marker = requirement.get("marker", "")
            path = ROOT / relative
            if not relative or not marker:
                errors.append(f"invalid_language_approval_routing_marker:{approval_id}")
            elif not path.is_file():
                errors.append(
                    f"missing_language_approval_routing_document:{approval_id}:{relative}"
                )
            elif marker not in path.read_text(errors="replace"):
                errors.append(
                    f"missing_language_approval_routing_marker:{approval_id}:{relative}:{marker}"
                )
        for relative in approval.get("required_artifacts", []):
            path = ROOT / relative
            if not path.is_file() or path.stat().st_size == 0:
                errors.append(
                    f"missing_language_approval_artifact:{approval_id}:{relative}"
                )
            elif path.name == "validation_summary.json":
                summary = json.loads(path.read_text())
                passed = summary.get("status") == "PASS" or summary.get(
                    "validation_pass"
                ) is True
                if not passed:
                    errors.append(
                        f"failed_language_approval_artifact:{approval_id}:{relative}"
                    )
    for relative in RESOURCE_AUTHORITIES:
        path = ROOT / relative
        if not path.is_file():
            continue
        lowered = path.read_text(errors="replace").lower()
        for term in exempt_terms:
            if term.lower() in approved_by_document.get(relative, set()):
                continue
            # Bare mentions inside a fenced code block or a path are unavoidable
            # when the document explains the exemption itself; require the term
            # to be absent from prose by checking word-ish boundaries.
            if re.search(rf"(?<![\w/-]){re.escape(term.lower())}(?![\w/-])", lowered):
                errors.append(f"exempt_term_in_resource_authority:{relative}:{term}")
        for scope in scope_paths:
            if scope.lower() in lowered:
                errors.append(f"exempt_scope_cited_by_resource:{relative}:{scope}")


def normalized(text: str) -> str:
    return re.sub(r"\s+", " ", text.replace("\n", " ")).strip()


def read(relative: str, errors: list[str]) -> str:
    path = ROOT / relative
    if not path.is_file() or path.stat().st_size == 0:
        errors.append(f"missing_or_empty:{relative}")
        return ""
    return path.read_text(errors="replace")


def require(text: str, needle: str, label: str, errors: list[str]) -> None:
    if needle not in text:
        errors.append(f"missing_marker:{label}:{needle}")


def require_order(text: str, tokens: list[str], label: str, errors: list[str]) -> None:
    cursor = -1
    lower = text.lower()
    for token in tokens:
        found = lower.find(token.lower(), cursor + 1)
        if found < 0:
            errors.append(f"ordered_marker_missing:{label}:{token}")
            return
        cursor = found


def main() -> int:
    errors: list[str] = []
    contract = json.loads(CONTRACT_PATH.read_text())
    if len(contract.get("main_figure_roles", [])) != 6:
        errors.append("main_figure_role_count_drift:expected_6")
    if contract.get("central_scientific_contribution") != (
        "inherited_genetic_variation_and_multicellular_disease_remodeling"
    ):
        errors.append("central_scientific_contribution_drift")
    expected_states = {
        "supported",
        "discordant",
        "tested_negative",
        "untestable",
        "not_applicable",
        "source_dependent",
        "indeterminate",
    }
    if set(contract.get("evidence_state_vocabulary", [])) != expected_states:
        errors.append("evidence_state_vocabulary_drift")
    required_decision_fields = {
        "claim",
        "assay_applicability",
        "biological_unit",
        "coverage_and_join_gate",
        "evidence_state",
        "state_reason",
        "unresolved_alternative",
        "next_discriminating_experiment_rule_id",
    }
    if set(contract.get("figure_decision_contract", [])) != required_decision_fields:
        errors.append("figure_decision_contract_drift")

    required_documents = [
        "docs/README.md",
        contract["authority_document"],
        contract["status_document"],
        contract["results_document"],
        contract["roadmap_document"],
        "docs/manuscript/README.md",
        "docs/manuscript/draft/title.md",
        "docs/manuscript/draft/abstract.md",
        "docs/manuscript/draft/01_intro.md",
        "docs/manuscript/draft/fig1.md",
        "docs/manuscript/draft/fig2.md",
        "docs/manuscript/draft/fig3.md",
        "docs/manuscript/draft/fig4.md",
        "docs/manuscript/draft/fig5.md",
        "docs/manuscript/draft/fig5_discussion.md",
        "docs/manuscript/METHODS.md",
        "docs/manuscript/05_figure_legends.md",
        "docs/technical/DATASETS_AND_PIPELINES.md",
        "docs/technical/RELEASE_PROVENANCE.md",
        "docs/technical/NUMBERS_HISTORY.md",
        "docs/literature/REVIEW.md",
        "docs/archive/INDEX.md",
        "figures/main/fig1_atlas_overview/README.md",
        "figures/main/fig2_genetics/README.md",
        "figures/main/INDEX.md",
        "figures/main/fig3_bulk_transcriptomics/README.md",
        "figures/main/fig3_bulk_transcriptomics/CANDIDATE_PANEL_INDEX.tsv",
        "figures/main/fig4_singlecell_programs/README.md",
        "figures/main/fig4_singlecell_programs/CANDIDATE_PANEL_INDEX.tsv",
        "figures/main/fig5_molecular_context/README.md",
        "figures/main/fig5_molecular_context/CANDIDATE_PANEL_INDEX.tsv",
        "figures/main/fig6_gene_catalog/README.md",
        "figures/main/fig6_gene_catalog/CANDIDATE_PANEL_INDEX.tsv",
        "masld-atlas-v2/README.md",
        "docs/manuscript/release/README.md",
    ]
    texts = {name: read(name, errors) for name in required_documents}

    reader_facing_documents = (
        "docs/README.md",
        contract["authority_document"],
        contract["status_document"],
        contract["results_document"],
        contract["roadmap_document"],
        "docs/manuscript/README.md",
        "docs/manuscript/METHODS.md",
        "docs/manuscript/05_figure_legends.md",
        "docs/manuscript/draft/title.md",
        "docs/manuscript/draft/abstract.md",
        "docs/manuscript/draft/01_intro.md",
        "docs/manuscript/draft/fig1.md",
        "docs/manuscript/draft/fig2.md",
        "docs/manuscript/draft/fig3.md",
        "docs/manuscript/draft/fig4.md",
        "docs/manuscript/draft/fig5.md",
        "docs/manuscript/draft/fig5_discussion.md",
    )
    for relative in reader_facing_documents:
        text = texts[relative] if relative in texts else read(relative, errors)
        for term in ("frozen", "progression", "cascade"):
            if re.search(rf"\b{term}\b", text, flags=re.IGNORECASE):
                errors.append(f"retired_reader_term:{relative}:{term}")

    paper = texts[contract["authority_document"]]
    if normalized(contract["ultimate_goal"]) not in normalized(paper):
        errors.append("ultimate_goal_drift:docs/PAPER.md")
    require(
        normalized(paper),
        contract["paper_title"],
        "paper_binding_title",
        errors,
    )
    for marker in (
        "## Ultimate goal",
        "## The central scientific contribution",
        "## Six-figure story",
        "## Claims we can make",
        "## Claims we cannot make",
        "## Two-paper firewall",
        "## Release acceptance",
    ):
        require(paper, marker, "paper", errors)
    require(
        normalized(paper),
        "The paper connects **inherited genetic variation** to **multicellular disease remodeling** in human MASLD.",
        "paper_central_scientific_contribution",
        errors,
    )
    require_order(
        paper,
        [
            "1. A multimodal human MASLD Resource",
            "2. Inherited genetic variation",
            "3. Multicohort disease remodeling",
            "4. Multicellular programs",
            "5. Molecular and physical tissue context",
            "6. MASLD Gene Catalog",
        ],
        "paper_figure_spine",
        errors,
    )

    figure_contract_markers = {
        "docs/manuscript/draft/fig1.md": "next discriminating experiment",
        "docs/manuscript/draft/fig2.md": "next discriminating experiment",
        "docs/manuscript/draft/fig3.md": "next discriminating experiment",
        "docs/manuscript/draft/fig4.md": "next discriminating experiment",
        "docs/manuscript/draft/fig5.md": "next discriminating experiment",
        "docs/manuscript/draft/fig5_discussion.md": "deterministic next-experiment rule",
        "docs/manuscript/05_figure_legends.md": "machine-readable evidence state",
        "figures/main/INDEX.md": "Shared Figure 1–6 decision contract",
    }
    for relative, marker in figure_contract_markers.items():
        require(
            normalized(texts[relative]),
            marker,
            f"figure_evidence_contract:{relative}",
            errors,
        )

    status = texts[contract["status_document"]]
    for marker in (
        "Validated update selected for the paper",
        # Updated 2026-08-17: the corrected COLOC rerun completed and was
        # promoted to canonical, so the marker pins the promoted state.
        "Corrected rerun PROMOTED to canonical 2026-08-17",
        "Five-cohort synchronized candidate",
        "117 program memberships",
        "Separate future paper",
    ):
        require(status, marker, "status", errors)

    results = texts[contract["results_document"]]
    for marker in (
        "| Pooled samples | 846 | 844 |",
        "| Genes tested | 27,638 | 23,370 |",
        # Canonical DEG gate migrated 2026-08-12 from TREAT (interval null at
        # lfc=0.25) to the conventional padj<0.05 & |log2FC|>0.50. Both rows are
        # pinned so neither the canonical nor the retained comparator can drift.
        "| Canonical DEGs | 1,853 |",
        "| TREAT DEGs (sensitivity arm) | 1,918 | 1,616 |",
        "complete 23,370-gene Benjamini–Hochberg families",
        "| Prespecified Hotspot programs | 117 |",
        # Updated 2026-08-17: recomputed on the promoted COLOC release
        # (437 jointly testable genetic genes, 20 overlapping, 14,920 background).
        "| Non-overlap | **95.4%** |",
        "| Primary SuSiE-nominated genes | **462** |",
        "COLOC rerun is COMPLETE and PROMOTED to canonical",
        "| CosMx complete-program observability | 12 / 117 programs |",
        "| ATAC `indeterminate` peaks | 38,066 / 39,914 |",
        "| Catalog strict `tested_negative` calls | 0 |",
    ):
        require(results, marker, "results", errors)

    roadmap = texts[contract["roadmap_document"]]
    for marker in (
        "## Critical path",
        "## Final release checklist",
        "Plans 40–44",
        "Plans 45–46C",
        "No Cas13 result",
        "user explicitly approves promotion",
        "every main figure must expose the claim",
    ):
        require(roadmap, marker, "roadmap", errors)

    removed_legacy_documents = (
        "docs/CODEBASE_CURRENT_STATE.md",
        "docs/RESOURCE_PAPER_SCOPE.md",
        "docs/paper_outline.md",
        "docs/paper_narrative.md",
        "docs/progress.md",
        "docs/competitor_paper_details.md",
        "docs/hmsma_access_request_2026-08-10.md",
        "docs/single_cell_analysis.md",
        "docs/manuscript/NUMBERS.md",
        "docs/manuscript/working",
        "docs/plans",
        "docs/progress",
        "docs/audits",
    )
    for relative in removed_legacy_documents:
        if (ROOT / relative).exists():
            errors.append(f"obsolete_active_document_reappeared:{relative}")

    guard_markers = {
        "scripts/figures/run_all_figures.sbatch": "ALLOW_LEGACY_FIGURE_REGEN",
        "scripts/figures/run_regen_fig2_fig5_panels.sh": "ALLOW_LEGACY_FIGURE_REGEN",
        "scripts/manuscript/run_manuscript_release.sh": "ALLOW_LEGACY_MANUSCRIPT_RELEASE",
        "scripts/portal/rebuild_web_data.sbatch": "ALLOW_LEGACY_PORTAL_REBUILD",
        "scripts/portal/stage_hf_dataset.py": "ALLOW_LEGACY_PORTAL_REBUILD",
        "scripts/portal/generate_convergence_evidence_json.py": "ALLOW_LEGACY_PORTAL_REBUILD",
        "masld-atlas-v2/scripts/preprocess_atlas_data.py": "ALLOW_LEGACY_PORTAL_REBUILD",
        "streamlit_convergence/app.py": "ALLOW_LEGACY_PORTAL_DEMO",
        "masld-atlas-v2/scripts/check-resource-release-gate.mjs": "ALLOW_LEGACY_PORTAL_BUILD",
        "RNA-seq/53_ncrna_landscape.R": "ALLOW_HISTORICAL_NCRNA_PIPELINE",
        "RNA-seq/55_ncrna_conservation.R": "ALLOW_HISTORICAL_NCRNA_PIPELINE",
        "RNA-seq/56_ncrna_epigenomic.R": "ALLOW_HISTORICAL_NCRNA_PIPELINE",
        "RNA-seq/57_ncrna_atlas_integration.R": "ALLOW_HISTORICAL_NCRNA_PIPELINE",
    }
    for relative, marker in guard_markers.items():
        require(read(relative, errors), marker, relative, errors)

    fig6_manifest = read("figures/main/fig6_gene_catalog/CANDIDATE_MAIN_PANELS.tsv", errors)
    require(fig6_manifest, "candidate_not_promoted", "figure6_manifest", errors)
    for retired in ("therapeutic_axes", "calibration", "top10", "top20"):
        if retired in fig6_manifest:
            errors.append(f"retired_figure6_panel_in_candidate_manifest:{retired}")

    require(
        texts["masld-atlas-v2/README.md"],
        "not publication-ready and not authorized for",
        "legacy_portal",
        errors,
    )
    require(
        read("masld-atlas-v2/src/app/layout.tsx", errors),
        "Legacy pre-Resource preview",
        "legacy_portal_runtime_banner",
        errors,
    )
    require(
        texts["docs/manuscript/release/README.md"],
        "not** a complete release",
        "historical_release",
        errors,
    )

    public_term_documents = (
        "docs/README.md",
        contract["authority_document"],
        contract["status_document"],
        contract["results_document"],
        contract["roadmap_document"],
        "docs/manuscript/README.md",
        "docs/manuscript/draft/title.md",
        "docs/manuscript/draft/abstract.md",
        "docs/manuscript/draft/01_intro.md",
        "docs/manuscript/draft/fig1.md",
        "docs/manuscript/draft/fig5_discussion.md",
        "docs/manuscript/05_figure_legends.md",
    )
    for relative in public_term_documents:
        text = read(relative, errors)
        if re.search(r"\bpassports?\b", text, flags=re.IGNORECASE):
            errors.append(f"retired_public_term:{relative}:passport")

    public_identity_documents = (
        "README.md",
        "masld-atlas-v2/README.md",
        "masld-atlas-v2/public/README.md",
        "masld-atlas-v2/src/app/layout.tsx",
        "masld-atlas-v2/src/app/downloads/page.tsx",
        "masld-atlas-v2/src/components/sidebar.tsx",
        "masld-atlas-v2/src/components/hero/hero.tsx",
        "figures/main/fig6_gene_catalog/README.md",
    )
    retired_product_names = (
        "MASLD Atlas",
        "Multi-Evidence Atlas",
        "evidence passport",
        "gene passport",
    )
    for relative in public_identity_documents:
        text = read(relative, errors)
        require(text, "MASLD Gene Catalog", relative, errors)
        for retired in retired_product_names:
            if retired.lower() in text.lower():
                errors.append(f"retired_public_product_name:{relative}:{retired}")

    check_claim_firewall(errors)

    if errors:
        print("RESOURCE_SCOPE_VALIDATION\tFAIL")
        for error in errors:
            print(f"ERROR\t{error}")
        return 1

    print("RESOURCE_SCOPE_VALIDATION\tPASS")
    print(f"CONTRACT\t{contract['contract_id']}")
    print(f"REQUIRED_DOCUMENTS\t{len(required_documents)}")
    print(f"MAIN_FIGURE_ROLES\t{len(contract['main_figure_roles'])}")
    print("PUBLIC_GENE_PRODUCT\tMASLD Gene Catalog")
    print("CAS13_FEEDBACK_TO_RESOURCE\tFALSE")
    print(f"PROMOTION_OWNER\t{contract['canonical_promotion_owner']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
