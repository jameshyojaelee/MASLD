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
        "figures/main/fig3_RNAseq/README.md",
        "figures/main/fig4_validation/README.md",
        "figures/main/fig5_convergence/README.md",
        "masld-atlas-v2/README.md",
        "docs/manuscript/release/README.md",
    ]
    texts = {name: read(name, errors) for name in required_documents}

    paper = texts[contract["authority_document"]]
    if normalized(contract["ultimate_goal"]) not in normalized(paper):
        errors.append("ultimate_goal_drift:docs/PAPER.md")
    for marker in (
        "## Ultimate goal",
        "## The central scientific contribution",
        "## Five-figure story",
        "## Claims we can make",
        "## Claims we cannot make",
        "## Two-paper firewall",
        "## Release acceptance",
    ):
        require(paper, marker, "paper", errors)
    require_order(
        paper,
        [
            "1. Resource design and observability",
            "2. Regulatory genetics",
            "3. Established disease state",
            "4. Molecular and physical context",
            "5. MASLD Gene Catalog",
        ],
        "paper_figure_spine",
        errors,
    )

    status = texts[contract["status_document"]]
    for marker in (
        "Validated update selected for the paper",
        "Corrected rerun active",
        "Fragment-native candidate validated",
        "117 program memberships",
        "Separate future paper",
    ):
        require(status, marker, "status", errors)

    results = texts[contract["results_document"]]
    for marker in (
        "| Pooled samples | 846 | 844 |",
        "| Genes tested | 27,638 | 23,370 |",
        "| TREAT DEGs | 1,918 | 1,616 |",
        # Stage contrast totals. Pinned per contrast since 2026-08-11 so that the
        # earlier-arm counts and cohorts-per-contrast are covered too; the old
        # single "313 / 361 / 306 / 149" row hid why 149 != 132 + 42.
        "| F1 vs F0 | 126 | 187 | 313 | 6 | 8 |",
        "| F2 vs F1 | 187 | 174 | 361 | 6 | 8 |",
        "| F3 vs F2 | 174 | 132 | 306 | 6 | 8 |",
        "| F4 vs F3 | **107** | 42 | 149 | **5** | 7 |",
        "| F1 vs F0 | 313 | 6 | 905 | 717 | 168 | 20 | 165 |",
        "| F4 vs F3 | 149 | 5 | 435 | 365 | 55 | 15 | 52 |",
        "| Frozen Hotspot programs | 117 |",
        "| Non-overlap | 92.4% |",
        "corrected 50-study × 22-chromosome COLOC rerun is in progress",
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

    fig5_manifest = read("figures/main/fig5_convergence/CANDIDATE_MAIN_PANELS.tsv", errors)
    require(fig5_manifest, "candidate_not_promoted", "figure5_manifest", errors)
    for retired in ("therapeutic_axes", "calibration", "top10", "top20"):
        if retired in fig5_manifest:
            errors.append(f"retired_figure5_panel_in_candidate_manifest:{retired}")

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
        "figures/main/fig5_convergence/README.md",
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
