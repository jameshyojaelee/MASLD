#!/usr/bin/env python3
"""Build the fixed five-section candidate manuscript and seal REL-02--04 products."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from rel01_snapshot_candidate import (
    materialize_immutable,
    materialize_immutable_json,
    materialize_immutable_tsv,
)
from rel03_render_candidate_panels import validate_rel02_state
from rel05_validate_candidate import validate_snapshot_base
from release_common import (
    CANDIDATE_ID,
    ReleaseContractError,
    is_relative_to,
    parse_nonnegative_int,
    read_tsv_exact,
    resolve_project_path,
    sha256_file,
)
from release_products import (
    ANALYSIS_MANIFEST_FIELDS,
    CLAIM_FIELDS,
    EXPECTED_RESULTS_SECTIONS,
    EXPECTED_SUPPLEMENTARY_SECTIONS,
    EXPECTED_SUPPLEMENTARY_PANELS,
    FIGURE_SOURCE_MANIFEST_FIELDS,
    JOURNAL_BRANCH,
    MANUSCRIPT_MANIFEST_FIELDS,
    MYOJIN_ROLE,
    NUMBERS_FIELDS,
    TRANSITION_PRODUCT_FIELDS,
    candidate_figure_root,
    ensure_candidate_output_root,
    load_frozen_input_index,
    project_relative_str,
    read_json_object,
)


SECTION_TITLES = {
    "resource_interface": (
        "The integrated human resource defines complementary, coverage-dependent evidence roles"
    ),
    "established_state_transcriptomics": (
        "Donor-resolved transcriptomics identifies reproducible established-state programs"
    ),
    "genetics_context": (
        "Phenotype provenance and biological context delimit regulatory genetics"
    ),
    "physical_context": (
        "Physical context, assay observability, and prespecified external challenges"
    ),
    "evidence_passports": (
        "The MASLD Gene Catalog connects heterogeneous results to discriminating experiments"
    ),
}
RETIRED_PHRASES = (
    "largest atlas",
    "broadest atlas",
    "most modalities",
    "causal gene",
    "reactive gene",
    "universal score",
    "validated biomarker",
    "observed progression",
    "longitudinal progression",
)
ABSTRACT_DATASET_NAMES = ("Yakubovsky", "GSE287826", "Myojin")


def validate_rel03_state(
    candidate: Path, spec_hash: str, fixture_mode: bool
) -> dict[str, object]:
    state = read_json_object(candidate / "manifests/rel03_state.json", "REL-03 state")
    expected = {
        "contract_version": 1,
        "candidate_id": CANDIDATE_ID,
        "state": "REL03_PANELS_FROZEN",
        "input_snapshot_spec_sha256": spec_hash,
        "n_main_figures": 5,
        "n_supplementary_panels": sum(
            len(panels) for panels in EXPECTED_SUPPLEMENTARY_PANELS.values()
        ),
        "normal_font_pt": 6,
        "format": "PDF",
        "myojin_role": MYOJIN_ROLE,
        "fixture_mode": fixture_mode,
        "canonical_promotion_status": "not_promoted",
    }
    for key, value in expected.items():
        if state.get(key) != value:
            raise ReleaseContractError(f"REL-03 state drift for {key}")
    if state.get("figure_source_manifest_sha256") != sha256_file(
        candidate / "manifests/figure_source_manifest.tsv"
    ):
        raise ReleaseContractError("REL-03 figure manifest hash drift")
    return state


def split_ids(value: str) -> list[str]:
    return [item for item in value.split(";") if item]


def build_manuscript_text(
    numbers: list[dict[str, str]], claims: list[dict[str, str]]
) -> tuple[str, set[str], set[str]]:
    number_by_id = {row["number_id"]: row for row in numbers}
    if len(number_by_id) != len(numbers):
        raise ReleaseContractError("numbers ledger contains duplicate number_id")
    claim_by_id = {row["claim_id"]: row for row in claims}
    if len(claim_by_id) != len(claims):
        raise ReleaseContractError("claim ledger contains duplicate claim_id")
    by_section: dict[str, list[dict[str, str]]] = defaultdict(list)
    for claim in claims:
        by_section[claim["section_id"]].append(claim)
    missing = set(EXPECTED_RESULTS_SECTIONS) - set(by_section)
    if missing:
        raise ReleaseContractError(f"manuscript lacks claims for sections: {missing}")
    used_numbers = set()
    used_claims = set()
    abstract = (
        "We organize inherited regulatory evidence and established disease-state "
        "programs as complementary, coverage-dependent evidence roles. Donor-resolved "
        "transcriptomic programs, phenotype-aware genetics, assay-native physical "
        "context, and the provenance-preserving MASLD Gene Catalog defines what is supported, "
        "indeterminate, or untestable. The resulting resource prioritizes falsifiable "
        "next experiments without collapsing modalities into a universal rank."
    )
    if any(name.lower() in abstract.lower() for name in ABSTRACT_DATASET_NAMES):
        raise ReleaseContractError("abstract names a processed external dataset")
    lines = [
        "# A context-aware multimodal resource for metabolic liver disease",
        "",
        f"**Candidate branch:** {JOURNAL_BRANCH}",
        "",
        "## Abstract",
        "",
        abstract,
        "",
        "## Introduction",
        "",
        (
            "Existing atlases motivate a different question: which conclusions remain "
            "supportable after phenotype provenance, assay observability, biological unit, "
            "and source reuse are made explicit? We distinguish inherited regulatory "
            "evidence from established-state programs and test their relationships without "
            "assuming that one evidence role substitutes for another."
        ),
        "",
        "## Results",
        "",
    ]
    for section_id in EXPECTED_RESULTS_SECTIONS:
        lines.extend([f"### {SECTION_TITLES[section_id]}", ""])
        section_claims = sorted(by_section[section_id], key=lambda row: row["claim_id"])
        included_any = False
        for claim in section_claims:
            number_ids = [
                number_id
                for number_id in split_ids(claim["number_ids"])
                if number_id in number_by_id
                and number_by_id[number_id]["manuscript_included"] == "true"
            ]
            if not number_ids:
                continue
            for phrase in split_ids(claim["prohibited_wording"].replace(",", ";")):
                if phrase.strip().lower() in claim["claim_text"].lower():
                    raise ReleaseContractError(
                        f"claim text contains its prohibited wording: {claim['claim_id']}"
                    )
            rendered_numbers = ", ".join(
                f"{number_by_id[number_id]['display_value']} "
                f"{number_by_id[number_id]['number_role'].replace('_', ' ')}"
                for number_id in number_ids
            )
            markers = " ".join(f"NUMBER:{number_id}" for number_id in number_ids)
            lines.extend(
                [
                    (
                        f"{claim['claim_text']} The frozen ledger records {rendered_numbers}. "
                        f"<!-- CLAIM:{claim['claim_id']} {markers} -->"
                    ),
                    "",
                ]
            )
            used_numbers.update(number_ids)
            used_claims.add(claim["claim_id"])
            included_any = True
        if not included_any:
            raise ReleaseContractError(
                f"Results section has no ledger-backed text: {section_id}"
            )
        lines.extend(
            [
                (
                    "The supported inference is bounded by the tested universe, biological "
                    "unit, and source-dependence class recorded in the claim ledger."
                ),
                "",
            ]
        )
    lines.extend(
        [
            "## Supplementary Results",
            "",
            (
                "The complete prespecified Myojin functional stress test remains "
                "supplementary under its frozen verdict. Statistically valid but "
                "nonconfirmatory effects are described as indeterminate unless an explicit "
                "adequate-negative or equivalence criterion was met."
            ),
            "",
            (
                "The k4 and k6 NMF products are retained only as continuous loading "
                "axes. Factor-level cross-seed stability does not make the fuzzy hard "
                "sample partition a reproducible patient subtype system."
            ),
            "",
        ]
    )
    for section_id in EXPECTED_SUPPLEMENTARY_SECTIONS:
        section_claims = sorted(
            by_section.get(section_id, []), key=lambda row: row["claim_id"]
        )
        if not section_claims:
            raise ReleaseContractError(
                f"supplementary section has no ledger-backed claims: {section_id}"
            )
        for claim in section_claims:
            number_ids = [
                number_id
                for number_id in split_ids(claim["number_ids"])
                if number_id in number_by_id
            ]
            rendered = ", ".join(
                f"{number_by_id[number_id]['display_value']} "
                f"{number_by_id[number_id]['number_role'].replace('_', ' ')}"
                for number_id in number_ids
            )
            markers = " ".join(f"NUMBER:{number_id}" for number_id in number_ids)
            lines.extend(
                [
                    f"{claim['claim_text']} The frozen ledger records {rendered}. "
                    f"<!-- CLAIM:{claim['claim_id']} {markers} -->",
                    "",
                ]
            )
            used_numbers.update(number_ids)
            used_claims.add(claim["claim_id"])
    lines.extend(
        [
            "## Discussion",
            "",
            (
                "The organizing principle is that evidence roles are complementary and "
                "their interpretation depends on observability and context. The strongest "
                "rival explanation is that apparent disease programs reflect ordinary "
                "tissue organization or source reuse. Assay-native external tests narrow "
                "that explanation for some contexts while leaving nonconfirmatory results "
                "indeterminate. Independent perturbation and longitudinal sampling are the "
                "next discriminating experiments."
            ),
            "",
            "## Translational boundaries",
            "",
            (
                "MASLD Gene Catalog entries report support, discordance, tested negatives only "
                "under an explicit adequacy criterion, indeterminate results, untestable "
                "assays, and next experiments. They are not probabilities, leaderboards, "
                "treatment recommendations, or prospectively validated clinical tools."
            ),
            "",
            "## Figure order",
            "",
            "1. Dataset interface and evidence observability.",
            "2. Established-state transcriptomics.",
            "3. Genetics and regulatory context.",
            "4. Physical context, assay observability, and prespecified external challenges.",
            "5. MASLD Gene Catalog and translational boundaries.",
            "",
        ]
    )
    manuscript = "\n".join(lines)
    lowered = manuscript.lower()
    for phrase in RETIRED_PHRASES:
        if phrase in lowered:
            raise ReleaseContractError(
                f"candidate manuscript contains retired phrase: {phrase}"
            )
    if "Figure 6" in manuscript or "Figure6" in manuscript:
        raise ReleaseContractError("five-figure manuscript contains Figure 6")
    return manuscript, used_numbers, used_claims


def add_product(
    rows: list[dict[str, object]],
    project: Path,
    spec_hash: str,
    product_id: str,
    phase: str,
    role: str,
    path: Path,
    producer: str,
) -> None:
    if not path.is_file() or path.is_symlink():
        raise ReleaseContractError(
            f"transition product is missing or symlinked: {path}"
        )
    rows.append(
        {
            "product_id": product_id,
            "phase": phase,
            "artifact_role": role,
            "project_relative_path": project_relative_str(project, path),
            "sha256": sha256_file(path),
            "bytes": path.stat().st_size,
            "producer_id": producer,
            "input_snapshot_spec_sha256": spec_hash,
        }
    )


def build_candidate_manuscript(
    project_root: Path, fixture_mode: bool = False
) -> dict[str, object]:
    base = validate_snapshot_base(project_root, fixture_mode=fixture_mode)
    project = base["project"]
    candidate = base["candidate"]
    if not isinstance(project, Path) or not isinstance(candidate, Path):
        raise ReleaseContractError("internal REL-04 path type drift")
    load_frozen_input_index(project, candidate, str(base["spec_hash"]), fixture_mode)
    validate_rel02_state(candidate, str(base["spec_hash"]), fixture_mode)
    rel03_state = validate_rel03_state(candidate, str(base["spec_hash"]), fixture_mode)
    analysis = read_tsv_exact(
        candidate / "manifests/analysis_release_manifest.tsv",
        ANALYSIS_MANIFEST_FIELDS,
    )
    for row in analysis:
        if row["producer_id"] != "rel02_build_candidate_tables":
            continue
        path = candidate / row["candidate_path"]
        expected_bytes = parse_nonnegative_int(
            row["bytes"], f"REL-02 {row['artifact_id']} bytes"
        )
        if (
            not path.is_file()
            or path.is_symlink()
            or path.stat().st_size != expected_bytes
            or sha256_file(path) != row["sha256"]
            or row["input_snapshot_spec_sha256"] != base["spec_hash"]
        ):
            raise ReleaseContractError(
                f"REL-02 product failed pre-manuscript verification: {row['artifact_id']}"
            )
    figure_manifest = read_tsv_exact(
        candidate / "manifests/figure_source_manifest.tsv",
        FIGURE_SOURCE_MANIFEST_FIELDS,
    )
    for row in figure_manifest:
        pdf = resolve_project_path(
            project,
            row["pdf_path"],
            f"{row['figure_id']} {row['panel_id']} candidate PDF",
        )
        expected_bytes = parse_nonnegative_int(
            row["pdf_bytes"], f"{row['figure_id']} {row['panel_id']} PDF bytes"
        )
        if (
            not is_relative_to(pdf, candidate_figure_root(project))
            or not pdf.is_file()
            or pdf.is_symlink()
            or pdf.stat().st_size != expected_bytes
            or sha256_file(pdf) != row["pdf_sha256"]
        ):
            raise ReleaseContractError(
                "candidate PDF failed pre-manuscript verification: "
                f"{row['figure_id']} {row['panel_id']}"
            )
    numbers = read_tsv_exact(candidate / "tables/numbers_ledger.tsv", NUMBERS_FIELDS)
    claims = read_tsv_exact(candidate / "claims/claim_ledger.tsv", CLAIM_FIELDS)
    manuscript, used_numbers, used_claims = build_manuscript_text(numbers, claims)
    manuscript_root = ensure_candidate_output_root(project, "manuscript")
    manuscript_path = manuscript_root / "manuscript.md"
    materialize_immutable(
        manuscript_path,
        manuscript.encode("utf-8"),
        manuscript_root,
        "candidate manuscript",
    )
    manuscript_manifest_rows = [
        {
            "document_id": "candidate_manuscript",
            "candidate_path": project_relative_str(project, manuscript_path),
            "sha256": sha256_file(manuscript_path),
            "bytes": manuscript_path.stat().st_size,
            "journal_branch": JOURNAL_BRANCH,
            "figure_count": "5",
            "myojin_role": MYOJIN_ROLE,
            "number_ids": ";".join(sorted(used_numbers)),
            "claim_ids": ";".join(sorted(used_claims)),
        }
    ]
    manuscript_manifest = candidate / "manifests/manuscript_manifest.tsv"
    materialize_immutable_tsv(
        manuscript_manifest,
        manuscript_manifest_rows,
        MANUSCRIPT_MANIFEST_FIELDS,
        candidate,
        "manuscript manifest",
    )

    products: list[dict[str, object]] = []
    for row in analysis:
        if row["producer_id"] != "rel02_build_candidate_tables":
            continue
        path = candidate / row["candidate_path"]
        add_product(
            products,
            project,
            str(base["spec_hash"]),
            f"REL02:{row['artifact_id']}",
            "REL02",
            row["artifact_role"],
            path,
            row["producer_id"],
        )
    for name, role in (
        ("rel02_state.json", "rel02_state"),
        ("figure_source_manifest.tsv", "figure_source_manifest"),
        ("rel03_state.json", "rel03_state"),
        ("manuscript_manifest.tsv", "manuscript_manifest"),
    ):
        path = candidate / "manifests" / name
        phase = (
            "REL02"
            if name.startswith("rel02")
            else "REL03"
            if "figure" in name or name.startswith("rel03")
            else "REL04"
        )
        add_product(
            products,
            project,
            str(base["spec_hash"]),
            f"{phase}:{name}",
            phase,
            role,
            path,
            "rel04_build_candidate_manuscript" if phase == "REL04" else phase.lower(),
        )
    for row in figure_manifest:
        pdf = project / row["pdf_path"]
        add_product(
            products,
            project,
            str(base["spec_hash"]),
            f"REL03:{row['figure_id']}:{row['panel_id']}",
            "REL03",
            "individual_pdf_panel",
            pdf,
            "rel03_render_candidate_panels",
        )
    add_product(
        products,
        project,
        str(base["spec_hash"]),
        "REL04:candidate_manuscript",
        "REL04",
        "candidate_manuscript",
        manuscript_path,
        "rel04_build_candidate_manuscript",
    )
    if len({str(row["product_id"]) for row in products}) != len(products):
        raise ReleaseContractError("duplicate REL-02--04 transition product ID")
    product_manifest = candidate / "manifests/rel02_04_product_manifest.tsv"
    materialize_immutable_tsv(
        product_manifest,
        sorted(products, key=lambda row: str(row["product_id"])),
        TRANSITION_PRODUCT_FIELDS,
        candidate,
        "REL-02--04 product manifest",
    )
    transition = {
        "contract_version": 1,
        "candidate_id": CANDIDATE_ID,
        "transition": "REL02_04_PRODUCTS_FROZEN",
        "from_candidate_state_sha256": sha256_file(
            candidate / "manifests/candidate_state.json"
        ),
        "input_snapshot_spec_sha256": base["spec_hash"],
        "base_input_transition_sha256": sha256_file(
            candidate / "manifests/base_input_transition.json"
        ),
        "rel02_state_sha256": sha256_file(candidate / "manifests/rel02_state.json"),
        "rel03_state_sha256": sha256_file(candidate / "manifests/rel03_state.json"),
        "manuscript_manifest_sha256": sha256_file(manuscript_manifest),
        "product_manifest_sha256": sha256_file(product_manifest),
        "journal_branch": JOURNAL_BRANCH,
        "figure_count": 5,
        "myojin_role": MYOJIN_ROLE,
        "n_products": len(products),
        "n_panels": rel03_state["n_panels"],
        "fixture_mode": fixture_mode,
        "canonical_promotion_status": "not_promoted",
        "full_release_pass": False,
    }
    transition_path = candidate / "manifests/rel02_04_transition.json"
    materialize_immutable_json(
        transition_path,
        transition,
        candidate,
        "REL-02--04 transition",
    )
    figure_root = candidate_figure_root(project)
    if any(
        path.is_file() and path.suffix.lower() in {".png", ".svg"}
        for path in figure_root.rglob("*")
    ):
        raise ReleaseContractError(
            "non-PDF raster/vector substitute entered candidate figures"
        )
    return {
        "status": "REL04_MANUSCRIPT_AND_PRODUCTS_FROZEN",
        "candidate_id": CANDIDATE_ID,
        "journal_branch": JOURNAL_BRANCH,
        "figure_count": 5,
        "myojin_role": MYOJIN_ROLE,
        "n_products": len(products),
        "n_manuscript_numbers": len(used_numbers),
        "transition_sha256": sha256_file(transition_path),
        "canonical_promotion_status": "not_promoted",
        "full_release_pass": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--fixture-mode", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        result = build_candidate_manuscript(args.project_root, args.fixture_mode)
    except ReleaseContractError as error:
        raise SystemExit(f"REL04_BLOCKED: {error}") from error
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
