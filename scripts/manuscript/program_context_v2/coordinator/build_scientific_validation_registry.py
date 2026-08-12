#!/usr/bin/env python3
"""Derive the ten fail-closed scientific reports required before real REL-05.

This is a post-REL-04 validator. It reads only the immutable candidate snapshot,
candidate-only ledgers/panels/manuscript, and protected frozen inputs. A report
row is written only after every predicate has passed in memory. Expected
limitations are validated boundaries, not warnings; an unexpected condition
blocks instead of being auto-adjudicated.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable, Mapping, Sequence


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
if str(PROGRAM_ROOT) not in sys.path:
    sys.path.insert(0, str(PROGRAM_ROOT))

from coordinator_contract import (  # noqa: E402
    COMPOSITION_READY_FIELDS,
    COMPOSITION_RESULT_FIELDS,
    COMPOSITION_TESTABILITY_FIELDS,
    CORE_PRODUCER_FIELDS,
    EXPECTED_COMPOSITION_COUNTS,
    FIXED_HERO_GENES,
    PLAN50_ANALYSIS_RELEASE_ID,
    CoordinatorContractError,
    build_myojin_rows,
    build_spatial_rows,
    parse_bool,
    read_tsv_flexible,
    validate_composition_contract,
    validate_hotspot_registry_projection,
    validate_passport_negative_semantics,
)
from compare_clean_rebuilds import revalidate_clean_rebuild_report  # noqa: E402
from rel01_snapshot_candidate import (  # noqa: E402
    canonical_tsv_bytes,
    materialize_immutable_json,
    materialize_immutable_tsv,
)
from rel02_build_candidate_tables import (  # noqa: E402
    build_ledgers,
    load_all_rows,
    validate_blueprint,
)
from rel03_render_candidate_panels import build_pdf  # noqa: E402
from rel04_build_candidate_manuscript import (  # noqa: E402
    ABSTRACT_DATASET_NAMES,
    RETIRED_PHRASES,
    build_manuscript_text,
)
from rel05_validate_candidate import (  # noqa: E402
    validate_rel02_04_products,
    validate_snapshot_base,
)
from release_common import (  # noqa: E402
    CANDIDATE_ID,
    REQUIRED_SCIENTIFIC_CHECKS,
    SCIENTIFIC_REPORT_PRODUCER_PATHS,
    SCIENTIFIC_REGISTRY_FIELDS,
    ReleaseContractError,
    read_tsv_exact,
    sha256_file,
)
from release_products import (  # noqa: E402
    CLAIM_FIELDS,
    EXPECTED_FIGURE_TITLES,
    FIGURE_SOURCE_MANIFEST_FIELDS,
    MANUSCRIPT_MANIFEST_FIELDS,
    NUMBERS_FIELDS,
    SOURCE_ROW_FIELDS,
    candidate_figure_root,
    load_frozen_input_index,
)


REPORT_CONTRACT = "plan60_real_scientific_validation_v2"
EXPECTED_MAIN_PANELS = {
    ("Figure1", "1A"),
    ("Figure1", "1B"),
    ("Figure2", "2A"),
    ("Figure2", "2B"),
    ("Figure2", "2C"),
    ("Figure3", "3A"),
    ("Figure3", "3B"),
    ("Figure4", "4A"),
    ("Figure4", "4B"),
    ("Figure5", "5A"),
    ("Figure5", "5B"),
}
EXPECTED_SUPPLEMENT_PANELS = {
    ("FigureS1", "S1A"),
    ("FigureS1", "S1B"),
    ("FigureS2", "S2A"),
    ("FigureS2", "S2B"),
}


class ScientificValidationError(ReleaseContractError):
    """Raised before any scientific report or registry row is written."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ScientificValidationError(message)


def close(left: object, right: object, context: str, tolerance: float = 1e-10) -> None:
    try:
        a = float(str(left))
        b = float(str(right))
    except ValueError as error:
        raise ScientificValidationError(
            f"{context} is not numeric: {left!r}, {right!r}"
        ) from error
    if (
        not math.isfinite(a)
        or not math.isfinite(b)
        or not math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)
    ):
        raise ScientificValidationError(f"{context} drift: {a!r} != {b!r}")


def bh_adjust(pvalues: Sequence[float]) -> list[float]:
    require(bool(pvalues), "BH family is empty")
    require(
        all(math.isfinite(value) and 0 <= value <= 1 for value in pvalues),
        "invalid BH p value",
    )
    order = sorted(range(len(pvalues)), key=lambda index: pvalues[index])
    adjusted = [1.0] * len(pvalues)
    running = 1.0
    total = len(pvalues)
    for rank_index in range(total - 1, -1, -1):
        index = order[rank_index]
        rank = rank_index + 1
        running = min(running, pvalues[index] * total / rank, 1.0)
        adjusted[index] = running
    return adjusted


def parquet_records(path: Path) -> list[dict[str, object]]:
    try:
        import pyarrow.parquet as pq
    except ImportError as error:
        raise ScientificValidationError(
            "scientific Gene Catalog validation requires the existing pyarrow environment"
        ) from error
    return pq.read_table(path).to_pylist()


def frozen_by_role(
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
    workstream: str,
    role: str,
) -> Path:
    matches = [
        row
        for row in source_index.values()
        if row["workstream_id"] == workstream and row["artifact_role"] == role
    ]
    require(
        len(matches) == 1,
        f"{workstream} requires one frozen {role}; found {len(matches)}",
    )
    row = matches[0]
    path = candidate / row["snapshot_path"]
    require(
        path.is_file() and not path.is_symlink(),
        f"frozen {workstream}:{role} is unsafe",
    )
    require(
        sha256_file(path) == row["sha256"], f"frozen {workstream}:{role} hash drift"
    )
    return path


def frozen_by_artifact(
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
    artifact_id: str,
) -> Path:
    matches = [
        row for row in source_index.values() if row["artifact_id"] == artifact_id
    ]
    require(
        len(matches) == 1,
        f"requires one frozen artifact {artifact_id}; found {len(matches)}",
    )
    path = candidate / matches[0]["snapshot_path"]
    require(
        path.is_file() and sha256_file(path) == matches[0]["sha256"],
        f"frozen {artifact_id} drift",
    )
    return path


def all_candidate_source_rows(candidate: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in sorted((candidate / "figure_sources").rglob("*.tsv")):
        rows.extend(read_tsv_exact(path, SOURCE_ROW_FIELDS))
    require(bool(rows), "candidate figure-source tree is empty")
    record_ids = [row["record_id"] for row in rows]
    require(
        len(record_ids) == len(set(record_ids)),
        "candidate figure-source record IDs are not unique",
    )
    return rows


def panel_rows(
    rows: Sequence[dict[str, str]], figure: str, panel: str
) -> list[dict[str, str]]:
    selected = [
        row for row in rows if row["figure_id"] == figure and row["panel_id"] == panel
    ]
    require(bool(selected), f"missing candidate source rows for {figure}/{panel}")
    return selected


def scientific_producer_binding(
    project: Path,
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
) -> dict[str, object]:
    manifest_path = frozen_by_role(candidate, source_index, "BASE", "producer_manifest")
    producers = read_tsv_exact(manifest_path, CORE_PRODUCER_FIELDS)
    by_path = {row["repository_path"]: row for row in producers}
    require(
        all(path in by_path for path in SCIENTIFIC_REPORT_PRODUCER_PATHS),
        "scientific validator producer set is absent from the frozen producer manifest",
    )
    bindings = []
    for relative in SCIENTIFIC_REPORT_PRODUCER_PATHS:
        row = by_path[relative]
        current = project / relative
        require(
            current.is_file()
            and not current.is_symlink()
            and current.stat().st_size == int(row["bytes"])
            and sha256_file(current) == row["sha256"],
            f"scientific validator producer differs from frozen snapshot: {relative}",
        )
        bindings.append(
            {
                "producer_id": row["producer_id"],
                "repository_path": relative,
                "sha256": row["sha256"],
                "bytes": int(row["bytes"]),
            }
        )
    return {
        "producer_manifest_snapshot_path": manifest_path.relative_to(
            candidate
        ).as_posix(),
        "producer_manifest_sha256": sha256_file(manifest_path),
        "producer_bindings": bindings,
    }


def report_payload(
    check_id: str,
    spec_hash: str,
    producer_binding: Mapping[str, object],
    metrics: Mapping[str, object],
) -> dict[str, object]:
    return {
        "contract": REPORT_CONTRACT,
        "candidate_id": CANDIDATE_ID,
        "check_id": check_id,
        "status": "pass",
        "input_snapshot_spec_sha256": spec_hash,
        "validator_producer_binding": dict(producer_binding),
        "metrics": dict(metrics),
        "warnings": [],
        "canonical_promotion_authorized": False,
    }


def check_hotspot(
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
    rows: Sequence[dict[str, str]],
) -> dict[str, object]:
    _, registry = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN20", "program_registry")
    )
    _, semantic = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN20", "semantic_adjudication")
    )
    _, figure = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN20", "figure_source")
    )
    require(
        len(registry) == len(semantic) == len(figure) == 117,
        "Hotspot frozen family is not 117 modules",
    )
    registry_by_uid = {row["program_uid"]: row for row in registry}
    semantic_by_uid = {row["program_uid"]: row for row in semantic}
    figure_by_uid = {row["program_uid"]: row for row in figure}
    require(
        len(registry_by_uid) == len(semantic_by_uid) == len(figure_by_uid) == 117
        and set(registry_by_uid) == set(semantic_by_uid) == set(figure_by_uid),
        "Hotspot frozen identities disagree",
    )
    validate_hotspot_registry_projection(registry, figure)
    pvalues = [float(row["pvalue"]) for row in figure]
    for row, expected_q in zip(figure, bh_adjust(pvalues), strict=True):
        close(row["qvalue"], expected_q, f"Hotspot BH {row['program_uid']}")
    selected = {
        row["program_uid"]: row
        for row in semantic
        if parse_bool(
            row["robust_display"], f"Hotspot {row['program_uid']} robust_display"
        )
    }
    require(
        len(selected) == 2,
        "Hotspot semantic registry does not select exactly two robust programs",
    )
    require(
        all(
            row["adjudicated_state"] == "supported_internal_stage_association"
            for row in selected.values()
        ),
        "Hotspot robust display includes a non-supported adjudicated state",
    )
    source = panel_rows(rows, "Figure2", "2C")
    boundary_rows = [
        row
        for row in source
        if row["number_role"] == "leave_gse244832_direction_agreement_count"
    ]
    require(
        len(boundary_rows) == 1
        and boundary_rows[0]["value"] == "2"
        and boundary_rows[0]["denominator"] == "2"
        and "Only GSE244832 spans" in boundary_rows[0]["claim_text"],
        "Figure2 lacks the exact stage-by-cohort/LODO boundary",
    )
    program_source = [row for row in source if row not in boundary_rows]
    grouped: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in program_source:
        grouped[row["group_id"]][row["number_role"]] = row
    require(
        set(grouped) == set(selected),
        "Figure2 Hotspot programs differ from semantic robust selection",
    )
    for uid, semantic_row in selected.items():
        observed = grouped[uid]
        require(
            set(observed)
            == {"stage_ordinal_beta", "stage_ordinal_se", "robust_display_flag"},
            f"Figure2 Hotspot source lacks beta/SE/selection flag for {uid}",
        )
        upstream = figure_by_uid[uid]
        require(
            observed["stage_ordinal_beta"]["label"] == semantic_row["module_name"],
            f"Hotspot label drift {uid}",
        )
        close(
            observed["stage_ordinal_beta"]["value"],
            upstream["beta"],
            f"Hotspot beta {uid}",
        )
        close(
            observed["stage_ordinal_se"]["value"], upstream["se"], f"Hotspot SE {uid}"
        )
        close(
            observed["stage_ordinal_beta"]["p_value"],
            upstream["pvalue"],
            f"Hotspot p {uid}",
        )
        close(
            observed["stage_ordinal_beta"]["q_value"],
            upstream["qvalue"],
            f"Hotspot q {uid}",
        )
        require(
            observed["robust_display_flag"]["value"] == "1",
            f"Hotspot robust flag drift {uid}",
        )
        require(
            observed["stage_ordinal_beta"]["denominator"] == upstream["n_donors"],
            f"Hotspot donor n drift {uid}",
        )
    _, design = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN20", "stage_design_audit")
    )
    stage_sets: dict[str, set[str]] = defaultdict(set)
    for row in design:
        if row["program_uid"] in selected:
            stage_sets[row["dataset"]].add(row["stage_ordinal"])
    require(
        {
            dataset
            for dataset, stages in stage_sets.items()
            if {"0", "1", "2"}.issubset(stages)
        }
        == {"GSE244832"},
        "stage-by-dataset audit no longer has GSE244832 as the sole H/S/SH-spanning cohort",
    )
    _, lodo = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN20", "cohort_lodo_effects")
    )
    selected_lodo = [
        row
        for row in lodo
        if row["program_uid"] in selected
        and row["analysis_type"] == "lodo"
        and row["held_out_dataset"] == "GSE244832"
    ]
    require(
        len(selected_lodo) == 2
        and all(
            row["estimable"] == "TRUE" and row["direction"] == "positive"
            for row in selected_lodo
        ),
        "selected-program GSE244832 LODO direction boundary drift",
    )
    require(
        not any(row["evidence_status"] == "tested_negative" for row in source),
        "legacy Hotspot negative semantics leaked into Figure2",
    )
    return {
        "n_programs": 117,
        "n_robust_display": 2,
        "n_figure2_rows": len(source),
        "bh_family_size": 117,
        "semantic_source": "post-freeze adjudication, not legacy tested_negative",
        "sole_h_s_sh_spanning_cohort": "GSE244832",
        "n_selected_positive_without_gse244832": 2,
    }


def check_spatial(
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
    rows: Sequence[dict[str, str]],
) -> dict[str, object]:
    _, matrix = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN13", "figure_source")
    )
    require(
        len(matrix) == 18,
        "Plan13 frozen Figure4 matrix is not the complete 18-row family",
    )
    states = Counter(row["evidence_state"] for row in matrix)
    require(
        states
        == {
            "robust": 2,
            "indeterminate": 8,
            "skipped": 4,
            "not_applicable": 2,
            "untestable": 2,
        },
        f"Plan13 state drift: {states}",
    )
    robust = [row for row in matrix if row["evidence_state"] == "robust"]
    require(
        len(robust) == 2 and len({row["program_uid"] for row in robust}) == 1,
        "Plan13 robust row count is being misread as two distinct robust programs",
    )
    gse192741 = [row for row in matrix if row["dataset"] == "GSE192741"]
    vu = [row for row in matrix if row["dataset"] == "Vu_et_al_2025"]
    require(
        len(gse192741) == 2
        and all(
            row["source_dependence"] == "independent"
            and row["biological_unit_resolution"] == "resolved"
            and row["n_biological"] == "4"
            for row in gse192741
        ),
        "GSE192741 four-donor independent spatial boundary drift",
    )
    require(
        len(vu) == 2
        and all(
            row["source_dependence"] == "source_dependent"
            and row["biological_unit_resolution"] == "unresolved"
            and not row["n_biological"]
            and row["n_technical"] == "10"
            for row in vu
        ),
        "Vu ten-section source-dependent/unresolved-donor boundary drift",
    )
    require(
        not any(row["negative_call_rule_id"] for row in matrix),
        "Plan13 unexpectedly carries a negative rule",
    )
    source = [row for row in rows if row["figure_id"] == "Figure4"]
    marks = {row["record_id"]: row for row in source if row["plot_role"] == "mark"}
    require(len(marks) == 18, "Figure4 does not retain all 18 prespecified matrix rows")
    status_map = {
        "robust": "context_supported",
        "indeterminate": "indeterminate",
        "skipped": "skipped",
        "not_applicable": "not_applicable",
        "untestable": "untestable",
    }
    uncertainty_columns = {
        "std_error": "assay_native_standard_error",
        "matched_null_sd": "matched_null_standard_deviation",
        "interval_low": "assay_native_interval_lower",
        "interval_high": "assay_native_interval_upper",
    }
    annotations = {
        row["record_id"]: row for row in source if row["plot_role"] == "annotation"
    }
    technical_annotations = [
        row for row in source if row["number_role"] == "technical_unit_count"
    ]
    require(
        len(technical_annotations) == 18
        and all(
            "not a biological denominator" in row["unit"]
            for row in technical_annotations
        ),
        "Figure4 technical-unit ledger is incomplete or presented as biological n",
    )
    for order, upstream in enumerate(matrix, 1):
        observed = marks.get(f"spatial_{order}")
        require(observed is not None, f"Figure4 lost matrix row {order}")
        expected_value = upstream["estimate"] or "1"
        close(observed["value"], expected_value, f"Figure4 estimate/gate row {order}")
        require(
            observed["evidence_status"] == status_map[upstream["evidence_state"]],
            f"Figure4 state drift row {order}",
        )
        expected_denominator = (
            upstream["n_biological"]
            if upstream["biological_unit_resolution"] == "resolved"
            else ""
        )
        require(
            observed["denominator"] == expected_denominator,
            f"Figure4 biological denominator drift row {order}",
        )
        if upstream["estimate"]:
            require(
                observed["effect_unit"] == upstream["effect_unit"],
                f"Figure4 native unit drift row {order}",
            )
            close(observed["p_value"], upstream["pvalue"], f"Figure4 p row {order}")
            close(observed["q_value"], upstream["padj"], f"Figure4 q row {order}")
        for column, number_role in uncertainty_columns.items():
            key = f"spatial_{order}_{number_role}"
            if upstream[column]:
                require(key in annotations, f"Figure4 lost {column} row {order}")
                close(
                    annotations[key]["value"],
                    upstream[column],
                    f"Figure4 {column} row {order}",
                )
            else:
                require(
                    key not in annotations, f"Figure4 fabricated {column} row {order}"
                )
    direction_rows = [
        row for row in source if row["record_id"].startswith("yak_module8_")
    ]
    direction_roles = {row["number_role"] for row in direction_rows}
    require(
        direction_roles
        == {
            "descriptive_median_slope_direction",
            "inferential_signed_stouffer_direction",
        },
        "Yakubovsky module 8 descriptive/inferential directions were collapsed",
    )
    require(
        len(direction_rows) == 2
        and all(row["group_id"] == "yak_module8" for row in direction_rows),
        "Yakubovsky module 8 direction annotation identity drift",
    )
    panel4b = panel_rows(rows, "Figure4", "4B")
    rendered4b = build_pdf("Zonation-adjusted lipid challenge", panel4b)
    for expected_text in (
        b"descriptive_median_slope_direction=positive",
        b"inferential_signed_stouffer_direction=negative",
    ):
        require(
            expected_text in rendered4b,
            f"Figure4B renderer silently drops {expected_text.decode('ascii')}",
        )
    require(
        not any("cross_assay" in row["number_role"].lower() for row in source),
        "Figure4 contains a cross-assay summary field",
    )
    return {
        "n_prespecified_matrix_rows": 18,
        "n_candidate_rows_including_uncertainty": len(source),
        "evidence_state_counts": dict(sorted(states.items())),
        "n_tested_negative": 0,
        "n_unique_robust_programs": 1,
        "gse192741_independent_donors": 4,
        "vu_technical_sections": 10,
        "vu_biological_donors": "unresolved",
        "technical_unit_annotations": len(technical_annotations),
        "yak_module8_direction_roles": sorted(direction_roles),
    }


def check_numbers_claims(
    candidate: Path, product_metrics: Mapping[str, object]
) -> dict[str, object]:
    numbers = read_tsv_exact(candidate / "tables/numbers_ledger.tsv", NUMBERS_FIELDS)
    claims = read_tsv_exact(candidate / "claims/claim_ledger.tsv", CLAIM_FIELDS)
    require(
        len(numbers) == int(product_metrics["n_numbers"]),
        "numbers-ledger/core validation count drift",
    )
    require(
        len(claims) == int(product_metrics["n_claims"]),
        "claim-ledger/core validation count drift",
    )
    require(
        len({row["number_id"] for row in numbers}) == len(numbers),
        "duplicate number IDs",
    )
    require(
        len({row["claim_id"] for row in claims}) == len(claims), "duplicate claim IDs"
    )
    require(
        all(row["source_artifact"].startswith("inputs/") for row in numbers),
        "a number reads a mutable source",
    )
    return {
        "n_numbers": len(numbers),
        "n_claims": len(claims),
        "n_manuscript_numbers": int(product_metrics["n_manuscript_numbers"]),
        "source_paths_candidate_snapshot_only": True,
    }


def check_panels(
    project: Path,
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
    rows: Sequence[dict[str, str]],
) -> dict[str, object]:
    manifest = read_tsv_exact(
        candidate / "manifests/figure_source_manifest.tsv",
        FIGURE_SOURCE_MANIFEST_FIELDS,
    )
    identities = {(row["figure_id"], row["panel_id"]) for row in manifest}
    require(
        identities == EXPECTED_MAIN_PANELS | EXPECTED_SUPPLEMENT_PANELS,
        f"panel manifest coverage drift: {identities}",
    )
    require(
        len(manifest) == 15,
        "panel manifest must contain 11 main plus 4 supplementary panels",
    )
    for row in manifest:
        pdf = project / row["pdf_path"]
        source = candidate / row["source_table_path"]
        require(
            pdf.is_file() and pdf.stat().st_size > 0 and pdf.suffix.lower() == ".pdf",
            f"missing PDF {row['figure_id']}/{row['panel_id']}",
        )
        require(
            sha256_file(pdf) == row["pdf_sha256"],
            f"PDF hash drift {row['figure_id']}/{row['panel_id']}",
        )
        require(
            source.is_file() and sha256_file(source) == row["source_table_sha256"],
            f"panel source drift {row['figure_id']}/{row['panel_id']}",
        )
    figure_tree = {
        path for path in candidate_figure_root(project).rglob("*") if path.is_file()
    }
    require(
        not any(path.suffix.lower() in {".png", ".svg"} for path in figure_tree),
        "candidate figure tree contains PNG/SVG panels",
    )
    myojin = [row for row in rows if row["figure_id"] == "FigureS1"]
    require(bool(myojin), "Myojin supplementary source rows are missing")
    require(
        all(
            row["evidence_status"] in {"indeterminate", "untestable"} for row in myojin
        ),
        "Myojin supplement overstates a call",
    )
    require(
        not any(row["manuscript_included"] == "true" for row in myojin),
        "Myojin supplementary rows leaked into main manuscript inclusion",
    )
    class_rows = [row for row in myojin if row["panel_id"] == "S1A"]
    omnibus = [
        row for row in class_rows if row["number_role"] == "class_omnibus_partial_f"
    ]
    pairwise = [
        row
        for row in class_rows
        if row["number_role"] == "class_pairwise_standardized_effect"
    ]
    require(
        len(class_rows) == 28
        and len(omnibus) == 7
        and len(pairwise) == 21
        and all(row["plot_role"] == "annotation" for row in omnibus)
        and all(row["plot_role"] == "mark" for row in pairwise),
        "Myojin FigureS1A mixes omnibus F statistics with standardized-effect marks",
    )
    rendered_s1a = build_pdf("Complete Myojin evidence-class stress test", class_rows)
    require(
        b"class_omnibus_partial_f=" in rendered_s1a,
        "Myojin omnibus annotation is absent from the rendered supplementary PDF",
    )
    nmf = [row for row in rows if row["figure_id"] == "FigureS2"]
    s2a = [row for row in nmf if row["panel_id"] == "S2A"]
    s2b = [row for row in nmf if row["panel_id"] == "S2B"]
    require(
        len(nmf) == 18
        and len(s2a) == 10
        and len(s2b) == 8
        and all(row["plot_role"] == "mark" for row in nmf)
        and all(row["section_id"] == "supplementary_nmf" for row in nmf)
        and all(row["manuscript_included"] == "false" for row in nmf),
        "FigureS2 must contain 10 continuous-axis marks and 8 stability marks",
    )
    require(
        b"median=" in build_pdf("Continuous k4/k6 NMF loading distributions", s2a)
        and b"mean Hungaria"
        in build_pdf("Three-seed factor stability and hard-partition boundary", s2b),
        "FigureS2 renderers omit the continuous-loading or stability evidence",
    )
    nmf_loadings = read_tsv_exact(
        frozen_by_role(candidate, source_index, "PLAN20", "supplement_loadings"),
        (
            "sample_id",
            "k",
            "program_code",
            "program_label",
            "continuous_loading",
            "source_sha256",
        ),
    )
    require(
        len(nmf_loadings) == 11040
        and len({row["sample_id"] for row in nmf_loadings}) == 1104
        and {row["k"] for row in nmf_loadings} == {"4", "6"},
        "FigureS2 frozen continuous-loading universe drift",
    )
    return {
        "n_main_panels": len(EXPECTED_MAIN_PANELS),
        "n_supplementary_panels": len(EXPECTED_SUPPLEMENT_PANELS),
        "myojin_class_omnibus_annotations": len(omnibus),
        "myojin_class_pairwise_effect_marks": len(pairwise),
        "n_pdf_panels": len(manifest),
        "myojin_supplement_rows": len(myojin),
        "nmf_supplement_rows": len(nmf),
        "nmf_continuous_loading_rows": len(nmf_loadings),
    }


def passport_inputs(
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    genes = parquet_records(
        frozen_by_role(candidate, source_index, "PLAN50", "gene_index")
    )
    evidence = parquet_records(
        frozen_by_role(candidate, source_index, "PLAN50", "gene_evidence")
    )
    coverage = parquet_records(
        frozen_by_role(candidate, source_index, "PLAN50", "gene_coverage")
    )
    contexts = parquet_records(
        frozen_by_role(candidate, source_index, "PLAN50", "program_context")
    )
    _, experiments = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN50", "next_experiment")
    )
    _, source_nodes = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN50", "source_nodes")
    )
    _, source_edges = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN50", "source_edges")
    )
    return genes, evidence, coverage, contexts, experiments, source_nodes, source_edges


def check_passports(
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
    rows: Sequence[dict[str, str]],
) -> dict[str, object]:
    (genes, evidence, coverage, contexts, experiments, source_nodes, source_edges) = (
        passport_inputs(candidate, source_index)
    )
    for name, records in (
        ("gene", genes),
        ("evidence", evidence),
        ("coverage", coverage),
        ("program_context", contexts),
        ("next_experiment", experiments),
    ):
        require(
            all(
                str(row.get("analysis_release_id", "")) == PLAN50_ANALYSIS_RELEASE_ID
                for row in records
            ),
            f"Gene Catalog {name} table has nested analysis-release drift",
        )
    validate_passport_negative_semantics(evidence)
    validate_passport_negative_semantics(contexts)
    nodes = {row["source_node_id"]: row for row in source_nodes}
    edges = {row["source_edge_id"]: row for row in source_edges}
    require(
        len(nodes) == len(source_nodes)
        and len(edges) == len(source_edges)
        and nodes
        and edges,
        "Gene Catalog source graph has duplicate or empty node/edge IDs",
    )
    for edge_id, edge in edges.items():
        require(
            edge["from_node_id"] in nodes and edge["to_node_id"] in nodes,
            f"Gene Catalog source edge references an unknown node: {edge_id}",
        )
    for result in evidence:
        result_id = str(result["evidence_result_id"])
        require(
            str(result.get("source_node_id") or "") in nodes,
            f"Gene Catalog entry lacks a graph-resolved source node: {result_id}",
        )
        require(
            any(
                edge.get("edge_type") == "tested_by"
                and edge.get("evidence_result_id") == result_id
                for edge in source_edges
            ),
            f"Gene Catalog entry lacks an exact tested_by edge: {result_id}",
        )
    source_releases = {
        str(node.get("source_release_id") or "") for node in source_nodes
    }
    require(
        all(
            str(row.get("source_release_id") or "") in source_releases
            for row in coverage
        ),
        "Gene Catalog coverage references a release absent from the source graph",
    )
    require(
        all(row["gene_call_expansion_authorized"] is False for row in contexts),
        "program-to-gene call expansion is authorized",
    )
    forbidden_columns = {
        "score",
        "rank",
        "priority_score",
        "combined_score",
        "overall_score",
        "modality_count",
        "evidence_count",
        "target_probability",
    }
    for name, records in (
        ("gene", genes),
        ("evidence", evidence),
        ("coverage", coverage),
        ("program_context", contexts),
    ):
        fields = set(records[0]) if records else set()
        require(
            not (fields & forbidden_columns),
            f"Gene Catalog {name} table contains score/rank columns: {fields & forbidden_columns}",
        )
    by_symbol = defaultdict(list)
    for row in genes:
        by_symbol[str(row["symbol"])].append(row)
    hero_passports = []
    expected_evidence = 0
    expected_coverage = 0
    for symbol in FIXED_HERO_GENES:
        require(
            len(by_symbol[symbol]) == 1, f"hero symbol identity is ambiguous: {symbol}"
        )
        gene = by_symbol[symbol][0]
        require(gene["symbol_collision"] is False, f"hero symbol collision: {symbol}")
        require(
            re.fullmatch(r"ENSG\d+(?:\.\d+)?", str(gene["ensembl_id"])) is not None,
            f"hero is not Ensembl-first: {symbol}",
        )
        passport_id = str(gene["passport_id"])
        hero_passports.append(passport_id)
        expected_evidence += sum(
            str(row["passport_id"]) == passport_id for row in evidence
        )
        expected_coverage += sum(
            str(row["passport_id"]) == passport_id for row in coverage
        )
        require(
            sum(row["passport_id"] == passport_id for row in experiments) == 1,
            f"hero lacks one next experiment: {symbol}",
        )
    figure5 = [row for row in rows if row["figure_id"] == "Figure5"]
    require(
        len(figure5)
        == 2 * len(FIXED_HERO_GENES) + expected_evidence + expected_coverage,
        "Figure5 does not preserve every hero evidence/coverage row",
    )
    require(
        not any(
            "score" in row["number_role"].lower()
            or "rank" in row["number_role"].lower()
            for row in figure5
        ),
        "Figure5 hides a score/rank number role",
    )
    graph_ids = set(nodes) | set(edges)
    for row in figure5:
        discovery = set(row["discovery_sources"].split(";"))
        evaluation = set(row["evaluation_sources"].split(";"))
        require(
            row["source_dependence"] == "source_dependent"
            and "figure5_vignette_selection_posthoc" in discovery
            and discovery - {"figure5_vignette_selection_posthoc"} == evaluation
            and evaluation
            and evaluation.issubset(graph_ids),
            f"Figure5 row is not bound to exact Plan50 source-graph IDs: {row['record_id']}",
        )
    return {
        "analysis_release_id": PLAN50_ANALYSIS_RELEASE_ID,
        "n_gene_passports": len(genes),
        "n_gene_evidence_rows": len(evidence),
        "n_coverage_rows": len(coverage),
        "n_program_context_rows": len(contexts),
        "hero_genes": list(FIXED_HERO_GENES),
        "hero_figure_rows": len(figure5),
        "program_to_gene_call_expansion_authorized": False,
        "source_graph_nodes": len(nodes),
        "source_graph_edges": len(edges),
        "figure5_source_graph_resolved": True,
    }


def check_determinism(
    project: Path,
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
    fixture_mode: bool,
) -> dict[str, object]:
    blueprint, _ = validate_blueprint(candidate, dict(source_index))
    first = load_all_rows(candidate, blueprint, dict(source_index))
    second = load_all_rows(candidate, blueprint, dict(source_index))
    require(
        canonical_tsv_bytes(first, SOURCE_ROW_FIELDS)
        == canonical_tsv_bytes(second, SOURCE_ROW_FIELDS),
        "two in-memory REL02 source rebuilds differ",
    )
    for rows in (first, second):
        require(
            len(rows) == len({str(row["record_id"]) for row in rows}),
            "determinism rebuild has duplicate records",
        )
    first_ledgers = build_ledgers(first)
    second_ledgers = build_ledgers(second)
    require(first_ledgers == second_ledgers, "two in-memory ledger rebuilds differ")
    manifest = read_tsv_exact(
        candidate / "manifests/figure_source_manifest.tsv",
        FIGURE_SOURCE_MANIFEST_FIELDS,
    )
    pdf_count = 0
    for row in manifest:
        source_rows = read_tsv_exact(
            candidate / row["source_table_path"], SOURCE_ROW_FIELDS
        )
        one = build_pdf(row["panel_title"], source_rows)
        two = build_pdf(row["panel_title"], source_rows)
        require(
            one == two, f"two PDF rebuilds differ: {row['figure_id']}/{row['panel_id']}"
        )
        pdf = project / row["pdf_path"]
        require(
            one == pdf.read_bytes(),
            f"recomputed PDF differs from frozen panel: {row['figure_id']}/{row['panel_id']}",
        )
        pdf_count += 1
    numbers = read_tsv_exact(candidate / "tables/numbers_ledger.tsv", NUMBERS_FIELDS)
    claims = read_tsv_exact(candidate / "claims/claim_ledger.tsv", CLAIM_FIELDS)
    manuscript_one, _, _ = build_manuscript_text(numbers, claims)
    manuscript_two, _, _ = build_manuscript_text(numbers, claims)
    require(manuscript_one == manuscript_two, "two manuscript rebuilds differ")
    manuscript_manifest = read_tsv_exact(
        candidate / "manifests/manuscript_manifest.tsv", MANUSCRIPT_MANIFEST_FIELDS
    )
    require(len(manuscript_manifest) == 1, "manuscript manifest is not singular")
    manuscript_path = project / manuscript_manifest[0]["candidate_path"]
    require(
        manuscript_one.encode("utf-8") == manuscript_path.read_bytes(),
        "recomputed manuscript differs from frozen candidate",
    )
    clean_report = candidate / "manifests/two_clean_rebuild_comparison.tsv"
    if fixture_mode:
        clean_metrics = {
            "retained_clean_rebuild_status": "fixture_mode_exempt",
            "retained_clean_rebuild_products": 0,
            "retained_clean_rebuild_report_sha256": "",
        }
    else:
        clean_rows = revalidate_clean_rebuild_report(project, clean_report)
        clean_metrics = {
            "retained_clean_rebuild_status": "authority_and_two_retained_builds_identical",
            "retained_clean_rebuild_products": len(clean_rows),
            "retained_clean_rebuild_report_sha256": sha256_file(clean_report),
            "retained_comparison_root": clean_rows[0]["comparison_root"],
            "build_a_project_root": clean_rows[0]["build_a_project_root"],
            "build_b_project_root": clean_rows[0]["build_b_project_root"],
        }
    return {
        "validation_scope": "in_process_repeatability_plus_retained_clean_rebuilds",
        "in_process_source_rebuilds": 2,
        "in_process_ledger_rebuilds": 2,
        "in_process_pdf_rebuilds_per_panel": 2,
        "n_pdf_panels": pdf_count,
        "in_process_manuscript_rebuilds": 2,
        "fixture_mode": fixture_mode,
        **clean_metrics,
    }


def check_null_compatibility(
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
    rows: Sequence[dict[str, str]],
) -> dict[str, object]:
    _, matrix = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN13", "figure_source")
    )
    synthetic_null = copy.deepcopy(matrix)
    for row in synthetic_null:
        if row["evidence_state"] in {"robust", "indeterminate"}:
            row["evidence_state"] = "indeterminate"
            row["robustness_pass"] = "FALSE"
            row["negative_call_rule_id"] = ""
    null_rows = build_spatial_rows(synthetic_null)
    require(
        not any(row["evidence_status"] == "tested_negative" for row in null_rows),
        "all-null spatial fixture fabricates tested-negative calls",
    )
    _, classes = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN40", "class_effects")
    )
    _, programs = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN40", "program_effects")
    )
    myojin = build_myojin_rows(classes, programs)
    require(
        all(
            row["evidence_status"] in {"indeterminate", "untestable"} for row in myojin
        ),
        "nonconfirmatory Myojin fixture fails",
    )
    figure5_main = [row for row in rows if row["figure_id"] == "Figure5"]
    require(
        bool(figure5_main)
        and not any(
            row["section_id"] == "supplementary_myojin" for row in figure5_main
        ),
        "null Myojin branch displaces the Gene Catalog in Figure 5",
    )
    return {
        "synthetic_all_nonpositive_spatial_rows": len(null_rows),
        "myojin_nonconfirmatory_rows": len(myojin),
        "myojin_main_figure_eligible": False,
        "passport_figure5_retained": True,
    }


def check_manuscript(
    project: Path, candidate: Path, rows: Sequence[dict[str, str]]
) -> dict[str, object]:
    manifest = read_tsv_exact(
        candidate / "manifests/manuscript_manifest.tsv", MANUSCRIPT_MANIFEST_FIELDS
    )
    require(len(manifest) == 1, "manuscript manifest is not singular")
    path = project / manifest[0]["candidate_path"]
    text = path.read_text(encoding="utf-8")
    lowered = text.lower()
    require(
        not any(phrase in lowered for phrase in RETIRED_PHRASES),
        "candidate manuscript contains a retired positive claim",
    )
    abstract = text.split("## Abstract", 1)[1].split("## Introduction", 1)[0]
    require(
        not any(name.lower() in abstract.lower() for name in ABSTRACT_DATASET_NAMES),
        "abstract names a processed external dataset",
    )
    require(
        EXPECTED_FIGURE_TITLES["Figure4"] in text,
        "Figure4 exact directional title is absent",
    )
    require(
        "Figure 6" not in text and "Figure6" not in text,
        "five-figure manuscript contains Figure6",
    )
    require(
        not any("broad validation" in row["claim_text"].lower() for row in rows),
        "figure source contains broad-validation overclaim",
    )
    return {
        "retired_phrases_checked": list(RETIRED_PHRASES),
        "abstract_dataset_names_excluded": list(ABSTRACT_DATASET_NAMES),
        "main_figure_count": 5,
        "figure4_title": EXPECTED_FIGURE_TITLES["Figure4"],
    }


def composition_inputs(
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
) -> tuple[
    list[dict[str, str]],
    list[dict[str, str]],
    dict[str, str],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    result_path = frozen_by_role(
        candidate, source_index, "PLAN20", "sample_composition"
    )
    testability_path = frozen_by_role(
        candidate, source_index, "PLAN20", "sample_composition_testability"
    )
    ready_path = frozen_by_role(
        candidate, source_index, "PLAN20", "composition_ready_seal"
    )
    audit_path = frozen_by_role(
        candidate, source_index, "PLAN20", "sample_composition_audit"
    )
    validation_path = frozen_by_role(
        candidate, source_index, "PLAN20", "sample_composition_validation"
    )
    results = read_tsv_exact(result_path, COMPOSITION_RESULT_FIELDS)
    testability = read_tsv_exact(testability_path, COMPOSITION_TESTABILITY_FIELDS)
    ready_rows = read_tsv_exact(ready_path, COMPOSITION_READY_FIELDS)
    require(len(ready_rows) == 1, "composition READY seal is not singular")
    ready = ready_rows[0]
    audit = read_tsv_exact(audit_path, ("audit_id", "status", "value", "detail"))
    validation = read_tsv_exact(
        validation_path, ("check_id", "status", "observed", "criterion")
    )
    for hash_field, path in (
        ("results_sha256", result_path),
        ("testability_sha256", testability_path),
        ("audit_sha256", audit_path),
        ("validation_sha256", validation_path),
    ):
        require(
            ready[hash_field] == sha256_file(path),
            f"composition READY {hash_field} differs from frozen bytes",
        )
    validate_composition_contract(results, testability, ready, audit, validation)
    return results, testability, ready, audit, validation


def check_composition_source(
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
    rows: Sequence[dict[str, str]],
) -> dict[str, object]:
    results, testability, ready, _, _ = composition_inputs(candidate, source_index)
    source = panel_rows(rows, "Figure2", "2B")
    marks = [row for row in source if row["plot_role"] == "mark"]
    summaries = [
        row
        for row in source
        if row["number_role"] == "bh_significant_testable_celltype_count"
    ]
    require(
        len(source) == 23 and len(marks) == 22 and len(summaries) == 1,
        "Figure2B must retain 22 cell types plus one family-summary annotation",
    )
    by_label = {row["label"]: row for row in marks}
    require(len(by_label) == 22, "Figure2B composition labels are not unique")
    testability_by_label = {row["celltype"]: row for row in testability}
    result_by_label = {row["celltype"]: row for row in results}
    require(
        set(by_label) == set(testability_by_label),
        "Figure2B differs from the frozen 22-cell-type testability universe",
    )
    for celltype, testability_row in testability_by_label.items():
        observed = by_label[celltype]
        require(
            observed["panel_title"] == "QC-passing sample-level cell composition",
            f"Figure2B panel title drift: {celltype}",
        )
        claim_language = " ".join(
            observed[field].lower()
            for field in (
                "panel_title",
                "allowed_wording",
                "claim_text",
                "biological_unit",
            )
        )
        require(
            "donor" not in claim_language, f"Figure2B donor wording leak: {celltype}"
        )
        if testability_row["testability"] == "testable":
            native = result_by_label[celltype]
            require(
                observed["number_role"] == "composition_effect_arcsin_sqrt"
                and observed["effect_unit"]
                == "arcsine-square-root proportion difference",
                f"Figure2B native effect unit drift: {celltype}",
            )
            close(
                observed["value"],
                native["effect_arcsin_sqrt"],
                f"composition effect {celltype}",
            )
            close(observed["p_value"], native["pvalue"], f"composition p {celltype}")
            close(observed["q_value"], native["padj"], f"composition q {celltype}")
            require(
                observed["denominator"] == native["n_analyzed"] == "1221",
                f"Figure2B sample denominator drift: {celltype}",
            )
            expected_status = (
                "state_associated" if float(native["padj"]) < 0.05 else "indeterminate"
            )
            require(
                observed["evidence_status"] == expected_status,
                f"Figure2B call-state drift: {celltype}",
            )
        else:
            require(
                celltype not in result_by_label
                and observed["number_role"] == "composition_testability_indicator"
                and observed["value"] == "1"
                and observed["evidence_status"] == "untestable"
                and not observed["p_value"]
                and not observed["q_value"]
                and observed["effect_unit"] == "not applicable",
                f"Figure2B structural-unavailability semantics drift: {celltype}",
            )
            require(
                "not an effect" in observed["unit"],
                f"Figure2B untestable indicator could be mistaken for an effect: {celltype}",
            )
    n_q_significant = sum(float(row["padj"]) < 0.05 for row in results)
    summary = summaries[0]
    require(
        summary["value"] == str(n_q_significant)
        and summary["numerator"] == str(n_q_significant)
        and summary["denominator"] == "16"
        and summary["manuscript_included"] == "true"
        and summary["evidence_status"] == "descriptive",
        "Figure2B manuscript family summary is absent or numerically inconsistent",
    )
    return {
        **EXPECTED_COMPOSITION_COUNTS,
        "n_figure2b_rows": len(source),
        "n_figure2b_celltype_marks": len(marks),
        "n_q_significant": n_q_significant,
        "model": ready["model"],
        "effect_unit": "arcsine-square-root proportion difference",
        "biological_unit": "human liver sample",
    }


def check_biological_units(
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
    rows: Sequence[dict[str, str]],
) -> dict[str, object]:
    composition = check_composition_source(candidate, source_index, rows)
    hotspot = [
        row
        for row in rows
        if row["figure_id"] == "Figure2"
        and row["panel_id"] == "2C"
        and row["number_role"] == "stage_ordinal_beta"
    ]
    require(
        len(hotspot) == 2
        and all(row["biological_unit"] == "human biological donor" for row in hotspot),
        "Hotspot replication unit is not donor",
    )
    _, matrix = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN13", "figure_source")
    )
    spatial_marks = {
        row["record_id"]: row
        for row in rows
        if row["figure_id"] == "Figure4" and row["plot_role"] == "mark"
    }
    for order, upstream in enumerate(matrix, 1):
        observed = spatial_marks[f"spatial_{order}"]
        if upstream["biological_unit_resolution"] == "unresolved":
            require(
                observed["denominator"] == "",
                f"unresolved Figure4 row exposes biological n: {order}",
            )
            require(
                "donor unknown" in observed["biological_unit"].lower(),
                f"unresolved Figure4 row lacks donor-unknown label: {order}",
            )
        else:
            require(
                observed["denominator"] == upstream["n_biological"],
                f"Figure4 donor/biological n drift: {order}",
            )
        require(
            not any(
                token in observed["biological_unit"].lower()
                for token in ("spot", "aoi", "section", "array")
            ),
            f"technical unit used as biological unit: {order}",
        )
    gse287826 = [
        row for row in spatial_marks.values() if row["label"].startswith("GSE287826:")
    ]
    require(
        len(gse287826) == 2
        and all(
            row["evidence_status"] == "skipped" and not row["denominator"]
            for row in gse287826
        ),
        "GSE287826 skip/donor-key boundary drift",
    )
    census = {
        row["number_role"]: row["value"]
        for row in rows
        if row["figure_id"] == "Figure1" and row["panel_id"] == "1A"
    }
    require(
        census
        == {
            "qc_samples": "1260",
            "qc_cohorts": "9",
            "pooled_samples": "846",
            "pooled_cohorts": "5",
        },
        f"Figure1 census drift: {census}",
    )
    bbj = [
        row
        for row in rows
        if row["number_role"] == "bbj_eas_gwas_with_eur_liver_eqtl_count"
    ]
    overlap = [row for row in rows if row["number_role"] == "fisher_overlap_odds_ratio"]
    require(
        len(bbj) == 1
        and bbj[0]["value"] == "3"
        and "EAS" in bbj[0]["allowed_wording"]
        and "EUR" in bbj[0]["allowed_wording"]
        and "cross-ancestry causality" in bbj[0]["prohibited_wording"],
        "BBJ EAS-GWAS/EUR-eQTL ancestry boundary drift",
    )
    require(len(overlap) == 1, "genetic/state overlap audit row is not singular")
    close(overlap[0]["value"], 0.889461813691, "genetic/state overlap OR")
    close(overlap[0]["p_value"], 0.603934955413, "genetic/state overlap Fisher p")
    require(
        overlap[0]["numerator"] == "34"
        and overlap[0]["denominator"] == "447"
        and overlap[0]["evidence_status"] == "indeterminate"
        and "non-enriched" in overlap[0]["allowed_wording"],
        "genetic/state overlap is presented as enriched or with wrong counts",
    )
    return {
        "sample_composition": composition,
        "hotspot_donor_programs": 2,
        "figure4_prespecified_rows": len(matrix),
        "gse287826_status": "skipped_no_donor_key",
        "qc_samples": 1260,
        "pooled_samples": 846,
        "bbj_eas_gwas_eur_liver_eqtl_strata": 3,
        "genetic_state_overlap": "34/447; OR=0.889461813691; Fisher p=0.603934955413",
    }


def check_multiplicity(
    candidate: Path,
    source_index: Mapping[str, Mapping[str, str]],
) -> dict[str, object]:
    _, hotspot = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN20", "figure_source")
    )
    pvalues = [float(row["pvalue"]) for row in hotspot]
    for row, expected in zip(hotspot, bh_adjust(pvalues), strict=True):
        close(row["qvalue"], expected, f"Hotspot multiplicity {row['program_uid']}")
    composition, _, _, _, _ = composition_inputs(candidate, source_index)
    composition_ordered = sorted(composition, key=lambda row: row["celltype"])
    composition_adjusted = bh_adjust(
        [float(row["pvalue"]) for row in composition_ordered]
    )
    for row, expected in zip(composition_ordered, composition_adjusted, strict=True):
        close(row["padj"], expected, f"composition multiplicity {row['celltype']}")
    _, matrix = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN13", "figure_source")
    )
    spatial_families: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in matrix:
        if row["pvalue"]:
            spatial_families[(row["dataset"], row["analysis_set_id"])].append(row)
    for family, family_rows in spatial_families.items():
        adjusted = bh_adjust([float(row["pvalue"]) for row in family_rows])
        for row, expected in zip(family_rows, adjusted, strict=True):
            close(
                row["padj"],
                expected,
                f"spatial multiplicity {family}/{row['program_uid']}",
            )
    _, program_rows = read_tsv_flexible(
        frozen_by_role(candidate, source_index, "PLAN40", "program_effects")
    )
    myojin_families: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in program_rows:
        if row["empirical_p"]:
            myojin_families[(row["analysis_variant"], row["weight_mode"])].append(row)
    for family, family_rows in myojin_families.items():
        adjusted = bh_adjust([float(row["empirical_p"]) for row in family_rows])
        for row, expected in zip(family_rows, adjusted, strict=True):
            close(
                row["BH_q"],
                expected,
                f"Myojin multiplicity {family}/{row['program_uid']}",
            )
    return {
        "hotspot_bh_family_size": len(hotspot),
        "composition_bh_family_size": len(composition_ordered),
        "spatial_bh_families": len(spatial_families),
        "spatial_bh_rows": sum(len(rows) for rows in spatial_families.values()),
        "myojin_program_bh_families": len(myojin_families),
        "myojin_program_testable_rows": sum(
            len(rows) for rows in myojin_families.values()
        ),
    }


def write_registry(
    candidate: Path,
    spec_hash: str,
    reports: Mapping[str, Mapping[str, object]],
) -> Path:
    require(
        tuple(reports) == tuple(REQUIRED_SCIENTIFIC_CHECKS),
        "scientific reports are not in the required exact order",
    )
    report_root = candidate / "logs/scientific_validation"
    prepared = []
    for check_id, payload in reports.items():
        require(
            payload.get("status") == "pass",
            f"refusing to register non-pass report: {check_id}",
        )
        path = report_root / f"{check_id}.json"
        prepared.append((check_id, path, payload))
    registry_rows = []
    for check_id, path, payload in prepared:
        materialize_immutable_json(
            path, payload, candidate, f"{check_id} scientific report"
        )
        registry_rows.append(
            {
                "check_id": check_id,
                "status": "pass",
                "report_path": path.relative_to(candidate).as_posix(),
                "report_sha256": sha256_file(path),
                "warnings_count": "0",
                "adjudication_path": "",
                "adjudication_sha256": "",
            }
        )
    registry = candidate / "manifests/scientific_validation_registry.tsv"
    materialize_immutable_tsv(
        registry,
        registry_rows,
        SCIENTIFIC_REGISTRY_FIELDS,
        candidate,
        "scientific validation registry",
    )
    return registry


def verify_registry(
    candidate: Path,
    reports: Mapping[str, Mapping[str, object]],
) -> Path:
    registry = candidate / "manifests/scientific_validation_registry.tsv"
    rows = read_tsv_exact(registry, SCIENTIFIC_REGISTRY_FIELDS)
    require(
        [row["check_id"] for row in rows] == list(REQUIRED_SCIENTIFIC_CHECKS),
        "scientific registry check order/membership drift",
    )
    for row in rows:
        check_id = row["check_id"]
        require(
            row["status"] == "pass" and row["warnings_count"] == "0",
            f"non-pass or warning-bearing registry row: {check_id}",
        )
        require(
            not row["adjudication_path"] and not row["adjudication_sha256"],
            f"unexpected zero-warning adjudication: {check_id}",
        )
        expected_path = f"logs/scientific_validation/{check_id}.json"
        require(
            row["report_path"] == expected_path,
            f"scientific report path drift: {check_id}",
        )
        path = candidate / expected_path
        require(
            path.is_file() and not path.is_symlink(),
            f"scientific report missing: {check_id}",
        )
        require(
            sha256_file(path) == row["report_sha256"],
            f"scientific report hash drift: {check_id}",
        )
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as error:
            raise ScientificValidationError(
                f"scientific report is invalid JSON: {check_id}"
            ) from error
        require(
            payload == reports[check_id],
            f"scientific report no longer rederives exactly: {check_id}",
        )
    return registry


def build_registry(
    project_root: Path,
    fixture_mode: bool = False,
    check_only: bool = False,
) -> dict[str, object]:
    base = validate_snapshot_base(project_root, fixture_mode=fixture_mode)
    project = base["project"]
    candidate = base["candidate"]
    require(
        isinstance(project, Path) and isinstance(candidate, Path),
        "candidate path type drift",
    )
    product_metrics = validate_rel02_04_products(
        project, candidate, str(base["spec_hash"]), fixture_mode
    )
    source_index = load_frozen_input_index(
        project, candidate, str(base["spec_hash"]), fixture_mode
    )
    source_rows = all_candidate_source_rows(candidate)
    producer_binding = scientific_producer_binding(project, candidate, source_index)

    checks: dict[str, Callable[[], dict[str, object]]] = {
        "hotspot_registry_consistency": lambda: check_hotspot(
            candidate, source_index, source_rows
        ),
        "spatial_assay_native_consistency": lambda: check_spatial(
            candidate, source_index, source_rows
        ),
        "numbers_claim_ledger": lambda: check_numbers_claims(
            candidate, product_metrics
        ),
        "figure_panel_manifest": lambda: check_panels(
            project, candidate, source_index, source_rows
        ),
        "passport_semantics": lambda: check_passports(
            candidate, source_index, source_rows
        ),
        "rebuild_determinism": lambda: check_determinism(
            project, candidate, source_index, fixture_mode
        ),
        "null_compatibility": lambda: check_null_compatibility(
            candidate, source_index, source_rows
        ),
        "manuscript_retired_claims": lambda: check_manuscript(
            project, candidate, source_rows
        ),
        "biological_unit_audit": lambda: check_biological_units(
            candidate, source_index, source_rows
        ),
        "multiplicity_universe_audit": lambda: check_multiplicity(
            candidate, source_index
        ),
    }
    require(
        tuple(checks) == tuple(REQUIRED_SCIENTIFIC_CHECKS),
        "scientific check implementation set drift",
    )

    # Run every predicate before the first report write. No synthetic PASS rows.
    reports = {
        check_id: report_payload(
            check_id, str(base["spec_hash"]), producer_binding, check()
        )
        for check_id, check in checks.items()
    }
    if check_only:
        registry = verify_registry(candidate, reports)
        status = "SCIENTIFIC_VALIDATION_REGISTRY_REDERIVED"
    else:
        registry = write_registry(candidate, str(base["spec_hash"]), reports)
        verify_registry(candidate, reports)
        status = "SCIENTIFIC_VALIDATION_REGISTRY_FROZEN"
    return {
        "status": status,
        "candidate_id": CANDIDATE_ID,
        "registry_path": registry.relative_to(project).as_posix(),
        "registry_sha256": sha256_file(registry),
        "n_checks": len(reports),
        "n_warnings": 0,
        "producer_manifest_sha256": producer_binding["producer_manifest_sha256"],
        "rel05_authorized_to_run": True,
        "canonical_promotion_authorized": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--fixture-mode", action="store_true")
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Rederive every predicate and verify an existing registry without writing.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        result = build_registry(args.project_root, args.fixture_mode, args.check_only)
    except (
        ScientificValidationError,
        CoordinatorContractError,
        ReleaseContractError,
    ) as error:
        print(f"PLAN60_SCIENTIFIC_VALIDATION_BLOCKED: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
