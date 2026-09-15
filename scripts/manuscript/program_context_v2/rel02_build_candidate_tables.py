#!/usr/bin/env python3
"""Build REL-02 ledgers and five-figure source tables from frozen snapshots only."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from rel01_snapshot_candidate import (
    materialize_immutable_json,
    materialize_immutable_tsv,
)
from rel05_validate_candidate import validate_snapshot_base
from release_common import (
    CANDIDATE_ID,
    ReleaseContractError,
    project_relative,
    sha256_file,
)
from release_products import (
    ACCEPTANCE_GATE_FIELDS,
    ANALYSIS_MANIFEST_FIELDS,
    CLAIM_FIELDS,
    EXPECTED_FIGURES,
    EXPECTED_FIGURE_TITLES,
    EXPECTED_RESULTS_SECTIONS,
    EXPECTED_SUPPLEMENTARY_SECTIONS,
    EXPECTED_SUPPLEMENTARY_FIGURES,
    EXPECTED_SUPPLEMENTARY_FIGURE_TITLES,
    EXPECTED_SUPPLEMENTARY_PANELS,
    JOURNAL_BRANCH,
    MYOJIN_ROLE,
    NUMBERS_FIELDS,
    PANEL_BLUEPRINT_FIELDS,
    SOURCE_DEPENDENCY_FIELDS,
    SOURCE_ROW_FIELDS,
    TRUE_KLEINER_TRANSITION_CENSUS,
    candidate_figure_root,
    load_frozen_input_index,
    parse_float,
    read_json_object,
    read_source_rows,
    read_table_flexible,
)


def validate_blueprint(
    candidate: Path,
    source_index: dict[str, dict[str, str]],
) -> tuple[dict[str, object], list[dict[str, str]], list[dict[str, str]]]:
    matches = [
        row
        for row in source_index.values()
        if row["artifact_role"] == "release_blueprint"
    ]
    if len(matches) != 1:
        raise ReleaseContractError(
            f"REL-02 requires exactly one frozen release_blueprint; observed {len(matches)}"
        )
    frozen = matches[0]
    blueprint_path = candidate / frozen["snapshot_path"]
    if sha256_file(blueprint_path) != frozen["sha256"]:
        raise ReleaseContractError("frozen release blueprint hash drift")
    blueprint = read_json_object(blueprint_path, "release blueprint")
    expected_scalars = {
        "contract_version": 1,
        "candidate_id": CANDIDATE_ID,
        "journal_branch": JOURNAL_BRANCH,
        "figure_count": 5,
        "myojin_role": MYOJIN_ROLE,
    }
    for key, expected in expected_scalars.items():
        if blueprint.get(key) != expected:
            raise ReleaseContractError(
                f"release blueprint drift for {key}: {blueprint.get(key)!r}"
            )
    figures = blueprint.get("figures")
    if not isinstance(figures, list):
        raise ReleaseContractError("release blueprint figures must be a list")
    observed_figures = []
    panel_rows = []
    panel_ids = set()
    for figure_order, figure in enumerate(figures, start=1):
        if not isinstance(figure, dict):
            raise ReleaseContractError("release blueprint figure must be an object")
        figure_id = str(figure.get("figure_id", ""))
        observed_figures.append(figure_id)
        if figure.get("title") != EXPECTED_FIGURE_TITLES.get(figure_id):
            raise ReleaseContractError(f"locked figure title drift: {figure_id}")
        panels = figure.get("panels")
        if not isinstance(panels, list) or not panels:
            raise ReleaseContractError(f"figure lacks panels: {figure_id}")
        for panel_order, panel in enumerate(panels, start=1):
            if not isinstance(panel, dict):
                raise ReleaseContractError(f"invalid panel object in {figure_id}")
            panel_id = str(panel.get("panel_id", ""))
            key = (figure_id, panel_id)
            if not panel_id or key in panel_ids:
                raise ReleaseContractError(f"blank or duplicate panel ID: {key}")
            panel_ids.add(key)
            panel_rows.append(
                {
                    "figure_id": figure_id,
                    "panel_id": panel_id,
                    "figure_order": str(figure_order),
                    "panel_order": str(panel_order),
                    "panel_title": str(panel.get("title", "")),
                    "panel_role": str(panel.get("role", "")),
                }
            )
    if tuple(observed_figures) != EXPECTED_FIGURES:
        raise ReleaseContractError(
            f"release blueprint must use locked five-figure order: {observed_figures}"
        )
    supplementary = blueprint.get("supplementary_figures")
    if not isinstance(supplementary, list):
        raise ReleaseContractError(
            "release blueprint supplementary_figures must be a list"
        )
    observed_supplementary = []
    supplementary_rows = []
    for figure_order, figure in enumerate(supplementary, start=1):
        if not isinstance(figure, dict):
            raise ReleaseContractError(
                "release blueprint supplementary figure must be an object"
            )
        figure_id = str(figure.get("figure_id", ""))
        observed_supplementary.append(figure_id)
        if figure.get("title") != EXPECTED_SUPPLEMENTARY_FIGURE_TITLES.get(figure_id):
            raise ReleaseContractError(
                f"locked supplementary figure title drift: {figure_id}"
            )
        panels = figure.get("panels")
        expected_panels = EXPECTED_SUPPLEMENTARY_PANELS.get(figure_id)
        if not isinstance(panels, list) or expected_panels is None:
            raise ReleaseContractError(
                f"unsupported supplementary figure or panel list: {figure_id}"
            )
        observed_panels = []
        for panel_order, panel in enumerate(panels, start=1):
            if not isinstance(panel, dict):
                raise ReleaseContractError(
                    f"invalid supplementary panel object in {figure_id}"
                )
            panel_id = str(panel.get("panel_id", ""))
            key = (figure_id, panel_id)
            if not panel_id or key in panel_ids:
                raise ReleaseContractError(
                    f"blank or duplicate supplementary panel ID: {key}"
                )
            panel_ids.add(key)
            observed_panels.append(
                (
                    panel_id,
                    str(panel.get("title", "")),
                    str(panel.get("role", "")),
                )
            )
            supplementary_rows.append(
                {
                    "figure_id": figure_id,
                    "panel_id": panel_id,
                    "figure_order": str(figure_order),
                    "panel_order": str(panel_order),
                    "panel_title": str(panel.get("title", "")),
                    "panel_role": str(panel.get("role", "")),
                }
            )
        if tuple(observed_panels) != expected_panels:
            raise ReleaseContractError(
                f"locked supplementary panel contract drift: {figure_id}"
            )
    if tuple(observed_supplementary) != EXPECTED_SUPPLEMENTARY_FIGURES:
        raise ReleaseContractError(
            "release blueprint must use the locked supplementary-figure order: "
            f"{observed_supplementary}"
        )
    sections = blueprint.get("results_sections")
    if sections != list(EXPECTED_RESULTS_SECTIONS):
        raise ReleaseContractError(
            "release blueprint Results sections differ from locked five-section branch"
        )
    sources = blueprint.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ReleaseContractError("release blueprint sources must be a nonempty list")
    seen_keys = set()
    for source in sources:
        if not isinstance(source, dict):
            raise ReleaseContractError("release blueprint source must be an object")
        source_key = str(source.get("source_key", ""))
        if not source_key or source_key in seen_keys:
            raise ReleaseContractError(
                f"blank or duplicate blueprint source: {source_key}"
            )
        seen_keys.add(source_key)
        if source_key not in source_index:
            raise ReleaseContractError(f"blueprint source is not frozen: {source_key}")
        if source.get("snapshot_path") != source_index[source_key]["snapshot_path"]:
            raise ReleaseContractError(
                f"blueprint source path differs from frozen manifest: {source_key}"
            )
        if source.get("kind") not in {
            "evidence_rows",
            "fibrosis_transition_raw",
        }:
            raise ReleaseContractError(
                f"unsupported blueprint source kind: {source.get('kind')}"
            )
    return blueprint, panel_rows, supplementary_rows


def row_template(**updates: object) -> dict[str, object]:
    row: dict[str, object] = {field: "" for field in SOURCE_ROW_FIELDS}
    row.update(updates)
    return row


def derive_fibrosis_rows(path: Path, source: dict[str, str]) -> list[dict[str, object]]:
    fields, raw = read_table_flexible(path)
    required = {"gene", "padj", "contrast", "n_cohorts", "n_samples"}
    if not required.issubset(fields):
        raise ReleaseContractError(
            f"fibrosis transition input lacks columns: {sorted(required - set(fields))}"
        )
    counts: Counter[str] = Counter()
    totals: Counter[str] = Counter()
    samples: dict[str, set[int]] = defaultdict(set)
    cohorts: dict[str, set[int]] = defaultdict(set)
    for row in raw:
        contrast = row["contrast"]
        if contrast not in TRUE_KLEINER_TRANSITION_CENSUS:
            continue
        totals[contrast] += 1
        padj = parse_float(row["padj"], f"{contrast} padj")
        if padj is not None and padj < 0.05:
            counts[contrast] += 1
        samples[contrast].add(int(row["n_samples"]))
        cohorts[contrast].add(int(row["n_cohorts"]))
    derived_census = {
        contrast: (
            next(iter(samples[contrast])) if len(samples[contrast]) == 1 else -1,
            next(iter(cohorts[contrast])) if len(cohorts[contrast]) == 1 else -1,
        )
        for contrast in TRUE_KLEINER_TRANSITION_CENSUS
    }
    if derived_census != TRUE_KLEINER_TRANSITION_CENSUS:
        raise ReleaseContractError(
            "true-Kleiner fibrosis transition census rederivation failed: "
            f"observed={derived_census}, expected={TRUE_KLEINER_TRANSITION_CENSUS}"
        )
    result = []
    for order, (contrast, (n_samples, n_cohorts)) in enumerate(
        TRUE_KLEINER_TRANSITION_CENSUS.items(), start=1
    ):
        deg_count = counts[contrast]
        label = contrast.replace("_vs_", " / ")
        claim_id = f"claim_fig2_{contrast}"
        shared = {
            "group_id": f"stage_{contrast}",
            "artifact_id": source["artifact_id"],
            "claim_id": claim_id,
            "section_id": "established_state_transcriptomics",
            "figure_id": "Figure2",
            "panel_id": "2A",
            "panel_title": "Cross-sectional stage-associated DEG counts",
            "panel_role": "transcriptomic_hero",
            "label": label,
            "category": "cross-sectional stage contrast",
            "biological_unit": "human biological sample",
            "model_contrast": contrast,
            "effect_unit": "count",
            "p_value": "",
            "q_value": "",
            "evidence_status": "stage_associated",
            "tested_universe": f"{totals[contrast]} tested genes",
            "source_dependence": "reused_source",
            "discovery_sources": "human_bulk_fibrosis_contrasts_v2",
            "evaluation_sources": "human_bulk_fibrosis_contrasts_v2",
            "reuse_detail": (
                "The displayed count is rederived from the same frozen cross-sectional "
                "contrast table used to define the stage-associated gene set."
            ),
            "independence_boundary": (
                "Descriptive rederivation only; this panel is not an independent replication."
            ),
            "allowed_wording": "cross-sectional stage-associated DEG count",
            "prohibited_wording": "longitudinal progression or transition rate",
            "claim_text": (
                f"The {contrast} cross-sectional contrast shows stage-associated "
                "transcriptomic remodeling."
            ),
            "next_experiment": "longitudinal paired-biopsy validation",
            "plot_order": str(order),
            "is_control": "false",
            "source_snapshot_path": source["snapshot_path"],
            "source_sha256": source["sha256"],
        }
        result.append(
            row_template(
                **shared,
                record_id=f"stage_{contrast}_deg_count",
                number_role="deg_count",
                plot_role="mark",
                value=str(deg_count),
                display_value=f"{deg_count:,}",
                numerator=str(deg_count),
                denominator=str(totals[contrast]),
                unit="genes at BH-adjusted p < 0.05",
                manuscript_included="true",
            )
        )
        result.append(
            row_template(
                **shared,
                record_id=f"stage_{contrast}_sample_count",
                number_role="sample_count",
                plot_role="annotation",
                value=str(n_samples),
                display_value=str(n_samples),
                numerator=str(n_samples),
                denominator="",
                unit="biological samples",
                manuscript_included="true",
            )
        )
        result.append(
            row_template(
                **shared,
                record_id=f"stage_{contrast}_cohort_count",
                number_role="cohort_count",
                plot_role="annotation",
                value=str(n_cohorts),
                display_value=str(n_cohorts),
                numerator=str(n_cohorts),
                denominator="",
                unit="cohorts",
                manuscript_included="true",
            )
        )
    return result


def load_all_rows(
    candidate: Path,
    blueprint: dict[str, object],
    source_index: dict[str, dict[str, str]],
) -> list[dict[str, object]]:
    all_rows: list[dict[str, object]] = []
    for declared in blueprint["sources"]:
        source = source_index[declared["source_key"]]
        path = candidate / source["snapshot_path"]
        if sha256_file(path) != source["sha256"]:
            raise ReleaseContractError(
                f"frozen source changed before REL-02: {declared['source_key']}"
            )
        kind = declared["kind"]
        if kind == "evidence_rows":
            rows = read_source_rows(path)
            for row in rows:
                if row["artifact_id"] != source["artifact_id"]:
                    raise ReleaseContractError(
                        f"source artifact identity drift in {declared['source_key']}"
                    )
                all_rows.append(
                    {
                        **row,
                        "source_snapshot_path": source["snapshot_path"],
                        "source_sha256": source["sha256"],
                    }
                )
        elif kind == "fibrosis_transition_raw":
            all_rows.extend(derive_fibrosis_rows(path, source))
    record_ids = [str(row["record_id"]) for row in all_rows]
    if len(record_ids) != len(set(record_ids)):
        duplicates = [item for item, count in Counter(record_ids).items() if count > 1]
        raise ReleaseContractError(f"duplicate REL-02 record IDs: {duplicates}")
    return sorted(
        all_rows,
        key=lambda row: (
            str(row["figure_id"]),
            str(row["panel_id"]),
            int(str(row["plot_order"])),
            str(row["plot_role"]),
            str(row["record_id"]),
        ),
    )


def build_ledgers(
    rows: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    numbers = []
    by_claim: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        number_id = f"N::{row['record_id']}"
        numbers.append(
            {
                "number_id": number_id,
                "claim_id": row["claim_id"],
                "number_role": row["number_role"],
                "raw_value": row["value"],
                "display_value": row["display_value"],
                "numerator": row["numerator"],
                "denominator": row["denominator"],
                "unit_of_replication": row["biological_unit"],
                "model_contrast": row["model_contrast"],
                "source_artifact": row["source_snapshot_path"],
                "source_row_filter": f"record_id={row['record_id']}",
                "release_id": CANDIDATE_ID,
                "allowed_wording": row["allowed_wording"],
                "prohibited_wording": row["prohibited_wording"],
                "figure_id": row["figure_id"],
                "panel_id": row["panel_id"],
                "manuscript_included": row["manuscript_included"],
            }
        )
        by_claim[str(row["claim_id"])].append(row)
    claims = []
    dependencies = []
    for claim_id in sorted(by_claim):
        claim_rows = by_claim[claim_id]
        invariant_fields = (
            "section_id",
            "claim_text",
            "allowed_wording",
            "prohibited_wording",
            "source_dependence",
            "discovery_sources",
            "evaluation_sources",
            "reuse_detail",
            "independence_boundary",
            "tested_universe",
            "biological_unit",
        )
        invariants = {
            field: {str(row[field]) for row in claim_rows} for field in invariant_fields
        }
        drift = {
            field: values for field, values in invariants.items() if len(values) != 1
        }
        if drift:
            raise ReleaseContractError(f"claim invariant drift for {claim_id}: {drift}")
        number_ids = sorted(f"N::{row['record_id']}" for row in claim_rows)
        artifact_ids = sorted({str(row["artifact_id"]) for row in claim_rows})
        statuses = sorted({str(row["evidence_status"]) for row in claim_rows})
        claims.append(
            {
                "claim_id": claim_id,
                "section_id": next(iter(invariants["section_id"])),
                "claim_text": next(iter(invariants["claim_text"])),
                "claim_status": statuses[0] if len(statuses) == 1 else "mixed",
                "number_ids": ";".join(number_ids),
                "source_artifact_ids": ";".join(artifact_ids),
                "allowed_wording": next(iter(invariants["allowed_wording"])),
                "prohibited_wording": next(iter(invariants["prohibited_wording"])),
                "source_dependence": next(iter(invariants["source_dependence"])),
                "tested_universe": next(iter(invariants["tested_universe"])),
                "biological_unit": next(iter(invariants["biological_unit"])),
            }
        )
        dependencies.append(
            {
                "claim_id": claim_id,
                "source_dependence_class": next(iter(invariants["source_dependence"])),
                "discovery_source": next(iter(invariants["discovery_sources"])),
                "evaluation_source": next(iter(invariants["evaluation_sources"])),
                "reuse_detail": next(iter(invariants["reuse_detail"])),
                "independence_boundary": next(
                    iter(invariants["independence_boundary"])
                ),
            }
        )
    return numbers, claims, dependencies


def build_candidate_tables(
    project_root: Path, fixture_mode: bool = False
) -> dict[str, object]:
    base = validate_snapshot_base(project_root, fixture_mode=fixture_mode)
    project = base["project"]
    candidate = base["candidate"]
    closure = base["closure"]
    if not isinstance(project, Path) or not isinstance(candidate, Path):
        raise ReleaseContractError("internal REL-02 path type drift")
    source_index = load_frozen_input_index(
        project, candidate, str(base["spec_hash"]), fixture_mode
    )
    blueprint, declared_panels, declared_supplementary_panels = validate_blueprint(
        candidate, source_index
    )
    plan40 = next(
        workstream
        for workstream in closure.workstreams
        if workstream.row["workstream_id"] == "PLAN40"
    )
    if (
        plan40.row["include_main"] != "false"
        or plan40.row["include_supplement"] != "true"
        or "supplement" not in plan40.row["figure_role"].lower()
    ):
        raise ReleaseContractError(
            "locked five-figure branch requires Myojin supplement-only verdict"
        )
    rows = load_all_rows(candidate, blueprint, source_index)
    observed_main_panels = {
        (str(row["figure_id"]), str(row["panel_id"]))
        for row in rows
        if str(row["figure_id"]) in EXPECTED_FIGURES
    }
    expected_main_panels = {
        (row["figure_id"], row["panel_id"]) for row in declared_panels
    }
    if observed_main_panels != expected_main_panels:
        raise ReleaseContractError(
            "main panel source coverage differs from locked blueprint: "
            f"missing={sorted(expected_main_panels - observed_main_panels)}, "
            f"extra={sorted(observed_main_panels - expected_main_panels)}"
        )
    observed_supplementary_panels = {
        (str(row["figure_id"]), str(row["panel_id"]))
        for row in rows
        if str(row["figure_id"]).startswith("FigureS")
    }
    expected_supplementary_panels = {
        (row["figure_id"], row["panel_id"]) for row in declared_supplementary_panels
    }
    if observed_supplementary_panels != expected_supplementary_panels:
        raise ReleaseContractError(
            "supplementary source coverage differs from locked blueprint: "
            f"missing={sorted(expected_supplementary_panels - observed_supplementary_panels)}, "
            f"extra={sorted(observed_supplementary_panels - expected_supplementary_panels)}"
        )
    allowed_figure_ids = set(EXPECTED_FIGURES) | set(EXPECTED_SUPPLEMENTARY_FIGURES)
    unexpected_figure_ids = {
        str(row["figure_id"])
        for row in rows
        if str(row["figure_id"]) not in allowed_figure_ids
    }
    if unexpected_figure_ids:
        raise ReleaseContractError(
            f"source rows contain figures outside the locked branch: {sorted(unexpected_figure_ids)}"
        )
    if any(str(row["figure_id"]) == "Figure6" for row in rows):
        raise ReleaseContractError("five-figure branch cannot contain Figure6")
    myojin_rows = [row for row in rows if row["section_id"] == "supplementary_myojin"]
    if not myojin_rows or any(
        str(row["figure_id"]) != "FigureS1" for row in myojin_rows
    ):
        raise ReleaseContractError(
            "complete Myojin fixture must remain in a supplementary figure"
        )
    nmf_rows = [row for row in rows if row["section_id"] == "supplementary_nmf"]
    if not nmf_rows or any(str(row["figure_id"]) != "FigureS2" for row in nmf_rows):
        raise ReleaseContractError(
            "continuous NMF axes and stability evidence must remain in FigureS2"
        )
    observed_supplementary_sections = {
        str(row["section_id"])
        for row in rows
        if str(row["section_id"]).startswith("supplementary_")
    }
    if observed_supplementary_sections != set(EXPECTED_SUPPLEMENTARY_SECTIONS):
        raise ReleaseContractError(
            "supplementary source sections differ from the locked branch: "
            f"{sorted(observed_supplementary_sections)}"
        )
    yak_module8 = [row for row in rows if str(row["group_id"]) == "yak_module8"]
    if yak_module8:
        observed_roles = {str(row["number_role"]) for row in yak_module8}
        required_roles = {
            "descriptive_median_slope_direction",
            "inferential_signed_stouffer_direction",
        }
        if not required_roles.issubset(observed_roles):
            raise ReleaseContractError(
                "Yakubovsky module 8 must keep descriptive median-slope and "
                "inferential signed-Stouffer directions separate"
            )
    # Complete every source/claim invariant check before the first REL-02
    # candidate product is materialized. A bad frozen contract must not leave a
    # partially extended read-only candidate.
    numbers, claims, dependencies = build_ledgers(rows)

    panel_table_paths: dict[tuple[str, str], Path] = {}
    for key in sorted({(str(row["figure_id"]), str(row["panel_id"])) for row in rows}):
        figure_id, panel_id = key
        panel_rows = [
            row
            for row in rows
            if str(row["figure_id"]) == figure_id and str(row["panel_id"]) == panel_id
        ]
        relative = Path("figure_sources") / figure_id.lower() / f"panel_{panel_id}.tsv"
        path = candidate / relative
        materialize_immutable_tsv(
            path,
            panel_rows,
            SOURCE_ROW_FIELDS,
            candidate,
            f"{figure_id} {panel_id} source table",
        )
        panel_table_paths[key] = path

    numbers_path = candidate / "tables/numbers_ledger.tsv"
    claims_path = candidate / "claims/claim_ledger.tsv"
    dependencies_path = candidate / "claims/source_dependency_ledger.tsv"
    materialize_immutable_tsv(
        numbers_path, numbers, NUMBERS_FIELDS, candidate, "numbers ledger"
    )
    materialize_immutable_tsv(
        claims_path, claims, CLAIM_FIELDS, candidate, "claim ledger"
    )
    materialize_immutable_tsv(
        dependencies_path,
        dependencies,
        SOURCE_DEPENDENCY_FIELDS,
        candidate,
        "source dependency ledger",
    )
    gate_rows = [
        {
            "gate_id": workstream.row["gate_id"],
            "workstream_id": workstream.row["workstream_id"],
            "terminal_state": workstream.row["terminal_state"],
            "gate_verdict": workstream.row["gate_verdict"],
            "include_main": workstream.row["include_main"],
            "include_supplement": workstream.row["include_supplement"],
            "figure_role": workstream.row["figure_role"],
            "closure_sha256": closure.closure_sha256,
        }
        for workstream in closure.workstreams
    ]
    gates_path = candidate / "manifests/acceptance_gates.tsv"
    materialize_immutable_tsv(
        gates_path,
        gate_rows,
        ACCEPTANCE_GATE_FIELDS,
        candidate,
        "acceptance gate ledger",
    )

    panel_blueprint_rows = []
    for declared in declared_panels:
        key = (declared["figure_id"], declared["panel_id"])
        source_path = panel_table_paths[key]
        pdf_path = (
            candidate_figure_root(project)
            / declared["figure_id"].lower()
            / f"panel_{declared['panel_id']}.pdf"
        )
        panel_blueprint_rows.append(
            {
                **declared,
                "source_table_path": source_path.relative_to(candidate).as_posix(),
                "source_table_sha256": sha256_file(source_path),
                "pdf_path": project_relative(project, pdf_path),
                "normal_font_pt": "6",
                "format": "PDF",
            }
        )
    for declared in declared_supplementary_panels:
        figure_id = declared["figure_id"]
        panel_id = declared["panel_id"]
        key = (figure_id, panel_id)
        source_path = panel_table_paths[key]
        pdf_path = (
            candidate_figure_root(project)
            / "supplementary"
            / figure_id.lower()
            / f"panel_{panel_id}.pdf"
        )
        panel_blueprint_rows.append(
            {
                "figure_id": figure_id,
                "panel_id": panel_id,
                "figure_order": "99",
                "panel_order": declared["panel_order"],
                "panel_title": declared["panel_title"],
                "panel_role": declared["panel_role"],
                "source_table_path": source_path.relative_to(candidate).as_posix(),
                "source_table_sha256": sha256_file(source_path),
                "pdf_path": project_relative(project, pdf_path),
                "normal_font_pt": "6",
                "format": "PDF",
            }
        )
    panel_blueprint_path = candidate / "manifests/figure_panel_blueprint.tsv"
    materialize_immutable_tsv(
        panel_blueprint_path,
        panel_blueprint_rows,
        PANEL_BLUEPRINT_FIELDS,
        candidate,
        "figure panel blueprint",
    )

    manifested = []
    for key, source in sorted(source_index.items()):
        path = candidate / source["snapshot_path"]
        manifested.append(
            {
                "artifact_id": key,
                "artifact_role": source["artifact_role"],
                "candidate_path": source["snapshot_path"],
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "producer_id": "REL01_or_BASE",
                "input_snapshot_spec_sha256": base["spec_hash"],
            }
        )
    for artifact_id, role, path in (
        ("numbers_ledger", "numbers_ledger", numbers_path),
        ("claim_ledger", "claim_ledger", claims_path),
        ("source_dependency_ledger", "source_dependency_ledger", dependencies_path),
        ("acceptance_gates", "acceptance_gates", gates_path),
        ("figure_panel_blueprint", "figure_panel_blueprint", panel_blueprint_path),
    ):
        manifested.append(
            {
                "artifact_id": artifact_id,
                "artifact_role": role,
                "candidate_path": path.relative_to(candidate).as_posix(),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "producer_id": "rel02_build_candidate_tables",
                "input_snapshot_spec_sha256": base["spec_hash"],
            }
        )
    for (figure_id, panel_id), path in sorted(panel_table_paths.items()):
        manifested.append(
            {
                "artifact_id": f"source_{figure_id}_{panel_id}",
                "artifact_role": "figure_source_table",
                "candidate_path": path.relative_to(candidate).as_posix(),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "producer_id": "rel02_build_candidate_tables",
                "input_snapshot_spec_sha256": base["spec_hash"],
            }
        )
    analysis_path = candidate / "manifests/analysis_release_manifest.tsv"
    materialize_immutable_tsv(
        analysis_path,
        manifested,
        ANALYSIS_MANIFEST_FIELDS,
        candidate,
        "analysis release manifest",
    )
    rel02_state = {
        "contract_version": 1,
        "candidate_id": CANDIDATE_ID,
        "state": "REL02_TABLES_FROZEN",
        "input_snapshot_spec_sha256": base["spec_hash"],
        "base_input_transition_sha256": sha256_file(
            candidate / "manifests/base_input_transition.json"
        ),
        "journal_branch": JOURNAL_BRANCH,
        "figure_count": 5,
        "myojin_role": MYOJIN_ROLE,
        "numbers_ledger_sha256": sha256_file(numbers_path),
        "claim_ledger_sha256": sha256_file(claims_path),
        "source_dependency_ledger_sha256": sha256_file(dependencies_path),
        "acceptance_gates_sha256": sha256_file(gates_path),
        "figure_panel_blueprint_sha256": sha256_file(panel_blueprint_path),
        "analysis_release_manifest_sha256": sha256_file(analysis_path),
        "n_numbers": len(numbers),
        "n_claims": len(claims),
        "n_main_panels": len(declared_panels),
        "n_supplementary_panels": len(declared_supplementary_panels),
        "fixture_mode": fixture_mode,
        "canonical_promotion_status": "not_promoted",
    }
    state_path = candidate / "manifests/rel02_state.json"
    materialize_immutable_json(state_path, rel02_state, candidate, "REL-02 state")
    return {
        "status": "REL02_TABLES_FROZEN",
        "candidate_id": CANDIDATE_ID,
        "journal_branch": JOURNAL_BRANCH,
        "figure_count": 5,
        "myojin_role": MYOJIN_ROLE,
        "n_numbers": len(numbers),
        "n_claims": len(claims),
        "n_main_panels": len(declared_panels),
        "n_supplementary_panels": len(declared_supplementary_panels),
        "rel02_state_sha256": sha256_file(state_path),
        "canonical_promotion_status": "not_promoted",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--fixture-mode", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        result = build_candidate_tables(args.project_root, args.fixture_mode)
    except ReleaseContractError as error:
        raise SystemExit(f"REL02_BLOCKED: {error}") from error
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
