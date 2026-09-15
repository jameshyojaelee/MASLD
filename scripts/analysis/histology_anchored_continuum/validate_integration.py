#!/usr/bin/env python3
"""Read-only integration checks for the continuum manuscript candidate."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
FIG4 = ROOT / "figures/main/fig4_singlecell_programs"
CORE = ROOT / (
    "RNA-seq/results/histology_anchored_continuum/candidates/"
    "hac-continuum-20260818T024923Z"
)
MOLECULAR = ROOT / (
    "RNA-seq/results/histology_anchored_continuum/molecular_layers/"
    "hac-molecular-layers-20260818T173348Z"
)
COMPARISON = ROOT / (
    "figures/candidates/"
    "histology-vs-continuum-comparison-20260818T192630Z"
)
FIVE_COHORT_PROMOTION = (
    FIG4 / "manifests/continuum_five_cohort_revision_20260824.tsv"
)
FIG4F_MINIMAL_PROMOTION = (
    FIG4 / "manifests/fig4f_minimal_promotion_20260825.tsv"
)
FIG4F_MINIMAL_SOURCE = (
    FIG4 / "source_tables/current_candidate/fig4f_minimal_20260825"
)
FIG4H_SYSTEM_PROMOTION = (
    FIG4 / "manifests/fig4h_system_stage_continuum_promotion_20260825.tsv"
)
FIG4H_SYSTEM_SOURCE = (
    FIG4 / "source_tables/current_candidate/fig4h_system_stage_continuum_20260825"
)
ATLAS_PROMOTION = (
    FIG4 / "manifests/molecular_systems_atlas_promotion_20260824.tsv"
)
ATLAS_SOURCE = (
    FIG4 / "source_tables/current_candidate/molecular_systems_atlas_20260824"
)
POOLED_SYSTEMS_PROMOTION = (
    FIG4 / "manifests/molecular_systems_pooled_promotion_20260824_v2.tsv"
)
POOLED_SYSTEMS_SOURCE = (
    FIG4 / "source_tables/current_candidate/molecular_systems_pooled_20260824_v2"
)
SYSTEMS_VISUAL_PROMOTION = FIG4 / (
    "manifests/molecular_systems_visual_revision_promotion_20260824.tsv"
)
SYSTEMS_VISUAL_SOURCE = FIG4 / (
    "source_tables/current_candidate/molecular_systems_visual_revision_20260824"
)
POOLED_NARROW_PROMOTION = FIG4 / (
    "manifests/molecular_systems_pooled_narrow_promotion_20260824.tsv"
)
POOLED_NARROW_SOURCE = FIG4 / (
    "source_tables/current_candidate/molecular_systems_pooled_narrow_20260824"
)
SYSTEMS_INFERENCE_PROMOTION = FIG4 / (
    "manifests/molecular_systems_inference_promotion_20260824.tsv"
)
SYSTEMS_INFERENCE_SOURCE = FIG4 / (
    "source_tables/current_candidate/molecular_systems_inference_20260824"
)
SYSTEMS_ASTERISK_PROMOTION = FIG4 / (
    "manifests/molecular_systems_inference_asterisk_promotion_20260825.tsv"
)
SYSTEMS_ASTERISK_SOURCE = FIG4 / (
    "source_tables/current_candidate/molecular_systems_inference_asterisk_20260825"
)
SYSTEMS_NO_DELTA_PROMOTION = FIG4 / (
    "manifests/molecular_systems_no_delta_promotion_20260825.tsv"
)
SYSTEMS_NO_DELTA_SOURCE = FIG4 / (
    "source_tables/current_candidate/molecular_systems_no_delta_20260825"
)
SYSTEMS_CAPTION_FREE_PROMOTION = FIG4 / (
    "manifests/molecular_systems_caption_free_promotion_20260825.tsv"
)
SYSTEMS_CAPTION_FREE_SOURCE = FIG4 / (
    "source_tables/current_candidate/molecular_systems_caption_free_20260825"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_tsv(path: Path) -> list[dict[str, str]]:
    require(path.is_file() and path.stat().st_size > 0, f"Missing TSV: {path}")
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    core_summary = json.loads(
        (CORE / "validation/validation_summary.json").read_text()
    )
    require(core_summary.get("status") == "PASS", "Core candidate did not pass")

    molecular_summary = json.loads(
        (MOLECULAR / "validation/validation_summary.json").read_text()
    )
    require(
        molecular_summary.get("validation_pass") is True
        and molecular_summary.get("n_failed") == 0,
        "Molecular-layer candidate did not pass",
    )

    index_rows = read_tsv(FIG4 / "CANDIDATE_PANEL_INDEX.tsv")
    row4f = [row for row in index_rows if row["callout"] == "4F"]
    require(len(row4f) == 1, "Figure 4 index must contain exactly one 4F row")
    indexed_4f = (FIG4 / row4f[0]["candidate_relative_path"]).resolve()
    five_cohort_rows = read_tsv(FIVE_COHORT_PROMOTION)
    promoted_4f = [row for row in five_cohort_rows if row["callout"] == "4F"]
    require(len(promoted_4f) == 1, "Five-cohort manifest lacks exactly one 4F row")
    require(
        indexed_4f.is_file()
        and indexed_4f.stat().st_size > 10_000
        and indexed_4f.read_bytes()[:4] == b"%PDF",
        "Figure 4F is missing or not a nonempty PDF",
    )
    require(
        sha256((ROOT / promoted_4f[0]["source"]).resolve())
        == promoted_4f[0]["new_sha256"],
        "Historical five-cohort Figure 4F source hash drift",
    )
    minimal_rows = read_tsv(FIG4F_MINIMAL_PROMOTION)
    require(
        len(minimal_rows) == 1 and minimal_rows[0]["callout"] == "4F",
        "Minimal Figure 4F promotion manifest drifted",
    )
    minimal_row = minimal_rows[0]
    require(
        minimal_row["old_sha256"] == promoted_4f[0]["new_sha256"],
        "Minimal Figure 4F does not supersede five-cohort context",
    )
    require(
        sha256(indexed_4f) == minimal_row["new_sha256"],
        "Figure 4F differs from its promoted minimal hash",
    )
    minimal_candidate = (ROOT / minimal_row["source"]).resolve()
    require(
        minimal_candidate.is_file() and sha256(minimal_candidate) == minimal_row["new_sha256"],
        "Minimal Figure 4F candidate hash drift",
    )
    require(
        minimal_row["trajectory_cohorts_displayed"] == "5"
        and minimal_row["forest_panel_displayed"] == "FALSE"
        and minimal_row["embedded_explanatory_text"] == "FALSE"
        and minimal_row["participant_size_legend_displayed"] == "FALSE"
        and minimal_row["source_overlap_wording_in_art"] == "FALSE"
        and minimal_row["window_values_changed"] == "FALSE"
        and minimal_row["statistical_model_changed"] == "FALSE",
        "Minimal Figure 4F display contract drifted",
    )
    for row in read_tsv(FIG4F_MINIMAL_SOURCE / "promoted_source_checksums.tsv"):
        source = FIG4F_MINIMAL_SOURCE / row["relative_path"]
        require(source.is_file() and sha256(source) == row["sha256"],
                f"Promoted minimal Figure 4F source drift: {source}")

    row4h = [row for row in index_rows if row["callout"] == "4H"]
    require(len(row4h) == 1, "Figure 4 index must contain exactly one 4H row")
    indexed_4h = (FIG4 / row4h[0]["candidate_relative_path"]).resolve()
    fig4h_rows = read_tsv(FIG4H_SYSTEM_PROMOTION)
    require(
        len(fig4h_rows) == 1 and fig4h_rows[0]["callout"] == "4H",
        "Figure 4H promotion manifest drifted",
    )
    fig4h_row = fig4h_rows[0]
    fig4h_candidate = (ROOT / fig4h_row["source"]).resolve()
    require(
        indexed_4h.is_file()
        and indexed_4h.stat().st_size > 10_000
        and indexed_4h.read_bytes()[:4] == b"%PDF"
        and sha256(indexed_4h) == fig4h_row["new_sha256"]
        and fig4h_candidate.is_file()
        and sha256(fig4h_candidate) == fig4h_row["new_sha256"],
        "Figure 4H active or candidate hash drift",
    )
    require(
        fig4h_row["complete_system_family"] == "43/43"
        and fig4h_row["strict_maxt_supported"] == "18/43"
        and fig4h_row["same_direction"] == "36/43"
        and abs(float(fig4h_row["spearman_rho"]) - 0.658) < 0.001
        and fig4h_row["fixed_outcome_blind_labels"] == "15/43"
        and fig4h_row["window_values_used"] == "FALSE"
        and fig4h_row["system_geometry_used"] == "FALSE"
        and fig4h_row["statistical_model_changed"] == "FALSE",
        "Figure 4H display or inference contract drifted",
    )
    fig4h_source = read_tsv(
        FIG4H_SYSTEM_SOURCE / "source_tables/fig4h_system_stage_continuum.tsv"
    )
    require(
        len(fig4h_source) == 43
        and sum(row["strict_maxt_supported"] == "TRUE" for row in fig4h_source) == 18
        and sum(row["fixed_outcome_blind_label"] == "TRUE" for row in fig4h_source) == 15,
        "Figure 4H complete source family drifted",
    )
    fig4h_audit = read_tsv(FIG4H_SYSTEM_SOURCE / "provenance/audit.tsv")
    require(
        fig4h_audit and all(row["passed"] == "TRUE" for row in fig4h_audit),
        "Figure 4H source audit failed",
    )
    for row in read_tsv(FIG4H_SYSTEM_SOURCE / "promoted_source_checksums.tsv"):
        source = FIG4H_SYSTEM_SOURCE / row["relative_path"]
        require(source.is_file() and sha256(source) == row["sha256"],
                f"Promoted Figure 4H source drift: {source}")

    atlas_rows = read_tsv(ATLAS_PROMOTION)
    require(
        {row["callout"] for row in atlas_rows} == {"4G", "S4Q", "S4R", "S4S", "S4T"},
        "Molecular-systems promotion callouts drifted",
    )
    indexed_atlas_rows = {
        row["callout"]: row for row in index_rows
        if row["callout"] in {"4G", "S4Q", "S4R", "S4S", "S4T"}
    }
    require(len(indexed_atlas_rows) == 5, "Figure 4 index lacks an atlas callout")

    for audit_name in ("outcome_blind_audit.tsv", "overlay_audit.tsv"):
        audit_rows = read_tsv(ATLAS_SOURCE / "provenance" / audit_name)
        require(
            audit_rows and all(row["passed"] == "TRUE" for row in audit_rows),
            f"Molecular-systems audit failed: {audit_name}",
        )
    atlas_summary = read_tsv(ATLAS_SOURCE / "overlays/atlas_summary.tsv")
    atlas_metrics = {row["metric"]: row["value"] for row in atlas_summary}
    require(atlas_metrics.get("membership_defined_systems") == "43",
            "Molecular-system count drifted")
    require(atlas_metrics.get("excluded_microcomponent_nodes") == "6",
            "Microcomponent-node count drifted")
    for row in read_tsv(ATLAS_SOURCE / "promoted_source_checksums.tsv"):
        path = ROOT / row["path"]
        require(path.is_file() and sha256(path) == row["sha256"],
                f"Promoted molecular-system source drift: {path}")

    pooled_rows = read_tsv(POOLED_SYSTEMS_PROMOTION)
    require(
        len(pooled_rows) == 1 and pooled_rows[0]["callout"] == "S4U",
        "Pooled molecular-systems promotion manifest drifted",
    )
    indexed_s4u = [row for row in index_rows if row["callout"] == "S4U"]
    require(len(indexed_s4u) == 1, "Figure 4 index lacks exactly one S4U row")
    pooled_panel = (FIG4 / indexed_s4u[0]["candidate_relative_path"]).resolve()
    manifested_panel = (ROOT / pooled_rows[0]["destination"]).resolve()
    require(pooled_panel == manifested_panel, "S4U index and manifest disagree")
    require(
        pooled_panel.is_file()
        and pooled_panel.stat().st_size > 10_000
        and pooled_panel.read_bytes()[:4] == b"%PDF",
        "S4U is missing or not a nonempty PDF",
    )
    pooled_audit = read_tsv(POOLED_SYSTEMS_SOURCE / "provenance/audit.tsv")
    require(
        pooled_audit and all(row["passed"] == "TRUE" for row in pooled_audit),
        "Pooled molecular-systems audit failed",
    )
    pooled_source = read_tsv(
        POOLED_SYSTEMS_SOURCE / "source_tables/pooled_molecular_system_windows.tsv"
    )
    require(len(pooled_source) == 43 * 9, "Pooled system-window family drifted")
    require(
        {row["n_cohorts"] for row in pooled_source} == {"5"},
        "Pooled system windows do not retain five cohorts",
    )
    for row in read_tsv(POOLED_SYSTEMS_SOURCE / "promoted_source_checksums.tsv"):
        path = ROOT / row["path"]
        require(path.is_file() and sha256(path) == row["sha256"],
                f"Promoted pooled-system source drift: {path}")

    visual_rows = read_tsv(SYSTEMS_VISUAL_PROMOTION)
    require(
        {row["callout"] for row in visual_rows}
        == {"4G", "S4Q", "S4R", "S4S", "S4T", "S4U"},
        "Molecular-systems visual-revision callouts drifted",
    )
    indexed_visual_rows = {
        row["callout"]: row for row in index_rows
        if row["callout"] in {"4G", "S4Q", "S4R", "S4S", "S4T", "S4U"}
    }
    require(len(indexed_visual_rows) == 6,
            "Figure 4 index lacks a visual-revision callout")
    for row in visual_rows:
        callout = row["callout"]
        indexed = (
            FIG4 / indexed_visual_rows[callout]["candidate_relative_path"]
        ).resolve()
        manifested = (ROOT / row["destination"]).resolve()
        require(indexed == manifested,
                f"{callout} visual-revision index and manifest disagree")
        require(
            indexed.is_file()
            and indexed.stat().st_size > 10_000
            and indexed.read_bytes()[:4] == b"%PDF",
            f"{callout} visual revision is missing or not a nonempty PDF",
        )
        visual_candidate = (ROOT / row["source"]).resolve()
        require(
            visual_candidate.is_file()
            and sha256(visual_candidate) == row["new_sha256"],
            f"Historical {callout} visual-revision hash drift",
        )
        require(row["data_values_changed"] == "FALSE",
                f"{callout} visual revision claims changed data values")
        require(row["geometry_changed"] == "FALSE",
                f"{callout} visual revision claims changed geometry")
    visual_audit = read_tsv(SYSTEMS_VISUAL_SOURCE / "provenance/audit.tsv")
    require(
        visual_audit and all(row["passed"] == "TRUE" for row in visual_audit),
        "Molecular-systems visual-revision audit failed",
    )
    direction_order = read_tsv(
        SYSTEMS_VISUAL_SOURCE
        / "source_tables/pooled_system_late_minus_early_order.tsv"
    )
    require(len(direction_order) == 43,
            "Pooled molecular-system display-order family drifted")
    require(
        [int(row["display_rank"]) for row in direction_order] == list(range(1, 44)),
        "Pooled molecular-system display ranks drifted",
    )
    display_values = read_tsv(
        SYSTEMS_VISUAL_SOURCE
        / "source_tables/pooled_system_heatmap_display_values.tsv"
    )
    require(len(display_values) == 43 * 10,
            "Pooled heatmap display-value family drifted")
    for row in read_tsv(SYSTEMS_VISUAL_SOURCE / "promoted_source_checksums.tsv"):
        path = ROOT / row["path"]
        require(path.is_file() and sha256(path) == row["sha256"],
                f"Promoted molecular-system visual source drift: {path}")

    narrow_rows = read_tsv(POOLED_NARROW_PROMOTION)
    require(
        len(narrow_rows) == 1 and narrow_rows[0]["callout"] == "S4U",
        "Narrow pooled-system promotion manifest drifted",
    )
    narrow_row = narrow_rows[0]
    narrow_panel = (ROOT / narrow_row["destination"]).resolve()
    require(narrow_panel == pooled_panel,
            "Narrow S4U manifest and Figure 4 index disagree")
    require(narrow_row["data_values_changed"] == "FALSE",
            "Narrow S4U claims changed data values")
    require(narrow_row["geometry_changed"] == "FALSE",
            "Narrow S4U claims changed system geometry")
    narrow_audit = read_tsv(POOLED_NARROW_SOURCE / "provenance/audit.tsv")
    require(
        narrow_audit and all(row["passed"] == "TRUE" for row in narrow_audit),
        "Narrow S4U source audit failed",
    )
    for row in read_tsv(POOLED_NARROW_SOURCE / "promoted_source_checksums.tsv"):
        path = ROOT / row["path"]
        require(path.is_file() and sha256(path) == row["sha256"],
                f"Promoted narrow S4U source drift: {path}")

    inference_rows = read_tsv(SYSTEMS_INFERENCE_PROMOTION)
    require(
        len(inference_rows) == 1 and inference_rows[0]["callout"] == "S4U",
        "System-inference promotion manifest drifted",
    )
    inference_row = inference_rows[0]
    inference_panel = (ROOT / inference_row["destination"]).resolve()
    require(inference_panel == pooled_panel,
            "System-inference S4U manifest and Figure 4 index disagree")
    require(inference_row["old_sha256"] == narrow_row["new_sha256"],
            "System-inference S4U does not supersede the narrow revision")
    inference_candidate_panel = (ROOT / inference_row["source"]).resolve()
    require(
        inference_candidate_panel.is_file()
        and sha256(inference_candidate_panel) == inference_row["new_sha256"],
        "Historical system-inference S4U candidate hash drift",
    )
    require(inference_row["primary_bh_supported"] == "33/43",
            "System-inference primary BH census drifted")
    require(inference_row["strict_maxT_supported"] == "18/43",
            "System-inference maxT census drifted")
    require(inference_row["window_values_changed"] == "FALSE",
            "System-inference S4U claims changed window values")
    require(inference_row["display_order_changed"] == "FALSE",
            "System-inference S4U claims changed display order")
    require(inference_row["windows_inferentially_tested"] == "FALSE",
            "System-inference S4U assigns inference to window means")
    inference_audit = read_tsv(SYSTEMS_INFERENCE_SOURCE / "provenance/audit.tsv")
    require(
        inference_audit and all(row["passed"] == "TRUE" for row in inference_audit),
        "System-inference S4U source audit failed",
    )
    inference_strip = read_tsv(
        SYSTEMS_INFERENCE_SOURCE / "source_tables/system_inference_strip.tsv"
    )
    require(len(inference_strip) == 43 * 2,
            "System-inference strip family drifted")
    require(
        sum(
            row["supported"] == "TRUE"
            for row in inference_strip if row["inference_column"] == "BH"
        ) == 33,
        "System-inference strip BH support count drifted",
    )
    require(
        sum(
            row["supported"] == "TRUE"
            for row in inference_strip if row["inference_column"] == "maxT"
        ) == 18,
        "System-inference strip maxT support count drifted",
    )
    for row in read_tsv(
        SYSTEMS_INFERENCE_SOURCE / "promoted_source_checksums.tsv"
    ):
        path = ROOT / row["path"]
        require(path.is_file() and sha256(path) == row["sha256"],
                f"Promoted system-inference S4U source drift: {path}")

    asterisk_rows = read_tsv(SYSTEMS_ASTERISK_PROMOTION)
    require(
        len(asterisk_rows) == 1 and asterisk_rows[0]["callout"] == "S4U",
        "Strict-support S4U promotion manifest drifted",
    )
    asterisk_row = asterisk_rows[0]
    asterisk_panel = (ROOT / asterisk_row["destination"]).resolve()
    require(asterisk_panel == pooled_panel,
            "Strict-support S4U manifest and Figure 4 index disagree")
    require(asterisk_row["old_sha256"] == inference_row["new_sha256"],
            "Strict-support S4U does not supersede the inference-strip revision")
    asterisk_candidate_panel = (ROOT / asterisk_row["source"]).resolve()
    require(
        asterisk_candidate_panel.is_file()
        and sha256(asterisk_candidate_panel) == asterisk_row["new_sha256"],
        "Historical strict-support S4U candidate hash drift",
    )
    require(asterisk_row["complete_system_family"] == "43/43",
            "Strict-support S4U complete-family census drifted")
    require(asterisk_row["primary_bh_supported"] == "33/43",
            "Strict-support S4U primary BH census drifted")
    require(asterisk_row["strict_maxT_supported"] == "18/43",
            "Strict-support S4U maxT census drifted")
    require(asterisk_row["rendered_strict_markers"] == "18/43",
            "Strict-support S4U marker census drifted")
    require(asterisk_row["window_values_changed"] == "FALSE",
            "Strict-support S4U claims changed window values")
    require(asterisk_row["display_order_changed"] == "FALSE",
            "Strict-support S4U claims changed display order")
    require(asterisk_row["windows_inferentially_tested"] == "FALSE",
            "Strict-support S4U assigns inference to window means")
    asterisk_audit = read_tsv(SYSTEMS_ASTERISK_SOURCE / "provenance/audit.tsv")
    require(
        asterisk_audit and all(row["passed"] == "TRUE" for row in asterisk_audit),
        "Strict-support S4U source audit failed",
    )
    asterisk_strip = read_tsv(
        SYSTEMS_ASTERISK_SOURCE / "source_tables/system_inference_strip.tsv"
    )
    require(len(asterisk_strip) == 43 * 2,
            "Strict-support S4U source family drifted")
    strict_markers = [
        row for row in asterisk_strip if row["rendered_in_panel"] == "TRUE"
    ]
    require(
        len(strict_markers) == 18
        and all(
            row["inference_column"] == "maxT" and row["render_symbol"] == "*"
            for row in strict_markers
        ),
        "Strict-support S4U marker registry drifted",
    )
    require(
        not any(
            row["rendered_in_panel"] == "TRUE"
            for row in asterisk_strip if row["inference_column"] == "BH"
        ),
        "Primary BH calls must remain table-only in strict-support S4U",
    )
    for row in read_tsv(
        SYSTEMS_ASTERISK_SOURCE / "promoted_source_checksums.tsv"
    ):
        path = ROOT / row["path"]
        require(path.is_file() and sha256(path) == row["sha256"],
                f"Promoted strict-support S4U source drift: {path}")

    no_delta_rows = read_tsv(SYSTEMS_NO_DELTA_PROMOTION)
    require(
        len(no_delta_rows) == 1 and no_delta_rows[0]["callout"] == "S4U",
        "No-Delta S4U promotion manifest drifted",
    )
    no_delta_row = no_delta_rows[0]
    no_delta_panel = (ROOT / no_delta_row["destination"]).resolve()
    require(no_delta_panel == pooled_panel,
            "No-Delta S4U manifest and Figure 4 index disagree")
    require(no_delta_row["old_sha256"] == asterisk_row["new_sha256"],
            "No-Delta S4U does not supersede the asterisk revision")
    no_delta_candidate_panel = (ROOT / no_delta_row["source"]).resolve()
    require(
        no_delta_candidate_panel.is_file()
        and sha256(no_delta_candidate_panel) == no_delta_row["new_sha256"],
        "Historical no-Delta S4U candidate hash drift",
    )
    require(no_delta_row["complete_system_family"] == "43/43",
            "No-Delta S4U complete-family census drifted")
    require(no_delta_row["strict_maxT_supported"] == "18/43",
            "No-Delta S4U maxT census drifted")
    require(no_delta_row["rendered_strict_markers"] == "18/43",
            "No-Delta S4U marker census drifted")
    require(no_delta_row["delta_tile_displayed"] == "FALSE",
            "No-Delta S4U still declares a Delta tile")
    require(no_delta_row["delta_order_retained"] == "TRUE",
            "No-Delta S4U lost the frozen row order")
    require(no_delta_row["window_values_changed"] == "FALSE",
            "No-Delta S4U claims changed window values")
    require(no_delta_row["display_order_changed"] == "FALSE",
            "No-Delta S4U claims changed display order")
    require(no_delta_row["windows_inferentially_tested"] == "FALSE",
            "No-Delta S4U assigns inference to window means")
    no_delta_audit = read_tsv(SYSTEMS_NO_DELTA_SOURCE / "provenance/audit.tsv")
    require(
        no_delta_audit and all(row["passed"] == "TRUE" for row in no_delta_audit),
        "No-Delta S4U source audit failed",
    )
    no_delta_figures = read_tsv(SYSTEMS_NO_DELTA_SOURCE / "figure_manifest.tsv")
    require(
        len(no_delta_figures) == 1
        and no_delta_figures[0]["delta_tile_displayed"] == "FALSE"
        and no_delta_figures[0]["delta_order_retained"] == "TRUE",
        "No-Delta S4U figure contract drifted",
    )
    no_delta_strip = read_tsv(
        SYSTEMS_NO_DELTA_SOURCE / "source_tables/system_inference_strip.tsv"
    )
    no_delta_markers = [
        row for row in no_delta_strip if row["rendered_in_panel"] == "TRUE"
    ]
    require(
        len(no_delta_strip) == 43 * 2
        and len(no_delta_markers) == 18
        and all(
            row["inference_column"] == "maxT" and row["render_symbol"] == "*"
            for row in no_delta_markers
        ),
        "No-Delta S4U strict-marker registry drifted",
    )
    for row in read_tsv(
        SYSTEMS_NO_DELTA_SOURCE / "promoted_source_checksums.tsv"
    ):
        path = ROOT / row["path"]
        require(path.is_file() and sha256(path) == row["sha256"],
                f"Promoted no-Delta S4U source drift: {path}")

    caption_rows = read_tsv(SYSTEMS_CAPTION_FREE_PROMOTION)
    require(
        {row["callout"] for row in caption_rows}
        == {"4G", "S4Q", "S4R", "S4S", "S4T", "S4U"},
        "Caption-free molecular-systems promotion callouts drifted",
    )
    visual_by_callout = {row["callout"]: row for row in visual_rows}
    for row in caption_rows:
        callout = row["callout"]
        indexed = (
            FIG4 / indexed_visual_rows[callout]["candidate_relative_path"]
        ).resolve()
        manifested = (ROOT / row["destination"]).resolve()
        require(indexed == manifested,
                f"{callout} caption-free index and manifest disagree")
        predecessor_sha = (
            no_delta_row["new_sha256"] if callout == "S4U"
            else visual_by_callout[callout]["new_sha256"]
        )
        require(row["old_sha256"] == predecessor_sha,
                f"{callout} caption-free predecessor hash drift")
        require(sha256(indexed) == row["new_sha256"],
                f"{callout} caption-free active hash drift")
        source = (ROOT / row["source"]).resolve()
        require(source.is_file() and sha256(source) == row["new_sha256"],
                f"{callout} caption-free candidate hash drift")
        require(row["embedded_explanatory_text"] == "FALSE",
                f"{callout} still declares embedded explanatory text")
        require(row["data_values_changed"] == "FALSE",
                f"{callout} caption-free revision claims changed values")
        require(row["geometry_changed"] == "FALSE",
                f"{callout} caption-free revision claims changed geometry")
        require(row["display_order_changed"] == "FALSE",
                f"{callout} caption-free revision claims changed order")
    caption_family_figures = read_tsv(
        SYSTEMS_CAPTION_FREE_SOURCE / "family/figure_manifest.tsv"
    )
    caption_s4u_figures = read_tsv(
        SYSTEMS_CAPTION_FREE_SOURCE / "s4u/figure_manifest.tsv"
    )
    require(
        len(caption_family_figures) == 6
        and len(caption_s4u_figures) == 1
        and all(
            row["embedded_explanatory_text"] == "FALSE"
            for row in caption_family_figures + caption_s4u_figures
        ),
        "Caption-free source figure contract drifted",
    )
    for row in read_tsv(
        SYSTEMS_CAPTION_FREE_SOURCE / "promoted_source_checksums.tsv"
    ):
        path = ROOT / row["path"]
        require(path.is_file() and sha256(path) == row["sha256"],
                f"Promoted caption-free source drift: {path}")

    molecular_figures = read_tsv(MOLECULAR / "figures/figure_manifest.tsv")
    require(molecular_figures, "Molecular figure manifest is empty")
    require(
        all(Path(row["path"]).is_file() for row in molecular_figures),
        "A molecular figure is missing",
    )
    require(
        all(row["current_figure3_modified"] == "FALSE" for row in molecular_figures),
        "A continuum figure claims to modify Figure 3",
    )

    comparison_figures = read_tsv(COMPARISON / "figure_manifest.tsv")
    require(len(comparison_figures) == 11, "Comparison figure family drifted")
    for row in comparison_figures:
        path = Path(row["path"])
        require(path.is_file() and path.stat().st_size > 0, f"Missing figure: {path}")
        require(sha256(path) == row["sha256"], f"Figure hash drift: {path}")
        require(
            row["claim_boundary"]
            == "same-substrate comparison; continuum windows descriptive only",
            f"Claim-boundary drift: {path}",
        )

    print("CONTINUUM_INTEGRATION\tPASS")
    print(f"CORE_VALIDATION\t{core_summary['status']}")
    print(f"MOLECULAR_VALIDATION\t{molecular_summary['validation_pass']}")
    print(f"FIGURE4F\t{indexed_4f}")
    print(f"FIGURE4H\t{indexed_4h}")
    print("FIGURE4H_SYSTEMS\t43_complete;18_maxT_outlines;15_fixed_labels;36_same_direction")
    print(f"FIGURE4G_S4Q_TO_S4U\t{len(visual_rows)}")
    print(f"FIGURE_S4U_POOLED_SYSTEMS\t{pooled_panel}")
    print("FIGURE_S4U_DONOR_SUPPORT\t43_complete;18_maxT_asterisks;33_BH_table_only;Delta_hidden")
    print("FIGURE4G_S4Q_TO_S4U_EMBEDDED_PROSE\tremoved;external_caption_retained")
    print(f"ATLAS_SYSTEMS\t{atlas_metrics['membership_defined_systems']}")
    print(f"MOLECULAR_FIGURES\t{len(molecular_figures)}")
    print(f"COMPARISON_FIGURES_HASHED\t{len(comparison_figures)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
