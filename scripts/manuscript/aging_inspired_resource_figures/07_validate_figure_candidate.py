#!/usr/bin/env python3
"""Independently validate the aging-inspired Figure 4–6 candidate."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
AMBIENT = (
    ROOT
    / "Analysis/SingleCell/candidates"
    / "ambient-program-recalibration-all117-candidate-2026-08-15-v2"
)
CROSS = (
    ROOT
    / "Analysis/SingleCell/candidates"
    / "cross-lineage-specificity-complete-atlas-candidate-2026-08-15-v2"
)
OUT = Path(os.environ["FIG_CAND_ROOT"])
SOURCE = OUT / "source_tables"
PANELS = OUT / "panels"
PROOFS = OUT / "proofs"
VALIDATION = OUT / "validation_report.tsv"
MANIFEST = OUT / "output_manifest.tsv"
VALIDATED = OUT / "VALIDATED"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def pdf_font_subtypes(path: Path) -> str:
    if shutil.which("pdffonts"):
        return run(["pdffonts", str(path)])
    require(shutil.which("gs") is not None, "pdffonts and gs are both unavailable")
    pdf_path = str(path.resolve()).replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")
    postscript = (
        f"({pdf_path}) (r) file runpdfbegin 1 pdfgetpage /Resources get exec "
        "/Font get exec { exch pop exec /Subtype get == } forall quit"
    )
    return run(["gs", "-q", "-dNOSAFER", "-dNODISPLAY", "-c", postscript])


def refuse(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")


def read_tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(SOURCE / name, sep="\t")


def run(command: list[str]) -> str:
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return result.stdout + result.stderr


def close(left: pd.Series, right: pd.Series) -> bool:
    return bool(
        np.allclose(
            pd.to_numeric(left, errors="coerce"),
            pd.to_numeric(right, errors="coerce"),
            equal_nan=True,
            atol=1e-12,
            rtol=1e-10,
        )
    )


def main() -> None:
    for path in (VALIDATION, MANIFEST, VALIDATED):
        refuse(path)
    rows: list[dict[str, object]] = []

    def add(check: str, passed: bool, detail: object) -> None:
        rows.append({"check": check, "passed": bool(passed), "detail": str(detail)})
        require(bool(passed), f"validation failed: {check}: {detail}")

    ambient_validation = pd.read_csv(AMBIENT / "results/validation_report.tsv", sep="\t")
    add("ambient_candidate_validated", (AMBIENT / "VALIDATED").exists(), AMBIENT / "VALIDATED")
    recompute = ambient_validation[ambient_validation["check"] == "all_234_effect_rows_recomputed"]
    add(
        "ambient_donor_level_effects_independently_recomputed",
        len(recompute) == 1 and bool(recompute.iloc[0]["passed"]),
        recompute.iloc[0]["detail"] if len(recompute) else "missing",
    )

    effects = pd.read_csv(AMBIENT / "results/ambient_program_effects.tsv", sep="\t")
    effects = effects[effects["analysis_universe"] == "complete_case_common_universe"].copy()
    skyline = read_tsv("fig4c_all117_disease_skyline.tsv")
    add("fig4c_rows_117", len(skyline) == 117, len(skyline))
    joined = skyline.merge(effects, on="program_uid", suffixes=("_plot", "_source"), validate="one_to_one")
    add("fig4c_beta_exact", close(joined["disease_beta"], joined["raw_beta"]), 117)
    add("fig4c_q_exact", close(joined["disease_qvalue"], joined["raw_qvalue"]), 117)
    expected_evidence = -np.log10(np.maximum(joined["raw_qvalue"].astype(float), np.finfo(float).tiny))
    add("fig4c_evidence_transform_exact", close(joined["family_evidence"], expected_evidence), 117)

    transport = read_tsv("fig4d_ambient_effect_transport.tsv")
    add("fig4d_rows_117", len(transport) == 117, len(transport))
    untestable_transport = transport[transport["corrected_beta"].isna()]
    add(
        "fig4d_11_tcell_programs_visible_as_untestable",
        len(untestable_transport) == 11
        and set(untestable_transport["cell_type"]) == {"tcells"}
        and set(untestable_transport["evidence_state"]) == {"untestable"},
        len(untestable_transport),
    )
    joined = transport.merge(effects, on="program_uid", suffixes=("_plot", "_source"), validate="one_to_one")
    for prefix in ("raw", "corrected", "delta"):
        for metric in ("beta", "ci_low", "ci_high", "qvalue", "hc3_qvalue"):
            add(
                f"fig4d_exact:{prefix}_{metric}",
                close(joined[f"{prefix}_{metric}_plot"], joined[f"{prefix}_{metric}_source"]),
                117,
            )
    status = read_tsv("fig4d_dataset_status_strip.tsv")
    add("fig4d_dataset_rows_7", len(status) == 7, len(status))
    ambient_qc = pd.read_csv(AMBIENT / "results/ambient_dataset_lineage_qc.tsv", sep="\t")
    expected_failed = set(
        ambient_qc.loc[
            ambient_qc["correction_status"] != "corrected", "dataset"
        ].astype(str)
    )
    plotted_failed = set(
        status.loc[status["correction_status"] != "corrected", "dataset"].astype(str)
    )
    add(
        "fig4d_every_failed_dataset_visible_with_reason",
        plotted_failed == expected_failed
        and status.loc[
            status["correction_status"] != "corrected", "failure_reason"
        ].fillna("").str.len().gt(0).all(),
        sorted(plotted_failed),
    )
    gse = status[status["dataset"] == "GSE189600"]
    add(
        "GSE189600_visible_passthrough",
        len(gse) == 1
        and gse.iloc[0]["correction_status"] == "uncorrected_passthrough"
        and "passthrough" in gse.iloc[0]["explicit_passthrough_label"],
        gse.to_dict("records"),
    )

    lineage = read_tsv("fig4e_same_atlas_lineage_specificity.tsv")
    add("fig4e_rows_two_by_six", len(lineage) == 12, len(lineage))
    contrasts = pd.read_csv(CROSS / "results/hero_lineage_contrasts.tsv", sep="\t")
    contrasts = contrasts[
        (contrasts["annotation_filter"] == "all_annotated_cells")
        & (contrasts["minimum_cells"] == 50)
    ]
    nonreference = lineage[lineage["comparison_lineage"] != "Hepatocytes"]
    joined = nonreference.merge(
        contrasts,
        on=["program_uid", "comparison_lineage"],
        suffixes=("_plot", "_source"),
        validate="one_to_one",
    )
    for metric in ("beta", "ci_low", "ci_high", "qvalue", "hc3_qvalue"):
        add(f"fig4e_exact:{metric}", close(joined[f"{metric}_plot"], joined[f"{metric}_source"]), 10)
    promotion = pd.read_csv(CROSS / "results/hero_lineage_promotion.tsv", sep="\t")
    destination = promotion["destination"].unique()
    add(
        "fig4e_destination_predeclared",
        len(destination) == 1 and set(lineage["destination"]) == {destination[0]},
        destination,
    )

    sc4f = read_tsv("fig4f_singlecell_estimates.tsv")
    add("fig4f_singlecell_rows_4", len(sc4f) == 4, len(sc4f))
    for estimate_type, group in sc4f.groupby("estimate_type"):
        source = effects[effects["program_uid"].isin(group["program_uid"])]
        joined = group.merge(source, on="program_uid", validate="one_to_one")
        for metric in ("beta", "ci_low", "ci_high", "qvalue"):
            add(
                f"fig4f_singlecell_exact:{estimate_type}_{metric}",
                close(joined[metric], joined[f"{estimate_type}_{metric}"]),
                len(joined),
            )
    bulk_plot = read_tsv("fig4f_bulk_stage_transport.tsv")
    bulk_source = pd.read_csv(
        ROOT
        / "figures/candidates/fig4-stage-terminology-corrected-2026-08-13-v5"
        / "source_tables/fig4e_bulk_projection.tsv",
        sep="\t",
    )
    add("fig4f_bulk_rows_8", len(bulk_plot) == 8, len(bulk_plot))
    joined = bulk_plot.merge(
        bulk_source,
        on=["program_uid", "stage"],
        suffixes=("_plot", "_source"),
        validate="one_to_one",
    )
    for metric in ("effect", "direction_agreement", "bh_weight_fraction"):
        add(f"fig4f_bulk_exact:{metric}", close(joined[f"{metric}_plot"], joined[f"{metric}_source"]), 8)

    spatial_plot = read_tsv("fig5f_spatial_map_and_matched_null.tsv")
    spatial_source = pd.read_csv(
        ROOT
        / "Analysis/Multimodal_Program_Projection/candidates"
        / "spatial-impact-figures-candidate-2026-08-11/data/fig4f_matched_null_discrimination.tsv",
        sep="\t",
    )
    add("fig5f_rows_4", len(spatial_plot) == 4, len(spatial_plot))
    joined = spatial_plot.merge(
        spatial_source,
        on=["program_id", "dataset"],
        suffixes=("_plot", "_source"),
        validate="one_to_one",
    )
    for metric in ("residual_moran_i", "primary_qvalue_spatial", "q050", "q500", "q950"):
        add(f"fig5f_exact:{metric}", close(joined[f"{metric}_plot"], joined[f"{metric}_source"]), 4)
    map_paths = set(spatial_plot["map_panel_path"])
    add(
        "fig5f_map_hash_exact",
        len(map_paths) == 1
        and Path(next(iter(map_paths))).exists()
        and set(spatial_plot["map_panel_sha256"]) == {sha256(Path(next(iter(map_paths))))},
        map_paths,
    )

    fingerprints = read_tsv("fig6c_evidence_fingerprints.tsv")
    nodes = read_tsv("fig6c_evidence_nodes.tsv")
    add("fig6c_three_examples", fingerprints["example"].nunique() == 3, fingerprints["example"].nunique())
    add("fig6c_nodes_15", len(nodes) == 15, len(nodes))
    add("fig6c_five_nodes_per_example", set(nodes.groupby("example").size()) == {5}, nodes.groupby("example").size().to_dict())
    add(
        "fig6c_genetics_blocked",
        set(
            fingerprints.loc[
                fingerprints["node"] == "Inherited shared signal", "state"
            ]
        )
        == {"blocked_pending_corrected_coloc"},
        fingerprints.loc[fingerprints["node"] == "Inherited shared signal", "state"].tolist(),
    )
    add(
        "fig6c_no_combined_score_or_rank",
        not any(
            token in column.lower()
            for column in fingerprints.columns
            for token in ("combined_score", "vote_count", "rank")
        ),
        list(fingerprints.columns),
    )
    source_hash_ok = True
    for path_string, expected in fingerprints[["source_path", "source_sha256"]].drop_duplicates().itertuples(index=False):
        path = Path(path_string)
        source_hash_ok &= path.exists() and sha256(path) == expected
    add("fig6c_all_source_hashes_exact", source_hash_ok, fingerprints["source_path"].nunique())
    add(
        "fig6c_native_units_complete",
        fingerprints["unit"].notna().all() and fingerprints["unit"].str.len().gt(0).all(),
        fingerprints["unit"].nunique(),
    )
    add(
        "fig6c_alternatives_and_experiments_complete",
        fingerprints["unresolved_alternative"].str.len().gt(0).all()
        and fingerprints["next_experiment"].str.len().gt(0).all(),
        len(fingerprints),
    )
    experiment_router = read_tsv("fig6d_experiment_router.tsv")
    expected_router = fingerprints[
        ["example", "node", "unresolved_alternative", "next_experiment"]
    ].drop_duplicates().sort_values(
        ["example", "node", "unresolved_alternative", "next_experiment"]
    ).reset_index(drop=True)
    observed_router = experiment_router.sort_values(
        ["example", "node", "unresolved_alternative", "next_experiment"]
    ).reset_index(drop=True)
    add(
        "fig6d_every_alternative_routes_to_discriminating_experiment",
        observed_router.equals(expected_router),
        len(observed_router),
    )

    supplement = read_tsv("supplement_asset_manifest.tsv")
    add("supplement_asset_rows_11", len(supplement) == 11, len(supplement))
    supplement_ok = True
    for row in supplement.itertuples(index=False):
        source_path = Path(row.source_path)
        supplement_destination = OUT / row.relative_destination
        supplement_ok &= (
            source_path.exists()
            and supplement_destination.exists()
            and sha256(source_path) == row.source_sha256
            and sha256(supplement_destination) == row.destination_sha256
            and row.source_sha256 == row.destination_sha256
        )
    add("supplement_assets_byte_identical", supplement_ok, len(supplement))

    expected_panel = (
        "fig4e_same_atlas_lineage_specificity.pdf"
        if destination[0] == "Figure_4E"
        else "figS4e_same_atlas_lineage_specificity.pdf"
    )
    expected_pdfs = {
        "fig4a_analysis_logic.pdf",
        "fig4c_all117_disease_skyline.pdf",
        "fig4d_ambient_effect_transport.pdf",
        expected_panel,
        "fig4f_tissue_state_transport.pdf",
        "fig5f_spatial_maps_and_matched_null.pdf",
        "fig6c_evidence_ribbons.pdf",
        "fig6d_experiment_router.pdf",
    }
    observed_pdfs = {path.name for path in PANELS.glob("*.pdf")}
    add("individual_panel_roster_exact", observed_pdfs == expected_pdfs, sorted(observed_pdfs))
    proof_pdfs = sorted(PROOFS.glob("*.pdf"))
    add("composite_proofs_two", len(proof_pdfs) == 2, [path.name for path in proof_pdfs])
    for pdf in sorted(PANELS.glob("*.pdf")) + proof_pdfs:
        info = run(["pdfinfo", str(pdf)])
        pages = [line for line in info.splitlines() if line.startswith("Pages:")]
        add(f"one_page:{pdf.name}", len(pages) == 1 and pages[0].split()[-1] == "1", pages)
        fonts = pdf_font_subtypes(pdf)
        add(f"no_type3:{pdf.name}", "Type 3" not in fonts and "Type3" not in fonts, fonts.replace("\n", " | ")[:400])
        run(
            [
                "gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage",
                str(pdf),
            ]
        )
        add(f"ghostscript:{pdf.name}", True, "pass")

    render_script = ROOT / "scripts/manuscript/aging_inspired_resource_figures/08_render_compact_panels.R"
    render_text = render_script.read_text(encoding="utf-8")
    add(
        "render_contract_six_point_plain",
        "GEOM_TEXT_6PT" in render_text
        and (
            "theme_masld()" in render_text
            or "theme_masld_compact()" in render_text
        )
        and "face = \"plain\"" in render_text,
        render_script,
    )
    add("useDingbats_false", "useDingbats = FALSE" in render_text, render_script)
    add(
        "compact_title_free_panel_contract",
        "plot.title = element_blank()" in render_text
        and 'fig4a_analysis_logic.pdf"), 3.0, 1.25' in render_text
        and 'fig4d_ambient_effect_transport.pdf"), 2.25, 2.45' in render_text
        and 'fig6c_evidence_ribbons.pdf"), 5.0, 2.55' in render_text,
        render_script,
    )

    validation = pd.DataFrame(rows)
    validation.to_csv(VALIDATION, sep="\t", index=False)
    release_files = sorted(
        path
        for path in OUT.rglob("*")
        if path.is_file() and path not in {MANIFEST, VALIDATED}
    )
    release = pd.DataFrame(
        [
            {
                "relative_path": str(path.relative_to(OUT)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in release_files
        ]
    )
    release.to_csv(MANIFEST, sep="\t", index=False)
    VALIDATED.write_text(
        "status\tvalidated_candidate\n"
        f"validation_report_sha256\t{sha256(VALIDATION)}\n"
        f"output_manifest_sha256\t{sha256(MANIFEST)}\n"
        f"ambient_validated_sha256\t{sha256(AMBIENT / 'VALIDATED')}\n",
        encoding="utf-8",
    )
    print(f"PASS {len(rows)} checks; {len(release)} candidate artifacts")


if __name__ == "__main__":
    main()
