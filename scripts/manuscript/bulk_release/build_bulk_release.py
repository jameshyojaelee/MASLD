#!/usr/bin/env python3
"""Assemble one immutable, bulk-only five-cohort manuscript release.

This tool packages prevalidated inputs; it never refits a model or falls back to
the historical integration directory.  It refuses an existing output path so a
named release cannot be overwritten.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


RELEASE_ID = "2026-09-09-bulk-r1"
N_GENES = 23_370
# Pooled counts of the adopted 2026-09-09 release. A corrected-control refit
# (scripts/manuscript/resource_f_five/refit_corrected_controls.R) ships
# expected_counts.json in its pooled root, which replaces these.
LEGACY_EXPECTED = {
    "n_samples": 844,
    "n_control": 157,
    "n_disease": 687,
    "n_genes": N_GENES,
    "canonical": {"n": 1347, "up": 1003, "down": 344},
    "treat": {"n": 1616, "up": 1144, "down": 472},
}


def load_expected(path: Path) -> dict:
    """Expected counts from expected_counts.json, or the adopted-release values."""
    if not path.is_file():
        return LEGACY_EXPECTED
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_table(path: Path):
    return gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")


def rows(path: Path):
    # deg_results.csv is comma-separated; every other source table is TSV.
    delimiter = "," if path.name.endswith(".csv") else "\t"
    with open_table(path) as handle:
        yield from csv.DictReader(handle, delimiter=delimiter)


def read_rows(path: Path) -> list[dict[str, str]]:
    return list(rows(path))


def write_rows(path: Path, fieldnames: list[str], data: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def finite(value: str) -> float:
    if value in {"", "NA", "NaN", "nan"}:
        raise ValueError("missing numeric value")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite numeric value: {value}")
    return result


def bh_max_error(data: list[dict[str, str]], p_col: str, q_col: str, n_genes: int = N_GENES) -> float:
    require(len(data) == n_genes, f"BH family has {len(data)} rows, expected {n_genes}")
    p_values = [finite(row[p_col]) for row in data]
    q_values = [finite(row[q_col]) for row in data]
    order = sorted(range(n_genes), key=p_values.__getitem__, reverse=True)
    running = 1.0
    expected = [0.0] * n_genes
    for rank_desc, index in enumerate(order, start=1):
        rank_asc = n_genes - rank_desc + 1
        running = min(running, n_genes * p_values[index] / rank_asc)
        expected[index] = running
    return max(abs(observed - calculated) for observed, calculated in zip(q_values, expected))


def copy_file(source: Path, destination: Path) -> dict[str, object]:
    require(source.is_file(), f"Required source is absent: {source}")
    require(not destination.exists(), f"Refusing to overwrite release artifact: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return {
        "release_path": str(destination),
        "source_path": str(source.resolve()),
        "size_bytes": destination.stat().st_size,
        "sha256": sha256(destination),
    }


def copy_named_files(source_dir: Path, destination_dir: Path, names: list[str]) -> list[dict[str, object]]:
    return [copy_file(source_dir / name, destination_dir / name) for name in names]


def validate_manifested_directory(source_dir: Path, manifest_name: str = "output_manifest.tsv") -> None:
    manifest = source_dir / manifest_name
    require(manifest.is_file(), f"Missing source output manifest: {manifest}")
    for row in rows(manifest):
        artifact = source_dir / row["relative_path"]
        require(artifact.is_file(), f"Manifested source artifact is absent: {artifact}")
        require(sha256(artifact) == row["sha256"], f"Manifest hash mismatch: {artifact}")


def validate_pooled(path: Path, expected: dict = LEGACY_EXPECTED) -> dict[str, object]:
    data = read_rows(path)
    n_genes = expected["n_genes"]
    require(len(data) == n_genes, f"Pooled table has {len(data)} rows, expected {n_genes}")
    genes = [row["gene"] for row in data]
    require(len(set(genes)) == n_genes, "Pooled table duplicates a versioned gene")
    primary = [row for row in data if finite(row["padj"]) < 0.05 and abs(finite(row["logFC"])) > 0.50]
    treat = [row for row in data if finite(row["treat_fdr"]) < 0.05]
    canonical, treat_expected = expected["canonical"], expected["treat"]
    require(len(primary) == canonical["n"], f"Primary pooled DEG count is {len(primary)}, expected {canonical['n']}")
    require(sum(finite(row["logFC"]) > 0 for row in primary) == canonical["up"], "Pooled up-DEG count drift")
    require(sum(finite(row["logFC"]) < 0 for row in primary) == canonical["down"], "Pooled down-DEG count drift")
    require(len(treat) == treat_expected["n"], f"TREAT pooled DEG count is {len(treat)}, expected {treat_expected['n']}")
    require(sum(finite(row["logFC"]) > 0 for row in treat) == treat_expected["up"], "TREAT up-DEG count drift")
    require(sum(finite(row["logFC"]) < 0 for row in treat) == treat_expected["down"], "TREAT down-DEG count drift")
    return {
        "pooled_point_null_max_abs_BH_error": bh_max_error(data, "P.Value", "padj", n_genes),
        "pooled_treat_max_abs_BH_error": bh_max_error(data, "treat_p", "treat_fdr", n_genes),
        "pooled_primary_degs": len(primary),
        "pooled_treat_degs": len(treat),
    }


def validate_stage(path: Path, n_genes: int = N_GENES) -> dict[str, object]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows(path):
        key = (row["axis"], row["contrast"])
        grouped[key].append(row)
    require(len(grouped) == 14, f"Stage table has {len(grouped)} fitted families, expected 14")
    required_axes = {"fibrosis_adjacent", "fibrosis", "NAS_adjacent", "NAS"}
    require({key[0] for key in grouped} == required_axes, "Stage axes drift")
    audits = []
    for (axis, contrast), family in sorted(grouped.items()):
        genes = [row["gene_id_versioned"] for row in family]
        require(len(set(genes)) == n_genes, f"Duplicate or missing genes in {axis}/{contrast}")
        require(all(finite(row["CI_low"]) <= finite(row["logFC"]) <= finite(row["CI_high"]) for row in family),
                f"Invalid pointwise interval in {axis}/{contrast}")
        error = bh_max_error(family, "P.Value", "FDR", n_genes)
        require(error < 1e-10, f"BH mismatch in {axis}/{contrast}: {error}")
        template = family[0]
        audits.append({
            "analysis": "stage",
            "axis": axis,
            "contrast": contrast,
            "n_genes": len(family),
            "n_reference": template["n_reference"],
            "n_comparison": template["n_comparison"],
            "n_cohorts": template["n_cohorts"],
            "bh_family_size": template["bh_family_size"],
            "max_abs_BH_error": f"{error:.3g}",
            "intervals": "finite_pointwise_95pct",
        })
    return {"stage_families": len(grouped), "stage_audits": audits}


def validate_f0(path: Path, summary_path: Path, expected: dict = LEGACY_EXPECTED) -> dict[str, object]:
    n_genes = expected["n_genes"]
    legacy = "f0_reference" not in expected
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows(path):
        grouped[(row["arm"], row["contrast"])].append(row)
    if legacy:
        require(len(grouped) == 23, f"F0 sensitivity has {len(grouped)} estimable fits, expected 23")
    audits = []
    for (arm, contrast), family in sorted(grouped.items()):
        require(len({row["gene_id_versioned"] for row in family}) == n_genes,
                f"F0 sensitivity family is incomplete: {arm}/{contrast}")
        require(all(finite(row["CI_low"]) <= finite(row["logFC"]) <= finite(row["CI_high"]) for row in family),
                f"Invalid F0 sensitivity interval: {arm}/{contrast}")
        error = bh_max_error(family, "P.Value", "FDR", n_genes)
        require(error < 1e-10, f"F0 sensitivity BH mismatch: {arm}/{contrast}")
        audits.append({"analysis": "F0_sensitivity", "axis": "fibrosis", "contrast": contrast,
                       "arm": arm, "n_genes": len(family), "max_abs_BH_error": f"{error:.3g}"})
    summary = {(row["arm"], row["contrast"]): row for row in rows(summary_path)}
    if legacy:
        reference = {
            ("A_all_F0", "F0_to_F1"): (1471, 31, 73),
            ("B_disease_only_F0", "F0_to_F1"): (7, 0, 73),
            ("C_all_F0_plus_term", "F0_to_F1"): (258, 31, 73),
        }
        for key, values in reference.items():
            require(key in summary, f"Missing F0 summary row: {key}")
            observed = summary[key]
            require((int(observed["n_BH_sig"]), int(observed["n_ref_control"]), int(observed["n_ref_disease"])) == values,
                    f"F0 sensitivity summary drift for {key}")
    else:
        # A corrected refit fixes the F0 reference composition from its labels;
        # the DEG counts are outcomes of this run, reported, not asserted.
        for arm, counts in expected["f0_reference"].items():
            key = (arm, "F0_to_F1")
            require(key in summary, f"Missing F0 summary row: {key}")
            observed = summary[key]
            require((int(observed["n_ref_control"]), int(observed["n_ref_disease"])) ==
                    (counts["n_ref_control"], counts["n_ref_disease"]),
                    f"F0 reference composition drift for {key}")
    require(("C_all_F0_plus_term", "F3_to_F4") not in grouped,
            "Non-estimable F3-to-F4 control-status sensitivity must remain unavailable")
    return {"f0_families": len(grouped), "f0_audits": audits}


def composition_summaries(composition: Path, stage_samples: Path) -> list[dict[str, object]]:
    comp = read_rows(composition)
    require(comp and {"sample_id", "dataset", "Hepatocytes"}.issubset(comp[0]), "Composition schema drift")
    fraction_columns = [column for column in comp[0] if column not in {"sample_id", "dataset"}]
    comp_by_sample: dict[str, dict[str, str]] = {}
    for row in comp:
        sample = row["sample_id"]
        require(sample not in comp_by_sample, f"Duplicate composition sample: {sample}")
        values = [finite(row[column]) for column in fraction_columns]
        require(all(0 <= value <= 1 for value in values) and abs(sum(values) - 1) <= 1e-6,
                f"Invalid composition fractions for {sample}")
        comp_by_sample[sample] = row
    groups: dict[tuple[str, str, str, str], list[tuple[float, float]]] = defaultdict(list)
    seen: set[tuple[str, str, str]] = set()
    for row in rows(stage_samples):
        identity = (row["axis"], row["contrast"], row["sample_id"])
        require(identity not in seen, f"Duplicate stage sample membership: {identity}")
        seen.add(identity)
        sample = row["sample_id"]
        require(sample in comp_by_sample, f"Stage sample missing accepted composition: {sample}")
        source = comp_by_sample[sample]
        require(source["dataset"] == row["dataset"], f"Composition cohort mismatch: {sample}")
        hepatocyte = finite(source["Hepatocytes"])
        npc = sum(finite(source[column]) for column in fraction_columns if column != "Hepatocytes")
        groups[(row["axis"], row["contrast"], row["dataset"], row["group"])].append((hepatocyte, npc))
    summary: list[dict[str, object]] = []
    for (axis, contrast, dataset, group), values in sorted(groups.items()):
        for compartment, index in (("hepatocyte", 0), ("non_parenchymal", 1)):
            measurements = [value[index] for value in values]
            n = len(measurements)
            se = 0.0 if n == 1 else math.sqrt(sum((value - sum(measurements) / n) ** 2 for value in measurements) / (n - 1)) / math.sqrt(n)
            summary.append({
                "summary_level": "cohort",
                "axis": axis,
                "contrast": contrast,
                "dataset": dataset,
                "group": group,
                "compartment": compartment,
                "mean": f"{sum(measurements) / n:.17g}",
                "SE": f"{se:.17g}",
                "n_participants": n,
                "estimand": "unadjusted participant mean plus SE",
            })
    return summary


def validate_pdf(path: Path) -> None:
    require(path.is_file(), f"Missing required rendered PDF: {path}")
    info = subprocess.run(["pdfinfo", str(path)], check=True, text=True, capture_output=True).stdout
    require("Pages:           1" in info, f"Expected one PDF page: {path}")
    raw = path.read_bytes()
    require(b"/Subtype /Type3" not in raw and b"/FontType 3" not in raw, f"Type 3 font: {path}")
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", str(path)], check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pooled", type=Path, required=True)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--f0", type=Path, required=True)
    parser.add_argument("--composition", type=Path, required=True)
    parser.add_argument("--lncrna", type=Path, required=True)
    parser.add_argument("--lncrna-biotype", type=Path, required=True)
    parser.add_argument("--figure-candidate", type=Path, required=True)
    parser.add_argument("--fig3g", type=Path, required=True)
    parser.add_argument("--fig4-candidate", type=Path, required=True)
    parser.add_argument("--release-id", default=RELEASE_ID)
    args = parser.parse_args()
    release_id = args.release_id
    expected_path = args.pooled / "expected_counts.json"
    expected = load_expected(expected_path)
    n_genes = expected["n_genes"]

    output = args.output.resolve()
    require(output.name == release_id, f"Bulk release ID must be {release_id}")
    require(not output.exists(), f"Refusing to overwrite named release: {output}")
    for source in (args.pooled, args.stage, args.f0, args.composition, args.lncrna,
                   args.lncrna_biotype, args.figure_candidate, args.fig3g, args.fig4_candidate):
        require(source.exists(), f"Required input is absent: {source}")

    pooled = args.pooled / "deg_results.csv"
    stage = args.stage / "stage_extension_all_gene_results.tsv"
    stage_samples = args.stage / "stage_extension_sample_manifest.tsv"
    f0 = args.f0 / "f0_arm_results.tsv.gz"
    f0_summary = args.f0 / "f0_arm_summary.tsv"
    validate_manifested_directory(args.pooled)
    validate_manifested_directory(args.stage)
    validate_manifested_directory(args.lncrna)
    validate_manifested_directory(args.lncrna_biotype)
    pooled_checks = validate_pooled(pooled, expected)
    stage_checks = validate_stage(stage, n_genes)
    f0_checks = validate_f0(f0, f0_summary, expected)
    comp_summary = composition_summaries(args.composition, stage_samples)

    figure_pdfs = [
        args.figure_candidate / "figure3/panels/fig3a_cohort_metadata_matrix.pdf",
        args.figure_candidate / "figure3/panels/fig3b_nas_fib_grid.pdf",
        args.figure_candidate / "figure3/panels/fig3c_cohort_alluvial.pdf",
        args.figure_candidate / "figure3/panels/fig3d_pca_fibrosis_gradient.pdf",
        args.figure_candidate / "figure3/panels/fig3e_stage_remodeling.pdf",
        args.figure_candidate / "supplementary/figureS3/panels/figs3c_cohort_robustness.pdf",
        args.figure_candidate / "supplementary/figureS3/panels/figs3_endpoint_deg_concordance.pdf",
        args.figure_candidate / "supplementary/figureS3/panels/figs3_stage_remodeling_full.pdf",
    ]
    for pdf in figure_pdfs:
        validate_pdf(pdf)

    output.mkdir(parents=True)
    inventory: list[dict[str, object]] = []
    def add(source: Path, relative: str, purpose: str, unit: str, state: str) -> None:
        record = copy_file(source, output / relative)
        record.update({"purpose": purpose, "biological_unit": unit, "state": state})
        record["release_path"] = str(Path(record["release_path"]).relative_to(output))
        inventory.append(record)

    for name in ["deg_results.csv", "treat_table.tsv.gz", "model_design.tsv", "model_input_manifest.tsv",
                 "VALIDATED.json", "validation/validation_report.json", "validation/validation_checks.tsv"] + (
                     ["expected_counts.json"] if expected_path.is_file() else []):
        add(args.pooled / name, f"artifacts/pooled/{Path(name).name}", "five-cohort pooled disease-versus-control", "one selected biopsy record", "adopted_bulk_only")
    for name in ["stage_extension_all_gene_results.tsv", "stage_all_gene_results.tsv", "stage_extension_design_audit.tsv",
                 "stage_extension_sample_manifest.tsv", "stage_sample_manifest.tsv", "five_cohort_sample_manifest.tsv",
                 "five_cohort_summary.tsv", "cohort_disease_all_gene_results.tsv", "prespecified_nmf_program_activity.tsv",
                 "input_manifest.tsv", "output_manifest.tsv", "sessionInfo.txt"]:
        add(args.stage / name, f"artifacts/stage/{name}", "cross-sectional stage or NAS remodeling", "one selected biopsy record", "adopted_bulk_only")
    for name in ["f0_arm_results.tsv.gz", "f0_arm_summary.tsv", "f0_composition_audit.tsv", "f0_arm_hero_genes.tsv",
                 "f0_arm_effect_agreement.tsv", "sessionInfo.txt"]:
        add(args.f0 / name, f"artifacts/f0_sensitivity/{name}", "F0 reference-composition sensitivity", "one selected biopsy record", "adopted_bulk_only")
    add(args.composition, "artifacts/composition/accepted_composition.tsv.gz", "accepted deconvolution fractions", "one selected biopsy record", "adopted_bulk_only")
    for source_dir, destination, purpose in [
        (args.lncrna, "artifacts/lncrna", "approved lncRNA disease and stage outputs"),
        (args.lncrna_biotype, "artifacts/lncrna_biotype", "lncRNA biotype sensitivity"),
    ]:
        for source in sorted(source_dir.iterdir()):
            if source.is_file():
                add(source, f"{destination}/{source.name}", purpose, "versioned gene", "adopted_bulk_only")
    for source, relative, purpose in [
        (figure_pdfs[0], "figures/figure3/fig3a_cohort_metadata_matrix.pdf", "Figure 3A"),
        (figure_pdfs[1], "figures/figure3/fig3b_nas_fib_grid.pdf", "Figure 3B"),
        (figure_pdfs[2], "figures/figure3/fig3c_cohort_alluvial.pdf", "Figure 3C"),
        (figure_pdfs[3], "figures/figure3/fig3d_pca_fibrosis_gradient.pdf", "Figure 3D endpoint sensitivity"),
        (figure_pdfs[4], "figures/figure3/fig3e_stage_remodeling.pdf", "Figure 3E"),
        (figure_pdfs[5], "figures/figureS3/figs3c_cohort_robustness.pdf", "Figure S3 cohort robustness"),
        (figure_pdfs[6], "figures/figureS3/figs3_endpoint_deg_concordance.pdf", "Figure S3 endpoint sensitivity"),
        (figure_pdfs[7], "figures/figureS3/figs3_stage_remodeling_full.pdf", "Figure S3 stage and composition"),
        (args.figure_candidate / "source_tables/fig3a_five_cohort_summary.tsv", "figures/source_tables/fig3a_five_cohort_summary.tsv", "Figure 3A source"),
        (args.figure_candidate / "source_tables/fig3b_fibrosis_nas_grid.tsv", "figures/source_tables/fig3b_fibrosis_nas_grid.tsv", "Figure 3B source"),
        (args.figure_candidate / "source_tables/fig3c_cohort_deg_replication.tsv", "figures/source_tables/fig3c_cohort_deg_replication.tsv", "Figure 3C source"),
        (args.figure_candidate / "source_tables/fig3e_deg_counts.tsv", "figures/source_tables/fig3e_deg_counts.tsv", "Figure 3E source"),
        (args.figure_candidate / "source_tables/fig3e_prespecified_nmf_activity.tsv", "figures/source_tables/fig3e_prespecified_nmf_activity.tsv", "Figure 3E source"),
        (args.figure_candidate / "source_tables/figs3_composition_summary.tsv", "figures/source_tables/figs3_composition_summary.tsv", "Figure S3 source"),
        (args.figure_candidate / "source_tables/figs3_composition_summary_by_cohort.tsv", "figures/source_tables/figs3_composition_summary_by_cohort.tsv", "Figure S3 cohort-stratified source"),
        (args.figure_candidate / "source_tables/figs3c_cohort_robustness.tsv", "figures/source_tables/figs3c_cohort_robustness.tsv", "Figure S3 robustness source"),
        (args.figure_candidate / "source_tables/figs3_endpoint_deg_concordance_metrics.tsv", "figures/source_tables/figs3_endpoint_deg_concordance_metrics.tsv", "Figure S3 endpoint source"),
        (args.fig3g / "figure3/panels/fig3g_genetic_state_complementarity.pdf", "figures/figure3/fig3g_genetic_state_complementarity.pdf", "retained Figure 3G"),
        (args.fig3g / "figureS3/panels/figs3g_expression_matched_pairs.pdf", "figures/figureS3/figs3g_expression_matched_pairs.pdf", "retained Figure S3G"),
        (args.fig3g / "source_tables/fig3g_raw_treat_interface.tsv", "figures/source_tables/fig3g_raw_treat_interface.tsv", "retained Figure 3G source"),
        (args.fig3g / "source_tables/fig3g_matched_treat_summary.tsv", "figures/source_tables/fig3g_matched_treat_summary.tsv", "retained Figure S3G source"),
        (args.fig4_candidate / "source_tables/fig4e_bulk_projection.tsv", "artifacts/fig4_bulk_inputs/fig4e_bulk_projection.tsv", "retained Figure 4 bulk-side input"),
        (args.fig4_candidate / "source_tables/fig4c_all_117_programs.tsv", "artifacts/fig4_bulk_inputs/fig4c_all_117_programs.tsv", "frozen 117-program reference"),
        (args.fig4_candidate / "source_tables/figs4_panel_t5_composition_coupling_source.tsv", "artifacts/fig4_bulk_inputs/figs4_panel_t5_composition_coupling_source.tsv", "retained Figure S4 composition input"),
    ]:
        state = "retained_external_genetics_dependency" if "3G" in purpose or "Fig3G" in purpose or "Figure 4" in purpose else "adopted_bulk_only"
        add(source, relative, purpose, "as documented by source artifact", state)

    write_rows(output / "derived/composition_summary_by_cohort.tsv",
               ["summary_level", "axis", "contrast", "dataset", "group", "compartment", "mean", "SE", "n_participants", "estimand"], comp_summary)
    inventory.append({
        "release_path": "derived/composition_summary_by_cohort.tsv", "source_path": "derived from accepted composition plus stage membership",
        "size_bytes": (output / "derived/composition_summary_by_cohort.tsv").stat().st_size,
        "sha256": sha256(output / "derived/composition_summary_by_cohort.tsv"),
        "purpose": "cohort-stratified unadjusted composition summaries", "biological_unit": "one selected biopsy record", "state": "adopted_bulk_only",
    })
    write_rows(output / "validation/statistical_family_checks.tsv",
               ["analysis", "axis", "contrast", "arm", "n_genes", "n_reference", "n_comparison", "n_cohorts", "bh_family_size", "max_abs_BH_error", "intervals"],
               [dict(row, arm="") for row in stage_checks["stage_audits"]] +
               [dict(row, n_reference="", n_comparison="", n_cohorts="", bh_family_size="", intervals="") for row in f0_checks["f0_audits"]])
    checks = {
        "release_id": release_id,
        "status": "PASS",
        "pooled": pooled_checks,
        "stage_families": stage_checks["stage_families"],
        "f0_estimable_families": f0_checks["f0_families"],
        "composition_summary_rows": len(comp_summary),
        "figure3F": "blocked: exact corrected-COLOC signal-pair export is not packaged by this bulk-only release",
        "resource_release": "pending: genetic memberships, Gene Catalog, portal, and final synchronized Resource adoption are out of scope",
    }
    (output / "validation/validation_report.json").parent.mkdir(parents=True, exist_ok=True)
    (output / "validation/validation_report.json").write_text(json.dumps(checks, indent=2) + "\n")
    manifest = [{
        "release_id": release_id,
        "release_state": "adopted_bulk_only",
        "scope": "five_cohort_bulk_dependents",
        "pooled_samples": str(expected["n_samples"]),
        "controls": str(expected["n_control"]),
        "disease": str(expected["n_disease"]),
        "tested_genes": str(n_genes),
        "resource_synchronized_release": "pending",
    }]
    write_rows(output / "release_manifest.tsv", list(manifest[0]), manifest)
    (output / "BULK_RELEASE.json").write_text(json.dumps({
        **manifest[0],
        "description": "Immutable five-cohort fragment-count bulk-only adoption. This is not a final Resource release.",
        "figure3F": checks["figure3F"],
    }, indent=2) + "\n")
    write_rows(output / "artifact_inventory.tsv",
               ["release_path", "source_path", "size_bytes", "sha256", "purpose", "biological_unit", "state"], inventory)
    print(f"BULK_RELEASE_COMPLETE\t{output}")


if __name__ == "__main__":
    main()
