#!/usr/bin/env python3
"""Independently validate the descriptive cell-context association candidate."""

from __future__ import annotations

import argparse
import importlib.util
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from spatial_resource_lib import sha256_file, write_tsv


HERE = Path(__file__).resolve().parent
PRODUCER_PATH = HERE / "35_build_cell_context_attribution.py"
SPEC = importlib.util.spec_from_file_location("cell_context_producer", PRODUCER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load producer: {PRODUCER_PATH}")
producer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(producer)


def spearman_rho(left: np.ndarray, right: np.ndarray) -> float:
    x = np.asarray(left, dtype=float)
    y = np.asarray(right, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    rx = rankdata(x[keep], method="average")
    ry = rankdata(y[keep], method="average")
    return float(np.corrcoef(rx, ry)[0, 1])


def command_output(command: list[str]) -> str:
    return subprocess.run(command, check=True, capture_output=True, text=True).stdout


def pdf_pages(path: Path) -> int:
    for line in command_output(["pdfinfo", str(path)]).splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise RuntimeError(f"pdfinfo did not report pages for {path}")


def reject_type3_fonts(path: Path) -> None:
    pdffonts = shutil.which("pdffonts")
    if pdffonts:
        output = command_output([pdffonts, str(path)])
        found = "Type 3" in output or "Type3" in output
    else:
        found = re.search(rb"/Subtype\s*/Type3\b", path.read_bytes()) is not None
    if found:
        raise RuntimeError(f"Type 3 font found: {path}")


def validate(output: Path) -> list[dict[str, object]]:
    checks: list[dict[str, object]] = []

    def check(condition: bool, check_id: str, observed: object, expected: object, detail: str) -> None:
        if not condition:
            raise RuntimeError(f"{check_id}: {detail}; observed={observed!r}, expected={expected!r}")
        checks.append({
            "candidate_release_id": producer.CANDIDATE_RELEASE_ID,
            "check_id": check_id,
            "pass": "TRUE",
            "observed": observed,
            "expected": expected,
            "detail": detail,
        })

    build_complete = output / "BUILD_COMPLETE"
    check(build_complete.is_file(), "CTX-BUILD-01", build_complete.is_file(), True, "producer completed atomically")
    check(not (output / "READY").exists(), "CTX-IMMUTABLE-01", (output / "READY").exists(), False, "candidate has not already been sealed")

    correlation = pd.read_csv(output / "data/per_physical_unit_cell_context.tsv", sep="\t")
    unit = pd.read_csv(output / "data/per_reporting_unit_cell_context.tsv", sep="\t")
    summary = pd.read_csv(output / "data/dataset_cell_context_summary.tsv", sep="\t")
    concordance = pd.read_csv(output / "data/cross_source_direction_concordance.tsv", sep="\t")
    scores = pd.read_parquet(output / "data/per_spot_residual_program_scores.parquet")
    abundance = pd.read_parquet(output / "data/per_spot_cell2location_q05.parquet")
    heatmap = pd.read_csv(output / "data/figS_cell_context_association.tsv", sep="\t")
    audit = pd.read_csv(output / "data/external_source_audit.tsv", sep="\t")

    check(len(correlation) == 480, "CTX-FAMILY-01", len(correlation), 480, "15 physical units x 2 programs x 16 factors")
    check(len(unit) == 448, "CTX-FAMILY-02", len(unit), 448, "4 GSE donors plus 10 Vu arrays x 2 x 16")
    check(len(summary) == 64, "CTX-FAMILY-03", len(summary), 64, "2 datasets x 2 programs x 16 factors")
    check(len(concordance) == 32, "CTX-FAMILY-04", len(concordance), 32, "2 programs x 16 cross-source factor comparisons")
    check(len(heatmap) == 64, "CTX-FIGURE-01", len(heatmap), 64, "heatmap contains the full dataset/program/factor family")
    check(set(correlation["cell2location_factor"]) == set(producer.FACTOR_ORDER), "CTX-FACTORS-01", correlation["cell2location_factor"].nunique(), 16, "all shared q05 factors retained")
    ordered = summary.sort_values(["dataset", "program_uid", "factor_order"])["cell2location_factor"].tolist()
    expected_order = producer.FACTOR_ORDER * 4
    check(ordered == expected_order, "CTX-FACTORS-02", ordered[:16], producer.FACTOR_ORDER, "factor order is fixed biologically, not outcome-sorted")
    check(set(correlation["spot_level_inferential_pvalue_authorized"].astype(str).str.upper()) == {"FALSE"}, "CTX-PVALUE-01", "FALSE", "FALSE", "spot-level inference is prohibited")
    all_columns = set(correlation) | set(unit) | set(summary) | set(concordance) | set(heatmap)
    banned = sorted(column for column in all_columns if column.lower() in {"pvalue", "qvalue", "padj", "fdr", "standard_error"})
    check(not banned, "CTX-PVALUE-02", banned, [], "no inferential fields were emitted")
    vu = correlation[correlation["dataset"] == "Vu_et_al_2025"]
    check(set(vu["source_dependence"]) == {"source_dependent"}, "CTX-VU-01", set(vu["source_dependence"]), {"source_dependent"}, "Vu remains source-dependent")
    check(vu["technical_id"].nunique() == 10 and vu["source_individual_label"].nunique() == 10, "CTX-VU-02", (vu["technical_id"].nunique(), vu["source_individual_label"].nunique()), (10, 10), "Vu arrays were not upgraded to donors")
    gse_unit = unit[unit["dataset"] == "GSE192741"]
    check(gse_unit["summary_unit_id"].nunique() == 4, "CTX-GSE-01", gse_unit["summary_unit_id"].nunique(), 4, "GSE sections collapse to four donors")
    h35 = gse_unit[gse_unit["summary_unit_id"] == "H35"]
    check(set(h35["n_technical_units_collapsed"]) == {2}, "CTX-GSE-02", set(h35["n_technical_units_collapsed"]), {2}, "H35 two-section collapse is explicit")

    score_keys = ["candidate_release_id", "program_release_id", "dataset", "technical_id", "source_individual_label", "spot_id"]
    abundance_keys = score_keys[:-1] + ["spot_id"]
    joined = scores.merge(abundance, on=abundance_keys, how="left", validate="many_to_one")
    check(not joined[producer.FACTOR_COLUMNS].isna().any().any(), "CTX-REDERIVE-01", int(joined[producer.FACTOR_COLUMNS].isna().sum().sum()), 0, "every stored score joins exactly to q05 abundance")
    factor_map = dict(zip(producer.FACTOR_ORDER, producer.FACTOR_COLUMNS, strict=True))
    rederived = []
    for keys, group in joined.groupby(["dataset", "technical_id", "program_uid"], observed=True):
        for factor in producer.FACTOR_ORDER:
            rederived.append((*keys, factor, spearman_rho(group["residual_program_score_z"], group[factor_map[factor]])))
    rederived = pd.DataFrame(rederived, columns=["dataset", "technical_id", "program_uid", "cell2location_factor", "rho_rederived"])
    merged = correlation.merge(rederived, on=["dataset", "technical_id", "program_uid", "cell2location_factor"], validate="one_to_one")
    max_delta = float(np.max(np.abs(merged["spearman_rho"] - merged["rho_rederived"])))
    check(max_delta <= 1e-12, "CTX-REDERIVE-02", max_delta, "<=1e-12", "all 480 correlations rederive from stored spot values")

    expected_unit = producer.build_unit_summary(correlation)
    keys = ["dataset", "summary_unit_id", "program_uid", "cell2location_factor"]
    unit_join = unit.merge(expected_unit[keys + ["spearman_rho"]], on=keys, suffixes=("", "_expected"), validate="one_to_one")
    max_unit_delta = float(np.max(np.abs(unit_join["spearman_rho"] - unit_join["spearman_rho_expected"])))
    check(max_unit_delta <= 1e-12, "CTX-REDERIVE-03", max_unit_delta, "<=1e-12", "donor/array summaries rederive")
    expected_summary = producer.build_dataset_summary(unit)
    keys = ["dataset", "program_uid", "cell2location_factor"]
    summary_join = summary.merge(expected_summary[keys + ["median_spearman_rho"]], on=keys, suffixes=("", "_expected"), validate="one_to_one")
    max_summary_delta = float(np.max(np.abs(summary_join["median_spearman_rho"] - summary_join["median_spearman_rho_expected"])))
    check(max_summary_delta <= 1e-12, "CTX-REDERIVE-04", max_summary_delta, "<=1e-12", "dataset medians rederive")

    heatmap_join = heatmap.merge(summary[keys + ["median_spearman_rho"]], on=keys, suffixes=("", "_expected"), validate="one_to_one")
    max_heatmap_delta = float(np.max(np.abs(heatmap_join["median_spearman_rho"] - heatmap_join["median_spearman_rho_expected"])))
    check(max_heatmap_delta <= 1e-12, "CTX-FIGURE-02", max_heatmap_delta, "<=1e-12", "heatmap source exactly matches dataset summary")
    panel = output / "panels/figS_cell_context_association.pdf"
    check(panel.is_file() and panel.stat().st_size > 0, "CTX-FIGURE-03", panel.stat().st_size if panel.exists() else 0, ">0", "supplementary panel exists")
    check(pdf_pages(panel) == 1, "CTX-FIGURE-04", pdf_pages(panel), 1, "supplementary PDF has one page")
    reject_type3_fonts(panel)
    check(True, "CTX-FIGURE-05", "no Type 3 fonts", "no Type 3 fonts", "PDF remains editable")

    check(audit["dataset"].tolist() == ["Govaere2026_CosMx", "Yakubovsky2026"], "CTX-EXTERNAL-01", audit["dataset"].tolist(), ["Govaere2026_CosMx", "Yakubovsky2026"], "both narrow external sources audited")
    check(set(audit["additional_analysis_decision"]) == {"no_new_analysis_or_panel"}, "CTX-EXTERNAL-02", set(audit["additional_analysis_decision"]), {"no_new_analysis_or_panel"}, "no source-owned mechanism was broadened")
    cosmx = audit[audit["dataset"] == "Govaere2026_CosMx"].iloc[0]
    check(int(cosmx["eligible_biological_or_physical_units"]) == 3 and int(cosmx["audit_value_1"]) == 3 and int(cosmx["audit_value_2"]) == 3, "CTX-EXTERNAL-03", (cosmx["eligible_biological_or_physical_units"], cosmx["audit_value_1"], cosmx["audit_value_2"]), (3, 3, 3), "CosMx preserves 3/3 MASH-array direction concordance")
    yak = audit[audit["dataset"] == "Yakubovsky2026"].iloc[0]
    check(int(yak["eligible_biological_or_physical_units"]) == 3 and int(yak["audit_value_1"]) == 2 and int(yak["audit_value_2"]) == 0, "CTX-EXTERNAL-04", (yak["eligible_biological_or_physical_units"], yak["audit_value_1"], yak["audit_value_2"]), (3, 2, 0), "Yakubovsky remains a three-donor indeterminate reference")

    superseded = output.parent / "spatial-cell-context-attribution-candidate-2026-08-11"
    if superseded.is_dir():
        unchanged_tables = [
            "per_physical_unit_cell_context.tsv",
            "per_reporting_unit_cell_context.tsv",
            "dataset_cell_context_summary.tsv",
            "cross_source_direction_concordance.tsv",
            "figS_cell_context_association.tsv",
            "external_source_audit.tsv",
        ]
        for filename in unchanged_tables:
            old = pd.read_csv(superseded / "data" / filename, sep="\t", dtype=str, keep_default_na=False)
            new = pd.read_csv(output / "data" / filename, sep="\t", dtype=str, keep_default_na=False)
            excluded = {"candidate_release_id", "interpretation"}
            columns = [column for column in old.columns if column not in excluded]
            same = columns == [column for column in new.columns if column not in excluded] and old[columns].equals(new[columns])
            check(same, f"CTX-VISUAL-R2-{filename}", "identical" if same else "drift", "identical", f"{filename} scientific values are unchanged from the superseded visual candidate")
    return checks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_root.resolve()
    try:
        checks = validate(output)
        write_tsv(output / "validation_report.tsv", tuple(checks[0]), checks)
        panel = output / "panels/figS_cell_context_association.pdf"
        source = output / "data/figS_cell_context_association.tsv"
        write_tsv(
            output / "panel_manifest.tsv",
            ("candidate_release_id", "panel", "relative_path", "bytes", "sha256", "pdf_pages", "status"),
            [{
                "candidate_release_id": producer.CANDIDATE_RELEASE_ID,
                "panel": "S-cell-context",
                "relative_path": panel.relative_to(output).as_posix(),
                "bytes": panel.stat().st_size,
                "sha256": sha256_file(panel),
                "pdf_pages": 1,
                "status": "candidate_not_promoted",
            }],
        )
        write_tsv(
            output / "source_table_manifest.tsv",
            ("candidate_release_id", "panel", "relative_path", "bytes", "sha256"),
            [{
                "candidate_release_id": producer.CANDIDATE_RELEASE_ID,
                "panel": "S-cell-context",
                "relative_path": source.relative_to(output).as_posix(),
                "bytes": source.stat().st_size,
                "sha256": sha256_file(source),
            }],
        )
        write_tsv(
            output / "READY",
            (
                "candidate_release_id", "program_release_id", "resource_release_id", "status", "n_checks",
                "n_failed_checks", "validation_report_sha256", "panel_manifest_sha256",
                "source_table_manifest_sha256", "canonical_output_written",
            ),
            [{
                "candidate_release_id": producer.CANDIDATE_RELEASE_ID,
                "program_release_id": producer.PROGRAM_RELEASE_ID,
                "resource_release_id": producer.RESOURCE_RELEASE_ID,
                "status": "validated_descriptive_cell_context_candidate_awaiting_adjudication",
                "n_checks": len(checks),
                "n_failed_checks": 0,
                "validation_report_sha256": sha256_file(output / "validation_report.tsv"),
                "panel_manifest_sha256": sha256_file(output / "panel_manifest.tsv"),
                "source_table_manifest_sha256": sha256_file(output / "source_table_manifest.tsv"),
                "canonical_output_written": "FALSE",
            }],
        )
        print(f"PASS {len(checks)} cell-context checks")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
