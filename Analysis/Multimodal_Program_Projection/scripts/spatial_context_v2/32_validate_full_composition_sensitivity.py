#!/usr/bin/env python3
"""Independently validate the all-cell-type Visium sensitivity candidate."""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from pathlib import Path

import numpy as np

from visium_rerun_lib import (
    N_NULL,
    N_SENSITIVITY_NULL,
    SEED,
    build_paths,
    read_tsv,
    sha256_file,
    universe_hash,
    utc_now,
    v2_universe,
    write_tsv,
)
RELEASE_ID = "spatial-full-composition-sensitivity-candidate-2026-08-11-r2"
UPSTREAM_READY_SHA256 = "fcff53888cd8a29adc817752be4f3798e7fb3bcd7903f02622a90af729089898"
EXPECTED_FACTOR_COUNT = 16
PROGRAM_RESULTS = "full_composition_program_results.tsv"
UNIT_RESULTS = "full_composition_unit_results.tsv"
NULL_SUMMARY = "full_composition_null_summary.tsv"
MATCHING_AUDIT = "full_composition_matching_audit.tsv"
GRAPH_AUDIT = "full_composition_graph_audit.tsv"


def default_output_root(project_root: Path) -> Path:
    return project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID


def upstream_root(paths) -> Path:
    return paths.native_root / "v2_candidate"


class Checks:
    def __init__(self) -> None:
        self.rows: list[dict[str, str]] = []

    def require(self, condition: bool, check_id: str, detail: str) -> None:
        self.rows.append(
            {"check_id": check_id, "status": "PASS" if condition else "FAIL", "detail": detail}
        )

    @property
    def passed(self) -> bool:
        return bool(self.rows) and all(row["status"] == "PASS" for row in self.rows)


def parse_bool(value) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def independent_bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def import_producer(script_dir: Path):
    """Import the numeric-prefixed producer without inventing a duplicate library."""
    import importlib.util

    path = script_dir / "31_run_full_composition_sensitivity.py"
    spec = importlib.util.spec_from_file_location("full_composition_producer", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import producer: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def compare_primary_passthrough(
    observed,
    primary,
    release_column: bool,
    checks: Checks,
    label: str,
    relaxed_float_columns: dict[str, float] | None = None,
) -> None:
    left = observed.drop(columns=["release_id"]) if release_column else observed.copy()
    right = primary.copy()
    checks.require(list(left.columns) == list(right.columns), f"{label}_schema", "upstream columns exact")
    if list(left.columns) != list(right.columns):
        return
    sort_columns = [column for column in ("dataset", "program_id", "sample_id", "gene_symbol") if column in left]
    if sort_columns:
        left = left.sort_values(sort_columns).reset_index(drop=True)
        right = right.sort_values(sort_columns).reset_index(drop=True)
    relaxed_float_columns = relaxed_float_columns or {}
    mismatches = []
    maxima = []
    for column in left.columns:
        if column in relaxed_float_columns:
            left_values = left[column].to_numpy(float)
            right_values = right[column].to_numpy(float)
            maximum = float(np.nanmax(np.abs(left_values - right_values)))
            maxima.append(f"{column}={maximum:.3g}")
            if not np.allclose(
                left_values,
                right_values,
                rtol=0,
                atol=relaxed_float_columns[column],
                equal_nan=True,
            ):
                mismatches.append(column)
        elif np.issubdtype(left[column].dtype, np.number) and np.issubdtype(right[column].dtype, np.number):
            if not np.allclose(
                left[column].to_numpy(float),
                right[column].to_numpy(float),
                rtol=1e-12,
                atol=1e-14,
                equal_nan=True,
            ):
                mismatches.append(column)
        elif not np.array_equal(
            left[column].fillna("").astype(str).to_numpy(),
            right[column].fillna("").astype(str).to_numpy(),
        ):
            mismatches.append(column)
    same = not mismatches
    detail = (
        f"rows={len(left)}; no key/bin/control-count drift; "
        + ("; ".join(maxima) if maxima else "strict numeric/string equality")
        if same
        else "mismatched columns=" + ",".join(mismatches)
    )
    checks.require(same, f"{label}_unchanged", detail)


def validate(project_root: Path | None, output_root: Path | None):
    import pandas as pd

    paths = build_paths(project_root)
    producer = import_producer(Path(__file__).resolve().parent)
    if any(
        getattr(producer, name) != value
        for name, value in (
            ("RELEASE_ID", RELEASE_ID),
            ("UPSTREAM_READY_SHA256", UPSTREAM_READY_SHA256),
            ("EXPECTED_FACTOR_COUNT", EXPECTED_FACTOR_COUNT),
        )
    ):
        raise RuntimeError("producer/validator release constants disagree")
    output = (output_root or default_output_root(paths.project_root)).resolve()
    if not output.is_dir():
        raise RuntimeError(f"candidate output is missing: {output}")
    if (output / "READY").exists() or (output / "validation_report.tsv").exists():
        raise RuntimeError("refusing to overwrite existing validation/READY artifacts")
    checks = Checks()

    candidates = (paths.project_root / "Analysis/Multimodal_Program_Projection/candidates").resolve()
    checks.require(output.parent == candidates, "isolated_candidate_root", str(output))
    checks.require(output.name == RELEASE_ID, "candidate_release_id", RELEASE_ID)

    manifest = output / "output_manifest.tsv"
    manifest_rows = read_tsv(manifest)
    checks.require(bool(manifest_rows), "output_manifest_nonempty", f"rows={len(manifest_rows)}")
    for row in manifest_rows:
        path = output / row["relative_path"]
        ok = (
            path.is_file()
            and path.stat().st_size == int(row["bytes"])
            and sha256_file(path) == row["sha256"]
        )
        checks.require(ok, f"manifest_{row['relative_path']}", "size and SHA256 rederived")
    complete = read_tsv(output / "COMPUTE_COMPLETE")
    checks.require(
        len(complete) == 1
        and complete[0].get("release_id") == RELEASE_ID
        and complete[0].get("status") == "compute_complete_pending_independent_validation"
        and complete[0].get("output_manifest_sha256") == sha256_file(manifest),
        "compute_complete_seal",
        "compute marker links release and output manifest",
    )

    for row in read_tsv(output / "source_manifest.tsv"):
        path = paths.project_root / row["relative_path"]
        checks.require(
            path.is_file()
            and path.stat().st_size == int(row["bytes"])
            and sha256_file(path) == row["sha256"],
            f"source_{row['source_role']}",
            "source size and SHA256 unchanged",
        )
    upstream = upstream_root(paths)
    checks.require(
        sha256_file(upstream / "READY") == UPSTREAM_READY_SHA256,
        "upstream_ready_gate",
        UPSTREAM_READY_SHA256,
    )

    contract = {row["parameter"]: row["value"] for row in read_tsv(output / "analysis_contract.tsv")}
    expected_universe = v2_universe(paths)
    checks.require(contract.get("release_id") == RELEASE_ID, "contract_release", RELEASE_ID)
    checks.require(
        contract.get("program_family_sha256") == universe_hash(expected_universe)
        and len(expected_universe) == 2,
        "frozen_two_program_family",
        universe_hash(expected_universe),
    )
    checks.require(
        contract.get("cell2location_factor_count") == str(EXPECTED_FACTOR_COUNT),
        "all_factor_count_contract",
        str(EXPECTED_FACTOR_COUNT),
    )
    checks.require(contract.get("n_null") == str(N_NULL), "null_count_contract", str(N_NULL))
    checks.require(
        contract.get("n_sensitivity_null") == str(N_SENSITIVITY_NULL),
        "sensitivity_null_contract",
        str(N_SENSITIVITY_NULL),
    )
    checks.require(contract.get("canonical_write_allowed") == "FALSE", "canonical_write_forbidden", "FALSE")

    designs, factor_names = producer.load_full_designs(paths)
    checks.require(
        contract.get("cell2location_factor_order") == "|".join(factor_names),
        "factor_order_rederived",
        f"n={len(factor_names)}",
    )
    design_audit = pd.read_csv(output / "full_composition_design_audit.tsv", sep="\t", dtype=str)
    checks.require(len(design_audit) == 2, "design_audit_rows", "one per dataset")
    for row in design_audit.itertuples(index=False):
        item = designs[row.dataset]
        checks.require(
            int(row.n_spots) == len(item["obs"])
            and int(row.n_cell2location_factors) == EXPECTED_FACTOR_COUNT
            and int(row.n_design_columns) == EXPECTED_FACTOR_COUNT + 3
            and int(row.design_rank) == item["rank"]
            and math.isclose(float(row.condition_number), item["condition_number"], rel_tol=1e-12),
            f"design_rederived_{row.dataset}",
            f"shape={item['design'].shape}; rank={item['rank']}; condition={item['condition_number']:.4g}",
        )
        checks.require(
            row.unit_gate == ("pass_donor_first" if row.dataset == "GSE192741" else "source_dependent"),
            f"unit_gate_{row.dataset}",
            row.unit_gate,
        )

    results = pd.read_csv(output / PROGRAM_RESULTS, sep="\t")
    expected_ids = {row["program_id"] for row in expected_universe}
    expected_keys = {(dataset, uid) for dataset in paths.datasets for uid in expected_ids}
    observed_keys = set(zip(results["dataset"], results["program_id"], strict=True))
    checks.require(len(results) == 4 and observed_keys == expected_keys, "complete_four_row_family", str(sorted(observed_keys)))
    checks.require(not results.duplicated(["dataset", "program_id"]).any(), "unique_result_keys", "no duplicates")
    checks.require(results["testable"].map(parse_bool).all(), "all_four_testable", "two programs x two datasets")
    checks.require(set(results["n_null"].astype(int)) == {N_NULL}, "result_null_count", str(N_NULL))
    checks.require(
        set(results.loc[results["dataset"] == "Vu_et_al_2025", "evidence_state"]) == {"source_dependent"}
        and set(results.loc[results["dataset"] == "Vu_et_al_2025", "biological_unit"]) == {"unresolved_physical_array"},
        "vu_source_dependence",
        "no donor upgrade",
    )

    for dataset, part in results.groupby("dataset", sort=False):
        expected_q = independent_bh(part["residual_pvalue"].to_numpy(float))
        checks.require(
            np.allclose(expected_q, part["residual_padj"].to_numpy(float), rtol=1e-12, atol=1e-14),
            f"bh_complete_family_{dataset}",
            "two-program family independently rederived",
        )
        centered = part["residual_moran_i"] - part["residual_null_mean"]
        equal = part["equal_residual_moran_i"] - part["equal_residual_null_mean"]
        leave = part["leave_top_residual_moran_i"] - part["leave_top_residual_null_mean"]
        expected_sign = (np.sign(centered) == np.sign(equal)) & (np.sign(centered) == np.sign(leave))
        expected_robust = (part["residual_padj"] < 0.05) & expected_sign
        checks.require(
            np.array_equal(expected_sign, part["sensitivity_sign_agree"].map(parse_bool)),
            f"sensitivity_sign_rule_{dataset}",
            "equal and leave-top centered directions",
        )
        checks.require(
            np.array_equal(expected_robust, part["robust"].map(parse_bool)),
            f"support_rule_{dataset}",
            "q<0.05 plus both sensitivity directions",
        )

    draws = pd.read_csv(output / "full_composition_null_draws.tsv.gz", sep="\t")
    checks.require(
        len(draws) == 4 * N_NULL
        and not draws.duplicated(["dataset", "program_id", "draw_id"]).any(),
        "complete_null_draws",
        f"rows={len(draws)}",
    )
    null_summary = pd.read_csv(output / NULL_SUMMARY, sep="\t")
    checks.require(len(null_summary) == 4, "null_summary_rows", "one residual null per result")
    for result in results.itertuples(index=False):
        part = draws[(draws["dataset"] == result.dataset) & (draws["program_id"] == result.program_id)]
        values = part.sort_values("draw_id")["residual_moran_i"].to_numpy(float)
        summary = null_summary[
            (null_summary["dataset"] == result.dataset)
            & (null_summary["program_id"] == result.program_id)
        ].iloc[0]
        pvalue = (1 + np.sum(values >= result.residual_moran_i)) / (1 + len(values))
        quantiles = np.quantile(values, [0.001, 0.01, 0.05, 0.50, 0.95, 0.99, 0.999])
        checks.require(
            len(values) == N_NULL
            and math.isclose(values.mean(), result.residual_null_mean, rel_tol=1e-12, abs_tol=1e-14)
            and math.isclose(values.std(ddof=1), result.residual_null_sd, rel_tol=1e-12, abs_tol=1e-14)
            and math.isclose(pvalue, result.residual_pvalue, rel_tol=0, abs_tol=1e-15),
            f"null_rederived_{result.dataset}_{result.program_id}",
            f"mean={values.mean():.6g}; sd={values.std(ddof=1):.6g}; p={pvalue:.6g}",
        )
        checks.require(
            np.allclose(
                quantiles,
                summary[["q001", "q010", "q050", "q500", "q950", "q990", "q999"]].to_numpy(float),
                rtol=1e-12,
                atol=1e-14,
            ),
            f"null_quantiles_{result.dataset}_{result.program_id}",
            "seven quantiles rederived",
        )

    reuse = pd.read_csv(output / "matched_set_reuse_audit.tsv", sep="\t", dtype=str)
    checks.require(len(reuse) == 4 and reuse["identical"].map(parse_bool).all(), "matched_sets_reused", "4/4 hashes identical")
    primary_results = pd.read_csv(upstream / "spatial_program_results.tsv", sep="\t")
    primary_hashes = primary_results.set_index(["dataset", "program_id"])["matched_set_sha256"].to_dict()
    checks.require(
        all(
            row.sensitivity_matched_set_sha256 == row.original_matched_set_sha256
            == primary_hashes[(row.dataset, row.program_id)]
            for row in reuse.itertuples(index=False)
        ),
        "matched_hashes_against_primary",
        "same matching inputs/seeds; residual design alone changed",
    )

    primary_graph = pd.read_csv(upstream / "spatial_graph_audit.tsv", sep="\t")
    observed_graph = pd.read_csv(output / GRAPH_AUDIT, sep="\t")
    compare_primary_passthrough(observed_graph, primary_graph, True, checks, "graph")
    primary_matching = pd.read_csv(upstream / "spatial_matching_audit.tsv", sep="\t")
    observed_matching = pd.read_csv(output / MATCHING_AUDIT, sep="\t")
    compare_primary_passthrough(
        observed_matching,
        primary_matching,
        True,
        checks,
        "matching",
        relaxed_float_columns={
            # Sparse matrix/BLAS threading can move recomputed Pearson
            # correlations by ~1e-8 without changing rank bins or controls.
            # The explicit 2e-8 absolute bound is paired with exact
            # bin/count/matched-set hashes, so it cannot bless control drift.
            "target_lineage_correlation": 2e-8,
            "mean_control_lineage_correlation": 2e-8,
        },
    )

    units = pd.read_csv(output / UNIT_RESULTS, sep="\t")
    checks.require(len(units) == 30, "unit_result_rows", "2 programs x (5 sections + 10 arrays)")
    checks.require(
        units.loc[units["dataset"] == "GSE192741", "sample_id"].nunique() == 5
        and units.loc[units["dataset"] == "GSE192741", "individual"].nunique() == 4,
        "gse_units",
        "five sections collapse to four donors",
    )
    checks.require(
        units.loc[units["dataset"] == "Vu_et_al_2025", "sample_id"].nunique() == 10
        and not units.loc[units["dataset"] == "Vu_et_al_2025", "population_inference_authorized"].map(parse_bool).any(),
        "vu_units",
        "ten arrays; population inference forbidden",
    )
    for result in results.itertuples(index=False):
        part = units[(units["dataset"] == result.dataset) & (units["program_id"] == result.program_id)]
        collapsed = part.groupby("individual", observed=False)["residual_moran_i"].mean().mean()
        checks.require(
            math.isclose(collapsed, result.residual_moran_i, rel_tol=1e-12, abs_tol=1e-14),
            f"unit_collapse_{result.dataset}_{result.program_id}",
            f"rederived={collapsed:.8g}",
        )

    comparison = pd.read_csv(output / "comparison_to_primary.tsv", sep="\t")
    checks.require(len(comparison) == 4, "primary_comparison_rows", "four paired effects")
    checks.require(comparison["matched_set_hash_identical"].map(parse_bool).all(), "comparison_hash_identity", "all paired")
    for row in comparison.itertuples(index=False):
        current = results[(results["dataset"] == row.dataset) & (results["program_id"] == row.program_id)].iloc[0]
        primary = primary_results[(primary_results["dataset"] == row.dataset) & (primary_results["program_id"] == row.program_id)].iloc[0]
        expected_primary = primary.residual_moran_i - primary.residual_null_mean
        expected_current = current.residual_moran_i - current.residual_null_mean
        checks.require(
            math.isclose(row.primary_centered_moran_i, expected_primary, rel_tol=1e-12)
            and math.isclose(row.full_composition_centered_moran_i, expected_current, rel_tol=1e-12)
            and math.isclose(row.centered_change, expected_current - expected_primary, rel_tol=1e-12, abs_tol=1e-14),
            f"comparison_rederived_{row.dataset}_{row.program_id}",
            "primary and sensitivity centered effects",
        )

    report = output / "validation_report.tsv"
    write_tsv(report, ("check_id", "status", "detail"), checks.rows)
    if not checks.passed:
        failures = [row["check_id"] for row in checks.rows if row["status"] == "FAIL"]
        raise RuntimeError("validation failed: " + ", ".join(failures))
    write_tsv(
        output / "READY",
        (
            "release_id", "status", "validation_report_sha256", "output_manifest_sha256",
            "upstream_ready_sha256", "canonical_written", "validated_utc",
        ),
        [{
            "release_id": RELEASE_ID,
            "status": "pass_full_composition_sensitivity",
            "validation_report_sha256": sha256_file(report),
            "output_manifest_sha256": sha256_file(manifest),
            "upstream_ready_sha256": UPSTREAM_READY_SHA256,
            "canonical_written": "FALSE",
            "validated_utc": utc_now(),
        }],
    )
    print(f"full-composition sensitivity validated: {len(checks.rows)} checks passed", flush=True)
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        validate(args.project_root, args.output_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
