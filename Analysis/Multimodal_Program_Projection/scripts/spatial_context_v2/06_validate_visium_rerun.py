#!/usr/bin/env python3
"""Independently validate SP-INT-03 v1 regression and frozen-v2 outputs."""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from pathlib import Path

from visium_rerun_lib import (
    CORE_OUTPUTS,
    N_NULL,
    RELEASE_ID,
    V1_ENGINE_SHA256,
    build_paths,
    read_tsv,
    sha256_file,
    universe_hash,
    utc_now,
    v1_universe,
    v2_universe,
    verify_hotspot_ready,
    verify_v1_anchors,
    write_tsv,
)


class Checks:
    def __init__(self) -> None:
        self.rows: list[dict[str, str]] = []

    def require(self, condition: bool, check_id: str, detail: str) -> None:
        self.rows.append({"check_id": check_id, "status": "PASS" if condition else "FAIL", "detail": detail})

    @property
    def passed(self) -> bool:
        return bool(self.rows) and all(row["status"] == "PASS" for row in self.rows)


def parse_bool(value) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def independent_bh(values):
    import numpy as np

    values = np.asarray(values, dtype=float)
    out = np.full(len(values), np.nan, dtype=float)
    valid = np.isfinite(values)
    p = values[valid]
    if not len(p):
        return out
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    out[np.flatnonzero(valid)] = restored
    return out


def verify_output_manifest(output: Path, checks: Checks) -> None:
    manifest = output / "output_manifest.tsv"
    rows = read_tsv(manifest)
    checks.require(bool(rows), "output_manifest_nonempty", f"rows={len(rows)}")
    for row in rows:
        path = output / row["relative_path"]
        ok = path.is_file() and path.stat().st_size == int(row["bytes"]) and sha256_file(path) == row["sha256"]
        checks.require(ok, f"manifest_{row['relative_path']}", "size and SHA256 rederived")
    complete = read_tsv(output / "COMPUTE_COMPLETE")
    checks.require(
        len(complete) == 1
        and complete[0].get("status") == "compute_complete_pending_independent_validation"
        and complete[0].get("output_manifest_sha256") == sha256_file(manifest),
        "compute_complete_seal",
        "compute marker links to output manifest",
    )


def compare_frames(observed, expected, name: str, checks: Checks) -> None:
    import numpy as np

    checks.require(list(observed.columns) == list(expected.columns), f"{name}_schema", "column order exact")
    checks.require(observed.shape == expected.shape, f"{name}_shape", f"observed={observed.shape}; expected={expected.shape}")
    if list(observed.columns) != list(expected.columns) or observed.shape != expected.shape:
        return
    mismatch = 0
    max_abs = 0.0
    for column in observed.columns:
        left = observed[column]
        right = expected[column]
        left_num = np.array([float(x) if str(x) not in {"", "nan", "NA", "NaN"} else np.nan for x in left], dtype=float) if all(
            str(x) in {"", "nan", "NA", "NaN"} or _is_float(str(x)) for x in left
        ) and all(str(x) in {"", "nan", "NA", "NaN"} or _is_float(str(x)) for x in right) else None
        if left_num is not None:
            right_num = np.array([float(x) if str(x) not in {"", "nan", "NA", "NaN"} else np.nan for x in right], dtype=float)
            valid = np.isfinite(left_num) & np.isfinite(right_num)
            if np.any(valid):
                max_abs = max(max_abs, float(np.max(np.abs(left_num[valid] - right_num[valid]))))
            mismatch += int(np.sum(~np.isclose(left_num, right_num, rtol=1e-10, atol=1e-12, equal_nan=True)))
        else:
            mismatch += sum(str(a) != str(b) for a, b in zip(left.fillna(""), right.fillna(""), strict=True))
    checks.require(mismatch == 0, f"{name}_semantic_regression", f"mismatches={mismatch}; max_abs_numeric={max_abs:.3g}")


def _is_float(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def validate(stage: str, paths) -> tuple[Checks, Path]:
    import numpy as np
    import pandas as pd

    checks = Checks()
    verify_v1_anchors(paths)
    if stage == "v2":
        verify_hotspot_ready(paths)
    output = paths.native_root / ("v1_regression" if stage == "v1" else "v2_candidate")
    if not output.is_dir():
        raise RuntimeError(f"missing {stage} candidate output: {output}")
    if (output / "READY").exists() or (output / "validation_report.tsv").exists():
        raise RuntimeError(f"refusing to overwrite existing validation artifacts in {output}")
    verify_output_manifest(output, checks)

    source_rows = read_tsv(output / "source_manifest.tsv")
    for row in source_rows:
        path = paths.project_root / row["relative_path"]
        checks.require(
            path.is_file() and path.stat().st_size == int(row["bytes"]) and sha256_file(path) == row["sha256"],
            f"source_{row['source_role']}",
            "input size and SHA256 unchanged since compute start",
        )
    execution = {row["parameter"]: row["value"] for row in read_tsv(output / "execution_manifest.tsv")}
    checks.require(execution.get("registry_version") == stage, "execution_registry_version", stage)
    checks.require(execution.get("v1_engine_sha256") == V1_ENGINE_SHA256, "engine_hash", V1_ENGINE_SHA256)
    checks.require(execution.get("n_null") == str(N_NULL), "primary_null_count_spec", str(N_NULL))
    checks.require(execution.get("n_sensitivity_null") == "999", "sensitivity_null_count_spec", "999")
    checks.require(execution.get("canonical_write_allowed") == "FALSE", "isolated_write_contract", "canonical write forbidden")

    universe = read_tsv(output / "tested_universe.tsv")
    expected_universe = v1_universe(paths) if stage == "v1" else v2_universe(paths)
    checks.require(
        universe_hash(universe) == universe_hash(expected_universe),
        "tested_family_hash",
        f"programs={len(universe)}; hash={universe_hash(universe)}",
    )
    checks.require(
        execution.get("tested_family_sha256") == universe_hash(expected_universe),
        "execution_family_link",
        "execution manifest links complete family",
    )

    frames = {
        name: pd.read_csv(output / name, sep="\t", dtype=str, keep_default_na=False)
        for name in CORE_OUTPUTS
    }
    results_text = frames["spatial_program_results.tsv"]
    expected_ids = {str(row["program_id"]) for row in expected_universe}
    expected_keys = {(dataset, pid) for dataset in paths.datasets for pid in expected_ids}
    observed_keys = set(zip(results_text["dataset"], results_text["program_id"], strict=True))
    checks.require(observed_keys == expected_keys, "complete_program_dataset_grid", f"rows={len(observed_keys)}")
    checks.require(not results_text.duplicated(["dataset", "program_id"]).any(), "unique_program_dataset_key", "no duplicates")
    checks.require(set(results_text["n_null"]) == {str(N_NULL)}, "result_null_count", str(N_NULL))

    samples = read_tsv(output / "native_sample_manifest.tsv")
    checks.require(len(samples) == 15, "native_sample_manifest_rows", "5 GSE192741 sections + 10 Vu arrays")
    sample_by_dataset = {}
    for dataset in paths.datasets:
        rows = [row for row in samples if row["dataset"] == dataset]
        sample_by_dataset[dataset] = rows
        checks.require(
            all(row["include_primary"] == "TRUE" and row["gate_state"] == "pass" for row in rows),
            f"sample_gate_{dataset}",
            "all frozen source units retained",
        )
    checks.require(
        len({row["technical_id"] for row in sample_by_dataset["GSE192741"]}) == 5
        and len({row["biological_id"] for row in sample_by_dataset["GSE192741"]}) == 4,
        "gse_sample_design",
        "5 sections -> 4 donors",
    )
    checks.require(
        len({row["technical_id"] for row in sample_by_dataset["Vu_et_al_2025"]}) == 10
        and len({row["biological_id"] for row in sample_by_dataset["Vu_et_al_2025"]}) == 10
        and all(row["biological_unit"] == "physical_array" for row in sample_by_dataset["Vu_et_al_2025"]),
        "vu_sample_design",
        "10 arrays remain 10 source-dependent reporting units",
    )
    design_rows = read_tsv(output / "native_design_audit.tsv")
    checks.require(len(design_rows) == 2, "native_design_rows", "one row per dataset")
    for row in design_rows:
        rows = sample_by_dataset[row["dataset"]]
        biological_ids = sorted({item["biological_id"] for item in rows})
        technical_ids = sorted({item["technical_id"] for item in rows})
        biological_sha = hashlib.sha256(("\n".join(biological_ids) + "\n").encode()).hexdigest()
        technical_sha = hashlib.sha256(("\n".join(technical_ids) + "\n").encode()).hexdigest()
        checks.require(
            int(row["n_biological"]) == len(biological_ids)
            and int(row["n_technical"]) == len(technical_ids)
            and row["biological_ids_sha256"] == biological_sha
            and row["technical_ids_sha256"] == technical_sha,
            f"design_rederived_{row['dataset']}",
            "counts and ID hashes rederived from sample manifest",
        )

    results = pd.read_csv(output / "spatial_program_results.tsv", sep="\t")
    for dataset, part in results.groupby("dataset", sort=False):
        for p_col, q_col in (
            ("residual_pvalue", "residual_padj"),
            ("raw_pvalue", "raw_padj"),
            ("zonation_pvalue", "zonation_padj"),
        ):
            expected_q = independent_bh(part[p_col].to_numpy(float))
            observed_q = part[q_col].to_numpy(float)
            checks.require(
                np.allclose(expected_q, observed_q, rtol=1e-12, atol=1e-14, equal_nan=True),
                f"bh_{dataset}_{q_col}",
                f"complete within-dataset family n={len(part)}",
            )
        expected_robust = (
            part["testable"].map(parse_bool)
            & (part["residual_padj"] < 0.05)
            & part["sensitivity_sign_agree"].map(parse_bool)
        )
        checks.require(
            np.array_equal(expected_robust.to_numpy(bool), part["robust"].map(parse_bool).to_numpy(bool)),
            f"robust_rule_{dataset}",
            "testable & residual BH<0.05 & both sensitivity directions",
        )

    graphs = frames["spatial_graph_audit.tsv"]
    graph_counts = graphs.groupby("dataset")["sample_id"].nunique().to_dict()
    checks.require(graph_counts == {"GSE192741": 5, "Vu_et_al_2025": 10}, "physical_array_counts", str(graph_counts))
    checks.require((pd.to_numeric(graphs["n_tissue_islands"]) >= 1).all(), "connected_island_graphs", "every physical array has >=1 retained island")
    checks.require(
        set(zip(graphs["dataset"], graphs["sample_id"], strict=True))
        == {(row["dataset"], row["technical_id"]) for row in samples},
        "graph_sample_manifest_link",
        "one graph audit row per frozen physical array/section",
    )
    sections = frames["spatial_section_results.tsv"]
    gse = sections.loc[sections["dataset"] == "GSE192741"]
    vu = sections.loc[sections["dataset"] == "Vu_et_al_2025"]
    checks.require(gse["sample_id"].nunique() == 5 and gse["individual"].nunique() == 4, "gse_donor_collapse_design", "5 sections -> 4 donors")
    checks.require(vu["sample_id"].nunique() == 10 and vu["individual"].nunique() == 10, "vu_array_resolution", "10 arrays remain 10 reporting units")

    testability = read_tsv(output / "spatial_testability_audit.tsv")
    expected_membership = {str(row["program_id"]): str(row["membership_sha256"]) for row in expected_universe}
    checks.require(len(testability) == len(results), "testability_audit_rows", "one row per dataset/program")
    testability_by_key = {(row["dataset"], row["program_id"]): row for row in testability}
    checks.require(set(testability_by_key) == observed_keys, "testability_audit_keys", "complete program grid")
    testability_ok = True
    for row in results.itertuples(index=False):
        audit = testability_by_key[(row.dataset, row.program_id)]
        testability_ok &= (
            audit["membership_sha256"] == expected_membership[row.program_id]
            and int(audit["n_genes_measured"]) == int(row.n_measured)
            and math.isclose(float(audit["retained_l1_weight"]), float(row.retained_l1_weight), rel_tol=1e-12, abs_tol=1e-14)
            and parse_bool(audit["testable"]) == parse_bool(row.testable)
        )
    checks.require(testability_ok, "testability_audit_values", "membership, coverage, L1 mass, and gate rederived")

    testable = results.loc[results["testable"].map(parse_bool)]
    nulls = frames["spatial_null_summary.tsv"]
    null_keys = set(zip(nulls["dataset"], nulls["program_id"], nulls["statistic"], strict=True))
    expected_null_keys = {
        (row.dataset, row.program_id, statistic)
        for row in testable.itertuples(index=False)
        for statistic in ("raw_moran_i", "residual_moran_i", "zonation_moran_i")
    }
    checks.require(null_keys == expected_null_keys, "complete_null_summary_family", f"rows={len(null_keys)}")
    checks.require(set(nulls["n_null"]) == {str(N_NULL)}, "finite_primary_nulls", f"all summaries n={N_NULL}")
    result_hash = {(row.dataset, row.program_id): row.matched_set_sha256 for row in testable.itertuples(index=False)}
    checks.require(
        all(row["matched_set_sha256"] == result_hash[(row["dataset"], row["program_id"])] for row in nulls.to_dict("records")),
        "matched_set_hash_linkage",
        "program result and all three null summaries agree",
    )
    matching = frames["spatial_matching_audit.tsv"]
    checks.require(
        len(matching) == int(testable["n_measured"].sum()),
        "matching_audit_completeness",
        f"rows={len(matching)}; expected={int(testable['n_measured'].sum())}",
    )
    sensitivities = read_tsv(output / "spatial_sensitivity_audit.tsv")
    expected_sensitivity_keys = {
        (row.dataset, row.program_id, sensitivity_id)
        for row in testable.itertuples(index=False)
        for sensitivity_id in (
            "primary_original_l1_weight",
            "equal_weight",
            "leave_highest_weight_gene_out",
        )
    }
    observed_sensitivity_keys = {
        (row["dataset"], row["program_id"], row["sensitivity_id"]) for row in sensitivities
    }
    checks.require(observed_sensitivity_keys == expected_sensitivity_keys, "sensitivity_audit_keys", f"rows={len(sensitivities)}")
    sensitivity_by_key = {
        (row["dataset"], row["program_id"], row["sensitivity_id"]): row for row in sensitivities
    }
    sensitivity_ok = True
    for row in testable.itertuples(index=False):
        centers = {
            "primary_original_l1_weight": float(row.residual_moran_i - row.residual_null_mean),
            "equal_weight": float(row.equal_residual_moran_i - row.equal_residual_null_mean),
            "leave_highest_weight_gene_out": float(row.leave_top_residual_moran_i - row.leave_top_residual_null_mean),
        }
        primary = centers["primary_original_l1_weight"]
        for sensitivity_id, center in centers.items():
            audit = sensitivity_by_key[(row.dataset, row.program_id, sensitivity_id)]
            sign_agree = (center == 0 and primary == 0) or center * primary > 0
            sensitivity_ok &= math.isclose(float(audit["centered_residual_moran_i"]), center, rel_tol=1e-12, abs_tol=1e-14)
            sensitivity_ok &= parse_bool(audit["sign_agree_with_primary"]) == sign_agree
        expected_overall = parse_bool(sensitivity_by_key[(row.dataset, row.program_id, "equal_weight")]["sign_agree_with_primary"]) and parse_bool(
            sensitivity_by_key[(row.dataset, row.program_id, "leave_highest_weight_gene_out")]["sign_agree_with_primary"]
        )
        sensitivity_ok &= expected_overall == parse_bool(row.sensitivity_sign_agree)
    checks.require(sensitivity_ok, "sensitivity_audit_values", "all centered effects and signs independently rederived")

    canonical_graph = pd.read_csv(paths.v1_results / "spatial_graph_audit.tsv", sep="\t", dtype=str, keep_default_na=False)
    compare_frames(graphs, canonical_graph, "graph_engine", checks)
    if stage == "v1":
        for name in CORE_OUTPUTS:
            observed = frames[name]
            expected = pd.read_csv(paths.v1_results / name, sep="\t", dtype=str, keep_default_na=False)
            compare_frames(observed, expected, name.removesuffix(".tsv"), checks)
    else:
        ready_path = paths.native_root / "v1_regression/READY"
        checks.require(ready_path.is_file(), "v1_regression_ready_present", str(ready_path))
        if ready_path.is_file():
            ready = read_tsv(ready_path)
            checks.require(len(ready) == 1 and ready[0].get("status") == "pass_v1_regression", "v1_regression_ready_valid", "v2 is downstream of passing regression")
            checks.require(
                execution.get("v1_regression_ready_sha256") == sha256_file(ready_path),
                "v1_regression_ready_link",
                "execution manifest pins regression READY",
            )
        checks.require(len(expected_ids) == 2 and all(pid.startswith("hotspot_hepatocytes_") for pid in expected_ids), "frozen_v2_family", f"ids={sorted(expected_ids)}")
        reuse = read_tsv(paths.candidate_root / "reuse_eligibility.tsv")
        checks.require(len(reuse) == 4 and all(row["reuse_eligible"] == "FALSE" for row in reuse), "rerun_not_reuse", "2 programs x 2 datasets explicitly ineligible")
    return checks, output


def write_validation(checks: Checks, output: Path, stage: str) -> None:
    report = output / "validation_report.tsv"
    write_tsv(report, ("check_id", "status", "detail"), checks.rows)
    if not checks.passed:
        raise RuntimeError(f"SP-INT-03 {stage} validation failed; see {report}")
    core_aggregate = "".join(
        f"{name}\0{(output / name).stat().st_size}\0{sha256_file(output / name)}\n" for name in CORE_OUTPUTS
    )
    import hashlib

    write_tsv(
        output / "READY",
        (
            "release_id",
            "registry_version",
            "status",
            "v1_engine_sha256",
            "validation_report_sha256",
            "core_outputs_aggregate_sha256",
            "validated_utc",
        ),
        [
            {
                "release_id": RELEASE_ID,
                "registry_version": stage,
                "status": "pass_v1_regression" if stage == "v1" else "pass_v2_candidate",
                "v1_engine_sha256": V1_ENGINE_SHA256,
                "validation_report_sha256": sha256_file(report),
                "core_outputs_aggregate_sha256": hashlib.sha256(core_aggregate.encode("utf-8")).hexdigest(),
                "validated_utc": utc_now(),
            }
        ],
    )
    print(f"SP-INT-03 {stage} validation passed ({len(checks.rows)} checks)", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("v1", "v2"))
    parser.add_argument("--project-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        checks, output = validate(args.stage, build_paths(args.project_root))
        write_validation(checks, output, args.stage)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
