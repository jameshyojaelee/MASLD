#!/usr/bin/env python3
"""Run the prespecified all-cell-type Visium composition sensitivity.

This is a falsification analysis for the two frozen display programs.  It
changes exactly one scientific component of the accepted Visium engine: the
spot-level nuisance design contains all 16 cell2location q05 abundance factors
instead of only the matching lineage factor.  Frozen program weights, measured
genes, six-neighbour disconnected-island graphs, matched-gene sets, seeds,
donor collapse, null counts, and complete within-dataset two-program BH
families are otherwise retained.

The result is written to an isolated candidate root and never changes the
validated spatial candidate or a canonical figure/release.
"""

from __future__ import annotations

import argparse
import gc
import os
import platform
import sys
from pathlib import Path

import numpy as np

from visium_rerun_lib import (
    N_NULL,
    N_SENSITIVITY_NULL,
    RELEASE_ID as UPSTREAM_RELEASE_ID,
    SEED,
    V1_ENGINE_SHA256,
    V2_MEMBERSHIP_SHA256,
    V2_REGISTRY_SHA256,
    build_paths,
    load_legacy_engine,
    prepare_engine_inputs,
    read_tsv,
    sha256_file,
    universe_hash,
    utc_now,
    verify_hotspot_ready,
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
    return (
        project_root
        / "Analysis/Multimodal_Program_Projection/candidates"
        / RELEASE_ID
    )


def upstream_root(paths) -> Path:
    return paths.native_root / "v2_candidate"


def verify_upstream(paths) -> Path:
    root = upstream_root(paths)
    ready = root / "READY"
    if sha256_file(ready) != UPSTREAM_READY_SHA256:
        raise RuntimeError("validated v2 candidate READY is missing or drifted")
    rows = read_tsv(ready)
    if len(rows) != 1 or rows[0].get("status") != "pass_v2_candidate":
        raise RuntimeError("validated v2 candidate does not carry pass_v2_candidate")
    manifest = root / "output_manifest.tsv"
    for row in read_tsv(manifest):
        path = root / row["relative_path"]
        if (
            not path.is_file()
            or path.stat().st_size != int(row["bytes"])
            or sha256_file(path) != row["sha256"]
        ):
            raise RuntimeError(f"upstream v2 output is missing or drifted: {path}")
    return root


def standardize_design(columns: list[np.ndarray]) -> np.ndarray:
    """Build the nuisance design using the pinned engine's convention."""
    if len(columns) < 2:
        raise ValueError("design must contain an intercept and nuisance columns")
    design = np.column_stack(columns).astype(float, copy=False)
    for column in range(1, design.shape[1]):
        mean = np.nanmean(design[:, column])
        sd = np.nanstd(design[:, column], ddof=1)
        design[:, column] = (
            (design[:, column] - mean) / sd
            if np.isfinite(sd) and sd > 0
            else 0.0
        )
    design[~np.isfinite(design)] = 0.0
    return design


def load_full_designs(paths):
    """Read only obs/obsm metadata and create a common 16-factor contract."""
    import anndata as ad
    import pandas as pd

    cached = {}
    factor_sets = []
    for dataset, path in paths.datasets.items():
        adata = ad.read_h5ad(path, backed="r")
        try:
            obs = adata.obs.copy()
            raw_names = [
                str(value).replace("means_per_cluster_mu_fg_", "")
                for value in adata.uns["mod"]["factor_names"]
            ]
            abundance = np.asarray(adata.obsm["q05_cell_abundance_w_sf"], dtype=float)
        finally:
            adata.file.close()
        if len(raw_names) != EXPECTED_FACTOR_COUNT or len(set(raw_names)) != EXPECTED_FACTOR_COUNT:
            raise RuntimeError(
                f"{dataset} must have exactly {EXPECTED_FACTOR_COUNT} unique cell2location factors; "
                f"observed {len(raw_names)} names/{len(set(raw_names))} unique"
            )
        if abundance.shape != (len(obs), EXPECTED_FACTOR_COUNT):
            raise RuntimeError(f"{dataset} q05 abundance matrix has unexpected shape {abundance.shape}")
        if not np.isfinite(abundance).all() or (abundance < 0).any():
            raise RuntimeError(f"{dataset} q05 abundance matrix is non-finite or negative")
        required = {"sample_id", "individual", "total_counts", "n_genes_by_counts"}
        missing = required.difference(obs.columns)
        if missing:
            raise RuntimeError(f"{dataset} lacks required design fields: {sorted(missing)}")
        factor_sets.append(set(raw_names))
        cached[dataset] = {
            "obs": obs,
            "raw_factor_names": raw_names,
            "abundance": abundance,
        }
    if factor_sets[0] != factor_sets[1]:
        raise RuntimeError("GSE192741 and Vu do not share the same 16 cell2location factors")
    canonical_factors = sorted(factor_sets[0])
    for dataset, item in cached.items():
        obs = item["obs"]
        reorder = [item["raw_factor_names"].index(name) for name in canonical_factors]
        abundance = item["abundance"][:, reorder]
        total = pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(float)
        detected = pd.to_numeric(obs["n_genes_by_counts"], errors="coerce").to_numpy(float)
        if (
            not np.isfinite(total).all()
            or not np.isfinite(detected).all()
            or (total < 0).any()
            or (detected < 0).any()
        ):
            raise RuntimeError(f"{dataset} has invalid library-size or detected-gene covariates")
        columns = [np.ones(len(obs), dtype=float)]
        columns.extend(abundance[:, column] for column in range(abundance.shape[1]))
        columns.extend([np.log1p(total), np.log1p(detected)])
        design = standardize_design(columns)
        item["factor_names"] = canonical_factors
        item["design"] = design
        item["rank"] = int(np.linalg.matrix_rank(design))
        item["condition_number"] = float(np.linalg.cond(design))
    return cached, canonical_factors


def install_design_override(engine, cached, dataset: str):
    """Replace only the nuisance-design constructor for one dataset call."""
    expected_obs = cached[dataset]["obs"]
    base_design = cached[dataset]["design"]

    def full_composition_design(obs, _single_lineage_abundance, zonation=None):
        if len(obs) != len(expected_obs) or not obs.index.equals(expected_obs.index):
            raise RuntimeError(f"{dataset} observation order drifted inside the pinned engine")
        if zonation is None:
            return base_design.copy()
        zonation = np.asarray(zonation, dtype=float)
        if zonation.shape != (len(obs),) or not np.isfinite(zonation).all():
            raise RuntimeError(f"{dataset} zonation sensitivity covariate is invalid")
        return standardize_design(
            [base_design[:, column] for column in range(base_design.shape[1])]
            + [zonation]
        )

    engine.design_matrix = full_composition_design


def source_manifest(paths, upstream: Path, producer: Path) -> list[dict[str, object]]:
    sources = [
        (producer, "candidate_producer"),
        (producer.with_name("visium_rerun_lib.py"), "candidate_engine_adapter"),
        (paths.v1_engine, "pinned_scientific_engine"),
        (paths.hotspot_root / "program_registry_v2.tsv", "frozen_program_registry"),
        (paths.hotspot_root / "program_membership_v2.tsv", "frozen_program_membership"),
        (paths.gene_metadata, "gene_biotype_metadata"),
        (upstream / "READY", "validated_primary_spatial_ready"),
        (upstream / "spatial_program_results.tsv", "primary_spatial_results"),
        (upstream / "spatial_null_summary.tsv", "primary_null_summary"),
        (upstream / "spatial_matching_audit.tsv", "primary_matching_audit"),
        (upstream / "spatial_graph_audit.tsv", "primary_graph_audit"),
    ]
    sources.extend((path, f"{dataset}_input_h5ad") for dataset, path in paths.datasets.items())
    rows = []
    for path, role in sources:
        if not path.is_file():
            raise RuntimeError(f"missing sensitivity source {role}: {path}")
        rows.append(
            {
                "relative_path": path.relative_to(paths.project_root).as_posix(),
                "source_role": role,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return rows


def write_output_manifest(staging: Path) -> None:
    excluded = {"output_manifest.tsv", "COMPUTE_COMPLETE", "READY", "validation_report.tsv"}
    rows = []
    for path in sorted(staging.iterdir()):
        if path.is_file() and path.name not in excluded:
            rows.append(
                {
                    "relative_path": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    write_tsv(staging / "output_manifest.tsv", ("relative_path", "bytes", "sha256"), rows)


def package_versions() -> dict[str, str]:
    import anndata
    import pandas
    import scipy

    return {
        "python": platform.python_version(),
        "anndata": anndata.__version__,
        "numpy": np.__version__,
        "pandas": pandas.__version__,
        "scipy": scipy.__version__,
    }


def run(project_root: Path | None, output_root: Path | None) -> Path:
    import pandas as pd

    paths = build_paths(project_root)
    verify_hotspot_ready(paths)
    upstream = verify_upstream(paths)
    final = (output_root or default_output_root(paths.project_root)).resolve()
    candidates = (paths.project_root / "Analysis/Multimodal_Program_Projection/candidates").resolve()
    if final.parent != candidates:
        raise RuntimeError(f"output must be a direct isolated candidate child of {candidates}")
    if final.exists():
        raise RuntimeError(f"refusing to overwrite candidate output: {final}")
    staging = final.parent / f".{final.name}.incomplete.{os.environ.get('SLURM_JOB_ID', os.getpid())}"
    if staging.exists():
        raise RuntimeError(f"refusing to reuse staging root: {staging}")
    staging.mkdir(parents=False)

    producer = Path(__file__).resolve()
    try:
        registry, membership, universe = prepare_engine_inputs(paths, "v2")
        if len(registry) != 2:
            raise RuntimeError(f"expected exactly two frozen programs, observed {len(registry)}")
        full_designs, factor_names = load_full_designs(paths)
        contract_rows = [
            ("release_id", RELEASE_ID),
            ("upstream_release_id", UPSTREAM_RELEASE_ID),
            ("analysis_role", "supplementary_falsification_sensitivity"),
            ("program_family", "two_frozen_robust_display_programs"),
            ("program_family_sha256", universe_hash(universe)),
            ("program_registry_sha256", V2_REGISTRY_SHA256),
            ("program_membership_sha256", V2_MEMBERSHIP_SHA256),
            ("score_weights", "frozen_positive_original_l1_renormalized_over_measured_genes"),
            ("residualization", "all_16_cell2location_q05_factors_plus_log1p_total_counts_and_log1p_detected_genes"),
            ("cell2location_factor_count", EXPECTED_FACTOR_COUNT),
            ("cell2location_factor_order", "|".join(factor_names)),
            ("graph", "six_neighbor_disconnected_tissue_islands"),
            ("gse_inferential_unit", "donor_first_5_sections_to_4_donors"),
            ("vu_inferential_unit", "10_unresolved_physical_arrays_source_dependent"),
            ("matched_set_policy", "reuse_original_expression_detection_biotype_mt_ribo_hepatocyte_lineage_match"),
            ("matched_set_seed_gse", SEED),
            ("matched_set_seed_vu", SEED + 1000),
            ("n_null", N_NULL),
            ("n_sensitivity_null", N_SENSITIVITY_NULL),
            ("multiplicity", "BH_complete_two_program_family_within_dataset"),
            ("support_rule", "q_lt_0.05_and_equal_weight_and_leave_top_centered_sign_agreement"),
            ("canonical_write_allowed", "FALSE"),
        ]
        write_tsv(
            staging / "analysis_contract.tsv",
            ("parameter", "value"),
            [{"parameter": key, "value": value} for key, value in contract_rows],
        )
        write_tsv(
            staging / "source_manifest.tsv",
            ("relative_path", "source_role", "bytes", "sha256"),
            source_manifest(paths, upstream, producer),
        )

        gene_meta = pd.read_csv(
            paths.gene_metadata,
            sep="\t",
            usecols=["gene_name", "gene_biotype"],
            compression="gzip",
        ).drop_duplicates("gene_name")
        biotype = dict(zip(gene_meta["gene_name"].astype(str), gene_meta["gene_biotype"].astype(str)))
        engine = load_legacy_engine(paths)
        original_standardized_against_null = engine.standardized_against_null
        results_list = []
        units_list = []
        matching_list = []
        graphs_list = []
        nulls_list = []
        null_draw_rows = []
        for dataset, path in paths.datasets.items():
            install_design_override(engine, full_designs, dataset)
            captured_nulls = []

            def capture_standardized_against_null(observed, null_values):
                captured_nulls.append(np.asarray(null_values, dtype=float).copy())
                return original_standardized_against_null(observed, null_values)

            engine.standardized_against_null = capture_standardized_against_null
            result, units, matching, graphs, nulls = engine.process_dataset(
                dataset, path, registry, membership, biotype
            )
            expected_calls = len(registry) * 3
            if len(captured_nulls) != expected_calls:
                raise RuntimeError(
                    f"{dataset} expected {expected_calls} raw/residual/zonation null calls; "
                    f"observed {len(captured_nulls)}"
                )
            for program_index, program_id in enumerate(registry["program_id"].astype(str)):
                residual_draws = captured_nulls[program_index * 3 + 1]
                if len(residual_draws) != N_NULL or not np.isfinite(residual_draws).all():
                    raise RuntimeError(f"{dataset}/{program_id} residual null draws are incomplete")
                null_draw_rows.extend(
                    {
                        "release_id": RELEASE_ID,
                        "dataset": dataset,
                        "program_id": program_id,
                        "draw_id": draw_id,
                        "residual_moran_i": value,
                    }
                    for draw_id, value in enumerate(residual_draws, start=1)
                )
            results_list.append(result)
            units_list.append(units)
            matching_list.append(matching)
            graphs_list.append(graphs)
            nulls_list.append(nulls.loc[nulls["statistic"] == "residual_moran_i"].copy())
            gc.collect()

        results = pd.concat(results_list, ignore_index=True)
        units = pd.concat(units_list, ignore_index=True)
        matching = pd.concat(matching_list, ignore_index=True)
        graphs = pd.concat(graphs_list, ignore_index=True)
        nulls = pd.concat(nulls_list, ignore_index=True)
        if len(results) != 4 or not results["testable"].all():
            raise RuntimeError("full-composition result must contain four testable dataset/program rows")

        output_results = results[
            [
                "program_id", "dataset", "display_order", "cell_type", "module",
                "program_name", "n_measured", "retained_l1_weight", "testable", "n_null",
                "residual_moran_i", "residual_null_mean", "residual_null_sd",
                "residual_moran_z", "residual_pvalue", "residual_padj",
                "equal_residual_moran_i", "equal_residual_null_mean",
                "leave_top_residual_moran_i", "leave_top_residual_null_mean",
                "matched_set_sha256", "target_lineage_corr_mean",
                "sensitivity_sign_agree", "robust",
            ]
        ].copy()
        output_results.insert(0, "release_id", RELEASE_ID)
        output_results["residualization"] = "all_16_cell2location_q05_plus_two_qc_covariates"
        output_results["biological_unit"] = np.where(
            output_results["dataset"] == "GSE192741", "donor", "unresolved_physical_array"
        )
        output_results["within_source_support"] = output_results["robust"]
        output_results["evidence_state"] = np.where(
            output_results["dataset"] == "Vu_et_al_2025",
            "source_dependent",
            np.where(output_results["robust"], "supported", "indeterminate"),
        )
        output_results.to_csv(staging / PROGRAM_RESULTS, sep="\t", index=False)

        unit_out = units[
            ["program_id", "dataset", "sample_id", "individual", "condition", "residual_moran_i"]
        ].copy()
        unit_out.insert(0, "release_id", RELEASE_ID)
        unit_out["unit_role"] = np.where(
            unit_out["dataset"] == "GSE192741", "physical_section", "unresolved_physical_array"
        )
        unit_out["population_inference_authorized"] = unit_out["dataset"] == "GSE192741"
        unit_out.to_csv(staging / UNIT_RESULTS, sep="\t", index=False)

        nulls.insert(0, "release_id", RELEASE_ID)
        nulls.to_csv(staging / NULL_SUMMARY, sep="\t", index=False)
        pd.DataFrame(null_draw_rows).to_csv(
            staging / "full_composition_null_draws.tsv.gz",
            sep="\t",
            index=False,
            compression={"method": "gzip", "mtime": 0},
        )
        matching.insert(0, "release_id", RELEASE_ID)
        matching.to_csv(staging / MATCHING_AUDIT, sep="\t", index=False)
        graphs.insert(0, "release_id", RELEASE_ID)
        graphs.to_csv(staging / GRAPH_AUDIT, sep="\t", index=False)

        design_rows = []
        for dataset, item in full_designs.items():
            obs = item["obs"]
            design_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "dataset": dataset,
                    "n_spots": len(obs),
                    "n_cell2location_factors": EXPECTED_FACTOR_COUNT,
                    "n_design_columns": item["design"].shape[1],
                    "design_rank": item["rank"],
                    "condition_number": item["condition_number"],
                    "factor_order": "|".join(factor_names),
                    "n_technical_units": obs["sample_id"].astype(str).nunique(),
                    "n_resolved_biological_units": (
                        obs["individual"].astype(str).nunique() if dataset == "GSE192741" else ""
                    ),
                    "unit_gate": "pass_donor_first" if dataset == "GSE192741" else "source_dependent",
                }
            )
        write_tsv(
            staging / "full_composition_design_audit.tsv",
            (
                "release_id", "dataset", "n_spots", "n_cell2location_factors",
                "n_design_columns", "design_rank", "condition_number", "factor_order",
                "n_technical_units", "n_resolved_biological_units", "unit_gate",
            ),
            design_rows,
        )

        primary = pd.read_csv(upstream / "spatial_program_results.tsv", sep="\t")
        primary = primary[[
            "program_id", "dataset", "residual_moran_i", "residual_null_mean",
            "residual_null_sd", "residual_pvalue", "residual_padj", "matched_set_sha256", "robust",
        ]].rename(columns=lambda value: f"primary_{value}" if value not in {"program_id", "dataset"} else value)
        compare = output_results.merge(primary, on=["program_id", "dataset"], validate="one_to_one")
        compare["primary_centered_moran_i"] = compare["primary_residual_moran_i"] - compare["primary_residual_null_mean"]
        compare["full_composition_centered_moran_i"] = compare["residual_moran_i"] - compare["residual_null_mean"]
        compare["centered_change"] = compare["full_composition_centered_moran_i"] - compare["primary_centered_moran_i"]
        compare["matched_set_hash_identical"] = compare["matched_set_sha256"] == compare["primary_matched_set_sha256"]
        compare[
            [
                "release_id", "program_id", "dataset", "primary_centered_moran_i",
                "full_composition_centered_moran_i", "centered_change",
                "primary_residual_pvalue", "primary_residual_padj", "residual_pvalue",
                "residual_padj", "primary_robust", "robust", "matched_set_sha256",
                "primary_matched_set_sha256", "matched_set_hash_identical", "evidence_state",
            ]
        ].to_csv(staging / "comparison_to_primary.tsv", sep="\t", index=False)

        reuse_rows = []
        for row in compare.itertuples(index=False):
            reuse_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "dataset": row.dataset,
                    "program_id": row.program_id,
                    "original_matched_set_sha256": row.primary_matched_set_sha256,
                    "sensitivity_matched_set_sha256": row.matched_set_sha256,
                    "identical": str(bool(row.matched_set_hash_identical)).upper(),
                    "scientific_basis": (
                        "same frozen genes/weights, dataset gene universe, expression/detection/biotype/"
                        "mt/ribo/hepatocyte-lineage matching, exclusion union, seed, and draw order; only residual design changed"
                    ),
                }
            )
        write_tsv(
            staging / "matched_set_reuse_audit.tsv",
            (
                "release_id", "dataset", "program_id", "original_matched_set_sha256",
                "sensitivity_matched_set_sha256", "identical", "scientific_basis",
            ),
            reuse_rows,
        )

        execution = [
            ("release_id", RELEASE_ID),
            ("upstream_ready_sha256", UPSTREAM_READY_SHA256),
            ("producer_sha256", sha256_file(producer)),
            ("v1_engine_sha256", V1_ENGINE_SHA256),
            ("program_family_sha256", universe_hash(universe)),
            ("canonical_write_allowed", "FALSE"),
            ("completed_utc", utc_now()),
            ("slurm_job_id", os.environ.get("SLURM_JOB_ID", "not_slurm")),
        ] + sorted(package_versions().items())
        write_tsv(
            staging / "execution_manifest.tsv",
            ("parameter", "value"),
            [{"parameter": key, "value": value} for key, value in execution],
        )
        write_output_manifest(staging)
        write_tsv(
            staging / "COMPUTE_COMPLETE",
            ("release_id", "status", "output_manifest_sha256", "completed_utc"),
            [{
                "release_id": RELEASE_ID,
                "status": "compute_complete_pending_independent_validation",
                "output_manifest_sha256": sha256_file(staging / "output_manifest.tsv"),
                "completed_utc": utc_now(),
            }],
        )
        os.replace(staging, final)
    except Exception:
        if staging.is_dir() and not (staging / "FAILED").exists():
            write_tsv(staging / "FAILED", ("status", "failed_utc"), [{"status": "failed", "failed_utc": utc_now()}])
        raise
    print(f"full-composition sensitivity compute complete: {final}", flush=True)
    return final


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        run(args.project_root, args.output_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
