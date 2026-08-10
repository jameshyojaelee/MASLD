#!/usr/bin/env python3
"""Fit prespecified donor-resolved binary lipid-zone models and combine donors."""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from yakubovsky_common import (
    DATASET,
    MAX_FAILED_BOOTSTRAP_FRACTION,
    MIN_DONOR_SPOTS,
    MIN_FINAL_BLOCKS,
    MODEL_STAGE_ARTIFACTS,
    N_BOOTSTRAP,
    PRIMARY_SPLINE_DF,
    RELEASE_ID,
    SPATIAL_BLOCK_AUDIT_COLUMNS,
    SPLINE_SENSITIVITY_DF,
    ContractError,
    atomic_write_frame,
    bh_adjust,
    block_bootstrap_lipid_slope,
    build_model_design,
    build_spatial_blocks,
    default_paths,
    deterministic_seed,
    direction_label,
    expected_sign,
    ols_lipid_fit,
    sha256_file,
    signed_stouffer,
    validate_registry_contract,
    validate_source_gate_summary,
    validate_stage_chain,
    write_stage_seal,
)


VARIANTS = (
    ("primary", "primary_score", PRIMARY_SPLINE_DF, False),
    ("equal_weight", "equal_weight_score", PRIMARY_SPLINE_DF, False),
    ("leave_top_gene", "leave_top_gene_score", PRIMARY_SPLINE_DF, False),
    ("zonation_df3", "primary_score", SPLINE_SENSITIVITY_DF[0], False),
    ("zonation_df5", "primary_score", SPLINE_SENSITIVITY_DF[1], False),
    ("zonation_landmark_overadjusted", "primary_score", PRIMARY_SPLINE_DF, True),
)
ROBUST_VARIANTS = ("equal_weight", "leave_top_gene", "zonation_df3", "zonation_df5")

DONOR_EFFECT_COLUMNS = (
    "release_id",
    "dataset",
    "donor",
    "program_uid",
    "legacy_module",
    "expected_direction",
    "variant",
    "score_column",
    "spline_df",
    "source_zonation_landmark_covariate",
    "estimable",
    "failure_reason",
    "n_complete_spots",
    "n_lipid_zone_spots",
    "n_non_lipid_zone_spots",
    "n_blocks",
    "observed_beta",
    "conventional_ols_se",
    "conventional_ols_statistic",
    "conventional_ols_p_descriptive",
    "bootstrap_ci_low",
    "bootstrap_ci_high",
    "bootstrap_pvalue",
    "bootstrap_sd",
    "n_bootstrap",
    "n_failed_bootstrap",
    "failed_bootstrap_fraction",
    "design_rank",
    "design_columns",
    "residual_df",
    "condition_number",
    "bootstrap_seed",
)
DONOR_COMBINATION_COLUMNS = (
    "release_id",
    "dataset",
    "program_uid",
    "combination",
    "omitted_donor",
    "donor",
    "n_donors",
    "combined_z",
    "combined_pvalue",
    "direction",
)
MODEL_DESIGN_AUDIT_COLUMNS = (
    "audit_stage",
    "dataset",
    "donor",
    "n_spots",
    "n_final_blocks",
    "block_gate_status",
    "source_label",
    "non_lipid_zone_is_lipid_free",
    "continuous_lipid_field_read",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def fit_task(task: dict[str, Any]) -> tuple[dict[str, Any], np.ndarray, pd.DataFrame | None]:
    donor = task["donor"]
    uid = task["program_uid"]
    variant = task["variant"]
    result: dict[str, Any] = {
        "release_id": RELEASE_ID,
        "dataset": DATASET,
        "donor": donor,
        "program_uid": uid,
        "legacy_module": task["legacy_module"],
        "expected_direction": task["expected_direction"],
        "variant": variant,
        "score_column": task["score_column"],
        "spline_df": task["spline_df"],
        "source_zonation_landmark_covariate": task[
            "zonation_landmark_adjusted"
        ],
        "estimable": False,
        "failure_reason": "",
        "n_complete_spots": 0,
        "n_lipid_zone_spots": 0,
        "n_non_lipid_zone_spots": 0,
        "n_blocks": 0,
        "observed_beta": np.nan,
        "conventional_ols_se": np.nan,
        "conventional_ols_statistic": np.nan,
        "conventional_ols_p_descriptive": np.nan,
        "bootstrap_ci_low": np.nan,
        "bootstrap_ci_high": np.nan,
        "bootstrap_pvalue": np.nan,
        "bootstrap_sd": np.nan,
        "n_bootstrap": task["n_bootstrap"],
        "n_failed_bootstrap": task["n_bootstrap"],
        "failed_bootstrap_fraction": 1.0,
        "design_rank": np.nan,
        "design_columns": np.nan,
        "residual_df": np.nan,
        "condition_number": np.nan,
        "bootstrap_seed": task["seed"],
    }
    frame = task["frame"].copy()
    needed = [
        task["score_column"],
        "lipid_zone",
        "zonation_eta",
        "background_corrected_spot_sum",
        "detected_genes",
        "block_id",
        "section",
    ]
    if task["zonation_landmark_adjusted"]:
        needed.append("source_zonation_landmark_expression")
    complete = frame[needed].notna().all(axis=1) & (frame["block_id"] >= 0)
    frame = frame.loc[complete].copy()
    result["n_complete_spots"] = len(frame)
    if len(frame) < MIN_DONOR_SPOTS:
        result["failure_reason"] = "fewer_than_100_complete_spots"
        return result, np.array([], dtype=float), None
    result["n_lipid_zone_spots"] = int(frame["lipid_zone"].sum())
    result["n_non_lipid_zone_spots"] = int((~frame["lipid_zone"]).sum())
    if frame["lipid_zone"].nunique() != 2:
        result["failure_reason"] = "binary_lipid_class_degenerate"
        return result, np.array([], dtype=float), None
    result["n_blocks"] = frame["block_id"].nunique()
    if result["n_blocks"] < MIN_FINAL_BLOCKS:
        result["failure_reason"] = "fewer_than_8_final_spatial_blocks"
        return result, np.array([], dtype=float), None
    try:
        design, columns = build_model_design(
            frame["lipid_zone"].astype(int).to_numpy(),
            frame["zonation_eta"].to_numpy(float),
            frame["background_corrected_spot_sum"].to_numpy(float),
            frame["detected_genes"].to_numpy(float),
            int(task["spline_df"]),
            (
                frame["source_zonation_landmark_expression"].to_numpy(float)
                if task["zonation_landmark_adjusted"]
                else None
            ),
        )
        outcome = frame[task["score_column"]].to_numpy(float)
        observed = ols_lipid_fit(design, outcome)
        bootstrap = block_bootstrap_lipid_slope(
            design,
            outcome,
            frame["block_id"].to_numpy(int),
            frame["section"].astype(str).to_numpy(),
            int(task["n_bootstrap"]),
            int(task["seed"]),
        )
    except (ContractError, np.linalg.LinAlgError, ValueError) as error:
        result["failure_reason"] = f"model_failure:{type(error).__name__}:{error}"
        return result, np.array([], dtype=float), None
    failed_fraction = bootstrap.n_failed / task["n_bootstrap"]
    result.update(
        {
            "observed_beta": observed.beta,
            "conventional_ols_se": observed.se,
            "conventional_ols_statistic": observed.statistic,
            "conventional_ols_p_descriptive": observed.pvalue,
            "bootstrap_ci_low": bootstrap.ci_low,
            "bootstrap_ci_high": bootstrap.ci_high,
            "bootstrap_pvalue": bootstrap.pvalue,
            "bootstrap_sd": float(np.nanstd(bootstrap.estimates, ddof=1)),
            "n_failed_bootstrap": bootstrap.n_failed,
            "failed_bootstrap_fraction": failed_fraction,
            "design_rank": observed.rank,
            "design_columns": observed.n_columns,
            "residual_df": observed.residual_df,
            "condition_number": observed.condition_number,
        }
    )
    if failed_fraction > MAX_FAILED_BOOTSTRAP_FRACTION:
        result["failure_reason"] = "more_than_5_percent_bootstrap_draws_failed"
        return result, bootstrap.estimates, None
    result["estimable"] = True
    reference = None
    if task["export_reference"]:
        reference = pd.DataFrame(design, columns=columns)
        reference.insert(0, "spot_id", frame["spot_id"].to_numpy())
        reference["outcome_primary_score"] = outcome
        reference["block_id"] = frame["block_id"].to_numpy(int)
        reference["section"] = frame["section"].astype(str).to_numpy()
        reference["donor"] = donor
        reference["program_uid"] = uid
        reference["variant"] = variant
    return result, bootstrap.estimates, reference


def combine_variant(part: pd.DataFrame) -> dict[str, Any]:
    valid = part[part["estimable"].astype(bool)].copy() if len(part) else part.copy()
    z, pvalue = signed_stouffer(valid["observed_beta"], valid["bootstrap_pvalue"])
    return {
        "n_valid_donors": len(valid),
        "combined_z": z,
        "combined_pvalue": pvalue,
        "median_donor_slope": float(valid["observed_beta"].median()) if len(valid) else np.nan,
        "minimum_donor_slope": float(valid["observed_beta"].min()) if len(valid) else np.nan,
        "maximum_donor_slope": float(valid["observed_beta"].max()) if len(valid) else np.nan,
        "n_positive_donor_slopes": int((valid["observed_beta"] > 0).sum()),
        "n_negative_donor_slopes": int((valid["observed_beta"] < 0).sum()),
        "valid_donors": ";".join(valid["donor"].astype(str)),
    }


def heterogeneity(part: pd.DataFrame) -> tuple[float, int, float]:
    numeric_beta = pd.to_numeric(part["observed_beta"], errors="coerce")
    numeric_sd = pd.to_numeric(part["bootstrap_sd"], errors="coerce")
    estimable = (
        part["estimable"]
        if part["estimable"].dtype == bool
        else part["estimable"].astype(str).str.lower().isin(["true", "1"])
    )
    valid = part[
        estimable
        & np.isfinite(numeric_beta)
        & np.isfinite(numeric_sd)
        & (numeric_sd > 0)
    ].copy()
    if len(valid) < 2:
        return np.nan, max(len(valid) - 1, 0), np.nan
    weights = 1 / np.square(valid["bootstrap_sd"].to_numpy(float))
    estimates = valid["observed_beta"].to_numpy(float)
    mean = float(np.sum(weights * estimates) / np.sum(weights))
    q = float(np.sum(weights * np.square(estimates - mean)))
    df = len(valid) - 1
    return q, df, float(stats.chi2.sf(q, df))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--workers", type=int, default=max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))))
    parser.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    parser.add_argument("--fixture-mode", action="store_true")
    args = parser.parse_args()
    if not args.fixture_mode and args.n_bootstrap != N_BOOTSTRAP:
        raise ContractError(f"Production bootstrap count is frozen at {N_BOOTSTRAP}")
    if args.fixture_mode and args.n_bootstrap < 19:
        raise ContractError("Fixture mode requires at least 19 bootstrap draws")
    paths = default_paths(args.base)
    output = (args.output_dir or paths["candidate"]).resolve()
    required = [
        output / "analysis_freeze_manifest.tsv",
        output / "sample_manifest.tsv",
        output / "program_testability.tsv",
        output / "per_sample_program_scores.tsv",
        output / "scoring_execution_manifest.tsv",
    ]
    if any(not path.is_file() for path in required):
        raise ContractError("Binary model prerequisites are incomplete")
    validate_stage_chain(output, "scoring")
    final_outputs = [
        output / "spatial_block_audit.tsv",
        output / "donor_program_effects.tsv",
        output / "spatial_bootstrap_summary.tsv",
        output / "donor_combination_audit.tsv",
        output / "program_effects.tsv",
        output / "sensitivity.tsv",
        output / "multiplicity_manifest.tsv",
        output / "model_sample_manifest.tsv",
        output / "model_design_reference.tsv",
        output / "bootstrap_reference_draws.tsv",
        output / "model_design_audit.tsv",
        output / "model_execution_manifest.tsv",
    ]
    existing = [str(path) for path in final_outputs if path.exists()]
    if existing:
        raise ContractError("Refusing to overwrite model output(s): " + ", ".join(existing))

    registry = validate_registry_contract(paths)
    gate = validate_source_gate_summary(paths["source_gate"])
    sample = pd.read_csv(output / "sample_manifest.tsv", sep="\t", dtype=str, keep_default_na=False)
    prohibited_columns = {
        "lipid_percentage",
        "tissue_percentage",
        "lipid_z",
        "continuous_lipid",
        "ordered_lipid_category",
    }
    overlap = prohibited_columns.intersection(sample.columns)
    if overlap:
        raise ContractError(f"Binary pipeline refuses continuous/ordered columns: {sorted(overlap)}")
    required_sample = {
        "spot_id",
        "donor",
        "section",
        "x_coordinate",
        "y_coordinate",
        "zonation_eta",
        "source_defined_binary_lipid_class",
        "binary_lipid_inference_eligible",
        "analysis_eligible",
        "background_corrected_spot_sum",
        "detected_genes",
        "source_zonation_landmark_expression",
    }
    if not required_sample.issubset(sample.columns):
        raise ContractError(f"Sample manifest lacks {sorted(required_sample)}")
    if sample["spot_id"].duplicated().any():
        raise ContractError("Sample manifest spot IDs are not unique")
    for column in (
        "x_coordinate",
        "y_coordinate",
        "zonation_eta",
        "background_corrected_spot_sum",
        "detected_genes",
        "source_zonation_landmark_expression",
    ):
        sample[column] = pd.to_numeric(sample[column], errors="coerce")
    sample["analysis_eligible"] = sample["analysis_eligible"].str.lower().isin(["true", "1"])
    sample["binary_lipid_inference_eligible"] = sample[
        "binary_lipid_inference_eligible"
    ].str.lower().isin(["true", "1"])
    eligible_donors = tuple(sorted(gate.passing_donors))
    observed_eligible = tuple(sorted(sample.loc[sample["binary_lipid_inference_eligible"], "donor"].unique()))
    if observed_eligible != eligible_donors:
        raise ContractError(
            f"Adapter binary-donor set {observed_eligible} != Plan 10 gate {eligible_donors}"
        )
    inference = sample[
        sample["analysis_eligible"] & sample["binary_lipid_inference_eligible"]
    ].copy()
    if not inference["source_defined_binary_lipid_class"].isin(
        ["lipid_zone", "non_lipid_zone"]
    ).all():
        raise ContractError("Inference rows contain an invalid binary lipid class")
    inference["lipid_zone"] = inference["source_defined_binary_lipid_class"].eq("lipid_zone")

    block_rows: list[pd.DataFrame] = []
    donor_block_status: dict[str, str] = {}
    inference["block_id"] = -1
    for donor, part in inference.groupby("donor", sort=True):
        try:
            built = build_spatial_blocks(
                donor,
                part["section"].astype(str).to_numpy(),
                part["x_coordinate"].to_numpy(float),
                part["y_coordinate"].to_numpy(float),
            )
            inference.loc[part.index, "block_id"] = built.labels
            block_rows.append(built.audit)
            donor_block_status[donor] = (
                "pass" if built.n_blocks >= MIN_FINAL_BLOCKS else "fewer_than_8_final_blocks"
            )
        except ContractError as error:
            donor_block_status[donor] = f"block_failure:{error}"
    block_audit = (
        pd.concat(block_rows, ignore_index=True)
        if block_rows
        else pd.DataFrame(columns=SPATIAL_BLOCK_AUDIT_COLUMNS)
    )
    valid_block_donors = [
        donor for donor, status in donor_block_status.items() if status == "pass"
    ]

    scores = pd.read_csv(
        output / "per_sample_program_scores.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    required_score = {
        "spot_id",
        "donor",
        "section",
        "program_uid",
        "primary_score",
        "equal_weight_score",
        "leave_top_gene_score",
    }
    if not required_score.issubset(scores.columns):
        raise ContractError(f"Score table lacks {sorted(required_score)}")
    for column in ("primary_score", "equal_weight_score", "leave_top_gene_score"):
        scores[column] = pd.to_numeric(scores[column], errors="coerce")
    if scores.duplicated(["spot_id", "program_uid"]).any():
        raise ContractError("Score table is not unique by spot/program")
    score_model_columns = [
        "spot_id",
        "donor",
        "section",
        "program_uid",
        "legacy_module",
        "primary_score",
        "equal_weight_score",
        "leave_top_gene_score",
    ]
    model_data = scores[score_model_columns].merge(
        inference[
            [
                "spot_id",
                "donor",
                "section",
                "lipid_zone",
                "zonation_eta",
                "background_corrected_spot_sum",
                "detected_genes",
                "source_zonation_landmark_expression",
                "block_id",
            ]
        ],
        on=["spot_id", "donor", "section"],
        how="inner",
        validate="many_to_one",
    )
    testability = pd.read_csv(
        output / "program_testability.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    testability["testable"] = testability["testable"].str.lower().isin(["true", "1"])
    family = registry.family.set_index("program_uid", drop=False)
    expected_programs = set(testability.loc[testability["testable"], "program_uid"])
    if not expected_programs.issubset(family.index):
        raise ContractError("A testable program is outside the sealed family")

    # Freeze the exact spot universe used by the primary models.  This is
    # separate from the source-native sample manifest because small disconnected
    # islands and incomplete model rows are valid adapter observations but not
    # inferential technical units.  Plan 13 consumes this table instead of
    # guessing spot inclusion from donor-level summaries.
    model_sample = inference[
        ["spot_id", "donor", "section", "block_id"]
    ].copy()
    model_sample["block_uid"] = [
        f"{donor}:{section}:block_{int(block_id)}" if int(block_id) >= 0 else ""
        for donor, section, block_id in model_sample[
            ["donor", "section", "block_id"]
        ].itertuples(index=False, name=None)
    ]
    score_presence = model_data[model_data["program_uid"].isin(expected_programs)].copy()
    score_presence["finite_primary_score"] = np.isfinite(score_presence["primary_score"])
    score_counts = score_presence.groupby("spot_id", sort=False).agg(
        n_testable_program_rows=("program_uid", "nunique"),
        n_finite_primary_scores=("finite_primary_score", "sum"),
    )
    model_sample = model_sample.merge(
        score_counts,
        left_on="spot_id",
        right_index=True,
        how="left",
        validate="one_to_one",
    )
    for column in ("n_testable_program_rows", "n_finite_primary_scores"):
        model_sample[column] = model_sample[column].fillna(0).astype(int)
    model_sample["n_expected_testable_programs"] = len(expected_programs)
    model_sample["primary_covariates_complete"] = inference[
        [
            "lipid_zone",
            "zonation_eta",
            "background_corrected_spot_sum",
            "detected_genes",
        ]
    ].notna().all(axis=1).to_numpy()
    model_sample["included_spatial_block"] = model_sample["block_id"] >= 0
    model_sample["primary_model_eligible"] = (
        bool(expected_programs)
        & model_sample["primary_covariates_complete"]
        & model_sample["included_spatial_block"]
        & (model_sample["n_testable_program_rows"] == len(expected_programs))
        & (model_sample["n_finite_primary_scores"] == len(expected_programs))
    )
    model_sample["primary_model_exclusion_reason"] = np.select(
        [
            ~model_sample["included_spatial_block"],
            ~model_sample["primary_covariates_complete"],
            model_sample["n_testable_program_rows"] != len(expected_programs),
            model_sample["n_finite_primary_scores"] != len(expected_programs),
        ],
        [
            "excluded_small_tissue_island",
            "incomplete_primary_covariates",
            "incomplete_frozen_program_family",
            "nonfinite_primary_program_score",
        ],
        default="included_primary_model" if expected_programs else "no_testable_program",
    )
    primary_spot_ids = set(
        model_sample.loc[model_sample["primary_model_eligible"], "spot_id"]
    )
    model_data = model_data[model_data["spot_id"].isin(primary_spot_ids)].copy()

    tasks: list[dict[str, Any]] = []
    for donor in sorted(valid_block_donors):
        donor_data = model_data[model_data["donor"] == donor].copy()
        for uid in sorted(expected_programs):
            program_data = donor_data[donor_data["program_uid"] == uid].copy()
            row = family.loc[uid]
            for (
                variant,
                score_column,
                spline_df,
                zonation_landmark_adjusted,
            ) in VARIANTS:
                tasks.append(
                    {
                        "donor": donor,
                        "program_uid": uid,
                        "legacy_module": row["module"],
                        "expected_direction": row["expected_direction"],
                        "variant": variant,
                        "score_column": score_column,
                        "spline_df": spline_df,
                        "zonation_landmark_adjusted": zonation_landmark_adjusted,
                        "n_bootstrap": args.n_bootstrap,
                        "seed": deterministic_seed(RELEASE_ID, donor, uid, variant),
                        "frame": program_data,
                        # Export every successful primary design to the parent
                        # process, then choose the lexicographically first one.
                        # This prevents a failed first task from silently
                        # disabling independent model/bootstrap validation.
                        "export_reference": variant == "primary",
                    }
                )

    effect_rows: list[dict[str, Any]] = []
    reference_design: pd.DataFrame | None = None
    reference_draws: pd.DataFrame | None = None
    reference_candidates: list[tuple[tuple[str, str], pd.DataFrame, np.ndarray, dict[str, Any]]] = []
    if tasks:
        with ProcessPoolExecutor(max_workers=max(1, min(args.workers, len(tasks)))) as executor:
            futures = {executor.submit(fit_task, task): task for task in tasks}
            for future in as_completed(futures):
                row, estimates, design = future.result()
                effect_rows.append(row)
                if design is not None:
                    reference_candidates.append(
                        ((str(row["donor"]), str(row["program_uid"])), design, estimates, row)
                    )
    if reference_candidates:
        _key, reference_design, reference_estimates, reference_row = sorted(
            reference_candidates, key=lambda item: item[0]
        )[0]
        reference_draws = pd.DataFrame(
            {
                "bootstrap_draw_1based": np.arange(1, len(reference_estimates) + 1),
                "lipid_zone_beta": reference_estimates,
                "bootstrap_seed": reference_row["bootstrap_seed"],
                "donor": reference_row["donor"],
                "program_uid": reference_row["program_uid"],
                "variant": reference_row["variant"],
            }
        )
    donor_effects = pd.DataFrame(effect_rows)
    if donor_effects.empty:
        donor_effects = pd.DataFrame(columns=DONOR_EFFECT_COLUMNS)
    else:
        donor_effects = donor_effects.sort_values(
            ["program_uid", "donor", "variant"]
        ).reset_index(drop=True)
    donor_effects["estimable"] = donor_effects["estimable"].astype(bool)

    sensitivity_rows: list[dict[str, Any]] = []
    combination_rows: list[dict[str, Any]] = []
    program_rows: list[dict[str, Any]] = []
    for uid, family_row in family.iterrows():
        test_row = testability[testability["program_uid"] == uid]
        if len(test_row) != 1:
            raise ContractError(f"Missing/duplicate testability row for {uid}")
        testable = bool(test_row.iloc[0]["testable"])
        variant_combined: dict[str, dict[str, Any]] = {}
        for variant, _, _, _ in VARIANTS:
            part = donor_effects[
                (donor_effects["program_uid"] == uid)
                & (donor_effects["variant"] == variant)
            ].copy()
            combined = combine_variant(part) if len(part) else combine_variant(
                pd.DataFrame(columns=["estimable", "observed_beta", "bootstrap_pvalue", "donor"])
            )
            variant_combined[variant] = combined
            sensitivity_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "dataset": DATASET,
                    "program_uid": uid,
                    "legacy_module": family_row["module"],
                    "variant": variant,
                    **combined,
                    "combined_qvalue": np.nan,
                    "expected_direction": family_row["expected_direction"],
                    "direction_matches_expected": (
                        direction_label(combined["combined_z"])
                        == family_row["expected_direction"]
                    ),
                }
            )
        primary_part = donor_effects[
            (donor_effects["program_uid"] == uid)
            & (donor_effects["variant"] == "primary")
            & donor_effects["estimable"]
        ].copy()
        primary = variant_combined["primary"]
        q_heterogeneity, q_df, q_p = heterogeneity(primary_part)
        lodo_directions: list[str] = []
        for donor in primary_part["donor"].astype(str):
            retained = primary_part[primary_part["donor"].astype(str) != donor]
            lodo_z, lodo_p = signed_stouffer(
                retained["observed_beta"], retained["bootstrap_pvalue"]
            )
            lodo_direction = direction_label(lodo_z)
            lodo_directions.append(lodo_direction)
            combination_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "dataset": DATASET,
                    "program_uid": uid,
                    "combination": "leave_one_donor_out",
                    "omitted_donor": donor,
                    "n_donors": len(retained),
                    "combined_z": lodo_z,
                    "combined_pvalue": lodo_p,
                    "direction": lodo_direction,
                }
            )
        for donor_row in primary_part.itertuples(index=False):
            donor_z = np.sign(donor_row.observed_beta) * stats.norm.isf(
                donor_row.bootstrap_pvalue / 2
            )
            combination_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "dataset": DATASET,
                    "program_uid": uid,
                    "combination": "donor_component",
                    "omitted_donor": "",
                    "donor": donor_row.donor,
                    "n_donors": 1,
                    "combined_z": donor_z,
                    "combined_pvalue": donor_row.bootstrap_pvalue,
                    "direction": direction_label(donor_row.observed_beta),
                }
            )
        program_rows.append(
            {
                "release_id": RELEASE_ID,
                "dataset": DATASET,
                "assay": "Visium",
                "biological_unit": "donor",
                "contrast_exposure": "source_defined_lipid_zone_vs_non_lipid_zone",
                "effect_unit": "source_lipid_zone_minus_non_lipid_zone_score",
                "program_uid": uid,
                "legacy_module": family_row["module"],
                "program_name": family_row["module_name"],
                "expected_direction": family_row["expected_direction"],
                "testable": testable,
                "testability_reason": test_row.iloc[0]["testability_reason"],
                "estimate": primary["median_donor_slope"],
                "se": np.nan,
                "ci_lower": primary["minimum_donor_slope"],
                "ci_upper": primary["maximum_donor_slope"],
                "interval_type": "donor_slope_range_not_confidence_interval",
                "combined_z": primary["combined_z"],
                "pvalue": primary["combined_pvalue"],
                "qvalue": np.nan,
                "n_donors": primary["n_valid_donors"],
                "n_technical_units": int(
                    inference[inference["donor"].isin(primary_part["donor"])][
                        ["donor", "section"]
                    ].drop_duplicates().shape[0]
                ),
                "n_spots": int(primary_part["n_complete_spots"].sum()) if len(primary_part) else 0,
                "n_mapped_genes": test_row.iloc[0]["n_mapped_genes"],
                "retained_l1_weight": test_row.iloc[0]["retained_original_l1_weight"],
                "n_expected_direction_donors": int(
                    (
                        np.sign(primary_part["observed_beta"].to_numpy(float))
                        == expected_sign(family_row["expected_direction"])
                    ).sum()
                ),
                "required_expected_direction_donors": math.ceil(
                    0.75 * len(gate.passing_donors)
                ),
                "primary_direction_matches_expected": (
                    direction_label(primary["combined_z"])
                    == family_row["expected_direction"]
                ),
                "combined_direction": direction_label(primary["combined_z"]),
                "median_donor_slope_direction": direction_label(
                    primary["median_donor_slope"]
                ),
                "combined_and_median_direction_agree": (
                    direction_label(primary["combined_z"])
                    == direction_label(primary["median_donor_slope"])
                ),
                "all_lodo_direction_match": bool(lodo_directions)
                and all(direction == family_row["expected_direction"] for direction in lodo_directions),
                "cochran_q_descriptive": q_heterogeneity,
                "cochran_q_df": q_df,
                "cochran_q_p_descriptive": q_p,
                "robust": False,
                "evidence_state": "pending_multiplicity",
                "robustness_status": "pending_multiplicity",
                "source_dependence_class": "source_defined_binary_loupe_annotation",
                "registry_sha256": registry.registry_sha256,
                "membership_sha256": family_row["membership_sha256"],
                "continuous_lipid_fields_read": False,
            }
        )

    sensitivity = pd.DataFrame(sensitivity_rows)
    programs = pd.DataFrame(program_rows)
    assay_testable_uids = set(
        testability.loc[testability["testable"].astype(bool), "program_uid"]
    )
    multiplicity_rows: list[dict[str, Any]] = []
    complete_variants: set[str] = set()
    if len(sensitivity):
        for variant in sensitivity["variant"].unique():
            family_mask = (sensitivity["variant"] == variant) & sensitivity[
                "program_uid"
            ].isin(assay_testable_uids)
            estimable = family_mask & (sensitivity["n_valid_donors"] >= 3) & np.isfinite(
                sensitivity["combined_pvalue"]
            )
            family_complete = bool(assay_testable_uids) and int(estimable.sum()) == len(
                assay_testable_uids
            )
            if family_complete:
                sensitivity.loc[estimable, "combined_qvalue"] = bh_adjust(
                    sensitivity.loc[estimable, "combined_pvalue"].to_numpy(float)
                )
                complete_variants.add(str(variant))
            multiplicity_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "dataset": DATASET,
                    "variant": variant,
                    "n_frozen_family": len(family),
                    "n_assay_testable_programs": len(assay_testable_uids),
                    "n_estimable_programs": int(estimable.sum()),
                    "assay_testable_program_uids": ";".join(sorted(assay_testable_uids)),
                    "estimable_program_uids": ";".join(
                        sorted(sensitivity.loc[estimable, "program_uid"])
                    ),
                    "family_complete": family_complete,
                    "bh_applied": family_complete,
                    "multiplicity_policy": "all_assay_testable_programs_or_no_BH",
                }
            )
    primary_family_complete = "primary" in complete_variants
    if primary_family_complete:
        valid_primary = (
            programs["program_uid"].isin(assay_testable_uids)
            & (programs["n_donors"] >= 3)
            & np.isfinite(programs["pvalue"])
        )
        programs.loc[valid_primary, "qvalue"] = bh_adjust(
            programs.loc[valid_primary, "pvalue"].to_numpy(float)
        )
    programs["primary_family_estimable"] = primary_family_complete
    for index, row in programs.iterrows():
        uid = row["program_uid"]
        if not bool(row["testable"]):
            programs.loc[index, ["evidence_state", "robustness_status"]] = [
                "untestable",
                "untestable",
            ]
            continue
        if not primary_family_complete:
            programs.loc[index, ["evidence_state", "robustness_status"]] = [
                "untestable",
                "incomplete_prespecified_primary_BH_family",
            ]
            continue
        sens = sensitivity[sensitivity["program_uid"] == uid].set_index("variant")
        sensitivity_direction_ok = all(
            variant in sens.index
            and variant in complete_variants
            and int(sens.loc[variant, "n_valid_donors"]) >= 3
            and bool(sens.loc[variant, "direction_matches_expected"])
            for variant in ROBUST_VARIANTS
        )
        donor_direction_ok = (
            int(row["n_expected_direction_donors"])
            >= int(row["required_expected_direction_donors"])
        )
        robust = bool(
            np.isfinite(row["qvalue"])
            and float(row["qvalue"]) < 0.05
            and bool(row["primary_direction_matches_expected"])
            and donor_direction_ok
            and sensitivity_direction_ok
            and bool(row["all_lodo_direction_match"])
            and bool(row["combined_and_median_direction_agree"])
        )
        programs.loc[index, "robust"] = robust
        programs.loc[index, "evidence_state"] = "robust" if robust else "tested_negative"
        programs.loc[index, "robustness_status"] = (
            "robust" if robust else "primary_nonsignificant" if float(row["qvalue"]) >= 0.05 else "failed_prespecified_robustness"
        )

    if reference_design is None:
        reference_design = pd.DataFrame(
            columns=["spot_id", "donor", "program_uid", "variant", "block_id", "section"]
        )
    if reference_draws is None:
        reference_draws = pd.DataFrame(
            columns=[
                "bootstrap_draw_1based",
                "lipid_zone_beta",
                "bootstrap_seed",
                "donor",
                "program_uid",
                "variant",
            ]
        )
    atomic_write_frame(output / "spatial_block_audit.tsv", block_audit)
    atomic_write_frame(output / "donor_program_effects.tsv", donor_effects)
    atomic_write_frame(output / "spatial_bootstrap_summary.tsv", donor_effects)
    atomic_write_frame(
        output / "donor_combination_audit.tsv",
        pd.DataFrame(combination_rows, columns=DONOR_COMBINATION_COLUMNS),
    )
    atomic_write_frame(output / "program_effects.tsv", programs)
    atomic_write_frame(output / "sensitivity.tsv", sensitivity)
    atomic_write_frame(output / "multiplicity_manifest.tsv", pd.DataFrame(multiplicity_rows))
    atomic_write_frame(output / "model_sample_manifest.tsv", model_sample)
    atomic_write_frame(output / "model_design_reference.tsv", reference_design)
    atomic_write_frame(output / "bootstrap_reference_draws.tsv", reference_draws)

    model_audit = pd.DataFrame(
        [
            {
                "audit_stage": "binary_lipid_model",
                "dataset": DATASET,
                "donor": donor,
                "n_spots": int((inference["donor"] == donor).sum()),
                "n_final_blocks": int(
                    block_audit.loc[
                        (block_audit["donor"] == donor)
                        & (block_audit["block_id"] >= 0),
                        "block_id",
                    ].nunique()
                )
                if len(block_audit)
                else 0,
                "block_gate_status": status,
                "source_label": "binary_lipid_zone_vs_non_lipid_zone",
                "non_lipid_zone_is_lipid_free": False,
                "continuous_lipid_field_read": False,
            }
            for donor, status in sorted(donor_block_status.items())
        ],
        columns=MODEL_DESIGN_AUDIT_COLUMNS,
    )
    atomic_write_frame(output / "model_design_audit.tsv", model_audit)
    execution = pd.DataFrame(
        [
            {
                "release_id": RELEASE_ID,
                "stage": "binary_lipid_spatial_bootstrap",
                "producer": str(Path(__file__).resolve()),
                "producer_sha256": sha256_file(Path(__file__).resolve()),
                "python": sys.version.replace("\n", " "),
                "platform": platform.platform(),
                "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
                "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
                "completed_utc": utc_now(),
                "workers": args.workers,
                "n_bootstrap": args.n_bootstrap,
                "fixture_mode": args.fixture_mode,
                "n_source_gate_donors": len(gate.passing_donors),
                "n_valid_block_donors": len(valid_block_donors),
                "n_frozen_family": len(family),
                "n_testable_programs": int(testability["testable"].sum()),
                "minimum_valid_primary_donors_across_testable_programs": int(
                    programs.loc[programs["testable"].astype(bool), "n_donors"].min()
                )
                if programs["testable"].astype(bool).any()
                else 0,
                "maximum_valid_primary_donors_across_testable_programs": int(
                    programs.loc[programs["testable"].astype(bool), "n_donors"].max()
                )
                if programs["testable"].astype(bool).any()
                else 0,
                "n_robust_programs": int(programs["robust"].astype(bool).sum()),
                "primary_prespecified_BH_family_complete": primary_family_complete,
                "binary_lipid_labels_read": True,
                "continuous_lipid_fields_read": False,
                "ordered_lipid_trend_tested": False,
                "image_reclassification_used": False,
                "exit_state": "pass",
            }
        ]
    )
    atomic_write_frame(output / "model_execution_manifest.tsv", execution)
    write_stage_seal(output, "model", MODEL_STAGE_ARTIFACTS)
    print(
        json.dumps(
            {
                "release_id": RELEASE_ID,
                "status": "binary_lipid_models_complete",
                "n_source_gate_donors": len(gate.passing_donors),
                "n_valid_block_donors": len(valid_block_donors),
                "n_testable_programs": int(testability["testable"].sum()),
                "n_robust_programs": int(programs["robust"].astype(bool).sum()),
                "continuous_lipid_fields_read": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
