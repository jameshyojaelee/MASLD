#!/usr/bin/env python3
"""Independently validate and seal the Plan 11 candidate bundle."""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
from scipy import stats

from yakubovsky_common import (
    ALLOWED_LIPID_CLASSES,
    DATASET,
    MIN_DONOR_SPOTS,
    MIN_PROGRAM_GENES,
    RELEASE_ID,
    ROBUST_RULE,
    TERMINAL_STAGE_ARTIFACTS,
    ContractError,
    atomic_write_frame,
    bulk_decode_matlab_cellstr,
    bool_value,
    decode_matlab_char,
    default_paths,
    direction_label,
    expected_sign,
    parse_h5ad_sparse,
    read_h5ad_dataframe_column,
    read_h5ad_dataframe_index,
    sha256_file,
    gene_axis_sha256,
    validate_gene_axis_manifest,
    validate_registry_contract,
    validate_source_gate_summary,
    validate_stage_chain,
    write_stage_seal,
)


TRUSTED_MAPPING_STATUS = "gencode_v49_unique_symbol_confirmed"


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def h5_scalar(dataset: h5py.Dataset) -> Any:
    value = dataset[()]
    if isinstance(value, np.ndarray) and value.shape == ():
        value = value.item()
    if isinstance(value, bytes):
        return value.decode()
    return value


def independent_bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1)
    return restored


def independent_stouffer(slopes: np.ndarray, pvalues: np.ndarray) -> tuple[float, float]:
    z = np.sign(slopes) * stats.norm.isf(np.clip(pvalues, 1e-300, 1) / 2)
    combined = float(np.sum(z) / math.sqrt(len(z)))
    return combined, float(2 * stats.norm.sf(abs(combined)))


def add_check(
    checks: list[dict[str, Any]],
    check_id: str,
    passed: bool,
    observed: Any,
    expected: Any,
    detail: str,
) -> None:
    checks.append(
        {
            "check_id": check_id,
            "pass": bool(passed),
            "observed": observed,
            "expected": expected,
            "detail": detail,
        }
    )
    if not passed:
        raise ContractError(f"Validation failed [{check_id}]: {detail}; observed={observed}")


def verify_v1(paths: dict[str, Path]) -> int:
    table = pd.read_csv(paths["v1_preservation"], sep="\t", dtype=str, keep_default_na=False)
    if len(table) != 44:
        raise ContractError(f"Expected 44 protected v1 files, found {len(table)}")
    for row in table.itertuples(index=False):
        target = paths["base"] / row.relative_path
        if row.unchanged != "TRUE" or sha256_file(target) != row.final_sha256:
            raise ContractError(f"Protected v1 drift: {target}")
    return len(table)


def source_abundance_fixture(
    paths: dict[str, Path], sample: pd.DataFrame, output: Path
) -> dict[str, float | int]:
    gene_audit = pd.read_csv(
        output / "gene_mapping_audit.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    source_to_symbol = gene_audit.set_index("source_gene")["canonical_symbol"].to_dict()
    spot_to_row = {spot_id: index for index, spot_id in enumerate(sample["spot_id"])}
    with h5py.File(output / "yakubovsky_human_adapter.h5ad", "r") as adapter:
        abundance = parse_h5ad_sparse(
            adapter["layers/source_background_corrected_abundance"]
        )
        transformed = parse_h5ad_sparse(
            adapter["layers/log1p_source_abundance_per_million"]
        )
        symbols = read_h5ad_dataframe_index(adapter["var"]).astype(str)
    symbol_index = {symbol: index for index, symbol in enumerate(symbols)}
    spot_sum_errors: list[float] = []
    abundance_errors: list[float] = []
    transform_errors: list[float] = []
    fractional_values = 0
    n_spot_fixtures = 0
    n_gene_fixtures = 0
    freeze = pd.read_csv(
        output / "analysis_freeze_manifest.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    v_mat_rows = freeze[
        (freeze["record_type"] == "file_input")
        & (freeze["name"] == "yakubovsky_v_mat")
    ]
    if len(v_mat_rows) != 1:
        raise ContractError("Freeze lacks one exact Yakubovsky v.mat input row")
    gene_axis = validate_gene_axis_manifest(
        paths["gene_axis_manifest"], v_mat_rows.iloc[0]["sha256"]
    )
    with h5py.File(paths["v_mat"], "r") as handle:
        references = np.asarray(handle["v"][...]).reshape(-1, order="F")
        fixture_groups: list[tuple[h5py.Group, str]] = []
        axes: dict[str, h5py.Dataset] = {}
        for reference, expected_axis in zip(
            references[:2], gene_axis.rows.itertuples(index=False), strict=True
        ):
            group = handle[reference]
            patient = decode_matlab_char(group["patient"])
            if patient != expected_axis.donor or group.name != expected_axis.group_path:
                raise ContractError("Source fixture donor/group order drift")
            axes[f"barcodes::{patient}"] = group["spot_name"]
            if patient == gene_axis.master_donor:
                axes["master_gene_axis"] = group["gene_name"]
            fixture_groups.append((group, patient))
        resolved = bulk_decode_matlab_cellstr(axes, handle, allow_null=False)
        source_genes = resolved.values["master_gene_axis"]
        if (
            len(source_genes) != gene_axis.n_genes
            or gene_axis_sha256(source_genes) != gene_axis.common_sha256
        ):
            raise ContractError("Source fixture master gene-axis certificate drift")
        for group, patient in fixture_groups:
            barcodes = resolved.values[f"barcodes::{patient}"]
            n_genes = len(source_genes)
            matrix = group["mat"]
            n_fixture_spots = min(5, len(barcodes))
            if matrix.shape == (len(barcodes), n_genes):
                source_block = np.asarray(
                    matrix[:n_fixture_spots, :], dtype=float
                )
            elif matrix.shape == (n_genes, len(barcodes)):
                source_block = np.asarray(
                    matrix[:, :n_fixture_spots], dtype=float
                ).T
            else:
                raise ContractError("Source fixture matrix axes do not align")
            for index in range(n_fixture_spots):
                source_values = source_block[index].reshape(-1)
                if not np.isfinite(source_values).all() or (source_values < 0).any():
                    raise ContractError("Deposited source-abundance fixture is invalid")
                fractional_values += int(
                    np.sum(
                        (source_values != 0)
                        & ~np.isclose(source_values, np.round(source_values), atol=1e-10, rtol=0)
                    )
                )
                spot_id = f"{patient}:{barcodes[index]}"
                if spot_id not in spot_to_row:
                    raise ContractError(f"Source fixture spot is absent from adapter: {spot_id}")
                row_index = spot_to_row[spot_id]
                source_sum = float(
                    np.sum(source_values, dtype=np.longdouble)
                )
                spot_sum_errors.append(
                    abs(
                        float(sample.loc[row_index, "background_corrected_spot_sum"])
                        - source_sum
                    )
                )
                expected_by_symbol: dict[str, float] = {}
                for source_gene, value in zip(source_genes, source_values, strict=True):
                    symbol = source_to_symbol.get(source_gene, "")
                    if symbol and symbol in symbol_index:
                        expected_by_symbol[symbol] = expected_by_symbol.get(symbol, 0.0) + float(value)
                for symbol in sorted(expected_by_symbol)[:5]:
                    expected = expected_by_symbol[symbol]
                    observed = float(abundance[row_index, symbol_index[symbol]])
                    abundance_errors.append(abs(observed - expected))
                    expected_transform = (
                        math.log1p(1e6 * expected / source_sum) if source_sum > 0 else 0.0
                    )
                    transform_errors.append(
                        abs(
                            float(transformed[row_index, symbol_index[symbol]])
                            - expected_transform
                        )
                    )
                    n_gene_fixtures += 1
                n_spot_fixtures += 1
    return {
        "n_spot_fixtures": n_spot_fixtures,
        "max_spot_sum_error": max(spot_sum_errors) if spot_sum_errors else np.nan,
        "n_gene_fixtures": n_gene_fixtures,
        "max_abundance_error": max(abundance_errors) if abundance_errors else np.nan,
        "max_transform_error": max(transform_errors) if transform_errors else np.nan,
        "n_fractional_source_values": fractional_values,
    }


def independent_score_checks(
    paths: dict[str, Path],
    output: Path,
    registry,
) -> tuple[int, float, int]:
    adapter = output / "yakubovsky_human_adapter.h5ad"
    with h5py.File(adapter, "r") as handle:
        norm = parse_h5ad_sparse(
            handle["layers/log1p_source_abundance_per_million"]
        )
        donors = read_h5ad_dataframe_column(handle["obs"], "donor").astype(str)
        eligible_raw = read_h5ad_dataframe_column(handle["obs"], "analysis_eligible")
        eligible = np.asarray([bool_value(value) for value in eligible_raw])
        spot_ids = read_h5ad_dataframe_column(handle["obs"], "spot_id").astype(str)
        genes = read_h5ad_dataframe_index(handle["var"]).astype(str)
    gene_index = {gene: index for index, gene in enumerate(genes)}
    score_table = pd.read_csv(
        output / "per_sample_program_scores.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    for column in ("primary_score", "equal_weight_score", "leave_top_gene_score"):
        score_table[column] = pd.to_numeric(score_table[column], errors="coerce")
    declared = pd.read_csv(
        output / "program_testability.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    declared_testable = declared["testable"].str.lower().isin(["true", "1"])
    program_ids = declared.loc[declared_testable, "program_uid"].iloc[:3].tolist()
    donor_ids = sorted(set(donors[eligible]))[:2]
    if not program_ids or not donor_ids:
        if len(score_table):
            raise ContractError(
                "Score table is nonempty despite a zero-testable or zero-eligible branch"
            )
        return 0, 0.0, len(donor_ids)
    max_error = 0.0
    n_comparisons = 0
    for donor in donor_ids:
        donor_index = np.flatnonzero((donors == donor) & eligible)
        for uid in program_ids:
            membership = registry.membership[
                (registry.membership["program_uid"] == uid)
                & registry.membership["mapped_symbol"].ne("")
                & registry.membership["mapped_symbol_status"].eq(TRUSTED_MAPPING_STATUS)
                & registry.membership["mapped_symbol"].isin(gene_index)
            ].copy()
            membership = membership.groupby("mapped_symbol", as_index=False)[
                "original_l1_weight"
            ].sum().sort_values("mapped_symbol")
            if len(membership) < MIN_PROGRAM_GENES:
                raise ContractError(f"Frozen score fixture program became untestable: {uid}")
            columns = np.asarray([gene_index[gene] for gene in membership["mapped_symbol"]])
            values = norm[donor_index][:, columns].toarray().astype(float)
            means = values.mean(axis=0)
            sd = values.std(axis=0, ddof=1)
            standardized = np.zeros_like(values)
            valid_sd = np.isfinite(sd) & (sd > 0)
            standardized[:, valid_sd] = (values[:, valid_sd] - means[valid_sd]) / sd[valid_sd]
            weights = membership["original_l1_weight"].to_numpy(float)
            primary = standardized @ (weights / weights.sum())
            equal = standardized.mean(axis=1)
            top = int(np.argmax(weights))
            keep = np.arange(len(weights)) != top
            leave = (
                standardized[:, keep] @ (weights[keep] / weights[keep].sum())
                if keep.sum() >= MIN_PROGRAM_GENES
                else np.repeat(np.nan, len(donor_index))
            )
            observed = score_table[
                (score_table["donor"] == donor) & (score_table["program_uid"] == uid)
            ].set_index("spot_id")
            ordered = observed.loc[spot_ids[donor_index]]
            for expected, column in (
                (primary, "primary_score"),
                (equal, "equal_weight_score"),
                (leave, "leave_top_gene_score"),
            ):
                observed_values = ordered[column].to_numpy(float)
                if np.isnan(expected).all() and np.isnan(observed_values).all():
                    error = 0.0
                elif not np.array_equal(np.isnan(expected), np.isnan(observed_values)):
                    error = np.inf
                else:
                    error = np.nanmax(np.abs(expected - observed_values))
                max_error = max(max_error, float(error))
                n_comparisons += 1
    return n_comparisons, max_error, len(donor_ids)


def independent_reference_model(output: Path) -> tuple[float, float, float]:
    design_table = pd.read_csv(
        output / "model_design_reference.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    draws = pd.read_csv(
        output / "bootstrap_reference_draws.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    if design_table.empty:
        return np.nan, np.nan, np.nan
    metadata = {
        "spot_id",
        "outcome_primary_score",
        "block_id",
        "section",
        "donor",
        "program_uid",
        "variant",
    }
    columns = [column for column in design_table.columns if column not in metadata]
    design = design_table[columns].apply(pd.to_numeric, errors="raise").to_numpy(float)
    outcome = pd.to_numeric(design_table["outcome_primary_score"], errors="raise").to_numpy(float)
    coefficient = np.linalg.lstsq(design, outcome, rcond=None)[0]
    residual = outcome - design @ coefficient
    residual_df = len(outcome) - design.shape[1]
    sigma2 = float(np.sum(residual**2) / residual_df)
    covariance = sigma2 * np.linalg.inv(design.T @ design)
    lipid_column = columns.index("lipid_zone")
    beta = float(coefficient[lipid_column])
    se = float(np.sqrt(covariance[lipid_column, lipid_column]))

    donor = design_table["donor"].iloc[0]
    uid = design_table["program_uid"].iloc[0]
    effects = pd.read_csv(
        output / "donor_program_effects.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    reference = effects[
        (effects["donor"] == donor)
        & (effects["program_uid"] == uid)
        & (effects["variant"] == "primary")
    ]
    if len(reference) != 1:
        raise ContractError("Reference donor-program model row is missing or duplicated")
    row = reference.iloc[0]
    beta_error = abs(beta - float(row["observed_beta"]))
    se_error = abs(se - float(row["conventional_ols_se"]))

    n_bootstrap = int(row["n_bootstrap"])
    seed = int(row["bootstrap_seed"])
    rng = np.random.default_rng(seed)
    blocks = pd.to_numeric(design_table["block_id"], errors="raise").to_numpy(int)
    sections = design_table["section"].astype(str).to_numpy()
    independent = np.full(n_bootstrap, np.nan)
    for draw in range(n_bootstrap):
        weights = np.zeros(len(outcome))
        for section in pd.unique(sections):
            section_blocks = np.unique(blocks[sections == section])
            counts = rng.multinomial(
                len(section_blocks), np.repeat(1 / len(section_blocks), len(section_blocks))
            )
            mask = sections == section
            for block, count in zip(section_blocks, counts, strict=True):
                if count:
                    weights[mask & (blocks == block)] = count
        positive = weights > 0
        if positive.sum() <= design.shape[1]:
            continue
        sqrt_w = np.sqrt(weights[positive])
        weighted_x = design[positive] * sqrt_w[:, None]
        if np.linalg.matrix_rank(weighted_x) != design.shape[1]:
            continue
        weighted_y = outcome[positive] * sqrt_w
        independent[draw] = np.linalg.lstsq(weighted_x, weighted_y, rcond=None)[0][
            lipid_column
        ]
    stored = pd.to_numeric(draws["lipid_zone_beta"], errors="coerce").to_numpy(float)
    if len(stored) != n_bootstrap:
        raise ContractError("Stored reference bootstrap draw count drift")
    bootstrap_error = float(np.nanmax(np.abs(independent - stored)))
    if not np.array_equal(np.isnan(independent), np.isnan(stored)):
        raise ContractError("Independent bootstrap failure pattern differs from stored pattern")
    return beta_error, se_error, bootstrap_error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    paths = default_paths(args.base)
    output = (args.output_dir or paths["candidate"]).resolve()
    terminal_outputs = [
        output / "validation_report.tsv",
        output / "execution_manifest.tsv",
        output / "candidate_release_manifest.tsv",
        output / "gate_status.tsv",
        output / "READY",
        output / "terminal_stage_seal.tsv",
    ]
    existing = [str(path) for path in terminal_outputs if path.exists()]
    if existing:
        raise ContractError("Refusing to overwrite terminal output(s): " + ", ".join(existing))
    required = [
        "analysis_freeze_manifest.tsv",
        "yakubovsky_human_adapter.h5ad",
        "sample_manifest.tsv",
        "gene_mapping_audit.tsv",
        "design_audit.tsv",
        "program_testability.tsv",
        "per_sample_program_scores.tsv",
        "per_donor_program_scores.tsv",
        "program_effects.tsv",
        "sensitivity.tsv",
        "multiplicity_manifest.tsv",
        "model_sample_manifest.tsv",
        "scoring_audit.tsv",
        "model_design_audit.tsv",
        "spatial_block_audit.tsv",
        "donor_program_effects.tsv",
        "donor_combination_audit.tsv",
        "zonation_reference.tsv",
        "model_design_reference.tsv",
        "bootstrap_reference_draws.tsv",
    ]
    missing = [name for name in required if not (output / name).is_file()]
    if missing:
        raise ContractError("Candidate bundle is incomplete: " + ", ".join(missing))
    validate_stage_chain(output, "model")

    checks: list[dict[str, Any]] = []
    registry = validate_registry_contract(paths)
    gate = validate_source_gate_summary(paths["source_gate"])
    add_check(
        checks,
        "VAL-REGISTRY-01",
        len(registry.family) == 2,
        len(registry.family),
        2,
        "Complete sealed external family is present",
    )
    add_check(
        checks,
        "VAL-SOURCE-01",
        len(gate.passing_donors) >= 3,
        len(gate.passing_donors),
        ">=3",
        "Plan 10 binary source gate remains valid",
    )

    freeze = pd.read_csv(
        output / "analysis_freeze_manifest.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    file_rows = freeze[freeze["record_type"] == "file_input"]
    drift = []
    for row in file_rows.itertuples(index=False):
        if sha256_file(Path(row.value)) != row.sha256:
            drift.append(row.name)
    add_check(
        checks,
        "VAL-FREEZE-01",
        not drift,
        ";".join(drift) or "none",
        "none",
        "Every pre-outcome frozen input hash rederives",
    )

    sample = pd.read_csv(
        output / "sample_manifest.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    add_check(
        checks,
        "VAL-BIOUNIT-01",
        not sample["spot_id"].duplicated().any(),
        sample["spot_id"].nunique(),
        len(sample),
        "Every retained spot maps once to donor/section/barcode",
    )
    observed_donors = sorted(
        sample.loc[
            sample["binary_lipid_inference_eligible"].str.lower().isin(["true", "1"]),
            "donor",
        ].unique()
    )
    add_check(
        checks,
        "VAL-BIOUNIT-02",
        observed_donors == sorted(gate.passing_donors),
        ";".join(observed_donors),
        ";".join(sorted(gate.passing_donors)),
        "Inference donors exactly match the Plan 10 gate",
    )

    adapter = output / "yakubovsky_human_adapter.h5ad"
    with h5py.File(adapter, "r") as handle:
        source = parse_h5ad_sparse(
            handle["layers/source_background_corrected_abundance"]
        )
        transformed = parse_h5ad_sparse(
            handle["layers/log1p_source_abundance_per_million"]
        )
        source_sparse = isinstance(
            handle["layers/source_background_corrected_abundance"], h5py.Group
        )
        source_valid = (
            (not len(source.data))
            or (np.isfinite(source.data).all() and (source.data >= 0).all())
        )
        x_absent = "X" not in handle
        h5ad_rows = source.shape[0]
        h5ad_columns = source.shape[1]
        authoritative_hepatocyte_eligibility = bool_value(
            h5_scalar(handle["uns/authoritative_hepatocyte_eligibility_available"])
        )
        zonation_landmark_interpretation = str(
            h5_scalar(handle["uns/source_zonation_landmark_interpretation"])
        )
        zonation_landmark_obs_present = (
            "source_zonation_landmark_expression" in handle["obs"]
        )
        prohibited_hepatocyte_proxy_present = (
            "source_hepatocyte_landmark_expression" in handle["obs"]
        )
        adapter_common_gene_axis_sha256 = str(
            h5_scalar(handle["uns/common_gene_axis_sha256"])
        )
        adapter_gene_axis_manifest_sha256 = str(
            h5_scalar(handle["uns/gene_axis_manifest_sha256"])
        )
        full_source_spot_universe_preserved = bool_value(
            h5_scalar(handle["uns/full_source_spot_universe_preserved"])
        )
        structural_join_required_for_analysis = bool_value(
            h5_scalar(handle["uns/structural_join_required_for_analysis"])
        )
    add_check(
        checks,
        "VAL-ADAPTER-01",
        source_sparse and source_valid and x_absent and transformed.shape == source.shape,
        (
            f"sparse={source_sparse};nonnegative_finite={source_valid};"
            f"X_absent={x_absent};transform_shape={transformed.shape}"
        ),
        "sparse=True;nonnegative_finite=True;X_absent=True;matching_transform_shape",
        "H5AD preserves source background-corrected abundance without a raw-count claim",
    )
    add_check(
        checks,
        "VAL-ADAPTER-02",
        h5ad_rows == len(sample) and h5ad_columns > 0,
        f"{h5ad_rows}x{h5ad_columns}",
        f"{len(sample)}x>0",
        "H5AD dimensions reconcile to the sample manifest",
    )
    eligibility_basis_ok = (
        "spot_eligibility_basis" in sample.columns
        and set(sample["spot_eligibility_basis"])
        == {
            "authoritative_structural_join_and_positive_source_abundance_no_hepatocyte_filter"
        }
    )
    structural_join = sample[
        "source_structural_expression_coordinate_zonation_join"
    ].map(bool_value)
    analysis_eligible = sample["analysis_eligible"].map(bool_value)
    positive_source = (
        pd.to_numeric(sample["background_corrected_spot_sum"], errors="raise") > 0
    )
    expected_analysis_eligible = structural_join & positive_source
    binary_eligible = sample["binary_lipid_inference_eligible"].map(bool_value)
    binary_gate_rederived = True
    for donor in gate.passing_donors:
        donor_mask = sample["donor"] == donor
        donor_binary = binary_eligible[donor_mask]
        donor_classes = set(
            sample.loc[donor_mask & binary_eligible, "source_defined_binary_lipid_class"]
        )
        binary_gate_rederived &= bool(
            donor_binary.sum() >= MIN_DONOR_SPOTS
            and donor_binary.mean() >= 0.90
            and set(ALLOWED_LIPID_CLASSES).issubset(donor_classes)
        )
    add_check(
        checks,
        "VAL-ADAPTER-02B",
        not authoritative_hepatocyte_eligibility
        and eligibility_basis_ok
        and zonation_landmark_obs_present
        and not prohibited_hepatocyte_proxy_present
        and full_source_spot_universe_preserved
        and structural_join_required_for_analysis
        and analysis_eligible.equals(expected_analysis_eligible)
        and structural_join[binary_eligible].all()
        and binary_gate_rederived
        and zonation_landmark_interpretation
        == "periportal_and_pericentral_zonation_landmark_expression_not_hepatocyte_identity_or_purity",
        (
            f"authoritative_hepatocyte_eligibility={authoritative_hepatocyte_eligibility};"
            f"basis_ok={eligibility_basis_ok};zonation_landmark_obs={zonation_landmark_obs_present};"
            f"prohibited_proxy={prohibited_hepatocyte_proxy_present};"
            f"full_source_universe={full_source_spot_universe_preserved};"
            f"eligibility_exact={analysis_eligible.equals(expected_analysis_eligible)};"
            f"binary_spot_structural_join={structural_join[binary_eligible].all()};"
            f"binary_gate_rederived={binary_gate_rederived};"
            f"interpretation={zonation_landmark_interpretation}"
        ),
        "full source universe;structural+positive analysis gate;binary spot joins and >=90% donor gate;zonation landmark only",
        "Reference-only spots remain represented but cannot enter scoring or lipid inference",
    )
    gene_axis_contract = validate_gene_axis_manifest(
        paths["gene_axis_manifest"],
        file_rows.loc[
            file_rows["name"] == "yakubovsky_v_mat", "sha256"
        ].iloc[0],
    )
    add_check(
        checks,
        "VAL-ADAPTER-02C",
        adapter_common_gene_axis_sha256 == gene_axis_contract.common_sha256
        and adapter_gene_axis_manifest_sha256 == gene_axis_contract.manifest_sha256,
        (
            f"axis={adapter_common_gene_axis_sha256};"
            f"manifest={adapter_gene_axis_manifest_sha256}"
        ),
        (
            f"axis={gene_axis_contract.common_sha256};"
            f"manifest={gene_axis_contract.manifest_sha256}"
        ),
        "Adapter reuse of one common gene axis is tied to the exact all-donor certificate",
    )
    adapter_design = pd.read_csv(
        output / "design_audit.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    optimized_chunks = pd.to_numeric(
        adapter_design["optimized_estimated_chunk_touches"], errors="raise"
    )
    physical_chunks = pd.to_numeric(
        adapter_design["n_physical_hdf5_chunks"], errors="raise"
    )
    legacy_gib = pd.to_numeric(
        adapter_design["legacy_256spot_estimated_logical_decompressed_gib"],
        errors="raise",
    )
    optimized_gib = pd.to_numeric(
        adapter_design["optimized_estimated_logical_decompressed_gib"],
        errors="raise",
    )
    aggregate_reduction = float(legacy_gib.sum() / optimized_gib.sum())
    add_check(
        checks,
        "VAL-ADAPTER-02D",
        set(adapter_design["matrix_traversal"])
        == {"hdf5_physical_chunk_aligned_exact_once"}
        and (optimized_chunks == physical_chunks).all()
        and set(adapter_design["physical_chunk_coverage_min"]) == {"1"}
        and set(adapter_design["physical_chunk_coverage_max"]) == {"1"}
        and aggregate_reduction > 10,
        (
            f"donors={len(adapter_design)};"
            f"optimized_touches={int(optimized_chunks.sum())};"
            f"physical_chunks={int(physical_chunks.sum())};"
            f"aggregate_estimated_reduction={aggregate_reduction:.6g}x"
        ),
        "all donor chunks covered exactly once;aggregate estimated reduction >10x",
        "Physical-chunk-aligned traversal and its estimated decompression reduction rederive",
    )
    fixture = source_abundance_fixture(paths, sample, output)
    add_check(
        checks,
        "VAL-ADAPTER-03",
        fixture["n_spot_fixtures"] == 10
        and fixture["n_gene_fixtures"] == 50
        and fixture["max_spot_sum_error"] <= 1e-10
        and fixture["max_abundance_error"] <= 1e-10
        and fixture["max_transform_error"] <= 1e-5,
        json.dumps(fixture, sort_keys=True),
        "10 spots;50 genes;source/sum<=1e-10;transform<=1e-5",
        "Fixed source values, spot sums, and log1p-per-million transform rederive",
    )
    add_check(
        checks,
        "VAL-ADAPTER-04",
        fixture["n_fractional_source_values"] > 0,
        fixture["n_fractional_source_values"],
        ">0",
        "Deposited background-corrected abundance demonstrates integer-count assertions are inapplicable",
    )

    comparisons, score_error, n_score_fixture_donors = independent_score_checks(
        paths, output, registry
    )
    declared_testable = pd.read_csv(
        output / "program_testability.tsv", sep="\t", dtype=str, keep_default_na=False
    )["testable"].str.lower().isin(["true", "1"]).sum()
    expected_comparisons = (
        min(3, int(declared_testable)) * n_score_fixture_donors * 3
    )
    add_check(
        checks,
        "VAL-SCORE-01",
        comparisons == expected_comparisons and score_error <= 1e-10,
        f"comparisons={comparisons};max_error={score_error:.3g}",
        f"comparisons={expected_comparisons};max_error<=1e-10",
        "Every testable frozen program (up to three) is hand-recomputed in up to two available donors; valid empty branches remain schema-safe",
    )

    blocks = pd.read_csv(
        output / "spatial_block_audit.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    model_samples = pd.read_csv(
        output / "model_sample_manifest.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    active_blocks = blocks[pd.to_numeric(blocks["block_id"], errors="coerce") >= 0]
    included_spots = model_samples[
        model_samples["included_spatial_block"].str.lower().isin(["true", "1"])
    ].copy()
    block_crosses = int(
        (
            included_spots.groupby("block_uid")[["donor", "section"]]
            .nunique()
            .max(axis=1)
            > 1
        ).sum()
    )
    audit_counts = pd.to_numeric(
        active_blocks.set_index("block_uid")["n_spots"], errors="raise"
    ).sort_index()
    spot_counts = included_spots.groupby("block_uid").size().sort_index()
    block_count_match = audit_counts.index.equals(spot_counts.index) and np.array_equal(
        audit_counts.to_numpy(int), spot_counts.to_numpy(int)
    )
    minimum_block_size = int(spot_counts.min()) if len(spot_counts) else 0
    add_check(
        checks,
        "VAL-BLOCK-01",
        block_crosses == 0
        and block_count_match
        and (not len(spot_counts) or minimum_block_size >= 5),
        (
            f"crosses={block_crosses};count_match={block_count_match};"
            f"minimum_spots={minimum_block_size}"
        ),
        "crosses=0;count_match=True;minimum_spots>=5 when any block is retained",
        "Spot-to-block membership reconciles and no block crosses donor/section boundaries",
    )

    beta_error, se_error, bootstrap_error = independent_reference_model(output)
    if np.isfinite(beta_error):
        add_check(
            checks,
            "VAL-MODEL-01",
            max(beta_error, se_error) <= 1e-10,
            f"beta_error={beta_error:.3g};se_error={se_error:.3g}",
            "<=1e-10",
            "Independent OLS refit recovers the exported reference model",
        )
        add_check(
            checks,
            "VAL-BOOTSTRAP-01",
            bootstrap_error <= 1e-8,
            bootstrap_error,
            "<=1e-8",
            "Independent whole-block replay recovers every reference bootstrap slope",
        )
    else:
        add_check(
            checks,
            "VAL-MODEL-01-SKIP",
            True,
            "no_estimable_reference_model",
            "valid_zonation_only_branch",
            "Reference refit is inapplicable when fewer than three spatially valid donors remain",
        )

    effects = pd.read_csv(
        output / "donor_program_effects.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    programs = pd.read_csv(
        output / "program_effects.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    programs["pvalue"] = pd.to_numeric(programs["pvalue"], errors="coerce")
    programs["qvalue"] = pd.to_numeric(programs["qvalue"], errors="coerce")
    programs["testable"] = programs["testable"].str.lower().isin(["true", "1"])
    add_check(
        checks,
        "VAL-FAMILY-01",
        len(programs) == len(registry.family)
        and set(programs["program_uid"]) == set(registry.family["program_uid"]),
        len(programs),
        len(registry.family),
        "All testable and untestable frozen-family members remain visible",
    )
    effects["estimable"] = effects["estimable"].str.lower().isin(["true", "1"])
    effects["observed_beta"] = pd.to_numeric(effects["observed_beta"], errors="coerce")
    effects["bootstrap_pvalue"] = pd.to_numeric(effects["bootstrap_pvalue"], errors="coerce")
    effects["n_complete_spots"] = pd.to_numeric(
        effects["n_complete_spots"], errors="coerce"
    )
    sensitivity = pd.read_csv(
        output / "sensitivity.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    for column in (
        "n_valid_donors",
        "combined_z",
        "combined_pvalue",
        "median_donor_slope",
        "combined_qvalue",
    ):
        sensitivity[column] = pd.to_numeric(sensitivity[column], errors="coerce")
    multiplicity = pd.read_csv(
        output / "multiplicity_manifest.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    combinations = pd.read_csv(
        output / "donor_combination_audit.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    for column in ("n_donors", "combined_z", "combined_pvalue"):
        combinations[column] = pd.to_numeric(combinations[column], errors="coerce")

    assay_testable_uids = set(programs.loc[programs["testable"], "program_uid"])
    multiplicity_errors: list[str] = []
    stouffer_errors: list[float] = []
    q_errors: list[float] = []
    complete_variants: set[str] = set()
    primary_expected_q_by_uid: dict[str, float] = {}
    for variant in sensitivity["variant"].unique():
        variant_rows = sensitivity[
            (sensitivity["variant"] == variant)
            & sensitivity["program_uid"].isin(assay_testable_uids)
        ].copy()
        estimable = (
            (variant_rows["n_valid_donors"] >= 3)
            & variant_rows["combined_pvalue"].notna()
        )
        expected_complete = bool(assay_testable_uids) and int(estimable.sum()) == len(
            assay_testable_uids
        )
        manifest_row = multiplicity[multiplicity["variant"] == variant]
        if len(manifest_row) != 1:
            multiplicity_errors.append(f"{variant}:missing_or_duplicate_manifest")
            continue
        manifest_row = manifest_row.iloc[0]
        observed_complete = bool_value(manifest_row["family_complete"])
        observed_applied = bool_value(manifest_row["bh_applied"])
        if (
            observed_complete != expected_complete
            or observed_applied != expected_complete
            or int(manifest_row["n_assay_testable_programs"]) != len(assay_testable_uids)
            or int(manifest_row["n_estimable_programs"]) != int(estimable.sum())
            or set(filter(None, manifest_row["assay_testable_program_uids"].split(";")))
            != assay_testable_uids
        ):
            multiplicity_errors.append(f"{variant}:family_manifest_drift")
        for row in variant_rows.loc[estimable].itertuples(index=False):
            donor = effects[
                (effects["program_uid"] == row.program_uid)
                & (effects["variant"] == variant)
                & effects["estimable"]
            ]
            z, pvalue = independent_stouffer(
                donor["observed_beta"].to_numpy(float),
                donor["bootstrap_pvalue"].to_numpy(float),
            )
            stouffer_errors.extend(
                [abs(z - float(row.combined_z)), abs(pvalue - float(row.combined_pvalue))]
            )
        if expected_complete:
            complete_variants.add(str(variant))
            expected_q = independent_bh(
                variant_rows.loc[estimable, "combined_pvalue"].to_numpy(float)
            )
            observed_q = variant_rows.loc[estimable, "combined_qvalue"].to_numpy(float)
            q_errors.extend(np.abs(expected_q - observed_q).tolist())
            if str(variant) == "primary":
                primary_expected_q_by_uid = {
                    str(uid): float(qvalue)
                    for uid, qvalue in zip(
                        variant_rows.loc[estimable, "program_uid"],
                        expected_q,
                        strict=True,
                    )
                }
        elif variant_rows["combined_qvalue"].notna().any():
            multiplicity_errors.append(f"{variant}:BH_present_for_incomplete_family")

    multiplicity_error = ";".join(multiplicity_errors) or "none"
    add_check(
        checks,
        "VAL-BH-01",
        not multiplicity_errors and (not q_errors or max(q_errors) <= 1e-12),
        f"manifest_errors={multiplicity_error};max_q_error={max(q_errors) if q_errors else 0}",
        "no manifest errors;max_q_error<=1e-12",
        "Every variant uses the complete assay-testable family or applies no BH",
    )

    primary_q_link_errors: list[str] = []
    primary_family_complete = "primary" in complete_variants
    if primary_family_complete:
        for uid in sorted(assay_testable_uids):
            program_row = programs[programs["program_uid"] == uid]
            sensitivity_row = sensitivity[
                (sensitivity["program_uid"] == uid)
                & (sensitivity["variant"] == "primary")
            ]
            if (
                len(program_row) != 1
                or len(sensitivity_row) != 1
                or uid not in primary_expected_q_by_uid
            ):
                primary_q_link_errors.append(f"{uid}:missing_primary_q_row")
                continue
            expected_q = primary_expected_q_by_uid[uid]
            program_q = float(program_row.iloc[0]["qvalue"])
            sensitivity_q = float(sensitivity_row.iloc[0]["combined_qvalue"])
            if (
                not np.isfinite(program_q)
                or not np.isfinite(sensitivity_q)
                or abs(program_q - expected_q) > 1e-12
                or abs(sensitivity_q - expected_q) > 1e-12
                or abs(program_q - sensitivity_q) > 1e-12
            ):
                primary_q_link_errors.append(f"{uid}:primary_q_mismatch")
        non_testable_q = programs.loc[
            ~programs["program_uid"].isin(assay_testable_uids), "qvalue"
        ]
        if non_testable_q.notna().any():
            primary_q_link_errors.append("untestable_program_has_primary_q")
    elif programs["qvalue"].notna().any():
        primary_q_link_errors.append("primary_q_present_for_incomplete_family")
    add_check(
        checks,
        "VAL-BH-PRIMARY-02",
        not primary_q_link_errors,
        ";".join(primary_q_link_errors) or "none",
        "none",
        "program_effects primary q exactly equals the independently rederived primary BH q and the corresponding sensitivity record",
    )
    stouffer_error = max(stouffer_errors) if stouffer_errors else 0.0
    add_check(
        checks,
        "VAL-STOUFFER-01",
        stouffer_error <= 1e-12,
        stouffer_error,
        "<=1e-12",
        "Every estimable primary and sensitivity signed-Stouffer result rederives",
    )

    registry_direction = registry.family.set_index("program_uid")[
        "expected_direction"
    ].to_dict()
    robust_errors: list[str] = []
    spot_count_errors: list[str] = []
    program_robust_bool = programs["robust"].str.lower().isin(["true", "1"])
    for row_index, row in programs.iterrows():
        uid = row["program_uid"]
        expected_direction = registry_direction[uid]
        primary_donor = effects[
            (effects["program_uid"] == uid)
            & (effects["variant"] == "primary")
            & effects["estimable"]
        ].copy()
        for donor_row in effects[
            (effects["program_uid"] == uid) & (effects["variant"] == "primary")
        ].itertuples(index=False):
            expected_spots = int(
                (
                    model_samples["primary_model_eligible"].str.lower().isin(["true", "1"])
                    & (model_samples["donor"] == donor_row.donor)
                ).sum()
            )
            if int(donor_row.n_complete_spots) != expected_spots:
                spot_count_errors.append(f"{uid}:{donor_row.donor}")

        if not bool(row["testable"]):
            expected_robust = False
            expected_state = "untestable"
        elif not primary_family_complete:
            expected_robust = False
            expected_state = "untestable"
        else:
            combined_z, combined_p = independent_stouffer(
                primary_donor["observed_beta"].to_numpy(float),
                primary_donor["bootstrap_pvalue"].to_numpy(float),
            )
            combined_direction = direction_label(combined_z)
            median_direction = direction_label(float(row["estimate"]))
            donor_direction_count = int(
                (
                    np.sign(primary_donor["observed_beta"].to_numpy(float))
                    == expected_sign(expected_direction)
                ).sum()
            )
            required_donor_directions = math.ceil(0.75 * len(gate.passing_donors))
            sensitivity_ok = True
            for variant in ("equal_weight", "leave_top_gene", "zonation_df3", "zonation_df5"):
                variant_row = sensitivity[
                    (sensitivity["program_uid"] == uid)
                    & (sensitivity["variant"] == variant)
                ]
                sensitivity_ok &= (
                    variant in complete_variants
                    and len(variant_row) == 1
                    and direction_label(float(variant_row.iloc[0]["combined_z"]))
                    == expected_direction
                )
            lodo = combinations[
                (combinations["program_uid"] == uid)
                & (combinations["combination"] == "leave_one_donor_out")
            ]
            lodo_ok = len(lodo) == len(primary_donor)
            for donor in primary_donor["donor"].astype(str):
                retained = primary_donor[primary_donor["donor"].astype(str) != donor]
                lodo_z, lodo_p = independent_stouffer(
                    retained["observed_beta"].to_numpy(float),
                    retained["bootstrap_pvalue"].to_numpy(float),
                )
                observed_lodo = lodo[lodo["omitted_donor"] == donor]
                lodo_ok &= (
                    len(observed_lodo) == 1
                    and abs(float(observed_lodo.iloc[0]["combined_z"]) - lodo_z) <= 1e-12
                    and abs(float(observed_lodo.iloc[0]["combined_pvalue"]) - lodo_p)
                    <= 1e-12
                    and direction_label(lodo_z) == expected_direction
                )
            expected_robust = bool(
                primary_expected_q_by_uid[uid] < 0.05
                and combined_direction == expected_direction
                and median_direction == combined_direction
                and donor_direction_count >= required_donor_directions
                and sensitivity_ok
                and lodo_ok
            )
            expected_state = "robust" if expected_robust else "tested_negative"
            if abs(float(row["combined_z"]) - combined_z) > 1e-12 or abs(
                float(row["pvalue"]) - combined_p
            ) > 1e-12:
                robust_errors.append(f"{uid}:primary_combination")
        if bool(program_robust_bool.iloc[row_index]) != expected_robust:
            robust_errors.append(f"{uid}:robust_flag")
        if row["evidence_state"] != expected_state:
            robust_errors.append(f"{uid}:state={row['evidence_state']}!={expected_state}")
    add_check(
        checks,
        "VAL-BIOUNIT-03",
        not spot_count_errors,
        ";".join(spot_count_errors) or "none",
        "none",
        "Every donor-program complete-spot count matches the exported spot inclusion map",
    )
    add_check(
        checks,
        "VAL-ROBUST-01",
        not robust_errors,
        ";".join(robust_errors) or "none",
        "none",
        f"Final robust/state calls independently rederive under: {ROBUST_RULE}",
    )

    stage_manifests = [
        output / "freeze_execution_manifest.tsv",
        output / "adapter_execution_manifest.tsv",
        output / "scoring_execution_manifest.tsv",
        output / "model_execution_manifest.tsv",
    ]
    stage_rows = [
        pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        for path in stage_manifests
    ]
    combined_execution = pd.concat(stage_rows, ignore_index=True, sort=False)
    continuous_read = any(
        "continuous_lipid_fields_read" in frame.columns
        and frame["continuous_lipid_fields_read"].str.lower().isin(["true", "1"]).any()
        for frame in stage_rows
    )
    cell2location_used = False
    count_likelihood_used = False
    adapter_execution = stage_rows[1]
    if "cell2location_used" in adapter_execution.columns:
        cell2location_used = adapter_execution["cell2location_used"].str.lower().isin(
            ["true", "1"]
        ).any()
    if "count_likelihood_used" in adapter_execution.columns:
        count_likelihood_used = adapter_execution["count_likelihood_used"].str.lower().isin(
            ["true", "1"]
        ).any()
    add_check(
        checks,
        "VAL-SCOPE-01",
        not continuous_read and not cell2location_used and not count_likelihood_used,
        (
            f"continuous_read={continuous_read};cell2location={cell2location_used};"
            f"count_likelihood={count_likelihood_used}"
        ),
        "False;False;False",
        "Binary-only adapter used no continuous lipid field, cell2location, or count likelihood",
    )
    n_v1 = verify_v1(paths)
    add_check(
        checks,
        "VAL-V1-01",
        n_v1 == 44,
        n_v1,
        44,
        "Every protected v1 file remains byte-identical",
    )

    validation = pd.DataFrame(checks)
    validation.insert(0, "release_id", RELEASE_ID)
    validation["validated_at_utc"] = utc_now()
    validation["validator"] = str(Path(__file__).resolve())
    validation["validator_sha256"] = sha256_file(Path(__file__).resolve())
    atomic_write_frame(output / "validation_report.tsv", validation)
    combined_execution = pd.concat(
        [
            combined_execution,
            pd.DataFrame(
                [
                    {
                        "release_id": RELEASE_ID,
                        "stage": "independent_validation",
                        "producer": str(Path(__file__).resolve()),
                        "producer_sha256": sha256_file(Path(__file__).resolve()),
                        "python": sys.version.replace("\n", " "),
                        "platform": platform.platform(),
                        "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
                        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
                        "completed_utc": utc_now(),
                        "n_checks": len(validation),
                        "n_failed_checks": 0,
                        "continuous_lipid_fields_read": False,
                        "exit_state": "pass",
                    }
                ]
            ),
        ],
        ignore_index=True,
        sort=False,
    )
    atomic_write_frame(output / "execution_manifest.tsv", combined_execution)

    manifest_rows = []
    terminal_names = {
        "candidate_release_manifest.tsv",
        "gate_status.tsv",
        "READY",
    }
    for path in sorted(output.iterdir()):
        if not path.is_file() or path.name in terminal_names or path.name.startswith("."):
            continue
        manifest_rows.append(
            {
                "release_id": RELEASE_ID,
                "role": "candidate_output",
                "relative_path": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    script_dir = Path(__file__).resolve().parent
    for path in sorted(script_dir.iterdir()):
        if path.is_file() and path.suffix in {".py", ".sh", ".sbatch", ".md"}:
            manifest_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "role": "producer_code",
                    "relative_path": str(path.relative_to(paths["base"])),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    atomic_write_frame(output / "candidate_release_manifest.tsv", pd.DataFrame(manifest_rows))

    model_execution = stage_rows[3].iloc[0]
    valid_block_donors = int(model_execution.get("n_valid_block_donors", "0"))
    minimum_valid_primary_donors = int(
        model_execution.get(
            "minimum_valid_primary_donors_across_testable_programs", "0"
        )
    )
    primary_family_complete = bool_value(
        model_execution.get("primary_prespecified_BH_family_complete", "False")
    )
    terminal_status = (
        "complete_binary_lipid_analysis"
        if minimum_valid_primary_donors >= 3 and primary_family_complete
        else "complete_zonation_reference_spatial_model_unestimable"
    )
    gate_row = pd.DataFrame(
        [
            {
                "release_id": RELEASE_ID,
                "dataset": "yakubovsky2026",
                "status": terminal_status,
                "source_gate_machine_label": gate.row["terminal_verdict"],
                "analysis_exposure_type": "binary_lipid_zone_vs_non_lipid_zone",
                "non_lipid_zone_is_lipid_free": False,
                "n_source_gate_donors": len(gate.passing_donors),
                "n_valid_spatial_model_donors": valid_block_donors,
                "minimum_valid_primary_donors_across_testable_programs": minimum_valid_primary_donors,
                "primary_prespecified_BH_family_complete": primary_family_complete,
                "n_frozen_family": len(registry.family),
                "n_testable_programs": int(programs["testable"].sum()),
                "n_robust_programs": int(
                    programs["robust"].str.lower().isin(["true", "1"]).sum()
                ),
                "biological_positivity_required": False,
                "continuous_lipid_fields_read": False,
                "validation_report_sha256": sha256_file(output / "validation_report.tsv"),
                "candidate_release_manifest_sha256": sha256_file(
                    output / "candidate_release_manifest.tsv"
                ),
                "registry_sha256": registry.registry_sha256,
                "membership_sha256": registry.membership_sha256,
                "completed_utc": utc_now(),
            }
        ]
    )
    atomic_write_frame(output / "gate_status.tsv", gate_row)
    ready = gate_row.copy()
    ready["status"] = "ready_for_plan13_candidate_integration"
    ready["gate_status_sha256"] = sha256_file(output / "gate_status.tsv")
    ready["canonical_promotion_authorized"] = False
    atomic_write_frame(output / "READY", ready)
    write_stage_seal(output, "terminal", TERMINAL_STAGE_ARTIFACTS)
    print(
        json.dumps(
            {
                "release_id": RELEASE_ID,
                "status": terminal_status,
                "n_validation_checks": len(validation),
                "n_valid_spatial_model_donors": valid_block_donors,
                "n_robust_programs": int(
                    programs["robust"].str.lower().isin(["true", "1"]).sum()
                ),
                "ready_sha256": sha256_file(output / "READY"),
                "canonical_promotion_authorized": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
