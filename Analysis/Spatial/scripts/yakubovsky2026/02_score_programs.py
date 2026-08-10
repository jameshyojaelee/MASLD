#!/usr/bin/env python3
"""Score the sealed hepatocyte family without reading binary lipid labels."""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats

from yakubovsky_common import (
    DATASET,
    MIN_PROGRAM_GENES,
    MIN_RETAINED_L1,
    RELEASE_ID,
    SCORING_STAGE_ARTIFACTS,
    ContractError,
    atomic_write_frame,
    bool_value,
    default_paths,
    parse_h5ad_sparse,
    read_h5ad_dataframe_column,
    read_h5ad_dataframe_index,
    score_frozen_program,
    sha256_file,
    standardize_columns,
    validate_registry_contract,
    validate_stage_chain,
    write_stage_seal,
)


SAFE_OBS_FIELDS = (
    "_index",
    "spot_id",
    "donor",
    "section",
    "barcode",
    "zonation_eta",
    "background_corrected_spot_sum",
    "detected_genes",
    "source_zonation_landmark_expression",
    "analysis_eligible",
)
FORBIDDEN_OUTCOME_TOKENS = ("lipid", "class", "zone_membership")
TRUSTED_MAPPING_STATUS = "gencode_v49_unique_symbol_confirmed"

TESTABILITY_COLUMNS = (
    "release_id",
    "dataset",
    "program_uid",
    "legacy_module",
    "program_name",
    "expected_direction",
    "n_original_membership_rows",
    "n_gencode_v49_unique_membership_rows",
    "n_mapped_genes",
    "source_membership_l1_weight_sum",
    "retained_original_l1_weight",
    "testable",
    "testability_reason",
    "leave_top_gene_testable",
    "scoring_transform",
    "registry_sha256",
    "membership_sha256",
)
SCORE_COLUMNS = (
    "release_id",
    "dataset",
    "spot_id",
    "donor",
    "section",
    "barcode",
    "program_uid",
    "legacy_module",
    "primary_score",
    "equal_weight_score",
    "leave_top_gene_score",
    "zonation_eta",
    "background_corrected_spot_sum",
    "detected_genes",
    "source_zonation_landmark_expression",
)
DONOR_SCORE_COLUMNS = (
    "release_id",
    "dataset",
    "donor",
    "program_uid",
    "score_definition",
    "n_spots",
    "mean_score",
    "sd_score",
    "median_score",
)
ZONATION_REFERENCE_COLUMNS = (
    "release_id",
    "dataset",
    "donor",
    "program_uid",
    "n_spots",
    "spearman_rho",
    "spearman_p_descriptive",
    "natural_spline_df",
    "spline_r_squared",
    "interpretation",
)
SCORING_AUDIT_COLUMNS = (
    "audit_stage",
    "dataset",
    "donor",
    "program_uid",
    "n_spots",
    "n_mapped_genes",
    "n_variable_genes_within_donor",
    "top_original_weight_gene",
    "top_original_l1_weight",
    "primary_weight_sum_after_renormalization",
    "lipid_label_read",
)


def assemble_scoring_output_frames(
    score_rows: list[dict[str, Any]],
    donor_summary_rows: list[dict[str, Any]],
    zonation_rows: list[dict[str, Any]],
    scoring_audit_rows: list[dict[str, Any]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Materialize fixed schemas, including every valid empty branch."""
    return (
        pd.DataFrame(score_rows, columns=SCORE_COLUMNS),
        pd.DataFrame(donor_summary_rows, columns=DONOR_SCORE_COLUMNS),
        pd.DataFrame(zonation_rows, columns=ZONATION_REFERENCE_COLUMNS),
        pd.DataFrame(scoring_audit_rows, columns=SCORING_AUDIT_COLUMNS),
    )


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def decode_scalar(dataset: h5py.Dataset) -> Any:
    value = dataset[()]
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray) and value.shape == ():
        value = value.item()
        return value.decode() if isinstance(value, bytes) else value
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    paths = default_paths(args.base)
    output = (args.output_dir or paths["candidate"]).resolve()
    adapter = output / "yakubovsky_human_adapter.h5ad"
    if not adapter.is_file():
        raise ContractError(f"Missing source-native adapter: {adapter}")
    required_freeze = [
        output / "analysis_freeze_manifest.tsv",
        output / "registry_manifest.tsv",
        output / "adapter_execution_manifest.tsv",
    ]
    if any(not path.is_file() for path in required_freeze):
        raise ContractError("Adapter/scoring prerequisites are incomplete")
    validate_stage_chain(output, "adapter")
    final_outputs = [
        output / "program_testability.tsv",
        output / "per_sample_program_scores.tsv",
        output / "per_donor_program_scores.tsv",
        output / "zonation_reference.tsv",
        output / "scoring_audit.tsv",
        output / "scoring_execution_manifest.tsv",
    ]
    existing = [str(path) for path in final_outputs if path.exists()]
    if existing:
        raise ContractError("Refusing to overwrite scoring output(s): " + ", ".join(existing))

    registry = validate_registry_contract(paths)
    family = registry.family.set_index("program_uid", drop=False)
    membership = registry.membership.copy()
    # Preserve every frozen membership row in the registry provenance, but only
    # score symbols that Plan 20 positively resolved against GENCODE v49.  A
    # non-empty convenience symbol is not, by itself, evidence of an
    # unambiguous mapping.
    membership = membership[
        membership["mapped_symbol"].ne("")
        & membership["mapped_symbol_status"].eq(TRUSTED_MAPPING_STATUS)
    ].copy()
    membership["original_l1_weight"] = pd.to_numeric(
        membership["original_l1_weight"], errors="raise"
    )

    accessed_h5ad_paths: list[str] = []
    with h5py.File(adapter, "r") as handle:
        if not isinstance(
            handle["layers/log1p_source_abundance_per_million"], h5py.Group
        ):
            raise ContractError("Adapter log1p-per-million abundance layer is not sparse")
        normalized = parse_h5ad_sparse(
            handle["layers/log1p_source_abundance_per_million"]
        )
        accessed_h5ad_paths.append("/layers/log1p_source_abundance_per_million")
        obs_group = handle["obs"]
        obs_values: dict[str, np.ndarray] = {}
        for field in SAFE_OBS_FIELDS:
            if field not in obs_group:
                if field == "_index" and "spot_id" in obs_group:
                    continue
                raise ContractError(f"Adapter lacks safe scoring field obs/{field}")
            obs_values[field] = read_h5ad_dataframe_column(obs_group, field)
            accessed_h5ad_paths.append(f"/obs/{field}")
        var_names = read_h5ad_dataframe_index(handle["var"]).astype(str)
        var_index_name = handle["var"].attrs.get("_index", "_index")
        if isinstance(var_index_name, bytes):
            var_index_name = var_index_name.decode()
        accessed_h5ad_paths.append(f"/var/{var_index_name}")
        if "uns" in handle:
            for field in ("registry_sha256", "membership_sha256", "continuous_lipid_fields_read"):
                if field not in handle["uns"]:
                    raise ContractError(f"Adapter lacks provenance field uns/{field}")
            adapter_registry = str(decode_scalar(handle["uns/registry_sha256"]))
            adapter_membership = str(decode_scalar(handle["uns/membership_sha256"]))
            continuous_read = decode_scalar(handle["uns/continuous_lipid_fields_read"])
            accessed_h5ad_paths.extend(
                [
                    "/uns/registry_sha256",
                    "/uns/membership_sha256",
                    "/uns/continuous_lipid_fields_read",
                ]
            )
        else:
            raise ContractError("Adapter lacks provenance metadata")

    if adapter_registry != registry.registry_sha256 or adapter_membership != registry.membership_sha256:
        raise ContractError("Adapter registry provenance does not match the frozen family")
    if bool_value(continuous_read):
        raise ContractError("Adapter reports that continuous lipid fields were read")
    if any(
        path.startswith("/obs/")
        and any(token in path.lower() for token in FORBIDDEN_OUTCOME_TOKENS)
        for path in accessed_h5ad_paths
    ):
        raise ContractError("Scoring accessed a forbidden lipid-outcome H5AD path")
    if normalized.shape != (len(obs_values["donor"]), len(var_names)):
        raise ContractError("Adapter normalized layer dimensions do not match obs/var")

    obs = pd.DataFrame(
        {
            "spot_id": obs_values["spot_id"].astype(str),
            "donor": obs_values["donor"].astype(str),
            "section": obs_values["section"].astype(str),
            "barcode": obs_values["barcode"].astype(str),
            "zonation_eta": pd.to_numeric(obs_values["zonation_eta"], errors="raise"),
            "background_corrected_spot_sum": pd.to_numeric(
                obs_values["background_corrected_spot_sum"], errors="raise"
            ),
            "detected_genes": pd.to_numeric(obs_values["detected_genes"], errors="raise"),
            "source_zonation_landmark_expression": pd.to_numeric(
                obs_values["source_zonation_landmark_expression"], errors="coerce"
            ),
            "analysis_eligible": [
                str(value).strip().lower() in {"true", "1", "yes"}
                if not isinstance(value, (bool, np.bool_))
                else bool(value)
                for value in obs_values["analysis_eligible"]
            ],
        }
    )
    if obs["spot_id"].duplicated().any():
        raise ContractError("Scoring input contains duplicate spot IDs")
    gene_to_index = {gene: index for index, gene in enumerate(var_names)}
    if len(gene_to_index) != len(var_names):
        raise ContractError("Adapter var names are not unique")

    testability_rows: list[dict[str, Any]] = []
    program_definitions: dict[str, pd.DataFrame] = {}
    for uid, family_row in family.iterrows():
        source_part = registry.membership[
            registry.membership["program_uid"] == uid
        ].copy()
        part = membership[membership["program_uid"] == uid].copy()
        part = (
            part.groupby("mapped_symbol", as_index=False, sort=True)["original_l1_weight"]
            .sum()
            .rename(columns={"mapped_symbol": "gene_symbol"})
        )
        source_weight_total = float(source_part["original_l1_weight"].astype(float).sum())
        part["measured"] = part["gene_symbol"].isin(gene_to_index)
        measured = part[part["measured"]].copy()
        retained = float(measured["original_l1_weight"].sum())
        testable = len(measured) >= MIN_PROGRAM_GENES and retained >= MIN_RETAINED_L1
        leave_top_testable = testable and len(measured) - 1 >= MIN_PROGRAM_GENES
        reason = (
            "testable"
            if testable
            else "fewer_than_8_mapped_genes"
            if len(measured) < MIN_PROGRAM_GENES
            else "retained_original_l1_weight_below_0.20"
        )
        testability_rows.append(
            {
                "release_id": RELEASE_ID,
                "dataset": DATASET,
                "program_uid": uid,
                "legacy_module": family_row["module"],
                "program_name": family_row["module_name"],
                "expected_direction": family_row["expected_direction"],
                "n_original_membership_rows": len(source_part),
                "n_gencode_v49_unique_membership_rows": len(part),
                "n_mapped_genes": len(measured),
                "source_membership_l1_weight_sum": source_weight_total,
                "retained_original_l1_weight": retained,
                "testable": testable,
                "testability_reason": reason,
                "leave_top_gene_testable": leave_top_testable,
                "scoring_transform": "within_donor_gene_z_on_log1p_source_abundance_per_million",
                "registry_sha256": registry.registry_sha256,
                "membership_sha256": family_row["membership_sha256"],
            }
        )
        if testable:
            program_definitions[uid] = measured.sort_values("gene_symbol").reset_index(drop=True)

    testability = pd.DataFrame(testability_rows, columns=TESTABILITY_COLUMNS)
    if len(testability) != len(family) or set(testability["program_uid"]) != set(family.index):
        raise ContractError("Program testability does not contain the complete frozen family")

    score_rows: list[dict[str, Any]] = []
    donor_summary_rows: list[dict[str, Any]] = []
    scoring_audit_rows: list[dict[str, Any]] = []
    zonation_rows: list[dict[str, Any]] = []
    for donor in pd.unique(obs["donor"]):
        donor_idx = np.flatnonzero((obs["donor"] == donor) & obs["analysis_eligible"])
        if not len(donor_idx):
            continue
        for uid, program in program_definitions.items():
            indices = np.asarray([gene_to_index[gene] for gene in program["gene_symbol"]], dtype=int)
            values = normalized[donor_idx][:, indices].toarray().astype(float)
            standardized = standardize_columns(values)
            variable = np.nanstd(values, axis=0, ddof=1) > 0
            weights = program["original_l1_weight"].to_numpy(float)
            primary, equal, leave, top = score_frozen_program(standardized, weights)
            family_row = family.loc[uid]
            for local_index, global_index in enumerate(donor_idx):
                score_rows.append(
                    {
                        "release_id": RELEASE_ID,
                        "dataset": DATASET,
                        "spot_id": obs.iloc[global_index]["spot_id"],
                        "donor": donor,
                        "section": obs.iloc[global_index]["section"],
                        "barcode": obs.iloc[global_index]["barcode"],
                        "program_uid": uid,
                        "legacy_module": family_row["module"],
                        "primary_score": primary[local_index],
                        "equal_weight_score": equal[local_index],
                        "leave_top_gene_score": (
                            leave[local_index] if leave is not None else np.nan
                        ),
                        "zonation_eta": obs.iloc[global_index]["zonation_eta"],
                        "background_corrected_spot_sum": obs.iloc[global_index][
                            "background_corrected_spot_sum"
                        ],
                        "detected_genes": obs.iloc[global_index]["detected_genes"],
                        "source_zonation_landmark_expression": obs.iloc[global_index][
                            "source_zonation_landmark_expression"
                        ],
                    }
                )
            for score_name, values_out in (
                ("primary", primary),
                ("equal_weight", equal),
                ("leave_top_gene", leave),
            ):
                if values_out is None:
                    continue
                donor_summary_rows.append(
                    {
                        "release_id": RELEASE_ID,
                        "dataset": DATASET,
                        "donor": donor,
                        "program_uid": uid,
                        "score_definition": score_name,
                        "n_spots": len(values_out),
                        "mean_score": float(np.mean(values_out)),
                        "sd_score": float(np.std(values_out, ddof=1)),
                        "median_score": float(np.median(values_out)),
                    }
                )
            zonation = obs.iloc[donor_idx]["zonation_eta"].to_numpy(float)
            rho, rho_p = stats.spearmanr(primary, zonation)
            spline = np.asarray(
                __import__("patsy").dmatrix(
                    "cr(z, df=4, constraints='center')",
                    {"z": zonation},
                    return_type="dataframe",
                )
            )
            fitted = spline @ np.linalg.lstsq(spline, primary, rcond=None)[0]
            total_ss = float(np.sum((primary - np.mean(primary)) ** 2))
            residual_ss = float(np.sum((primary - fitted) ** 2))
            zonation_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "dataset": DATASET,
                    "donor": donor,
                    "program_uid": uid,
                    "n_spots": len(primary),
                    "spearman_rho": rho,
                    "spearman_p_descriptive": rho_p,
                    "natural_spline_df": 4,
                    "spline_r_squared": 1 - residual_ss / total_ss if total_ss > 0 else np.nan,
                    "interpretation": "descriptive_normal_zonation_reference_not_lipid_inference",
                }
            )
            scoring_audit_rows.append(
                {
                    "audit_stage": "outcome_blind_scoring",
                    "dataset": DATASET,
                    "donor": donor,
                    "program_uid": uid,
                    "n_spots": len(donor_idx),
                    "n_mapped_genes": len(program),
                    "n_variable_genes_within_donor": int(variable.sum()),
                    "top_original_weight_gene": program.iloc[top]["gene_symbol"],
                    "top_original_l1_weight": weights[top],
                    "primary_weight_sum_after_renormalization": 1.0,
                    "lipid_label_read": False,
                }
            )

    scores, donor_scores, zonation_reference, scoring_audit = (
        assemble_scoring_output_frames(
            score_rows,
            donor_summary_rows,
            zonation_rows,
            scoring_audit_rows,
        )
    )
    expected_score_rows = int(
        sum(obs["analysis_eligible"]) * testability["testable"].astype(bool).sum()
    )
    if len(scores) != expected_score_rows:
        raise ContractError(
            f"Score row count mismatch: observed {len(scores)}, expected {expected_score_rows}"
        )
    if scores.duplicated(["spot_id", "program_uid"]).any():
        raise ContractError("Per-spot program scores are not unique")
    if any("lipid" in column.lower() for column in scores.columns):
        raise ContractError("Outcome-blind score table contains a lipid label column")

    atomic_write_frame(output / "program_testability.tsv", testability)
    atomic_write_frame(output / "per_sample_program_scores.tsv", scores)
    atomic_write_frame(output / "per_donor_program_scores.tsv", donor_scores)
    atomic_write_frame(output / "zonation_reference.tsv", zonation_reference)

    atomic_write_frame(output / "scoring_audit.tsv", scoring_audit)
    execution = pd.DataFrame(
        [
            {
                "release_id": RELEASE_ID,
                "stage": "outcome_blind_program_scoring",
                "producer": str(Path(__file__).resolve()),
                "producer_sha256": sha256_file(Path(__file__).resolve()),
                "python": sys.version.replace("\n", " "),
                "platform": platform.platform(),
                "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
                "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
                "completed_utc": utc_now(),
                "n_frozen_family": len(family),
                "n_testable_programs": int(testability["testable"].astype(bool).sum()),
                "n_score_rows": len(scores),
                "h5ad_paths_accessed": ";".join(accessed_h5ad_paths),
                "binary_lipid_labels_read": False,
                "continuous_lipid_fields_read": False,
                "external_outcomes_used_for_selection": False,
                "exit_state": "pass",
            }
        ]
    )
    atomic_write_frame(output / "scoring_execution_manifest.tsv", execution)
    write_stage_seal(output, "scoring", SCORING_STAGE_ARTIFACTS)
    print(
        json.dumps(
            {
                "release_id": RELEASE_ID,
                "status": "outcome_blind_scoring_complete",
                "n_frozen_family": len(family),
                "n_testable_programs": int(testability["testable"].astype(bool).sum()),
                "n_score_rows": len(scores),
                "binary_lipid_labels_read": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
