#!/usr/bin/env python3
"""Independent numerical and release-contract validation for the candidate."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from common import (
    MEMBERSHIP,
    READY,
    REGISTRY,
    bh_adjust,
    load_donor_metadata,
    refuse_existing,
    require,
    sha256,
)


CANDIDATE = Path(os.environ["CAND_ROOT"])
RESULTS = CANDIDATE / "results"
CROSS_CANDIDATE = Path(os.environ["CROSS_CAND_ROOT"])
CROSS_RESULTS = CROSS_CANDIDATE / "results"


def independent_fit(
    frame: pd.DataFrame, value_column: str
) -> dict[str, float | int]:
    data = frame[
        (~frame["exclude"])
        & frame["stage_ordinal"].notna()
        & frame[value_column].notna()
    ].copy()
    if len(data) < 20 or data["dataset"].nunique() < 2:
        raise RuntimeError("independent model is not estimable")
    dummies = pd.get_dummies(data["dataset"], drop_first=True, dtype=float)
    design = np.column_stack(
        [np.ones(len(data)), data["stage_ordinal"].to_numpy(float), dummies.to_numpy(float)]
    )
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise RuntimeError("independent model is rank deficient")
    y = data[value_column].to_numpy(float)
    inverse = np.linalg.inv(design.T @ design)
    coefficients = inverse @ design.T @ y
    residuals = y - design @ coefficients
    n, k = design.shape
    df = n - k
    variance = float(residuals @ residuals) / df * inverse
    se = float(np.sqrt(variance[1, 1]))
    leverage = np.clip(np.diag(design @ inverse @ design.T), 0, 1 - 1e-12)
    omega = (residuals / (1 - leverage)) ** 2
    hc3 = inverse @ (design.T @ (omega[:, None] * design)) @ inverse
    hc3_se = float(np.sqrt(hc3[1, 1]))
    beta = float(coefficients[1])
    return {
        "beta": beta,
        "se": se,
        "pvalue": float(2 * stats.t.sf(abs(beta / se), df)) if se else 1.0,
        "hc3_se": hc3_se,
        "hc3_pvalue": float(2 * stats.t.sf(abs(beta / hc3_se), df)) if hc3_se else 1.0,
        "n_donors": int(data["donor"].nunique()),
        "n_datasets": int(data["dataset"].nunique()),
    }


def independent_paired_fit(
    scores: pd.DataFrame, uid: str, comparison_lineage: str
) -> dict[str, float | int]:
    subset = scores[scores["lineage"].isin(["Hepatocytes", comparison_lineage])]
    wide = subset.pivot_table(
        index=["dataset", "donor"], columns="lineage", values=uid, aggfunc="first"
    ).reset_index()
    if "Hepatocytes" not in wide or comparison_lineage not in wide:
        raise RuntimeError("independent paired columns unavailable")
    wide = wide.dropna(subset=["Hepatocytes", comparison_lineage]).copy()
    wide["difference"] = wide[comparison_lineage] - wide["Hepatocytes"]
    counts = wide["dataset"].value_counts()
    eligible = sorted(counts[counts >= 3].index)
    wide = wide[wide["dataset"].isin(eligible)].copy()
    if len(wide) < 6 or len(eligible) < 2:
        raise RuntimeError("independent paired model is not estimable")
    design = pd.get_dummies(wide["dataset"], drop_first=False, dtype=float).to_numpy(float)
    y = wide["difference"].to_numpy(float)
    inverse = np.linalg.inv(design.T @ design)
    coefficients = inverse @ design.T @ y
    contrast = np.repeat(1 / design.shape[1], design.shape[1])
    beta = float(contrast @ coefficients)
    residuals = y - design @ coefficients
    n, k = design.shape
    df = n - k
    variance = float(residuals @ residuals) / df * inverse
    se = float(np.sqrt(contrast @ variance @ contrast))
    leverage = np.clip(np.diag(design @ inverse @ design.T), 0, 1 - 1e-12)
    omega = (residuals / (1 - leverage)) ** 2
    hc3 = inverse @ (design.T @ (omega[:, None] * design)) @ inverse
    hc3_se = float(np.sqrt(contrast @ hc3 @ contrast))
    return {
        "beta": beta,
        "se": se,
        "pvalue": float(2 * stats.t.sf(abs(beta / se), df)) if se else 1.0,
        "hc3_se": hc3_se,
        "hc3_pvalue": float(2 * stats.t.sf(abs(beta / hc3_se), df)) if hc3_se else 1.0,
        "n_paired_donors": int(wide["donor"].nunique()),
        "n_datasets": int(len(eligible)),
    }


def check_close(observed: float, expected: float, label: str, tolerance: float = 1e-8) -> None:
    if np.isnan(observed) and np.isnan(expected):
        return
    require(
        np.isclose(observed, expected, atol=tolerance, rtol=tolerance),
        f"numeric validation failed for {label}: {observed} != {expected}",
    )


def main() -> None:
    validation_path = RESULTS / "validation_report.tsv"
    release_manifest_path = RESULTS / "release_manifest.tsv"
    retirement_path = RESULTS / "retired_five_lineage_cross_prototype.tsv"
    validated_path = CANDIDATE / "VALIDATED"
    for path in (validation_path, release_manifest_path, retirement_path, validated_path):
        refuse_existing(path)

    require((CROSS_CANDIDATE / "VALIDATED").exists(), "complete six-lineage candidate is not validated")

    rows: list[dict[str, object]] = []

    def add(check: str, passed: bool, detail: object) -> None:
        rows.append({"check": check, "passed": bool(passed), "detail": str(detail)})
        require(bool(passed), f"validation failed: {check}: {detail}")

    ready = pd.read_csv(READY, sep="\t", dtype=str)
    registry = pd.read_csv(REGISTRY, sep="\t")
    add("registry_family_117", len(registry) == 117, len(registry))
    add(
        "registry_hash_matches_ready",
        sha256(REGISTRY) == ready.loc[0, "registry_sha256"],
        sha256(REGISTRY),
    )
    add(
        "membership_hash_matches_ready",
        sha256(MEMBERSHIP) == ready.loc[0, "membership_table_sha256"],
        sha256(MEMBERSHIP),
    )

    effects = pd.read_csv(RESULTS / "ambient_program_effects.tsv", sep="\t")
    reproduction = pd.read_csv(RESULTS / "raw_score_reproduction.tsv", sep="\t")
    donor_scores = pd.read_csv(RESULTS / "donor_program_scores.tsv.gz", sep="\t")
    donor_roster = pd.read_csv(
        RESULTS / "current_source_diagnosis_donor_roster.tsv", sep="\t"
    )
    dataset_qc = pd.read_csv(RESULTS / "ambient_dataset_lineage_qc.tsv", sep="\t")
    burden = pd.read_csv(RESULTS / "program_dataset_ambient_burden.tsv", sep="\t")
    add("effect_rows_117_x_2", len(effects) == 234, len(effects))
    add("reproduction_rows_117", len(reproduction) == 117, len(reproduction))
    add(
        "minimum_raw_reproduction_r",
        float(reproduction["reproduction_r"].min()) >= 0.995,
        reproduction["reproduction_r"].min(),
    )
    transport_pass = (
        reproduction["transport_gate_pass"].astype(str).str.lower().eq("true")
    )
    add(
        "transport_reproduction_gate_has_106_testable_and_11_tcell_untestable_programs",
        int(transport_pass.sum()) == 106
        and set(reproduction.loc[~transport_pass, "cell_type"]) == {"tcells"},
        reproduction.groupby(["cell_type", "transport_gate_pass"]).size().to_dict(),
    )
    add(
        "program_keys_complete",
        set(effects["program_uid"]) == set(registry["program_uid"]),
        effects["program_uid"].nunique(),
    )
    add(
        "donor_score_key_unique",
        not donor_scores[["program_uid", "donor"]].duplicated().any(),
        len(donor_scores),
    )
    current_roster = load_donor_metadata()[
        [
            "donor", "dataset", "disease_stage_coarse", "exclude_stage_analysis",
            "exclude", "stage_ordinal",
        ]
    ].sort_values("donor").reset_index(drop=True)
    donor_roster = donor_roster.sort_values("donor").reset_index(drop=True)
    for frame in (current_roster, donor_roster):
        for column in ("exclude_stage_analysis", "exclude"):
            frame[column] = (
                frame[column].astype(str).str.strip().str.lower().map(
                    {"true": True, "false": False}
                )
            )
    add(
        "biological_donor_and_diagnosis_roster_matches_current_release",
        donor_roster.astype(str).equals(current_roster.astype(str)),
        {"rows": len(donor_roster), "datasets": donor_roster["dataset"].nunique()},
    )
    add("dataset_lineage_qc_7_x_5", len(dataset_qc) == 35, len(dataset_qc))
    required_qc_columns = {
        "correction_status", "failure_reason", "failure_detail", "n_cells",
        "n_biological_donors",
        "n_ambient_batches", "n_pooled_small_samples",
        "n_cells_with_contamination_estimate", "contamination_median",
        "contamination_mean", "contamination_q10", "contamination_q90",
    }
    add(
        "dataset_lineage_qc_contract_complete",
        required_qc_columns <= set(dataset_qc.columns),
        sorted(required_qc_columns),
    )
    corrected_qc = dataset_qc[dataset_qc["correction_status"] == "corrected"]
    failed_qc = dataset_qc[dataset_qc["correction_status"] != "corrected"]
    add(
        "corrected_datasets_use_unpooled_sequencing_sample_batches",
        (corrected_qc["n_pooled_small_samples"].astype(int) == 0).all(),
        sorted(corrected_qc["dataset"].unique()),
    )
    add(
        "corrected_lineages_have_contamination_summaries_when_cells_present",
        (
            corrected_qc.loc[corrected_qc["n_cells"] > 0, "n_cells_with_contamination_estimate"]
            .astype(int)
            > 0
        ).all(),
        int((corrected_qc["n_cells"] > 0).sum()),
    )
    add(
        "failed_datasets_visible_with_specific_reason",
        not failed_qc.empty
        and failed_qc["failure_reason"].fillna("").str.len().gt(0).all()
        and failed_qc["failure_detail"].fillna("").str.len().gt(0).all(),
        sorted(failed_qc["dataset"].unique()),
    )
    add("ambient_burden_117_x_7", len(burden) == 819, len(burden))
    primary_states = set(
        effects.loc[
            (effects["analysis_universe"] == "complete_case_common_universe")
            & effects["primary_selected_frozen"].astype(bool),
            "evidence_state",
        ]
    )
    add(
        "primary_evidence_states_predeclared",
        primary_states
        <= {
            "retained",
            "ambient_sensitive_attenuation",
            "attenuation_indeterminate",
            "reversal_supported",
            "untestable",
        },
        sorted(primary_states),
    )
    nonselected = effects[
        (effects["analysis_universe"] == "complete_case_common_universe")
        & ~effects["primary_selected_frozen"].astype(bool)
    ]
    add(
        "nonselected_state_is_not_applicable_or_untestable",
        set(nonselected["evidence_state"])
        <= {"not_applicable_not_frozen_disease_association", "untestable"},
        sorted(nonselected["evidence_state"].unique()),
    )
    tcell_primary = effects[
        (effects["analysis_universe"] == "complete_case_common_universe")
        & (effects["cell_type"] == "tcells")
    ]
    add(
        "tcell_ambient_transport_is_explicitly_untestable",
        len(tcell_primary) == 11
        and tcell_primary["raw_beta"].notna().all()
        and tcell_primary["corrected_beta"].isna().all()
        and tcell_primary["delta_beta"].isna().all()
        and set(tcell_primary["evidence_state"]) == {"untestable"},
        len(tcell_primary),
    )
    untestable = effects[
        (effects["analysis_universe"] == "complete_case_common_universe")
        & (effects["evidence_state"] == "untestable")
    ]
    failure_columns = [
        "raw_failure_reason", "corrected_failure_reason", "delta_failure_reason",
        "score_reproduction_failure_reason",
    ]
    add(
        "untestable_requires_model_failure",
        untestable.empty
        or untestable[failure_columns]
        .fillna("")
        .astype(str)
        .apply(lambda row: any(bool(value) for value in row), axis=1)
        .all(),
        len(untestable),
    )
    gse189600 = dataset_qc[dataset_qc["dataset"] == "GSE189600"]
    add("GSE189600_explicit_five_lineages", len(gse189600) == 5, len(gse189600))
    add(
        "GSE189600_passthrough_explicit",
        set(gse189600["correction_status"]) == {"uncorrected_passthrough"},
        sorted(gse189600["correction_status"].unique()),
    )

    for universe, group in effects.groupby("analysis_universe"):
        add(f"family_size:{universe}", len(group) == 117, len(group))
        for prefix in ("raw", "corrected", "delta"):
            q = bh_adjust(group[f"{prefix}_pvalue"], family_size=117)
            hc3_q = bh_adjust(group[f"{prefix}_hc3_pvalue"], family_size=117)
            add(
                f"BH_exact:{universe}:{prefix}",
                np.allclose(q, group[f"{prefix}_qvalue"], equal_nan=True, atol=1e-12),
                117,
            )
            add(
                f"HC3_BH_exact:{universe}:{prefix}",
                np.allclose(
                    hc3_q, group[f"{prefix}_hc3_qvalue"], equal_nan=True, atol=1e-12
                ),
                117,
            )

    failed = set(
        dataset_qc.loc[
            dataset_qc["correction_status"] != "corrected", "dataset"
        ].astype(str)
    )
    metadata = load_donor_metadata()
    for result in effects.itertuples(index=False):
        frame = donor_scores[donor_scores["program_uid"] == result.program_uid][
            ["donor", "raw", "corrected"]
        ].merge(metadata, on="donor", how="inner", validate="one_to_one")
        if result.analysis_universe == "complete_case_common_universe":
            frame = frame[~frame["dataset"].isin(failed)].copy()
        frame["delta"] = frame["corrected"] - frame["raw"]
        for prefix in ("raw", "corrected", "delta"):
            observed_beta = float(getattr(result, f"{prefix}_beta"))
            if not np.isfinite(observed_beta):
                require(
                    bool(str(getattr(result, f"{prefix}_failure_reason")))
                    or bool(str(result.score_reproduction_failure_reason)),
                    f"missing estimate lacks failure reason: {result.program_uid}:{prefix}",
                )
                continue
            fit = independent_fit(frame, prefix)
            for metric in ("beta", "se", "pvalue", "hc3_se", "hc3_pvalue"):
                check_close(
                    float(getattr(result, f"{prefix}_{metric}")),
                    float(fit[metric]),
                    f"{result.program_uid}:{result.analysis_universe}:{prefix}:{metric}",
                )
    add("all_234_effect_rows_recomputed", True, 234)

    full = effects[
        effects["analysis_universe"] == "full_universe_passthrough_sensitivity"
    ].merge(
        registry[["program_uid", "primary_beta", "primary_se", "primary_hc3_se"]],
        on="program_uid",
        validate="one_to_one",
    )
    add(
        "full_raw_beta_matches_registry",
        np.allclose(full["raw_beta"], full["primary_beta"], atol=1e-8, rtol=1e-8),
        float(np.max(np.abs(full["raw_beta"] - full["primary_beta"]))),
    )
    add(
        "full_raw_se_matches_registry",
        np.allclose(full["raw_se"], full["primary_se"], atol=1e-8, rtol=1e-8),
        float(np.max(np.abs(full["raw_se"] - full["primary_se"]))),
    )
    add(
        "full_raw_hc3_matches_registry",
        np.allclose(full["raw_hc3_se"], full["primary_hc3_se"], atol=1e-8, rtol=1e-8),
        float(np.max(np.abs(full["raw_hc3_se"] - full["primary_hc3_se"]))),
    )

    contrasts = pd.read_csv(CROSS_RESULTS / "hero_lineage_contrasts.tsv", sep="\t")
    hero_scores = pd.read_csv(CROSS_RESULTS / "hero_lineage_scores.tsv.gz", sep="\t")
    promotion = pd.read_csv(CROSS_RESULTS / "hero_lineage_promotion.tsv", sep="\t")
    add(
        "complete_six_lineage_candidate_validated",
        (CROSS_CANDIDATE / "VALIDATED").exists(),
        sha256(CROSS_CANDIDATE / "VALIDATED"),
    )
    add("hero_contrast_rows_40", len(contrasts) == 40, len(contrasts))
    primary_tcell = contrasts[
        (contrasts["annotation_filter"] == "all_annotated_cells")
        & (contrasts["minimum_cells"] == 50)
        & (contrasts["comparison_lineage"] == "T cells")
    ]
    add(
        "primary_tcell_rows_estimable",
        len(primary_tcell) == 2 and primary_tcell["beta"].notna().all(),
        primary_tcell[["program_name", "beta"]].to_dict("records"),
    )
    for family, group in contrasts.groupby(["annotation_filter", "minimum_cells"]):
        add(f"hero_family_size:{family}", len(group) == 10, len(group))
        add(
            f"hero_BH_exact:{family}",
            np.allclose(
                bh_adjust(group["pvalue"]), group["qvalue"], equal_nan=True, atol=1e-12
            ),
            10,
        )
        add(
            f"hero_HC3_BH_exact:{family}",
            np.allclose(
                bh_adjust(group["hc3_pvalue"]),
                group["hc3_qvalue"],
                equal_nan=True,
                atol=1e-12,
            ),
            10,
        )
    for result in contrasts.itertuples(index=False):
        subset = hero_scores[
            (hero_scores["annotation_filter"] == result.annotation_filter)
            & (hero_scores["minimum_cells"] == result.minimum_cells)
        ]
        try:
            fit = independent_paired_fit(subset, result.program_uid, result.comparison_lineage)
        except RuntimeError:
            add(
                f"hero_unestimable_matches:{result.program_uid}:{result.comparison}",
                pd.isna(result.beta),
                result.beta,
            )
            continue
        for metric in ("beta", "se", "pvalue", "hc3_se", "hc3_pvalue"):
            check_close(
                float(getattr(result, metric)),
                float(fit[metric]),
                f"hero:{result.program_uid}:{result.annotation_filter}:{result.minimum_cells}:{result.comparison_lineage}:{metric}",
            )
    add("all_40_hero_contrasts_recomputed", True, 40)
    add("promotion_rows_two", len(promotion) == 2, len(promotion))
    family_eligible = bool(promotion["program_pass"].all())
    add(
        "promotion_destination_consistent",
        set(promotion["main_figure_family_eligible"].astype(bool)) == {family_eligible}
        and set(promotion["destination"])
        == ({"Figure_4E"} if family_eligible else {"Figure_S4"}),
        family_eligible,
    )

    pd.DataFrame(
        [
            {
                "artifact_family": "hero_lineage_* in this ambient candidate",
                "status": "retired_not_for_claims",
                "reason": "prototype atlas omitted T cells",
                "replacement_candidate": str(CROSS_CANDIDATE.resolve()),
                "replacement_validated_sha256": sha256(CROSS_CANDIDATE / "VALIDATED"),
            }
        ]
    ).to_csv(retirement_path, sep="\t", index=False)
    validation = pd.DataFrame(rows)
    validation.to_csv(validation_path, sep="\t", index=False)
    candidate_files = sorted(
        path for path in RESULTS.rglob("*") if path.is_file() and path != release_manifest_path
    )
    release = pd.DataFrame(
        [
            {
                "relative_path": str(path.relative_to(CANDIDATE)),
                "bytes": int(path.stat().st_size),
                "sha256": sha256(path),
            }
            for path in candidate_files
        ]
    )
    release.to_csv(release_manifest_path, sep="\t", index=False)
    validated_path.write_text(
        "status\tvalidated\n"
        f"registry_sha256\t{sha256(REGISTRY)}\n"
        f"membership_sha256\t{sha256(MEMBERSHIP)}\n"
        f"validation_report_sha256\t{sha256(validation_path)}\n"
        f"release_manifest_sha256\t{sha256(release_manifest_path)}\n"
        f"complete_cross_lineage_validated_sha256\t{sha256(CROSS_CANDIDATE / 'VALIDATED')}\n",
        encoding="utf-8",
    )
    print(f"PASS {len(rows)} validation checks; {len(release)} release artifacts")


if __name__ == "__main__":
    main()
