#!/usr/bin/env python3
"""Frozen statistical primitives for the one-pass Myojin phase-C analysis."""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.contingency_tables import StratifiedTable


CLASS_ORDER = ("neither", "genetic_only", "disease_state_only", "convergent")
PAIRWISE = (
    ("genetic_only", "genetic_only_vs_neither"),
    ("disease_state_only", "disease_state_only_vs_neither"),
    ("convergent", "convergent_vs_neither"),
)


def bh_adjust(values: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg adjustment with deterministic monotonic correction."""

    p = np.asarray(values, dtype=float)
    if (
        p.ndim != 1
        or len(p) == 0
        or np.any(~np.isfinite(p))
        or np.any((p < 0) | (p > 1))
    ):
        raise ValueError("BH inputs must be finite probabilities")
    order = np.argsort(p, kind="stable")
    ranked = p[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(p) / np.arange(1, len(p) + 1))[::-1]
    )[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result.tolist()


def rank_first_bins(
    values: Sequence[float], symbols: Sequence[str], n_bins: int = 5
) -> np.ndarray:
    """Assign deterministic rank-first bins using value then symbol ordering."""

    numeric = np.asarray(values, dtype=float)
    names = np.asarray(symbols, dtype=str)
    if len(numeric) != len(names) or len(numeric) == 0 or np.any(~np.isfinite(numeric)):
        raise ValueError("Rank-bin inputs must be nonempty, aligned, and finite")
    order = np.lexsort((names, numeric))
    bins = np.empty(len(numeric), dtype=np.int8)
    bins[order] = np.minimum(
        n_bins, (np.arange(len(numeric)) * n_bins // len(numeric)) + 1
    )
    return bins


def guide_bin(values: Sequence[int]) -> np.ndarray:
    numeric = np.asarray(values, dtype=int)
    if np.any(numeric < 4):
        raise ValueError("Guide-bin input contains a value below the primary floor")
    return np.minimum(numeric, 8)


def add_frozen_strata(frame: pd.DataFrame) -> pd.DataFrame:
    """Add strata once, or preserve and validate strata frozen upstream.

    Sensitivity subsets must inherit the bins defined on their prespecified
    parent universe.  Re-ranking after removing known hits or raising the guide
    floor would make the null outcome-dependent.  The executor therefore adds
    bins to each parent universe before subsetting; this function is also safe
    for standalone fixtures that do not yet carry frozen bins.
    """

    result = frame.copy()
    frozen = {
        "guide_bin",
        "expression_quintile",
        "chronos_missing",
        "chronos_quintile",
        "permutation_stratum",
    }
    present = frozen.intersection(result.columns)
    if present and present != frozen:
        raise ValueError(f"Partial frozen-stratum columns: {sorted(present)}")
    if present == frozen:
        expected_guide = guide_bin(result["guide_count_min"].to_numpy())
        if not np.array_equal(result["guide_bin"].to_numpy(dtype=int), expected_guide):
            raise ValueError("Frozen guide bins do not match guide counts")
        missing = result["HLF_Chronos_missing"].to_numpy(dtype=int)
        if not np.array_equal(result["chronos_missing"].to_numpy(dtype=int), missing):
            raise ValueError("Frozen Chronos-missingness stratum is inconsistent")
    else:
        result["guide_bin"] = guide_bin(result["guide_count_min"].to_numpy())
        result["expression_quintile"] = rank_first_bins(
            result["log1p_HLF_TPM"].to_numpy(),
            result["gene_symbol"].astype(str).to_numpy(),
        )
        result["chronos_missing"] = result["HLF_Chronos_missing"].astype(int)
        nonmissing = result["chronos_missing"].to_numpy() == 0
        chronos_quintile = np.zeros(len(result), dtype=np.int8)
        if np.any(nonmissing):
            chronos_quintile[nonmissing] = rank_first_bins(
                result.loc[nonmissing, "HLF_Chronos"].to_numpy(),
                result.loc[nonmissing, "gene_symbol"].astype(str).to_numpy(),
            )
        result["chronos_quintile"] = chronos_quintile
        result["permutation_stratum"] = [
            f"g{g}_e{e}_m{m}"
            for g, e, m in zip(
                result["guide_bin"],
                result["expression_quintile"],
                result["chronos_missing"],
            )
        ]
    for column in (
        "guide_bin",
        "expression_quintile",
        "chronos_missing",
        "chronos_quintile",
    ):
        numeric = result[column].to_numpy(dtype=int)
        if column in {"expression_quintile", "chronos_quintile"}:
            upper = 5
            lower = 0 if column == "chronos_quintile" else 1
            if np.any((numeric < lower) | (numeric > upper)):
                raise ValueError(f"Invalid frozen values in {column}")
        result[column] = numeric
    expected_labels = np.asarray(
        [
            f"g{g}_e{e}_m{m}"
            for g, e, m in zip(
                result["guide_bin"],
                result["expression_quintile"],
                result["chronos_missing"],
            )
        ],
        dtype=str,
    )
    if not np.array_equal(
        result["permutation_stratum"].astype(str).to_numpy(), expected_labels
    ):
        raise ValueError("Frozen permutation-stratum labels are inconsistent")
    return result


def _design(
    frame: pd.DataFrame, include_missing_indicator: bool
) -> tuple[np.ndarray, np.ndarray]:
    class_columns = [
        (frame["primary_evidence_class"].to_numpy() == class_name).astype(float)
        for class_name in CLASS_ORDER[1:]
    ]
    covariates = [
        np.ones(len(frame), dtype=float),
        *class_columns,
        frame["log1p_HLF_TPM"].to_numpy(dtype=float),
        frame["guide_count_min"].to_numpy(dtype=float),
        frame["HLF_Chronos"].to_numpy(dtype=float),
    ]
    if include_missing_indicator:
        covariates.append(frame["HLF_Chronos_missing"].to_numpy(dtype=float))
    full = np.column_stack(covariates)
    reduced = np.column_stack([covariates[0], *covariates[4:]])
    return full, reduced


def _stratum_indices(frame: pd.DataFrame) -> list[np.ndarray]:
    groups: dict[str, list[int]] = defaultdict(list)
    for index, label in enumerate(frame["permutation_stratum"].astype(str)):
        groups[label].append(index)
    return [np.asarray(indices, dtype=int) for _, indices in sorted(groups.items())]


def fit_class_continuous(
    frame: pd.DataFrame,
    outcome_column: str,
    analysis_variant: str,
    seed: int,
    accepted_draws: int = 100_000,
    include_missing_indicator: bool = False,
    batch_size: int = 128,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Fit the frozen OLS/HC3 model and restricted Freedman-Lane null."""

    if len(frame) == 0:
        raise ValueError(f"No rows for class analysis {analysis_variant}")
    frame = add_frozen_strata(frame.reset_index(drop=True))
    y = frame[outcome_column].to_numpy(dtype=float)
    if np.any(~np.isfinite(y)):
        raise ValueError(f"Non-finite outcome in {analysis_variant}")
    full, reduced = _design(frame, include_missing_indicator)
    rank_full = int(np.linalg.matrix_rank(full))
    rank_reduced = int(np.linalg.matrix_rank(reduced))
    if rank_full != full.shape[1] or rank_reduced != reduced.shape[1]:
        raise ValueError(
            f"Design rank failure in {analysis_variant}: "
            f"full {rank_full}/{full.shape[1]}, reduced {rank_reduced}/{reduced.shape[1]}"
        )
    if rank_full - rank_reduced != 3:
        raise ValueError("Class design does not add exactly three indicators")

    full_fit = sm.OLS(y, full).fit(cov_type="HC3")
    beta = np.asarray(full_fit.params)
    hc3_se = np.asarray(full_fit.bse)
    hc3_ci = np.asarray(full_fit.conf_int(alpha=0.05))
    pinv_full = np.linalg.pinv(full)
    q_full, _ = np.linalg.qr(full, mode="reduced")
    q_reduced, _ = np.linalg.qr(reduced, mode="reduced")
    reduced_beta = np.linalg.lstsq(reduced, y, rcond=None)[0]
    reduced_fitted = reduced @ reduced_beta
    reduced_residual = y - reduced_fitted
    rss_full = float(np.sum((y - full @ np.linalg.lstsq(full, y, rcond=None)[0]) ** 2))
    rss_reduced = float(np.sum(reduced_residual**2))
    df_extra = rank_full - rank_reduced
    df_residual = len(y) - rank_full
    if df_residual <= 0 or rss_full <= 0:
        raise ValueError(f"Invalid residual degrees/variance in {analysis_variant}")
    observed_f = max(0.0, (rss_reduced - rss_full) / df_extra) / (
        rss_full / df_residual
    )
    residual_sd = math.sqrt(rss_full / df_residual)

    groups = _stratum_indices(frame)
    rng = np.random.default_rng(seed)
    accepted = 0
    attempted = 0
    omnibus_extreme = 0
    pair_extreme = np.zeros(3, dtype=np.int64)
    null_beta_sum = np.zeros(3, dtype=float)
    null_beta_sumsq = np.zeros(3, dtype=float)
    null_f_sum = 0.0
    observed_pair = beta[1:4].copy()

    while accepted < accepted_draws:
        batch = min(batch_size, accepted_draws - accepted)
        permuted_y = np.empty((len(y), batch), dtype=float)
        for draw_index in range(batch):
            permuted_residual = reduced_residual.copy()
            for indices in groups:
                if len(indices) > 1:
                    permuted_residual[indices] = rng.permutation(
                        reduced_residual[indices]
                    )
            permuted_y[:, draw_index] = reduced_fitted + permuted_residual
        attempted += batch
        sum_y2 = np.sum(permuted_y * permuted_y, axis=0)
        rss_r = np.maximum(
            0.0, sum_y2 - np.sum((q_reduced.T @ permuted_y) ** 2, axis=0)
        )
        rss_f = np.maximum(0.0, sum_y2 - np.sum((q_full.T @ permuted_y) ** 2, axis=0))
        null_f = np.maximum(0.0, (rss_r - rss_f) / df_extra) / (rss_f / df_residual)
        null_beta = (pinv_full @ permuted_y)[1:4, :]
        finite = np.isfinite(null_f) & np.all(np.isfinite(null_beta), axis=0)
        if not np.any(finite):
            continue
        null_f = null_f[finite]
        null_beta = null_beta[:, finite]
        keep = min(len(null_f), accepted_draws - accepted)
        null_f = null_f[:keep]
        null_beta = null_beta[:, :keep]
        accepted += keep
        omnibus_extreme += int(np.sum(null_f >= observed_f))
        pair_extreme += np.sum(
            np.abs(null_beta) >= np.abs(observed_pair[:, None]), axis=1
        )
        null_beta_sum += np.sum(null_beta, axis=1)
        null_beta_sumsq += np.sum(null_beta * null_beta, axis=1)
        null_f_sum += float(np.sum(null_f))

    omnibus_p = (1 + omnibus_extreme) / (accepted + 1)
    pair_p = ((1 + pair_extreme) / (accepted + 1)).astype(float)
    pair_q = bh_adjust(pair_p.tolist())
    class_counts = frame["primary_evidence_class"].value_counts().to_dict()
    rows: list[dict[str, object]] = [
        {
            "analysis_variant": analysis_variant,
            "outcome": outcome_column,
            "contrast": "omnibus_evidence_class",
            "n_total": len(frame),
            "n_class": len(frame),
            "n_reference": 0,
            "estimate": "",
            "HC3_SE": "",
            "CI95_low": "",
            "CI95_high": "",
            "standardized_effect": "",
            "statistic": observed_f,
            "statistic_type": "partial_F",
            "permutation_p": omnibus_p,
            "BH_q": "",
            "accepted_draws": accepted,
            "attempted_draws": attempted,
            "seed": seed,
            "df_numerator": df_extra,
            "df_residual": df_residual,
            "residual_SD": residual_sd,
        }
    ]
    for index, (class_name, contrast) in enumerate(PAIRWISE):
        rows.append(
            {
                "analysis_variant": analysis_variant,
                "outcome": outcome_column,
                "contrast": contrast,
                "n_total": len(frame),
                "n_class": int(class_counts.get(class_name, 0)),
                "n_reference": int(class_counts.get("neither", 0)),
                "estimate": float(beta[index + 1]),
                "HC3_SE": float(hc3_se[index + 1]),
                "CI95_low": float(hc3_ci[index + 1, 0]),
                "CI95_high": float(hc3_ci[index + 1, 1]),
                "standardized_effect": float(beta[index + 1] / residual_sd),
                "statistic": float(beta[index + 1]),
                "statistic_type": "absolute_OLS_coefficient",
                "permutation_p": float(pair_p[index]),
                "BH_q": float(pair_q[index]),
                "accepted_draws": accepted,
                "attempted_draws": attempted,
                "seed": seed,
                "df_numerator": "",
                "df_residual": df_residual,
                "residual_SD": residual_sd,
            }
        )
    null_beta_mean = null_beta_sum / accepted
    null_beta_variance = np.maximum(
        0.0, (null_beta_sumsq - accepted * null_beta_mean**2) / max(1, accepted - 1)
    )
    audit = {
        "analysis_variant": analysis_variant,
        "outcome": outcome_column,
        "accepted_draws": accepted,
        "attempted_draws": attempted,
        "seed": seed,
        "n_rows": len(frame),
        "n_strata": len(groups),
        "n_singleton_strata": sum(len(indices) == 1 for indices in groups),
        "rank_full": rank_full,
        "rank_reduced": rank_reduced,
        "observed_omnibus_F": observed_f,
        "null_omnibus_F_mean": null_f_sum / accepted,
        "null_genetic_only_beta_mean": null_beta_mean[0],
        "null_disease_state_only_beta_mean": null_beta_mean[1],
        "null_convergent_beta_mean": null_beta_mean[2],
        "null_genetic_only_beta_SD": math.sqrt(null_beta_variance[0]),
        "null_disease_state_only_beta_SD": math.sqrt(null_beta_variance[1]),
        "null_convergent_beta_SD": math.sqrt(null_beta_variance[2]),
    }
    return rows, audit


def _risk_difference(strata: list[dict[str, int]]) -> tuple[float, float]:
    total_weight = sum(row["n_class"] + row["n_reference"] for row in strata)
    if total_weight == 0:
        return float("nan"), float("nan")
    estimate = 0.0
    variance = 0.0
    for row in strata:
        weight = (row["n_class"] + row["n_reference"]) / total_weight
        p_class = row["hits_class"] / row["n_class"]
        p_reference = row["hits_reference"] / row["n_reference"]
        estimate += weight * (p_class - p_reference)
        variance += weight**2 * (
            p_class * (1 - p_class) / row["n_class"]
            + p_reference * (1 - p_reference) / row["n_reference"]
        )
    return estimate, math.sqrt(max(0.0, variance))


def fit_class_binary(
    frame: pd.DataFrame,
    event_column: str,
    analysis_variant: str,
    seed: int,
    accepted_draws: int = 100_000,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    frame = add_frozen_strata(frame.reset_index(drop=True))
    seed_sequences = np.random.SeedSequence(seed).spawn(len(PAIRWISE))
    rows: list[dict[str, object]] = []
    audits: list[dict[str, object]] = []
    for pair_index, ((class_name, contrast), seed_sequence) in enumerate(
        zip(PAIRWISE, seed_sequences)
    ):
        subset = frame[
            frame["primary_evidence_class"].isin(["neither", class_name])
        ].copy()
        shared: list[dict[str, int]] = []
        tables = []
        for _, group in subset.groupby("permutation_stratum", sort=True):
            class_group = group[group["primary_evidence_class"] == class_name]
            reference_group = group[group["primary_evidence_class"] == "neither"]
            if len(class_group) == 0 or len(reference_group) == 0:
                continue
            hits_class = int(class_group[event_column].astype(bool).sum())
            hits_reference = int(reference_group[event_column].astype(bool).sum())
            row = {
                "n_class": len(class_group),
                "n_reference": len(reference_group),
                "hits_class": hits_class,
                "hits_reference": hits_reference,
            }
            shared.append(row)
            tables.append(
                np.asarray(
                    [
                        [hits_class, hits_reference],
                        [
                            len(class_group) - hits_class,
                            len(reference_group) - hits_reference,
                        ],
                    ],
                    dtype=float,
                )
            )
        observed_rd, rd_se = _risk_difference(shared)
        if not math.isfinite(observed_rd):
            raise ValueError(
                f"No shared binary strata for {analysis_variant} {contrast}"
            )
        rng = np.random.default_rng(seed_sequence)
        total_weight = sum(row["n_class"] + row["n_reference"] for row in shared)
        null_rd = np.zeros(accepted_draws, dtype=float)
        for stratum in shared:
            n_class = stratum["n_class"]
            n_reference = stratum["n_reference"]
            hits = stratum["hits_class"] + stratum["hits_reference"]
            n_total = n_class + n_reference
            permuted_class_hits = rng.hypergeometric(
                hits, n_total - hits, n_class, size=accepted_draws
            )
            permuted_reference_hits = hits - permuted_class_hits
            weight = n_total / total_weight
            null_rd += weight * (
                permuted_class_hits / n_class - permuted_reference_hits / n_reference
            )
        p_value = (1 + int(np.sum(np.abs(null_rd) >= abs(observed_rd)))) / (
            accepted_draws + 1
        )
        odds_ratio = float("nan")
        or_low = float("nan")
        or_high = float("nan")
        if tables:
            try:
                stratified = StratifiedTable(
                    np.stack(tables, axis=2), shift_zeros=False
                )
                odds_ratio = float(stratified.oddsratio_pooled)
                or_low, or_high = map(
                    float, stratified.oddsratio_pooled_confint(alpha=0.05)
                )
            except Exception:
                pass
        rows.append(
            {
                "analysis_variant": analysis_variant,
                "event": event_column,
                "contrast": contrast,
                "n_class": int(sum(row["n_class"] for row in shared)),
                "n_reference": int(sum(row["n_reference"] for row in shared)),
                "hits_class": int(sum(row["hits_class"] for row in shared)),
                "hits_reference": int(sum(row["hits_reference"] for row in shared)),
                "risk_difference": observed_rd,
                "RD_CI95_low": observed_rd - 1.959963984540054 * rd_se,
                "RD_CI95_high": observed_rd + 1.959963984540054 * rd_se,
                "odds_ratio_MH": odds_ratio,
                "OR_CI95_low": or_low,
                "OR_CI95_high": or_high,
                "permutation_p": p_value,
                "BH_q": "",
                "accepted_draws": accepted_draws,
                "attempted_draws": accepted_draws,
                "seed": seed,
                "rng_substream_index": pair_index,
                "n_shared_strata": len(shared),
            }
        )
        audits.append(
            {
                "analysis_variant": analysis_variant,
                "event": event_column,
                "contrast": contrast,
                "accepted_draws": accepted_draws,
                "attempted_draws": accepted_draws,
                "seed": seed,
                "rng_substream_index": pair_index,
                "n_shared_strata": len(shared),
                "observed_risk_difference": observed_rd,
                "null_risk_difference_mean": float(np.mean(null_rd)),
                "null_risk_difference_SD": float(np.std(null_rd, ddof=1)),
            }
        )
    q_values = bh_adjust([float(row["permutation_p"]) for row in rows])
    for row, q_value in zip(rows, q_values):
        row["BH_q"] = q_value
    return rows, audits


def weighted_score(values: Sequence[float], weights: Sequence[float]) -> float:
    numeric = np.asarray(values, dtype=float)
    weight = np.asarray(weights, dtype=float)
    if len(numeric) == 0 or len(numeric) != len(weight):
        raise ValueError("Weighted-score inputs are empty or misaligned")
    if (
        np.any(~np.isfinite(numeric))
        or np.any(~np.isfinite(weight))
        or np.any(weight <= 0)
    ):
        raise ValueError("Weighted-score inputs must be finite with positive weights")
    weight = weight / np.sum(weight)
    return float(np.sum(numeric * weight))


def matched_program_test(
    frame: pd.DataFrame,
    members: list[dict[str, object]],
    program_member_symbols: set[str],
    stage_direction: int,
    program_uid: str,
    analysis_variant: str,
    weight_mode: str,
    seed_sequence: np.random.SeedSequence,
    draws: int = 10_000,
) -> tuple[dict[str, object], dict[str, object]]:
    if stage_direction not in {-1, 1}:
        raise ValueError("Program stage direction must be -1 or 1")
    if frame["gene_symbol"].duplicated().any():
        raise ValueError("Program test frame contains duplicate gene symbols")
    member_names = [str(member["gene_symbol"]) for member in members]
    if len(member_names) != len(set(member_names)):
        raise ValueError(f"Program {program_uid} contains duplicate frozen members")
    if set(member_names) != set(program_member_symbols):
        raise ValueError(f"Program {program_uid} member-exclusion set is inconsistent")
    indexed = frame.set_index("gene_symbol", drop=False)
    eligible = [
        member for member in members if str(member["gene_symbol"]) in indexed.index
    ]
    retained_weight = sum(float(member["original_l1_weight"]) for member in eligible)
    if weight_mode == "leave_highest_weight_out" and eligible:
        drop = sorted(
            eligible,
            key=lambda row: (
                -float(row["original_l1_weight"]),
                str(row["gene_symbol"]),
            ),
        )[0]
        eligible = [member for member in eligible if member is not drop]
        retained_weight = sum(
            float(member["original_l1_weight"]) for member in eligible
        )
    testable = len(eligible) >= 8 and retained_weight >= 0.20
    reason = (
        "testable" if testable else "fewer_than_8_genes_or_retained_weight_below_0.20"
    )
    base_result: dict[str, object] = {
        "analysis_variant": analysis_variant,
        "weight_mode": weight_mode,
        "program_uid": program_uid,
        "n_eligible_genes": len(eligible),
        "retained_original_L1_weight": retained_weight,
        "testable": str(testable).upper(),
        "testability_reason": reason,
        "stage_direction": stage_direction,
        "observed_signed_score": "",
        "matched_null_mean": "",
        "matched_null_SD": "",
        "signed_effect": "",
        "empirical_p": "",
        "BH_q": "",
        "draws": 0,
        "seed_entropy": int(seed_sequence.entropy),
    }
    audit: dict[str, object] = {
        "analysis_variant": analysis_variant,
        "weight_mode": weight_mode,
        "program_uid": program_uid,
        "testable": str(testable).upper(),
        "draws": 0,
        "n_match_cells": 0,
        "minimum_control_pool_per_cell": 0,
        "null_signed_mean": "",
        "null_signed_SD": "",
        "control_reuse_within_draw": "FALSE",
        "control_reuse_between_draws": "TRUE",
        "failure_reason": "" if testable else reason,
    }
    if not testable:
        return base_result, audit

    symbols = [str(member["gene_symbol"]) for member in eligible]
    if weight_mode == "equal_weight":
        weights = np.repeat(1.0 / len(eligible), len(eligible))
    else:
        weights = np.asarray(
            [float(member["original_l1_weight"]) for member in eligible], dtype=float
        )
        weights /= np.sum(weights)
    target = indexed.loc[symbols]
    observed = stage_direction * weighted_score(target["outcome_z"].to_numpy(), weights)
    cell_columns = [
        "guide_bin",
        "expression_quintile",
        "chronos_quintile",
        "chronos_missing",
    ]
    target_cells: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for position, cell in enumerate(
        target[cell_columns].itertuples(index=False, name=None)
    ):
        target_cells[tuple(map(int, cell))].append(position)
    # The sealed specification excludes the members of the program currently
    # being tested, not members of every program in the external family.
    control_frame = frame[~frame["gene_symbol"].isin(program_member_symbols)]
    control_cells: dict[tuple[int, ...], np.ndarray] = {}
    for cell, group in control_frame.groupby(cell_columns, sort=True):
        control_cells[tuple(map(int, cell if isinstance(cell, tuple) else (cell,)))] = (
            group["outcome_z"].to_numpy(dtype=float)
        )
    minimum_pool = math.inf
    for cell, positions in target_cells.items():
        pool = control_cells.get(cell, np.asarray([], dtype=float))
        minimum_pool = min(minimum_pool, len(pool))
        if len(pool) < len(positions):
            base_result["testable"] = "FALSE"
            base_result["testability_reason"] = (
                "insufficient_exact_match_cell_no_relaxation"
            )
            audit["testable"] = "FALSE"
            audit["failure_reason"] = "insufficient_exact_match_cell_no_relaxation"
            audit["minimum_control_pool_per_cell"] = int(minimum_pool)
            return base_result, audit

    rng = np.random.default_rng(seed_sequence)
    null_scores = np.zeros(draws, dtype=float)
    for cell, positions in sorted(target_cells.items()):
        pool = control_cells[cell]
        positions_array = np.asarray(positions, dtype=int)
        cell_weights = weights[positions_array]
        if len(positions) == 1:
            chosen = rng.choice(pool, size=draws, replace=True)
            null_scores += cell_weights[0] * chosen
        else:
            for draw_index in range(draws):
                chosen = rng.choice(pool, size=len(positions), replace=False)
                null_scores[draw_index] += float(np.sum(cell_weights * chosen))
    null_signed = stage_direction * null_scores
    null_mean = float(np.mean(null_signed))
    null_sd = float(np.std(null_signed, ddof=1))
    effect = observed - null_mean
    p_value = (1 + int(np.sum(null_signed >= observed))) / (draws + 1)
    base_result.update(
        {
            "observed_signed_score": observed,
            "matched_null_mean": null_mean,
            "matched_null_SD": null_sd,
            "signed_effect": effect,
            "empirical_p": p_value,
            "draws": draws,
        }
    )
    audit.update(
        {
            "draws": draws,
            "n_match_cells": len(target_cells),
            "minimum_control_pool_per_cell": int(minimum_pool),
            "null_signed_mean": null_mean,
            "null_signed_SD": null_sd,
        }
    )
    return base_result, audit


def apply_program_bh(rows: list[dict[str, object]]) -> None:
    families: dict[tuple[str, str], list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        if row["testable"] == "TRUE" and row["empirical_p"] != "":
            families[(str(row["analysis_variant"]), str(row["weight_mode"]))].append(
                index
            )
    for indices in families.values():
        adjusted = bh_adjust([float(rows[index]["empirical_p"]) for index in indices])
        for index, q_value in zip(indices, adjusted):
            rows[index]["BH_q"] = q_value


def mechanical_figure_verdict(
    source_gate_pass: bool,
    class_rows: Iterable[dict[str, object]],
    program_rows: Iterable[dict[str, object]],
    alpha: float = 0.05,
    class_magnitude: float = 0.20,
    program_magnitude: float = 0.20,
) -> dict[str, object]:
    classes = list(class_rows)
    programs = list(program_rows)
    primary_omnibus = next(
        row
        for row in classes
        if row["analysis_variant"] == "primary"
        and row["contrast"] == "omnibus_evidence_class"
    )
    primary_pairs = [
        row
        for row in classes
        if row["analysis_variant"] == "primary"
        and row["contrast"] != "omnibus_evidence_class"
    ]
    known_omnibus = next(
        row
        for row in classes
        if row["analysis_variant"] == "remove_all_known_hits"
        and row["contrast"] == "omnibus_evidence_class"
    )
    known_pairs = {
        str(row["contrast"]): row
        for row in classes
        if row["analysis_variant"] == "remove_all_known_hits"
        and row["contrast"] != "omnibus_evidence_class"
    }
    qualifying_class = []
    for row in primary_pairs:
        known = known_pairs.get(str(row["contrast"]))
        if (
            float(row["BH_q"]) < alpha
            and float(row["estimate"]) > 0
            and float(row["standardized_effect"]) >= class_magnitude
            and known is not None
            and float(known["BH_q"]) < alpha
            and float(known["estimate"]) > 0
            and float(known["standardized_effect"]) >= class_magnitude
        ):
            qualifying_class.append(row)
    known_class_omnibus_pass = float(known_omnibus["permutation_p"]) < alpha
    class_branch = (
        float(primary_omnibus["permutation_p"]) < alpha
        and known_class_omnibus_pass
        and bool(qualifying_class)
    )

    primary_programs = {
        str(row["program_uid"]): row
        for row in programs
        if row["analysis_variant"] == "primary"
        and row["weight_mode"] == "original_weight"
    }
    equal_programs = {
        str(row["program_uid"]): row
        for row in programs
        if row["analysis_variant"] == "primary" and row["weight_mode"] == "equal_weight"
    }
    leave_programs = {
        str(row["program_uid"]): row
        for row in programs
        if row["analysis_variant"] == "primary"
        and row["weight_mode"] == "leave_highest_weight_out"
    }
    known_programs = {
        str(row["program_uid"]): row
        for row in programs
        if row["analysis_variant"] == "remove_all_known_hits"
        and row["weight_mode"] == "original_weight"
    }
    qualifying_program = []
    for uid, row in primary_programs.items():
        equal = equal_programs.get(uid)
        leave = leave_programs.get(uid)
        known = known_programs.get(uid)
        if (
            row["testable"] == "TRUE"
            and row["BH_q"] != ""
            and float(row["BH_q"]) < alpha
            and float(row["signed_effect"]) >= program_magnitude
            and equal is not None
            and equal["testable"] == "TRUE"
            and float(equal["signed_effect"]) > 0
            and leave is not None
            and leave["testable"] == "TRUE"
            and float(leave["signed_effect"]) > 0
            and known is not None
            and known["testable"] == "TRUE"
            and known["BH_q"] != ""
            and float(known["BH_q"]) < alpha
            and float(known["signed_effect"]) >= program_magnitude
        ):
            qualifying_program.append(row)
    program_branch = bool(qualifying_program)
    known_hit_survival = bool(
        (known_class_omnibus_pass and qualifying_class) or qualifying_program
    )
    main = source_gate_pass and (class_branch or program_branch)
    reasons = []
    if not source_gate_pass:
        reasons.append("source_gate_failed")
    if not class_branch:
        reasons.append("class_branch_failed")
    if not program_branch:
        reasons.append("program_branch_failed")
    if main:
        reasons.append("at_least_one_prespecified_branch_passed")
    qualifying_effects = [float(row["standardized_effect"]) for row in qualifying_class]
    qualifying_effects += [float(row["signed_effect"]) for row in qualifying_program]
    return {
        "main_figure_eligible": str(main).upper(),
        "source_gate_pass": str(source_gate_pass).upper(),
        "class_omnibus_pass": str(
            float(primary_omnibus["permutation_p"]) < alpha
        ).upper(),
        "known_hit_class_omnibus_pass": str(known_class_omnibus_pass).upper(),
        "class_pairwise_branch_pass": str(class_branch).upper(),
        "program_branch_pass": str(program_branch).upper(),
        "known_hit_survival_pass": str(known_hit_survival).upper(),
        "class_minimum_standardized_effect": class_magnitude,
        "program_minimum_signed_effect": program_magnitude,
        "observed_qualifying_effect": max(qualifying_effects)
        if qualifying_effects
        else "",
        "qualifying_class_contrasts": ";".join(
            str(row["contrast"]) for row in qualifying_class
        ),
        "qualifying_programs": ";".join(
            str(row["program_uid"]) for row in qualifying_program
        ),
        "verdict_reason_codes": ";".join(reasons),
    }
