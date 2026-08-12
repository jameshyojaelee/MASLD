#!/usr/bin/env python3
"""Candidate-only lncRNA versus protein-coding observability analysis.

The analysis consumes assay-native observability calls prepared by accepted
source adapters.  It never converts statistical significance, testability, or
an evidence state into observability.  Chromatin can observe a regulatory
locus but not an RNA molecule; proteomics is not applicable to lncRNA
molecules.
"""

from __future__ import annotations

import hashlib
import json
import math
from bisect import bisect_left
from dataclasses import dataclass
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression


MATCH_SEED = 20260811
BOOTSTRAP_REPLICATES = 10_000
CALIPER_SD = 0.2
MIN_MATCH_RATE = 0.70
MAX_ABS_SMD = 0.10

GENE_REQUIRED_COLUMNS = (
    "ensembl_id",
    "gene_biotype",
    "bulk_abundance",
    "expression_variability",
    "gene_length_bp",
    "cohort_detection_count",
)

COVERAGE_REQUIRED_COLUMNS = (
    "ensembl_id",
    "assay_id",
    "assay_family",
    "molecular_scope",
    "applicability_state",
    "observable",
)

MATCH_COVARIATES = (
    "bulk_abundance",
    "expression_variability",
    "log_gene_length",
    "cohort_detection_count",
)

BIOLOGICAL_BIOTYPES = {"lncRNA", "protein_coding"}
APPLICABILITY_STATES = {
    "applicable",
    "untestable",
    "not_applicable",
    "source_dependent",
    "indeterminate",
}
OBSERVABILITY_DENOMINATOR_STATES = {"applicable", "source_dependent"}
MOLECULAR_SCOPES = {"rna_molecule", "protein_molecule", "regulatory_locus"}
ASSAY_FAMILIES = {
    "bulk_rna",
    "liver_eqtl",
    "single_cell_rna",
    "spatial_rna",
    "chromatin",
    "proteomics",
}


class ObservabilityContractError(ValueError):
    """Fail-closed input or analysis error with a stable code."""

    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


@dataclass(frozen=True)
class MatchResult:
    pairs: pd.DataFrame
    balance: pd.DataFrame
    eligible_genes: pd.DataFrame
    caliper: float
    match_rate: float
    max_abs_matched_smd: float
    gate_status: str
    gate_reason: str


def _require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ObservabilityContractError("MISSING_COLUMNS", f"{label}: {missing}")


def _stable_assay_seed(assay_id: str, seed: int) -> int:
    digest = hashlib.sha256(f"{seed}|{assay_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % (2**32)


def _as_nullable_boolean(series: pd.Series, label: str) -> pd.Series:
    if str(series.dtype) == "boolean":
        return series
    mapping = {
        True: True,
        False: False,
        "true": True,
        "false": False,
        "TRUE": True,
        "FALSE": False,
        "1": True,
        "0": False,
        "": pd.NA,
        None: pd.NA,
    }
    values = []
    for value in series.tolist():
        if pd.isna(value):
            values.append(pd.NA)
        elif value in mapping:
            values.append(mapping[value])
        else:
            raise ObservabilityContractError("BOOLEAN", f"{label}: {value!r}")
    return pd.Series(values, index=series.index, dtype="boolean")


def prepare_gene_covariates(genes: pd.DataFrame) -> pd.DataFrame:
    """Validate, transform, and deterministically order the matching universe."""

    _require_columns(genes, GENE_REQUIRED_COLUMNS, "genes")
    frame = genes.loc[:, GENE_REQUIRED_COLUMNS].copy()
    frame["ensembl_id"] = frame["ensembl_id"].astype("string").str.strip()
    frame["gene_biotype"] = frame["gene_biotype"].astype("string").str.strip()
    if frame["ensembl_id"].isna().any() or (frame["ensembl_id"] == "").any():
        raise ObservabilityContractError(
            "ENSEMBL_ID", "stable Ensembl IDs are required"
        )
    if frame["ensembl_id"].duplicated().any():
        duplicates = frame.loc[
            frame["ensembl_id"].duplicated(False), "ensembl_id"
        ].tolist()
        raise ObservabilityContractError("DUPLICATE_GENE", f"{duplicates[:5]}")
    unexpected = sorted(set(frame["gene_biotype"].dropna()) - BIOLOGICAL_BIOTYPES)
    if unexpected:
        raise ObservabilityContractError("BIOTYPE", f"unexpected values: {unexpected}")

    numeric = (
        "bulk_abundance",
        "expression_variability",
        "gene_length_bp",
        "cohort_detection_count",
    )
    for column in numeric:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    finite = np.isfinite(frame[list(numeric)].to_numpy(dtype=float)).all(axis=1)
    frame = frame.loc[finite].copy()
    if frame.empty:
        raise ObservabilityContractError(
            "NO_COMPLETE_GENES", "no complete matching records"
        )
    if (frame["gene_length_bp"] <= 0).any():
        raise ObservabilityContractError(
            "GENE_LENGTH", "gene_length_bp must be positive"
        )
    if (frame["expression_variability"] < 0).any():
        raise ObservabilityContractError(
            "EXPRESSION_VARIABILITY", "variability cannot be negative"
        )
    if (frame["cohort_detection_count"] < 0).any():
        raise ObservabilityContractError("COHORT_DETECTION", "count cannot be negative")
    if not np.allclose(
        frame["cohort_detection_count"], np.round(frame["cohort_detection_count"])
    ):
        raise ObservabilityContractError(
            "COHORT_DETECTION", "count must be integer-valued"
        )
    if set(frame["gene_biotype"]) != BIOLOGICAL_BIOTYPES:
        raise ObservabilityContractError(
            "BIOTYPE_ARMS", "both lncRNA and protein_coding are required"
        )

    frame["log_gene_length"] = np.log(frame["gene_length_bp"].astype(float))
    frame = frame.sort_values("ensembl_id", kind="mergesort").reset_index(drop=True)
    return frame


def _smd(left: pd.Series, right: pd.Series) -> float:
    left_values = left.to_numpy(dtype=float)
    right_values = right.to_numpy(dtype=float)
    variance = (np.var(left_values, ddof=1) + np.var(right_values, ddof=1)) / 2.0
    difference = float(np.mean(left_values) - np.mean(right_values))
    if not np.isfinite(variance) or variance < 0:
        return math.nan
    if variance == 0:
        return 0.0 if difference == 0 else math.copysign(math.inf, difference)
    return difference / math.sqrt(variance)


def _balance_rows(
    genes: pd.DataFrame, pairs: pd.DataFrame | None, stage: str
) -> list[dict[str, Any]]:
    if pairs is None:
        left = genes.loc[genes["gene_biotype"] == "lncRNA"].set_index("ensembl_id")
        right = genes.loc[genes["gene_biotype"] == "protein_coding"].set_index(
            "ensembl_id"
        )
    else:
        indexed = genes.set_index("ensembl_id")
        left = indexed.loc[pairs["lncrna_ensembl_id"]]
        right = indexed.loc[pairs["protein_coding_ensembl_id"]]
    rows = []
    for covariate in MATCH_COVARIATES:
        smd = _smd(left[covariate], right[covariate]) if len(left) else math.nan
        rows.append(
            {
                "stage": stage,
                "covariate": covariate,
                "lncrna_mean": float(left[covariate].mean()) if len(left) else math.nan,
                "protein_coding_mean": (
                    float(right[covariate].mean()) if len(right) else math.nan
                ),
                "smd": smd,
                "abs_smd": abs(smd),
            }
        )
    return rows


def match_genes(
    genes: pd.DataFrame,
) -> MatchResult:
    """Fit propensity scores and perform deterministic 1:1 matching."""

    frame = prepare_gene_covariates(genes)
    covariates = frame[list(MATCH_COVARIATES)].to_numpy(dtype=float)
    means = covariates.mean(axis=0)
    scales = covariates.std(axis=0, ddof=0)
    if np.any(~np.isfinite(scales)):
        raise ObservabilityContractError("COVARIATE_SCALE", "nonfinite covariate scale")
    # A constant covariate remains in balance reporting but contributes zero to
    # the propensity model. This is common when every gene is detected in all
    # input cohorts and is not a source-gate failure.
    scales[scales == 0] = 1.0
    design = (covariates - means) / scales
    outcome = (frame["gene_biotype"] == "lncRNA").astype(int).to_numpy()
    model = LogisticRegression(
        penalty="l2",
        solver="lbfgs",
        max_iter=10_000,
        tol=1e-10,
        random_state=MATCH_SEED,
    )
    model.fit(design, outcome)
    if int(model.n_iter_[0]) >= int(model.max_iter):
        raise ObservabilityContractError(
            "PROPENSITY_CONVERGENCE", f"reached max_iter={model.max_iter}"
        )
    propensity = np.clip(model.predict_proba(design)[:, 1], 1e-12, 1 - 1e-12)
    logits = np.log(propensity / (1 - propensity))
    logit_sd = float(np.std(logits, ddof=1))
    if not np.isfinite(logit_sd):
        raise ObservabilityContractError("PROPENSITY_SD", "logit SD is nonfinite")
    caliper = CALIPER_SD * logit_sd
    frame["propensity"] = propensity
    frame["propensity_logit"] = logits

    treated = frame.loc[frame["gene_biotype"] == "lncRNA"].copy()
    controls = frame.loc[frame["gene_biotype"] == "protein_coding"].copy()
    sorted_control_logits = np.sort(controls["propensity_logit"].to_numpy(dtype=float))
    treated["n_controls_in_caliper"] = [
        int(
            np.searchsorted(sorted_control_logits, value + caliper, side="right")
            - np.searchsorted(sorted_control_logits, value - caliper, side="left")
        )
        for value in treated["propensity_logit"].to_numpy(dtype=float)
    ]
    treated = treated.sort_values(
        ["n_controls_in_caliper", "ensembl_id"], kind="mergesort"
    )
    available = sorted(
        zip(
            controls["propensity_logit"].to_numpy(dtype=float),
            controls["ensembl_id"].astype(str),
        )
    )
    control_propensity = controls.set_index("ensembl_id")["propensity"].to_dict()
    pair_rows: list[dict[str, Any]] = []
    for row in treated.itertuples(index=False):
        if not available:
            break
        target = float(row.propensity_logit)
        position = bisect_left(available, (target, ""))
        candidate_positions = {
            index for index in (position - 1, position) if 0 <= index < len(available)
        }
        if not candidate_positions:
            continue
        best_index = min(
            candidate_positions,
            key=lambda index: (
                abs(available[index][0] - target),
                available[index][1],
            ),
        )
        control_logit, control_id = available[best_index]
        distance = abs(control_logit - target)
        if distance > caliper + 1e-12:
            continue
        available.pop(best_index)
        pair_rows.append(
            {
                "pair_id": f"pair:{len(pair_rows) + 1:06d}",
                "lncrna_ensembl_id": str(row.ensembl_id),
                "protein_coding_ensembl_id": control_id,
                "lncrna_propensity": float(row.propensity),
                "protein_coding_propensity": float(control_propensity[control_id]),
                "absolute_logit_distance": distance,
                "caliper": caliper,
            }
        )
    pairs = pd.DataFrame(
        pair_rows,
        columns=(
            "pair_id",
            "lncrna_ensembl_id",
            "protein_coding_ensembl_id",
            "lncrna_propensity",
            "protein_coding_propensity",
            "absolute_logit_distance",
            "caliper",
        ),
    )
    if (
        pairs["lncrna_ensembl_id"].duplicated().any()
        or pairs["protein_coding_ensembl_id"].duplicated().any()
    ):
        raise ObservabilityContractError(
            "MATCH_REUSE", "matching must be without replacement"
        )
    match_rate = len(pairs) / len(treated)
    balance = pd.DataFrame(
        _balance_rows(frame, None, "before") + _balance_rows(frame, pairs, "matched")
    )
    matched_smd = balance.loc[balance["stage"] == "matched", "abs_smd"]
    max_observed_smd = float(matched_smd.max()) if len(pairs) else math.inf
    failures = []
    if match_rate < MIN_MATCH_RATE:
        failures.append(f"match_rate={match_rate:.6f}<{MIN_MATCH_RATE:.6f}")
    if not np.isfinite(max_observed_smd) or max_observed_smd > MAX_ABS_SMD:
        failures.append(f"max_abs_smd={max_observed_smd:.6f}>{MAX_ABS_SMD:.6f}")
    gate_status = "pass" if not failures else "fail"
    return MatchResult(
        pairs=pairs,
        balance=balance,
        eligible_genes=frame,
        caliper=caliper,
        match_rate=match_rate,
        max_abs_matched_smd=max_observed_smd,
        gate_status=gate_status,
        gate_reason="all_matching_gates_passed" if not failures else ";".join(failures),
    )


def validate_coverage(coverage: pd.DataFrame) -> pd.DataFrame:
    """Validate the source-adapter observability contract and assay semantics."""

    _require_columns(coverage, COVERAGE_REQUIRED_COLUMNS, "coverage")
    if coverage.empty:
        raise ObservabilityContractError("EMPTY_COVERAGE", "coverage table is empty")
    frame = coverage.loc[:, COVERAGE_REQUIRED_COLUMNS].copy()
    for column in (
        "ensembl_id",
        "assay_id",
        "assay_family",
        "molecular_scope",
        "applicability_state",
    ):
        frame[column] = frame[column].astype("string").str.strip()
    if frame[["ensembl_id", "assay_id"]].isna().any().any():
        raise ObservabilityContractError("COVERAGE_KEY", "gene and assay are required")
    if frame.duplicated(["ensembl_id", "assay_id"]).any():
        raise ObservabilityContractError(
            "DUPLICATE_COVERAGE", "one row per gene and assay"
        )
    unexpected_families = sorted(set(frame["assay_family"]) - ASSAY_FAMILIES)
    unexpected_scopes = sorted(set(frame["molecular_scope"]) - MOLECULAR_SCOPES)
    unexpected_states = sorted(set(frame["applicability_state"]) - APPLICABILITY_STATES)
    if unexpected_families:
        raise ObservabilityContractError("ASSAY_FAMILY", str(unexpected_families))
    if unexpected_scopes:
        raise ObservabilityContractError("MOLECULAR_SCOPE", str(unexpected_scopes))
    if unexpected_states:
        raise ObservabilityContractError("APPLICABILITY", str(unexpected_states))
    frame["observable"] = _as_nullable_boolean(frame["observable"], "observable")
    should_be_null = ~frame["applicability_state"].isin(
        OBSERVABILITY_DENOMINATOR_STATES
    )
    if frame.loc[should_be_null, "observable"].notna().any():
        raise ObservabilityContractError(
            "OBSERVABILITY_WITHOUT_APPLICABILITY",
            "untestable, indeterminate, and not_applicable rows require null observable",
        )
    if frame.loc[~should_be_null, "observable"].isna().any():
        raise ObservabilityContractError(
            "MISSING_OBSERVABILITY",
            "applicable/source-dependent rows require observable",
        )
    if (
        (frame["assay_family"] == "proteomics")
        & (frame["molecular_scope"] != "protein_molecule")
    ).any():
        raise ObservabilityContractError(
            "PROTEOMICS_SCOPE", "proteomics must declare protein_molecule"
        )
    if (
        (frame["assay_family"] == "chromatin")
        & (frame["molecular_scope"] != "regulatory_locus")
    ).any():
        raise ObservabilityContractError(
            "CHROMATIN_SCOPE",
            "chromatin must declare regulatory_locus, not RNA observation",
        )
    return frame.sort_values(["assay_id", "ensembl_id"], kind="mergesort").reset_index(
        drop=True
    )


def _percentage(values: pd.Series) -> float:
    return 100.0 * float(values.astype(float).mean())


def _bootstrap_pair_difference(
    differences: np.ndarray, assay_id: str, replicates: int, seed: int
) -> tuple[float, float]:
    if replicates < 100:
        raise ObservabilityContractError(
            "BOOTSTRAP_REPLICATES", "at least 100 required"
        )
    rng = np.random.default_rng(_stable_assay_seed(assay_id, seed))
    n_pairs = len(differences)
    estimates = np.empty(replicates, dtype=float)
    chunk_size = max(1, min(256, 2_000_000 // n_pairs))
    for start in range(0, replicates, chunk_size):
        stop = min(start + chunk_size, replicates)
        draws = rng.integers(0, n_pairs, size=(stop - start, n_pairs))
        estimates[start:stop] = differences[draws].mean(axis=1) * 100.0
    lower, upper = np.quantile(estimates, [0.025, 0.975])
    return float(lower), float(upper)


def compare_observability(
    genes: pd.DataFrame,
    coverage: pd.DataFrame,
    *,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = MATCH_SEED,
) -> tuple[MatchResult, pd.DataFrame, pd.DataFrame]:
    """Return matching, assay comparisons, and complete status counts."""

    match = match_genes(genes)
    coverage_frame = validate_coverage(coverage)
    gene_identity = genes[["ensembl_id", "gene_biotype"]].copy()
    gene_identity["ensembl_id"] = (
        gene_identity["ensembl_id"].astype("string").str.strip()
    )
    gene_identity["gene_biotype"] = (
        gene_identity["gene_biotype"].astype("string").str.strip()
    )
    known_genes = set(gene_identity["ensembl_id"].astype(str))
    unknown = sorted(set(coverage_frame["ensembl_id"].astype(str)) - known_genes)
    if unknown:
        raise ObservabilityContractError(
            "COVERAGE_GENE", f"unknown genes: {unknown[:5]}"
        )
    assay_ids = set(coverage_frame["assay_id"].astype(str))
    expected_keys = {(gene, assay) for gene in known_genes for assay in assay_ids}
    observed_keys = set(
        coverage_frame[["ensembl_id", "assay_id"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if observed_keys != expected_keys:
        missing = sorted(expected_keys - observed_keys)
        raise ObservabilityContractError(
            "INCOMPLETE_COVERAGE_MATRIX",
            f"missing {len(missing)} gene-assay rows; first={missing[:3]}",
        )
    merged = coverage_frame.merge(
        gene_identity,
        on="ensembl_id",
        how="inner",
        validate="many_to_one",
    )
    lnc_proteomics = merged[
        (merged["gene_biotype"] == "lncRNA") & (merged["assay_family"] == "proteomics")
    ]
    if not lnc_proteomics.empty and set(lnc_proteomics["applicability_state"]) != {
        "not_applicable"
    }:
        raise ObservabilityContractError(
            "LNCRNA_PROTEOMICS",
            "every lncRNA proteomics row must be not_applicable",
        )

    status_counts = (
        merged.groupby(
            [
                "assay_id",
                "assay_family",
                "molecular_scope",
                "gene_biotype",
                "applicability_state",
            ],
            observed=True,
            dropna=False,
        )
        .size()
        .rename("n_genes")
        .reset_index()
        .sort_values(
            ["assay_id", "gene_biotype", "applicability_state"], kind="mergesort"
        )
        .reset_index(drop=True)
    )

    pair_long = pd.concat(
        [
            match.pairs[["pair_id", "lncrna_ensembl_id"]]
            .rename(columns={"lncrna_ensembl_id": "ensembl_id"})
            .assign(gene_biotype="lncRNA"),
            match.pairs[["pair_id", "protein_coding_ensembl_id"]]
            .rename(columns={"protein_coding_ensembl_id": "ensembl_id"})
            .assign(gene_biotype="protein_coding"),
        ],
        ignore_index=True,
    )
    comparison_rows = []
    for assay_id, assay in merged.groupby("assay_id", sort=True):
        families = set(assay["assay_family"])
        scopes = set(assay["molecular_scope"])
        if len(families) != 1 or len(scopes) != 1:
            raise ObservabilityContractError("ASSAY_METADATA", str(assay_id))
        family = next(iter(families))
        scope = next(iter(scopes))
        raw = assay.loc[
            assay["applicability_state"].isin(OBSERVABILITY_DENOMINATOR_STATES)
        ]
        raw_lnc = raw.loc[raw["gene_biotype"] == "lncRNA", "observable"]
        raw_pc = raw.loc[raw["gene_biotype"] == "protein_coding", "observable"]
        joined_pairs = pair_long.merge(
            assay[["ensembl_id", "gene_biotype", "applicability_state", "observable"]],
            on=["ensembl_id", "gene_biotype"],
            how="left",
            validate="many_to_one",
        )
        applicable = joined_pairs["applicability_state"].isin(
            OBSERVABILITY_DENOMINATOR_STATES
        )
        paired = joined_pairs.loc[applicable].pivot(
            index="pair_id", columns="gene_biotype", values="observable"
        )
        if {"lncRNA", "protein_coding"}.issubset(paired.columns):
            paired = paired.dropna(subset=["lncRNA", "protein_coding"])
        else:
            paired = paired.iloc[0:0]

        if len(raw_lnc) == 0:
            status = "not_applicable" if family == "proteomics" else "indeterminate"
            raw_lnc_percent = math.nan
            raw_pc_percent = _percentage(raw_pc) if len(raw_pc) else math.nan
            raw_difference = math.nan
            matched_lnc_percent = matched_pc_percent = matched_difference = math.nan
            ci_lower = ci_upper = math.nan
        elif len(raw_pc) == 0:
            status = "indeterminate"
            raw_lnc_percent = _percentage(raw_lnc)
            raw_pc_percent = math.nan
            raw_difference = raw_lnc_percent - raw_pc_percent
            matched_lnc_percent = matched_pc_percent = matched_difference = math.nan
            ci_lower = ci_upper = math.nan
        elif match.gate_status != "pass":
            status = "matching_gate_failed"
            raw_lnc_percent = _percentage(raw_lnc)
            raw_pc_percent = _percentage(raw_pc)
            raw_difference = raw_lnc_percent - raw_pc_percent
            matched_lnc_percent = matched_pc_percent = matched_difference = math.nan
            ci_lower = ci_upper = math.nan
        elif paired.empty:
            status = "indeterminate"
            raw_lnc_percent = _percentage(raw_lnc)
            raw_pc_percent = _percentage(raw_pc)
            raw_difference = raw_lnc_percent - raw_pc_percent
            matched_lnc_percent = matched_pc_percent = matched_difference = math.nan
            ci_lower = ci_upper = math.nan
        else:
            status = "matched_comparison"
            raw_lnc_percent = _percentage(raw_lnc)
            raw_pc_percent = _percentage(raw_pc)
            raw_difference = raw_lnc_percent - raw_pc_percent
            matched_lnc_percent = _percentage(paired["lncRNA"])
            matched_pc_percent = _percentage(paired["protein_coding"])
            matched_difference = matched_lnc_percent - matched_pc_percent
            differences = (
                paired["lncRNA"].astype(float).to_numpy()
                - paired["protein_coding"].astype(float).to_numpy()
            )
            ci_lower, ci_upper = _bootstrap_pair_difference(
                differences, str(assay_id), bootstrap_replicates, seed
            )
        comparison_rows.append(
            {
                "assay_id": str(assay_id),
                "assay_family": family,
                "molecular_scope": scope,
                "comparison_status": status,
                "n_raw_lncrna": len(raw_lnc),
                "n_raw_protein_coding": len(raw_pc),
                "n_matched_pairs_assay_eligible": len(paired),
                "raw_lncrna_percent": raw_lnc_percent,
                "raw_protein_coding_percent": raw_pc_percent,
                "raw_difference_percentage_points": raw_difference,
                "matched_lncrna_percent": matched_lnc_percent,
                "matched_protein_coding_percent": matched_pc_percent,
                "matched_difference_percentage_points": matched_difference,
                "bootstrap_ci_lower_percentage_points": ci_lower,
                "bootstrap_ci_upper_percentage_points": ci_upper,
                "bootstrap_replicates": bootstrap_replicates,
                "bootstrap_seed": _stable_assay_seed(str(assay_id), seed),
                "matching_gate_status": match.gate_status,
                "interpretation": (
                    "locus_observability_not_lncrna_molecule"
                    if family == "chromatin"
                    else "lncrna_molecule_not_applicable"
                    if family == "proteomics"
                    else "assay_native_observability"
                ),
            }
        )
    comparisons = (
        pd.DataFrame(comparison_rows)
        .sort_values("assay_id", kind="mergesort")
        .reset_index(drop=True)
    )
    return match, comparisons, status_counts


def matching_gate_row(result: MatchResult) -> dict[str, Any]:
    """Produce the machine-readable promotion gate for a candidate bundle."""

    return {
        "gate_id": "lncrna_protein_coding_observability_matching",
        "status": result.gate_status,
        "promotion_allowed": result.gate_status == "pass",
        "n_eligible_lncrna": int(
            (result.eligible_genes["gene_biotype"] == "lncRNA").sum()
        ),
        "n_matched_pairs": len(result.pairs),
        "match_rate": result.match_rate,
        "required_match_rate": MIN_MATCH_RATE,
        "max_abs_matched_smd": result.max_abs_matched_smd,
        "required_max_abs_smd": MAX_ABS_SMD,
        "caliper_logit": result.caliper,
        "caliper_definition": "0.2 times pooled propensity-logit sample SD",
        "matching_method": (
            "one-to-one deterministic nearest-neighbor without replacement; "
            "fewest eligible controls first; Ensembl tie-break"
        ),
        "reason": result.gate_reason,
    }


def contract_json() -> str:
    """Stable contract summary for candidate manifests and tests."""

    payload: Mapping[str, Any] = {
        "gene_required_columns": GENE_REQUIRED_COLUMNS,
        "coverage_required_columns": COVERAGE_REQUIRED_COLUMNS,
        "match_covariates": MATCH_COVARIATES,
        "gene_length_transform": "natural_log",
        "caliper_sd": CALIPER_SD,
        "min_match_rate": MIN_MATCH_RATE,
        "max_abs_smd": MAX_ABS_SMD,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "seed": MATCH_SEED,
        "observability_denominator_states": sorted(OBSERVABILITY_DENOMINATOR_STATES),
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))
