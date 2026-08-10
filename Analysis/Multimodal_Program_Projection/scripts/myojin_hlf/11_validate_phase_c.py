#!/usr/bin/env python3
"""Independently validate the completed noncanonical Myojin phase-C release."""

from __future__ import annotations

import argparse
import math
import platform
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.api as sm

from myojin_firewall_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    RELEASE_ID,
    atomic_write_text,
    parse_bool,
    read_tsv,
    require_within,
    sha256_file,
    stable_bundle_sha256,
    write_tsv,
)
from phase_c_stats import mechanical_figure_verdict


SPECIFICATION_SHA256 = (
    "1f99b7a7a95d9aebbbb911bc2317860f4a3fada0d1da2a4932d5c1ab190c4cdd"
)
RAW_SHA256 = "2663a811676a1abe2778375e6d28f6397eaded1f895d274f3c241733669c1c66"
EXPECTED_VERSIONS = {
    "python": "3.10",
    "numpy": "2.2.6",
    "pandas": "2.3.3",
    "scipy": "1.15.2",
    "statsmodels": "0.14.5",
}
CONTINUOUS_VARIANTS = {
    "primary",
    "guide_floor_5",
    "remove_all_known_hits",
    "remove_ACSL3",
    "missing_Chronos_indicator",
    "PA_D21_vs_D0",
    "absolute_PA_vs_vehicle_LFC",
}
BINARY_VARIANTS = {
    "primary_protective_hit",
    "guide_floor_5",
    "remove_all_known_hits",
    "remove_ACSL3",
    "missing_Chronos_indicator",
    "positive_selection_FDR_without_LFC_floor",
}
PROGRAM_FAMILIES = {
    ("primary", "original_weight"),
    ("primary", "equal_weight"),
    ("primary", "leave_highest_weight_out"),
    ("guide_floor_5", "original_weight"),
    ("remove_all_known_hits", "original_weight"),
    ("remove_ACSL3", "original_weight"),
}
EXPECTED_RELEASE_FILES = {
    "unseal_record.tsv",
    "UNSEALED",
    "phase_c_preflight_audit.tsv",
    "phase_c_code_manifest.tsv",
    "PHASE_C_CODE_READY",
    "screen_results.tsv",
    "class_effects.tsv",
    "class_protective_hits.tsv",
    "program_effects.tsv",
    "permutation_audit.tsv",
    "matched_null_audit.tsv",
    "sensitivity.tsv",
    "fig5_verdict.tsv",
}


class ValidationError(RuntimeError):
    pass


def exact_bh(values: list[float]) -> list[float]:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p, kind="stable")
    ranked = p[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(p) / np.arange(1, len(p) + 1))[::-1]
    )[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result.tolist()


def close(left: object, right: object, tolerance: float = 1e-10) -> bool:
    return math.isclose(float(left), float(right), rel_tol=tolerance, abs_tol=tolerance)


def add_check(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append(
        {
            "check": name,
            "passed": str(passed).upper(),
            "detail": detail,
            "specification_sha256": SPECIFICATION_SHA256,
            "release_id": RELEASE_ID,
        }
    )
    if not passed:
        raise ValidationError(f"{name}::{detail}")


def one_row(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    if len(rows) != 1:
        raise ValidationError(f"expected_one_row::{path.name}::{len(rows)}")
    return rows[0]


def validate_code_freeze(outdir: Path, checks: list[dict[str, object]]) -> str:
    marker = one_row(outdir / "PHASE_C_CODE_READY")
    manifest_path = outdir / "phase_c_code_manifest.tsv"
    add_check(
        checks,
        "code_manifest_hash_unchanged",
        sha256_file(manifest_path) == marker["manifest_sha256"],
        marker["manifest_sha256"],
    )
    manifest = read_tsv(manifest_path)
    code_paths = []
    files_ok = True
    environment = {}
    for row in manifest:
        if row["role"] == "environment_version":
            environment[row["relative_path"].split("::", 1)[1]] = row["sha256"]
            continue
        path = PROJECT_ROOT / row["relative_path"]
        files_ok &= (
            path.is_file()
            and str(path.stat().st_size) == row["bytes"]
            and sha256_file(path) == row["sha256"]
        )
        if row["role"] == "phase_c_code":
            code_paths.append(path)
    add_check(
        checks, "frozen_code_and_inputs_unchanged", files_ok, f"rows={len(manifest)}"
    )
    bundle_hash, _ = stable_bundle_sha256(code_paths, PROJECT_ROOT)
    add_check(
        checks,
        "phase_c_code_bundle_unchanged",
        bundle_hash == marker["code_bundle_sha256"],
        bundle_hash,
    )
    observed_environment = {
        "python": ".".join(platform.python_version().split(".")[:2]),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "statsmodels": statsmodels.__version__,
    }
    add_check(
        checks,
        "execution_environment_unchanged",
        observed_environment == environment == EXPECTED_VERSIONS,
        str(observed_environment),
    )
    return bundle_hash


def validate_screen(
    outdir: Path, checks: list[dict[str, object]]
) -> tuple[list[dict[str, str]], pd.DataFrame]:
    screen = read_tsv(outdir / "screen_results.tsv")
    source_rows = [int(row["source_row"]) for row in screen]
    add_check(
        checks,
        "full_screen_universe_released_once",
        len(screen) == 18_343
        and source_rows == list(range(4, 18_347))
        and all(row["outcome_values_read"] == "TRUE" for row in screen),
        f"rows={len(screen)}",
    )
    mechanical_hits = True
    direction_policy = True
    for row in screen:
        lfc_text = row["PA_vs_vehicle_positive_selection_LFC"]
        fdr_text = row["PA_vs_vehicle_positive_selection_FDR"]
        lfc = None if lfc_text == "" else float(lfc_text)
        fdr = None if fdr_text == "" else float(fdr_text)
        protective = lfc is not None and fdr is not None and fdr < 0.05 and lfc >= 0.25
        fdr_only = fdr is not None and fdr < 0.05
        mechanical_hits &= (
            parse_bool(row["protective_hit_FDR_lt_0.05_LFC_ge_0.25"]) == protective
            and parse_bool(row["positive_selection_FDR_hit_without_LFC_floor"])
            == fdr_only
        )
        expected_policy = (
            "no_direct_result"
            if lfc is None
            else "continuous_negative_only_not_source_significant_sensitizing"
            if lfc < 0
            else "positive_selection_direction_continuous"
        )
        direction_policy &= row["source_significance_interpretation"] == expected_policy
    add_check(
        checks, "protective_hits_mechanical", mechanical_hits, "FDR<0.05 and LFC>=0.25"
    )
    add_check(
        checks,
        "negative_LFC_not_promoted_to_source_significant_sensitizing",
        direction_policy,
        "continuous direction only",
    )
    primary_rows = [row for row in screen if parse_bool(row["primary_eligible"])]
    bins_complete = all(
        row["analysis_guide_bin"] != ""
        and row["analysis_expression_quintile"] != ""
        and row["analysis_chronos_quintile"] != ""
        for row in primary_rows
    )
    add_check(
        checks,
        "frozen_primary_universe_and_bins",
        len(primary_rows) == 7_762 and bins_complete,
        f"eligible={len(primary_rows)}",
    )
    frame = pd.DataFrame(primary_rows)
    for column in (
        "guide_count_min",
        "log1p_HLF_TPM",
        "HLF_Chronos",
        "PA_vs_vehicle_positive_selection_LFC",
    ):
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    return screen, frame


def validate_observed_primary_model(
    frame: pd.DataFrame,
    class_rows: list[dict[str, str]],
    checks: list[dict[str, object]],
) -> None:
    y = frame["PA_vs_vehicle_positive_selection_LFC"].to_numpy(dtype=float)
    class_names = frame["primary_evidence_class"].to_numpy()
    full = np.column_stack(
        [
            np.ones(len(frame)),
            (class_names == "genetic_only").astype(float),
            (class_names == "disease_state_only").astype(float),
            (class_names == "convergent").astype(float),
            frame["log1p_HLF_TPM"].to_numpy(dtype=float),
            frame["guide_count_min"].to_numpy(dtype=float),
            frame["HLF_Chronos"].to_numpy(dtype=float),
        ]
    )
    reduced = full[:, [0, 4, 5, 6]]
    fit = sm.OLS(y, full).fit(cov_type="HC3")
    rss_full = float(np.sum((y - full @ np.linalg.lstsq(full, y, rcond=None)[0]) ** 2))
    rss_reduced = float(
        np.sum((y - reduced @ np.linalg.lstsq(reduced, y, rcond=None)[0]) ** 2)
    )
    residual_df = len(y) - np.linalg.matrix_rank(full)
    residual_sd = math.sqrt(rss_full / residual_df)
    observed_f = ((rss_reduced - rss_full) / 3) / (rss_full / residual_df)
    by_contrast = {
        row["contrast"]: row
        for row in class_rows
        if row["analysis_variant"] == "primary"
    }
    comparison = close(by_contrast["omnibus_evidence_class"]["statistic"], observed_f)
    for index, contrast in enumerate(
        (
            "genetic_only_vs_neither",
            "disease_state_only_vs_neither",
            "convergent_vs_neither",
        ),
        start=1,
    ):
        row = by_contrast[contrast]
        comparison &= (
            close(row["estimate"], fit.params[index])
            and close(row["HC3_SE"], fit.bse[index])
            and close(row["standardized_effect"], fit.params[index] / residual_sd)
            and close(row["residual_SD"], residual_sd)
        )
    add_check(
        checks,
        "primary_observed_model_independently_rederived",
        comparison,
        "OLS/HC3/F/SD",
    )


def validate_class_results(
    outdir: Path,
    primary: pd.DataFrame,
    checks: list[dict[str, object]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    continuous = read_tsv(outdir / "class_effects.tsv")
    binary = read_tsv(outdir / "class_protective_hits.tsv")
    continuous_ok = (
        len(continuous) == 28
        and {row["analysis_variant"] for row in continuous} == CONTINUOUS_VARIANTS
        and all(int(row["accepted_draws"]) == 100_000 for row in continuous)
        and all(int(row["attempted_draws"]) >= 100_000 for row in continuous)
    )
    binary_ok = (
        len(binary) == 18
        and {row["analysis_variant"] for row in binary} == BINARY_VARIANTS
        and all(int(row["accepted_draws"]) == 100_000 for row in binary)
    )
    add_check(
        checks, "continuous_families_complete", continuous_ok, f"rows={len(continuous)}"
    )
    add_check(checks, "binary_families_complete", binary_ok, f"rows={len(binary)}")
    bh_ok = True
    for variant in CONTINUOUS_VARIANTS:
        pair = [
            row
            for row in continuous
            if row["analysis_variant"] == variant
            and row["contrast"] != "omnibus_evidence_class"
        ]
        pair.sort(key=lambda row: row["contrast"])
        expected = exact_bh([float(row["permutation_p"]) for row in pair])
        bh_ok &= all(close(row["BH_q"], q) for row, q in zip(pair, expected))
    for variant in BINARY_VARIANTS:
        pair = sorted(
            [row for row in binary if row["analysis_variant"] == variant],
            key=lambda row: row["contrast"],
        )
        expected = exact_bh([float(row["permutation_p"]) for row in pair])
        bh_ok &= all(close(row["BH_q"], q) for row, q in zip(pair, expected))
    add_check(checks, "class_BH_independently_rederived", bh_ok, "three-test families")
    validate_observed_primary_model(primary, continuous, checks)
    primary_n = {
        int(row["n_total"])
        for row in continuous
        if row["analysis_variant"] == "primary"
    }
    add_check(
        checks,
        "primary_class_n_is_frozen_universe",
        primary_n == {7_762},
        str(primary_n),
    )
    return continuous, binary


def validate_program_results(
    outdir: Path,
    primary: pd.DataFrame,
    checks: list[dict[str, object]],
) -> list[dict[str, str]]:
    rows = read_tsv(outdir / "program_effects.tsv")
    family_counts = defaultdict(int)
    for row in rows:
        family_counts[(row["analysis_variant"], row["weight_mode"])] += 1
    structure_ok = (
        len(rows) == 12
        and set(family_counts) == PROGRAM_FAMILIES
        and all(count == 2 for count in family_counts.values())
        and all(
            (
                int(row["draws"]) == 10_000
                if row["testable"] == "TRUE"
                else int(row["draws"]) == 0
            )
            for row in rows
        )
    )
    add_check(checks, "program_families_complete", structure_ok, f"rows={len(rows)}")
    bh_ok = True
    for family in PROGRAM_FAMILIES:
        family_rows = [
            row
            for row in rows
            if (row["analysis_variant"], row["weight_mode"]) == family
            and row["testable"] == "TRUE"
        ]
        if family_rows:
            family_rows.sort(key=lambda row: row["program_uid"])
            expected = exact_bh([float(row["empirical_p"]) for row in family_rows])
            bh_ok &= all(close(row["BH_q"], q) for row, q in zip(family_rows, expected))
    add_check(
        checks, "program_BH_independently_rederived", bh_ok, "complete testable family"
    )

    prediction = read_tsv(outdir / "prediction_manifest.tsv")
    memberships: dict[str, list[dict[str, str]]] = defaultdict(list)
    effects: dict[str, set[float]] = defaultdict(set)
    for row in prediction:
        if row["object_type"] == "hepatocyte_program_gene":
            memberships[row["object_id"]].append(row)
            effects[row["object_id"]].add(float(row["frozen_source_effect"]))
    values = primary["PA_vs_vehicle_positive_selection_LFC"].to_numpy(dtype=float)
    z_by_symbol = dict(
        zip(
            primary["gene_symbol"],
            (values - np.mean(values)) / np.std(values, ddof=1),
        )
    )
    observed_ok = True
    for uid, members in memberships.items():
        eligible = [row for row in members if row["gene_symbol"] in z_by_symbol]
        direction = 1 if next(iter(effects[uid])) > 0 else -1
        for mode in ("original_weight", "equal_weight", "leave_highest_weight_out"):
            use = eligible.copy()
            if mode == "leave_highest_weight_out" and use:
                drop = sorted(
                    use,
                    key=lambda row: (-float(row["frozen_weight"]), row["gene_symbol"]),
                )[0]
                use.remove(drop)
            weights = np.asarray(
                [
                    1.0 if mode == "equal_weight" else float(row["frozen_weight"])
                    for row in use
                ]
            )
            weights /= weights.sum()
            observed = direction * float(
                np.sum(
                    weights
                    * np.asarray([z_by_symbol[row["gene_symbol"]] for row in use])
                )
            )
            result = one_row_from(
                rows,
                analysis_variant="primary",
                weight_mode=mode,
                program_uid=uid,
            )
            observed_ok &= int(result["n_eligible_genes"]) == len(use)
            if result["testable"] == "TRUE":
                observed_ok &= close(result["observed_signed_score"], observed)
            else:
                observed_ok &= result["observed_signed_score"] == ""
    add_check(
        checks,
        "primary_program_scores_independently_rederived",
        observed_ok,
        "weights/z/direction",
    )
    return rows


def one_row_from(rows: list[dict[str, str]], **criteria: str) -> dict[str, str]:
    selected = [
        row for row in rows if all(row[key] == value for key, value in criteria.items())
    ]
    if len(selected) != 1:
        raise ValidationError(
            f"Expected one row for {criteria}; observed {len(selected)}"
        )
    return selected[0]


def validate_sensitivities(outdir: Path, checks: list[dict[str, object]]) -> None:
    rows = read_tsv(outdir / "sensitivity.tsv")
    observed_cont = {
        row["analysis_variant"]
        for row in rows
        if row["analysis_family"] == "class_continuous"
    }
    observed_binary = {
        row["analysis_variant"]
        for row in rows
        if row["analysis_family"] == "class_binary"
    }
    observed_program = {
        (row["analysis_variant"], row["outcome"])
        for row in rows
        if row["analysis_family"] == "program_matched_null"
    }
    valid = (
        len(rows) == 19
        and observed_cont == CONTINUOUS_VARIANTS
        and observed_binary == BINARY_VARIANTS
        and observed_program == PROGRAM_FAMILIES
        and all(row["status"].startswith("completed") for row in rows)
    )
    add_check(checks, "mandatory_sensitivities_complete", valid, f"rows={len(rows)}")


def validate_release_manifest(outdir: Path, checks: list[dict[str, object]]) -> None:
    manifest_path = outdir / "phase_c_release_manifest.tsv"
    rows = read_tsv(manifest_path)
    names = {row["relative_path"] for row in rows}
    valid = len(rows) == len(EXPECTED_RELEASE_FILES) and names == EXPECTED_RELEASE_FILES
    for row in rows:
        path = outdir / row["relative_path"]
        valid &= (
            path.is_file()
            and path.stat().st_size == int(row["bytes"])
            and sha256_file(path) == row["sha256"]
            and sum(1 for _ in path.open(encoding="utf-8")) == int(row["line_count"])
            and row["canonical"] == "FALSE"
            and row["release_id"] == RELEASE_ID
        )
    complete = one_row(outdir / "PHASE_C_COMPLETE")
    valid &= (
        complete["release_manifest_sha256"] == sha256_file(manifest_path)
        and complete["biological_positivity_required"] == "FALSE"
        and complete["external_outcomes_read"] == "TRUE"
    )
    add_check(
        checks,
        "release_manifest_and_complete_marker_rederived",
        valid,
        f"files={len(rows)}",
    )


def run_validation(outdir: Path, checks: list[dict[str, object]]) -> str:
    required = [
        "UNSEALED",
        "PHASE_C_COMPLETE",
        "phase_c_code_manifest.tsv",
        "PHASE_C_CODE_READY",
        "unseal_record.tsv",
        "screen_results.tsv",
        "class_effects.tsv",
        "class_protective_hits.tsv",
        "program_effects.tsv",
        "permutation_audit.tsv",
        "matched_null_audit.tsv",
        "sensitivity.tsv",
        "fig5_verdict.tsv",
        "phase_c_release_manifest.tsv",
    ]
    missing = [name for name in required if not (outdir / name).is_file()]
    add_check(checks, "all_phase_c_outputs_present", not missing, f"missing={missing}")
    add_check(
        checks,
        "no_terminal_execution_failure",
        not (outdir / "phase_c_terminal_blocker.tsv").exists()
        and not (outdir / "phase_c_terminal_failure.tsv").exists(),
        "no blocker/failure marker",
    )
    code_hash = validate_code_freeze(outdir, checks)
    unseal = one_row(outdir / "unseal_record.tsv")
    unsealed = one_row(outdir / "UNSEALED")
    add_check(
        checks,
        "one_pass_unseal_record_valid",
        unseal["specification_sha256"] == SPECIFICATION_SHA256
        and unseal["raw_outcome_source_sha256"] == RAW_SHA256
        and unseal["one_pass_release_id"] == RELEASE_ID
        and unsealed["phase_c_code_bundle_sha256"] == code_hash
        and unsealed["external_outcomes_read"] == "TRUE",
        "one immutable unseal row",
    )
    preflight = read_tsv(outdir / "phase_c_preflight_audit.tsv")
    add_check(
        checks,
        "preflight_all_passed_before_unseal",
        bool(preflight)
        and all(
            row["passed"] == "TRUE" and row["external_outcomes_read"] == "FALSE"
            for row in preflight
        ),
        f"checks={len(preflight)}",
    )
    screen, primary = validate_screen(outdir, checks)
    continuous, _ = validate_class_results(outdir, primary, checks)
    programs = validate_program_results(outdir, primary, checks)
    validate_sensitivities(outdir, checks)

    permutation = read_tsv(outdir / "permutation_audit.tsv")
    matched = read_tsv(outdir / "matched_null_audit.tsv")
    nulls_ok = (
        len(permutation) == 25
        and len(matched) == 12
        and all(int(row["accepted_draws"]) == 100_000 for row in permutation)
        and all(
            (
                int(row["draws"]) == 10_000
                if row["testable"] == "TRUE"
                else int(row["draws"]) == 0
            )
            for row in matched
        )
    )
    add_check(
        checks,
        "all_prespecified_null_audits_present",
        nulls_ok,
        "25 permutation; 12 matched",
    )

    verdict = one_row(outdir / "fig5_verdict.tsv")
    rederived = mechanical_figure_verdict(True, continuous, programs)
    verdict_ok = all(
        verdict.get(key, "") == str(value) for key, value in rederived.items()
    )
    verdict_ok &= (
        verdict["biological_positivity_required_for_software_acceptance"] == "FALSE"
    )
    add_check(
        checks,
        "figure_gate_mechanically_rederived",
        verdict_ok,
        "zero discoveries remain valid software output",
    )
    validate_release_manifest(outdir, checks)
    return code_hash


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=CANDIDATE_ROOT)
    args = parser.parse_args()
    outdir = require_within(args.outdir, CANDIDATE_ROOT)
    checks: list[dict[str, object]] = []
    try:
        code_hash = run_validation(outdir, checks)
    except Exception as error:
        checks.append(
            {
                "check": "phase_c_validation_terminal",
                "passed": "FALSE",
                "detail": str(error),
                "specification_sha256": SPECIFICATION_SHA256,
                "release_id": RELEASE_ID,
            }
        )
        write_tsv(
            outdir / "phase_c_validation_status.tsv",
            checks,
            ["check", "passed", "detail", "specification_sha256", "release_id"],
        )
        atomic_write_text(
            outdir / "PHASE_C_VALIDATION_FAILED",
            "release_id\tstatus\tspecification_sha256\texternal_outcomes_read\n"
            f"{RELEASE_ID}\tterminal_validation_failure\t{SPECIFICATION_SHA256}\tTRUE\n",
        )
        raise
    write_tsv(
        outdir / "phase_c_validation_status.tsv",
        checks,
        ["check", "passed", "detail", "specification_sha256", "release_id"],
    )
    atomic_write_text(
        outdir / "PHASE_C_VALIDATED",
        "release_id\tstatus\tspecification_sha256\tphase_c_code_bundle_sha256\t"
        "validation_sha256\texternal_outcomes_read\tbiological_positivity_required\n"
        f"{RELEASE_ID}\tphase_c_release_validated\t{SPECIFICATION_SHA256}\t{code_hash}\t"
        f"{sha256_file(outdir / 'phase_c_validation_status.tsv')}\tTRUE\tFALSE\n",
    )
    print(
        "PHASE_C_VALIDATION_PASS "
        f"checks={len(checks)} specification_sha256={SPECIFICATION_SHA256} "
        "biological_positivity_not_required=TRUE"
    )


if __name__ == "__main__":
    main()
