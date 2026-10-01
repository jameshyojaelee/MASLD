#!/usr/bin/env python3
"""Descriptive source-quality and precision strata for fixed predictions.

Operational MAF/R2 bins were fixed after aggregate development results were
inspected, before this report. They are not original preregistration or
biological thresholds. No model is refitted and no weak-effect row is removed.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[3]
LABELS = ROOT/"GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz"
SOURCE_SE = ROOT/"GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/endpoints-revised-21772544/caqtl_source_leads.tsv.gz"
POPULATIONS = ("native_full_existing_leads", "fixed_representation_common_rows", "adaptation_single_split", "specialist_eligible_common_rows")
BIN_NAMES = {
    "imputation_r2": ("R2_lt0.3", "R2_0.3_to_lt0.8", "R2_ge0.8", "missing_or_invalid"),
    "maf": ("MAF_lt0.01", "MAF_0.01_to_lt0.05", "MAF_ge0.05", "missing_or_invalid"),
    "source_se": ("SE_training_Q1", "SE_training_Q2", "SE_training_Q3", "SE_training_Q4", "missing_or_unverified_SE"),
}


def quality_bins(values, name):
    x = np.asarray(values, dtype=float)
    result = np.full(len(x), "missing_or_invalid", dtype=object)
    a, b, ceiling = (0.3, 0.8, 1.0) if name == "imputation_r2" else (0.01, 0.05, 0.5)
    okay = np.isfinite(x) & (x >= 0) & (x <= ceiling)
    result[okay & (x < a)] = BIN_NAMES[name][0]
    result[okay & (x >= a) & (x < b)] = BIN_NAMES[name][1]
    result[okay & (x >= b)] = BIN_NAMES[name][2]
    return result


def add_source_se(labels, path):
    labels = labels.copy()
    labels["source_se"] = np.nan
    labels["source_se_status"] = "source_not_available"
    if not path.is_file():
        return labels, dict(status="not_run_source_missing", path=str(path))
    source = pd.read_csv(path, sep="\t", usecols=["variant_id", "target_id", "beta_nominal", "varbeta", "effect_se", "uncertainty_state"])
    duplicated = source.duplicated(["variant_id", "target_id"], keep=False)
    unique = source.loc[~duplicated].rename(columns={"variant_id": "lead_variant_id", "target_id": "peak_id", "effect_se": "deposited_se"})
    joined = labels.drop(columns=["source_se", "source_se_status"]).merge(unique, on=["lead_variant_id", "peak_id"], how="left", validate="one_to_one", indicator=True)
    # SEs stay in the same association-beta units. No p-value inversion,
    # population-count inference, allele sign change, or target substitution.
    same_beta = np.isclose(joined.beta_source, joined.beta_nominal, rtol=1e-5, atol=1e-7)
    native = joined.uncertainty_state.eq("sqrt_source_nominal_varbeta_matched_by_variant_and_peak")
    positive = np.isfinite(joined.varbeta) & (joined.varbeta > 0)
    rooted = np.sqrt(joined.varbeta.where(positive))
    same_se = np.isclose(rooted, joined.deposited_se, rtol=1e-8, atol=1e-12)
    eligible = same_beta & native & positive & same_se & joined._merge.eq("both")
    joined["source_se"] = rooted.where(eligible)
    joined["source_se_status"] = np.select(
        [joined._merge.ne("both"), ~same_beta, ~native, ~positive, ~same_se],
        ["missing_or_nonunique_exact_variant_peak", "source_beta_mismatch", "uncertainty_definition_unverified", "nonpositive_or_missing_varbeta", "deposited_SE_disagrees_with_varbeta"],
        default="verified_exact_variant_peak_native_varbeta")
    return joined.drop(columns="_merge"), dict(status="exact_variant_peak_join_checked", path=str(path),
        source_rows=len(source), nonunique_source_rows=int(duplicated.sum()), admitted_se=int(eligible.sum()),
        reasons=joined.source_se_status.value_counts().to_dict(),
        units="same_source_association_beta_units; sign_invariant_standard_error")


def se_strata(source, adaptation):
    """Reference cutoffs use source annotations in permitted training folds.

    These are source-precision quartiles, not model-specific training exposure.
    Equal cutoff values may leave empty strata, which are retained explicitly.
    """
    result = pd.Series("missing_or_unverified_SE", index=source.index, dtype=object)
    receipts = []
    for held in ([1] if adaptation else range(5)):
        train = ~source.heldout_fold.isin([0, 1]) if adaptation else source.heldout_fold.ne(held)
        train_values = source.loc[train, "source_se"].dropna().to_numpy()
        test = source.heldout_fold.eq(held) & source.source_se.notna()
        if not len(train_values):
            receipts.append(dict(held_fold=held, training_se_n=0, status="no_verified_training_SE"))
            continue
        cuts = np.quantile(train_values, [0.25, 0.5, 0.75])
        code = np.searchsorted(cuts, source.loc[test, "source_se"], side="right")
        result.loc[test] = [BIN_NAMES["source_se"][i] for i in code]
        receipts.append(dict(held_fold=held, training_se_n=len(train_values), q25=cuts[0], q50=cuts[1], q75=cuts[2],
            training_folds="2,3,4" if adaptation else ",".join(str(f) for f in range(5) if f != held),
            status="training_source_annotations_only"))
    return result, receipts


def metrics(y, predicted):
    if not len(y):
        return np.nan, np.nan, "empty_stratum"
    rmse = float(np.sqrt(np.mean((predicted-y)**2)))
    if len(y) < 3:
        return rmse, np.nan, "fewer_than_three_variants"
    if np.ptp(y) == 0 or np.ptp(predicted) == 0:
        return rmse, np.nan, "constant_measurement_or_prediction"
    return rmse, float(spearmanr(y, predicted).statistic), "defined"


def self_tests():
    assert list(quality_bins([0, .299, .3, .799, .8, 1, np.nan, -1, 1.1], "imputation_r2")) == [
        "R2_lt0.3", "R2_lt0.3", "R2_0.3_to_lt0.8", "R2_0.3_to_lt0.8", "R2_ge0.8", "R2_ge0.8", *["missing_or_invalid"]*3]
    assert list(quality_bins([.009, .01, .049, .05, .5, .6], "maf")) == [
        "MAF_lt0.01", "MAF_0.01_to_lt0.05", "MAF_0.01_to_lt0.05", "MAF_ge0.05", "MAF_ge0.05", "missing_or_invalid"]
    test = pd.DataFrame(dict(heldout_fold=np.repeat(np.arange(5), 5), source_se=np.arange(1, 26, dtype=float)))
    _, before = se_strata(test, True)
    test.loc[test.heldout_fold.isin([0, 1]), "source_se"] *= 1000
    _, after = se_strata(test, True)
    assert before == after, "Held-fold precision altered training cutoffs"
    assert metrics(np.array([0., .00001, -.00001]), np.zeros(3))[0] > 0
    assert metrics(np.zeros(3), np.zeros(3))[2] == "constant_measurement_or_prediction"


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run scientific checks and workload on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    self_tests()
    labels = pd.read_csv(LABELS, sep="\t", usecols=["lead_variant_id", "peak_id", "beta_source", "beta_alt", "maf", "imputation_r2", "heldout_fold", "block_1mb"])
    assert len(labels) == 32322 and not labels.lead_variant_id.duplicated().any()
    labels, se_receipt = add_source_se(labels, args.source_se)
    labels["key"] = labels.lead_variant_id.str.removeprefix("chr")
    labels = labels.set_index("key", drop=False)
    for quality in ("imputation_r2", "maf"):
        labels[quality+"_stratum"] = quality_bins(labels[quality], quality)
    missingness = []
    outputs, cutoffs, population_receipts = [], [], []
    source_hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (LABELS, args.source_se) if p.is_file()}
    for population in POPULATIONS:
        path = args.final/(population+"_matched_predictions.tsv.gz")
        performance_path = args.final/(population+"_performance.tsv")
        frame = pd.read_csv(path, sep="\t").set_index("key")
        arms = pd.read_csv(performance_path, sep="\t").arm.tolist()
        incomplete_path = args.final/(population+"_incomplete_arms.json")
        incomplete = json.loads(incomplete_path.read_text())
        failed_arms = [row["arm"] for row in incomplete]
        assert not frame.index.duplicated().any() and len(arms) == len(set(arms))
        source = labels.loc[frame.index]
        np.testing.assert_allclose(frame.beta_alt, source.beta_alt, atol=1e-7, rtol=1e-6)
        np.testing.assert_array_equal(frame.heldout_fold, source.heldout_fold)
        np.testing.assert_array_equal(frame.block_1mb, source.block_1mb)
        assert np.isfinite(frame[["beta_alt", *arms]].to_numpy(dtype=float)).all()
        adaptation = population == "adaptation_single_split"
        eligible = labels.loc[labels.heldout_fold.eq(1)] if adaptation else labels
        se_codes, cuts = se_strata(labels, adaptation)
        labels["source_se_stratum"] = se_codes
        source = labels.loc[frame.index]
        eligible = labels.loc[eligible.index]
        cutoffs.extend([dict(population=population, **r) for r in cuts])
        for quality, bins in BIN_NAMES.items():
            for label in bins:
                source_mask = eligible[quality+"_stratum"].eq(label)
                matched = source[quality+"_stratum"].eq(label)
                n_source, n_matched = int(source_mask.sum()), int(matched.sum())
                missingness.append(dict(population=population, source_variable=quality, stratum=label,
                    source_eligible_variants=n_source, matched_variants=n_matched,
                    unavailable_prediction_variants=n_source-n_matched,
                    source_eligible_locus_blocks=eligible.loc[source_mask, "block_1mb"].nunique(),
                    matched_locus_blocks=source.loc[matched, "block_1mb"].nunique()))
                y = frame.loc[matched, "beta_alt"].to_numpy(dtype=float)
                for arm in [*arms, *failed_arms]:
                    prediction = frame.loc[matched, arm].to_numpy(dtype=float)
                    available = int(np.isfinite(prediction).sum())
                    error, rho, rank_state = (np.nan, np.nan, "incomplete_fit_no_subset_performance") if arm in failed_arms else metrics(y, prediction)
                    native_se = source.loc[matched, "source_se"].dropna()
                    outputs.append(dict(population=population, model=arm, source_variable=quality, stratum=label,
                        matched_population_variants=len(y), prediction_available_variants=available,
                        matched_population_locus_blocks=source.loc[matched, "block_1mb"].nunique(),
                        source_eligible_variants=n_source, coverage_fraction=available/n_source if n_source else np.nan,
                        model_state="incomplete_fit" if arm in failed_arms else "complete_matched_fit",
                        RMSE_source_beta_units=error, signed_spearman=rho, correlation_state=rank_state,
                        source_SE_available_variants=len(native_se), source_SE_median=native_se.median(),
                        effective_participant_n="unknown_source_variant_specific_n_not_available",
                        evidence_state="descriptive_development;no_stratum_selection_or_inferential_test"))
            assert sum(r["matched_variants"] for r in missingness if r["population"] == population and r["source_variable"] == quality) == len(frame)
        source_hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        source_hashes[str(performance_path)] = hashlib.sha256(performance_path.read_bytes()).hexdigest()
        source_hashes[str(incomplete_path)] = hashlib.sha256(incomplete_path.read_bytes()).hexdigest()
        population_receipts.append(dict(population=population, source_eligible_n=len(eligible), matched_n=len(frame), complete_models=len(arms), incomplete_fits=incomplete))
    pd.DataFrame(outputs).to_csv(args.out/"stratified_performance.tsv", sep="\t", index=False)
    pd.DataFrame(missingness).to_csv(args.out/"source_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(cutoffs).to_csv(args.out/"source_SE_training_cutoffs.tsv", sep="\t", index=False)
    labels[["lead_variant_id", "peak_id", "imputation_r2", "imputation_r2_stratum", "maf", "maf_stratum", "source_se", "source_se_status"]].to_csv(args.out/"source_quality_annotations.tsv.gz", sep="\t", index=False)
    report = dict(status="pass", invariants="bin_boundaries_missingness_weak_rows_and_training_only_SE_cutoffs_passed",
        annotation_source=str(LABELS), source_SE=se_receipt, populations=population_receipts,
        bin_timing="operational_cuts_fixed_after_aggregate_development_inspection_before_this_stratified_report; not_original_preregistration",
        bin_meaning="source_quality_descriptions_not_biological_thresholds",
        source_SE_strata="quartiles_of_verified_SE_in_permitted_source_training_folds; native_cutoffs_recorded; ties_use_upper_bin",
        participant_n="unknown_effective_participants_per_variant; variant_and_locus_counts_are_not_participant_n",
        outcomes="stored_prediction_and_C2_beta_unchanged; no_allele_sign_reinterpretation",
        limits=["existing_significance_selected_deduplicated_Currin_leads_not_all_tested_variants",
                "descriptive_strata_no_bootstrap_p_values_or_multiple_testing_claim",
                "no_refitting_or_choice_of_favorable_stratum; all_nonmissing_and_missing_strata_reported",
                "native_1Mb_and_short_window_coverage_and_legacy_training_exposure_differ",
                "SE_stratum_reference_source_training_folds_not_model_specific_training_coverage",
                "low_MAF_or_imputation_bins_may_be_empty_or_sparse_due_to_historical_source_eligibility"],
        protected_outcomes_read=False, source_sha256=source_hashes)
    (args.out/"analysis.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final", type=Path, required=True)
    parser.add_argument("--source-se", type=Path, default=SOURCE_SE)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
