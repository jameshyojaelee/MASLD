#!/usr/bin/env python3
"""Independent saved-prediction comparison of the twenty Gaussian-loss runs."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from review_adaptation_precision import ADAPTERS, CONTROLS, bh, contrasts, metrics, paired_draws, scientific_checks, sha256

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/"GWAS/finemapping/results/alphagenome_campaign/week1-20260915"
OLD = BASE/"model/sweep_21772726"
ARMS = ADAPTERS+CONTROLS
KINDS = ("source_se", "constant_se")


def read_json(path):
    return json.loads(path.read_text())


def verify_source_se(run):
    labels = pd.read_csv(run/"labels_with_exact_source_se.tsv.gz", sep="\t")
    assert len(labels) == 32322 and labels.lead_variant_id.is_unique
    columns = ["variant_id", "target_id", "beta_nominal", "varbeta"]
    nominal = pd.read_csv(BASE/"data/endpoints-revised-21772544/caqtl_source_leads.tsv.gz", sep="\t", usecols=columns)
    assert not nominal.duplicated(["variant_id", "target_id"]).any()
    matched = labels.merge(nominal, left_on=["lead_variant_id", "peak_id"],
                           right_on=["variant_id", "target_id"], how="left", validate="one_to_one", suffixes=("", "_verified_again"))
    assert matched.varbeta_verified_again.gt(0).all()
    np.testing.assert_allclose(matched.beta_source, matched.beta_nominal_verified_again, rtol=1e-5, atol=1e-7)
    np.testing.assert_allclose(labels.verified_se, np.sqrt(matched.varbeta_verified_again), rtol=0, atol=1e-12)
    train = labels.heldout_fold.isin([2, 3, 4])
    assert train.sum() == 18914
    median = float(np.median(labels.loc[train, "verified_se"]))
    np.testing.assert_allclose(labels.constant_se, median, rtol=0, atol=1e-15)
    np.testing.assert_allclose(float(read_json(run/"recipe.json")["constant_se"]), median, rtol=0, atol=1e-15)
    return labels, median


def load_run(path, expected_steps, expected_variants):
    frame = pd.read_csv(path/"validation_predictions.tsv", sep="\t").sort_values("variant_id").reset_index(drop=True)
    assert len(frame) == 6834 and frame.variant_id.is_unique and frame.validation_fold.eq(1).all()
    np.testing.assert_array_equal(frame.variant_id, expected_variants)
    assert np.isfinite(frame[["observed_beta", "predicted_beta"]].to_numpy(float)).all()
    split = read_json(path/"split.json")
    assert split["training_n"] == 18914 and split["validation_n"] == 6834
    assert split["held_fold"] == 0 and split["validation_fold"] == 1
    assert split["held_outcomes_evaluated"] is False and split["training_native_slope"] == 0
    assert split["seed"] == 1103 and split["loss_units"] == "source_FastQTL_ALT_dosage_beta"
    feasibility = read_json(path/"feasibility.json")
    assert feasibility["completed_steps"] == feasibility["requested_steps"] == expected_steps
    assert feasibility["microbatch"] == 4 and feasibility["validation_examples"] == 6834
    for flag in ("frozen_parameters_unchanged", "running_statistics_unchanged", "save_reload_agreement"):
        assert feasibility[flag] is True
    return frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=2000)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run numerical review on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    completion = read_json(args.run/"complete.json")
    records = read_json(args.run/"runs.json")
    pd.DataFrame(records).to_csv(args.out/"producer_run_dispositions.tsv", sep="\t", index=False)
    if completion.get("completed_full_arms") != 20 or completion["status"] != "completed_fixed_development_comparison":
        receipt = dict(status="unresolved_incomplete_twenty_arm_comparison", producer=completion,
                       full_runs_required=20, statistical_comparisons_run=False,
                       reason="retain_failed_and_deferred_runs; no_restricted_success_only_family")
        (args.out/"checks.json").write_text(json.dumps(receipt, indent=2)+"\n")
        print(json.dumps(receipt, indent=2))
        return
    scientific_checks()
    labels, median = verify_source_se(args.run)
    expected = labels.loc[labels.heldout_fold.eq(1)].sort_values("lead_variant_id")
    expected_variants = expected.lead_variant_id.to_numpy()
    expected_y = expected.beta_alt.to_numpy(np.float32)
    groups = expected.block_1mb.to_numpy()
    chromosomes = expected.chr.to_numpy()
    names, predictions, nuisance, drift, hashes = [], [], [], [], {}
    run_records = {r["output_name"]: r for r in records if r.get("output_name") != "probe100"}
    assert len(run_records) == 20
    old_configs = {}
    for arm in ARMS:
        source = OLD/arm
        frame = load_run(source, 5000, expected_variants)
        np.testing.assert_array_equal(frame.observed_beta.to_numpy(np.float32), expected_y)
        names.append(arm+"__mse")
        predictions.append(frame.predicted_beta.to_numpy(float))
        old_configs[arm] = read_json(source/"config.json")
        hashes[str(source/"validation_predictions.tsv")] = sha256(source/"validation_predictions.tsv")
    for kind in KINDS:
        for arm in ARMS:
            name = arm+"__gaussian_"+kind
            source = args.run/name
            record = run_records[name]
            assert record["status"] == "complete" and record["exit_code"] == 0
            assert record["complete_validation_and_frozen_state_checked"] is True
            assert record["steps_requested"] == record["steps_completed"] == 5000
            assert record["base_arm"] == arm and record["kind"] == kind
            assert record["se_column"] == ("verified_se" if kind == "source_se" else "constant_se")
            frame = load_run(source, 5000, expected_variants)
            np.testing.assert_array_equal(frame.observed_beta.to_numpy(np.float32), expected_y)
            config = read_json(source/"config.json")
            assert config["id"] == name and config["loss"] == "gaussian_effect"
            for key, value in old_configs[arm].items():
                if key not in ("id", "loss", "steps", "status"):
                    assert config[key] == value
            with np.load(source/"weights.npz", allow_pickle=False) as weights:
                log_variance = float(weights["log_residual_variance"])
            assert np.isfinite(log_variance)
            variance = float(np.logaddexp(0, log_variance))
            assert record["initial_log_residual_variance"] == -2
            assert record["final_log_residual_variance"] == log_variance
            assert record["final_residual_variance"] == variance
            nuisance.append(dict(arm=name, initial_log_residual_variance=-2,
                final_log_residual_variance=log_variance, final_residual_variance_beta_squared=variance,
                interpretation="shared_training_nuisance_variance_not_validated_individual_prediction_uncertainty"))
            diagnostic = read_json(source/"native_function_drift.json")
            assert diagnostic["measured_native_accuracy"] == "not_available"
            for assay in ("atac", "dnase"):
                assert np.isfinite(list(diagnostic[assay].values())).all()
                drift.append(dict(arm=name, assay=assay, **diagnostic[assay], interpretation="single_region_prediction_drift_only"))
            names.append(name)
            predictions.append(frame.predicted_beta.to_numpy(float))
            hashes[str(source/"validation_predictions.tsv")] = sha256(source/"validation_predictions.tsv")
    pred = np.column_stack(predictions)
    y = expected_y.astype(float)
    point = metrics(y, pred)
    assert np.isfinite(np.column_stack(point[:2])).all()
    points = pd.DataFrame(dict(arm=names, rows=6834, RMSE_beta_units=point[1], signed_spearman=point[2]))
    producer_points = pd.read_csv(args.run/"point_performance.tsv", sep="\t").set_index("arm")
    assert len(producer_points) == 20
    for row in points.iloc[10:].itertuples():
        np.testing.assert_allclose([row.RMSE_beta_units, row.signed_spearman],
            producer_points.loc[row.arm, ["RMSE_beta_units", "signed_spearman"]].to_numpy(float), rtol=0, atol=1e-12, equal_nan=True)
    points.to_csv(args.out/"point_performance.tsv", sep="\t", index=False)
    pair_specs = []
    for family, comparator in (("source_gaussian_vs_MSE", "mse"), ("source_gaussian_vs_constant_gaussian", "gaussian_constant_se")):
        for arm in ARMS:
            pair_specs.append((names.index(arm+"__gaussian_source_se"), names.index(arm+"__"+comparator), family))
    pairs = [(a, b) for a, b, _ in pair_specs]
    improvements = contrasts(point, pairs)
    output = []
    for unit, labels_for_unit, seed in (("historical_1mb_bins", groups, 20260915), ("whole_chromosomes", chromosomes, 20260916)):
        draws, units, _ = paired_draws(y, pred, labels_for_unit, pairs, args.draws, seed)
        assert np.isfinite(draws[:, :, :2]).all()
        limits = np.full((2, len(pairs), 3), np.nan)
        limits[:, :, :2] = np.quantile(draws[:, :, :2], [.025, .975], axis=0)
        ranking_defined = np.isfinite(draws[:, :, 2]).all(axis=0) & np.isfinite(improvements[:, 2])
        for j in np.flatnonzero(ranking_defined):
            limits[:, j, 2] = np.quantile(draws[:, j, 2], [.025, .975])
        nominal = (1+np.sum(np.abs(draws[:, :, 0]-improvements[:, 0]) >= np.abs(improvements[:, 0]), axis=0))/(args.draws+1)
        q = np.concatenate((bh(nominal[:10]), bh(nominal[10:])))
        for j, (a, b, family) in enumerate(pair_specs):
            row = dict(arm=names[a], comparator=names[b], family=family, uncertainty_unit=unit,
                rows=6834, resampling_units=len(units), draws=args.draws, resampling_seed=seed,
                planned_family_contrasts=10, nominal_centered_bootstrap_MSE_p=nominal[j], exploratory_BH_q=q[j],
                ranking_interval_status="defined" if ranking_defined[j] else "undefined_for_constant_predictions_or_resampled_predictions; all_MSE_rows_and_tests_retained")
            for k, quantity in enumerate(("MSE_improvement_beta_squared", "RMSE_improvement_beta", "signed_spearman_improvement")):
                row[quantity] = improvements[j, k]
                row[quantity+"_low95"] = limits[0, j, k]
                row[quantity+"_high95"] = limits[1, j, k]
            output.append(row)
    pd.DataFrame(output).to_csv(args.out/"fixed_loss_comparisons.tsv", sep="\t", index=False)
    pd.DataFrame(nuisance).to_csv(args.out/"training_nuisance_variance.tsv", sep="\t", index=False)
    pd.DataFrame(drift).to_csv(args.out/"native_prediction_drift.tsv", sep="\t", index=False)
    receipt = dict(status="pass", complete_Gaussian_runs=20, MSE_comparators=10, rows=6834,
        native_source_SE_exact_variant_peak_verified=True, constant_training_median_SE=median,
        fixed_comparison_families=2, planned_contrasts_per_family=10, resampling_sensitivities=2,
        training_variance_reloaded_and_matches_receipts=True, all_producer_point_metrics_reconstructed=True,
        uncertainty="conditional_paired_bin_or_chromosome_bootstrap_of_fixed_predictions; no_refitting_seed_or_source_label_uncertainty",
        limits=["Only_five_validation_chromosomes; historical_bins_not_verified_LD_units.",
            "Significance_selected_marginal_associations; no_winners_curse_or_LD_correction.",
            "Source_SE_vs_MSE_changes_whole_likelihood_recipe; constant_SE_comparator_controls_Gaussian_recipe.",
            "Native_prediction_drift_is_not_measured_native_accuracy.",
            "One_seed_split; already_inspected_development; no_protected_or_external_generalization.",
            "Shared_residual_variance_is_not_validated_individual_prediction_uncertainty."],
        protected_outcomes_read=False, refitting=False, hashes=hashes)
    (args.out/"checks.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({k: v for k, v in receipt.items() if k != "hashes"}, indent=2))


if __name__ == "__main__":
    main()
