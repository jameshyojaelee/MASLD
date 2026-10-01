#!/usr/bin/env python3
"""Post-pilot paired uncertainty from saved adaptation validation predictions.

No fitting, model selection, new sequence inference, or fold-0 outcomes.
The ten contrasts and two resampling units were fixed after pilot inspection.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata


ADAPTERS = [f"adapter_r{rank}_last{blocks}_{pool}"
            for rank in (4, 16) for blocks in (3, 5)
            for pool in ("variant", "symmetric")]
PARTIALS = ["partial_qv_last3_lr1e-05", "partial_qv_last3_lr3e-06"]
CONTROLS = [f"frozen_2048_{pool}_shared_mlp64" for pool in ("variant", "symmetric")]
ARMS = ADAPTERS + PARTIALS + CONTROLS


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bh(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    adjusted = np.minimum.accumulate((values[order] * len(values) /
                                     np.arange(1, len(values)+1))[::-1])[::-1]
    result = np.empty_like(values)
    result[order] = np.minimum(1, adjusted)
    return result


def metrics(y, predictions):
    mse = np.mean(np.square(predictions-y[:, None]), axis=0)
    ranks = rankdata(np.column_stack((y, predictions)), axis=0)
    ranks -= ranks.mean(axis=0)
    denominator = np.sqrt(np.sum(ranks[:, :1]**2)*np.sum(ranks[:, 1:]**2, axis=0))
    rho = np.divide(np.sum(ranks[:, :1]*ranks[:, 1:], axis=0), denominator,
                    out=np.full(predictions.shape[1], np.nan), where=denominator > 0)
    return mse, np.sqrt(mse), rho


def contrasts(values, pairs):
    mse, rmse, rho = values
    return np.array([[mse[control]-mse[arm], rmse[control]-rmse[arm],
                      rho[arm]-rho[control]] for arm, control in pairs])


def paired_draws(y, predictions, units, pairs, draws, seed):
    unique, inverse = np.unique(np.asarray(units, dtype=str), return_inverse=True)
    if len(unique) < 2:
        raise ValueError("At least two resampling units are required")
    members = [np.flatnonzero(inverse == j) for j in range(len(unique))]
    rng = np.random.default_rng(seed)
    sampled = np.empty((draws, len(pairs), 3), dtype=float)
    counts = np.empty(draws, dtype=int)
    for iteration in range(draws):
        indices = np.concatenate([members[j] for j in rng.integers(0, len(unique), len(unique))])
        counts[iteration] = len(indices)
        # Ranks are recomputed after resampling, including repeated whole blocks.
        sampled[iteration] = contrasts(metrics(y[indices], predictions[indices]), pairs)
    return sampled, unique, counts


def scientific_checks():
    y = np.arange(12, dtype=float)/3
    pred = np.column_stack((y, y+1, y+1))
    values, units, _ = paired_draws(y, pred, np.repeat(["a", "b", "c"], 4),
                                   [(0, 1), (2, 1)], 40, 710)
    np.testing.assert_allclose(values[:, 0, :2], 1, rtol=0, atol=1e-14)
    np.testing.assert_allclose(values[:, 1], 0, rtol=0, atol=0)
    np.testing.assert_allclose(bh([.01, .04, 1]), [.03, .06, 1], rtol=0, atol=1e-14)
    assert len(units) == 3
    return dict(identical_prediction_zero_difference=True,
                constant_error_benefit_interval=True, complete_family_BH=True)


def load_inputs(sweep, manifest):
    metadata = pd.read_csv(manifest, sep="\t", dtype={"key": str, "chr": str})
    assert metadata.key.is_unique and not metadata.key.isna().any()
    metadata = metadata.set_index("key")
    runs = {row["name"]: row for row in read_json(sweep/"runs.json")}
    predictions, exposures, drift, hashes = [], [], [], {}
    authority = None
    for arm in ARMS:
        source = sweep/arm
        frame = pd.read_csv(source/"validation_predictions.tsv", sep="\t")
        assert frame.variant_id.is_unique and len(frame) == 6834
        assert frame.validation_fold.eq(1).all()
        frame["key"] = frame.variant_id.str.replace(r"^chr", "", regex=True)
        frame = frame.sort_values("key").reset_index(drop=True)
        if authority is None:
            authority = frame[["key", "variant_id", "observed_beta", "validation_fold"]].copy()
        else:
            pd.testing.assert_frame_equal(frame[authority.columns], authority, check_exact=True)
        predictions.append(frame.predicted_beta.to_numpy(dtype=float))
        config = read_json(source/"config.json")
        split = read_json(source/"split.json")
        run = runs[arm]
        feasibility = read_json(source/"feasibility.json")
        assert config["id"] == arm and config["length"] == 2048 and config["seed"] == 1103
        assert config["hidden"] == 64 and config["loss"] == "mse"
        assert split["held_fold"] == 0 and split["validation_fold"] == 1
        assert split["training_n"] == 18914 and split["validation_n"] == 6834
        assert split["held_outcomes_evaluated"] is False and split["training_native_slope"] == 0
        assert split["loss_units"] == "source_FastQTL_ALT_dosage_beta"
        assert run["status"] == "completed_fixed_exposure" and run["exit_code"] == 0
        assert run["steps_completed"] == run["steps_requested"] == 5000
        assert feasibility["completed_steps"] == feasibility["requested_steps"] == 5000
        assert feasibility["microbatch"] == 4 and feasibility["validation_examples"] == 6834
        for check in ("frozen_parameters_unchanged", "running_statistics_unchanged", "save_reload_agreement"):
            assert feasibility[check] is True
        if arm in ADAPTERS:
            assert config["mode"] == "adapter"
            assert arm == f"adapter_r{config['rank']}_last{config['last_blocks']}_{config['pooling']}"
        elif arm in PARTIALS:
            assert config["mode"] == "partial" and config["pooling"] == "variant"
            assert config["last_blocks"] == 3
            assert config["trainable_scope"] == "full_query_value_matrices_final_three_attention_blocks"
            assert config["backbone_lr"] in (3e-6, 1e-5)
        else:
            assert config["mode"] == "frozen" and config["head"] == "shared_mlp64"
            assert arm == f"frozen_2048_{config['pooling']}_shared_mlp64"
        exposures.append(dict(arm=arm, mode=config["mode"], pooling=config["pooling"],
            training_rows=18914, validation_rows=6834, training_folds="2,3,4",
            validation_fold=1, closed_fold=0, seed=1103, completed_steps=5000,
            microbatch=4, length_bp=2048, config_default_steps=config["steps"],
            config_steps_note="actual_CLI_exposure_from_runs_and_feasibility_overrides_config_default"))
        diagnostic = read_json(source/"native_function_drift.json")
        assert diagnostic["measured_native_accuracy"] == "not_available"
        assert diagnostic["scope"] == "separate_validation_region_predicted_native_ATAC_DNase_track_drift"
        for assay in ("atac", "dnase"):
            values = diagnostic[assay]
            assert np.isfinite([values["rms_log2_track_change"], values["max_abs_log2_track_change"]]).all()
            assert 0 <= values["rms_log2_track_change"] <= values["max_abs_log2_track_change"]
            drift.append(dict(arm=arm, assay=assay, **values, measured_native_accuracy="not_available",
                              scope="one_validation_region_prediction_drift_only"))
        for filename in ("validation_predictions.tsv", "config.json", "split.json", "feasibility.json", "native_function_drift.json"):
            hashes[str(source/filename)] = sha256(source/filename)
    identities = metadata.loc[authority.key].reset_index()
    assert identities.heldout_fold.eq(1).all()
    np.testing.assert_array_equal(identities.key, authority.key)
    assert identities.chr.str.fullmatch(r"chr[0-9]+").all()
    for row in identities.itertuples():
        assert row.block_1mb == f"{row.chr}:{row.pos_hg38//1_000_000}"
    predictions = np.column_stack(predictions)
    assert np.isfinite(predictions).all() and np.isfinite(authority.observed_beta).all()
    hashes[str(manifest)] = sha256(manifest)
    return authority, identities, predictions, pd.DataFrame(exposures), pd.DataFrame(drift), hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=2000)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run numerical analyses on a compute node")
    if args.draws < 1000:
        raise ValueError("At least 1,000 draws required for this report")
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    checks = scientific_checks()
    authority, identities, predictions, exposures, drift, hashes = load_inputs(args.sweep, args.manifest)
    y = authority.observed_beta.to_numpy(dtype=float)
    pairs = [(ARMS.index(arm), ARMS.index(f"frozen_2048_{exposures.set_index('arm').loc[arm, 'pooling']}_shared_mlp64"))
             for arm in ADAPTERS+PARTIALS]
    assert len(pairs) == 10
    point = metrics(y, predictions)
    improvement = contrasts(point, pairs)
    points = pd.DataFrame(dict(arm=ARMS, rows=len(y), source_units="FastQTL_ALT_dosage_beta",
        RMSE_beta_units=point[1], signed_spearman=point[2], coverage=1.0,
        historical_bins=identities.block_1mb.nunique(), chromosomes=identities.chr.nunique()))
    points.to_csv(args.out/"point_performance.tsv", sep="\t", index=False)
    output, draw_summary = [], []
    for unit, column, seed in (("historical_1mb_bins", "block_1mb", 20260915),
                                ("whole_chromosomes", "chr", 20260916)):
        sampled, units, rowcounts = paired_draws(y, predictions, identities[column], pairs, args.draws, seed)
        if not np.isfinite(sampled).all():
            raise ValueError("Nonfinite bootstrap contrasts; cannot silently omit an estimand")
        interval = np.quantile(sampled, [.025, .975], axis=0)
        # Two-sided centered bootstrap on the native-unit mean squared error contrast.
        nominal = (1+np.sum(np.abs(sampled[:, :, 0]-improvement[:, 0]) >=
                            np.abs(improvement[:, 0]), axis=0))/(args.draws+1)
        adjusted = bh(nominal)
        for j, (arm, control) in enumerate(pairs):
            row = dict(arm=ARMS[arm], comparator=ARMS[control], uncertainty_unit=unit,
                rows=len(y), units=len(units), validation_fold=1, seed=1103,
                resampling_seed=seed, bootstrap_draws=args.draws,
                native_RMSE=point[1][arm], comparator_native_RMSE=point[1][control],
                signed_spearman=point[2][arm], comparator_signed_spearman=point[2][control],
                nominal_centered_bootstrap_MSE_p=nominal[j], exploratory_BH_q=adjusted[j],
                planned_family_contrasts=10, estimable_family_contrasts=10,
                hypothesis_family=f"post_pilot_ten_adaptation_vs_pool_matched_control_{unit}")
            for k, name in enumerate(("MSE_improvement_beta_squared", "RMSE_improvement_beta", "signed_spearman_improvement")):
                row[name] = improvement[j, k]
                row[name+"_low95"] = interval[0, j, k]
                row[name+"_high95"] = interval[1, j, k]
            output.append(row)
        draw_summary.append(dict(unit=unit, unique_units=units.tolist(), draws=args.draws,
            minimum_resampled_rows=int(rowcounts.min()), maximum_resampled_rows=int(rowcounts.max()), seed=seed))
    pd.DataFrame(output).to_csv(args.out/"ten_paired_contrasts.tsv", sep="\t", index=False)
    identities.to_csv(args.out/"matched_molecular_identities.tsv.gz", sep="\t", index=False)
    exposures.to_csv(args.out/"training_exposure.tsv", sep="\t", index=False)
    drift.to_csv(args.out/"diagnostic_native_function_drift.tsv", sep="\t", index=False)
    report = dict(status="pass", scientific_checks=checks, analysis_timing="post_pilot_precision_analysis_after_aggregate_development_inspection_not_preregistration",
        rows=6834, arms=12, contrasts=10, families="same_fixed_ten_contrasts_separately_adjusted_with_BH_per_resampling_unit",
        population="existing_significance_selected_Currin_leads_validation_fold_1_with_identical_complete_predictions",
        source_units="FastQTL_ALT_dosage_beta; MSE_beta_squared; RMSE_beta; Spearman_unitless",
        uncertainty="paired_percentile_cluster_bootstrap_conditional_on_observed_labels_fixed_fitted_models_and_single_split_seed",
        resampling=draw_summary, inference_limits=[
            "Historical_1Mb_bins_are_not_verified_LD_independent_and_overlapping_sequence_windows_cross_bins.",
            "Whole_chromosome_resampling_is_a_sensitivity_with_few_validation_chromosomes_not_participant_replication.",
            "No_refitting_seed_variation_training_selection_uncertainty_or_source_beta_estimation_error_is_propagated.",
            "Five_seed_nested_selection_and_external_generalization_are_not_established.",
            "Validation_fold_1_has_already_been_inspected_during_pilot_development.",
            "Native_function_drift_is_prediction_change_at_one_region_without_measured_native_accuracy.",
            "No_native_1Mb_predictor_or_partial_native_coverage_is_compared.",
            "Pretraining_exposure_is_not_resolved_by_chromosome_holdout.",
            "Variant_or_chromosome_counts_are_not_effective_participant_n; per_variant_effective_n_is_unknown."],
        protected_outcomes_read=False, fitting_or_selection=False,
        hashes=hashes, environment=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
        elapsed_seconds=time.monotonic()-started)
    (args.out/"checks_and_interpretation.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k: v for k, v in report.items() if k != "hashes"}, indent=2))


if __name__ == "__main__":
    main()
