#!/usr/bin/env python3
"""Independently reconstruct saved regional ATAC retention, without model inference."""
import argparse
import hashlib
import json
import os
import platform
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915"
FOLDS = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/fold_manifest.tsv"
PEAKS = ROOT / "GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1/supplementalData1_liver_ATAC_peaks.bed.gz"
MODELS = ("original",) + tuple(
    f"adapter_r{rank}_last{blocks}_{pool}" for rank in (4, 16)
    for blocks in (3, 5) for pool in ("variant", "symmetric")
) + ("partial_qv_last3_lr3e-06", "partial_qv_last3_lr1e-05")
ALL_MODELS = MODELS + ("training_mean_target",)
TRACK_INDICES = (129, 130, 145)
TRACK_NAMES = ("UBERON:0001114 ATAC-seq", "UBERON:0001115 ATAC-seq", "UBERON:0002107 ATAC-seq")
SEED = 20260915
# Independent weighted algebra uses float64; producer centers float32 predictors.
FLOAT_TOLERANCE = 64 * np.finfo(np.float32).eps


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(2 ** 20), b""):
            digest.update(block)
    return digest.hexdigest()


def emit(path, value):
    path.write_text(json.dumps(value, indent=2, default=str, allow_nan=False) + "\n")


def complete_status(completion):
    return (completion.get("status") == "measured_regional_native_ATAC_comparison_complete"
            and completion.get("models") == 11 and completion.get("regions_per_model") == 1024
            and completion.get("all_required_models_complete") is True)


def require_tracks(values, names, peaks, panel):
    if tuple(names) != MODELS or len(set(names)) != 11:
        raise ValueError("Incomplete or substituted native-backbone family")
    np.testing.assert_array_equal(peaks, panel.peak_id)
    if values.shape != (11, 1024, 3) or not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("Missing, nonfinite, or negative native track predictions")


def weighted_fit(x, y, weights):
    """Affine fit from multiplicity weights, independent of expanded-row producer."""
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    if not np.isfinite(x).all() or not np.isfinite(y).all() or w.sum() <= 0 or np.any(w < 0):
        raise ValueError("Invalid weighted calibration inputs")
    xmean = np.einsum("i,ij->j", w, x) / w.sum()
    ymean = np.dot(w, y) / w.sum()
    centered = x - xmean
    denominator = np.einsum("i,ij,ij->j", w, centered, centered)
    numerator = np.einsum("i,ij,i->j", w, centered, y - ymean)
    slope = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 1e-20)
    intercept = ymean - slope * xmean
    return slope, intercept, ymean


def weighted_errors(x, y, training_weights, validation_weights):
    slope, intercept, mean = weighted_fit(x, y, training_weights)
    prediction = np.column_stack((np.asarray(x, float) * slope + intercept, np.full(len(y), mean)))
    error = np.einsum("i,ij,ij->j", validation_weights, prediction - y[:, None], prediction - y[:, None])
    return error / validation_weights.sum(), prediction, slope, intercept


def rho_columns(x, y):
    return np.array([spearmanr(x[:, j], y).statistic if np.ptp(x[:, j]) > 1e-12 else np.nan
                     for j in range(x.shape[1])])


def family_bh(pvalues):
    p = np.asarray(pvalues, float)
    order = np.argsort(p)
    ordered = p[order] * len(p) / (np.arange(len(p)) + 1)
    corrected = np.minimum.accumulate(ordered[::-1])[::-1]
    answer = np.empty(len(p))
    answer[order] = np.minimum(corrected, 1)
    return answer


def compare(observed, expected, label, differences, tight=False):
    observed, expected = np.asarray(observed, float), np.asarray(expected, float)
    if observed.shape != expected.shape or not np.array_equal(np.isfinite(observed), np.isfinite(expected)):
        raise ValueError(label + ": shape or missingness differs")
    finite = np.isfinite(expected)
    delta = np.abs(observed[finite] - expected[finite])
    bound = np.full(delta.shape, 1e-12) if tight else FLOAT_TOLERANCE * np.maximum(1, np.abs(expected[finite]))
    if np.any(delta > bound):
        raise ValueError(label + ": exceeds declared numerical tolerance; max=" + str(delta.max()))
    differences[label] = {"max_absolute_difference": float(delta.max()) if len(delta) else 0.,
                          "maximum_tolerance": float(bound.max()) if len(bound) else 1e-12}


def check_panel(fixture):
    panel = pd.read_csv(fixture / "panel.tsv", sep="\t", dtype={"peak_id": str})
    target_receipt = read_json(fixture / "targets_receipt.json")
    recipe = read_json(fixture / "panel_recipe.json")
    assert target_receipt["status"] == "native_profile_targets_prepared"
    assert target_receipt["donors"] == 138 and target_receipt["full_peak_library_regions"] == 349685
    assert target_receipt["source_md5"] == "6f3af5bac0509909d015e8fd99cc1d5e"
    assert target_receipt["source_VST_or_normalization_factors_used"] is False
    assert target_receipt["panel_sha256"] == sha256(fixture / "panel.tsv")
    assert recipe["count_values_accessed_by_panel_selector"] is False
    assert recipe["fold0_per_region_outcomes_exported"] is False
    assert len(panel) == 1024 and panel.peak_id.is_unique
    assert not panel.duplicated(["chrom", "start0", "end0"]).any()
    assert set(panel.role) == {"training", "validation"} and not panel.fold.eq(0).any()
    assert panel.role.eq("training").sum() == panel.role.eq("validation").sum() == 512
    assert set(panel.loc[panel.role.eq("training"), "fold"]) == {2, 3, 4}
    assert set(panel.loc[panel.role.eq("validation"), "fold"]) == {1}
    assert panel.window_end0.sub(panel.window_start0).eq(2048).all()
    assert panel.start0.ge(panel.window_start0).all() and panel.end0.le(panel.window_end0).all()
    assert panel.end0.gt(panel.start0).all()
    # Rebuild coordinate selection from public BED, exact source count axis and saved sequence exclusions.
    assert recipe["fold_manifest_sha256"] == sha256(FOLDS)
    assert recipe["peaks_sha256"] == sha256(PEAKS)
    axis_path = fixture / "source_peak_axis.tsv"
    assert recipe["source_axis_sha256"] == sha256(axis_path)
    axis = pd.read_csv(axis_path, sep="\t", dtype=str).peak_id
    assert len(axis) == 349685 and axis.is_unique
    fold_map = {}
    for row in pd.read_csv(FOLDS, sep="\t").itertuples(index=False):
        for chrom in row.test_chromosomes.split(","):
            assert chrom not in fold_map
            fold_map[chrom] = int(row.fold)
    assert set(fold_map) == {f"chr{i}" for i in range(1, 23)}
    np.testing.assert_array_equal(panel.chrom.map(fold_map), panel.fold)
    peaks = pd.read_csv(PEAKS, sep="\t", dtype={"peakID": str})
    assert peaks.peakID.is_unique
    source = peaks.loc[peaks.peakID.isin(axis), ["#chr", "start", "end", "peakID"]].copy()
    assert len(source) == 349685
    source.columns = ["chrom", "start0", "end0", "peak_id"]
    source["fold"] = source.chrom.map(fold_map)
    assert source.fold.notna().all()
    source = source.loc[source.fold.isin([1, 2, 3, 4]) & source.end0.sub(source.start0).between(1, 2048)].copy()
    source["coordinate_hash"] = [hashlib.sha256(f"native-retention|{SEED}|{c}:{s}-{e}|{p}".encode()).hexdigest()
                                 for c, s, e, p in source[["chrom", "start0", "end0", "peak_id"]].itertuples(index=False, name=None)]
    rejected = pd.read_csv(fixture / "sequence_rejections.tsv", sep="\t", dtype={"peak_id": str})
    assert rejected.peak_id.is_unique
    assert set(rejected.reason) <= {"incomplete_reference_window", "non_ACGT_reference_window"}
    assert set(rejected.peak_id) <= set(source.peak_id) and not set(rejected.peak_id) & set(panel.peak_id)
    for role, allowed in (("training", [2, 3, 4]), ("validation", [1])):
        candidates = source.loc[source.fold.isin(allowed)].sort_values(["coordinate_hash", "peak_id"])
        chosen = candidates.loc[~candidates.peak_id.isin(rejected.peak_id)].head(512)
        observed = panel.loc[panel.role.eq(role)]
        for column in ("chrom", "start0", "end0", "peak_id", "fold", "coordinate_hash"):
            np.testing.assert_array_equal(chosen[column], observed[column])
    # Independent interval-component partition; numeric component labels need not coincide.
    memberships = np.full(len(panel), -1, int)
    component = -1
    for chrom in sorted(panel.chrom.unique()):
        right = -1
        for index in panel.index[panel.chrom.eq(chrom)].to_numpy()[np.argsort(panel.loc[panel.chrom.eq(chrom), "window_start0"].to_numpy(), kind="stable")]:
            row = panel.loc[index]
            if row.window_start0 >= right:
                component += 1
            memberships[index] = component
            right = max(right, int(row.window_end0))
    partition = pd.DataFrame({"independent": memberships, "producer": panel.input_component, "role": panel.role})
    assert partition.groupby("independent").producer.nunique().eq(1).all()
    assert partition.groupby("producer").independent.nunique().eq(1).all()
    assert partition.groupby("independent").role.nunique().eq(1).all()
    return panel, target_receipt


def count_targets(fixture, panel):
    root = fixture / "private_targets"
    raw = pd.read_csv(root / "raw_counts.tsv.gz", sep="\t", dtype={"peak_id": str})
    stored = pd.read_csv(root / "log2cpm.tsv.gz", sep="\t", dtype={"peak_id": str})
    totals = pd.read_csv(root / "donor_library_totals.tsv", sep="\t", dtype={"donor_id": str})
    assert len(totals) == 138 and totals.donor_id.is_unique
    for frame in (raw, stored):
        np.testing.assert_array_equal(frame.peak_id, panel.peak_id)
        assert list(frame.columns[1:]) == list(totals.donor_id)
    counts = raw.iloc[:, 1:].to_numpy(float)
    library = totals.full_peak_library_count_total.to_numpy(float)
    assert counts.shape == (1024, 138) and np.isfinite(counts).all() and np.isfinite(library).all()
    assert (counts >= 0).all() and (counts == np.floor(counts)).all() and (library > 0).all()
    assert (counts.sum(axis=0) <= library).all()
    target = np.log1p(counts / library * 1e6) / np.log(2)
    np.testing.assert_allclose(target, stored.iloc[:, 1:].to_numpy(float), rtol=0, atol=1e-12)
    return target


def verify_models(run, fixture, panel):
    design = read_json(run / "design.json")
    assert tuple(design["models"]) == MODELS and design["sequence_length"] == 2048
    assert design["orientation"] == "reference_forward_only"
    assert design["all_ten_MSE_adaptations_retained"] is True
    assert design["source_GC_correction_VST_offsets_used"] is False
    assert design["no_DNase_or_whole_native_function_claim"] is True
    assert design["panel_sha256"] == sha256(fixture / "panel.tsv")
    assert design["targets_receipt_sha256"] == sha256(fixture / "targets_receipt.json")
    completed = read_json(run / "completed_models.json")
    assert tuple(row["model"] for row in completed) == MODELS
    assert all(row["regions"] == 1024 and row["repeat_exact"] is True for row in completed)
    assert tuple(row["name"] for row in design["configuration_records"]) == MODELS[1:]
    for row in design["configuration_records"]:
        weights = Path(row["weights"])
        assert row["actual_steps"] == 5000 and sha256(weights) == row["weights_sha256"]
        saved = read_json(Path(str(weights) + ".json"))
        config = saved["config"]
        assert config["id"] == row["name"] and config["loss"] == "mse" and config["length"] == 2048
        assert saved["checkpoint"] == design["checkpoint"] and saved["paths"] == row["paths"]
        assert config["mode"] == row["mode"]
        feasibility = read_json(weights.parent / "feasibility.json")
        assert feasibility["completed_steps"] == 5000
        assert feasibility["running_statistics_unchanged"] is True and feasibility["frozen_parameters_unchanged"] is True
    tracks = pd.read_csv(run / "native_liver_ATAC_tracks.tsv", sep="\t")
    assert tuple(tracks.native_channel_index) == TRACK_INDICES and tuple(tracks.name) == TRACK_NAMES
    for column, value in (("biosample_life_stage", "adult"), ("biosample_type", "tissue"),
                          ("data_source", "encode"), ("strand", "."), ("Assay title", "ATAC-seq")):
        assert tracks[column].eq(value).all()
    with np.load(run / "native_track_predictions.npz", allow_pickle=False) as stored:
        values = stored["log2_track_sums"]
        require_tracks(values, stored["model_names"], stored["peak_id"], panel)
        assert tuple(stored["native_track_indices"]) == TRACK_INDICES
    return values.mean(axis=2).T, int(design["bootstrap_draws"])


def group_plan(panel, mask, unit, stratify):
    plans = []
    for fold in sorted(panel.loc[mask, "fold"].unique()) if stratify else [None]:
        use = mask if fold is None else mask & panel.fold.eq(fold).to_numpy()
        indices = np.flatnonzero(use)
        unique, inverse = np.unique(panel.loc[use, unit].to_numpy(), return_inverse=True)
        plans.append((indices, unique, inverse))
    return plans


def draw_weights(plan, rng, n):
    weights = np.zeros(n, float)
    for indices, unique, inverse in plan:
        sampled = rng.choice(unique, len(unique), replace=True)
        multiplicities = np.bincount(np.searchsorted(unique, sampled), minlength=len(unique))
        weights[indices] = multiplicities[inverse]
    return weights


def weighted_bootstrap(x, target, panel, unit, draws):
    rng = np.random.default_rng(SEED)
    training = group_plan(panel, panel.role.eq("training").to_numpy(), unit, True)
    validation = group_plan(panel, panel.role.eq("validation").to_numpy(), unit, False)
    result = np.empty((draws, 11))
    for draw in range(draws):
        donor_ids = rng.choice(target.shape[1], target.shape[1], replace=True)
        donor_weights = np.bincount(donor_ids, minlength=target.shape[1]) / target.shape[1]
        average = target @ donor_weights
        train_weights = draw_weights(training, rng, len(panel))
        validation_weights = draw_weights(validation, rng, len(panel))
        errors, _, _, _ = weighted_errors(x, average, train_weights, validation_weights)
        result[draw] = errors[0] - errors[1:]
    assert np.isfinite(result).all()
    return result


def invariants():
    rng = np.random.default_rng(1103)
    x = rng.normal(size=(14, 3)); y = rng.normal(size=14)
    tr = np.array([2, 0, 3, 1, 0, 2, 1] + [0] * 7, float)
    va = np.array([0] * 7 + [1, 2, 0, 1, 3, 1, 1], float)
    errors, prediction, slopes, intercepts = weighted_errors(x, y, tr, va)
    expanded = np.repeat(np.arange(len(y)), tr.astype(int))
    held = np.repeat(np.arange(len(y)), va.astype(int))
    for j in range(3):
        coefficients = np.linalg.lstsq(np.column_stack((x[expanded, j], np.ones(len(expanded)))), y[expanded], rcond=None)[0]
        np.testing.assert_allclose([slopes[j], intercepts[j]], coefficients, rtol=0, atol=1e-12)
    np.testing.assert_allclose(errors, ((prediction[held] - y[held, None]) ** 2).mean(0), rtol=0, atol=1e-12)
    changed = y.copy(); changed[7:] += 1000000
    for first, second in zip(weighted_fit(x, y, tr), weighted_fit(x, changed, tr)):
        np.testing.assert_allclose(first, second, rtol=0, atol=1e-12)
    constant = np.ones((len(y), 1))
    cerror, cpred, cslope, _ = weighted_errors(constant, y, tr, va)
    assert cslope[0] == 0 and np.isfinite(cerror).all() and np.ptp(cpred) == 0
    assert not complete_status({"status": "measured_regional_native_ATAC_comparison_complete", "models": 10})
    fake = pd.DataFrame({"peak_id": np.arange(1024).astype(str)})
    for values, names in ((np.full((11, 1024, 3), np.nan), MODELS), (np.ones((11, 1024, 3)), MODELS[:-1])):
        try:
            require_tracks(values, names, fake.peak_id.to_numpy(), fake)
        except ValueError:
            pass
        else:
            raise AssertionError("Missing/nonfinite/incomplete-family invariant failed")
    np.testing.assert_array_equal(family_bh(np.ones(11)), np.ones(11))
    return ["weighted_affine_matches_independent_lstsq", "weighted_errors_match_duplicate_rows", "held_targets_cannot_change_calibration",
            "constant_predictor_preserves_MSE", "missing_model_and_nonfinite_outputs_rejected", "incomplete_completion_rejected", "complete_family_retained"]


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run numerical review on an authorized compute node")
    args.out.mkdir(parents=True, exist_ok=False, mode=0o700)
    shutil.copy2(__file__, args.out / "executed_review_native_retention.py")
    started = time.monotonic()
    checks = invariants()
    completion_path = args.run / "completion.json"
    completion = read_json(completion_path) if completion_path.is_file() else {}
    if not complete_status(completion):
        emit(args.out / "checks.json", {"status": "not_evaluated_incomplete_native_comparison", "producer_completion": completion,
             "all_eleven_backbones_required": True, "statistical_comparisons_run": False, "synthetic_invariants": checks})
        return
    differences = {}
    panel, target_receipt = check_panel(args.fixture)
    target = count_targets(args.fixture, panel)
    x, draws = verify_models(args.run, args.fixture, panel)
    assert draws == 1000
    training = panel.role.eq("training").to_numpy()
    validation = panel.role.eq("validation").to_numpy()
    average = target.mean(axis=1)
    errors, predicted, slopes, intercepts = weighted_errors(x, average, training.astype(float), validation.astype(float))
    raw = np.column_stack((x, np.zeros(len(x))))
    native_rho = rho_columns(raw[validation], average[validation])
    calibrated_rho = rho_columns(predicted[validation], average[validation])
    points = pd.DataFrame({"model": ALL_MODELS, "RMSE_donor_mean_log2_peak_CPM": np.sqrt(errors),
        "raw_predicted_scale_spearman": native_rho, "calibrated_spearman": calibrated_rho,
        "training_slope": np.r_[slopes, 0], "training_intercept": np.r_[intercepts, predicted[0, -1]],
        "input_order_reversed_by_calibration": np.r_[slopes < 0, False]})
    deposited = pd.read_csv(args.run / "measured_native_performance.tsv", sep="\t")
    assert tuple(deposited.model) == ALL_MODELS and deposited.regions.eq(512).all()
    assert deposited.source_described_donors.eq(138).all() and deposited.native_context_bp.eq(2048).all()
    for column in points.columns[1:-1]:
        compare(points[column], deposited[column], column, differences, tight="spearman" in column)
    np.testing.assert_array_equal(points.input_order_reversed_by_calibration, deposited.input_order_reversed_by_calibration)
    regional = pd.read_csv(args.run / "regional_predictions.tsv.gz", sep="\t", dtype={"peak_id": str})
    np.testing.assert_array_equal(regional.peak_id, panel.peak_id)
    compare(average, regional.observed_donor_mean_log2_peak_CPM, "regional_measured_donor_means", differences, tight=True)
    for i, name in enumerate(ALL_MODELS):
        compare(predicted[:, i], regional[name + "__calibrated"], name + "_regional_calibration", differences)
        compare(raw[:, i], regional[name + "__raw_log_track_sum"], name + "_raw_log_track_sum", differences)
    changed = average.copy(); changed[validation] += 1e6
    for first, second in zip(weighted_fit(x, average, training), weighted_fit(x, changed, training)):
        np.testing.assert_allclose(first, second, rtol=0, atol=1e-12)
    producer = pd.read_csv(args.run / "measured_native_paired_uncertainty.tsv", sep="\t")
    assert len(producer) == 66 and not producer.duplicated(["model", "endpoint", "resampling_unit"]).any()
    assert set(producer.endpoint) == {"MSE_improvement", "calibrated_spearman_improvement", "raw_spearman_improvement"}
    assert set(producer.resampling_unit) == {"input_component", "chrom"}
    assert producer.planned_family_n.eq(11).all() and producer.draws.eq(draws).all()
    assert producer.paired_donor_resampling.all() and producer.calibration_refit_each_draw.all()
    assert producer.comparator.eq("original").all()
    assert np.isfinite(producer[["nominal_p", "BH_q_complete_exploratory_family"]].to_numpy()).all()
    assert producer.nominal_p.between(0, 1).all() and producer.BH_q_complete_exploratory_family.between(0, 1).all()
    for family_key, family in producer.groupby(["endpoint", "resampling_unit"]):
        assert set(family.model) == set(ALL_MODELS[1:])
        compare(family_bh(family.nominal_p), family.BH_q_complete_exploratory_family, "family_accounting_" + str(family_key), differences, tight=True)
    point = errors[0] - errors[1:]
    uncertainty_rows = []
    interpretation_flags = []
    for unit in ("input_component", "chrom"):
        sampled = weighted_bootstrap(x, target, panel, unit, draws)
        low, high = np.quantile(sampled, [0.025, 0.975], axis=0)
        pvalues = (1 + np.sum(np.abs(sampled - point[None]) >= np.abs(point[None]), axis=0)) / (draws + 1)
        qvalues = family_bh(pvalues)
        observed = producer.loc[producer.resampling_unit.eq(unit) & producer.endpoint.eq("MSE_improvement")].set_index("model").loc[list(ALL_MODELS[1:])]
        assert observed.estimable_draws.eq(draws).all()
        for column, values in (("point", point), ("low95", low), ("high95", high), ("nominal_p", pvalues), ("BH_q_complete_exploratory_family", qvalues)):
            compare(values, observed[column], unit + "_MSE_" + column, differences, tight=column in ("nominal_p", "BH_q_complete_exploratory_family"))
        for i, name in enumerate(ALL_MODELS[1:]):
            numeric_bound = FLOAT_TOLERANCE * max(1, abs(point[i]), abs(low[i]), abs(high[i]))
            decisive_distance = min(abs(point[i]), abs(low[i]), abs(high[i]))
            small = decisive_distance <= 10 * numeric_bound
            if small:
                interpretation_flags.append({"model": name, "resampling_unit": unit,
                    "reason": "point_or_CI_endpoint_within_ten_declared_float32_bounds_of_zero; no_stable_direction_claim_from_this_check"})
            uncertainty_rows.append({"model": name, "comparator": "original", "resampling_unit": unit,
                "endpoint": "MSE_improvement", "point": point[i], "low95": low[i], "high95": high[i],
                "nominal_p": pvalues[i], "BH_q_complete_exploratory_family": qvalues[i], "planned_family_n": 11,
                "draws": draws, "independent_weighted_algebra": True, "numerically_small_relative_to_tolerance": small})
    points.to_csv(args.out / "independent_point_performance.tsv", sep="\t", index=False)
    pd.DataFrame(uncertainty_rows).to_csv(args.out / "independent_MSE_uncertainty.tsv", sep="\t", index=False)
    emit(args.out / "checks.json", {"status": "independent_native_retention_reconstruction_passed", "job_id": os.environ["SLURM_JOB_ID"],
        "synthetic_invariants": checks, "complete_models": 11, "point_rows_including_mean": 12, "MSE_contrasts": 22,
        "calibration_regions": 512, "validation_regions": 512, "source_described_donors": 138,
        "validation_components": int(panel.loc[validation, "input_component"].nunique()),
        "validation_chromosomes": int(panel.loc[validation, "chrom"].nunique()), "draws_per_scheme": draws,
        "statistical_comparisons_run": True, "differences": differences, "interpretation_flags": interpretation_flags,
        "tolerance": "64*float32_epsilon*max(1,abs(producer_value)) for independent float64 affine/error algebra; 1e-12 ranking points and discrete p/BH",
        "bootstrap_reconstruction": "donor multiplicities and stratified region-component multiplicities; independent weighted affine fit; shared donor draw in both partitions",
        "ranking_CI_independently_reconstructed": False, "ranking_CI_scope": "family membership and BH accounting checked; ranking bootstrap draws and intervals not independently rederived",
        "inference_scope": "same_Currin_source_fixed_adapted_models_regional_ATAC_2kb; neither independent_cohort_nor_global_native_function_accuracy",
        "uncertainty_scope": "actual_window_components_not_claimed_LD_independent; only_five_validation_chromosomes; adaptation_training_and_selection_not_resampled",
        "source_scope": "private_extracted_raw_counts_and_full_deposited_peak_library_totals; full_RData_extraction_and_GPU_reference_sequence_calls_not_reexecuted",
        "backbone_scope": "all_saved_backbone_identifiers_weight_hashes_training_receipts_verified; no_new_model_forward_call_or_parameter_tree_hash_reconstruction",
        "participant_level_outputs_written": False, "fold0_regional_outcomes_read": False, "protected_outcomes_read": False,
        "source_files_sha256": {str(path): sha256(path) for path in (args.fixture / "panel.tsv", args.fixture / "targets_receipt.json",
            args.fixture / "private_targets/raw_counts.tsv.gz", args.fixture / "private_targets/donor_library_totals.tsv",
            args.run / "native_track_predictions.npz", args.run / "design.json")},
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "seconds": time.monotonic() - started})
    assert sum(path.stat().st_size for path in args.out.iterdir() if path.is_file()) < 100000000


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        main(args)
    except Exception as exc:
        if args.out.is_dir():
            emit(args.out / "checks.json", {"status": "independent_native_retention_check_failed", "error": repr(exc),
                 "successful_comparison_claim": False, "failure_requires_resolution_not_success_only_subset" : True})
        raise
