#!/usr/bin/env python3
"""Independent source/error/bootstrap reconstruction of the stronger donor study."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import review_donor_baselines as B

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915"
SOURCE = ROOT / "Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture"
FIXTURE = BASE / "model/donor_fixture_21773589"
FEATURES = BASE / "model/donor_features_21773588/features.npz"
SEED = 20260915
ARMS = ("all_rna", "annotated_all_rna", "chromosome", "random_chromosome")
BASELINES = ("training_baseline", "sequence_only", "rna_global", "additive_linear", "product_kernel")
NEURAL = tuple(architecture + "_" + objective for architecture in ("film", "additive_nn") for objective in ("profile", "residual", "pairwise"))
REGION_RNA = ("full_rna_ridge", "pcr", "pls_svd_ridge", "selected_rna_linear")


def read_json(path):
    return json.loads(path.read_text())


def stable(value):
    return int(hashlib.sha256(value.encode()).hexdigest()[:15], 16)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(2 ** 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def emit(path, value):
    path.write_text(json.dumps(value, indent=2, default=str, allow_nan=False) + "\n")


def close(first, second, label, differences, tolerance=1e-9):
    first, second = np.asarray(first, float), np.asarray(second, float)
    np.testing.assert_allclose(first, second, rtol=0, atol=tolerance, equal_nan=True, err_msg=label)
    finite = np.isfinite(first) & np.isfinite(second)
    delta = float(np.max(np.abs(first[finite] - second[finite]))) if finite.any() else 0.
    differences[label] = max(delta, differences.get(label, 0.))


def load_source():
    reg = pd.read_csv(FIXTURE / "regions.tsv", sep="\t")
    with np.load(FIXTURE / "paired.npz", allow_pickle=False) as stored:
        ids = stored["participant_ids"]; fold = stored["donor_fold"]
        x = stored["rna_log2cpm_full_library"]; y = stored["target_log2cpm_full_library"]
        genes = stored["gene_ids"]; chrom = stored["gene_chrom"]
    participant_axis = pd.read_csv(SOURCE / "molecular/participant_axis.tsv", sep="\t")
    np.testing.assert_array_equal(ids, participant_axis.participant_id.astype(str))
    assert len(ids) == len(set(ids)) == 99 and x.shape == (99, 42163) and y.shape == (99, 128)
    authority = pd.read_csv(SOURCE / "folds/participant_outer_folds.tsv", sep="\t", dtype={"participant_id": str})
    np.testing.assert_array_equal(fold, authority.set_index("participant_id").loc[ids, "outer_fold"].to_numpy(int))
    assert sorted(np.bincount(fold)) == [17, 19, 21, 21, 21]
    rna = np.load(SOURCE / "molecular/rna_values.npy", mmap_mode="r")
    h3 = np.load(SOURCE / "molecular/h3k27ac_counts.npy", mmap_mode="r")
    assert rna.shape == (99, 42163) and h3.shape == (99, 96460)
    assert np.load(SOURCE / "molecular/rna_observed_mask.npy").all()
    assert np.load(SOURCE / "molecular/h3k27ac_observed_mask.npy").all()
    assert np.issubdtype(h3.dtype, np.integer) and np.isfinite(rna).all() and np.isfinite(h3).all()
    assert (rna >= 0).all() and (h3 >= 0).all()
    rna_library = np.asarray(rna, float).sum(axis=1)
    h3_library = np.asarray(h3, float).sum(axis=1)
    assert (rna_library > 0).all() and (h3_library > 0).all()
    np.testing.assert_allclose(x, np.log1p(np.asarray(rna, float) / rna_library[:, None] * 1e6) / np.log(2), rtol=0, atol=1e-12)
    np.testing.assert_allclose(y, np.log1p(np.asarray(h3[:, reg.region_index.to_numpy()], float) / h3_library[:, None] * 1e6) / np.log(2), rtol=0, atol=1e-12)
    gene_axis = pd.read_csv(SOURCE / "molecular/rna_feature_axis.tsv", sep="\t")
    np.testing.assert_array_equal(genes, gene_axis.stable_gene_id.astype(str))
    h3_axis = pd.read_csv(SOURCE / "molecular/h3k27ac_feature_axis.tsv", sep="\t")
    np.testing.assert_array_equal(reg.region_key, h3_axis.iloc[reg.region_index].opaque_source_feature_key)
    assert (chrom != "unknown").sum() == 42126 and (chrom == "unknown").sum() == 37
    old = pd.read_csv(BASE / "review/b2-21771990/regions.tsv", sep="\t")
    for column in old.columns:
        np.testing.assert_array_equal(reg[column], old[column])
    assert reg.region_key.is_unique and reg.region_index.is_unique
    assert reg.groupby("tile_index").size().eq(2).all() and reg.tile_index.nunique() == 64
    assert reg.region_role.eq("train").sum() == reg.region_role.eq("held").sum() == 64
    assert reg.groupby("input_component").region_role.nunique().eq(1).all()
    assert reg.groupby("input_component").inner_region_role.nunique().eq(1).all()
    assert reg.loc[reg.region_role.eq("held"), "inner_region_role"].eq("outer_held").all()
    mids = (reg.start0.to_numpy() + reg.end0.to_numpy()) // 2
    for i in range(len(reg)):
        for j in range(i):
            if reg.chrom.iloc[i] == reg.chrom.iloc[j] and abs(mids[i] - mids[j]) < 16384:
                assert reg.input_component.iloc[i] == reg.input_component.iloc[j]
    with np.load(FEATURES, allow_pickle=False) as features:
        np.testing.assert_array_equal(features["region_key"], reg.region_key)
        seq = features["mean_strands"]
    assert seq.shape[0] == 128 and np.isfinite(seq).all()
    feature_receipt = read_json(FEATURES.parent / "receipt.json")
    assert feature_receipt["whole_counted_interval_weighted128bp_pooling"] is True and feature_receipt["outcomes_read"] is False
    return x, y, fold, reg, genes, chrom, seq


def splits(fold, reg, held, mode, stage):
    donors = np.flatnonzero((fold != held) & ((fold != (held + 1) % 5) if stage == "inner" else True))
    if mode == "trained_regions":
        regions = np.arange(len(reg))
    else:
        regions = np.flatnonzero(reg.inner_region_role.eq("inner_train") if stage == "inner" else reg.region_role.eq("train"))
    return donors, regions


def verify_training_and_masks(run, x, y, fold, reg, genes, chrom):
    fits = read_json(run / "completed_fits.json")
    metadata = pd.read_csv(run / "gene_masks.tsv", sep="\t", dtype={"group": str})
    assert not metadata.duplicated(["mode", "held_fold", "stage", "arm", "group"]).any()
    selected = pd.read_csv(run / "inner_selection.tsv", sep="\t", dtype={"group": str})
    actual_fit_keys = {(r["mode"], r["held_fold"], r["arm"], r["group"], r["model"].split("|", 1)[1]) for r in fits}
    assert len(actual_fit_keys) == len(fits)
    expected_fits, expected_masks = set(), set()
    for mode in ("trained_regions", "held_regions"):
        target_regions = reg if mode == "trained_regions" else reg.loc[reg.region_role.eq("held")]
        assignments = [("all_rna", "all"), ("annotated_all_rna", "all")]
        assignments += [(arm, c) for c in sorted(target_regions.chrom.unique()) for arm in ("chromosome", "random_chromosome")]
        for held in range(5):
            for arm, group in assignments:
                expected_fits.update((mode, held, arm, group, model) for model in NEURAL)
                expected_masks.update((mode, held, stage, arm, group) for stage in ("inner", "outer"))
    assert actual_fit_keys == expected_fits
    assert set(metadata[["mode", "held_fold", "stage", "arm", "group"]].itertuples(index=False, name=None)) == expected_masks
    for row in fits:
        objective = row["model"].rsplit("_", 1)[1]
        assert objective in ("profile", "residual", "pairwise") and row["selected_steps"] in (100, 300, 1000)
        for stage in ("inner", "outer"):
            donors, regions = splits(fold, reg, row["held_fold"], row["mode"], stage)
            values = y[np.ix_(donors, regions)]
            cost = row[stage]
            assert cost["training_donors"] == len(donors) and cost["training_regions"] == len(regions)
            assert cost["steps"] == (1000 if stage == "inner" else row["selected_steps"])
            assert cost["objective"] == objective and cost["interaction"] == row["model"].split("|", 1)[1].startswith("film_")
            assert cost["normalization"] == "training_targets_only" and cost["pair_partners_per_donor_per_epoch"] == 2
            level = float(values.mean()) if objective == "profile" else 0.
            scale = max(float(values.std() if objective == "profile" else (values - values.mean(axis=0)).std()), 1e-6)
            np.testing.assert_allclose([cost["global_level"], cost["target_scale"]], [level, scale], rtol=0, atol=1e-12)
    with np.load(run / "gene_mask_membership.npz", allow_pickle=False) as source:
        np.testing.assert_array_equal(source["gene_ids"], genes)
        np.testing.assert_array_equal(source["gene_chrom"], chrom)
        admitted = source["admitted"]; removed = source["explicit_removed"]
    assert admitted.shape == removed.shape == (len(metadata), 42163)
    known = chrom != "unknown"
    masks, flags = [], []
    variance_cache = {}
    for row in metadata.itertuples(index=False):
        donors, _ = splits(fold, reg, row.held_fold, row.mode, row.stage)
        key = tuple(donors)
        if key not in variance_cache:
            variance_cache[key] = (x[donors].var(axis=0) > 1e-10, x[donors].std(axis=0) > 1e-10)
        variable, eligible = variance_cache[key]
        removal = np.zeros(len(genes), bool)
        if row.arm in ("chromosome", "random_chromosome"):
            removal = (chrom == row.group) & variable
            if row.arm == "random_chromosome":
                count = int(removal.sum()); removal[:] = False
                rng = np.random.default_rng(SEED + stable(row.arm + "|" + row.group))
                removal[rng.choice(np.flatnonzero(known & variable), count, replace=False)] = True
        keep = np.ones(len(genes), bool) if row.arm == "all_rna" else known & ~removal
        wanted = keep & eligible
        np.testing.assert_array_equal(admitted[row.mask_index], wanted)
        np.testing.assert_array_equal(removed[row.mask_index], removal)
        assert row.admitted_genes == wanted.sum() and row.removed_explicit_genes == removal.sum() and row.unknown_coordinate_genes == 37
        assert row.actual_NN_PCA_rank == min(32, row.actual_RNA_training_rank) and row.actual_RNA_training_rank <= len(donors)
        surviving = int((wanted & (chrom == row.group)).sum()) if row.arm == "chromosome" else 0
        if surviving:
            flags.append({"mode": row.mode, "held_fold": int(row.held_fold), "stage": row.stage, "group": row.group,
                          "surviving_target_chromosome_genes": surviving, "reason": "variance_removal_and_SD_admission_thresholds_differ"})
        masks.append({"mode": row.mode, "held_fold": row.held_fold, "stage": row.stage, "arm": row.arm, "group": row.group,
                      "admitted_genes": int(wanted.sum()), "removed_genes": int(removal.sum()), "surviving_target_chromosome_genes": surviving})
    check = pd.DataFrame(masks)
    left = check.loc[check.arm.eq("chromosome")].set_index(["mode", "held_fold", "stage", "group"])
    right = check.loc[check.arm.eq("random_chromosome")].set_index(["mode", "held_fold", "stage", "group"]).loc[left.index]
    np.testing.assert_array_equal(left.admitted_genes, right.admitted_genes)
    np.testing.assert_array_equal(left.removed_genes, right.removed_genes)
    assert set(selected["mode"]) == {"trained_regions", "held_regions"}
    assert set(selected.loc[selected["mode"].eq("held_regions"), "model"]) <= {"rna_global", "additive_linear", "product_kernel"}
    assert selected.held_fold.isin(range(5)).all() and selected.arm.isin(ARMS).all()
    return fits, selected, check, flags


def difference_mse(error, fold, weights=None):
    weights = np.ones(len(fold)) if weights is None else weights
    result = np.zeros(error.shape[1])
    for f in range(5):
        use = fold == f; w = weights[use]; n = w.sum()
        total = np.einsum("i,ij->j", w, error[use])
        squares = np.einsum("i,ij,ij->j", w, error[use], error[use])
        result += 2 * n / weights.sum() * (squares - total * total / n) / (n - 1)
    return result


def bootstrap(ea, eb, fold, components, rng, draws=1000):
    group_ids, inverse = np.unique(components, return_inverse=True)
    donor_groups = [np.flatnonzero(fold == f) for f in range(5)]
    squared_difference = ea ** 2 - eb ** 2
    values = np.empty((draws, 2))
    for i in range(draws):
        donors = np.concatenate([rng.choice(group, len(group), replace=True) for group in donor_groups])
        donor_weight = np.bincount(donors, minlength=len(fold))
        groups = rng.choice(group_ids, len(group_ids), replace=True)
        component_counts = np.bincount(np.searchsorted(group_ids, groups), minlength=len(group_ids))
        region_weight = component_counts[inverse]
        denominator = donor_weight.sum() * region_weight.sum()
        values[i, 0] = np.einsum("i,ij,j->", donor_weight, squared_difference, region_weight) / denominator
        diff = difference_mse(ea, fold, donor_weight) - difference_mse(eb, fold, donor_weight)
        values[i, 1] = np.dot(diff, region_weight) / region_weight.sum()
    assert np.isfinite(values).all()
    return values


def bh(p):
    p = np.asarray(p, float); order = np.argsort(p)
    corrected = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(p); result[order] = np.minimum(corrected, 1)
    return result


def contrasts(keys, mode):
    result = []
    anchor_kind = "selected_rna_linear" if mode == "trained_regions" else "product_kernel"
    for key in keys:
        arm, model = key.split("|", 1)
        anchor = arm + "|" + anchor_kind
        if key != anchor:
            result.append((key, anchor, "architecture_MSE"))
        if arm == "annotated_all_rna":
            result.append((key, "all_rna|" + model, "annotation_coverage_MSE"))
    for key in keys:
        arm, model = key.split("|", 1)
        if arm == "chromosome":
            result += [(key, "annotated_all_rna|" + model, "context_exclusion_MSE"),
                       (key, "random_chromosome|" + model, "context_exclusion_MSE")]
    return result


def self_checks():
    rng = np.random.default_rng(42)
    fold = np.repeat(np.arange(5), 4); error = rng.normal(size=(20, 5))
    expanded = np.concatenate([rng.choice(np.flatnonzero(fold == f), 4, replace=True) for f in range(5)])
    weights = np.bincount(expanded, minlength=20)
    expected = np.zeros(5)
    for f in range(5):
        e = error[expanded[fold[expanded] == f]]
        # Enumerate ordered distinct sample-position pairs, including identical resampled donors.
        pair = [(e[i] - e[j]) ** 2 for i in range(4) for j in range(4) if i != j]
        expected += np.mean(pair, axis=0) / 5
    np.testing.assert_allclose(difference_mse(error, fold, weights), expected, rtol=0, atol=1e-12)
    np.testing.assert_array_equal(bootstrap(error, error, fold, np.arange(5), np.random.default_rng(SEED), 3), np.zeros((3, 2)))
    assert len(contrasts([a + "|" + m for a in ARMS for m in BASELINES + REGION_RNA + NEURAL], "trained_regions")) == 101
    assert len(contrasts([a + "|" + m for a in ARMS for m in BASELINES + NEURAL], "held_regions")) == 73
    return ["weighted_donor_difference_equals_explicit_distinct_pairs", "identical_predictors_zero_bootstrap", "full_101_and73_contrast_families"]


def main(args):
    assert os.environ.get("SLURM_JOB_ID"), "Compute allocation required"
    args.out.mkdir(parents=True, exist_ok=False, mode=0o700)
    shutil.copy2(__file__, args.out / "executed_review_donor_context.py")
    shutil.copy2(Path(B.__file__), args.out / "executed_review_donor_baselines.py")
    started = time.monotonic(); tests = self_checks()
    completion = read_json(args.run / "complete.json") if (args.run / "complete.json").is_file() else {}
    if completion.get("status") != "completed_nested_development":
        emit(args.out / "checks.json", {"status": "not_evaluated_incomplete_donor_comparison", "producer_completion": completion,
             "statistical_comparisons_run": False, "all_planned_families_retained_no_success_only_comparison": True})
        return
    invariant = read_json(BASE / "model/donor_fit_invariants_21773708.json")
    assert sum(value is True for value in invariant.values()) == 10
    assert completion["no_evaluation_fitted_scaling"] is True and completion["no_H3_covariates"] is True
    recipe = read_json(args.run / "recipe.json")
    assert recipe["planned_paired_contrasts"] == 174 and recipe["seed"] == SEED and recipe["unknown_coordinate_gene_count"] == 37
    assert recipe["exclusion_levels"] == ["chromosome"] and recipe["planned_maximum_NN_fits"] == 5040
    precheck_path = BASE / "review/donor-mask-precheck-21774187/checks.json"
    precheck = read_json(precheck_path)
    assert precheck["status"] == "no_actual_mask_threshold_mismatch" and precheck["affected_partitions"] == 0
    assert precheck["partitions_checked"] == 390 and precheck["all_effective_random_dimensions_match"] is True
    x, y, fold, reg, genes, chrom, seq = load_source()
    fits, selected, masks, flags = verify_training_and_masks(args.run, x, y, fold, reg, genes, chrom)
    assert completion["completed_NN_fits"] == len(fits)
    predictions = {"trained_regions": {}, "held_regions": {}}
    with np.load(args.run / "internal_predictions.npz", allow_pickle=False) as saved:
        for key in saved.files:
            mode, model = key.split("|", 1)
            predictions[mode][model] = saved[key]
    differences, metrics, rows = {}, [], []
    for mode, values in predictions.items():
        expected_models = {a + "|" + m for a in ARMS for m in BASELINES + NEURAL + (REGION_RNA if mode == "trained_regions" else ())}
        assert set(values) == expected_models
        regions = np.arange(len(reg)) if mode == "trained_regions" else np.flatnonzero(reg.region_role.eq("held"))
        target = y[:, regions]; blocks = reg.input_component.to_numpy()[regions]
        deposited_metrics = pd.read_csv(args.run / f"{mode}_metrics.tsv", sep="\t").set_index("model")
        assert set(deposited_metrics.index) == expected_models
        for key, matrix in values.items():
            assert matrix.shape == y.shape and np.isfinite(matrix[:, regions]).all(), "Completed model has missing/nonfinite evaluation values"
            if mode == "held_regions":
                assert np.isnan(matrix[:, reg.region_role.eq("train")]).all()
            pred = matrix[:, regions]; error = pred - target
            rho = [spearmanr(pred[:, j], target[:, j]).statistic for j in range(len(regions)) if np.ptp(pred[:, j]) > 1e-10]
            points = [np.sqrt(np.mean(error ** 2)), difference_mse(error, fold).mean(), np.mean(rho) if rho else np.nan]
            original = deposited_metrics.loc[key]
            assert original.donors == 99 and original.regions == len(regions) and original.genomic_components == len(set(blocks))
            close(points, original[["RMSE_log2CPM", "donor_difference_MSE", "mean_region_donor_spearman"]].to_numpy(float), "all_point_metrics", differences)
            anchor = values[key.split("|")[0] + "|training_baseline"][:, regions] - target
            metrics.append({"evaluation": mode, "model": key, "donors": 99, "regions": len(regions), "RMSE_log2CPM": points[0],
                "donor_difference_MSE": points[1], "mean_region_donor_spearman": points[2],
                "MSE_skill_vs_training_baseline": 1 - points[0] ** 2 / np.mean(anchor ** 2),
                "donor_difference_skill_vs_training_baseline": 1 - points[1] / difference_mse(anchor, fold).mean()})
        pairs = contrasts(list(values), mode)
        expected_counts = {"architecture_MSE": 56 if mode == "trained_regions" else 40,
                           "annotation_coverage_MSE": 15 if mode == "trained_regions" else 11,
                           "context_exclusion_MSE": 30 if mode == "trained_regions" else 22}
        original = pd.read_csv(args.run / f"{mode}_paired_uncertainty.tsv", sep="\t")
        assert [(r.model, r.comparator, r.family) for r in original.itertuples()] == pairs
        assert original.status.eq("development_fixed_predictions").all()
        rng = np.random.default_rng(SEED)
        begin = len(rows)
        for key, comparator, family in pairs:
            ea = values[comparator][:, regions] - target; eb = values[key][:, regions] - target
            point = np.mean(ea ** 2 - eb ** 2)
            ddpoint = (difference_mse(ea, fold) - difference_mse(eb, fold)).mean()
            draws = bootstrap(ea, eb, fold, blocks, rng)
            low, high = np.quantile(draws, [.025, .975], axis=0)
            p = (1 + np.sum(np.abs(draws[:, 0] - point) >= abs(point))) / 1001
            rows.append({"evaluation": mode, "model": key, "comparator": comparator, "family": family,
                         "MSE_improvement": point, "low95": low[0], "high95": high[0],
                         "donor_difference_MSE_improvement": ddpoint, "difference_low95": low[1], "difference_high95": high[1], "p_family_accounting": p})
        for family, count in expected_counts.items():
            indices = [i for i in range(begin, len(rows)) if rows[i]["family"] == family]
            assert len(indices) == count
            q = bh([rows[i]["p_family_accounting"] for i in indices])
            for i, adjusted in zip(indices, q):
                rows[i]["BH_q_complete_mode_family"] = adjusted; rows[i]["planned_family_n"] = count
        rebuilt = pd.DataFrame(rows[begin:])
        for column in ("MSE_improvement", "low95", "high95", "donor_difference_MSE_improvement", "difference_low95", "difference_high95", "p_family_accounting", "BH_q_complete_mode_family", "planned_family_n"):
            close(rebuilt[column], original[column], column, differences, tolerance=1e-10)
    assert len(rows) == 174 and len(metrics) == 104
    pd.DataFrame(metrics).to_csv(args.out / "independent_point_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(args.out / "independent_paired_uncertainty.tsv", sep="\t", index=False)
    masks.to_csv(args.out / "independent_mask_counts.tsv", sep="\t", index=False)
    analytic = B.reconstruct(x, seq, y, fold, reg, predictions, selected, fits, chrom, seconds_cap=3000)
    emit(args.out / "selected_analytical_baseline_checks.json", analytic)
    emit(args.out / "checks.json", {"status": "independent_donor_context_reconstruction_passed" if not flags else "donor_metrics_reconstructed_exclusion_interpretation_limited",
        "job_id": os.environ["SLURM_JOB_ID"], "checks": tests, "source_scientific_invariants_true": 10,
        "source_donors": 99, "trained_output_regions": 128, "crossed_held_output_regions": 64,
        "point_models": 104, "paired_contrasts": 174, "complete_family_counts": {"trained_regions": [56, 15, 30], "held_regions": [40, 11, 22]},
        "differences": differences, "mask_flags": flags, "unrestricted_gene_universe": 42163, "annotated_gene_universe": 42126,
        "actual_input_mask_precheck_sha256": sha256(precheck_path), "actual_input_mask_precheck_partitions": 390,
        "fit_record_count_inner_outer_pairs": len(fits), "individual_NN_fits": 2 * len(fits), "planned_maximum_individual_NN_fits": 5040,
        "producer_completed_NN_fits_field_semantics": "counts_inner_outer_pairs; multiply_by_two_for_individual_network_fits",
        "analytical_baseline_check_status": analytic["status"], "NN_or_inner_selection_refit_performed": False,
        "bootstrap": "Independent donor multiplicity and region-component weights; donor-difference sum-of-squares; same RNG strata and full families",
        "uncertainty_scope": "fixed_predictions_only; no_fitting_or_selection_uncertainty; genomic_input_components_not_verified_LD_independent",
        "scientific_scope": "same99_source_participants_small128_region_development; no_trans_causality_or_independent_cohort_claim",
        "normalization_scope": "full_deposited_RNA_and_H3_libraries; fractionalRNA_continuous; excludedgenes_aggregate_denominator_retained",
        "protected_outcomes_read": False, "participant_level_values_written": False, "seconds": time.monotonic() - started,
        "source_sha256": {str(path): sha256(path) for path in (FIXTURE / "paired.npz", FIXTURE / "regions.tsv", FEATURES, args.run / "internal_predictions.npz", args.run / "recipe.json")}})
    assert sum(path.stat().st_size for path in args.out.iterdir()) < 100000000


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        main(args)
    except Exception as exc:
        if args.out.is_dir():
            emit(args.out / "checks.json", {"status": "independent_donor_context_check_failed", "error": repr(exc), "successful_comparison_claim": False})
        raise
