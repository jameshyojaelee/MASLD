#!/usr/bin/env python3
"""Independent denominator and participant-level checks of chromatin candidates.

No predictive model is fitted. The response-only concentration adjustment is
reconstructed to verify that the residual experiment and saved target agree.
Intervals condition on the existing fitted models and observed region panels.
"""
import argparse
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT / "Analysis/MASLD_Model_Benchmark"
FIX = BENCH / "executions/model-data-064-21079902/fixture"
CANDIDATE = ROOT / "GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT"
PAIRED = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model/donor_curve_fixtures_21784771/size_8192"
OOF = BENCH / "executions/chromatin-stable-rrr-forms-20260908T192024Z/out/stable_oof.npz"
WEIGHTS = BENCH / "release/masld-liver-chromatin-state-v1.2/weights/chromatin_state_v1_2.npz"
SEQUENCE = CANDIDATE / "sequence_rna_residual_21849236"
ACQUISITION = CANDIDATE / "chromatin_acquisition_21849197"
SEED = 20260924
DRAWS = 20000


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def concentration(raw):
    """Reconstruct the declared four sample concentration descriptors."""
    library = raw.sum(axis=1)
    fraction = raw / library[:, None]
    positive_log = np.log(np.maximum(fraction, np.finfo(float).tiny))
    entropy = -(fraction * positive_log).sum(axis=1)
    sorted_values = np.sort(raw, axis=1)
    width = raw.shape[1]
    rank_weight = 2 * np.arange(1, width + 1) - width - 1
    gini = sorted_values @ rank_weight / (width * library)
    iqr = np.array([np.diff(np.quantile(np.log2(row[row > 0]), [.25, .75]))[0]
                    for row in raw])
    top_fraction = sorted_values[:, -round(.05 * width):].sum(axis=1) / library
    return np.column_stack([gini, entropy, iqr, top_fraction])


def interval(draws):
    return [float(x) for x in np.quantile(draws, [.025, .975])]


def holm(values):
    values = np.asarray(values, float)
    order = np.argsort(values)
    adjusted = np.empty_like(values)
    corrected = np.maximum.accumulate(values[order] * np.arange(len(values), 0, -1))
    adjusted[order] = np.minimum(corrected, 1.0)
    return adjusted


def sequence_check(axis, fold, sampled, out):
    regions = pd.read_csv(PAIRED / "regions.tsv", sep="\t")
    subset = regions.loc[regions.training_1024 & regions.region_role.eq("train")]
    if len(subset) != 1024 or subset.region_index.duplicated().any():
        raise ValueError("Sequence region census differs")
    selected = subset.region_index.to_numpy(int)
    with np.load(PAIRED / "paired.npz", allow_pickle=True) as source:
        if not np.array_equal(source["participant_ids"].astype(str), axis):
            raise ValueError("Sequence fixture and raw participant axis differ")
        if not np.array_equal(source["donor_fold"].astype(int), fold):
            raise ValueError("Sequence fixture and saved held-participant folds differ")
    with np.load(WEIGHTS, allow_pickle=True) as weights:
        if not np.array_equal(weights["prof_region_key"].astype(str)[selected],
                              subset.region_key.astype(str).to_numpy()):
            raise ValueError("Sequence region identity differs from profile region axis")
    with np.load(PAIRED / "features.npz", allow_pickle=False) as features:
        if not np.array_equal(features["region_key"].astype(str), regions.region_key.astype(str)):
            raise ValueError("Frozen sequence feature region order differs")
    with np.load(OOF, allow_pickle=False) as saved:
        if not np.array_equal(saved["fold"].astype(int), fold):
            raise ValueError("Stored response folds differ")
        measured = saved["Cres_oof"][:, selected].astype(float)

    raw = np.asarray(np.load(FIX / "molecular/h3k27ac_counts.npy", mmap_mode="r"), float)
    if raw.shape != (99, 96460):
        raise ValueError("Measured chromatin axis differs")
    response = np.log2(raw[:, selected] / raw.sum(axis=1)[:, None] * 1e6 + 1)
    descriptor = concentration(raw)
    design = np.column_stack([np.ones(len(axis)), descriptor])
    reference = np.zeros(len(axis))
    reconstruction_error, largest_training_mean = 0.0, 0.0
    for outer in sorted(set(fold)):
        train, held = fold != outer, fold == outer
        coefficient = np.linalg.lstsq(design[train], response[train], rcond=None)[0]
        train_residual = response[train] - design[train] @ coefficient
        held_residual = response[held] - design[held] @ coefficient
        reconstruction_error = max(reconstruction_error,
                                   float(np.max(np.abs(held_residual - measured[held]))))
        train_mean = train_residual.mean(axis=0)
        largest_training_mean = max(largest_training_mean, float(np.max(np.abs(train_mean))))
        reference[held] = np.sum((measured[held] - train_mean) ** 2, axis=1)
    if reconstruction_error > 4e-5 or largest_training_mean > 1e-8:
        raise ValueError("Original concentration-residual response was not reproduced")

    error = pd.read_csv(SEQUENCE / "held_participant_errors.tsv", sep="\t")
    if error.participant_id.duplicated().any() or set(error.participant_id) != set(axis):
        raise ValueError("Sequence error-table participant census differs")
    error = error.set_index("participant_id").loc[axis].reset_index()
    if not np.array_equal(error.fold.to_numpy(int), fold):
        raise ValueError("Sequence error-table held folds differ")
    if not np.isfinite(error.filter(like="_SSE").to_numpy()).all():
        raise ValueError("Nonfinite sequence errors")
    numerator = error.additive_SSE.to_numpy() - error.interaction_SSE.to_numpy()
    residual_reference = error.baseline_SSE.to_numpy()
    numerator_draws = numerator[sampled].sum(axis=1)
    original_draws = numerator_draws / reference[sampled].sum(axis=1)
    residual_draws = numerator_draws / residual_reference[sampled].sum(axis=1)
    error["original_target_reference_SSE"] = reference
    error["additive_minus_interaction_SSE"] = numerator
    error.to_csv(out / "sequence_participant_denominators.tsv", sep="\t", index=False)
    region_scope = subset[["region_index", "region_key", "coordinate_semantics", "coordinate_residual_bp"]].copy()
    region_scope["called_peak_boundary_status"] = "unresolved"
    region_scope.to_csv(out / "sequence_regions.tsv", sep="\t", index=False)
    models = [c.removesuffix("_SSE") for c in error if c.endswith("_SSE")
              and c not in ("original_target_reference_SSE", "additive_minus_interaction_SSE")]
    return {
        "n_participants": len(axis), "n_regions": len(selected),
        "target_reconstruction_max_absolute_difference": reconstruction_error,
        "training_target_mean_max_absolute_value": largest_training_mean,
        "residual_baseline_SSE": float(residual_reference.sum()),
        "original_target_reference_SSE": float(reference.sum()),
        "original_skill": {model: float(1 - error[model + "_SSE"].sum() / reference.sum())
                           for model in models},
        "interaction_minus_additive_original_skill": float(numerator.sum() / reference.sum()),
        "interaction_minus_additive_original_skill_ci95": interval(original_draws),
        "interaction_minus_additive_residual_skill": float(numerator.sum() / residual_reference.sum()),
        "interaction_minus_additive_residual_skill_ci95": interval(residual_draws),
        "original_skill_complexity_margin": .01,
        "original_skill_complexity_margin_met": bool(numerator.sum() / reference.sum() >= .01),
        "interpretation": "Both gains have the identical paired SSE numerator; only the reference denominator differs. The negative increment does not qualify for added complexity.",
    }


def acquisition_check(axis, fold, sampled, out):
    all_errors = pd.read_csv(ACQUISITION / "acquisition_per_participant.tsv.gz", sep="\t")
    errors = all_errors.loc[np.isclose(all_errors.budget_fraction, .2)].copy()
    expected = {"random", "histology", "rna_diversity", "global_local_disagreement"}
    if set(errors.policy) != expected or set(errors.participant_id) != set(axis):
        raise ValueError("Acquisition participant or policy census differs")
    if errors.duplicated(["participant_id", "policy", "random_repeat"]).any():
        raise ValueError("Repeated acquisition error rows")
    for participant, rows in errors.groupby("participant_id"):
        if len(rows) != 6 or set(rows.loc[rows.policy.eq("random"), "random_repeat"]) != {0, 1, 2}:
            raise ValueError(f"Unexpected policy/repeat census: {participant}")
        if rows.outer_fold.nunique() != 1 or not np.allclose(rows.base_SSE, rows.base_SSE.iloc[0], atol=1e-10, rtol=0):
            raise ValueError("Comparisons differ in held fold or base prediction")
    observed_fold = errors.groupby("participant_id").outer_fold.first().loc[axis].to_numpy(int)
    if not np.array_equal(observed_fold, fold):
        raise ValueError("Acquisition held folds differ from fixture")
    nregions = len(pd.read_csv(ACQUISITION / "fixed_region_panel.tsv", sep="\t"))
    if nregions != 1000:
        raise ValueError("Acquisition region-panel census differs")
    averaged = errors.groupby(["participant_id", "policy"]).agg(
        refit_SSE=("refit_SSE", "mean"), base_SSE=("base_SSE", "mean")).reset_index()
    refit = averaged.pivot(index="participant_id", columns="policy", values="refit_SSE").loc[axis]
    base = averaged.groupby("participant_id").base_SSE.first().loc[axis]
    if not np.isfinite(refit.to_numpy()).all() or not np.isfinite(base.to_numpy()).all():
        raise ValueError("Nonfinite acquisition errors")
    participant = refit.rename(columns=lambda name: name + "_refit_SSE").copy()
    participant["outer_fold"] = fold
    participant["base_SSE"] = base
    rows = []
    for policy in ("histology", "rna_diversity", "global_local_disagreement"):
        gain = (refit["random"].to_numpy() - refit[policy].to_numpy()) / nregions
        observed = float(gain.mean())
        draws = gain[sampled].mean(axis=1)
        nominal_p = float((1 + np.count_nonzero(np.abs(draws - observed) >= abs(observed))) / (len(draws) + 1))
        participant[policy + "_MSE_gain_vs_random"] = gain
        rows.append({"policy": policy, "n_participants": len(axis),
                     "participant_pooled_MSE_gain_vs_random": observed,
                     "ci95_low": interval(draws)[0], "ci95_high": interval(draws)[1],
                     "p_approx_centered_participant_bootstrap": nominal_p})
    table = pd.DataFrame(rows)
    table["p_holm_three_policy_family"] = holm(table.p_approx_centered_participant_bootstrap)
    table.to_csv(out / "acquisition_primary20_paired_participant.tsv", sep="\t", index=False)
    participant.to_csv(out / "acquisition_primary20_participant_errors.tsv", sep="\t", index_label="participant_id")
    per_fold = pd.read_csv(ACQUISITION / "acquisition_per_fold.tsv", sep="\t")
    census = per_fold.loc[np.isclose(per_fold.budget_fraction, .2),
                         ["outer_fold", "n_base", "n_pool", "n_acquired", "n_evaluation"]].drop_duplicates()
    if len(census) != 5 or not (census.n_acquired == 20).all():
        raise ValueError("Primary acquisition budget census differs")
    return {
        "n_participants": len(axis), "n_regions": nregions, "n_random_repeats_averaged_per_participant": 3,
        "budget": "20 of 99 total participants; 20 of 33 to 37 selectable pool participants per fold",
        "fold_census": census.to_dict("records"),
        "participant_pooled_base_MSE": float(base.mean() / nregions),
        "participant_pooled_MSE_reduction": {policy: float((base - refit[policy]).mean() / nregions)
                                              for policy in sorted(expected)},
        "contrasts": table.to_dict("records"),
        "weighting": "Each biological participant has equal weight; differs from unweighted means of five unequal-size evaluation folds.",
    }


def main(out):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    if out.exists():
        raise FileExistsError(out)
    axis = pd.read_csv(FIX / "molecular/participant_axis.tsv", sep="\t").participant_id.astype(str).to_numpy()
    if len(axis) != 99 or len(set(axis)) != 99:
        raise ValueError("Expected 99 distinct participants")
    fold = pd.read_csv(FIX / "folds/participant_outer_folds.tsv", sep="\t").set_index("participant_id").loc[axis, "outer_fold"].to_numpy(int)
    sampled = np.random.default_rng(SEED).integers(0, len(axis), size=(DRAWS, len(axis)))
    out.mkdir(parents=True)
    report = {
        "status": "independent numerical follow-up of inspected single-cohort development",
        "sequence": sequence_check(axis, fold, sampled, out),
        "acquisition": acquisition_check(axis, fold, sampled, out),
        "uncertainty": "Paired resampling of 99 distinct participants, conditional on the existing fitted predictions and observed fixed region panels. Shared training sets induce fitted-model dependence; retraining, region-population and random-selection Monte Carlo uncertainty are not included. Centered-bootstrap P values are approximate conditional sensitivities, not exact inferential tests. Earlier fold sign-flips also remain dependence sensitivities.",
        "bootstrap_seed": SEED, "bootstrap_draws": DRAWS,
        "source_sha256": {str(path): sha256(path) for path in (
            Path(__file__), SEQUENCE / "held_participant_errors.tsv",
            ACQUISITION / "acquisition_per_participant.tsv.gz", OOF,
            PAIRED / "regions.tsv", FIX / "molecular/participant_axis.tsv",
            FIX / "folds/participant_outer_folds.tsv")},
        "environment": {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
                        "platform": platform.platform(), "slurm_job_id": os.environ["SLURM_JOB_ID"]},
    }
    (out / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"sequence": report["sequence"], "acquisition": report["acquisition"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=CANDIDATE / "chromatin_secondary_review")
    main(parser.parse_args().out)
