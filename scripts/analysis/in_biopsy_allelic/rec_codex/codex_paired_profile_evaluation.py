"""Fixed donor-level RNA-to-H3 profile comparison; no fit or input adapter.

The command runs synthetic checks only. Receiving evaluation still requires
the declared complete measured RNA/counting contract and frozen B99 labels.
"""
import json
import os
from pathlib import Path
import sys

import numpy as np

from profile_error_decomposition import profile_error_decomposition
from export_training_h3_target_transform import FIX, read_axis, require, sha

DONORS = ("B1", "B7", "B8", "B9", "B10", "B12", "B15", "B19", "B21",
          "B22", "B24", "B26", "B36", "B38", "B41", "B46", "B47")
REGIONS = 96460
REGION_AXIS_SHA = "cbe35aeb188e531283a1bcd636ba14de4cfbf215dfc544b99e3ebb97f9831986"
BOOTSTRAPS = 10000
SEED = 20260930


def paired_loss_summary(target, global_prediction, local_prediction, training_mean):
    """One paired-donor primary contrast, conditional on frozen predictions.

    Arrays must already use identical donor and region order. This function
    does not establish their identity or measurement comparability.
    """
    target = np.asarray(target, dtype=np.float64)
    global_prediction = np.asarray(global_prediction, dtype=np.float64)
    local_prediction = np.asarray(local_prediction, dtype=np.float64)
    training_mean = np.asarray(training_mean, dtype=np.float64)
    require(target.ndim == 2 and target.shape[0] >= 2 and target.shape[1] >= 1,
            "Require at least two donors and one region")
    require(global_prediction.shape == local_prediction.shape == target.shape,
            "Both predictions must preserve the entire common target axis")
    require(training_mean.shape == (target.shape[1],), "Training mean axis differs")
    require(all(np.isfinite(x).all() for x in (target, global_prediction, local_prediction, training_mean)),
            "Nonfinite values: no model-specific or common complete-case omission")
    decompositions = {"global": profile_error_decomposition(target, global_prediction),
                      "local": profile_error_decomposition(target, local_prediction)}
    losses = np.column_stack([decompositions["global"]["donor_mse"],
                              decompositions["local"]["donor_mse"],
                              np.mean((target - training_mean) ** 2, axis=1),
                              np.mean(target ** 2, axis=1)])
    require(np.isfinite(losses).all(), "Squared error overflow")
    delta = losses[:, 0] - losses[:, 1]
    # Resample biological donors once per draw, sharing indices across arms.
    draws = np.random.default_rng(SEED).integers(0, len(delta), size=(BOOTSTRAPS, len(delta)))
    boot_delta = delta[draws].mean(axis=1)
    interval = np.quantile(boot_delta, [0.025, 0.975], method="linear")
    mse = losses.mean(axis=0)
    baseline = float(mse[2])
    skill = {name: None if baseline == 0 else float(1 - mse[j] / baseline)
             for j, name in enumerate(("global", "local"))}
    return dict(biological_n=len(delta), regions=target.shape[1],
                donor_mse=dict(zip(("global", "local", "fixed_training_mean", "zero_secondary"), losses.T.tolist())),
                donor_paired_error_reduction=delta.tolist(),
                mean_paired_error_reduction=float(delta.mean()),
                paired_error_reduction_percentile95=interval.tolist(),
                cohort_mse=dict(zip(("global", "local", "fixed_training_mean", "zero_secondary"), mse.tolist())),
                descriptive_skill_vs_fixed_training_mean=skill,
                skill_definition="1 - mean donor model MSE / mean donor fixed-training-mean MSE; undefined when baseline MSE is zero",
                error_decomposition=decompositions,
                bootstrap=dict(draws=BOOTSTRAPS, seed=SEED, unit="biological donor",
                               paired_indices_shared=True, method="percentile, linear quantiles",
                               uncertainty="approximate interval conditional on frozen models/B99/counting and measured labels; excludes refitting, technical error and assay shift"),
                hypothesis_family="one prespecified primary paired mean-error contrast; no region tests or p-values",
                target_or_prediction_recalibrated=False, model_selected=False,
                independently_validated=False)


def evaluate_fixed_arrays(target, global_prediction, local_prediction, training_mean,
                          target_donors, global_donors, local_donors,
                          target_regions, global_regions, local_regions, training_mean_regions):
    """Enforce the fixed external axes without acquiring or opening any files.

    Caller must independently establish measurement admission and construct
    the complete fixed B99 residual target. Axes here cannot establish that.
    """
    require(os.environ.get("SLURM_JOB_ID"), "Receiving numerical work requires compute")
    require(all(tuple(axis) == DONORS for axis in (target_donors, global_donors, local_donors)),
            "Require the same exact seventeen donors and order in all arrays")
    keys = tuple(target_regions)
    require(len(keys) == len(set(keys)) == REGIONS, "Require all unique fixed region keys")
    require(all(tuple(axis) == keys for axis in (global_regions, local_regions, training_mean_regions)),
            "Prediction/baseline region identity or ordering differs")
    axis_file = FIX / "h3k27ac_feature_axis.tsv"
    require(sha(axis_file) == REGION_AXIS_SHA, "Frozen region identity changed")
    expected_keys = tuple(row["opaque_source_feature_key"] for row in read_axis(axis_file))
    require(keys == expected_keys, "Supplied region axis differs from fixed counted intervals")
    require(np.shape(target) == (len(DONORS), REGIONS), "Fixed target shape changed")
    result = paired_loss_summary(target, global_prediction, local_prediction, training_mean)
    result["donor_id"] = list(DONORS)
    result["fixed_region_axis_sha256"] = REGION_AXIS_SHA
    return result


def check_invariants():
    require(os.environ.get("SLURM_JOB_ID"), "Execute synthetic checks on compute")
    # Independent arithmetic: global donor losses .5,6.5,20.5; local loss 1.
    # Mean reduction 49/6 and skill 49/55, with equal donor/region weighting.
    target = np.array([[0., 1.], [2., 3.], [4., 5.]])
    global_prediction = np.zeros_like(target)
    local_prediction = target - 1
    baseline = np.zeros(2)
    report = paired_loss_summary(target, global_prediction, local_prediction, baseline)
    require(np.isclose(report["mean_paired_error_reduction"], 49 / 6, rtol=0, atol=1e-12), "Hand-calculated paired loss differs")
    require(np.isclose(report["descriptive_skill_vs_fixed_training_mean"]["local"], 49 / 55, rtol=0, atol=1e-12), "Hand-calculated skill differs")
    # Identical regional duplication cannot manufacture more biological n or precision.
    repeated = paired_loss_summary(np.tile(target, (1, 4)), np.tile(global_prediction, (1, 4)),
                                   np.tile(local_prediction, (1, 4)), np.tile(baseline, 4))
    require(report["biological_n"] == repeated["biological_n"] == 3, "Region duplication changed donor n")
    for key in ("mean_paired_error_reduction", "paired_error_reduction_percentile95"):
        require(np.array_equal(report[key], repeated[key]), "Region duplication changed paired uncertainty")
    require(report == paired_loss_summary(target, global_prediction, local_prediction, baseline), "Seeded repetition differs")
    zeros = np.zeros((17, 2))
    perfect = paired_loss_summary(zeros, zeros, zeros, np.zeros(2))
    require(perfect["paired_error_reduction_percentile95"] == [0., 0.], "Identical predictors generated uncertainty")
    require(perfect["descriptive_skill_vs_fixed_training_mean"] == {"global": None, "local": None}, "Zero baseline acquired an invented skill")
    swapped = paired_loss_summary(target, local_prediction, global_prediction, baseline)
    require(np.allclose(swapped["paired_error_reduction_percentile95"], -np.array(report["paired_error_reduction_percentile95"])[::-1], rtol=0, atol=1e-12), "Swapped paired interval fails sign symmetry")
    for bad_global in (global_prediction[:-1], np.full_like(global_prediction, np.nan)):
        try:
            paired_loss_summary(target, bad_global, local_prediction, baseline)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid common observations were accepted")
    # Identity refusal must precede shape/data work; no receiving arrays used.
    for donor_axis in (DONORS[::-1], DONORS[:-1] + (DONORS[0],)):
        try:
            evaluate_fixed_arrays(target, global_prediction, local_prediction, baseline,
                                  donor_axis, DONORS, DONORS, (), (), (), ())
        except ValueError as error:
            require("seventeen donors" in str(error), "Donor identity was not checked first")
        else:
            raise AssertionError("Incorrect donor identity was accepted")
    # Exercise the full declared identity path with synthetic values. Only
    # permitted region metadata is read; no molecular/count array is opened.
    region_file = FIX / "h3k27ac_feature_axis.tsv"
    require(sha(region_file) == REGION_AXIS_SHA, "Frozen synthetic-check metadata changed")
    keys = tuple(row["opaque_source_feature_key"] for row in read_axis(region_file))
    full_zero = np.broadcast_to(0., (len(DONORS), REGIONS))
    full_one = np.broadcast_to(1., (len(DONORS), REGIONS))
    fixed = evaluate_fixed_arrays(full_zero, full_zero, full_one, np.zeros(REGIONS),
                                  DONORS, DONORS, DONORS, keys, keys, keys, keys)
    require(fixed["donor_id"] == list(DONORS), "Per-donor errors lack their identities")
    require(fixed["mean_paired_error_reduction"] == -1 and
            fixed["paired_error_reduction_percentile95"] == [-1., -1.],
            "Hand-calculated full-axis paired comparison differs")
    for bad_keys in (keys[:-1] + (keys[0],), keys[::-1]):
        try:
            evaluate_fixed_arrays(full_zero, full_zero, full_one, np.zeros(REGIONS),
                                  DONORS, DONORS, DONORS, keys, keys, keys, bad_keys)
        except ValueError:
            pass
        else:
            raise AssertionError("Changed/duplicate baseline region axis was accepted")
    output = Path(os.environ["CODEX_REC_OUTPUT"]) / ("codex_paired_profile_checks_" + os.environ["SLURM_JOB_ID"])
    output.mkdir(exist_ok=False)
    summary = dict(hand_calculated_paired_loss=True, seeded_donor_bootstrap_repeat=True,
                   region_duplication_preserves_n_and_uncertainty=True,
                   swapped_arm_sign_symmetry=True, identical_predictions_zero_interval=True,
                   zero_baseline_skill_undefined=True, donor_identity_mismatch_refused=True,
                   nonfinite_or_partial_common_axis_refused=True,
                   full_fixed_axis_synthetic_hand_calculation=True,
                   returned_donor_ids_verified=True, changed_baseline_region_axis_refused=True,
                   receiving_data_read=False, model_fitted=False, scientific_validation=False,
                   test_examples_biological_n=[3, 17], fixed_receiving_donors=list(DONORS),
                   fixed_receiving_regions=REGIONS, primary_example=report,
                   script_sha256=sha(Path(__file__)),
                   decomposition_helper_sha256=sha(Path(__file__).with_name("profile_error_decomposition.py")),
                   imported_helper_sha256=sha(Path(__file__).with_name("export_training_h3_target_transform.py")),
                   launcher_sha256=sha(Path(__file__).with_name("run_codex_paired_profile_checks.sbatch")),
                   environment=dict(python=sys.version, numpy=np.__version__, slurm_job_id=os.environ["SLURM_JOB_ID"]))
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k:v for k,v in summary.items() if isinstance(v, bool)}, sort_keys=True))


if __name__ == "__main__":
    check_invariants()
