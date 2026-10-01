#!/usr/bin/env python3
"""Compare source-SE and constant-SE Gaussian scalar-effect training.

Twenty fixed-exposure development runs reuse ten existing MSE comparators.
SEs are source marginal-association uncertainty, not causal-effect uncertainty.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from catalog_build import sha
from catalog_currin import BASE, LABELS, SOURCE_SE, matched_uncertainty

ARMS = [f"adapter_r{rank}_last{blocks}_{pool}" for rank in (4,16)
        for blocks in (3,5) for pool in ("variant","symmetric")]
ARMS += [f"frozen_2048_{pool}_shared_mlp64" for pool in ("variant","symmetric")]
OLD = BASE / "model/sweep_21772726"


def emit(path, value):
    path.write_text(json.dumps(value,indent=2)+"\n")


def constant_training_se(labels):
    return float(labels.loc[labels.heldout_fold.isin([2,3,4]),"verified_se"].median())


def prepare(out):
    labels = pd.read_csv(LABELS,sep="\t")
    source = pd.read_csv(SOURCE_SE,sep="\t",usecols=["variant_id","target_id","beta_nominal","varbeta","effect_se","uncertainty_state"])
    joined = matched_uncertainty(labels,source)
    assert len(joined) == 32322 and joined.lead_variant_id.is_unique
    pd.testing.assert_frame_equal(joined[labels.columns],labels,check_exact=True)
    assert joined.se_status.eq("verified_exact_variant_peak_native_varbeta").all()
    assert np.isfinite(joined.verified_se).all() and joined.verified_se.gt(0).all()
    constant = constant_training_se(joined)
    changed = joined.copy()
    changed.loc[changed.heldout_fold.isin([0,1]),"verified_se"] *= 1e6
    assert constant_training_se(changed) == constant
    joined["constant_se"] = constant
    path = out/"labels_with_exact_source_se.tsv.gz"
    joined.to_csv(path,sep="\t",index=False)
    reread = pd.read_csv(path,sep="\t")
    np.testing.assert_array_equal(reread.beta_alt.to_numpy(np.float32),labels.beta_alt.to_numpy(np.float32))
    np.testing.assert_array_equal(reread.verified_se.to_numpy(np.float32),joined.verified_se.to_numpy(np.float32))
    config_dir = out/"configs"
    config_dir.mkdir()
    configs = []
    for arm in ARMS:
        old_config = OLD/arm/"config.json"
        original = json.loads(old_config.read_text())
        assert original["loss"] == "mse" and original["seed"] == 1103 and original["length"] == 2048
        for kind, column in (("source_se","verified_se"),("constant_se","constant_se")):
            config = {**original,"id":arm+"__gaussian_"+kind,"loss":"gaussian_effect","steps":5000,
                "status":"fixed_development_loss_comparison_before_new_fits"}
            filename = config_dir/(config["id"]+".json")
            emit(filename,config)
            configs.append(dict(base_arm=arm,kind=kind,se_column=column,path=str(filename),
                name=config["id"],MSE_comparator=str(OLD/arm),original_config_sha256=sha(old_config)))
    recipe = dict(population="32322_existing_significance_selected_Currin_leads; no_new_eligibility_selection",
        training_folds=[2,3,4],validation_fold=1,closed_fold=0,training_rows=18914,validation_rows=6834,
        seed=1103,steps=5000,microbatch=4,configuration_pairs=configs,
        source_se="sqrt_native_nominal_varbeta_matched_by_exact_variant_peak_and_beta",
        constant_se=constant,constant_se_definition="median_verified_SE_in_training_folds_2_3_4_only",
        loss="0.5*mean((prediction-beta)^2/(SE^2+softplus(log_residual_variance)+1e-8)+log(total_variance))",
        initial_log_residual_variance=-2.,initial_residual_variance=float(np.logaddexp(0,-2.)),
        optimization="unchanged_AdamW_3e-4_head_adapter_weight_decay1e-4_global_clip1; no_validation_updates",
        comparisons="source_SE_Gaussian_vs_MSE_and_source_SE_Gaussian_vs_constant_SE_Gaussian; fixed_ten_per_family",
        inference="paired_historical_bins_and_chromosome_sensitivity_after_complete_runs; no_finalist_selection",
        assumptions="approximate_Gaussian_marginal_beta_estimates_plus_shared_training_fitted_residual_variance",
        limits=["No correction for LD, significance selection, winner's curse, or causal interpretation.",
            "MSE comparison changes the whole likelihood recipe; constant-SE comparison isolates heterogeneous SE input within Gaussian recipes.",
            "Fitted nuisance variance is not validated individual prediction uncertainty.",
            "One split/seed; no native-model superiority or protected evaluation."],
        input_hashes={str(p):sha(p) for p in (LABELS,SOURCE_SE,path)},protected_outcomes_read=False)
    emit(out/"recipe.json",recipe)
    return joined,path,configs


def likelihood_checks(out):
    import jax
    import jax.numpy as jnp
    import model_scalar as M
    predicted=jnp.array([.3,-.4],dtype=jnp.float32)
    observed=jnp.array([.1,-.1],dtype=jnp.float32)
    se=jnp.array([.05,.2],dtype=jnp.float32)
    parameter=jnp.asarray(-2.,jnp.float32)
    expected_variance=np.asarray(se,dtype=float)**2+np.logaddexp(0,-2.)+1e-8
    expected=.5*np.mean((np.asarray(predicted,dtype=float)-np.asarray(observed,dtype=float))**2/expected_variance+np.log(expected_variance))
    actual=M.effect_loss(predicted,observed,loss="gaussian_effect",se=se,log_residual_variance=parameter)
    np.testing.assert_allclose(float(actual),expected,rtol=1e-6,atol=1e-7)
    loss=lambda pred,p:M.effect_loss(pred,observed,loss="gaussian_effect",se=se,log_residual_variance=p)
    gp,gv=jax.grad(loss,argnums=(0,1))(predicted,parameter)
    np.testing.assert_allclose(np.asarray(gp),(np.asarray(predicted)-np.asarray(observed))/expected_variance/2,rtol=1e-6,atol=1e-7)
    assert np.isfinite(gv) and float(gv) != 0
    emit(out/"likelihood_checks.json",dict(Gaussian_loss_matches_independent_formula=True,
        prediction_gradient_matches_SE_weighted_formula=True,finite_nonzero_nuisance_gradient=True,
        constant_SE_ignores_closed_fold_annotations=True,precision="float32"))


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("GPU compute-node execution required")
    args.out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    labels,path,configs=prepare(args.out)
    likelihood_checks(args.out)
    trainer=Path(__file__).with_name("model_train.py")
    emit(args.out/"source_code.json",{str(p):sha(p) for p in
        (Path(__file__),trainer,Path(__file__).with_name("model_scalar.py"),Path(__file__).with_name("catalog_currin.py"))})
    records=[]
    def run(config,steps,limit,name,max_hours):
        command=[sys.executable,str(trainer),"--config",config["path"],"--labels",str(path),
            "--se-column",config["se_column"],"--out",str(args.out/name),"--steps",str(steps),
            "--microbatch","4","--validation-limit",str(limit),"--max-hours",str(max_hours)]
        tick=time.monotonic()
        result=subprocess.run(command,check=False)
        receipt=args.out/name/"feasibility.json"
        completed=json.loads(receipt.read_text()) if receipt.exists() else {}
        record={**config,"output_name":name,"exit_code":result.returncode,"seconds":time.monotonic()-tick,
            "steps_requested":steps,"steps_completed":completed.get("completed_steps",0),
            "status":"complete" if result.returncode == 0 and completed.get("completed_steps") == steps else "failed_or_incomplete"}
        if record["status"] == "complete":
            try:
                expected_n=limit or 6834
                assert completed["requested_steps"] == completed["completed_steps"] == steps
                assert completed["microbatch"] == 4 and completed["validation_examples"] == expected_n
                for flag in ("frozen_parameters_unchanged","running_statistics_unchanged","save_reload_agreement"):
                    assert completed[flag] is True
                split=json.loads((args.out/name/"split.json").read_text())
                assert split["held_fold"] == 0 and split["validation_fold"] == 1
                assert split["training_n"] == 18914 and split["validation_n"] == expected_n
                assert split["held_outcomes_evaluated"] is False and split["training_native_slope"] == 0
                frame=pd.read_csv(args.out/name/"validation_predictions.tsv",sep="\t")
                assert len(frame) == expected_n and frame.variant_id.is_unique and frame.validation_fold.eq(1).all()
                assert np.isfinite(frame[["predicted_beta","observed_beta"]].to_numpy(dtype=float)).all()
                valid=labels.loc[labels.heldout_fold.eq(1)]
                chosen=valid.iloc[np.linspace(0,len(valid)-1,expected_n,dtype=int)].set_index("lead_variant_id")
                assert set(frame.variant_id) == set(chosen.index)
                np.testing.assert_array_equal(frame.observed_beta.to_numpy(np.float32),chosen.loc[frame.variant_id,"beta_alt"].to_numpy(np.float32))
                record["complete_validation_and_frozen_state_checked"]=True
            except (AssertionError,KeyError,ValueError) as error:
                record["status"]="failed_validation_guard"
                record["validation_error"]=repr(error)
        weights=args.out/name/"weights.npz"
        if weights.exists():
            with np.load(weights,allow_pickle=False) as archive:
                parameter=float(archive["log_residual_variance"])
            record.update(initial_log_residual_variance=-2.,final_log_residual_variance=parameter,
                final_residual_variance=float(np.logaddexp(0,parameter)))
            if not np.isfinite(parameter):
                raise ValueError("Nonfinite learned training residual variance")
        records.append(record)
        emit(args.out/"runs.json",records)
        return record,completed
    probe,feasibility=run(configs[0],100,32,"probe100",min(1.,args.max_hours))
    if probe["status"] != "complete":
        raise RuntimeError("Gaussian probe failed; full comparison not admitted")
    overhead=max(0,probe["seconds"]-feasibility["training_wall_seconds_including_batching_and_progress_writes"]-feasibility["validation_wall_seconds"])
    train_cost=feasibility["iteration_wall_seconds_median"]*5000
    eval_cost=feasibility["validation_first_batch_seconds"]+feasibility["validation_warm_seconds_per_example"]*(6834-4)
    per_arm=2*(train_cost+eval_cost+overhead)
    remaining=args.max_hours*3600-(time.monotonic()-started)
    projection=dict(probe_seconds=probe["seconds"],arms=20,training_seconds_per_arm=train_cost,
        evaluation_seconds_per_arm=eval_cost,overhead_per_arm=overhead,safety_multiplier=2,
        projected_20_arm_seconds=20*per_arm,remaining_seconds=remaining,admitted=20*per_arm < remaining-600)
    emit(args.out/"projection.json",projection)
    if not projection["admitted"]:
        emit(args.out/"complete.json",dict(status="resource_deferred_after_probe",projection=projection))
        return
    for config in configs:
        remaining=args.max_hours*3600-(time.monotonic()-started)
        if remaining < per_arm+600:
            records.append({**config,"status":"resource_deferred_before_fit"})
            emit(args.out/"runs.json",records)
            continue
        run(config,5000,0,config["name"],min(remaining-120,per_arm+900)/3600)
    # Descriptive point errors only here; paired comparison has a separate reviewer.
    points=[]
    expected=labels.loc[labels.heldout_fold.eq(1)].set_index("lead_variant_id")
    for record in records[1:]:
        if record["status"] != "complete":
            continue
        frame=pd.read_csv(args.out/record["output_name"]/"validation_predictions.tsv",sep="\t").set_index("variant_id")
        if len(frame) != 6834 or not frame.index.is_unique or set(frame.index) != set(expected.index):
            raise ValueError("Full identical validation coverage required")
        np.testing.assert_allclose(frame.observed_beta,expected.loc[frame.index,"beta_alt"],rtol=1e-6,atol=1e-7)
        np.testing.assert_array_equal(frame.validation_fold,np.ones(6834,dtype=int))
        assert np.isfinite(frame[["predicted_beta","observed_beta"]].to_numpy(dtype=float)).all()
        mse=float(np.mean((frame.predicted_beta-frame.observed_beta)**2))
        points.append(dict(arm=record["name"],base_arm=record["base_arm"],SE_input=record["kind"],
            rows=6834,RMSE_beta_units=np.sqrt(mse),signed_spearman=spearmanr(frame.observed_beta,frame.predicted_beta).statistic,
            final_training_residual_variance=record["final_residual_variance"],seed=1103,steps=5000))
    pd.DataFrame(points).to_csv(args.out/"point_performance.tsv",sep="\t",index=False)
    emit(args.out/"complete.json",dict(status="completed_fixed_development_comparison" if len(points) == 20 else "partial_results_only",
        completed_full_arms=len(points),requested_full_arms=20,seconds=time.monotonic()-started,
        training_variance_only=True,protected_outcomes_read=False,selection_or_adoption=False))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--max-hours",type=float,default=3.8)
    main(parser.parse_args())
