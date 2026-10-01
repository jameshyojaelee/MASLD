#!/usr/bin/env python3
"""Matched development tables and paired locus uncertainty after the GPU chain.

Draws resample historical 1-Mb bins and, separately, whole chromosomes within
folds, preserving matched rows across all predictors. The bins are not verified
independent LD units. CIs concern held-prediction sampling, not retraining or
model-selection uncertainty. Complete exploratory MSE-contrast families use BH;
these are not confirmatory protected tests or five-seed nested finalists.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

from model_summarize import main as summarize_sweep
PROJ=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(PROJ/"scripts/analysis/alphagenome_program"))
import i1_common as C


def bh(p):
    p=np.asarray(p,dtype=float);order=np.argsort(p)
    q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q


def describe(frame,arms,anchors,name,out,resamples,coverage,block_column="block_1mb",save_predictions=True):
    frame=frame.copy()
    planned_contrasts=sum(a != b for a in arms for b in anchors)
    if frame.groupby(block_column).heldout_fold.nunique().max() != 1:
        raise ValueError(f"A locus block spans folds in {name}")
    incomplete=[a for a in arms if not np.isfinite(frame[a].to_numpy(dtype=float)).all()]
    (out/(name+"_incomplete_arms.json")).write_text(json.dumps(
        [{"arm":a,"finite_predictions":int(np.isfinite(frame[a]).sum()),"required_predictions":len(frame),
          "disposition":"incomplete_fit_not_entered_in_paired_population"} for a in incomplete],indent=2)+"\n")
    arms=[a for a in arms if a not in incomplete]
    anchors=[a for a in anchors if a not in incomplete]
    if not arms or not anchors:
        raise ValueError(f"No complete models or anchors in {name}")
    keep=np.isfinite(frame[["beta_alt",*arms]].to_numpy(dtype=float)).all(axis=1)
    frame=frame.loc[keep].copy()
    if frame.key.duplicated().any():
        raise ValueError(f"Duplicate identity in {name}")
    y=frame.beta_alt.to_numpy(dtype=float);fold=frame.heldout_fold.to_numpy(dtype=int)
    if len(set(fold)) != 1 and set(fold) != {0,1,2,3,4}:
        raise ValueError(f"Unsupported partial fold population {name}")
    blocks=frame[block_column].to_numpy(dtype=str)
    uncertainty_unit=("paired_whole_chromosomes_within_fold" if block_column == "uncertainty_chromosome"
                      else "paired_historical_1Mb_bins_within_fold_not_verified_LD_units")
    if not len(frame):
        raise ValueError(f"Empty matched population {name}")
    predicted={a:frame[a].to_numpy(dtype=float) for a in arms}
    if len(set(fold)) == 5:
        macro_draws,rho_draws=C.block_bootstrap_within_fold(y,predicted,fold,blocks,resamples=resamples,seed=20260915)
    else:
        macro_draws={a:np.full(resamples,np.nan) for a in arms}
        rho_draws=C.block_bootstrap_plain(y,predicted,blocks,resamples=resamples,seed=20260915)
    # The error bootstrap uses exactly the same fold/block resampling scheme.
    # Error sums avoid constructing repeated large prediction tables.
    unique=np.array(sorted(set(blocks)))
    lookup={b:i for i,b in enumerate(unique)}
    code=np.array([lookup[b] for b in blocks]);counts=np.bincount(code,minlength=len(unique))
    squared=np.stack([(predicted[a]-y)**2 for a in arms],axis=1)
    sums=np.zeros((len(unique),len(arms)),dtype=float);np.add.at(sums,code,squared)
    groups=[]
    for f in sorted(set(fold)):
        ids=np.unique(code[fold == f]);groups.append(ids)
    rng=np.random.default_rng(20260915);error_draws=np.empty((resamples,len(arms)))
    for draw in range(resamples):
        chosen=np.concatenate([rng.choice(ids,len(ids),replace=True) for ids in groups])
        error_draws[draw]=sums[chosen].sum(axis=0)/counts[chosen].sum()
    metrics=[]
    for i,arm in enumerate(arms):
        p=predicted[arm];rho=C.fast_spearman(y,p)
        lo,hi,_=C.interval_of(rho_draws[arm])
        mlo,mhi,_=C.interval_of(macro_draws[arm])
        rmse=np.sqrt(np.mean((p-y)**2));rl,rh=np.quantile(np.sqrt(error_draws[:,i]),(0.025,0.975))
        metrics.append({"population":name,"arm":arm,"n":len(y),"source_eligible_n":coverage,
            "coverage_fraction":len(y)/coverage,"blocks":len(unique),"folds":len(set(fold)),
            "pooled_signed_spearman":rho,"rho_low95":lo,"rho_high95":hi,
            "macro_spearman":C.macro_over_folds(y,p,fold) if len(set(fold)) == 5 else np.nan,
            "macro_low95":mlo,"macro_high95":mhi,"RMSE_beta_units":rmse,"RMSE_low95":rl,"RMSE_high95":rh,
            "MAE_beta_units":np.mean(np.abs(p-y)),"uncertainty_unit":uncertainty_unit})
    contrasts=[]
    for arm in arms:
        for anchor in anchors:
            if arm == anchor:
                continue
            i,j=arms.index(arm),arms.index(anchor)
            improvement=np.mean((predicted[anchor]-y)**2)-np.mean((predicted[arm]-y)**2)
            draws=error_draws[:,j]-error_draws[:,i]
            lo,hi=np.quantile(draws,(0.025,0.975))
            p=(1+np.sum(np.abs(draws-improvement)>=abs(improvement)))/(resamples+1)
            rd=rho_draws[arm]-rho_draws[anchor];rlo,rhi,_=C.interval_of(rd)
            contrasts.append({"family":name,"arm":arm,"comparator":anchor,"n":len(y),
                "uncertainty_unit":uncertainty_unit,"resampling_units":len(unique),
                "MSE_improvement":improvement,"MSE_improvement_low95":lo,"MSE_improvement_high95":hi,
                "rho_improvement":C.fast_spearman(y,predicted[arm])-C.fast_spearman(y,predicted[anchor]),
                "rho_improvement_low95":rlo,"rho_improvement_high95":rhi,
                "p_nominal_centered_block_bootstrap_MSE":p})
    if contrasts:
        correction_p=[r["p_nominal_centered_block_bootstrap_MSE"] for r in contrasts]
        correction_p += [1.0]*(planned_contrasts-len(contrasts))
        q=bh(correction_p)[:len(contrasts)]
        for row,value in zip(contrasts,q):
            row["BH_q_complete_exploratory_MSE_family"]=value
            row["planned_family_contrasts"]=planned_contrasts
            row["estimable_family_contrasts"]=len(contrasts)
            row["untestable_contrasts_in_correction"]="p1_accounting_only_not_observed_pvalues"
    pd.DataFrame(metrics).to_csv(out/(name+"_performance.tsv"),sep="\t",index=False)
    pd.DataFrame(contrasts).to_csv(out/(name+"_paired_contrasts.tsv"),sep="\t",index=False)
    if save_predictions:
        frame.to_csv(out/(name+"_matched_predictions.tsv.gz"),sep="\t",index=False)
    print(f"{name}: {len(frame)} rows, {len(unique)} blocks, {len(arms)} arms",flush=True)


def describe_both(frame,arms,anchors,name,out,resamples,coverage):
    frame=frame.copy()
    frame["uncertainty_chromosome"]=frame.key.str.split(":").str[0].str.removeprefix("chr")
    if not frame.uncertainty_chromosome.isin([str(i) for i in range(1,23)]).all():
        raise ValueError("Unexpected chromosome in Currin uncertainty comparison")
    describe(frame,arms,anchors,name,out,resamples,coverage)
    describe(frame,arms,anchors,name+"_chromosome_sensitivity",out,resamples,coverage,
             block_column="uncertainty_chromosome",save_predictions=False)


def calibrate_common(frame,raw_columns):
    if set(frame.heldout_fold) != {0,1,2,3,4}:
        raise ValueError("Common-population calibration requires the five established folds")
    receipts=[];calibrated=[]
    for raw in raw_columns:
        name=raw+"__training_calibrated";frame[name]=np.nan
        for f in sorted(frame.heldout_fold.unique()):
            tr=frame.heldout_fold != f;te=~tr
            x=frame.loc[tr,raw].to_numpy();y=frame.loc[tr,"beta_alt"].to_numpy()
            slope=float(np.dot(x,y)/max(np.dot(x,x),1e-20))
            frame.loc[te,name]=frame.loc[te,raw]*slope
            receipts.append({"raw":raw,"heldout_fold":int(f),"training_n":int(tr.sum()),"slope":slope,"intercept":0})
        calibrated.append(name)
    return calibrated,receipts


def verified_specialists(out):
    specialist=pd.read_csv(C.TIER_A4,sep="\t")
    specialist=specialist.loc[specialist.model_peak_overlap == 1].copy()
    labels=pd.read_csv(C.C2_LABELS,sep="\t")
    labels["key"]=[C.key_of(*r) for r in labels[["chr","pos_hg38","ref","alt"]].itertuples(index=False,name=None)]
    authority=labels[["key","peak_id","beta_alt","heldout_fold"]]
    check=specialist[["key","peak_id","beta_alt","heldout_fold"]].merge(
        authority,on="key",how="left",validate="one_to_one",suffixes=("_specialist","_C2"),indicator=True)
    check["molecular_target_matches"]=check.peak_id_specialist == check.peak_id_C2
    check["effect_matches"]=np.isclose(check.beta_alt_specialist,check.beta_alt_C2,rtol=0,atol=1e-12)
    check["fold_matches"]=check.heldout_fold_specialist == check.heldout_fold_C2
    common_target=(check._merge == "both")&check.molecular_target_matches
    okay=common_target&check.effect_matches&check.fold_matches
    check.loc[~common_target].to_csv(out/"specialist_other_target_or_population_exclusions.tsv",sep="\t",index=False)
    check.loc[common_target&~okay].to_csv(out/"specialist_source_disagreements.tsv",sep="\t",index=False)
    if (common_target&~okay).any():
        raise ValueError("Specialist effect/fold differs for the same C2 variant-target; see specialist_source_disagreements.tsv")
    (out/"specialist_source_identity.json").write_text(json.dumps({"model_peak_overlap_rows":len(check),
        "exact_variant_target_eligible_rows":int(okay.sum()),"other_target_or_population_exclusions":int((~common_target).sum()),
        "variant_target_effect_fold_agreement_on_common_population":True,"beta_absolute_roundoff_tolerance":1e-12,
        "eligibility":"exact_variant_and_measured_peak_identity; different_source_target_is_not_a_mislabeled_effect",
        "source":str(C.TIER_A4),"label_authority":str(C.C2_LABELS)},indent=2)+"\n")
    return specialist.loc[specialist.key.isin(check.loc[okay,"key"])].copy()


def load_adaptation(path, authority):
    # The sweep writes variant identity as its first (index) column. Historical
    # native intersections lost that index's name; never infer identity by row order.
    adaptation=pd.read_csv(path,sep="\t",index_col=0).rename_axis("variant_id").reset_index()
    adaptation["key"]=adaptation.variant_id.str.removeprefix("chr")
    adaptation=adaptation.merge(authority[["key","beta_alt","heldout_fold","block_1mb"]],
                                on="key",how="left",validate="one_to_one",indicator=True)
    if not (adaptation._merge == "both").all() or not (adaptation.heldout_fold == 1).all():
        raise ValueError("Adaptation identity missing from native population or outside validation fold 1")
    np.testing.assert_allclose(adaptation.observed_beta,adaptation.beta_alt,rtol=1e-7,atol=1e-7,
                               err_msg="Adaptation effect differs from native label authority")
    return adaptation.drop(columns="_merge")


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    native=args.model/"native"
    completion=json.loads((native/"completion.json").read_text())
    if completion["unresolved"]:
        raise ValueError("Native population incomplete")
    specialist=verified_specialists(args.out)
    frame=pd.read_csv(native/"matched_predictions.tsv.gz",sep="\t")
    summarize_sweep(SimpleNamespace(sweep=args.sweep,out=args.out/"sweep_reconstruction",native=native))
    adaptation=load_adaptation(args.out/"sweep_reconstruction/native_matched_predictions.tsv.gz",frame)
    native_arms=["local_atac_liver__train_calibrated","local_dnase_liver__train_calibrated"]
    old=[c for c in frame if c.endswith("__ensemble")]
    describe_both(frame,[*native_arms,*old],native_arms,"native_full_existing_leads",args.out,args.resamples,completion["eligible"])
    frozen=pd.read_csv(args.model/"frozen/comparisons/predictions.tsv.gz",sep="\t")
    fa=[c for c in frozen if c.startswith("frozen_")]
    frozen_training_population=len(frozen)
    frozen=frozen.merge(frame[["key","local_atac_liver","local_dnase_liver",*old]],on="key",validate="one_to_one")
    fixed_anchors,fixed_cal=calibrate_common(frozen,["local_atac_liver","local_dnase_liver"])
    pd.DataFrame(fixed_cal).to_csv(args.out/"fixed_representation_native_calibration.tsv",sep="\t",index=False)
    describe_both(frozen,[*fa,*fixed_anchors,*old],fixed_anchors,"fixed_representation_common_rows",args.out,args.resamples,completion["eligible"])
    aa=[c for c in adaptation if c.startswith(("adapter_","partial_","frozen_"))]+["allele_identity","zero_effect","native_atac_calibrated","native_dnase_calibrated"]
    describe_both(adaptation,aa,["native_atac_calibrated","native_dnase_calibrated"],"adaptation_single_split",args.out,args.resamples,int((frame.heldout_fold == 1).sum()))
    raw=["chrombpnet_adult_hep","chrombpnet_adult_hep_gse281367","borzoi_atac","borzoi_dnase"]
    specialist=specialist[["key","peak_id",*raw]].merge(frame[["key","peak_id","beta_alt","heldout_fold","block_1mb","local_atac_liver","local_dnase_liver"]],on=["key","peak_id"],validate="one_to_one")
    raw += ["local_atac_liver","local_dnase_liver"]
    specialist=specialist.loc[np.isfinite(specialist[raw]).all(axis=1)].copy()
    sa,cal=calibrate_common(specialist,raw)
    pd.DataFrame(cal).to_csv(args.out/"specialist_training_calibration.tsv",sep="\t",index=False)
    describe_both(specialist,sa,sa[-2:],"specialist_eligible_common_rows",args.out,args.resamples,completion["eligible"])
    (args.out/"analysis.json").write_text(json.dumps({"status":"development_not_protected_confirmation","bootstrap_draws":args.resamples,
        "seed":20260915,"MSE_multiple_testing":"BH_within_each_complete_named_exploratory_population_family",
        "p_value_method":"two_sided_centered_paired_block_bootstrap_with_plus_one",
        "macro_spearman_definition":"tanh(mean(arctanh(per_fold_spearman))); equal_weight_folds; established_Fisher_z_summary",
        "inferential_unit":"participant_based_source_QTL_estimates; paired_prediction_resampling_by_historical_1Mb_bins_and_separately_whole_chromosomes; per_variant_effective_donor_n_not_established_here",
        "chromosome_sensitivity":"amended_after_coordinate_audit_before_complete_native_results; identical_models_rows_and_point_estimates; separate_BH_family_per_population_and_unit; few_chromosomes_limit_precision; no_choice_of_unit_by_significance",
        "historical_bin_limit":"1Mb_bins_are_not_verified_independent_LD_units; overlapping_native_windows_cross_bin_boundaries",
        "uncertainty_excludes":"retraining_and_selection_uncertainty","finalists_five_seeds_nested":"not_run",
        "raw_native_error":"inapplicable_different_effect_units; only_train_calibrated_errors_compared",
        "specialist_eligibility":"historicalTierA4_model_peak_overlap1_exact_C2_variant_and_measured_target_complete_predictions; allrawscores_calibratedoncommon_trainingfoldrows; other_targets_excluded_explicitly",
        "legacy_DNA_baselines":"existing_five_seed_ensembles; adaptationusesone_seed_andreservesfold0; training_exposure_differs",
        "fixed_feature_training_population":frozen_training_population,
        "fixed_feature_native_comparable_population":len(frozen),
        "remaining_training_coverage_asymmetry":"fixed_heads_include_short_window_eligible_rows_without_1Mb_native_scores; nativecalibration_uses_available_common_rows",
        "genetic_membership_or_Resource_adoption_changed":False},indent=2)+"\n")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model",type=Path,required=True)
    parser.add_argument("--sweep",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--resamples",type=int,default=1000)
    main(parser.parse_args())
