#!/usr/bin/env python3
"""Target-aware endogenous specialist versus a jointly trained shared/private model.

Separate assay-native outputs are normalized using training-only RMS scales.
Targets are weighted so each variant contributes one unit per assay. Shared
variant features and private measured-target/construct features are trained
jointly, rather than adding a prior reporter prediction as an extra column.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge

ALPHAS=(0.0001,0.01,1.0)
RATIOS=(0.1,1.0,10.0)
SEED=1103


def delta(alleles):
    return 0.5*((alleles[:,1]-alleles[:,0])+(alleles[:,3]-alleles[:,2]))


def compress(raw,train,dimension):
    scale=raw[train].std(axis=0);scale[scale<1e-6]=1
    model=TruncatedSVD(n_components=min(dimension,int(train.sum())-1),random_state=SEED)
    model.fit(raw[train]/scale)
    return model.transform(raw/scale)


def variant_weights(ids):
    counts=pd.Series(ids).value_counts()
    return np.array([1.0/counts[x] for x in ids])


def mse(y,p,w):
    return float(np.average((y-p)**2,weights=w))


def fit_prediction(v,t,c,pi,yp,yr,wp,trainp,trainr,testp,testr,alpha,ratio=None):
    ep=np.concatenate([v[pi],t],axis=1)
    er=np.concatenate([v,c],axis=1)
    sp=max(float(np.sqrt(np.average(yp[trainp]**2,weights=wp[trainp]))),1e-8)
    sr=max(float(np.sqrt(np.mean(yr[trainr]**2))),1e-8)
    if alpha is None:
        return np.zeros(testp.sum()),np.zeros(testr.sum()),sp,sr
    if ratio is None:
        model=Ridge(alpha=alpha,fit_intercept=False)
        weights=wp[trainp]/wp[trainp].sum()
        model.fit(ep[trainp],yp[trainp]/sp,sample_weight=weights)
        return model.predict(ep[testp])*sp,np.full(testr.sum(),np.nan),sp,sr
    # Shared variant coefficients; each assay additionally has private variant
    # and target/construct coefficients. Private variant penalty = alpha*ratio.
    dv=v.shape[1];dt=t.shape[1];dc=c.shape[1]
    endo=np.concatenate([v[pi],v[pi]/np.sqrt(ratio),t,np.zeros((len(pi),dv+dc))],axis=1)
    reporter=np.concatenate([v,np.zeros((len(v),dv+dt)),v/np.sqrt(ratio),c],axis=1)
    design=np.concatenate([endo[trainp],reporter[trainr]],axis=0)
    outcomes=np.concatenate([yp[trainp]/sp,yr[trainr]/sr])
    weights=np.concatenate([wp[trainp]/wp[trainp].sum(),np.full(trainr.sum(),1/trainr.sum())])
    model=Ridge(alpha=alpha,fit_intercept=False)
    model.fit(design,outcomes,sample_weight=weights)
    return model.predict(endo[testp])*sp,model.predict(reporter[testr])*sr,sp,sr


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    spec={"primary_comparison":"learned_endogenous_specialist_vs_learned_joint_shared_private",
        "development_amendment":"after_original_both_selected_zero; retain_null_selection_and_expose_best_learned_candidates_from_same_grid; no_stratum_or_fold_changes",
        "operational_margin":"5_percent_training_zero_predictor_MSE_within_each_declared_stratum",
        "margin_interpretation":"operational_precision_target_not_a_biological_minimum",
        "strata":"fixed_prior_transfer_strata_spec; no outcome-based revisions",
        "assays":"normalized_accessibility_ALT_dosage_beta_and_LX2_control_log2_reporter_allele_contrast",
        "folds":"operative_source_long_range_blocks_all_targets_conditions_assays_together",
        "weighting":"equal_variants_within_assay; equal_total_training_loss_per_assay",
        "representations":"16kb_local_variant_and_measured_target_or_construct; training_SVD64_variant_32_private",
        "alphas":ALPHAS,"private_variant_penalty_ratios":RATIOS,"seed":SEED,
        "no_protected_outcomes":True,"source_QTL_biological_n":"not_inferred_from_variants_or_peak_counts"}
    (args.out/"recipe.json").write_text(json.dumps(spec,indent=2)+"\n")
    complete=json.loads((args.features/"complete.json").read_text())
    source_spec=json.loads(args.strata_spec.read_text())
    if source_spec["strata"]["center_distance_bp"] != [0,500,1000,5000,10000,100000,"larger"]:
        raise ValueError("Source distance strata differ from the fixed recipe")
    manifest=pd.read_csv(args.features/"reporter_manifest.tsv",sep="\t")
    genomic_identity=["contig","variant_pos1","genomic_ref","genomic_alt"]
    if manifest.element_id.duplicated().any() or manifest.duplicated(genomic_identity).any():
        raise ValueError("Element weighting requires exactly one source element per genomic variant")
    pairs=pd.read_csv(args.features/"target_labels.tsv.gz",sep="\t")
    with np.load(args.features/"reporters_16384.npz",allow_pickle=False) as a:
        if not np.array_equal(a["element_ids"],manifest.element_id):
            raise ValueError("Reporter feature identity differs")
        if list(a["allele_order"]) != ["REF","ALT","REF_RC","ALT_RC"]:
            raise ValueError("Reporter feature allele order differs")
        for name in ("pooled_centre_bin","pooled_construct"):
            if not np.isfinite(a[name][a["extracted"]]).all():
                raise ValueError("An extracted reporter representation is nonfinite")
        validr=a["extracted"].copy();vraw=delta(a["pooled_centre_bin"]);craw=delta(a["pooled_construct"])
    with np.load(args.features/"target_pools_16384.npz",allow_pickle=False) as a:
        if not np.array_equal(a["pair_index"],pairs.pair_index):
            raise ValueError("Target feature identity differs")
        if not np.isfinite(a["pooled_target"][a["extracted"]]).all():
            raise ValueError("An extracted target representation is nonfinite")
        validp=a["extracted"].copy();traw=delta(a["pooled_target"])
    index={e:i for i,e in enumerate(manifest.element_id)}
    pi=np.array([index[e] for e in pairs.element_id])
    validp &= validr[pi]
    pairs=pairs.loc[validp].reset_index(drop=True);traw=traw[validp];pi=pi[validp]
    foldr=manifest.outer_fold.to_numpy(dtype=int);foldp=foldr[pi]
    if not np.array_equal(foldp,pairs.outer_fold.to_numpy(dtype=int)):
        raise ValueError("An allele/target/assay differs in source fold")
    if manifest.groupby("long_range_block_id").outer_fold.nunique().max() != 1:
        raise ValueError("A source locus spans folds")
    if set(foldp) != set(range(5)):
        raise ValueError("No complete five-fold endogenous panel")
    source=pd.read_csv(args.reporters,sep="\t")
    source=source.loc[(source.cell_line == "LX2")&source.eligible.astype(str).str.lower().eq("true")]
    rep=source.pivot(index="element_id",columns="experimental_replicate",values="effect_control").reindex(manifest.element_id)
    if rep.shape[1] != 4:
        raise ValueError("Four source control replicate indices required")
    validr &= np.isfinite(rep.to_numpy()).all(axis=1)
    yr=rep.to_numpy().mean(axis=1);yp=pairs.effect.to_numpy(dtype=float)
    wp=variant_weights(pairs.element_id.to_numpy())
    # Invalid feature rows are never training examples; finite zero storage
    # placeholders permit projection of arrays while masks exclude those rows.
    vraw=np.nan_to_num(vraw);craw=np.nan_to_num(craw)
    preds={name:np.full(len(pairs),np.nan) for name in ("endogenous_specialist","shared_private","endogenous_selected_with_null","shared_private_selected_with_null")}
    reporter_pred=np.full(len(manifest),np.nan)
    selections=[]
    # Definitions are fixed independently of new transfer fit outcomes.
    strata={"all_fully_observed":np.ones(len(pairs),dtype=bool),
        "variant_inside_peak":pairs.physical_variant_peak_overlap.to_numpy(dtype=bool),
        "construct_overlaps_peak":pairs.physical_construct_peak_overlap_bp.to_numpy()>0}
    for value in ("0-500","501-1000","1001-5000","5001-10000","10001-100000",">100000"):
        strata["distance_"+value]=pairs.distance_stratum.astype(str).to_numpy() == value
    for value in ("RNA_over_DNA_le1","RNA_over_DNA_gt1","missing"):
        strata["activity_"+str(value)]=pairs.activity_stratum.to_numpy() == value
    stratum_weights={name:np.where(mask,1.0,0.0) for name,mask in strata.items()}
    for name,mask in strata.items():
        stratum_weights[name][mask]=variant_weights(pairs.loc[mask,"element_id"].to_numpy())
    margins={name:np.full(5,np.nan) for name in strata}
    for held in range(5):
        inner=(held+1)%5
        outerp=foldp != held;outerr=(foldr != held)&validr
        innerp=outerp&(foldp != inner);innerr=outerr&(foldr != inner)
        testp=foldp == held;testr=(foldr == held)&validr
        valp=foldp == inner;valr=(foldr == inner)&validr
        for name,mask in strata.items():
            train=outerp&mask
            if train.any():
                margins[name][held]=0.05*float(np.average(yp[train]**2,weights=stratum_weights[name][train]))
        vi=compress(vraw,innerr,64);ci=compress(craw,innerr,32);ti=compress(traw,innerp,32)
        choices=[]
        for alpha in (*ALPHAS,None):
            p,_,_,_=fit_prediction(vi,ti,ci,pi,yp,yr,wp,innerp,innerr,valp,valr,alpha)
            choices.append((mse(yp[valp],p,wp[valp]),alpha))
        _,selected_special_alpha=min(choices,key=lambda z:z[0])
        _,special_alpha=min((z for z in choices if z[1] is not None),key=lambda z:z[0])
        choices=[]
        for alpha,ratio in [(a,r) for a in ALPHAS for r in RATIOS]+[(None,1.0)]:
            p,_,_,_=fit_prediction(vi,ti,ci,pi,yp,yr,wp,innerp,innerr,valp,valr,alpha,ratio)
            choices.append((mse(yp[valp],p,wp[valp]),alpha,ratio))
        _,selected_joint_alpha,selected_joint_ratio=min(choices,key=lambda z:z[0])
        _,joint_alpha,joint_ratio=min((z for z in choices if z[1] is not None),key=lambda z:z[0])
        vo=compress(vraw,outerr,64);co=compress(craw,outerr,32);to=compress(traw,outerp,32)
        preds["endogenous_specialist"][testp],_,sp,sr=fit_prediction(vo,to,co,pi,yp,yr,wp,outerp,outerr,testp,testr,special_alpha)
        preds["shared_private"][testp],reporter_pred[testr],_,_=fit_prediction(vo,to,co,pi,yp,yr,wp,outerp,outerr,testp,testr,joint_alpha,joint_ratio)
        preds["endogenous_selected_with_null"][testp],_,_,_=fit_prediction(vo,to,co,pi,yp,yr,wp,outerp,outerr,testp,testr,selected_special_alpha)
        preds["shared_private_selected_with_null"][testp],_,_,_=fit_prediction(vo,to,co,pi,yp,yr,wp,outerp,outerr,testp,testr,selected_joint_alpha,selected_joint_ratio)
        selections.append({"held_fold":held,"inner_fold":inner,"specialist_alpha":special_alpha,
            "joint_alpha":joint_alpha,"private_penalty_ratio":joint_ratio,"training_caQTL_RMS":sp,
            "selected_with_null_specialist_alpha":selected_special_alpha,
            "selected_with_null_joint_alpha":selected_joint_alpha,"selected_with_null_joint_ratio":selected_joint_ratio,
            "training_reporter_RMS":sr,"training_pairs":int(outerp.sum()),"training_reporters":int(outerr.sum())})
    for prediction in preds.values():
        if not np.isfinite(prediction).all():
            raise ValueError("Incomplete held predictions")
    pairs["fold"]=foldp;pairs["variant_weight"]=wp
    for name,prediction in preds.items():pairs[name]=prediction
    pairs.to_csv(args.out/"endogenous_oof.tsv.gz",sep="\t",index=False)
    pd.DataFrame({"element_id":manifest.element_id,"fold":foldr,"measured_reporter_log2_effect":yr,
                  "predicted_shared_private_log2_effect":reporter_pred,"eligible":validr}).to_csv(args.out/"reporter_oof.tsv.gz",sep="\t",index=False)
    pd.DataFrame(selections).to_csv(args.out/"inner_selection.tsv",sep="\t",index=False)
    pd.DataFrame(margins).rename_axis("held_fold").to_csv(args.out/"training_only_operational_margins.tsv",sep="\t")
    results=[];selected_rows=[];rng=np.random.default_rng(20260915)
    block=pairs.borzoi_long_range_group_id.to_numpy()
    if pairs.groupby("borzoi_long_range_group_id").fold.nunique().max() != 1:
        raise ValueError("A source association block spans folds")
    error_a=(preds["endogenous_specialist"]-yp)**2;error_b=(preds["shared_private"]-yp)**2
    for name,mask in strata.items():
        if mask.any():
            selected_a=preds["endogenous_selected_with_null"][mask]
            selected_b=preds["shared_private_selected_with_null"][mask]
            selected_rows.append({"stratum":name,"pairs":int(mask.sum()),
                "specialist_selected_RMSE_beta":np.sqrt(mse(yp[mask],selected_a,stratum_weights[name][mask])),
                "shared_private_selected_RMSE_beta":np.sqrt(mse(yp[mask],selected_b,stratum_weights[name][mask])),
                "both_all_zero":bool(np.all(selected_a == 0)&np.all(selected_b == 0)),
                "interpretation":"conditional_identical_zero_predictions_have_degenerate_difference; not_an_information_limit_or_useful_transfer_exclusion"})
        support={f:len(set(block[mask&(foldp == f)])) for f in range(5)}
        row={"stratum":name,"pairs":int(mask.sum()),"variants":pairs.loc[mask,"element_id"].nunique(),
             "blocks":len(set(block[mask])),"minimum_blocks_per_fold":min(support.values())}
        if not mask.any() or min(support.values())<2:
            row["status"]="insufficient_fold_support_or_outside_observable_targets";row["p_for_family_accounting"]=1
            results.append(row);continue
        weights=stratum_weights[name]
        point=float(np.average(error_a[mask]-error_b[mask],weights=weights[mask]))
        margin=float(np.average(margins[name][foldp[mask]],weights=weights[mask]))
        unique=np.unique(block[mask]);by_block={b:np.flatnonzero(mask&(block == b)) for b in unique}
        by_fold={f:np.unique(block[mask&(foldp == f)]) for f in range(5)}
        draws=[]
        for _ in range(args.bootstrap):
            chosen=np.concatenate([rng.choice(v,len(v),replace=True) for v in by_fold.values()])
            ix=np.concatenate([by_block[b] for b in chosen])
            draws.append(np.average(error_a[ix]-error_b[ix],weights=weights[ix]))
        draws=np.array(draws);lo,hi=np.quantile(draws,(0.025,0.975))
        p=(1+np.sum(np.abs(draws-point)>=abs(point)))/(len(draws)+1)
        row.update({"MSE_improvement":point,"low95":lo,"high95":hi,"training_only_operational_margin":margin,
            "specialist_RMSE_beta":np.sqrt(np.average(error_a[mask],weights=weights[mask])),
            "shared_private_RMSE_beta":np.sqrt(np.average(error_b[mask],weights=weights[mask])),
            "p_nominal_centered_block_bootstrap":p,"p_for_family_accounting":p,
            "status":"learned_candidates_development_estimate","precision_disposition":
                "identical_fixed_predictions_do_not_establish_information_limit" if np.array_equal(preds["endogenous_specialist"][mask],preds["shared_private"][mask]) else
                "conditional_fixed_prediction_upper_bound_below_operational_margin; excludes_neither_other_models_nor_retraining_gain" if hi<margin else
                "operational_useful_gain_not_excluded"})
        results.append(row)
    from model_final import bh
    q=bh([r["p_for_family_accounting"] for r in results])
    for row,value in zip(results,q):row["BH_q_complete_fixed_stratum_family"]=value
    pd.DataFrame(results).to_csv(args.out/"performance_by_fixed_stratum.tsv",sep="\t",index=False)
    pd.DataFrame(selected_rows).to_csv(args.out/"selected_with_null_disposition.tsv",sep="\t",index=False)
    (args.out/"complete.json").write_text(json.dumps({"status":"nested_development_comparison",
        "seconds":time.monotonic()-started,"eligible_pairs":len(pairs),"eligible_reporters":int(validr.sum()),
        "replicate_control_measurement":"four_paired_LX2_experimental_replicates",
        "caQTL_unit":"source_participant_based_association; variant_count_not_biological_n",
        "source_strata_spec":str(args.strata_spec),"source_strata_sha256":hashlib.sha256(args.strata_spec.read_bytes()).hexdigest(),
        "feature_receipt":complete,"planned_stratum_contrasts":len(results),
        "estimable_stratum_contrasts":sum(r["status"] == "learned_candidates_development_estimate" for r in results),
        "uncertainty":"paired_source_locus_blocks; conditional_on_source_beta_estimates_and_four_replicate_reporter_means; source_effect_SE_retained_but_not_propagated; reporter_measurement_and_retraining_uncertainty_not_propagated",
        "bootstrap":args.bootstrap,"bootstrap_seed":20260915,"source_expansion":"no_new_protected_data",
        "target_link":"not_supplied; no_nearest_gene_substitution"},indent=2)+"\n")
    print(pd.DataFrame(results).to_string(index=False))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features",type=Path,required=True)
    parser.add_argument("--reporters",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--strata-spec",type=Path,required=True)
    parser.add_argument("--bootstrap",type=int,default=2000)
    main(parser.parse_args())
