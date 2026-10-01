#!/usr/bin/env python3
"""Nested strong RNA and nonlinear shared-context comparisons on fixed regions."""
import argparse,gzip,hashlib,importlib.metadata,json,os,platform,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import model_donor_common as C
import model_donor_nn as N
import review_donor as D

CHECKPOINTS=(100,300,1000)
OBJECTIVES=("profile","residual","pairwise")
ARCHITECTURES=("film","additive_nn")


def selection_indices(fold,reg,held,mode):
    tr=np.flatnonzero(fold != held);te=np.flatnonzero(fold == held)
    inner=(held+1)%5;it=np.flatnonzero((fold != held)&(fold != inner));iv=np.flatnonzero(fold == inner)
    if mode == "trained_regions":rt=np.arange(len(reg));rv=rt;irt=rt;irv=rt
    else:
        rt=np.flatnonzero(reg.region_role == "train");rv=np.flatnonzero(reg.region_role == "held")
        irt=np.flatnonzero(reg.inner_region_role == "inner_train");irv=np.flatnonzero(reg.inner_region_role == "inner_valid")
    return tr,te,it,iv,rt,rv,irt,irv


def baseline_pair(y,seq,it,iv,irt,irv,tr,rt,mode):
    si=C.sequence_projection(seq,irt);so=C.sequence_projection(seq,rt)
    inner_mean=y[np.ix_(it,irt)].mean(0);outer_mean=y[np.ix_(tr,rt)].mean(0)
    if mode == "trained_regions":return inner_mean,outer_mean,si,so,None
    candidates={float(c):C.sequence_mean(si,inner_mean,c) for c in C.RIDGE_C}
    best=min(candidates,key=lambda c:np.mean((candidates[c][irv][None]-y[np.ix_(iv,irv)])**2))
    return candidates[best],C.sequence_mean(so,outer_mean,best),si,so,best


def baseline_models(y,ri,ro,si,so,bi,bo,it,iv,irt,irv,tr,rt,mode):
    predictions={"training_baseline":np.broadcast_to(bo,(len(y),len(bo))).copy()};chosen=[]
    inner_y=y[np.ix_(it,irt)];outer_y=y[np.ix_(tr,rt)]
    predictions["sequence_only"]=np.broadcast_to(C.sequence_mean(so,outer_y.mean(0),.01),(len(y),len(bo))).copy()
    for kind in ("rna_global","additive_linear","product_kernel"):
        additive=kind != "product_kernel"
        inner_base=np.full_like(bi,inner_y.mean()) if kind == "rna_global" else bi
        outer_base=np.full_like(bo,outer_y.mean()) if kind == "rna_global" else bo
        candidates={float(c):C.product_prediction(ri,si,inner_y,c,inner_base,additive) for c in C.RIDGE_C}
        best=min(candidates,key=lambda c:np.mean((candidates[c][np.ix_(iv,irv)]-y[np.ix_(iv,irv)])**2))
        predictions[kind]=C.product_prediction(ro,so,outer_y,best,outer_base,additive)
        chosen.append({"model":kind,"c":best})
    if mode == "trained_regions":
        inner=C.rna_candidates(ri,inner_y);outer=C.rna_candidates(ro,outer_y)
        keys=list(inner);errors=np.stack([np.mean((inner[k][iv]-y[iv])**2,axis=0) for k in keys])
        for family in ("full_rna_ridge","pcr","pls_svd_ridge","selected_rna_linear"):
            allowed=np.arange(len(keys)) if family == "selected_rna_linear" else np.array([i for i,k in enumerate(keys) if k.startswith(family+":")])
            best=allowed[np.argmin(errors[allowed],axis=0)]
            predictions[family]=np.column_stack([outer[keys[int(k)]][:,j] for j,k in enumerate(best)])
            for j,k in enumerate(best):chosen.append({"model":family,"region_index_in_panel":j,"recipe":keys[int(k)]})
    return predictions,chosen


def bh(values):
    p=np.asarray(values,float);order=np.argsort(p);q=np.empty(len(p))
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q


def contrast_pairs(keys,mode):
    pairs=[]
    for key in keys:
        arm,model=key.split("|",1)
        anchor=arm+"|"+("selected_rna_linear" if mode == "trained_regions" else "product_kernel")
        if key != anchor:pairs.append((key,anchor,"architecture_MSE"))
        if arm == "annotated_all_rna":pairs.append((key,"all_rna|"+model,"annotation_coverage_MSE"))
    for level in sorted({k.split("|")[0] for k in keys if k.startswith("random_")}):
        exclusion=level.removeprefix("random_")
        for key in keys:
            arm,model=key.split("|",1)
            if arm == exclusion:
                pairs.extend([(key,"annotated_all_rna|"+model,"context_exclusion_MSE"),(key,level+"|"+model,"context_exclusion_MSE")])
    return pairs


def bootstrap_errors(ea,eb,fold,blocks,draws,rng):
    by_fold=[np.flatnonzero(fold == f) for f in range(5)]
    unique=np.unique(blocks);by_block={b:np.flatnonzero(blocks == b) for b in unique}
    sampled=[];sampled_diff=[]
    for _ in range(draws):
        di=np.concatenate([rng.choice(d,len(d),replace=True) for d in by_fold])
        ri=np.concatenate([by_block[b] for b in rng.choice(unique,len(unique),replace=True)])
        a=ea[np.ix_(di,ri)];b=eb[np.ix_(di,ri)]
        sampled.append(np.mean(a*a-b*b));sampled_diff.append((D.donor_difference_loss(a,fold[di])-D.donor_difference_loss(b,fold[di])).mean())
    return np.asarray(sampled),np.asarray(sampled_diff)


def summarize(predictions,y,fold,reg,mode,out,draws):
    regions=np.arange(len(reg)) if mode == "trained_regions" else np.flatnonzero(reg.region_role == "held")
    target=y[:,regions];blocks=reg.input_component.to_numpy()[regions];metrics=[];comparisons=[]
    valid={k:v[:,regions] for k,v in predictions.items() if np.isfinite(v[:,regions]).all()}
    for key,pred in valid.items():
        corr=[spearmanr(pred[:,j],target[:,j]).statistic for j in range(len(regions)) if np.ptp(pred[:,j])>1e-10]
        metrics.append({"evaluation":mode,"model":key,"donors":99,"regions":len(regions),"genomic_components":len(set(blocks)),
            "RMSE_log2CPM":float(np.sqrt(np.mean((pred-target)**2))),
            "donor_difference_MSE":float(D.donor_difference_loss(pred-target,fold).mean()),
            "mean_region_donor_spearman":float(np.mean(corr)) if corr else None})
    rng=np.random.default_rng(C.SEED);pairs=contrast_pairs(predictions,mode)
    for key,anchor,family in pairs:
        if key not in valid or anchor not in valid:
            comparisons.append({"evaluation":mode,"family":family,"model":key,"comparator":anchor,"status":"incomplete_not_tested","p_family_accounting":1.});continue
        ea=valid[anchor]-target;eb=valid[key]-target
        point=float(np.mean(ea**2-eb**2));dd=float((D.donor_difference_loss(ea,fold)-D.donor_difference_loss(eb,fold)).mean())
        sampled,sampled_diff=bootstrap_errors(ea,eb,fold,blocks,draws,rng)
        lo,hi=np.quantile(sampled,(.025,.975));dl,dh=np.quantile(sampled_diff,(.025,.975))
        p=(1+np.sum(np.abs(np.asarray(sampled)-point)>=abs(point)))/(draws+1)
        comparisons.append({"evaluation":mode,"family":family,"model":key,"comparator":anchor,"status":"development_fixed_predictions",
            "MSE_improvement":point,"low95":lo,"high95":hi,"donor_difference_MSE_improvement":dd,"difference_low95":dl,"difference_high95":dh,
            "p_family_accounting":p,"uncertainty":"donors_within_fold_and_genomic_components; fitting_selection_uncertainty_not_included"})
    for family in sorted({r["family"] for r in comparisons}):
        rows=[r for r in comparisons if r["family"] == family];q=bh([r["p_family_accounting"] for r in rows])
        for row,value in zip(rows,q):row["BH_q_complete_mode_family"]=value;row["planned_family_n"]=len(rows)
    pd.DataFrame(metrics).to_csv(out/f"{mode}_metrics.tsv",sep="\t",index=False)
    pd.DataFrame(comparisons).to_csv(out/f"{mode}_paired_uncertainty.tsv",sep="\t",index=False)


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):raise ValueError("Compute allocation required")
    args.out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    with np.load(args.fixture/"paired.npz",allow_pickle=False) as f:
        x=f["rna_log2cpm_full_library"];y=f["target_log2cpm_full_library"];fold=f["donor_fold"]
        genes=f["gene_ids"];chrom=f["gene_chrom"];gs=f["gene_start"];ge=f["gene_end"];seq=f["dinucleotide_features"]
    reg=pd.read_csv(args.fixture/"regions.tsv",sep="\t")
    if args.sequence_features:
        with np.load(args.sequence_features,allow_pickle=False) as f:
            if not np.array_equal(f["region_key"],reg.region_key):raise ValueError("Sequence region identities differ")
            seq=f["mean_strands"]
    if x.shape != (99,42163) or y.shape != (99,args.expect_regions) or len(reg) != args.expect_regions:
        raise ValueError(f"Paired source population differs: expected (99,{args.expect_regions})")
    if not np.isfinite(x).all() or not np.isfinite(y).all() or not np.isfinite(seq).all():raise ValueError("Nonfinite input")
    assignments=C.groups(reg,args.exclusion_levels);mask_meta=[];admitted_masks=[];removed_masks=[];selection=[];fit_records=[]
    projection_cache={};all_predictions={};timing_receipt=None
    planned_fits=len(assignments)*2*5*len(OBJECTIVES)*len(ARCHITECTURES)*2
    baseline_names=["training_baseline","sequence_only","rna_global","additive_linear","product_kernel"]
    neural_names=[a+"_"+o for a in ARCHITECTURES for o in OBJECTIVES]
    planned_keys={mode:[arm+"|"+name for arm in sorted({a[0] for a in assignments}) for name in
        baseline_names+neural_names+(["full_rna_ridge","pcr","pls_svd_ridge","selected_rna_linear"] if mode == "trained_regions" else [])]
        for mode in ("trained_regions","held_regions")}
    planned_contrasts=sum(len(contrast_pairs(keys,mode)) for mode,keys in planned_keys.items())
    recipe={"status":"development_recipe_fixed_before_new_fits","donors":99,"regions":len(reg),"sequence_features":str(args.sequence_features or "counted_interval_dinucleotides"),
        "RNA":"full42163_source_gene_features; real_training_PCA32_for_NN; fullRNA_kernel_and_PCR_PLS_SVD_baselines",
        "objectives":OBJECTIVES,"architectures":ARCHITECTURES,"checkpoints":CHECKPOINTS,"learning_rate":1e-3,"weight_decay":1e-3,"gradient_clip":1.,
        "full_library_normalization_before_exclusion":True,"excluded_genes_aggregate_library_denominator_retained":True,
        "exclusion_levels":args.exclusion_levels,"mask_strategy":"one_shared_model_per_target_chromosome_or_interval; same_mask_for_all_training_regions; evaluate_target_group_only",
        "random_exclusions":"training_variable_annotated_gene_universe; exact_removed_count; no_outcomes",
        "unknown_coordinates":"retained_in_unrestricted_allRNA; excluded_from_annotated_allRNA_and_both_exclusion_random_controls",
        "unknown_coordinate_gene_count":int(np.sum(chrom == "unknown")),
        "exclusion_contrast_anchor":"annotated_all_rna; same_known_gene_universe_before_explicit_removal",
        "annotation_coverage_contrast":"annotated_all_rna_vs_unrestricted_all_rna; separate_exploratory_family",
        "inner_donor_fold":"next established outer fold","inner_region_split":"training_tile_and16kb_overlap_components; outerheldregions_closed_for_selection",
        "strong_RNA_comparator_crossed":"not_eligible_not_run","held_region_mean":"predicted_from_training_regions_only",
        "planned_maximum_NN_fits":planned_fits,"planned_paired_contrasts":planned_contrasts,"max_core_hours":args.max_core_hours,"seed":C.SEED,
        "environment":{"python":sys.version,"platform":platform.platform(),"SLURM_JOB_ID":os.environ.get("SLURM_JOB_ID"),
            "versions":{p:importlib.metadata.version(p) for p in ("numpy","scipy","pandas","jax","jaxlib","optax")}},
        "no_external_or_protected_outcomes":True,"counted_interval_residual_bp":1}
    (args.out/"recipe.json").write_text(json.dumps(recipe,indent=2)+"\n")
    def projection(donors,arm,group,targets,mode,held,stage):
        keep,removed=C.gene_mask(x,donors,chrom,gs,ge,reg,arm,group,targets)
        key=(tuple(donors),hashlib.sha256(keep.tobytes()).hexdigest())
        if key not in projection_cache:projection_cache[key]=C.kernel_projection(x,donors,keep)
        p=projection_cache[key];admitted=np.zeros(len(genes),bool);admitted[p["genes"]]=True
        mask_meta.append({"mode":mode,"held_fold":held,"stage":stage,"arm":arm,"group":group,
            "admitted_genes":int(admitted.sum()),"removed_explicit_genes":int(removed.sum()),"unknown_coordinate_genes":int(np.sum(chrom == "unknown")),"mask_index":len(admitted_masks)})
        mask_meta[-1]["actual_RNA_training_rank"]=p["scores"].shape[1]
        mask_meta[-1]["actual_NN_PCA_rank"]=min(32,p["scores"].shape[1])
        admitted_masks.append(admitted);removed_masks.append(removed);return p
    def throughput(measured):
        nonlocal timing_receipt
        if timing_receipt is not None:return
        overhead=time.monotonic()-started-measured["first100_seconds"]-measured["full_panel_inference_seconds"]
        probe_tick=time.monotonic()
        bootstrap_errors(y,np.zeros_like(y),fold,reg.input_component.to_numpy(),32,np.random.default_rng(C.SEED))
        bootstrap_per_draw=(time.monotonic()-probe_tick)/32
        metric_tick=time.monotonic()
        for j in range(y.shape[1]):spearmanr(y[:,j],y[:,j])
        metric_per_model=time.monotonic()-metric_tick
        components={"NN_training_host_inclusive_seconds":planned_fits*max(CHECKPOINTS)*measured["warm_iteration_p90_seconds"],
            "nested_RNA_sequence_projections_and_linear_grids_seconds":overhead*len(assignments)*10,
            "checkpoint_inference_seconds":planned_fits*len(CHECKPOINTS)*measured["full_panel_inference_seconds"],
            "objective_architecture_compilation_allowance_seconds":6*measured["first100_seconds"],
            "all_planned_paired_bootstraps_seconds":planned_contrasts*args.bootstrap*bootstrap_per_draw,
            "all_planned_metrics_seconds":sum(map(len,planned_keys.values()))*metric_per_model}
        projection=2*sum(components.values())
        timing_receipt={**measured,"elapsed_all_overhead_seconds":time.monotonic()-started,"projected_full_seconds_with_factor2":projection,
            "projection_components_before_factor2":components,"planned_mask_assignments":len(assignments),
            "planned_maximum_NN_fits":planned_fits,"planned_paired_contrasts":planned_contrasts,
            "bootstrap_probe_draws":32,"bootstrap_seconds_per_draw":bootstrap_per_draw,
            "projected_CPU_core_hours":projection*args.cpus/3600,"admitted":projection*args.cpus/3600 <= args.max_core_hours}
        (args.out/"throughput_100steps.json").write_text(json.dumps(timing_receipt,indent=2)+"\n")
        if not timing_receipt["admitted"]:raise RuntimeError("Measured donor-comparison projection exceeds reserved core hours")
    try:
        for mode in ("trained_regions","held_regions"):
            pred={};all_predictions[mode]=pred
            for held in range(5):
                tr,te,it,iv,rt,rv,irt,irv=selection_indices(fold,reg,held,mode)
                bi,bo,si,so,base_c=baseline_pair(y,seq,it,iv,irt,irv,tr,rt,mode)
                for arm,group,targets in assignments:
                    destination=np.intersect1d(targets,rv)
                    if not len(destination):continue
                    if (time.monotonic()-started)*args.cpus/3600 > args.max_core_hours-.1:raise RuntimeError("Reserved donor core-hour wall guard")
                    ri=projection(it,arm,group,targets,mode,held,"inner");ro=projection(tr,arm,group,targets,mode,held,"outer")
                    baseline,choices=baseline_models(y,ri,ro,si,so,bi,bo,it,iv,irt,irv,tr,rt,mode)
                    for name,values in baseline.items():pred.setdefault(arm+"|"+name,np.full_like(y,np.nan))[np.ix_(te,destination)]=values[np.ix_(te,destination)]
                    for choice in choices:selection.append({"mode":mode,"held_fold":held,"arm":arm,"group":group,**choice})
                    for architecture in ARCHITECTURES:
                        for objective in OBJECTIVES:
                            inner,cost=N.train(ri["nn"],si["nn"],y,it,irt,bi,objective,architecture == "film",CHECKPOINTS,throughput)
                            best=min(inner,key=lambda epoch:np.mean((inner[epoch][np.ix_(iv,irv)]-y[np.ix_(iv,irv)])**2))
                            outer,outer_cost=N.train(ro["nn"],so["nn"],y,tr,rt,bo,objective,architecture == "film",(best,))
                            name=arm+"|"+architecture+"_"+objective
                            pred.setdefault(name,np.full_like(y,np.nan))[np.ix_(te,destination)]=outer[best][np.ix_(te,destination)]
                            fit_records.append({"mode":mode,"held_fold":held,"arm":arm,"group":group,"model":name,
                                "selected_steps":best,"inner":cost,"outer":outer_cost,"sequence_mean_c":base_c})
                    print(json.dumps({"mode":mode,"fold":held,"arm":arm,"group":group,"elapsed_seconds":time.monotonic()-started}),flush=True)
                    (args.out/"completed_fits.json").write_text(json.dumps(fit_records,indent=2)+"\n")
            summarize(pred,y,fold,reg,mode,args.out,args.bootstrap)
        status="completed_nested_development"
    except Exception as exc:
        status="failed_or_resource_deferred";error=repr(exc)
        (args.out/"failure.json").write_text(json.dumps({"error":error,"seconds":time.monotonic()-started,"completed_NN_fits":len(fit_records)},indent=2)+"\n")
        raise
    finally:
        np.savez_compressed(args.out/"internal_predictions.npz",**{mode+"|"+k:v for mode,pred in all_predictions.items() for k,v in pred.items()})
        np.savez_compressed(args.out/"gene_mask_membership.npz",gene_ids=genes,gene_chrom=chrom,admitted=np.asarray(admitted_masks),explicit_removed=np.asarray(removed_masks))
        pd.DataFrame(mask_meta).to_csv(args.out/"gene_masks.tsv",sep="\t",index=False)
        pd.DataFrame(selection).to_csv(args.out/"inner_selection.tsv",sep="\t",index=False)
    (args.out/"complete.json").write_text(json.dumps({"status":status,"seconds":time.monotonic()-started,"donors":99,"regions":len(reg),
        "completed_NN_fits":len(fit_records),"projection_cache_entries":len(projection_cache),"biological_unit":"participant",
        "source_RNA_fractional_values":"continuous_not_integer_counts","no_evaluation_fitted_scaling":True,"no_H3_covariates":True,
        "uncertainty":"paired_donors_within_fold_and_input_overlap_components; no_retraining_uncertainty; exploratory_BH_per_mode",
        "trained_and_held_region_comparators_separate":True,"exclusion_survival_not_trans_regulation_or_causality":True},indent=2)+"\n")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--fixture",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    p.add_argument("--sequence-features",type=Path);p.add_argument("--exclusion-levels",nargs="+",default=["chromosome"])
    p.add_argument("--max-core-hours",type=float,default=80);p.add_argument("--cpus",type=int,default=4);p.add_argument("--bootstrap",type=int,default=1000)
    p.add_argument("--expect-regions",type=int,default=128,help="Region count the fixture must hold; 128 is the original pilot.")
    main(p.parse_args())
