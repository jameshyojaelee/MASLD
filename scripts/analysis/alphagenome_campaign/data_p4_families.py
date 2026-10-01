#!/usr/bin/env python3
"""Retain original P4 questions, minimum groups, and explicit missing dependencies.

F1-F3 use the existing registered methods, with historical and amended groups
distinguished. F4 retains original pooled, supported, and source-dependent groups.
F5 preserves the original gene universe; lack of model coverage cannot redefine it.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, spearmanr
from scipy.io import mmread

import data_p4 as p4

spec=importlib.util.spec_from_file_location("old_p4_tests",p4.ATLAS / "44_saturation_tests.py")
old=importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
CTX=p4.ROOT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
PROGRAMS=p4.ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
MIN_GROUP=200


def corrected(rows,family_name):
    available=[r for r in rows if r.get("p_nominal") is not None]
    if available:
        qs=false_discovery_control([r["p_nominal"] for r in available],method="bh")
        for row,q in zip(available,qs):
            row["BH_q"]=float(q);row["hypothesis_family"]=family_name
    return rows


def measured_contrast(df,value,group,match_promoter=True):
    valid=np.isfinite(pd.to_numeric(df[value],errors="coerce"))
    d=df.loc[valid].reset_index(drop=True)
    g=np.asarray(group)[valid]
    result={"n_group":int(g.sum()),"n_other":int((~g).sum()),"n_missing_feature":int((~valid).sum())}
    if min(result["n_group"],result["n_other"])<MIN_GROUP:
        return {**result,"status":"indeterminate_minimum_group","p_nominal":None}
    strata=old.stratum_labels(d,use_promoter=match_promoter)
    n_matched_group,n_matched_control=matching_support(strata,g)
    result.update({"matched_group_support":n_matched_group,"distinct_matched_control_support":n_matched_control})
    if min(n_matched_group,n_matched_control)<MIN_GROUP:
        return {**result,"status":"indeterminate_minimum_after_matching","p_nominal":None}
    draw=old.matched_draw_null(d[value].to_numpy(),strata,g,old.N_DRAWS,old.SEED)
    if draw["n_group"]<MIN_GROUP:
        return {**result,"status":"indeterminate_minimum_after_matching","matched_group":draw["n_group"],"p_nominal":None}
    return {**result,"status":"computed_development_comparison",**draw,"p_nominal":draw["exceedance_p"]}


def matching_support(strata,group):
    group=np.asarray(group,dtype=bool)
    shared=np.intersect1d(np.unique(strata[group]),np.unique(strata[~group]))
    usable=np.isin(strata,shared)
    return int((usable&group).sum()),int((usable&~group).sum())


def u2_signal_covariates(cohort,lineage,region_keys):
    """Donor logCPM covariates use all source peaks for each donor's library size."""
    path=CTX / "counts" / cohort / lineage
    peaks=pd.read_csv(str(path)+".peaks.tsv",sep="\t")
    donors=pd.read_csv(str(path)+".donors.tsv",sep="\t")
    matrix=mmread(str(path)+".mtx.gz").tocsr()
    if matrix.shape!=(len(donors),len(peaks)) or peaks.peak_coordinate.duplicated().any():
        raise ValueError("U2 donor/peak count axes disagree")
    library=np.asarray(matrix.sum(axis=1)).ravel()
    if not np.array_equal(library,donors.total_counted_fragments.to_numpy()):
        raise ValueError("full peak sums disagree with source counted fragments")
    eligible=donors.contrast_eligible.astype(str).str.upper().eq("TRUE").to_numpy() & (library>0)
    chosen=peaks.peak_coordinate.isin(region_keys).to_numpy()
    data=matrix[eligible][:,chosen].toarray().astype(float)
    logcpm=np.log2(1+data/library[eligible,None]*1_000_000)
    return pd.DataFrame({"region_key":peaks.loc[chosen,"peak_coordinate"].to_numpy(),"signal_mean":logcpm.mean(axis=0),"signal_sd":logcpm.std(axis=0,ddof=1),"n_donors":int(eligible.sum()),"source_count_columns":matrix.shape[1],"library_scope":"all_source_peaks_before_U2_subset"})


def u2_families(feats,tf,universe,out):
    import pysam
    fasta=Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")
    genome=pysam.FastaFile(str(fasta))
    da=pd.read_csv(CTX / "da/da_peak_results.tsv.gz",sep="\t",keep_default_na=False)
    u2=universe.loc[universe.universe.eq("U2_snatac")]
    results=[];compositions=[];covariates=[]
    for lineage in ("hepatocyte","stellate"):
        selected=u2.loc[u2.lineage.eq(lineage)].drop_duplicates("region_key")
        annotation=da.loc[da.lineage.eq(lineage),["peak_coordinate","evidence_state","promoter_genes"]]
        selected=selected.merge(annotation,left_on="region_key",right_on="peak_coordinate",how="left",validate="one_to_one")
        selected=selected.merge(feats.loc[feats.universe.eq("U2_snatac")],on="region_key",suffixes=("","_features"),validate="one_to_one")
        selected["gc"]=[(seq.count("G")+seq.count("C"))/len(seq) for seq in (genome.fetch(r.chrom,int(r.start0),int(r.end)).upper() for r in selected.itertuples())]
        selected["promoter"]=selected.promoter_genes.fillna("").apply(lambda v:int(v not in ("","NA","none","None","nan")))
        for cohort in ("GSE244832","GSE281367"):
            signals=u2_signal_covariates(cohort,lineage,set(selected.region_key))
            d=selected.merge(signals,on="region_key",validate="one_to_one").reset_index(drop=True)
            covariates.append(d[["region_key","gc","width","promoter","signal_mean","signal_sd","n_donors","source_count_columns","library_scope"]].assign(cohort=cohort,lineage=lineage))
            for label,states in (("original_nonindeterminate",("supported","source_dependent","discordant")),("supported",("supported",)),("source_dependent",("source_dependent",))):
                # Every named group is compared with source-indeterminate regions; another
                # nonindeterminate evidence class cannot silently enter the controls.
                eligible=d.evidence_state.isin(states+("indeterminate",))
                sub=d.loc[eligible].reset_index(drop=True)
                group=sub.evidence_state.isin(states).to_numpy()
                for trackgroup in ("liver_atac","liver_atac_adult_verified","liver_dnase_adult_verified","liver_h3k27ac_adult_verified"):
                    result=measured_contrast(sub,f"{trackgroup}_rel_mean",group)
                    results.append({"cohort":cohort,"lineage":lineage,"evidence_group":label,"track_group":trackgroup,**result})
                strata=old.stratum_labels(sub)
                if min(matching_support(strata,group))<MIN_GROUP:
                    compositions.append({"cohort":cohort,"lineage":lineage,"evidence_group":label,"status":"indeterminate_minimum_after_matching","p_nominal":None})
                    continue
                index={k:i for i,k in enumerate(sub.region_key)}
                relevant=tf.loc[tf.region_key.isin(index)&tf.group.isin(("chip_tf_liver","chip_tf_liver_adult_verified"))]
                for (trackgroup,track),rows in relevant.groupby(["group","track_name"]):
                    present=np.zeros(len(sub));present[[index[k] for k in rows.region_key]]=1
                    if present.sum()<50:
                        continue
                    cmh=old.cmh_test(present,group,strata)
                    draw=old.matched_draw_null(present,strata,group,1000,old.SEED,statistic="mean")
                    compositions.append({"cohort":cohort,"lineage":lineage,"evidence_group":label,"group":trackgroup,"track":track,"transcription_factor":rows.transcription_factor.iloc[0],"status":"computed_development_comparison","share_group":draw["observed_group"],"share_matched_null":draw["null_median"],"p_nominal":cmh["p"]})
    genome.close()
    pd.concat(covariates,ignore_index=True).to_csv(out / "F4_cohort_native_covariates.tsv.gz",sep="\t",index=False)
    corrected(results,"F4_all_cohort_lineage_evidence_track_contrasts")
    corrected(compositions,"F4_all_cohort_lineage_evidence_TF_contrasts")
    pd.DataFrame(results).to_csv(out / "F4_DA_sensitivity.tsv",sep="\t",index=False)
    pd.DataFrame(compositions).to_csv(out / "F4_DA_TF_composition.tsv",sep="\t",index=False)
    return {"contrasts":len(results),"minimum_group_indeterminate":sum(r["status"].startswith("indeterminate") for r in results),"covariates":"cohort_native_donor_logCPM_mean_sd;all_peak_library_size;source_eligible_donors","matching_amendment":"cohort_native_signal_covariates_chosen_before_full_F4_outcomes;historical_pooled_evidence_group_retained"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    check=json.loads((args.features / "reduction_summary.json").read_text())
    if check["status"]!="full_reduction_complete" or check["regions"]!=103390 or check["failures"]:
        raise SystemExit("full successful103390-region reduction required")
    args.out.mkdir(parents=True,exist_ok=False)
    feats=pd.read_csv(args.features / "saturation_region_features.tsv.gz",sep="\t",low_memory=False)
    cov=pd.read_csv(old.COVARIATES,sep="\t",low_memory=False)
    u1=feats.loc[feats.universe.eq("U1_h3k27ac")].merge(cov,on="region_key",suffixes=("","_cov"),how="left",validate="one_to_one")
    needed=["gc","width","promoter","signal_mean","signal_sd","reliable_shipped_form"]
    if any(c not in u1 for c in needed):
        raise ValueError("original U1 matching covariates unavailable")
    missing=u1[needed].isna().any(axis=1)
    d=u1.loc[~missing].reset_index(drop=True)
    group=d.reliable_shipped_form.astype(str).str.lower().eq("true").to_numpy()
    historical=list(old.GROUPS_PRIMARY)
    amended=[f"{g}_{stage}" for g in historical for stage in p4.STAGES]
    f1=[];f2=[]
    for name in historical+amended:
        for scale in p4.OLD.SCALES:
            f1.append({"group":name,"scale":scale,"interpretation":"historical" if name in historical else "life_stage_amendment",**measured_contrast(d,f"{name}_{scale}_mean",group)})
            f2.append({"group":name,"scale":scale,"interpretation":"historical" if name in historical else "life_stage_amendment",**measured_contrast(d,f"{name}_{scale}_share_best50",d.promoter.eq(1).to_numpy(),False)})
    corrected(f1,"F1_all_historical_and_life_stage_by_scale_development")
    corrected(f2,"F2_all_historical_and_life_stage_by_scale_development")
    pd.DataFrame(f1).to_csv(args.out / "F1_predictability.tsv",sep="\t",index=False)
    pd.DataFrame(f2).to_csv(args.out / "F2_localisation.tsv",sep="\t",index=False)
    # The second registered F1 form is descriptive within-stratum skill correlation.
    skill="skill_shipped_mean" if "skill_shipped_mean" in d else None
    second=[]
    for name in historical+amended:
        for scale in p4.OLD.SCALES:
            second.append({"group":name,"scale":scale,"status":"computed" if skill else "missing_source_skill_column",**(old.spearman_within_strata(d[skill],d[f"{name}_{scale}_mean"],old.stratum_labels(d)) if skill else {})})
    pd.DataFrame(second).to_csv(args.out / "F1_skill_correlation.tsv",sep="\t",index=False)
    signal=[]
    for name in historical+amended:
        valid=np.isfinite(d.signal_mean)&np.isfinite(d[f"{name}_rel_mean"])
        n=int(valid.sum())
        signal.append({"group":name,"n_regions":n,"rho":float(spearmanr(d.loc[valid,"signal_mean"],d.loc[valid,f"{name}_rel_mean"]).statistic) if n>=3 else None,"status":"descriptive_P5_signal_abundance"})
    pd.DataFrame(signal).to_csv(args.out / "P5_signal_mean.tsv",sep="\t",index=False)
    tf=pd.read_csv(args.features / "saturation_region_tf_top5.tsv.gz",sep="\t")
    all_tf=tf
    tf=tf.loc[tf.region_key.isin(d.region_key)]
    index={k:i for i,k in enumerate(d.region_key)}
    strata=old.stratum_labels(d)
    f3=[]
    for (name,track),sub in tf.groupby(["group","track_name"]):
        present=np.zeros(len(d))
        present[[index[k] for k in sub.region_key]]=1
        if present.sum()<50:
            f3.append({"group":name,"track":track,"status":"original_under50_top5_occurrences","p_nominal":None})
            continue
        if min(matching_support(strata,group))<MIN_GROUP:
            f3.append({"group":name,"track":track,"status":"indeterminate_minimum_group","p_nominal":None})
            continue
        draw=old.matched_draw_null(present,strata,group,1000,old.SEED,statistic="mean")
        cmh=old.cmh_test(present,group,strata)
        f3.append({"group":name,"track":track,"transcription_factor":sub.transcription_factor.iloc[0],"status":"computed_development_comparison","n_regions_top5":int(present.sum()),"share_predictable":draw["observed_group"],"share_matched_null":draw["null_median"],"matched_null_lo":draw["null_lo"],"matched_null_hi":draw["null_hi"],"matched_draws":1000,"p_nominal":cmh["p"]})
    corrected(f3,"F3_all_eligible_track_groups_original_occurrence_rule")
    pd.DataFrame(f3).to_csv(args.out / "F3_TF_composition.tsv",sep="\t",index=False)
    # Original U2 labels are source-defined; no new genetic group is manufactured.
    universe=pd.read_csv(p4.SAT / "tables/saturation_universe.tsv.gz",sep="\t",keep_default_na=False)
    u2=universe.loc[universe.universe.eq("U2_snatac")].copy()
    f4=u2_families(feats,all_tf,universe,args.out)
    da=pd.read_csv(CTX / "da/da_peak_results.tsv.gz",sep="\t",keep_default_na=False)
    manifest=pd.read_csv(CTX / "consensus_peak_manifest.tsv",sep="\t",dtype=str)
    manifest=manifest.loc[~manifest.blacklist_overlap.str.upper().eq("TRUE")]
    eligible_peak_ids=set(manifest.lineage+"|"+manifest.chrom+":"+manifest.start0+"-"+manifest.end)
    da=da.loc[(da.lineage+"|"+da.peak_coordinate).isin(eligible_peak_ids)]
    links={}
    for row in da.itertuples():
        for gene in str(row.promoter_genes).replace(";",",").split(","):
            if gene and gene not in ("NA","nan","none","None"):
                links.setdefault(gene,set()).add(row.peak_coordinate)
    available=set(feats.loc[feats.universe.eq("U2_snatac"),"region_key"])
    missing_genes={g:len(regions-available) for g,regions in links.items() if regions-available}
    pd.DataFrame([{"canonical_gene":gene,"source_promoter_regions":len(regions),"scored_regions":len(regions&available),"missing_saturation_regions":len(regions-available)} for gene,regions in sorted(links.items())]).to_csv(args.out / "F5_original_background_coverage.tsv",sep="\t",index=False)
    pd.DataFrame([{"canonical_gene":gene,"region_key":region,"required_inference":"same_hosted_P4_scorers_and_archive_terms;not_acquired_by_this_continuation"} for gene,regions in sorted(links.items()) for region in sorted(regions-available)]).to_csv(args.out / "F5_original_background_missing_inference.tsv.gz",sep="\t",index=False)
    programs=pd.read_csv(PROGRAMS,sep="\t")
    f5=[]
    for uid,part in programs.groupby("program_uid"):
        genes=set(part.canonical_gene.dropna())
        regions=set().union(*(links.get(g,set()) for g in genes))
        n=len(regions & available)
        reason="original_random_gene_set_background_contains_unscored_promoter_regions" if missing_genes else "original_size_matched_gene_null_not_yet_executed"
        if n<MIN_GROUP:
            reason="minimum200_profile_regions_not_met;"+reason
        f5.append({"program_uid":uid,"n_source_genes":len(genes),"n_promoter_regions":len(regions),"n_scored_regions":n,"n_missing_regions":len(regions-available),"status":"indeterminate" if n<MIN_GROUP else "not_run_original_null","reason":reason})
    pd.DataFrame(f5).to_csv(args.out / "F5_program_dispositions.tsv",sep="\t",index=False)
    summary={"full_reduction":True,"U1_regions":len(u1),"U1_missing_matching_covariates":int(missing.sum()),"original_minimum_regions_each_side":MIN_GROUP,"F1_primary_comparison":"historical_liver_h3k27ac_relative_scale;amended_adult_and_other_channels_development","F1_F2_BH":"complete reported development comparison families; native matched-draw p retained","F3_BH":"all eligible TF track groups passing original50-occurrence rule","F4":f4,"F5_complete":False,"F5_original_eligible_source_genes":len(links),"F5_background_genes_with_unscored_promoter_regions":len(missing_genes),"F5_background_not_restricted_to_archive":True,"U2_original_genetic_overlap":0,"seed":old.SEED,"matched_draws":old.N_DRAWS,"U1_coordinate_limit":check["U1_coordinate_limit"],"scientific_F1_F5_complete":False}
    (args.out / "family_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":
    main()
