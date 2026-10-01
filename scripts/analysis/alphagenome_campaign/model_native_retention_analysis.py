#!/usr/bin/env python3
"""Measured regional liver-ATAC retention with paired donor/region uncertainty."""
import json,time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SEED=20260915
ADAPTATION_NAMES=tuple(f"adapter_r{r}_last{b}_{p}" for r in (4,16) for b in (3,5) for p in ("variant","symmetric"))+("partial_qv_last3_lr3e-06","partial_qv_last3_lr1e-05")


def calibrate(x,y,training):
    """Separate affine calibration for each fixed predictor; training rows only."""
    xx=x[training];yy=y[training];xm=xx.mean(0);ym=yy.mean()
    centered=xx-xm;denominator=np.sum(centered**2,axis=0)
    slope=np.divide(np.sum(centered*(yy-ym)[:,None],axis=0),denominator,
        out=np.zeros(x.shape[1]),where=denominator>1e-20)
    intercept=ym-slope*xm
    return x*slope+intercept,slope,intercept,np.full(len(y),ym)


def correlations(x,y):
    return np.array([spearmanr(x[:,j],y).statistic if np.ptp(x[:,j])>1e-12 else np.nan for j in range(x.shape[1])])


def resample_regions(panel,indices,unit,rng,stratified):
    parts=[]
    strata=sorted(panel.iloc[indices].fold.unique()) if stratified else [None]
    for fold in strata:
        subset=indices if fold is None else indices[panel.iloc[indices].fold.to_numpy()==fold]
        block=panel.iloc[subset][unit].to_numpy();unique=np.unique(block)
        groups={b:subset[block==b] for b in unique}
        parts.extend(groups[b] for b in rng.choice(unique,len(unique),replace=True))
    return np.concatenate(parts)


def bh(values):
    p=np.asarray(values,float);order=np.argsort(p);q=np.empty(len(p))
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q


def targets(fixture,panel):
    root=fixture/"private_targets"
    values=pd.read_csv(root/"log2cpm.tsv.gz",sep="\t",dtype={"peak_id":str})
    raw=pd.read_csv(root/"raw_counts.tsv.gz",sep="\t",dtype={"peak_id":str})
    lib=pd.read_csv(root/"donor_library_totals.tsv",sep="\t",dtype={"donor_id":str})
    if not np.array_equal(values.peak_id,panel.peak_id.astype(str)) or not np.array_equal(raw.peak_id,values.peak_id):raise ValueError("Measured peak identities/order differ")
    if list(values.columns[1:])!=list(lib.donor_id) or list(raw.columns[1:])!=list(lib.donor_id):raise ValueError("Donor axes/order differ")
    if len(lib)!=138 or lib.donor_id.duplicated().any():raise ValueError("Source donor population differs")
    y=values.iloc[:,1:].to_numpy(float);counts=raw.iloc[:,1:].to_numpy(float);total=lib.full_peak_library_count_total.to_numpy(float)
    if not np.isfinite(counts).all() or not np.isfinite(y).all() or not np.isfinite(total).all() or np.any(counts<0) or np.any(counts!=np.floor(counts)) or np.any(total<=0):raise ValueError("Raw count/target/denominator contract differs")
    np.testing.assert_allclose(y,np.log2(1+counts/total[None]*1e6),rtol=0,atol=1e-12)
    return y,lib.donor_id.to_numpy()


def bootstrap(x,y,panel,training,validation,unit,draws):
    """One common donor draw changes both calibration and validation targets."""
    rng=np.random.default_rng(SEED);mse=[];rho=[];rawrho=[]
    for _ in range(draws):
        donor=rng.choice(y.shape[1],y.shape[1],replace=True);average=y[:,donor].mean(1)
        tr=resample_regions(panel,training,unit,rng,True);va=resample_regions(panel,validation,unit,rng,False)
        pred,_,_,mean=calibrate(x,average,tr);pred=np.column_stack([pred,mean])
        errors=np.mean((pred[va]-average[va,None])**2,axis=0)
        calibrated=correlations(pred[va],average[va]);native=correlations(np.column_stack([x,np.zeros(len(x))])[va],average[va])
        mse.append(errors[0]-errors[1:]);rho.append(calibrated[1:]-calibrated[0]);rawrho.append(native[1:]-native[0])
    return {"MSE_improvement":np.asarray(mse),"calibrated_spearman_improvement":np.asarray(rho),"raw_spearman_improvement":np.asarray(rawrho)}


def run(fixture,panel,raw_tracks,names,out,draws=1000):
    started=time.monotonic();y,donors=targets(fixture,panel);average=y.mean(1)
    if tuple(names)!=("original",)+ADAPTATION_NAMES or len(set(names))!=len(names):raise ValueError("Native model identities/order differ from complete prespecified set")
    training=np.flatnonzero(panel.role=="training");validation=np.flatnonzero(panel.role=="validation")
    if len(training)!=512 or len(validation)!=512 or np.any(panel.fold==0):raise ValueError("Fixed panel role counts differ")
    if not np.isfinite(raw_tracks).all() or raw_tracks.shape!=(11,1024,3):raise ValueError("Incomplete native-model/track predictions")
    x=raw_tracks.mean(2).T
    predicted,slopes,intercepts,mean=calibrate(x,average,training)
    keys=list(names)+["training_mean_target"];predicted=np.column_stack([predicted,mean])
    raw=np.column_stack([x,np.zeros(len(x))]);errors=np.mean((predicted[validation]-average[validation,None])**2,axis=0)
    native_rho=correlations(raw[validation],average[validation]);cal_rho=correlations(predicted[validation],average[validation])
    rows=[]
    for j,name in enumerate(keys):
        rows.append({"model":name,"regions":len(validation),"source_described_donors":len(donors),
            "RMSE_donor_mean_log2_peak_CPM":np.sqrt(errors[j]),"raw_predicted_scale_spearman":native_rho[j],
            "calibrated_spearman":cal_rho[j],"training_slope":slopes[j] if j<len(names) else 0.,
            "training_intercept":intercepts[j] if j<len(names) else mean[0],
            "input_order_reversed_by_calibration":bool(slopes[j]<0) if j<len(names) else False,
            "native_context_bp":2048,"assessment":"same_Currin_source_regional_liver_ATAC_retention_only"})
    pd.DataFrame(rows).to_csv(out/"measured_native_performance.tsv",sep="\t",index=False)
    table=panel.copy();table["observed_donor_mean_log2_peak_CPM"]=average
    for j,name in enumerate(keys):table[name+"__calibrated"]=predicted[:,j];table[name+"__raw_log_track_sum"]=raw[:,j]
    table.to_csv(out/"regional_predictions.tsv.gz",sep="\t",index=False)
    comparisons=[]
    point={"MSE_improvement":errors[0]-errors[1:],"calibrated_spearman_improvement":cal_rho[1:]-cal_rho[0],
        "raw_spearman_improvement":native_rho[1:]-native_rho[0]}
    for unit in ("input_component","chrom"):
        samples=bootstrap(x,y,panel,training,validation,unit,draws)
        for endpoint,sampled in samples.items():
            for j,name in enumerate(keys[1:]):
                valid=sampled[:,j][np.isfinite(sampled[:,j])];value=point[endpoint][j]
                # Constant mean predictor has no defined rank; retain it as p=1 in complete family.
                low,high=np.quantile(valid,(.025,.975)) if len(valid)==draws else (np.nan,np.nan)
                p=(1+np.sum(np.abs(valid-value)>=abs(value)))/(len(valid)+1) if np.isfinite(value) and len(valid)==draws else 1.
                comparisons.append({"model":name,"comparator":"original","endpoint":endpoint,"point":value,"low95":low,"high95":high,
                    "resampling_unit":unit,"paired_donor_resampling":True,"calibration_refit_each_draw":True,
                    "draws":draws,"estimable_draws":len(valid),"nominal_p":p,
                    "status":"source_development_fixed_adapted_models" if np.isfinite(value) and len(valid)==draws else "undefined_rank_or_bootstrap_CI_unavailable"})
    result=pd.DataFrame(comparisons)
    for (_, _),indices in result.groupby(["resampling_unit","endpoint"]).groups.items():
        result.loc[indices,"BH_q_complete_exploratory_family"]=bh(result.loc[indices,"nominal_p"])
        result.loc[indices,"planned_family_n"]=len(indices)
    result.to_csv(out/"measured_native_paired_uncertainty.tsv",sep="\t",index=False)
    receipt={"status":"measured_regional_ATAC_retention_complete","seconds":time.monotonic()-started,"models":11,
        "calibration_regions":512,"validation_regions":512,"donors":len(donors),"training_chromosomes":sorted(panel.iloc[training].chrom.unique()),
        "validation_chromosomes":sorted(panel.iloc[validation].chrom.unique()),"validation_input_components":int(panel.iloc[validation].input_component.nunique()),
        "full_peak_library_denominator":"includes_all_deposited_peak_rows_not_offpeak_reads; no_library_size_predictor",
        "bootstrap":"same_donor_draw_in_training_and_validation_and_all_models; independently_resampled_training_and_validation_regions; training_only_affine_refit",
        "uncertainty_limit":"fixed_adapted_models; adaptation_training_and_selection_uncertainty_not_propagated; actual2kb_components_not_claimed_LD_independent",
        "native_scope":"regional_adult_liver_ATAC_2kb_reference_inputs; no_DNase_or_whole_model_function_claim",
        "source_scope":"same_Currin_development_family; not_independent_cohort; foundation_exposure_unresolved",
        "no_donor_prediction_or_personal_risk_score":True,"source_GC_correction_and_VST_not_used":True}
    (out/"measured_native_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    return receipt
