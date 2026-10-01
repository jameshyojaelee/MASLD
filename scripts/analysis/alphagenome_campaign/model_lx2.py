#!/usr/bin/env python3
"""Replicate-preserving treatment interaction comparison in LX2 and HepG2.

Uses local frozen AlphaGenome 2,048-bp variant-bin features, averaged across
forward/reverse-complement ALT-minus-REF differences. Labels remain paired
within experimental replicate. Inference concerns reporter constructs in
these cell lines, not participant fibrosis or treatment response.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge

PROJ = Path(__file__).resolve().parents[3]
ROW_AUTHORITY = PROJ/"Analysis/MASLD_Model_Benchmark/executions/model-check-220-21088696/contract/row_universe.tsv"
EMBEDDINGS = PROJ/"GWAS/finemapping/results/alphagenome_program/ad-arm2-mpra-20260915T114407Z/embeddings/alphagenome_frozen_trunk_2048bp.npz"
ALPHAS = (1.0, 10.0, 100.0, 1000.0, 10000.0)
ARMS = ("direct_interaction", "separate_conditions", "shared_condition_residual", "zero_interaction")


def projection(raw, train):
    scale = raw[train].std(axis=0)
    scale[scale < 1e-6] = 1.0
    # A linear, uncentered compression preserves allele antisymmetry.
    svd = TruncatedSVD(n_components=min(64, int(train.sum())-1, raw.shape[1]-1), random_state=1103)
    svd.fit(raw[train]/scale)
    return svd.transform(raw/scale)


def ridge_predict(x, y, train, test, alpha):
    model = Ridge(alpha=alpha, fit_intercept=False)
    model.fit(x[train], y[train])
    return model.predict(x[test])


def select_alpha(x, y, train, valid):
    return min(ALPHAS, key=lambda a:(np.mean((ridge_predict(x,y,train,valid,a)-y[valid])**2), -a))


def shared_predict(x, control, treated, train, test, alpha, ratio):
    # Control = x*w; treated = x*w+x*v. v has penalty alpha*ratio.
    # Scaling the residual design by sqrt(ratio) expresses unequal penalties
    # using ordinary ridge, without changing native outcome units.
    z = x[train]
    design = np.concatenate([np.concatenate([z,np.zeros_like(z)],axis=1),
                             np.concatenate([z,z/np.sqrt(ratio)],axis=1)],axis=0)
    model = Ridge(alpha=alpha, fit_intercept=False)
    model.fit(design, np.concatenate([control[train], treated[train]]))
    return x[test] @ (model.coef_[x.shape[1]:]/np.sqrt(ratio))


def uncertainty(replicates, predictions, blocks, resamples):
    """Paired locus resampling plus one shared replicate resampling per draw.

    The four experimental replicate indices are resampled jointly across all
    constructs and both conditions. Thus shared experimental perturbations and
    shared-control covariance are not replaced by independent per-row noise.
    Predictors are fixed OOF fits; this interval does not include fitting/model
    selection uncertainty and only four experimental replicates are available.
    """
    rng = np.random.default_rng(20260915)
    unique = np.asarray(sorted(set(blocks)))
    lookup = {b:np.flatnonzero(blocks == b) for b in unique}
    draws = {name:[] for name in ARMS}
    for _ in range(resamples):
        idx = np.concatenate([lookup[b] for b in rng.choice(unique, len(unique), replace=True)])
        repidx = rng.integers(0, replicates.shape[1], replicates.shape[1])
        observed = replicates[idx][:,repidx].mean(axis=1)
        for name in ARMS:
            draws[name].append(float(np.mean((predictions[name][idx]-observed)**2)))
    rows = []
    for name in ARMS[:-1]:
        improvement = np.asarray(draws["zero_interaction"])-np.asarray(draws[name])
        ci = np.quantile(improvement,(0.025,0.975))
        rows.append({"arm":name, "MSE_improvement_vs_zero_low95":float(ci[0]),
                     "MSE_improvement_vs_zero_high95":float(ci[1]),
                     "resamples":resamples, "seed":20260915,
                     "unit":"long_range_locus_and_joint_experimental_replicate",
                     "fit_selection_uncertainty":"not_included", "multiplicity":"exploratory_CI_no_confirmatory_test"})
    return rows


def run(args):
    args.out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    rows = pd.read_csv(args.row_authority, sep="\t", usecols=["element_id","outer_fold","long_range_block_id"]).drop_duplicates()
    if rows.element_id.duplicated().any():
        raise ValueError("Operative row authority has inconsistent split identities")
    rows["fold"] = rows.outer_fold.astype(str).str.replace("fold-","",regex=False).astype(int)
    if rows.groupby("long_range_block_id").fold.nunique().max() != 1:
        raise ValueError("A locus spans folds")
    with np.load(args.embeddings, allow_pickle=False) as archive:
        ids = archive["element_ids"].astype(str)
        if list(archive["allele_order"]) != ["REF","ALT","REF_RC","ALT_RC"]:
            raise ValueError("Unexpected allele order")
        e = archive["pooled_centre_bin"].astype(np.float32)
        raw = 0.5*((e[:,1]-e[:,0])+(e[:,3]-e[:,2]))
        extracted = archive["extracted"]
    if set(ids) != set(rows.element_id):
        raise ValueError("Representations and operative source-fold authority differ")
    rowmap = rows.set_index("element_id").loc[ids]
    frame = pd.read_csv(args.endpoints, sep="\t")
    metrics, selection, excluded, intervals, census = [], [], [], [], []
    for cell in ("LX2","HepG2"):
        source = frame.loc[frame.cell_line == cell].copy()
        if source.duplicated(["element_id","experimental_replicate"]).any():
            raise ValueError("Repeated source construct/replicate")
        eligible = source.eligible.astype(str).str.lower().eq("true")
        source = source.loc[eligible]
        per_rep = {col:source.pivot(index="element_id",columns="experimental_replicate",values=col).reindex(ids)
                   for col in ("effect_control","effect_treated","effect")}
        for col in per_rep:
            if not per_rep[col].columns.equals(per_rep["effect"].columns):
                raise ValueError("Condition and interaction replicate indices differ")
        # All four pairs are required for the primary equal-replication panel;
        # incomplete pairs are reported, and never imputed as zero interaction.
        okay = extracted & np.isfinite(raw).all(axis=1)
        for col in per_rep:
            if per_rep[col].shape[1] != 4:
                raise ValueError("Expected four source experimental replicate indices")
            okay &= np.isfinite(per_rep[col].to_numpy()).all(axis=1)
        for element in ids[~okay]:
            excluded.append({"cell_line":cell,"element_id":element,"reason":"incomplete_four_pairs_or_representation"})
        names = ids[okay]
        census.append({"cell_line":cell, "full_source_constructs":frame.loc[frame.cell_line == cell,"element_id"].nunique(),
            "source_constructs_with_eligible_replicates":source.element_id.nunique(),
            "representation_population":len(ids), "matched_four_pair_constructs":len(names),
            "representation_population_rule":args.population_rule,
            "pairing_assumption":"source_replicate_index_alignment; no_independent_donor_inference"})
        xraw = raw[okay]
        folds = rowmap.fold.to_numpy()[okay]
        blocks = rowmap.long_range_block_id.to_numpy()[okay]
        control_rep, treated_rep, diff_rep = (per_rep[c].to_numpy(dtype=float)[okay] for c in ("effect_control","effect_treated","effect"))
        np.testing.assert_allclose(treated_rep-control_rep,diff_rep,rtol=0,atol=1e-10)
        control, treated, y = control_rep.mean(axis=1), treated_rep.mean(axis=1), diff_rep.mean(axis=1)
        covariance = np.sum((control_rep-control[:,None])*(treated_rep-treated[:,None]),axis=1)/3
        var_mean = diff_rep.var(axis=1,ddof=1)/4
        reconstructed = (control_rep.var(axis=1,ddof=1)+treated_rep.var(axis=1,ddof=1)-2*covariance)/4
        np.testing.assert_allclose(var_mean,reconstructed,atol=1e-12)
        pred = {name:np.zeros(len(y)) if name == "zero_interaction" else np.full(len(y),np.nan) for name in ARMS}
        for held in range(5):
            inner = (held+1)%5
            train, test = folds != held, folds == held
            inner_train, inner_valid = train & (folds != inner), folds == inner
            xi = projection(xraw,inner_train)
            xo = projection(xraw,train)
            direct_alpha = select_alpha(xi,y,inner_train,inner_valid)
            ca = select_alpha(xi,control,inner_train,inner_valid)
            ta = select_alpha(xi,treated,inner_train,inner_valid)
            shared_alpha,ratio = min(((a,r) for a in ALPHAS for r in (0.1,1.0,10.0)),
                key=lambda ar:np.mean((shared_predict(xi,control,treated,inner_train,inner_valid,*ar)-y[inner_valid])**2))
            pred["direct_interaction"][test] = ridge_predict(xo,y,train,test,direct_alpha)
            pred["separate_conditions"][test] = ridge_predict(xo,treated,train,test,ta)-ridge_predict(xo,control,train,test,ca)
            pred["shared_condition_residual"][test] = shared_predict(xo,control,treated,train,test,shared_alpha,ratio)
            selection.append({"cell_line":cell,"held_fold":held,"inner_fold":inner,
                "direct_alpha":direct_alpha,"control_alpha":ca,"treated_alpha":ta,
                "shared_alpha":shared_alpha,"residual_penalty_ratio":ratio,
                "train_n":int(train.sum()),"test_n":int(test.sum()),"preprocessing":"inner_and_outer_train_only"})
        table = pd.DataFrame({"element_id":names,"fold":folds,"long_range_block_id":blocks,
            "observed_interaction":y,"paired_mean_variance":var_mean,"control_treated_sample_covariance":covariance,
            "experimental_replicates":4,**pred})
        table.to_csv(args.out/f"{cell}_oof.tsv",sep="\t",index=False)
        for arm in ARMS:
            if not np.isfinite(pred[arm]).all():
                raise ValueError("Nonfinite prediction")
            metrics.append({"cell_line":cell,"arm":arm,"constructs":len(y),"locus_blocks":len(set(blocks)),
                "experimental_replicates":4,"RMSE_log2_activity_interaction":float(np.sqrt(np.mean((pred[arm]-y)**2))),
                "signed_spearman":float(spearmanr(y,pred[arm]).statistic) if np.std(pred[arm]) else np.nan,
                "evidence":"nested_development_held_locus", "measurement_limit":"cell_line_reporter_not_patient_response"})
        ci = uncertainty(diff_rep,pred,blocks,args.bootstrap)
        intervals.extend([{"cell_line":cell,**r} for r in ci])
    pd.DataFrame(metrics).to_csv(args.out/"performance.tsv",sep="\t",index=False)
    pd.DataFrame(census).to_csv(args.out/"population_coverage.tsv",sep="\t",index=False)
    pd.DataFrame(selection).to_csv(args.out/"inner_selection.tsv",sep="\t",index=False)
    pd.DataFrame(excluded,columns=["cell_line","element_id","reason"]).to_csv(args.out/"excluded.tsv",sep="\t",index=False)
    pd.DataFrame(intervals).to_csv(args.out/"paired_intervals.tsv",sep="\t",index=False)
    (args.out/"analysis.json").write_text(json.dumps({"source":str(args.endpoints),"representations":str(args.embeddings),
        "split_authority":str(args.row_authority),"historical_source_outer_fold_used":False,
        "representation_population_rule":args.population_rule,
        "elapsed_seconds":time.monotonic()-start,"shared_control_covariance_preserved":True,
        "protected_outcomes_read":False,"seed":1103,"bootstrap_seed":20260915,
        "biological_n":"four_experimental_replicates_per_construct_per_condition; no_participants",
        "measurement_unit":"difference_of_genomic_ALT_minus_REF_log2_RNA_over_DNA_contrasts",
        "source_counts":"existing_pseudocount0.5_log_ratio_not_negative_binomial_model",
        "status":"development_no_confirmatory_claim"},indent=2)+"\n")
    print(pd.DataFrame(metrics).to_string(index=False))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoints",type=Path,required=True)
    parser.add_argument("--embeddings",type=Path,default=EMBEDDINGS)
    parser.add_argument("--row-authority",type=Path,default=ROW_AUTHORITY)
    parser.add_argument("--population-rule",default="existing_1033_comparator_sequence_fixture; not_full4359_reporter_population")
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--bootstrap",type=int,default=2000)
    run(parser.parse_args())
