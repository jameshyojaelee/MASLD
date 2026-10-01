"""Exact/grouped histology sensitivity, frozen scores and complete families.

Run on a compute node. No discovery, weight changes, or adopted-file writes.
"""
import sys
import platform
import numpy as np
import pandas as pd
from scipy.linalg import helmert
from common import ROOT, HAC, ML, AXES, FOCAL, SEED, PROGRAM, read, write, bh, design, fit, meta, snapshot, plot_setup


def main(out):
    out.mkdir(parents=True, exist_ok=False)
    manifest = ROOT / "figures/candidates/bulk-release-corrected-controls-20260925T175215Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv"
    system_path = HAC / "molecular_systems_inference/hac-system-inference-20260824T232135Z/tables/participant_system_scores.tsv.gz"
    program_path = ML / "programs/hotspot/hotspot_participant_scores.tsv.gz"
    m, a = read(manifest), read(AXES)
    m["source_control"] = m.group_binary.eq("Control")
    m["fibrosis_stage_exact"] = m.fibrosis_stage
    m["fibrosis_group_reported"] = m.fibrosis_stage.map(lambda v: f"F{int(v)}" if pd.notna(v) else "Not recorded")
    grouped = m.dataset.eq("GSE213621")
    m.loc[grouped, "fibrosis_group_reported"] = m.loc[grouped, "fibrosis_stage"].map({0:"Control",1:"F0–F1",2:"F2",3:"F3–F4"})
    m.loc[grouped, "fibrosis_stage_exact"] = np.where(m.loc[grouped,"fibrosis_stage"].eq(2), 2, np.nan)
    m["stage_resolution"] = np.where(grouped, "grouped", np.where(m.fibrosis_stage_exact.notna(), "exact", "not_recorded"))
    m["stage_group"] = m.fibrosis_group_reported
    m["stage_coarsened"] = m.fibrosis_stage_exact.map({0:"F0–F1",1:"F0–F1",2:"F2",3:"F3–F4",4:"F3–F4"})
    m.loc[grouped, "stage_coarsened"] = m.loc[grouped,"fibrosis_group_reported"]
    assert len(m)==844 and not m.sample_id.duplicated().any()
    assert m.loc[grouped & m.fibrosis_stage.eq(1),"fibrosis_stage_exact"].isna().all()
    write(m, out/"participant_histology.tsv")
    write(m.groupby(["dataset","fibrosis_group_reported","source_control"], dropna=False).size().rename("n_participants").reset_index(), out/"histology_census.tsv")
    scores = read(program_path).rename(columns={"program_uid":"feature_id"})
    scores["family"] = "programs"; scores["family_size"] = 117
    systems = read(system_path).query("aggregation_method == 'collection_balanced'").rename(columns={"community_id":"feature_id","system_score_z":"outcome_z"})
    systems["family"] = "systems"; systems["family_size"] = 43
    s = pd.concat([scores[["sample_id","dataset","feature_id","outcome_z","family","family_size"]], systems[["sample_id","dataset","feature_id","outcome_z","family","family_size"]]], ignore_index=True)
    d = s.merge(m, on=["sample_id","dataset"], validate="many_to_one").merge(a[["sample_id","fixed_projection_raw","signature_pc1_raw"]], on="sample_id", validate="many_to_one")
    fractions = []
    inputs = [manifest, AXES, system_path, program_path, PROGRAM/"program_membership_v2.tsv"]
    for cohort in ["GSE162694","GSE213621"]:
        p = ROOT/f"Analysis/Deconvolution/results/{cohort}/{cohort}_bayesprism_proportions.tsv"
        f = pd.read_csv(p, sep="\t", index_col=0)
        assert f.shape[1]==16 and np.allclose(f.sum(axis=1),1,atol=1e-5)
        # Same orthonormal coordinates as the established R contr.helmert design.
        ilr = np.log(np.maximum(f.to_numpy(float),1e-6)) @ (-helmert(16).T)
        z = pd.DataFrame(ilr,columns=[f"ilr_{i:02d}" for i in range(1,16)])
        z["sample_id"] = f.index; fractions.append(z); inputs.append(p)
    d = d.merge(pd.concat(fractions),on="sample_id",how="left",validate="many_to_one")
    snapshot(inputs,out/"inputs.tsv")
    rows, invariants, corr = [], [], []
    for cohort in ["GSE162694","GSE213621"]:
        c = d[d.dataset.eq(cohort) & d.fibrosis_stage.notna() & d.inferred_sex.notna()].copy()
        unique = c.drop_duplicates("sample_id")
        corr.append(dict(dataset=cohort,n=len(unique),pearson=unique.fixed_projection_raw.corr(unique.signature_pc1_raw)))
        for arm in ["original_categories","disease_only","exact_stage","coarsened_exact"]:
            if cohort=="GSE213621" and arm in ["exact_stage","coarsened_exact"]: continue
            z = c if arm=="original_categories" else c[~c.source_control]
            stage = "stage_coarsened" if arm=="coarsened_exact" else "stage_group"
            for axis in ["fixed_projection","signature_pc1"]:
                for (family,uid), g in z.groupby(["family","feature_id"],sort=False):
                    g = g[np.isfinite(g.outcome_z) & np.isfinite(g[axis+"_raw"])].copy()
                    ilr_cols = list(g.filter(regex="^ilr_"))
                    g = g.dropna(subset=ilr_cols)
                    if len(g):g["axis_z"]=(g[axis+"_raw"]-g[axis+"_raw"].mean())/g[axis+"_raw"].std(ddof=1)
                    for composition in [False,True]:
                        r = fit(g.outcome_z,design(g,stage,composition)) if len(g)>=20 else dict(n=len(g),estimable=False,beta=np.nan,se=np.nan,p=np.nan,hc3_p=np.nan,hc3_se=np.nan,failure_reason="fewer_than_20_participants")
                        r.update(dataset=cohort,feature_id=uid,family=family,family_size=117 if family=="programs" else 43,axis=axis,arm=arm,composition=composition)
                        rows.append(r)
                        if arm=="original_categories" and not composition and r["estimable"]:
                            old=fit(g.outcome_z,design(g,"fibrosis_stage",False))
                            invariants.append(dict(dataset=cohort,feature_id=uid,axis=axis,absolute_difference=abs(r["beta"]-old["beta"])))
    fits=pd.DataFrame(rows)
    keys=["dataset","family","axis","arm","composition"]
    for _,idx in fits.groupby(keys).groups.items():
        n=int(fits.loc[idx,"family_size"].iloc[0])
        fits.loc[idx,"q"]=bh(fits.loc[idx,"p"],n)
        fits.loc[idx,"hc3_q"]=bh(fits.loc[idx,"hc3_p"],n)
    write(fits,out/"continuum_cohort_effects.tsv")
    mets=[]
    for key,g in fits[fits.arm.isin(["original_categories","disease_only"])].groupby(["family","feature_id","axis","arm","composition"]):
        for robust in [False,True]:
            r=meta(g,robust=robust);r.update(dict(zip(["family","feature_id","axis","arm","composition"],key)));r["robust"]=robust;mets.append(r)
    mets=pd.DataFrame(mets)
    for key,idx in mets.groupby(["family","axis","arm","composition","robust"]).groups.items():mets.loc[idx,"q"]=bh(mets.loc[idx,"p"],117 if key[0]=="programs" else 43)
    write(mets,out/"continuum_meta_effects.tsv")
    iv=pd.DataFrame(invariants);assert iv.absolute_difference.max()<1e-10
    write(iv,out/"categorical_relabel_invariance.tsv");write(pd.DataFrame(corr),out/"axis_correlations.tsv")
    # Reproduce the stored cohort and meta estimates, not merely code execution.
    old=read(HAC/"translation/hac-translation-20260906T230500Z/program_cohort_effects.tsv")
    chk=fits.query("family=='programs' and arm=='original_categories'").merge(old,left_on=["dataset","feature_id","axis","composition"],right_on=["dataset","program_uid","axis_id",old.model.eq("stage_sex_composition")],suffixes=("_new","_old"))
    chk["absolute_beta_difference"]=(chk.beta_new-chk.beta_old).abs()
    assert chk.absolute_beta_difference.max()<1e-8
    write(chk[["dataset","feature_id","axis","composition","absolute_beta_difference"]],out/"stored_effect_reproduction.tsv")
    # Paired donor bootstrap; all 43 systems and both focal programs share each draw.
    rng=np.random.default_rng(SEED); bootrows=[]
    c=d[d.dataset.eq("GSE162694") & ~d.source_control & d.fibrosis_stage.notna() & (d.family.eq("systems")|d.feature_id.isin(FOCAL))].copy()
    wide=c.pivot(index="sample_id",columns="feature_id",values="outcome_z")
    md=c.drop_duplicates("sample_id").set_index("sample_id").loc[wide.index].copy()
    good=np.isfinite(wide).all(axis=1)&md.filter(regex="^ilr_").notna().all(axis=1)
    wide=wide.loc[good];md=md.loc[good];y=wide.to_numpy(float)
    for axis in ["fixed_projection","signature_pc1"]:
        md["axis_z"]=(md[axis+"_raw"]-md[axis+"_raw"].mean())/md[axis+"_raw"].std(ddof=1)
        for comp in [False,True]:
            x=design(md,"stage_group",comp);xc=design(md,"stage_coarsened",comp)
            observed=np.linalg.lstsq(xc,y,rcond=None)[0][1]-np.linalg.lstsq(x,y,rcond=None)[0][1]
            draws=[]
            for _ in range(2000):
                ii=rng.integers(0,len(md),len(md));xx=x[ii];cc=xc[ii]
                if np.linalg.matrix_rank(xx)!=xx.shape[1] or np.linalg.matrix_rank(cc)!=cc.shape[1]:continue
                draws.append(np.linalg.lstsq(cc,y[ii],rcond=None)[0][1]-np.linalg.lstsq(xx,y[ii],rcond=None)[0][1])
            limits=np.quantile(draws,[.025,.975],axis=0)
            for j,uid in enumerate(wide.columns):bootrows.append(dict(feature_id=uid,axis=axis,composition=comp,n=len(md),delta_coarse_minus_exact=observed[j],low=limits[0,j],high=limits[1,j],n_bootstrap=len(draws),n_attempted=2000,seed=SEED,interval="pointwise_percentile_conditional_on_frozen_scores"))
    write(pd.DataFrame(bootrows),out/"coarsening_paired_bootstrap.tsv")
    # KEY MESSAGE: the source cohort records grouped histology, not exact stages.
    plt=plot_setup()
    for family,data,ids,prefix in [("programs",scores,list(FOCAL),"S4G"),("nmf",None,[],"S4N")]:
        if data is None: continue
        plot=data[data.feature_id.isin(ids)].merge(m[["sample_id","fibrosis_group_reported","stage_resolution"]],on="sample_id",validate="many_to_one")
        plot=plot[plot.dataset.isin(["GSE162694","GSE213621"]) & plot.fibrosis_group_reported.ne("Not recorded")]
        sums=plot.groupby(["feature_id","dataset","fibrosis_group_reported"]).outcome_z.agg(["mean","std","count"]).reset_index()
        write(sums,out/f"{prefix}_stage_display.tsv")
        for uid in ids:
            fig,axs=plt.subplots(1,2,figsize=(4.6,1.9),sharey=True)
            for ax,cohort in zip(axs,["GSE162694","GSE213621"]):
                zz=sums[(sums.feature_id==uid)&(sums.dataset==cohort)]
                order=["F0","F1","F2","F3","F4"] if cohort=="GSE162694" else ["Control","F0–F1","F2","F3–F4"]
                zz=zz.set_index("fibrosis_group_reported").reindex(order)
                err=1.96*zz["std"]/np.sqrt(zz["count"])
                ax.errorbar(range(len(order)),zz["mean"],yerr=err,fmt="o",ms=3,color="#C9265E",lw=.6)
                if cohort=="GSE213621":ax.plot(0,zz["mean"].iloc[0],"o",color="#9E9E9E",ms=3)
                ax.set_xticks(range(len(order)),order);ax.set_title(cohort);ax.axhline(0,color="#9E9E9E",lw=.5)
            axs[0].set_ylabel("Program score (cohort SD)")
            fig.tight_layout();fig.savefig(out/f"{prefix}_{FOCAL[uid].replace(' ','_')}.pdf");plt.close(fig)
    focal=fits[fits.feature_id.isin(FOCAL)]
    write(focal,out/"focal_effects.tsv")
    (out/"environment.txt").write_text(f"Python {platform.python_version()}\nnumpy {np.__version__}\npandas {pd.__version__}\nseed={SEED}\n")
    (out/"COMPLETE").write_text("Frozen effects reproduced; categorical relabeling invariant; candidate only.\n")


if __name__=="__main__":
    from pathlib import Path
    main(Path(sys.argv[1]))
