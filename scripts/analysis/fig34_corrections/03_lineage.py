"""Fixed-program follow-up within content lineages; stream counts on CPU.

Corrected counts use each cell's original full-transcriptome library exposure.
Raw gene centering/scaling is frozen for the paired corrected comparison.
This is a projected weighted expression score, not a native Hotspot score.
"""
import sys
import json
import platform
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
from anndata.io import read_elem
from common import ROOT, PROGRAM, FOCAL, SEED, read, write, fit, bh, meta, plot_setup, snapshot

AMBIENT=ROOT/"Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2"
RAW=ROOT/"Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/work"
TARGET={list(FOCAL)[0]:"Fibroblasts",list(FOCAL)[1]:"Cholangiocytes"}
LINEAGES=["Hepatocytes","Fibroblasts","Cholangiocytes","Endothelial cells","Macrophages","T cells"]


def aggregate(selected,codes,confidence,nlevels):
    result={}
    for filt in ["all","confidence_0.90"]:
        keep=(codes>=0)&((confidence>=.90) if filt!="all" else True)
        ii=np.flatnonzero(keep);cc=codes[keep]
        ind=sp.csr_matrix((np.ones(len(ii)),(cc,np.arange(len(ii)))),shape=(nlevels,len(ii)))
        result[filt]=((ind@selected[:,ii].T).toarray(),np.bincount(cc,minlength=nlevels))
    return result


def main(out):
    out.mkdir(parents=True,exist_ok=False)
    member=read(PROGRAM/"program_membership_v2.tsv")
    weights={uid:member[member.program_uid.eq(uid)].groupby("source_gene").original_l1_weight.sum() for uid in FOCAL}
    genes=sorted(set().union(*(w.index for w in weights.values())))
    rawgenes=(RAW/"genes.txt").read_text().splitlines()
    scoregenes=(AMBIENT/"work/score_genes.txt").read_text().splitlines()
    rawpos=[rawgenes.index(g) for g in genes];decpos=[scoregenes.index(g) for g in genes]
    assert all(np.isclose(w.sum(),1) for w in weights.values())
    with h5py.File(ROOT/"Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad") as h:
        conf=pd.Series(np.asarray(read_elem(h["obs/cell_type_conf"]),float),index=read_elem(h["obs/_index"]).astype(str))
    diagnoses=read(PROGRAM.parent/"metadata/donor_metadata_stage_corrected_donor.tsv").set_index("sample")
    outputs=[];coverage=[];pairchanges=[]
    for cohort in ["GSE202379","GSE244832"]:
        print("Streaming",cohort,flush=True)
        pairing=pd.read_csv(ROOT/f"data/{cohort}/metadata/donor_pairing.csv",dtype=str)
        donor={s.strip():cohort+"_"+r.donor_id for r in pairing.itertuples() for s in r.rna_srrs.split(";") if s.strip()}
        md=pd.read_csv(AMBIENT/f"work/{cohort}/meta.csv.gz",dtype=str)
        md["donor"]=md["sample"].map(donor);assert md.donor.notna().all()
        md["confidence"]=conf.reindex(md.cell_id).to_numpy();assert md.confidence.notna().all()
        keys=md.donor+"||"+md.cell_type
        levels=pd.Index(sorted(keys[md.cell_type.isin(LINEAGES)].unique()))
        codes=levels.get_indexer(keys);nlevels=len(levels)
        frame=pd.DataFrame([k.split("||") for k in levels],columns=["donor","lineage"])
        frame["dataset"]=cohort;frame["disease"]=frame.donor.map(diagnoses.disease_stage_coarse)
        assert frame.disease.notna().all()
        cells=pd.read_csv(AMBIENT/f"work/{cohort}/dec_cells.csv.gz",dtype=str)
        available=md.cell_id.isin(cells.cell_id).to_numpy()
        acc={f:[np.zeros((nlevels,len(genes))),np.zeros(nlevels,int)] for f in ["all","confidence_0.90"]}
        matched_acc={f:[np.zeros((nlevels,len(genes))),np.zeros(nlevels,int)] for f in acc}
        totals=np.zeros(len(md));dims=json.loads((RAW/cohort/"dims.json").read_text());offset=0
        for chunk in dims["chunks"]:
            prefix=RAW/cohort/f"chunk{chunk['chunk']}"
            ix=np.fromfile(str(prefix)+"_indices.bin",dtype="<i4")
            ip=np.fromfile(str(prefix)+"_indptr.bin",dtype="<i4")
            val=np.fromfile(str(prefix)+"_data.bin",dtype="<i4")
            x=sp.csc_matrix((val,ix,ip),shape=(len(rawgenes),int(chunk["ncells"])))
            end=offset+x.shape[1];tt=np.asarray(x.sum(axis=0)).ravel();totals[offset:end]=tt
            small=x[rawpos].astype(float)@sp.diags(np.divide(1e4,tt,out=np.zeros_like(tt,dtype=float),where=tt>0));small.data=np.log1p(small.data)
            for f,(s,n) in aggregate(small,codes[offset:end],md.confidence.to_numpy()[offset:end],nlevels).items():acc[f][0]+=s;acc[f][1]+=n
            matched_codes=np.where(available[offset:end],codes[offset:end],-1)
            for f,(s,n) in aggregate(small,matched_codes,md.confidence.to_numpy()[offset:end],nlevels).items():matched_acc[f][0]+=s;matched_acc[f][1]+=n
            offset=end
        assert offset==len(md)
        dd=json.loads((AMBIENT/f"work/{cohort}/dec_dims.json").read_text());assert dd["correction_status"]=="corrected"
        cells=pd.read_csv(AMBIENT/f"work/{cohort}/dec_cells.csv.gz",dtype=str)
        positions=pd.Index(md.cell_id).get_indexer(cells.cell_id);assert (positions>=0).all()
        dec_acc={f:[np.zeros((nlevels,len(genes))),np.zeros(nlevels,int)] for f in acc}
        dp=AMBIENT/f"work/{cohort}"
        ip=np.memmap(dp/"dec_indptr.bin",dtype="<i4",mode="r")
        ix=np.memmap(dp/"dec_indices.bin",dtype="<i4",mode="r")
        values=np.memmap(dp/"dec_data.bin",dtype="<f8",mode="r")
        for start in range(0,len(cells),5000):
            end=min(start+5000,len(cells));lo,hi=int(ip[start]),int(ip[end])
            x=sp.csc_matrix((np.asarray(values[lo:hi]),np.asarray(ix[lo:hi]),np.asarray(ip[start:end+1])-lo),shape=(len(scoregenes),end-start))
            pos=positions[start:end];tt=totals[pos]
            small=x[decpos]@sp.diags(np.divide(1e4,tt,out=np.zeros_like(tt),where=tt>0));small.data=np.log1p(small.data)
            for f,(s,n) in aggregate(small,codes[pos],md.confidence.to_numpy()[pos],nlevels).items():dec_acc[f][0]+=s;dec_acc[f][1]+=n
        for filt in acc:
            for minimum in [50,20]:
                sums,n=acc[filt];dsums,dn=dec_acc[filt]
                keep=n>=minimum;expr=sums[keep]/n[keep,None]
                center=expr.mean(axis=0);scale=expr.std(axis=0,ddof=1)
                assert np.isfinite(scale).all() and (scale>0).all()
                rawz=(expr-center)/scale
                corrected=np.divide(dsums[keep],dn[keep,None],out=np.full_like(dsums[keep],np.nan),where=dn[keep,None]>0)
                corrz=(corrected-center)/scale
                ff=frame.loc[keep].copy();ff["n_cells"]=n[keep];ff["n_corrected_cells"]=dn[keep]
                ff["minimum_cells"]=minimum;ff["annotation_filter"]=filt
                # Compare the identical available cells; never substitute raw counts
                # for missing corrected cells. Incomplete correction remains explicit.
                msums,mn=matched_acc[filt];assert np.array_equal(mn,dn)
                common=ff.n_corrected_cells.ge(minimum).to_numpy()
                matched=np.divide(msums[keep],mn[keep,None],out=np.full_like(msums[keep],np.nan),where=mn[keep,None]>0)
                matchedz=(matched-center)/scale
                ff["correction_cell_coverage"]=ff.n_corrected_cells/ff.n_cells
                ff["correction_complete"]=ff.n_cells.eq(ff.n_corrected_cells)
                total_by_donor=pd.Series(n,index=frame.donor).groupby(level=0).sum()
                ff["captured_lineage_fraction"]=ff.n_cells/ff.donor.map(total_by_donor)
                for uid,w in weights.items():
                    ww=w.reindex(genes,fill_value=0).to_numpy()
                    assert len(w)>=8 and np.isclose(ww.sum(),1)
                    coverage.append(dict(dataset=cohort,program_uid=uid,minimum_cells=minimum,annotation_filter=filt,n_genes=len(w),retained_weight=ww.sum()))
                    for arm,zz in [("raw_full",rawz),("raw",matchedz),("ambient_corrected_original_exposure",corrz)]:
                        z=ff.copy();z["program_uid"]=uid;z["arm"]=arm;z["score"]=zz@ww
                        z["common_corrected_cells"]=common
                        outputs.append(z)
                # Store both weighted expressions and disjoint membership scores for the paired-lineage sensitivity.
                shared=set(weights[list(FOCAL)[0]].index)&set(weights[list(FOCAL)[1]].index)
                for uid,w in weights.items():
                    dis=w.drop(list(shared),errors="ignore");retained=dis.sum()
                    z=ff.copy();z["program_uid"]=uid;z["retained_disjoint_weight"]=retained
                    z["score_disjoint"]=rawz@dis.reindex(genes,fill_value=0).to_numpy()/retained if retained>=.8 else np.nan
                    pairchanges.append(z)
    s=pd.concat(outputs,ignore_index=True);write(s,out/"lineage_scores.tsv.gz")
    write(pd.DataFrame(coverage),out/"weight_coverage.tsv")
    dis=pd.concat(pairchanges,ignore_index=True);write(dis,out/"disjoint_lineage_scores.tsv.gz")
    # Exact reproduction of the stored raw localization scores, on the same cohort and donor-lineage gate.
    stored=read(ROOT/"Analysis/SingleCell/candidates/cross-lineage-specificity-complete-atlas-candidate-2026-08-15-v2/results/hero_lineage_scores.tsv.gz")
    stored["annotation_filter"]=stored.annotation_filter.replace({"all_annotated_cells":"all","cell_type_conf_ge_0.90":"confidence_0.90"})
    checks=[]
    for uid in FOCAL:
        z=s[s.arm.eq("raw_full")&s.program_uid.eq(uid)].merge(stored,on=["donor","lineage","dataset","minimum_cells","annotation_filter"])
        checks.append(dict(program_uid=uid,n=len(z),max_abs_difference=(z.score-z[uid]).abs().max()))
    check=pd.DataFrame(checks);write(check,out/"raw_score_reproduction.tsv");assert check.max_abs_difference.max()<1e-8
    rows=[]
    for (uid,cohort,minimum,filt,arm),g in s.groupby(["program_uid","dataset","minimum_cells","annotation_filter","arm"]):
        g=g[g.lineage.eq(TARGET[uid]) & (g.common_corrected_cells | (arm=="raw_full"))]
        for contrast,ref,test in [("primary_steatohepatitis_vs_steatosis","Steatosis","Steatohepatitis"),("secondary_steatohepatitis_vs_healthy","Healthy","Steatohepatitis")]:
            z=g[g.disease.isin([ref,test])];counts=z.disease.value_counts()
            if counts.get(ref,0)>=3 and counts.get(test,0)>=3:
                x=np.column_stack([np.ones(len(z)),z.disease.eq(test).to_numpy(float)])
                r=fit(z.score,x);ar=fit(np.log((z.n_cells+.5)/(z.n_cells/z.captured_lineage_fraction-z.n_cells+.5)),x)
            else:r=dict(n=len(z),estimable=False,beta=np.nan,se=np.nan,p=np.nan,hc3_se=np.nan,hc3_p=np.nan,failure_reason="fewer_than_three_donors_in_compared_group");ar={}
            r.update(program_uid=uid,program=FOCAL[uid],dataset=cohort,lineage=TARGET[uid],minimum_cells=minimum,annotation_filter=filt,arm=arm,contrast=contrast,n_reference=counts.get(ref,0),n_disease=counts.get(test,0),abundance_beta=ar.get("beta",np.nan),abundance_hc3_p=ar.get("hc3_p",np.nan))
            rows.append(r)
    fits=pd.DataFrame(rows)
    for _,idx in fits.groupby(["dataset","minimum_cells","annotation_filter","arm","contrast"]).groups.items():fits.loc[idx,"hc3_q"]=bh(fits.loc[idx,"hc3_p"],2)
    write(fits,out/"lineage_cohort_effects.tsv")
    metas=[]
    for keys,g in fits.groupby(["program_uid","minimum_cells","annotation_filter","arm","contrast"]):
        r=meta(g,robust=True);r.update(dict(zip(["program_uid","minimum_cells","annotation_filter","arm","contrast"],keys)));metas.append(r)
    metas=pd.DataFrame(metas)
    for _,idx in metas.groupby(["minimum_cells","annotation_filter","arm","contrast"]).groups.items():metas.loc[idx,"q_two_programs"]=bh(metas.loc[idx,"p"],2)
    write(metas,out/"lineage_meta_effects.tsv")
    # Same-donor cross-lineage relationship; residualize categorical diagnosis within cohort.
    rng=np.random.default_rng(SEED);paired=[]
    for (cohort,minimum,filt),g in dis.groupby(["dataset","minimum_cells","annotation_filter"]):
        a=g[g.program_uid.eq(list(FOCAL)[0])&g.lineage.eq("Fibroblasts")]
        b=g[g.program_uid.eq(list(FOCAL)[1])&g.lineage.eq("Cholangiocytes")]
        z=a.merge(b,on="donor",suffixes=("_ecm","_duct"));z=z[z.disease_ecm.isin(["Healthy","Steatosis","Steatohepatitis"])]
        if len(z)<10 or not np.isfinite(z[["score_disjoint_ecm","score_disjoint_duct"]]).all().all():continue
        y=z[["score_disjoint_ecm","score_disjoint_duct"]].to_numpy(float)
        x=np.column_stack([np.ones(len(z)),pd.get_dummies(z.disease_ecm,drop_first=True).to_numpy(float)])
        def correlation(ii):
            xx=x[ii];yy=y[ii]
            if np.linalg.matrix_rank(xx)<xx.shape[1]:return np.nan
            rr=yy-xx@np.linalg.lstsq(xx,yy,rcond=None)[0]
            return np.corrcoef(rr.T)[0,1]
        draws=[correlation(rng.integers(0,len(z),len(z))) for _ in range(2000)]
        lo,hi=np.nanquantile(draws,[.025,.975]);paired.append(dict(dataset=cohort,minimum_cells=minimum,annotation_filter=filt,n=len(z),r=correlation(np.arange(len(z))),low=lo,high=hi,n_bootstrap=int(np.isfinite(draws).sum()),retrospective=True))
    write(pd.DataFrame(paired),out/"paired_lineage_covariance.tsv")
    # KEY MESSAGE: disease contrasts are estimated within the content lineage, by cohort.
    plt=plot_setup();z=fits[(fits.minimum_cells==50)&fits.annotation_filter.eq("all")&fits.contrast.str.startswith("primary")&fits.estimable]
    fig,ax=plt.subplots(figsize=(4.6,2.4))
    for j,(_,r) in enumerate(z.iterrows()):
        ax.errorbar(r.beta,j,xerr=[[r.beta-r.hc3_low],[r.hc3_high-r.beta]],fmt="o",ms=3,color="#C9265E" if r.arm!="raw" else "#9E9E9E",lw=.5)
    ax.set_yticks(range(len(z)),[f"{r.program} / {r.dataset} / {r.arm.split('_')[0]} (n={r.n})" for r in z.itertuples()]);ax.axvline(0,color="#9E9E9E",lw=.5);ax.set_xlabel("Steatohepatitis − steatosis projected score (HC3 95% CI)")
    fig.tight_layout();fig.savefig(out/"S4_content_lineage_disease_effects.pdf");plt.close(fig)
    snapshot([PROGRAM/"program_membership_v2.tsv",PROGRAM/"program_registry_v2.tsv",AMBIENT/"work/manifest.json"],out/"frozen_inputs.tsv")
    (out/"analysis_definition.txt").write_text("Primary: corrected counts per original full-transcriptome exposure; raw donor-lineage gene z calibration held fixed.\nRaw/corrected models use identical cells and donors; no claims of complete cell autonomy.\nBH: two programs per contrast, arm and sensitivity; no selection across sensitivity arms.\nPaired lineage bootstrap is retrospective; pointwise intervals.\nSeed=20260929\n")
    (out/"environment.txt").write_text(f"Python {platform.python_version()}\nnumpy {np.__version__}\npandas {pd.__version__}\n")
    (out/"COMPLETE").write_text("Frozen raw scores reconstructed; targeted follow-up, not new discovery.\n")


if __name__=="__main__":
    from pathlib import Path
    main(Path(sys.argv[1]))
