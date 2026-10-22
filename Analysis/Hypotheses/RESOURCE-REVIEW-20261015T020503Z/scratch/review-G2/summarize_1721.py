import pandas as pd, numpy as np, glob, os
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN=ROOT+"/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z"
J=RUN+"/jobs/1721/fits"; S=RUN+"/scratch/review-G2"
ext=pd.read_csv(ROOT+"/figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv",sep="\t",
    usecols=["gene_id_versioned","gene_id_base","gene_name","logFC","FDR","contrast","axis"])
ext=ext[ext.axis=="fibrosis"]
adj=pd.read_csv(ROOT+"/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_all_gene_results.tsv",sep="\t",
    usecols=["gene_id_versioned","gene_id_base","gene_name","logFC","FDR","contrast"])
rel=pd.concat([ext.drop(columns="axis"),adj])
ver2base=dict(zip(rel.gene_id_versioned,rel.gene_id_base)); ver2name=dict(zip(rel.gene_id_versioned,rel.gene_name))
relsig={c:set(g.loc[g.FDR<0.05,"gene_id_versioned"]) for c,g in rel.groupby("contrast")}
relfc={c:g.set_index("gene_id_versioned").logFC for c,g in rel.groupby("contrast")}
rows=[]; fits={}
for f in sorted(glob.glob(J+"/*.genes.tsv.gz")):
    lab,con=os.path.basename(f).replace(".genes.tsv.gz","").split("__")
    g=pd.read_csv(f,sep="\t"); fits[(lab,con)]=g
    sig=set(g.loc[g.FDR<0.05,"gene_id_versioned"]); R=relsig[con]
    fc=g.set_index("gene_id_versioned").logFC
    both=list(sig&R)
    rows.append(dict(label=lab,contrast=con,released=len(R),refit=len(sig),retained=len(sig&R),
        retained_frac=round(len(sig&R)/len(R),4) if R else np.nan,new=len(sig-R),
        lost=len(R-sig),sign_agree_retained=int((np.sign(fc[both])==np.sign(relfc[con][both])).sum()),
        logFC_pearson_all=round(np.corrcoef(fc.loc[relfc[con].index],relfc[con])[0,1],4),
        logFC_slope_released_sig=round(np.polyfit(relfc[con][list(R)],fc[list(R)],1)[0],4) if R else np.nan))
out=pd.DataFrame(rows); out.to_csv(S+"/refit_vs_released_counts.tsv",sep="\t",index=False)
pd.set_option("display.width",250); print(out.to_string(index=False))
# Fig 3F reselection
g=pd.read_csv(ROOT+"/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
pos=set(g.loc[g.coloc_best_susie_pp4>0.5,"ensembl"]); assert len(pos)==436
tr=["F0_to_F1","F1_to_F2","F2_to_F3","F3_to_F4"]
def select(tabs):
    # tabs: contrast -> df(gene_id_versioned, FDR, logFC)
    sel=set()
    for c in tr:
        t=tabs[c]; b=t.gene_id_versioned.map(ver2base)
        sel|=set(b[(t.FDR<0.05)&b.isin(pos)])
    marks={c:int(((tabs[c].FDR<0.05)&tabs[c].gene_id_versioned.map(ver2base).isin(sel)).sum()) for c in tr}
    return sel,marks
reltabs={c:adj[adj.contrast==c] for c in tr}
rsel,rmarks=select(reltabs); print("released 3F",len(rsel),4*len(rsel),sum(rmarks.values()),rmarks)
res=[]
for lab in ["joint","lib","gini","nf"]:
    need={"F0_to_F1":(lab,"F1_vs_F0"),"F1_to_F2":(lab,"F1_to_F2"),"F2_to_F3":(lab,"F2_to_F3"),"F3_to_F4":(lab,"F3_to_F4")}
    if not all(k in fits for k in need.values()): print(lab,"3F incomplete:",[c for c,k in need.items() if k not in fits]); continue
    tabs={c:fits[k] for c,k in need.items()}
    sel,marks=select(tabs)
    res.append(dict(label=lab,genes=len(sel),effects=4*len(sel),marks=sum(marks.values()),**{f"marks_{c}":marks[c] for c in tr},
        kept_of_54=len(sel&rsel),new=len(sel-rsel),dropped=";".join(sorted(ver2name.get(next(v for v in ver2base if ver2base[v]==x),x) for x in rsel-sel)),
        added=";".join(sorted(ver2name.get(next(v for v in ver2base if ver2base[v]==x),x) for x in sel-rsel))))
    # per-cell mark agreement on released 54 genes
    agree=0; lostm=0; gained=0
    for c in tr:
        a=reltabs[c].assign(b=reltabs[c].gene_id_versioned.map(ver2base)); a=a[a.b.isin(rsel)].set_index("gene_id_versioned").FDR<0.05
        t=tabs[c].set_index("gene_id_versioned").FDR.loc[a.index]<0.05
        agree+=int((a&t).sum()); lostm+=int((a&~t).sum()); gained+=int((~a&t).sum())
    res[-1].update(released54_marks_kept=agree,released54_marks_lost=lostm,released54_marks_gained=gained)
pd.DataFrame(res).to_csv(S+"/fig3f_reselection.tsv",sep="\t",index=False)
print(pd.DataFrame(res).to_string(index=False))
# hero genes used in the text
hero=["GNMT","MAT1A","CYP2C19"]
hv=[v for v,n in ver2name.items() if n in hero]
hr=[]
for (lab,con),gtab in fits.items():
    x=gtab[gtab.gene_id_versioned.isin(hv)]
    for _,r in x.iterrows():
        rr=rel[(rel.contrast==con)&(rel.gene_id_versioned==r.gene_id_versioned)].iloc[0]
        hr.append(dict(gene=ver2name[r.gene_id_versioned],label=lab,contrast=con,logFC=round(r.logFC,4),FDR=r.FDR,released_logFC=round(rr.logFC,4),released_FDR=rr.FDR))
hr=pd.DataFrame(hr).sort_values(["gene","contrast","label"]); hr.to_csv(S+"/hero_genes_refit.tsv",sep="\t",index=False)
print(hr[hr.contrast.isin(["F1_vs_F0","F1_to_F2"])].to_string(index=False))
# F0->F1 marks among the released 54 genes under each available F1_vs_F0 refit
a=reltabs["F0_to_F1"].assign(b=lambda d:d.gene_id_versioned.map(ver2base)); a=a[a.b.isin(rsel)]
relm=set(a.loc[a.FDR<0.05,"gene_id_versioned"])
for (lab,con),gtab in fits.items():
    if con!="F1_vs_F0": continue
    t=gtab[gtab.gene_id_versioned.isin(a.gene_id_versioned)]
    m=set(t.loc[t.FDR<0.05,"gene_id_versioned"])
    print("F0->F1 marks in released 54:",lab,"refit",len(m),"of released",len(relm),"kept",len(m&relm))
