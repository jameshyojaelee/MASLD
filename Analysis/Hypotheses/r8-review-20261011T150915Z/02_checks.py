import pandas as pd, numpy as np, random
pd.set_option("display.width",250); pd.set_option("display.max_columns",80); pd.set_option("display.max_colwidth",60); pd.set_option("display.max_rows",200)
C="RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/"
g=pd.read_parquet(C+"passport_gene_index.parquet"); e=pd.read_parquet(C+"passport_evidence_long.parquet")
deg=pd.read_csv("RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv")
deg["ens"]=deg.gene.str.replace(r"\..*","",regex=True)
col=pd.read_csv("GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
ec=pd.read_csv("RNA-seq/results/manuscript_release/candidates/evidence-classes-corrected-r8-20261009T142027Z/evidence_class_table.tsv",sep="\t")
gen=pd.read_csv(C+"preflight/gen_identity_adjudication_v1/gen_frozen_classes_ensembl_adjudicated.tsv",sep="\t",usecols=lambda c: c not in("source_rows_json",))
q=pd.read_csv(C+"passport_identity_quarantine.tsv",sep="\t")
print("== sizes: genes",len(g),"evidence",len(e),"deg",len(deg),"deg unique ens",deg.ens.nunique(),"coloc",len(col),"coloc uniq",col.ensembl.nunique(),"evclass",len(ec),"gen",len(gen),"quarantine",len(q))
print(g.primary_evidence_class.value_counts().to_string())
print("call_state:",e.call_state.value_counts().to_dict()); print("neg rule non-null:",e.negative_call_rule_id.notna().sum(), e.negative_call_passed.notna().sum())
w=e.pivot(index="ensembl_id",columns="evidence_domain",values=["call_state","testability_state","estimate","q_value","phenotype","dataset_id"])
w.columns=[a+"_"+b for a,b in w.columns]; w=w.reset_index().merge(g[["ensembl_id","symbol","primary_evidence_class"]],on="ensembl_id")
print("\n== class x (genetics call/testability, tx call/testability)")
print(w.groupby(["primary_evidence_class","call_state_genetics","testability_state_genetics","call_state_transcriptomics","testability_state_transcriptomics"]).size().to_string())
# bulk calls
deg["call"]=(deg.padj<0.05)&(deg.logFC.abs()>0.5)
print("\n== bulk: rows",len(deg),"calls",deg.call.sum(),"up",((deg.logFC>0)&deg.call).sum(),"down",((deg.logFC<0)&deg.call).sum())
sup_tx=set(w.loc[w.call_state_transcriptomics=="supported","ensembl_id"])
degcall=set(deg.loc[deg.call,"ens"])
print("catalog tx supported",len(sup_tx),"deg calls",len(degcall),"in deg not catalog-supported",sorted(degcall-sup_tx),"catalog-supported not deg call",sorted(sup_tx-degcall))
miss=deg[deg.ens.isin(degcall-sup_tx)][["gene","symbol","logFC","padj"]]; print(miss.to_string())
print("in catalog at all?", g[g.ensembl_id.isin(degcall-sup_tx)][["ensembl_id","symbol","primary_evidence_class"]].to_string())
print("quarantine cols",list(q.columns)); 
for c in q.columns:
    hit=q[q[c].astype(str).str.contains("|".join(list(degcall-sup_tx))+"|"+"|".join(miss.symbol.astype(str)) if len(miss) else "ZZZZ",na=False)]
    if len(hit): print("quarantine hit via",c); print(hit.iloc[:, :8].to_string()); break
print("deg genes not in catalog:",len(set(deg.ens)-set(g.ensembl_id)),"catalog tx testable:",(w.testability_state_transcriptomics=="testable").sum())
# estimates agree?
m=w.merge(deg[["ens","logFC","padj","call"]],left_on="ensembl_id",right_on="ens",how="left")
tt=m[m.testability_state_transcriptomics=="testable"]
print("tx testable rows w/o deg match",tt.logFC.isna().sum(),"max |logFC diff|",(tt.estimate_transcriptomics.astype(float)-tt.logFC).abs().max(),"max |q diff|",(tt.q_value_transcriptomics.astype(float)-tt.padj).abs().max(),"call mismatch",((tt.call_state_transcriptomics=="supported")!=tt.call.fillna(False)).sum())
# coloc
print("\n== coloc cols"); print(col.susie_state_t12.value_counts(dropna=False).to_string()); print(col.driving_tier.value_counts(dropna=False).to_string())
for thr in (0.8,):
    s=col[col.coloc_best_susie_pp4>=thr]; print("susie pp4>=",thr,len(s),"uniq",s.ensembl.nunique(),"abf>=0.8",(col.coloc_best_abf_pp4>=0.8).sum(),"abf>=.8 & not susie",((col.coloc_best_abf_pp4>=0.8)&~(col.coloc_best_susie_pp4>=0.8)).sum())
sus=set(col.loc[col.coloc_best_susie_pp4>=0.8,"ensembl"]); sup_g=set(w.loc[w.call_state_genetics=="supported","ensembl_id"])
print("catalog genetics supported",len(sup_g),"susie>=.8 not catalog",sorted(sus-sup_g),"catalog not susie",sorted(sup_g-sus))
print("class of genetics-supported:",w[w.call_state_genetics=="supported"].primary_evidence_class.value_counts().to_dict())
print("class of tx-supported:",w[w.call_state_transcriptomics=="supported"].primary_evidence_class.value_counts().to_dict())
print("genetics supported & untested class:"); print(w[(w.call_state_genetics=="supported")&(w.primary_evidence_class=="untested")][["symbol","ensembl_id","estimate_genetics","phenotype_genetics","call_state_transcriptomics","testability_state_transcriptomics"]].to_string())
print("coloc genes not in catalog",len(set(col.ensembl)-set(g.ensembl_id)), "catalog genetics testable",(w.testability_state_genetics=="testable").sum(), "of which in coloc", w[w.testability_state_genetics=="testable"].ensembl_id.isin(col.ensembl).sum())
mc=w.merge(col,left_on="ensembl_id",right_on="ensembl",how="left")
sg=mc[mc.call_state_genetics=="supported"]; print("supported: max|pp4 diff|",(sg.estimate_genetics.astype(float)-sg.coloc_best_susie_pp4).abs().max(),"dataset mismatch",(sg.dataset_id_genetics!=sg.susie_driving_gwas).sum(), "min pp4",sg.estimate_genetics.min())
ig=mc[(mc.call_state_genetics=="indeterminate")]; print("indeterminate genetics n",len(ig),"estimate NA",ig.estimate_genetics.isna().sum(),"max susie",ig.coloc_best_susie_pp4.max(),"max estimate",ig.estimate_genetics.max(), "estimate==susie",np.isclose(ig.estimate_genetics.astype(float),ig.coloc_best_susie_pp4).sum(),"estimate==abf",np.isclose(ig.estimate_genetics.astype(float),ig.coloc_best_abf_pp4).sum())
# heroes
print("\n== heroes")
H=["MTARC1","HKDC1","GLP1R","THRB"]
print(w[w.symbol.isin(H)].drop(columns=["ensembl_id"]).T.to_string())
print(deg[deg.symbol.isin(H)][["gene","symbol","logFC","SE","P.Value","padj","treat_fdr","AveExpr","call"]].to_string())
he=g[g.symbol.isin(H)].ensembl_id
print(col[col.ensembl.isin(he)].T.to_string())
print(ec[ec.symbol.isin(H)].T.to_string())
print(gen[gen.gene_symbol.isin(H)][["gene_symbol","bulk_tested","primary_genetic_map_tested","joint_testable","primary_genetic","established_state_associated","static_class","bulk_logFC","bulk_padj","coloc_best_susie_pp4","coloc_best_abf_pp4","driving_gwas","driving_trait","abf_driving_gwas","genetic_untestable_shared_posterior"]].T.to_string())
# overview cards
print(w[w.symbol.isin(["7SK","A1BG","A2M-AS1","A2ML1","AADACP1","AADAT","A1CF","AACS"])][["symbol","primary_evidence_class","call_state_genetics","testability_state_genetics","call_state_transcriptomics","testability_state_transcriptomics"]].to_string())
# evidence class table vs catalog
ec["ens"]=ec.ensembl_bulk.fillna(ec.ensembl_atlas)
print("\n== evclass table classes:",ec.primary_evidence_class.value_counts().to_dict())
MAP={"convergent":"concordant","genetic_only":"genetically_anchored","disease_state_only":"established_state_associated","neither":"unresolved","not_jointly_testable":"untested","genetic_untestable_shared_posterior":"untested"}
x=g.merge(ec[["ens","primary_evidence_class","symbol"]].drop_duplicates("ens"),left_on="ensembl_id",right_on="ens",how="left",suffixes=("","_ec"))
print("catalog genes without evclass row",x.primary_evidence_class_ec.isna().sum(),"class disagreements",(x.primary_evidence_class_ec.map(MAP)!=x.primary_evidence_class)[x.primary_evidence_class_ec.notna()].sum())
print(pd.crosstab(x.primary_evidence_class,x.primary_evidence_class_ec.fillna("NO_ROW")).to_string())
# spot check
random.seed(20260923)
rows=[]
for cls in ["genetically_anchored","concordant","established_state_associated","unresolved","untested"]:
    ids=sorted(g.loc[g.primary_evidence_class==cls,"ensembl_id"]); rows+=random.sample(ids,4)
dd=deg.set_index("ens"); cc=col.set_index("ensembl")
out=[]
for ens in rows:
    r=g[g.ensembl_id==ens].iloc[0]
    bt=ens in dd.index; lf=dd.loc[ens,"logFC"] if bt else np.nan; pa=dd.loc[ens,"padj"] if bt else np.nan
    bcall=bool(bt and pa<0.05 and abs(lf)>0.5)
    gt=ens in cc.index; sp=cc.loc[ens,"coloc_best_susie_pp4"] if gt else np.nan; ab=cc.loc[ens,"coloc_best_abf_pp4"] if gt else np.nan; st=cc.loc[ens,"susie_state_t12"] if gt else ""
    gcall=bool(gt and sp>=0.8)
    if not(bt and gt): rc="untested"
    elif gcall and bcall: rc="concordant"
    elif gcall: rc="genetically_anchored"
    elif st not in ("tested","no_susie_result") and not gcall and str(st).startswith("untest"): rc="untested(shared_posterior)"
    elif bcall: rc="established_state_associated"
    else: rc="unresolved"
    out.append(dict(symbol=r.symbol,ens=ens,catalog=r.primary_evidence_class,recomputed=rc,bulk_in_fit=bt,logFC=lf,padj=pa,coloc_in_table=gt,susie_pp4=sp,abf_pp4=ab,susie_state=st,ok=rc.split("(")[0]==r.primary_evidence_class))
o=pd.DataFrame(out); print("\n== spot check"); print(o.to_string()); print("agree",o.ok.sum(),"/",len(o))
o.to_csv("Analysis/Hypotheses/r8-review-20261011T150915Z/spot_check_20.tsv",sep="\t",index=False)
# whole-catalog recompute using same rule
