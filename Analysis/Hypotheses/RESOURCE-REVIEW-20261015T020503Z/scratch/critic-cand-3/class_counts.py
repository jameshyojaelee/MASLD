import pandas as pd, numpy as np
C='figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
g=pd.read_parquet(C+'passport_gene_index.parquet')[['ensembl_id','symbol','primary_evidence_class']]
ev=pd.read_parquet(C+'passport_evidence_long.parquet')
gen=ev[(ev.evidence_domain=='genetics')&(ev.call_state=='supported')][['ensembl_id']].assign(gen_sup=True)
g=g.merge(gen,on='ensembl_id',how='left'); g['gen_sup']=g.gen_sup.fillna(False)
print('genetic supported by class'); print(pd.crosstab(g.primary_evidence_class,g.gen_sup))
d=pd.read_csv('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/critic-cand-3/r2_gene_direction_collapse.tsv',sep='\t')
d.columns=['ensembl_id','dir']
g=g.merge(d,on='ensembl_id',how='left')
E='figures/candidates/fig3e-s3-774-20261011T151616Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv'
st=pd.read_csv(E,sep='\t',usecols=['gene_id_base','axis','contrast','logFC','FDR'])
fib=st[st.axis=='fibrosis']
anyfib=fib.groupby('gene_id_base').FDR.min().lt(0.05).rename('stage_any_fib_vs_F0')
adj=st[st.axis=='fibrosis_adjacent'].groupby('gene_id_base').FDR.min().lt(0.05).rename('stage_any_adjacent')
f4=fib[fib.contrast=='F4_vs_F0'].set_index('gene_id_base')[['logFC','FDR']].rename(columns={'logFC':'f4_lfc','FDR':'f4_fdr'})
g=g.merge(anyfib,left_on='ensembl_id',right_index=True,how='left').merge(adj,left_on='ensembl_id',right_index=True,how='left').merge(f4,left_on='ensembl_id',right_index=True,how='left')
deg=pd.read_csv('RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv',keep_default_na=False,na_values=[''])
deg['ensembl_id']=deg.gene.str.replace(r'\..*$','',regex=True)
sym=deg[deg.symbol.astype(str).str.len()>0]; u=sym.groupby('symbol').size(); u=u[u==1].index
prot=pd.read_csv('figures/candidates/bulk-context-774-20261014T140904Z_22236311/proteomics/protein_de_adjusted.tsv',sep='\t',keep_default_na=False,na_values=[''])
pm=prot.merge(sym[sym.symbol.isin(u)][['symbol','ensembl_id','logFC']].rename(columns={'symbol':'gene','logFC':'bulk_lfc'}),on='gene')
pm=pm[['ensembl_id','logFC','padj','bulk_lfc']].rename(columns={'logFC':'prot_lfc','padj':'prot_q'})
g=g.merge(pm,on='ensembl_id',how='left')
g['prot_meas']=g.prot_lfc.notna()
g['prot_sig']=g.prot_q<0.05
g['prot_conc']=g.prot_sig & (np.sign(g.prot_lfc)==np.sign(g.bulk_lfc))
g['dir_call']=g.dir.isin(['call_up','call_down'])
rows=[]
for cl in ['genetically_anchored','concordant','established_state_associated','unresolved','untested']:
    x=g[g.primary_evidence_class==cl]
    rows.append(dict(cls=cl,n=len(x),gen_sup=int(x.gen_sup.sum()),dir_call=int(x.dir_call.sum()),
      stage_tested=int(x.stage_any_fib_vs_F0.notna().sum()),stage_any_fib=int(x.stage_any_fib_vs_F0.fillna(False).sum()),stage_any_adj=int(x.stage_any_adjacent.fillna(False).sum()),
      f4_sig=int((x.f4_fdr<0.05).sum()),prot_meas=int(x.prot_meas.sum()),prot_sig=int(x.prot_sig.sum()),prot_conc=int(x.prot_conc.sum())))
print(pd.DataFrame(rows).to_string())
x=g[g.gen_sup]
print('436 SuSiE: n',len(x),'dir',x.dir.value_counts().to_dict(),'stage any fib',int(x.stage_any_fib_vs_F0.fillna(False).sum()),'prot meas',int(x.prot_meas.sum()),'prot conc',int(x.prot_conc.sum()))
# risk-allele direction vs F4 stage direction among anchored with call and F4 sig
y=x[x.dir_call & (x.f4_fdr<0.05)]
conc=((y.dir=='call_up')&(y.f4_lfc>0))|((y.dir=='call_down')&(y.f4_lfc<0))
print('dir call & F4 sig:',len(y),'concordant',int(conc.sum()))
y2=x[x.dir_call & x.stage_any_fib_vs_F0.fillna(False)]
print('dir call & any fib stage sig:',len(y2))
g.to_csv('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/critic-cand-3/catalog_annotation_preview.tsv',sep='\t',index=False)
