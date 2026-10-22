import pandas as pd, numpy as np
CR='RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10'
ann=pd.read_csv(CR+'/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz',sep='\t',usecols=['gene_id','gene_name'])
ann['key']=ann.gene_id.str.replace(r'\..*$','',regex=True); ann['sym']=ann.gene_name.str.strip().str.upper()
ann=ann.dropna(subset=['sym']).drop_duplicates(['key','sym'])
mem=pd.read_csv('Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv',sep='\t')
comm=pd.read_csv(CR+'/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction/program_communities.tsv',sep='\t')
feats=comm.feature_id.tolist()
mem['sym']=mem.mapped_symbol.astype(str).str.strip().str.upper()
mem=mem[mem.mapped_symbol.notna()&(mem.sym!='')&mem.program_uid.isin(feats)]
w=mem.groupby(['program_uid','sym'],as_index=False).original_l1_weight.sum()
print('features',len(feats),'panel symbols',w.sym.nunique(),'program rows',len(w))
def cov(counts, cols, label):
    X=counts[cols]
    lib=X.sum(0)
    lc=np.log2((X+0.5*lib/lib.mean()*0+1)/ (lib+2)*1e6)  # approx edgeR prior.count=1
    keys=X.index.str.replace(r'\..*$','',regex=True)
    sym=ann.set_index('key').sym.reindex(keys)
    ok=sym.notna().values
    lc=lc[ok]; lc.index=sym[ok].values
    S=lc.groupby(level=0).mean()
    sd=S.std(axis=1)
    var=set(S.index[sd>0])
    tot=w.groupby('program_uid').original_l1_weight.sum()
    ret=w[w.sym.isin(var)].groupby('program_uid').original_l1_weight.sum().reindex(tot.index).fillna(0)
    frac=ret/tot
    sym_cov=w.sym[w.sym.isin(var)].nunique()
    print(f"{label}: n={len(cols)} panel symbols variable {sym_cov}/{w.sym.nunique()}; programs testable (>=0.80 retained weight) {(frac>=0.8).sum()}/{len(feats)}; min frac {frac.min():.3f}; programs <0.95: {(frac<0.95).sum()}")
    return frac
a=pd.read_csv('data/GSE281797/metadata/sample_phenotypes_resolved.tsv',sep='\t')
m=pd.read_csv('data/GSE281797/GSE281797_counts_ensembl_stable.tsv.gz',sep='\t',index_col=0)
c91=a[a.sex_concordant==True]
cov(m,a.matrix_column.tolist(),'GSE281797 all94')
cov(m,c91.matrix_column.tolist(),'GSE281797 n91')
cov(m,c91[c91.fibrosis<=1].matrix_column.tolist(),'GSE281797 F0-F1 n91')
r=pd.read_csv('data/GSE225740/metadata/sample_phenotypes_resolved.tsv',sep='\t')
m2=pd.read_csv('data/GSE225740/GSE225740_counts_ensembl_stable.tsv.gz',sep='\t',index_col=0)
d=r[r.fibrosis.notna()&r.nas.notna()]
cov(m2,d.matrix_column.astype(str).tolist(),'GSE225740 87')
f2=cov(m2,d[d.fibrosis<=1].matrix_column.astype(str).tolist(),'GSE225740 F0-F1')
print(f2.sort_values().head(8).round(3).to_string())
