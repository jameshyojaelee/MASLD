import pandas as pd, numpy as np
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
p=pd.read_csv(B+'Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv',sep='\t')
print(len(p), p.primary_signal_pair.value_counts().to_dict(), p.evidence_state.value_counts().to_dict())
p=p[p.primary_signal_pair.astype(str).str.upper()=='TRUE'] if p.primary_signal_pair.dtype==object else p[p.primary_signal_pair]
pairs=p[['gwas_name','ensembl','trait_class']].drop_duplicates()
print('pairs',len(pairs), pairs.trait_class.value_counts().to_dict())
genes=pairs.drop_duplicates(['ensembl','trait_class'])
print('genes by class',genes.trait_class.value_counts().to_dict(),'unique genes',pairs.ensembl.nunique(),
      'both',len(set(genes[genes.trait_class=='direct_MASLD'].ensembl)&set(genes[genes.trait_class=='liver_enzyme'].ensembl)))
# installed pair-level panel check
inst=pd.read_csv(B+'figures/main/fig5_molecular_context/panels/data/fig5d_variant_consistent_accessibility.tsv',sep='\t')
pc=p.groupby(['trait_class','lineage','evidence_state']).size().rename('n').reset_index()
mm=inst.merge(pc,on=['trait_class','lineage','evidence_state'],how='left').fillna({'n':0})
print('installed vs recomputed pair counts mismatches',(mm.n_pairs!=mm.n).sum())
p['mx']=p[['gse244832_accessible_mass','gse281367_accessible_mass']].max(axis=1)
def st(r):
    if r.evidence_state in('replicated_accessible','source_dependent'): return r.evidence_state
    if r.evidence_state=='indeterminate': return 'indet'
    return 'partial' if r.mx>=0.05 else 'under5'
p['s']=p.apply(st,axis=1)
rk={'replicated_accessible':0,'source_dependent':1,'partial':2,'under5':3,'indet':4}
p['rk']=p.s.map(rk)
# sanity: is evidence_state consistent with >=0.5 rule?
print(p.groupby('evidence_state').mx.agg(['min','max']))
g=p.groupby(['trait_class','lineage','ensembl']).rk.min().reset_index()
g['s']=g.rk.map({v:k for k,v in rk.items()})
g.loc[g.s=='indet','s']='under5'
tab=g.groupby(['trait_class','lineage','s']).size().unstack(fill_value=0)
print(tab)
rep=pd.read_csv(B+'figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/repair/fig5d_gene_level_display_counts.tsv',sep='\t')
rep['disp']=rep.disp.replace({'partial_ge_floor':'partial','under_5pct':'under5'})
mine=tab.stack().rename('mine').reset_index().rename(columns={'s':'disp'})
c=rep.merge(mine,on=['trait_class','lineage','disp'],how='outer').fillna(0)
print('display count mismatches', (c.N!=c.mine).sum(), len(c))
# genes under 5% in every lineage
ge=g.groupby(['trait_class','ensembl']).s.apply(lambda x:(x=='under5').all()).groupby('trait_class').agg(['sum','size'])
print('genes under5 every lineage\n',ge)
# pair-level partial rows and under 1%
part=p[p.evidence_state=='partial']
print('partial rows',len(part),'under1pct',(part.mx<0.01).sum())
pm=p.groupby(['gwas_name','ensembl','trait_class']).mx.max().reset_index()
print('pairs under 1% every lineage',(pm.mx<0.01).sum(), 'direct',((pm.mx<0.01)&(pm.trait_class=='direct_MASLD')).sum())
print('pairs under 5% every lineage',(pm.mx<0.05).sum())
print('no_overlap rows', (p.evidence_state=='no_overlap').sum(), 'partial rows with mx<0.05',(part.mx<0.05).sum())
