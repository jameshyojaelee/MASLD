# Collapse COLOC r2 gene x trait direction rows to one call per gene (H26 rule without the GTEx and coding flags).
import pandas as pd
D='Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z/coloc_direction/full-20260930T120450Z/direction_by_trait.tsv'
t=pd.read_csv(D,sep='\t',low_memory=False,keep_default_na=False,na_values=[''])
r=t[t.release=='candidate_r2_20260924T053244Z']
def collapse(s):
    s=s.dropna(); s=s[s!=0]
    if len(s)==0: return 'untestable'
    return ('call_up' if s.iloc[0]>0 else 'call_down') if s.nunique()==1 else 'ambiguous'
gs=r.groupby('ensembl').direction_sign.apply(collapse)
print('rows',len(r),'genes',r.ensembl.nunique()); print(gs.value_counts().to_dict())
gs.rename('dir').to_csv('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/critic-cand-3/r2_gene_direction_collapse.tsv',sep='\t')
