"""H78: DepMap 24Q4 coverage counts (header and list sizes only; no route-level flag fractions)."""
import pandas as pd
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
D=R+'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf/source/depmap_24q4/'
W=R+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H78/'
cols=pd.read_csv(D+'CRISPRGeneEffect.csv',nrows=0).columns[1:]
dm=pd.DataFrame({'col':cols}); dm['sym']=dm.col.str.split(' \\(').str[0]; dm['entrez']=dm.col.str.extract(r'\((\d+)\)')[0].astype(int)
ce=pd.read_csv(W+'depmap24q4/CRISPRInferredCommonEssentials.csv'); ce['entrez']=ce.iloc[:,0].str.extract(r'\((\d+)\)')[0].astype(int)
ctl=pd.read_csv(W+'depmap24q4/AchillesCommonEssentialControls.csv'); ctl['entrez']=ctl.iloc[:,0].str.extract(r'\((\d+)\)')[0].astype(int)
print('24Q4 genes',len(dm),'unique entrez',dm.entrez.nunique(),'| inferred CE',len(ce),'in matrix',ce.entrez.isin(dm.entrez).sum(),'| CE controls',len(ctl),'overlap with inferred CE',ctl.entrez.isin(ce.entrez).sum())
g=pd.read_parquet(R+'figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'); g['ens']=g.ensembl_id.str.split('.').str[0]
gn=pd.read_csv(R+'data/external/gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv',sep='\t',usecols=['gene','gene_id','mane_select'],low_memory=False)
gn=gn[gn.mane_select.astype(str).str.lower().eq('true')]
e=gn[gn.gene_id.str.startswith('ENSG')][['gene','gene_id']]; z=gn[gn.gene_id.str.fullmatch(r'\d+')][['gene','gene_id']].rename(columns={'gene_id':'entrez'})
m=e.merge(z,on='gene'); m['entrez']=m.entrez.astype(int)
g=g.merge(m.rename(columns={'gene_id':'ens'})[['ens','entrez']],on='ens',how='left')
g['in24']=g.entrez.isin(dm.entrez)
t=g.groupby('next_experiment_rule_id').in24.agg(['size','sum']); t.loc['ALL']=[len(g),g.in24.sum()]
print(t.to_string())
print('inferred CE genes reachable from Catalog by Entrez',ce.entrez.isin(g.entrez).sum())
for nm in ['CEGv2','NEGv1']:
    h=pd.read_csv(R+f'data/public_r4/Hart2017_CEG2/{nm}.txt',sep='\t'); h['ENTREZ_ID']=h.ENTREZ_ID.astype('Int64')
    print(nm,'in 24Q4 matrix by Entrez',h.ENTREZ_ID.isin(dm.entrez).sum())
