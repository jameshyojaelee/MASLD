# Synthetic inputs for testing 05_fit_h90.R. Random scores only; no GSE289173 expression.
import numpy as np, pandas as pd, sys
T=sys.argv[1]; R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
rng=np.random.default_rng(1)
d=pd.read_csv(R+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gse289173/write/donor_design.tsv',sep='\t').sort_values('gsm')
reg=pd.read_csv(R+'figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv',sep='\t')[['cell_type','module']]
d['n_hep']=rng.integers(40,6000,len(d)); d['median_umi_hep']=rng.integers(800,3000,len(d))
d[['gsm','n_hep','median_umi_hep','title','group','stage','sex','series']].to_csv(T+'/donor_covariates_hep.tsv',sep='\t',index=False)
rows=[];ln=[]
for ct,md in reg.itertuples(index=False):
    keep = d.gsm if ct=='hepatocytes' else d.gsm[rng.random(len(d))>0.2]
    for g in d.gsm: ln.append(dict(cell_type=ct,gsm=g,n_nuclei=int(rng.integers(30,500)),eligible=True))
    for g in keep: rows.append(dict(cell_type=ct,module=md,gsm=g,score=rng.normal()+0.2*float(d.stage[d.gsm==g].iloc[0])*(md==20)))
pd.DataFrame(rows).to_csv(T+'/donor_scores.tsv',sep='\t',index=False)
ln=pd.DataFrame(ln).drop_duplicates(['cell_type','gsm']); ln.to_csv(T+'/lineage_donor_nuclei.tsv',sep='\t',index=False)
pd.DataFrame([dict(score='projection',dataset='n-weighted',n=55,m1=-0.2,m2=0.8)]).to_csv(T+'/D4_bridge.tsv',sep='\t',index=False)
pd.DataFrame([dict(cell_type='hepatocytes',module=8,n_scored=40,retained_weight=0.9)]).to_csv(T+'/D3_coverage.tsv',sep='\t',index=False)
pd.DataFrame(rng.normal(size=(1000,len(d))),columns=d.gsm).to_csv(T+'/NC3_set_scores.tsv.gz',sep='\t',index=False)
genes=[f'G{i}' for i in range(21)]
e6=pd.DataFrame(rng.poisson(50,(len(d),21)),columns=genes); e6.insert(0,'lib_size',rng.integers(1e6,1e7,len(d))); e6.insert(0,'gsm',d.gsm.values); e6.to_csv(T+'/E6_pseudobulk_counts.tsv',sep='\t',index=False)
pd.DataFrame(dict(cell_type='hepatocytes',module=20,gene=genes,weight=np.ones(21)/21)).to_csv(T+'/scored_genes_weights.tsv',sep='\t',index=False)
pd.DataFrame(dict(gsm=d.gsm,score=rng.normal(size=len(d)),n_hep_kept=d.n_hep-5)).to_csv(T+'/E10_scores.tsv',sep='\t',index=False)
d.assign(included=True).to_csv(T+'/donor_qc.tsv',sep='\t',index=False)
pd.DataFrame([dict(check='G4c_spearman_within_dataset (qualifier)',value=0.7,expected='',**{'pass':True})]).to_csv(T+'/gate_G4_bridge.tsv',sep='\t',index=False)
pd.Series({'x':1}).to_csv(T+'/gate_G1_G3.tsv',sep='\t',header=['value'])
