# Discovery-side effect sizes for the GSE289173 feasibility gate (no GSE289173 expression read).
import pandas as pd, numpy as np, gzip
np.random.seed(20260923)
R='Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/'
G='Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gse289173/gate/'
s=pd.read_csv(R+'donor_program_scores_primary.tsv',sep='\t',usecols=['cell_type','module','donor','score'])
m=pd.read_csv('Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv',sep='\t')
m=m[(~m.exclude_stage_analysis.astype(str).str.upper().eq('TRUE'))&m.disease_stage_coarse.isin(['Healthy','Steatosis','Steatohepatitis'])]
m=m[m.dataset!='GSE189600']
d=s.merge(m[['sample','dataset','disease_stage_coarse']],left_on='donor',right_on='sample')
d['stage']=d.disease_stage_coarse.map({'Healthy':0,'Steatosis':1,'Steatohepatitis':2})
print('donors',d.donor.nunique(), d.drop_duplicates('donor').groupby(['dataset','stage']).size().to_dict())
rows=[]
for (ct,mod),g in d.groupby(['cell_type','module']):
  for ds,h in g.groupby('dataset'):
    if h.stage.nunique()<2: continue
    def cd(a,b):
      if len(a)<2 or len(b)<2: return np.nan
      sp=np.sqrt(((len(a)-1)*a.var()+(len(b)-1)*b.var())/(len(a)+len(b)-2)); return (b.mean()-a.mean())/sp
    rows.append(dict(cell_type=ct,module=mod,dataset=ds,
      d_dis_vs_ctrl=cd(h.score[h.stage==0],h.score[h.stage>0]),
      d_mash_vs_masl=cd(h.score[h.stage==1],h.score[h.stage==2]),
      n0=(h.stage==0).sum(),n1=(h.stage==1).sum(),n2=(h.stage==2).sum()))
e=pd.DataFrame(rows); e.to_csv(G+'discovery_standardized_effects.tsv',sep='\t',index=False)
for ct,mod,name in [('hepatocytes',20,'Ductular injury'),('hepatocytes',8,'Stromal ECM')]:
  print(name); print(e[(e.cell_type==ct)&(e.module==mod)].round(3).to_string(index=False))
# median |d| across 117 programs, n-weighted average of two datasets
e['w']=e.n0+e.n1+e.n2
a=e.groupby(['cell_type','module']).apply(lambda x: pd.Series({'d_dc':np.nansum(x.d_dis_vs_ctrl*x.w)/x.w[x.d_dis_vs_ctrl.notna()].sum(),'d_mm':np.nansum(x.d_mash_vs_masl*x.w)/x.w[x.d_mash_vs_masl.notna()].sum()}))
print('117 programs |d| disease-vs-ctrl quantiles', a.d_dc.abs().quantile([.5,.75,.9,1]).round(2).to_dict())
print('117 programs |d| MASH-vs-MASL quantiles', a.d_mm.abs().quantile([.5,.75,.9,1]).round(2).to_dict())
a.to_csv(G+'discovery_pooled_d.tsv',sep='\t')
# gene coverage in GSE289173 10x 2020-A features
feat=set(l.split('\t')[1] for l in gzip.open('data/public_r4/GSE289173/features/GSM8785686_FL514_features.tsv.gz','rt'))
mem=pd.read_csv(R+'program_membership_v2.tsv',sep='\t')
cov=mem.assign(present=mem.mapped_symbol.isin(feat)).groupby(['cell_type','module']).apply(lambda x: pd.Series({'n':len(x),'present':x.present.sum(),'weight_present':x.original_l1_weight[x.present].sum()}))
cov.to_csv(G+'program_gene_coverage_10x2020A.tsv',sep='\t')
print('weight coverage min/median', cov.weight_present.min().round(3), cov.weight_present.median().round(3))
print(cov.loc[[('hepatocytes',20),('hepatocytes',8)]])
