# NG04 feasibility counts. Does NOT compute the primary estimand (sign agreement with liver RNA).
import pandas as pd, numpy as np, json, sys
from scipy import stats
from statsmodels.stats.power import TTestIndPower
import statsmodels.api as sm
np.random.seed(20260923)
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
O=R+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG04/'
m=pd.read_csv(R+'data/PXD051911/meta_data.txt',sep='\t'); p=pd.read_csv(R+'data/PXD051911/plasma_protein_quant.txt',sep='\t')
i=m[(m.sample_group=='initial_sample')&m.kleiner_fibrosis_grade.notna()&m.plasma_proteomics_filename.isin(p.columns)].copy()
out={'n_participants':len(i),'unique_patients':int(i.patient_name.nunique()),'saf':i.saf_diagnosis.value_counts().to_dict(),
     'sex':i.gender.value_counts().to_dict(),'age_na':int(i.alder.isna().sum()),'bmi_na':int(i.bmi.isna().sum())}
X=np.log2(p[i.plasma_proteomics_filename.tolist()].replace(0,np.nan).astype(float))
keep=(X.notna().mean(1)>=0.5)&p.Genes.notna()&~p.Genes.fillna('').str.contains(';')
X=X[keep]; genes=p.Genes[keep].values; X.index=genes
out['n_single_gene_det50']=int(keep.sum()); out['dup_genes']=int(pd.Series(genes).duplicated().sum())
g=pd.read_parquet(R+'RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/passport_gene_index.parquet')
gi=g.drop_duplicates('symbol').set_index('symbol').primary_evidence_class
cls=pd.Series(genes,index=genes).map(gi).fillna('not_in_catalog')
out['class_counts']=cls.value_counts().to_dict()
e=pd.read_parquet(R+'RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/passport_evidence_long.parquet')
t=e[e.evidence_domain=='transcriptomics'].drop_duplicates('symbol').set_index('symbol')
state=[s for s in genes if cls[s] in ('established_state_associated','concordant')]
st=t.loc[[s for s in state if s in t.index],['estimate','q_value','direction','call_state']]
out['n_state_direction']=len(state); out['state_rna_direction']=st.direction.value_counts().to_dict(); out['state_call']=st.call_state.value_counts().to_dict()
out['state_rna_abs_log2fc_median']=float(st.estimate.abs().median())
st.to_csv(O+'state_proteins_rna.tsv',sep='\t')
# detection per state protein
out['state_detect_frac_min']=float(X.loc[state].notna().mean(1).min())
# unsigned plasma signal strength, all 448, MASLD vs No_MASLD, age+sex adjusted OLS (direction relative to RNA NOT computed)
y=(i.saf_diagnosis!='No_MASLD').astype(int).values
cov=pd.DataFrame({'d':y,'age':i.alder.values,'sex':(i.gender=='Male').astype(int).values})
res=[]
for gname,row in X.iterrows():
    v=row.values; ok=~np.isnan(v)
    zz=(v[ok]-v[ok].mean())/v[ok].std(ddof=1)
    f=sm.OLS(zz,sm.add_constant(cov[ok].astype(float))).fit()
    res.append((gname,cls[gname],abs(f.params['d']),f.pvalues['d']))
r=pd.DataFrame(res,columns=['gene','cls','abs_d','p'])
out['all448_abs_d_median']=float(r.abs_d.median()); out['all448_frac_p05']=float((r.p<0.05).mean())
rs=r[r.gene.isin(state)]
out['state45_abs_d_median']=float(rs.abs_d.median()); out['state45_frac_p05']=float((rs.p<0.05).mean()); out['state45_n_p05']=int((rs.p<0.05).sum())
# expected agreement if every state protein were truly concordant with true |d| = observed |d| (upper bound) and with half that
se=np.sqrt(1/100+1/36)
for scale in (1.0,0.5):
    out[f'expected_agreement_if_all_concordant_dscale{scale}']=float(np.mean(stats.norm.cdf(scale*rs.abs_d.values/se)))
# MDE
pw=TTestIndPower()
out['mde_d_per_protein_a05']=pw.solve_power(nobs1=100,ratio=0.36,alpha=0.05,power=0.8)
out['mde_d_per_protein_a05_over56']=pw.solve_power(nobs1=100,ratio=0.36,alpha=0.05/56,power=0.8)
k=int(stats.binom.ppf(0.95,len(state),0.5))+1
out['binom_crit']=k
for pa in (0.6,0.65,0.7,0.75): out[f'binom_power_{pa}']=float(stats.binom.sf(k-1,len(state),pa))
# effective number of independent proteins among state set (eigenvalue-based, Li & Ji)
C=X.loc[state].T.corr(min_periods=30).fillna(0).values; ev=np.linalg.eigvalsh(C); ev=np.abs(ev)
out['state_Meff_LiJi']=float(np.sum((ev>=1)+(ev-np.floor(ev))))
# positive control IGF1
ig=X.loc['IGF1'].values; kl=i.kleiner_fibrosis_grade.str[1].astype(int).values; ok=~np.isnan(ig)
rho,pp=stats.spearmanr(ig[ok],kl[ok]); out['IGF1_detect']=float(ok.mean()); out['IGF1_spearman_kleiner']=float(rho); out['IGF1_p_two_sided']=float(pp)
out['IGF1_rna']=t.loc['IGF1',['estimate','direction']].astype(str).to_dict()
json.dump(out,open(O+'gate_counts.json','w'),indent=1,default=str); print(json.dumps(out,indent=1,default=str))
