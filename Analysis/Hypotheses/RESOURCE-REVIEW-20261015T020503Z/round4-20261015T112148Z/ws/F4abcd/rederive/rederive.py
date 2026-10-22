"""F4abcd login-node re-derivations from saved tables (no new model families)."""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
B = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
OUT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F4abcd/rederive')
S = B/'figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014'
H = B/'figures/candidates/source-diagnosis-programs-20261014T161500Z_22237556/hotspot'
STG = {'Healthy':0,'Steatosis':1,'Steatohepatitis':2}

pe = pd.read_csv(S/'primary_effects.tsv', sep='\t')
print('primary_effects cols', pe.columns.tolist()[:30])
sc = pd.read_csv(H/'donor_program_scores_primary.tsv', sep='\t')
sc = sc[(sc.exclude_stage_analysis.astype(str).str.upper()=='FALSE') & sc.disease_stage_coarse.isin(STG)].copy()
print('donors 62?', sc.donor.nunique(), 'datasets', sorted(sc.dataset.unique()))

def ols(d):
    X = np.column_stack([np.ones(len(d)), d.stage_ordinal.values.astype(float), pd.get_dummies(d.dataset, drop_first=True, dtype=float).values])
    y = d.score.values; inv = np.linalg.pinv(X.T@X); b = inv@X.T@y; r = y-X@b; df = len(y)-np.linalg.matrix_rank(X)
    se = np.sqrt(r@r/df*inv[1,1]); return b[1], se, 2*stats.t.sf(abs(b[1]/se), df), df

# R13-04 / R04-01: family counts
q = pe.qvalue if 'qvalue' in pe else pe.primary_qvalue
print('R13-04: q<0.05', int((pe.qvalue<0.05).sum()), 'hc3 q<0.05', int((pe.hc3_qvalue<0.05).sum()), 'min hc3 q', pe.hc3_qvalue.min())
r20 = pe[(pe.cell_type=='hepatocytes')&(pe.module==20)].iloc[0]; r8 = pe[(pe.cell_type=='hepatocytes')&(pe.module==8)].iloc[0]
print('Ductular', r20.beta, r20.pvalue, r20.qvalue, r20.hc3_qvalue, '| ECM', r8.beta, r8.qvalue, r8.hc3_qvalue)

# R04-02: within-dataset slopes, refit and LODO family
for mod,lab in [(20,'Ductular'),(8,'ECM')]:
    d = sc[(sc.cell_type=='hepatocytes')&(sc.module==mod)]
    b,se,p,df = ols(d); print(lab,'62 OLS', round(b,6), round(se,6), p)
    for ds in ['GSE202379','GSE244832']:
        x = d[d.dataset==ds]; s = stats.linregress(x.stage_ordinal, x.score); print('  ',ds,'n',len(x),'slope',round(s.slope,4),'se',round(s.stderr,4),'p',round(s.pvalue,4))
fam = []
for drop in [None,'GSE189600','GSE244832','GSE202379']:
    res = []
    for (ct,m),d in sc.groupby(['cell_type','module']):
        dd = d if drop is None else d[d.dataset!=drop]
        b,se,p,df = ols(dd); res.append((ct,m,b,p,len(dd)))
    R = pd.DataFrame(res, columns=['ct','m','b','p','n']); R['q'] = stats.false_discovery_control(R.p) if hasattr(stats,'false_discovery_control') else np.nan
    h = R[(R.ct=='hepatocytes')&(R.m==20)].iloc[0]
    print('drop',drop,'calls',int((R.q<0.05).sum()),'Ductular b',round(h.b,6),'q',round(h.q,6),'n',h.n,'minq',round(R.q.min(),4))
    fam.append(dict(drop=drop, calls=int((R.q<0.05).sum()), ductular_beta=h.b, ductular_q=h.q, n=h.n, min_q=R.q.min()))
pd.DataFrame(fam).to_csv(OUT/'R04-01_02_lodo_family.tsv', sep='\t', index=False)

# R04-03: informative donors by lineage (datasets with >1 stage)
inf = []
for ct,d in sc.drop_duplicates(['cell_type','donor']).groupby('cell_type'):
    nst = d.groupby('dataset').stage_ordinal.nunique(); infds = nst[nst>1].index
    inf.append(dict(cell_type=ct, donors=len(d), datasets=d.dataset.nunique(), informative_donors=int(d.dataset.isin(infds).sum()), informative_datasets=len(infds), informative_list=';'.join(infds),
                    steatosis_in_informative=int(((d.stage_ordinal==1)&d.dataset.isin(infds)).sum())))
inf = pd.DataFrame(inf); print(inf.to_string(index=False)); inf.to_csv(OUT/'R04-03_informative_donors.tsv', sep='\t', index=False)
design = sc.drop_duplicates(['cell_type','donor']).groupby(['cell_type','dataset','disease_stage_coarse']).size().unstack(fill_value=0)
print(design); design.to_csv(OUT/'R04-03_design_by_lineage_dataset_stage.tsv', sep='\t')

# R04-04: ambient two-program BH
am = pd.read_csv(S/'ambient_program_effects.tsv', sep='\t')
f = am[(am.cell_type=='hepatocytes')&am.module.isin([8,20])&(am.analysis_universe=='complete_case_common_universe')]
print(f[['module','n_model_donors','raw_beta','corrected_beta','delta_beta','delta_pvalue','delta_hc3_pvalue','delta_qvalue','raw_qvalue','corrected_qvalue','evidence_state']].to_string())
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); r=np.empty(n); r[o]=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]; return np.minimum(r,1)
print('two-program BH delta p', dict(zip(f.module, bh(f.delta_pvalue))), 'HC3', dict(zip(f.module, bh(f.delta_hc3_pvalue))))
print('pct change ECM', (f[f.module==8].corrected_beta.iloc[0]/f[f.module==8].raw_beta.iloc[0]-1)*100)
cu = am[am.analysis_universe=='complete_case_common_universe']
print('rank of ECM delta p among 117:', int((cu.delta_pvalue<=f[f.module==8].delta_pvalue.iloc[0]).sum()), 'of', cu.delta_pvalue.notna().sum())
f.to_csv(OUT/'R04-04_focal_ambient_rows.tsv', sep='\t', index=False)
