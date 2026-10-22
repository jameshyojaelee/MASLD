"""Covariate balance after 10x10x10 stratification (no outcome), observability-column counts, CMH MDE."""
import glob, numpy as np, pandas as pd
from scipy.stats import norm
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
S=R+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L03-01/'
t=pd.read_csv(S+'gate_gene_table.tsv',sep='\t',index_col=0)
t['ctrl']=~t.case
# balance: case minus stratum control mean, averaged over cases (strata with >=1 control)
for v in ['AveExpr','len','n_snps']:
    x=np.log10(t[v]) if v=='len' else t[v]
    cm=x[t.ctrl].groupby(t.stratum[t.ctrl]).mean()
    cases=t[t.case]; xc=x[t.case]
    diff=(xc-cases.stratum.map(cm)).dropna()
    sd=np.sqrt((x[t.case].var()+x[t.ctrl].var())/2)
    raw=(x[t.case].mean()-x[t.ctrl].mean())/sd
    print(f'{v}{" (log10)" if v=="len" else ""}: raw SMD {raw:.3f}; within-stratum SMD {diff.mean()/sd:.3f} over {len(diff)} cases')
lo=t[t.case & (t.stratum//100==0)]
print('cases in lowest AveExpr decile:',len(lo))
# CMH-style MDE with all controls, effective n per stratum (sum of harmonic terms)
st=t.groupby('stratum').agg(n1=('case','sum'),n0=('ctrl','sum'))
st=st[(st.n1>0)&(st.n0>0)]
p0=t[t.ctrl].groupby(t.stratum).apply(lambda s: s.egene_cur.mean()) if False else None
# outcome-blind: use overall control eGene rate within strata is outcome information on controls only
pc=t[t.ctrl].groupby('stratum').egene_cur.mean().reindex(st.index)
w=st.n1*st.n0/(st.n1+st.n0)
var=(w*pc*(1-pc)).sum()
z=norm.ppf(0.975)+norm.ppf(0.8)
d=z*np.sqrt(var)/w.sum()
pbar=(w*pc).sum()/w.sum()
print(f'CMH all controls: cases used {int(st.n1.sum())}, controls used {int(st.n0.sum())}, weighted control eGene rate {pbar:.3f}, MDE abs diff {d:.4f}, OR {((pbar-d)/(1-pbar+d))/(pbar/(1-pbar)):.3f} / {((pbar+d)/(1-pbar-d))/(pbar/(1-pbar)):.3f}')
print('control eGene rate by AveExpr decile:',t[t.ctrl].groupby(t.stratum//100).egene_cur.mean().round(3).to_dict())
# Observability column counts over the whole catalog (rule as written, current fit set)
P=R+'figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
e=pd.read_parquet(P+'passport_evidence_long.parquet')
e['ens']=e.ensembl_id.astype(str).str.split('.').str[0]
tr=e[e.evidence_domain=='transcriptomics'].set_index('ens').call_state
ge=e[e.evidence_domain=='genetics'].set_index('ens').testability_reason
def load(d):
    x=pd.concat([pd.read_csv(f,sep='\t') for f in sorted(glob.glob(R+f'GWAS/finemapping/results/{d}/chr*_summary.tsv'))])
    x['ens']=x.ensg.astype(str).str.split('.').str[0]; return x.drop_duplicates('ens',keep='last').set_index('ens')
cur=load('eqtl_susie_polyfun'); rf=load('eqtl_susie_refit')
cur.loc[rf.index.intersection(cur.index),rf.columns.intersection(cur.columns)]=rf.loc[rf.index.intersection(cur.index),rf.columns.intersection(cur.columns)]
genes=pd.Index(sorted(set(tr.index)|set(ge.index)))
def col(g):
    if tr.get(g)=='untestable': return 'not_detected_bulk'
    if ge.get(g)=='underpowered_source': return 'not_in_eqtl_panel'
    if g not in cur.index: return 'testable_but_absent_from_fit_table'
    r=cur.loc[g]
    if r.status!='ok': return 'fit_failed_'+str(r.status)
    if str(r.converged)!='True': return 'not_converged(n_cs>0)' if r.n_cs>0 else 'not_converged(n_cs=0)'
    return 'egene' if r.n_cs>0 else 'in_panel_no_credible_set'
c=pd.Series({g:col(g) for g in genes})
print('catalog genes',len(genes)); print(c.value_counts().to_string())
