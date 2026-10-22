# Independent refit of the Fig 4H x-axis (feature-level fibrosis-stage trend -> fixed-effect meta -> system aggregation)
# under three rosters. Read-only inputs. Usage: step2_refit_feature_stage.py OUTDIR [collections comma list]
import sys, os, numpy as np, pandas as pd
from scipy.stats import spearmanr, norm
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT=sys.argv[1]; os.makedirs(OUT, exist_ok=True)
COLLS=sys.argv[2].split(',') if len(sys.argv)>2 else ['hallmark','kegg','reactome','go_bp','go_mf','go_cc']
ML=ROOT+'RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261011T152106Z/'
MAN=ROOT+'figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv'
CAND=ROOT+'figures/candidates/continuum-774-panels-20261011T152106Z/'
EVAL=['GSE162694','GSE213621']
man=pd.read_csv(MAN,sep='\t')
man=man[man.dataset.isin(EVAL)][['sample_id','dataset','inferred_sex','fibrosis_stage','group_binary']]
print(pd.crosstab([man.dataset,man.group_binary],man.fibrosis_stage.fillna(-1)).to_string())
nodes=pd.read_csv(CAND+'molecular_systems_atlas/overlays/all_node_stage_continuum_overlays.tsv.gz',sep='\t')

def standardize(x):
    x=np.asarray(x,float); ok=np.isfinite(x); out=np.full(len(x),np.nan)
    if ok.sum()>=2 and np.std(x[ok],ddof=1)>0: out[ok]=(x[ok]-x[ok].mean())/np.std(x[ok],ddof=1)
    return out

def fit_cohort(wide, meta):
    """wide: samples x features outcome matrix (rows aligned with meta). OLS outcome ~ stage_z + sex per feature."""
    st=standardize(meta.fibrosis_stage.values)
    sex=(meta.inferred_sex.values=='M').astype(float)
    X=np.column_stack([np.ones(len(st)),st,sex])
    Y=wide.values
    res=[]
    if np.isfinite(Y).all():
        XtXi=np.linalg.inv(X.T@X); B=XtXi@X.T@Y; R=Y-X@B
        df=len(st)-3; s2=(R**2).sum(0)/df; se=np.sqrt(s2*XtXi[1,1])
        return pd.DataFrame({'feature_id':wide.columns,'beta':B[1],'se':se,'n':len(st)})
    for j,f in enumerate(wide.columns):
        y=Y[:,j]; ok=np.isfinite(y)
        Xo=X[ok]; XtXi=np.linalg.inv(Xo.T@Xo); b=XtXi@Xo.T@y[ok]; r=y[ok]-Xo@b
        se=np.sqrt((r**2).sum()/(ok.sum()-3)*XtXi[1,1]); res.append((f,b[1],se,ok.sum()))
    return pd.DataFrame(res,columns=['feature_id','beta','se','n'])

def meta(a,b):
    m=a.merge(b,on='feature_id',suffixes=('_a','_b'))
    wa=1/m.se_a**2; wb=1/m.se_b**2
    m['beta']=(wa*m.beta_a+wb*m.beta_b)/(wa+wb); m['se']=np.sqrt(1/(wa+wb))
    m['p']=2*norm.sf(np.abs(m.beta/m.se)); m['conc']=np.sign(m.beta_a)==np.sign(m.beta_b)
    return m

def run_layer(long, collection, restandardize_from_published=True):
    # long: feature_id, sample_id, dataset, outcome_z (already per published rule), fibrosis_stage, inferred_sex, group_binary
    out={}
    for roster in ['published','disease_only','disease_only_restd']:
        fits=[]
        for ds in EVAL:
            d=long[(long.dataset==ds)&long.fibrosis_stage.notna()&long.inferred_sex.notna()&np.isfinite(long.outcome_z)]
            if roster!='published' and ds=='GSE213621': d=d[d.group_binary=='Disease']
            d=d.copy()
            if roster=='disease_only_restd' and ds=='GSE213621':
                d['outcome_z']=d.groupby('feature_id').outcome_z.transform(lambda v:(v-v.mean())/v.std(ddof=1))
            wide=d.pivot(index='sample_id',columns='feature_id',values='outcome_z')
            md=d.drop_duplicates('sample_id').set_index('sample_id').loc[wide.index]
            f=fit_cohort(wide,md); f['dataset']=ds; fits.append(f)
        m=meta(fits[0],fits[1]); m['collection']=collection; m['roster']=roster
        m['n_GSE162694']=fits[0].n.iloc[0]; m['n_GSE213621']=fits[1].n.iloc[0]
        out[roster]=m
    return pd.concat(out.values())

allres=[]
# Hotspot: outcome_z as stored; stage and sex from the manifest (checked equal to stored)
hs=pd.read_csv(ML+'programs/hotspot/hotspot_participant_scores.tsv.gz',sep='\t')
hs=hs[hs.dataset.isin(EVAL)].rename(columns={'program_uid':'feature_id'})
chk=hs[['sample_id','fibrosis_stage','inferred_sex']].drop_duplicates().merge(man,on='sample_id',suffixes=('_h','_m'))
print('hotspot stage/sex agree with manifest:',((chk.fibrosis_stage_h.fillna(-9)==chk.fibrosis_stage_m.fillna(-9))&(chk.inferred_sex_h==chk.inferred_sex_m)).all(), len(chk))
hs=hs[['feature_id','sample_id','outcome_z']].merge(man,on='sample_id')
allres.append(run_layer(hs,'hotspot'))
for coll in COLLS:
    keep=set(nodes.feature_id[nodes.collection==coll])
    sc=pd.read_csv(ML+f'pathway_tf/pathways/{coll}/donor_scores.tsv.gz',sep='\t')
    sc=sc[sc.dataset.isin(EVAL)&sc.set_id.isin(keep)].rename(columns={'set_id':'feature_id'})
    sc['outcome_z']=sc.groupby(['dataset','feature_id']).pathway_score.transform(lambda v:standardize(v.values))
    sc=sc[['feature_id','sample_id','outcome_z']].merge(man,on='sample_id')
    allres.append(run_layer(sc,coll)); print(coll,'done',len(keep),flush=True)
res=pd.concat(allres)
res.to_csv(OUT+'/feature_stage_meta_three_rosters.tsv.gz',sep='\t',index=False)

# Check reproduction of published feature-level stage beta
pub=nodes[['feature_id','collection','community_id','stage_beta','continuum_beta']]
p=res[res.roster=='published'].merge(pub,on=['feature_id','collection'])
print('feature reproduction: n',len(p),'max|diff| stage beta',np.abs(p.beta-p.stage_beta).max())

def agg(df,col):
    d=df[df.community_id.notna()]
    return d.groupby(['community_id','collection'])[col].median().dropna().groupby('community_id').mean()
cont=agg(pub,'continuum_beta')
rows=[]
for roster in ['published','disease_only','disease_only_restd']:
    r=res[res.roster==roster].merge(pub[['feature_id','collection','community_id','continuum_beta']],on=['feature_id','collection'])
    s=agg(r,'beta'); s=s.loc[cont.index]
    rows.append(dict(level='system_43',roster=roster,n_GSE162694=int(r.n_GSE162694.iloc[0]),n_GSE213621=int(r.n_GSE213621.iloc[0]),
        spearman_rho=spearmanr(s,cont)[0],same_direction=int((np.sign(s)==np.sign(cont)).sum()),n=len(s),median_stage_beta=float(s.median())))
    pd.DataFrame({'community_id':s.index,'stage_beta':s.values,'continuum_beta_published':cont.values}).to_csv(OUT+f'/system_stage_{roster}.tsv',sep='\t',index=False)
    # Hotspot program-level (S4J-type) comparison with the published continuum beta
    hh=r[r.collection=='hotspot']
    hpub=pd.read_csv(CAND+'stage_continuum_comparison/source_tables/hotspot_stage_vs_continuum.tsv',sep='\t').rename(columns={'program_uid':'feature_id'})
    hh=hpub[['feature_id','continuum_beta']].merge(res[(res.roster==roster)&(res.collection=='hotspot')][['feature_id','beta']],on='feature_id').dropna()
    rows.append(dict(level='hotspot_programs',roster=roster,spearman_rho=spearmanr(hh.beta,hh.continuum_beta)[0],same_direction=int((np.sign(hh.beta)==np.sign(hh.continuum_beta)).sum()),n=len(hh)))
summ=pd.DataFrame(rows); summ.to_csv(OUT+'/summary.tsv',sep='\t',index=False)
pd.set_option('display.width',250); print(summ.to_string())
