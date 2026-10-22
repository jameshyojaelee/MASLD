import numpy as np, pandas as pd
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
S=pd.read_csv(root+'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv',sep='\t',usecols=['program_uid','donor','score','dataset'])
M=pd.read_csv(root+'figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/qualified_donor_metadata.tsv',sep='\t')
M=M[(M.exclude_stage_analysis.astype(str).str.upper()=='FALSE')&(M.stage_ordinal.isin([0,1,2]))]
import statsmodels.formula.api as smf
for name,u in [('Ductular','hotspot_hepatocytes_48f39dd4d817a10e'),('ECM','hotspot_hepatocytes_f05c535ae5bbc0b9')]:
    d=S[S.program_uid==u].merge(M[['donor','stage_ordinal']],on='donor')
    for lab,dd in [('all',d),('drop single-stage datasets',d[d.groupby('dataset').stage_ordinal.transform('nunique')>1])]:
        f=smf.ols('score ~ C(dataset) + stage_ordinal',dd).fit()
        h=f.get_robustcov_results('HC3')
        print(name,lab,'n',len(dd),'datasets',dd.dataset.nunique(),'beta %.6f'%f.params['stage_ordinal'],'p %.3g'%f.pvalues['stage_ordinal'],'HC3 p %.3g'%h.pvalues[list(f.params.index).index('stage_ordinal')])
    # per-dataset slopes
    for ds,g in d.groupby('dataset'):
        if g.stage_ordinal.nunique()>1 and len(g)>3:
            f=smf.ols('score ~ stage_ordinal',g).fit(); print('   ',ds,len(g),'slope %.3f p %.3g'%(f.params['stage_ordinal'],f.pvalues['stage_ordinal']))
print('--- LODO and interaction')
for name,u in [('Ductular','hotspot_hepatocytes_48f39dd4d817a10e'),('ECM','hotspot_hepatocytes_f05c535ae5bbc0b9')]:
    d=S[S.program_uid==u].merge(M[['donor','stage_ordinal']],on='donor')
    for ds in ['GSE202379','GSE244832','GSE189600']:
        dd=d[d.dataset!=ds]
        f=smf.ols('score ~ C(dataset) + stage_ordinal',dd).fit()
        print(name,'drop',ds,'n',len(dd),'beta %.3f p %.3g'%(f.params['stage_ordinal'],f.pvalues['stage_ordinal']), 'BH117-equivalent q if rank1: %.3g'%(f.pvalues['stage_ordinal']*117))
    dd=d[d.dataset.isin(['GSE202379','GSE244832'])]
    f0=smf.ols('score ~ C(dataset) + stage_ordinal',dd).fit(); f1=smf.ols('score ~ C(dataset) * stage_ordinal',dd).fit()
    import statsmodels.api as sm
    print(name,'stage x dataset interaction (GSE202379 vs GSE244832) F-test', sm.stats.anova_lm(f0,f1).iloc[1][['F','Pr(>F)']].to_dict())
    print(name,'stage counts', d.groupby(['dataset','stage_ordinal']).size().to_dict())
