import pandas as pd, numpy as np
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
F=R+'figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/'
D=R+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F4abcd/run-20261015T_recheck/'
H=R+'Analysis/Hypotheses/H59-20261015T112432Z/out/'
fam=pd.read_csv(H+'outcome/results_primary_family.tsv',sep='\t')
nt=pd.read_csv(H+'null_gls_t.tsv.gz',sep='\t')
crit=nt.abs().quantile(1-0.05/117)  # pandas linear = R type 7
cr=pd.read_csv(H+'null_critical_abs_t.tsv',sep='\t').set_index('program').null_abs_t_crit
print('max|crit diff|', (crit-cr.reindex(crit.index)).abs().max(), 'hep20', cr['hepatocytes_20'], 'hep8', cr['hepatocytes_8'])
reg=pd.read_csv(F+'repair/source_tables/fig4c_h59_registry_117_with_familywise_interval.tsv',sep='\t')
reg=reg.merge(fam[['program','gls_beta','gls_se','gls_ci_low','gls_ci_high','perm_q']],on='program')
reg['my_fw_lo']=reg.gls_beta-cr.reindex(reg.program).values*reg.gls_se
reg['my_fw_hi']=reg.gls_beta+cr.reindex(reg.program).values*reg.gls_se
print('fw diff', (reg.my_fw_lo-reg.fw_lo).abs().max(), (reg.my_fw_hi-reg.fw_hi).abs().max(), 'beta diff', (reg.beta-reg.gls_beta).abs().max(), 'q diff',(reg.q-reg.perm_q).abs().max())
print('117: nominal excl0', ((reg.gls_ci_low>0)|(reg.gls_ci_high<0)).sum(), 'fw excl0', ((reg.my_fw_lo>0)|(reg.my_fw_hi<0)).sum())
# 24-row selection: non-T, per lineage 3 up + 3 down smallest perm q (tie-break perm_p then |t|?)
sel=[]
for ct,g in fam[fam.cell_type!='tcells'].groupby('cell_type'):
    for sgn in ['up','down']:
        gg=g[g.gls_beta>0] if sgn=='up' else g[g.gls_beta<0]
        gg=gg.sort_values(['perm_q','perm_p'])
        sel.append(gg.head(3).assign(dir=sgn))
        print(ct,sgn,[ (m,round(q,4),round(p,4)) for m,q,p in zip(gg.module,gg.perm_q,gg.perm_p)][:5])
sel=pd.concat(sel)
m24=pd.read_csv(F+'source_tables/fig4c_h59_main_24rows.tsv',sep='\t')
print('24 rows match set:', set(zip(sel.cell_type,sel.module))==set(zip(m24.cell_type,m24.module)), len(m24))
print(sorted(set(zip(sel.cell_type,sel.module))^set(zip(m24.cell_type,m24.module))))
m24=m24.merge(reg[['cell_type','module','my_fw_lo','my_fw_hi']],on=['cell_type','module'])
print('24: nominal excl0', ((m24.lo>0)|(m24.hi<0)).sum(), 'fw excl0', ((m24.my_fw_lo>0)|(m24.my_fw_hi<0)).sum(), 'selected', m24.selected_H59.astype(str).str.lower().eq('true').sum())
hm=pd.read_csv(D+'heatmap_medians_refit.tsv',sep='\t')
x=m24.merge(hm,on=['cell_type','module'],suffixes=('','_my'))
print('heatmap max diff', max((x[s]-x[s+'_my']).abs().max() for s in ['Healthy','Steatosis','Steatohepatitis']))
print(m24[['cell_type','module','name','beta','lo','hi','q']].to_string())
