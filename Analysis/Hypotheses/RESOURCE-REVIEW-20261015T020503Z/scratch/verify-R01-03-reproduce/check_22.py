import pandas as pd, numpy as np
w=pd.read_csv('Analysis/Multimodal_Program_Projection/results/program_context_wide.tsv',sep='\t')
p=pd.read_csv('figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/primary_effects.tsv',sep='\t')
p['program_id']=p.cell_type+'::'+p.module.astype(str)
print('primary rows',len(p),'estimable',p.estimable.sum(),'q<0.05',(p.qvalue<0.05).sum(), 'hc3 q<0.05',(p.hc3_qvalue<0.05).sum())
# re-derive BH within estimable
e=p[p.estimable].copy()
pv=e.pvalue.values; n=len(pv); o=np.argsort(pv); r=np.empty(n); 
q=pv[o]*n/np.arange(1,n+1); q=np.minimum.accumulate(q[::-1])[::-1]; q=np.minimum(q,1); qq=np.empty(n); qq[o]=q
e['q_rederived']=qq
print('max |q_rederived - qvalue|', np.nanmax(np.abs(e.q_rederived-e.qvalue)), 'family n', n)
sel=pd.read_csv('figures/main/fig5_molecular_context/panels/data/fig5e_multimodal_program_summary_selection.tsv',sep='\t')
m=w[['program_id','program_name','fig3_beta','fig3_q']].merge(p[['program_id','estimable','beta','n_donors','qvalue','hc3_qvalue']],on='program_id',how='left')
m['displayed_5E']=m.program_id.isin(sel.program_id)
pd.set_option('display.width',250)
print(m.to_string())
print('22 universe: corrected q<0.05 =', (m.qvalue<0.05).sum(), list(m.loc[m.qvalue<0.05,'program_name']))
print('displayed 13: corrected q<0.05 =', ((m.qvalue<0.05)&m.displayed_5E).sum(), 'displayed n', m.displayed_5E.sum())
print('displayed 13 min corrected q', m.loc[m.displayed_5E,'qvalue'].min())
print('sign flips beta vs fig3_beta among 22:', (np.sign(m.beta)!=np.sign(m.fig3_beta)).sum())
