"""Focused independent numerical checks and existing Figure4C rendering."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from scipy import stats

base=Path(os.environ['MASLD_PROJECT_ROOT']); out=Path(os.environ['DIAGNOSIS_REFIT_ROOT'])
old=base/'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot'
checks=[]
for name in ['donor_program_scores_primary.tsv','donor_program_scores_equal_run.tsv']:
    source=pd.read_csv(old/name,sep='\t'); new=pd.read_csv(out/'hotspot'/name,sep='\t')
    cols=['program_uid','donor','score']
    expected=source[~source.donor.isin(['SRR14602057','SRR14602058'])][cols].sort_values(cols[:2]).reset_index(drop=True)
    actual=new[cols].sort_values(cols[:2]).reset_index(drop=True)
    assert expected.equals(actual),name
    checks.append(name+': scores and keys exactly unchanged except fixed-eligibility exclusions')
assert (old/'program_membership_v2.tsv').read_bytes()==(out/'hotspot/program_membership_v2.tsv').read_bytes()
checks.append('117 memberships and source weights byte-identical')
p=pd.read_csv(out/'hotspot/donor_program_scores_primary.tsv',sep='\t')
results=pd.read_csv(out/'hotspot/primary_effects.tsv',sep='\t')
errors=[]
for _,row in results.iterrows():
    d=p[p.cell_type.eq(row.cell_type)&p.module.eq(row.module)]
    if not row.estimable: continue
    x=np.column_stack([np.ones(len(d)),d.stage_ordinal,pd.get_dummies(d.dataset,drop_first=True,dtype=float)])
    y=d.score.to_numpy(); inv=np.linalg.inv(x.T@x); coef=inv@x.T@y; residual=y-x@coef; df=len(y)-np.linalg.matrix_rank(x)
    se=np.sqrt((residual@residual)/df*inv[1,1]); pv=2*stats.t.sf(abs(coef[1]/se),df)
    leverage=np.einsum('ij,jk,ik->i',x,inv,x)
    hc=inv@(x.T@((residual/(1-leverage))[:,None]**2*x))@inv
    hp=2*stats.t.sf(abs(coef[1]/np.sqrt(hc[1,1])),df)
    assert np.allclose([coef[1],se,pv,hp],[row.beta,row.se,row.pvalue,row.hc3_pvalue],rtol=1e-8,atol=1e-10)
checks.append('Independent NumPy OLS/classical/HC3 reconstruction passes all estimable programs')
reg=pd.read_csv(out/'hotspot/program_registry_v2.tsv',sep='\t')
assert len(reg)==117 and reg.primary_selected.sum()==1
checks.append('Complete117-family; exactly BICC1 primary selected')
(out/'verification.json').write_text(json.dumps(checks,indent=2)+'\n')
code=(base/'scripts/figures/singlecell_module_heatmap.R').read_text()
code=code[:code.index('draw_modules(W,')]
code += '\ndraw_modules(W, OUT_SUPP, block_titles=TRUE)\n'
code += '\ndraw_modules(W[rk<=3 & cell_type!="tcells"],OUT_MAIN,block_titles=FALSE,page_width_in=255/72,page_height_in=261/72,row_mm=2.0,overhead_in=.56)\n'
(out/'source/render_qualified_4c.R').write_text(code)
subprocess.run([os.environ['RSCRIPT'],'--vanilla',str(out/'source/render_qualified_4c.R')],check=True)
