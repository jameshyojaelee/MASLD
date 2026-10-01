"""Shared B-FIX P1/P2/P3 under v1.7, exactly fifty supplied draws per test."""
import json
import os
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from bfix_v1 import ROOT,FIXTURE,load,table
from likelihood_v1 import Model,conditional_draw
from tests_v1 import mc_result,p2_kappa,p2_tau


def serial(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


output=Path(os.environ['CODEX_REC_OUTPUT'])/('bfix_tests_'+os.environ['SLURM_JOB_ID'])
output.mkdir(exist_ok=False)
def save(name,value):
    (output/(name+'.json')).write_text(json.dumps(value,default=serial,indent=2)+'\n')
    print(name,'saved',flush=True)

design,settings=load();model=design.model();B=settings['B'];assert B==50
prior=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/bfix_fits_21984981'
old_null=json.loads((prior/'null.json').read_text())['fit']
assert old_null['converged']
start=time.monotonic()
null=model.fit(np.asarray(old_null['parameters']))
save('observed_null',null);assert null['converged']
null_p=np.asarray(null['parameters'])
marginal=model.fit_marginal(settings['tau_max_fixture'],null_p,null_fit=null)
save('observed_marginal',marginal);assert marginal['converged']
se=model.sandwich(null_p,design.clusters)
p1_exceed=[];p1_failed=[];p1_z=[]
if np.isfinite(se) and se>0:
    observed_z=float(null_p[2*model.G]/se)
    stage_rows=table('draws/p1_Sstar.tsv');assert len(stage_rows)==B
    for b,row in enumerate(stage_rows):
        assert int(row['draw'])==b
        sample=design.with_stage([int(row[p.individual]) for p in design.participants]);fit=None
        try:
            fitted_model=sample.model();fit=fitted_model.fit(null_p)
            sd=fitted_model.sandwich(fit['parameters'],sample.clusters)
            z=fit['parameters'][2*model.G]/sd
            okay=fit['converged'] and np.isfinite(z)
        except (ValueError,FloatingPointError,np.linalg.LinAlgError):okay=False
        if not okay:z=float('nan');p1_failed.append(b)
        p1_z.append(float(z));p1_exceed.append(not okay or abs(z)>=abs(observed_z)-1e-6)
        save(f'p1_{b:03d}',{'b':b,'fit':fit,'z':z,'exceed':p1_exceed[-1]})
    p1={**mc_result(p1_exceed),'computed':True,'observed_z':observed_z,'observed_se':se,
        'draw_z':p1_z,'failed':p1_failed}
else:p1={'computed':False,'reason':'observed cluster sandwich SE invalid','se':se}
save('p1_result',p1)

source_names=('GSE193066','PXD051911_A1')
lambda_table=table('draws/lambda_draws.tsv')
lambda_map={(int(r['draw']),r['source']):float(r['lambda']) for r in lambda_table}
lambda_draws=np.asarray([[lambda_map[b,s] for s in source_names] for b in range(B)])
groups=design.cohort_clusters()
group_index={c:{design.clusters[members[0]]:j for j,members in enumerate(v)} for c,v in groups.items()}
indices={b:{c:[] for c in groups} for b in range(B)}
for r in table('draws/p2_resample.tsv'):
    b=int(r['draw']);c=r['stratum']
    assert int(r['slot'])==len(indices[b][c])
    indices[b][c].append(group_index[c][r['cluster_id']])
kappas=[];p2_failed=[];genes={r.gene for r in design.rows}
for b in range(B):
    sample=design.resample(indices[b]);fit=None
    try:
        fit=sample.model().fit(null_p);okay=fit['converged'];k=fit['parameters'][2*model.G]
    except (ValueError,FloatingPointError,np.linalg.LinAlgError):okay=False
    if not okay:k=float('nan');p2_failed.append(b)
    kappas.append(float(k))
    save(f'p2_{b:03d}',{'b':b,'fit':fit,'kappa':k,'lambda':lambda_draws[b]})
kappa_result=p2_kappa(kappas,lambda_draws,settings['Delta'])
def profile_progress(event,payload):save('p2_'+event+'_'+str(payload['source']),payload)
tau_result=p2_tau(model,marginal,settings['tau_max_fixture'],lambda_draws,profile_progress)
p2={'computed':True,'p':max(kappa_result['p_kappa'],tau_result['p_tau']),
    'kappa':kappa_result,'tau':tau_result,'kappa_draws':kappas,'failed':p2_failed,'sources':source_names}
save('p2_result',p2)

ordered=[(r['individual_id'],r['gene_id']) for r in table('rows.tsv')]
assert ordered==[(r.individual,settings['gene_order'][r.gene]) for r in model.rows]
uniforms=np.load(FIXTURE/'draws/p3_uniforms.npy');assert uniforms.shape==(B,len(model.rows),2)
observed_lr=0. if marginal['v']==0. else max(0.,2*(marginal['loglik']-null['loglik']))
p3_exceed=[];p3_failed=[];p3_lr=[]
for b,draw in enumerate(uniforms):
    fit=None
    try:
        rows=[]
        for row,(u1,u2) in zip(model.rows,draw):
            eta,r,_,_=model.predictors(row,null_p)
            rows.append(conditional_draw(row,eta,r,u1,u2))
        fit=Model(rows).fit_marginal(settings['tau_max_fixture'],null_p,section4_from_null=True)
        okay=fit['converged'] and fit['null']['converged']
        lr=0. if fit['v']==0. else max(0.,2*(fit['loglik']-fit['null']['loglik']))
    except (ValueError,FloatingPointError,np.linalg.LinAlgError):okay=False
    if not okay:lr=float('nan');p3_failed.append(b)
    p3_lr.append(float(lr));p3_exceed.append(not okay or lr>=observed_lr-1e-6)
    save(f'p3_{b:03d}',{'b':b,'fit':fit,'LR':lr,'exceed':p3_exceed[-1]})
p3={**mc_result(p3_exceed),'computed':True,'observed_LR':observed_lr,'draw_LR':p3_lr,'failed':p3_failed}
save('p3_result',p3)
save('result',{'fixture':str(FIXTURE),'spec':'v1.7','B':B,'P1':p1,'P2':p2,'P3':p3,
               'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
               'elapsed_seconds':time.monotonic()-start,'note':'Synthetic agreement fixture, not calibration.'})
