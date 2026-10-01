"""Independent C8 from approved development-only z moments and public errors.

Primary uncertainty holds observed adjusted R² fixed (v1.3). It therefore
does not include uncertainty in that quantity. Headcount weights are planning
weights; fixed Gate A information weights are still needed for P2.
"""
import csv
import gzip
import json
import os
import platform
from pathlib import Path
import numpy as np
import scipy
from lambda_v1 import marginal_reliability,target_distribution


ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
SPEC=os.environ.get('CODEX_C8_SPEC','v1.3')
assert SPEC in ('v1.3','v1.4')
DESIGN=ROOT/'Analysis/MASLD_Model_Benchmark/executions'/('model-a-gateA-design-20260930T121812Z'
    if SPEC=='v1.4' else 'model-a-gateA-design-20260930T115256Z')
MOMENTS=ROOT/'Analysis/MASLD_Model_Benchmark/executions'/('model-a-gateA-r2-20260930T121832Z'
    if SPEC=='v1.4' else 'model-a-gateA-r2-20260930T120017Z')
PRIOR=ROOT/'docs/technical/agent_exchange/codex_lambda_20260929'
B=300


def read_table(path):
    with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


moments=json.loads((MOMENTS/'r2_moments_observed.json').read_text())
people=[r for r in read_table(DESIGN/'design_individuals.tsv') if r['set']=='development']
assert len(people)==558 and len({r['individual_id'] for r in people})==558
assert len({r['cluster_id'] for r in people})==558
by_id={r['individual_id']:r for r in people}
cross=read_table(ROOT/'Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/seal-20260929T141204Z/frozen_crosswalk.tsv')
run_by_id={r['individual_id']:r['run'] for r in cross if r['run_role']=='development' and r['unit_library']=='True' and r['S']}
with (ROOT/'RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv').open() as handle:
    meta={r['sample_id']:r for r in csv.DictReader(handle)}
labels={i:int(r['S']) if r['cohort']=='GSE213621' else int(float(meta[run_by_id[i]]['fibrosis_stage']))
        for i,r in by_id.items()}
cohorts=sorted(moments)
assert set(cohorts)=={r['cohort'] for r in people}
audit={};r2={};raw_r2={};weights={}
for cohort in cohorts:
    m=moments[cohort];ppl=sorted((r for r in people if r['cohort']==cohort),key=lambda r:r['cluster_id'])
    n=len(ppl);weights[cohort]=n;assert n==m['n']
    y=np.asarray([int(r['S']) for r in ppl],float);y-=y.mean()
    np.testing.assert_allclose(y@y,m['StS'],atol=1e-10)
    non_axes=[j for j,name in enumerate(m['columns']) if not name.startswith('axis_')]
    columns=[];constant_non_axes=[]
    for j in non_axes:
        name=m['columns'][j]
        base=name.removesuffix('_missing')
        values=np.asarray([float(r[base]) if r[base] not in ('NA','') else np.nan for r in ppl])
        missing=~np.isfinite(values)
        values=missing.astype(float) if name.endswith('_missing') else np.where(missing,np.nanmean(values),values)
        if np.ptp(values)==0:constant_non_axes.append(j)
        columns.append(values-values.mean())
    x=np.column_stack(columns)
    gram=np.asarray(m['XtX']);cross_y=np.asarray(m['XtS'])
    np.testing.assert_allclose(x.T@x,gram[np.ix_(non_axes,non_axes)],atol=1e-9,rtol=1e-9)
    np.testing.assert_allclose(x.T@y,cross_y[non_axes],atol=1e-9,rtol=1e-9)
    # Scale the Gram system so the small FFPE column does not lose rank.
    active=[j for j in range(len(gram)) if j not in constant_non_axes and gram[j,j]>0]
    reduced=gram[np.ix_(active,active)]
    scale=np.sqrt(np.diag(reduced));assert np.all(scale>0)
    correlation=reduced/np.outer(scale,scale);score=cross_y[active]/scale
    rank=int(np.linalg.matrix_rank(correlation))
    assert rank==m['rank_p'],(cohort,rank,m['rank_p'])
    raw=float(score@np.linalg.solve(correlation,score)/m['StS'])
    adjusted=max(0.,1-(1-raw)*(n-1)/(n-rank-1))
    np.testing.assert_allclose([raw,adjusted],[m['r2_raw'],m['r2_adjusted']],atol=1e-10,rtol=1e-9)
    raw_r2[cohort]=raw;r2[cohort]=adjusted
    audit[cohort]={'n':n,'rank':rank,'raw_r_squared':raw,'adjusted_r_squared':adjusted,
                   'non_ancestry_crossproducts':'match independently centred individual export'}

previous=json.loads((PRIOR/'lambda_558.json').read_text())
sources=('GSE193066','PXD051911_A1')
points={}
for source in sources:
    old=previous['sources'][source];spec=old['conservative_spec_by_headcount']
    marginal=old['specifications'][spec]['per_cohort']
    per={c:(marginal[c]-r2[c])/(1-r2[c]) for c in cohorts}
    raw_per={c:(marginal[c]-raw_r2[c])/(1-raw_r2[c]) for c in cohorts}
    points[source]={'carried_spec':spec,'per_cohort':per,
                    'lambda_headcount':sum(weights[c]*per[c] for c in cohorts)/558,
                    'raw_r_squared_sensitivity':sum(weights[c]*raw_per[c] for c in cohorts)/558}

draw_indices={b:{c:[] for c in cohorts} for b in range(B)}
with gzip.open(MOMENTS/'p2_resample_development.tsv.gz','rt') as handle:
    for row in csv.DictReader(handle,delimiter='\t'):
        b=int(row['draw'])
        if b>=B:break
        draw_indices[b][row['stratum']].append(row['cluster_id'])
resampled_r2={b:{} for b in range(B)}
with gzip.open(MOMENTS/'r2_by_resample.tsv.gz','rt') as handle:
    for row in csv.DictReader(handle,delimiter='\t'):
        b=int(row['draw'])
        if b>=B:break
        raw=float(row['r2_raw']);n=int(row['n']);p=int(row['rank_p'])
        adjusted=max(0.,1-(1-raw)*(n-1)/(n-p-1))
        np.testing.assert_allclose(adjusted,float(row['r2_adjusted']),atol=1e-10)
        resampled_r2[b][row['cohort']]=adjusted
error_draws={s:json.loads((PRIOR/f'lambda_558_{s}_draws.json').read_text()) for s in sources}
cluster_people={r['cluster_id']:r['individual_id'] for r in people}
draws={s:[] for s in sources}
for b in range(B):
    rng=np.random.default_rng(np.random.SeedSequence([20260923,2,b]));params={}
    for cohort in cohorts:
        original=sorted(r['cluster_id'] for r in people if r['cohort']==cohort)
        expected=[original[j] for j in rng.integers(0,len(original),len(original))]
        assert draw_indices[b][cohort]==expected,(b,cohort)
        selected=[labels[cluster_people[k]] for k in expected]
        params[cohort]=target_distribution(selected,collapsed=cohort=='GSE213621')
    for source in sources:
        entry=error_draws[source][b];assert entry['b']==b
        marginal={c:marginal_reliability(*params[c],entry['var_e']) for c in cohorts}
        per={c:(marginal[c]-r2[c])/(1-r2[c]) for c in cohorts}
        sensitive={c:(marginal[c]-resampled_r2[b][c])/(1-resampled_r2[b][c]) for c in cohorts}
        draws[source].append({'b':b,'per_cohort':per,
                              'lambda_headcount':sum(weights[c]*per[c] for c in cohorts)/558,
                              'recomputed_r_squared_sensitivity':sum(weights[c]*sensitive[c] for c in cohorts)/558})
for source in sources:
    values=np.asarray([r['lambda_headcount'] for r in draws[source]])
    assert np.isfinite(values).all()
    points[source]['ci95_fixed_r_squared']=np.quantile(values,[.025,.975]).tolist()
    points[source]['ci95_recomputed_r_squared_sensitivity']=np.quantile(
        [r['recomputed_r_squared_sensitivity'] for r in draws[source]],[.025,.975]).tolist()
    points[source]['invalid_lambda_draws']=int(np.sum(values<=.05))
out={'spec':SPEC,'design':str(DESIGN),'moments':str(MOMENTS),'n':558,'bootstrap':B,'r_squared_audit':audit,'sources':points,'draws':draws,
     'target_stream':'SeedSequence([20260923,2,b]); matches approved fixed resample file',
     'public_pair_error_seed':20260929,'python':platform.python_version(),
     'numpy':np.__version__,'scipy':scipy.__version__,
     'limits':'Headcount weights. Primary holds observed adjusted R² fixed; no R² estimation uncertainty. Not a P2 result.'}
output=Path(os.environ['CODEX_REC_OUTPUT'])/('lambda_558_conditional_'+SPEC.replace('.','')+'.json')
output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'r_squared':audit,'sources':points}),flush=True)
