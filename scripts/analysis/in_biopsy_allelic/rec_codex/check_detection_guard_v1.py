"""Outcome-free feasibility of spec §4's n-star detection guard."""
import csv
import json
import os
from pathlib import Path

import numpy as np
from scipy.special import betainc
from scipy.stats import betabinom

from likelihood_v1 import allowed_counts


root=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
rule_path=root/'Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-t2b-20260929T220019Z/cohort_rules.tsv'
with rule_path.open() as handle:
    rules=list(csv.DictReader(handle,delimiter='\t'))
triples=sorted({(int(r['min_dp']),int(r['min_minor_reads']),int(r['min_minor_pct'])) for r in rules})
mean=.2;rho=.02;t=(1-rho)/rho;alpha=mean*t;beta=(1-mean)*t
result={str(rule):{'rule':rule,'n_star':None,'max_probability':0.,'depth_at_max':None,
                   'large_depth_limit':float(betainc(alpha,beta,1-rule[2]/100)-betainc(alpha,beta,rule[2]/100))}
        for rule in triples}
mass=np.array([1.])
cap=10000
for n in range(1,cap+1):
    k=np.arange(n)
    updated=np.zeros(n+1)
    updated[:-1]+=mass*(beta+n-1-k)/(t+n-1)
    updated[1:]+=mass*(alpha+k)/(t+n-1)
    mass=updated
    if n in (10,100,1000,6821,10000):
        assert abs(mass.sum()-1)<1e-10
        checkpoints=np.unique(np.r_[0,n//5,n//2,n-1,n])
        np.testing.assert_allclose(mass[checkpoints],betabinom.pmf(checkpoints,n,alpha,beta),atol=2e-12,rtol=1e-8)
    for rule in triples:
        if n<rule[0]:
            continue
        lo=max(rule[1],(rule[2]*n+99)//100)
        probability=float(mass[lo:n-lo+1].sum())
        record=result[str(rule)]
        if probability>record['max_probability']:
            record.update(max_probability=probability,depth_at_max=n)
        if probability>=.99 and record['n_star'] is None:
            record['n_star']=n
out={'mean':mean,'rho':rho,'all_integer_depths_checked':[1,cap],
     'observed_maximum_n_in_spec':6821,'rules':list(result.values()),
     'cohort_rules':[{k:r[k] for k in ('cohort','set','min_dp','min_minor_reads','min_minor_pct')} for r in rules]}
(Path(os.environ['CODEX_REC_OUTPUT'])/'detection_guard_v1.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
