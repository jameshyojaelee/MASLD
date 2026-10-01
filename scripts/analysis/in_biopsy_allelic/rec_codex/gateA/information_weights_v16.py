"""Expected information for the selected RNA beta-binomial mixture, outcome-free.

The row weights are sufficient for the owning builder to assemble the full
expected Fisher matrix with its ancestry covariates. No realized allele counts
or per-person ancestry values are read here.
"""
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import sys
from collections import Counter
from functools import lru_cache
from pathlib import Path
import numpy as np
import scipy
from scipy.special import expit,logsumexp
from scipy.stats import betabinom

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from likelihood_v1 import EPS,bb_terms,fixed_bb

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
BASE=ROOT/'Analysis/MASLD_Model_Benchmark/executions'
DESIGN=BASE/'model-a-gateA-design-20260930T121812Z'
RULES=BASE/'in-biopsy-allelic-t2b-20260930T123844Z'
QC=BASE/'in-biopsy-allelic-obs-20260926T155651Z/qc/r1c-20260930T123747Z'
GTEX=ROOT/'data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz'


def table(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt') as fh:return list(csv.DictReader(fh,delimiter='\t'))


def r_coordinate(rho):
    p=(rho-EPS)/(1-2*EPS)
    return math.log(p/(1-p))


@lru_cache(maxsize=65536)
def information(n,k,F,eta,r,ph,plow,phigh,e,phi):
    lo=max(k,(F*n+99)//100)
    if lo>n-lo:return (0.,0.,0.,0.)
    counts=np.arange(lo,n-lo+1)
    hmass,hscore,_=bb_terms(counts,n,eta,r)
    priors=np.asarray([ph,plow,phigh],float)
    with np.errstate(divide='ignore'):
        parts=np.stack((hmass,fixed_bb(counts,n,e,phi),fixed_bb(counts,n,1-e,phi)))+np.log(priors)[:,None]
    mass=logsumexp(parts,axis=0);den=logsumexp(mass)
    probability=np.exp(mass-den);included=float(np.exp(den))
    score=np.exp(parts[0]-mass)[:,None]*hscore
    score-=probability@score
    fisher=included*np.einsum('n,ni,nj->ij',probability,score,score)
    assert np.isfinite(fisher).all() and 0<=included<=1+1e-9
    np.testing.assert_allclose(probability@score,0.,atol=1e-10)
    assert np.linalg.eigvalsh(fisher)[0]>=-1e-10
    return (included,float(fisher[0,0]),float(fisher[0,1]),float(fisher[1,1]))


def check_numeric_hessian():
    """Independent scipy probability calculation and finite-difference log Hessian."""
    n=35;k=2;F=10;eta=.7;r=r_coordinate(.02);e=.025;phi=.04
    prior=np.asarray([.4,.35,.25]);counts=np.arange(4,32)
    def logq(x):
        mean=expit(x[0]);rho=EPS+(1-2*EPS)*expit(x[1]);t=(1-rho)/rho
        raw=np.stack((betabinom.logpmf(counts,n,mean*t,(1-mean)*t),
            betabinom.logpmf(counts,n,e*(1-phi)/phi,(1-e)*(1-phi)/phi),
            betabinom.logpmf(counts,n,(1-e)*(1-phi)/phi,e*(1-phi)/phi)))+np.log(prior)[:,None]
        mass=logsumexp(raw,axis=0);den=logsumexp(mass)
        return mass-den,float(np.exp(den))
    point=np.array([eta,r]);lp,included=logq(point);probability=np.exp(lp);step=1e-4
    numeric=np.zeros((2,2))
    for j in range(2):
        a=np.eye(2)[j]*step
        numeric[j,j]=-included*float(probability@(logq(point+a)[0]-2*lp+logq(point-a)[0]))/step**2
    a=np.array([step,0]);b=np.array([0,step])
    numeric[0,1]=numeric[1,0]=-included*float(probability@(logq(point+a+b)[0]-logq(point+a-b)[0]
                                    -logq(point-a+b)[0]+logq(point-a-b)[0]))/(4*step**2)
    value=information(n,k,F,eta,r,*prior,e,phi)
    analytic=np.array([[value[1],value[2]],[value[2],value[3]]])
    np.testing.assert_allclose(value[0],included,atol=1e-12)
    np.testing.assert_allclose(analytic,numeric,atol=2e-5,rtol=2e-5)
    # Common prior scaling changes expected row count, but not conditional information.
    scaled=information(n,k,F,eta,r,*(prior*.8),e,phi)
    np.testing.assert_allclose(scaled,np.asarray(value)*.8,atol=1e-12)
    return {'analytic':analytic.tolist(),'numeric_hessian':numeric.tolist(),'inclusion':included}


check=check_numeric_hessian()
people={r['individual_id']:r for r in table(DESIGN/'design_individuals.tsv') if r['set']=='development'}
assert len(people)==558
cohorts=sorted({r['cohort'] for r in people.values()})
rules={r['cohort']:tuple(int(r[k]) for k in ('min_dp','min_minor_reads','min_minor_pct'))
       for r in table(RULES/'cohort_rules.tsv') if r['set']=='development' and r['cohort'] in cohorts}
inputs={r['cohort']:r for r in table(RULES/'gateA_cohort_inputs.tsv') if r['set']=='development' and r['cohort'] in cohorts}
cal={r['cohort']:r for r in table(QC/'cohort_model_inputs.tsv') if r['rule_set']=='development' and r['cohort'] in cohorts}
leads={r['gene_id']:r['lead_variant_id'] for r in table(RULES/'tag_table/tag_table.tsv') if r['kept_primary']=='True'}
afc={r['gene_id']:r for r in table(GTEX) if r['gene_id'] in leads}
assert set(afc)==set(leads) and all(afc[g]['variant_id']==leads[g] for g in leads)
alpha={g:float(np.clip(float(afc[g]['log2_aFC'])*math.log(2),-3,3)) for g in leads}
rho_specs={'rho_0.02':.02,'geuvadis_rho':float(inputs[cohorts[0]]['rho_het_geuvadis'])}
out=Path(os.environ['CODEX_REC_OUTPUT'])/('information_weights_v16_'+os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)
totals=Counter();written=0;seen=set()
with gzip.open(out/'weights_for_full_z.tsv.gz','wt') as fh,gzip.open(DESIGN/'design_rows.tsv.gz','rt') as source:
    writer=csv.writer(fh,delimiter='\t');writer.writerow(['rho_spec','cohort','individual_id','gene_id','alpha_plant',
        's_t','S','n','inclusion_probability','W_eta_eta','W_eta_r','W_r_r'])
    for row in csv.DictReader(source,delimiter='\t'):
        if row['set']!='development':continue
        c=row['cohort'];g=row['gene_id'];person=people[row['individual_id']];n=int(row['n'])
        assert person['cohort']==c
        key=(row['individual_id'],g);assert key not in seen;seen.add(key)
        d,k,F=rules[c]
        if n<d:continue
        prior=np.array([float(row[v]) for v in ('P_H','P_R','P_A')])
        assert np.isfinite(prior).all() and (prior>=0).all() and prior.sum()<=1+1e-6
        if prior.sum()>1:prior/=prior.sum()
        s=int(row['s_t']);assert s in (-1,1)
        if s==1:prior=prior[[0,2,1]]
        e=max(math.exp(float(person['z_b_log_e']))*float(cal[c]['e_multiplier_applied']),1e-4)
        eta=alpha[g]+s*float(inputs[c]['omega_c']);phi=float(inputs[c]['phi'])
        for name,rho in rho_specs.items():
            value=information(n,k,F,eta,r_coordinate(rho),*prior,e,phi)
            writer.writerow([name,c,row['individual_id'],g,alpha[g],s,person['S'],n,*value])
            totals[name]+=value[0];written+=1
expected_path=BASE/'codex-rec-20260929T142434Z/expected_counts_v16_21986201/summary.json'
previous=json.loads(expected_path.read_text())
for name,v in totals.items():np.testing.assert_allclose(v,previous['summaries'][name]['total_expected_included'],rtol=1e-10)
paths=[DESIGN/'design_individuals.tsv',DESIGN/'design_rows.tsv.gz',RULES/'cohort_rules.tsv',
       RULES/'gateA_cohort_inputs.tsv',RULES/'tag_table/tag_table.tsv',QC/'cohort_model_inputs.tsv',GTEX,expected_path]
summary={'spec':'v1.6','plant':'kappa=delta=omega_S=dispersion_slopes=tau=0; alpha from capped signed GTEx',
    'biological_n':558,'written_rows':written,'total_inclusion':dict(totals),'rho_specs':rho_specs,
    'kernel_check':check,'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
    'inputs_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    'weight_definition':'W_ab = P(included) * E[score_a score_b | included], a,b = eta,r',
    'builder_jacobian':'eta derivatives: alpha_g:1; kappa:alpha_g*S_centered; delta:alpha_g*z; '
        'omega_S:s_t*S_centered. r derivatives: r_g:1; rho_S:S_centered; rho_b:b_centered; rho_d:d_centered.',
    'limits':'Approved cohort-mixture priors and assumed GTEx transport. No outcome counts or ancestry values. '
        'Full adjusted information and power require the full covariate Jacobian and nuisance projection.'}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({'output':str(out),'written_rows':written,'total_inclusion':dict(totals),'kernel_check':check}),flush=True)
