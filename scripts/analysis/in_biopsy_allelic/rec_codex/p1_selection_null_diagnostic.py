"""Exact synthetic RNA-mixture diagnostic; not P1 size or power calibration."""
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import scipy
from scipy.special import expit, logit
from scipy.stats import betabinom

from likelihood_v1 import EPS, Model, Row, mixture_terms


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Requires a compute allocation.')
    out=Path(os.environ['CODEX_REC_OUTPUT'])/('p1_selection_null_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    alpha=np.array([.7,-.6,1.1,-.3,.5,-.9,.8,-.4])
    orientations=np.array([1,-1,1,-1,-1,1,-1,1])
    z=.8
    cumulative=expit(np.array([-.7,.7])-z)
    q_pre=np.r_[cumulative[0],np.diff(cumulative),1-cumulative[-1]]
    rows=[Row(g,'synthetic','synthetic',80,40,.7,.3,0.,.002,.02,.01,
              int(orientations[g]),0.,(z,),0.,0.,20,3,15) for g in range(8)]
    model=Model(rows)
    base=np.zeros(model.size)
    base[:8]=alpha
    base[8:16]=logit((.02-EPS)/(1-2*EPS))
    patterns={'all_selected':np.ones(8,dtype=bool),
              'alternating_selected':np.array([1,0,1,0,1,0,1,0],dtype=bool),
              'none_selected':np.zeros(8,dtype=bool)}
    result=[]
    for name,omega_s,rho_s in [('nuisance_off',0.,0.),('nuisance_on',.02,.10)]:
        p=base.copy()
        p[2*model.G+1+model.Z]=omega_s
        p[-3]=rho_s
        include=np.zeros((8,3))
        conditioned=np.zeros((8,3,57))
        means=np.zeros((8,3))
        for g,row in enumerate(rows):
            for s in range(3):
                current=replace(row,stage=float(s-1))
                eta,r,_,_=model.predictors(current,p)
                counts=np.arange(current.n+1)
                logmass,_,_=mixture_terms(current,counts,eta,r,need_hessian=False)
                mass=np.exp(logmass)
                # Independent, assay-native beta-binomial formula for all three classes.
                rho=EPS+(1-2*EPS)*expit(r)
                t=(1-rho)/rho
                m=expit(eta-current.orientation*current.omega)
                hom_t=(1-current.phi)/current.phi
                independently=(current.prior[0]*betabinom.pmf(counts,current.n,m*t,(1-m)*t)
                    +current.prior[1]*betabinom.pmf(counts,current.n,current.e*hom_t,(1-current.e)*hom_t)
                    +current.prior[2]*betabinom.pmf(counts,current.n,(1-current.e)*hom_t,current.e*hom_t))
                np.testing.assert_allclose(mass,independently,atol=2e-12,rtol=2e-10)
                np.testing.assert_allclose(mass.sum(),1.,atol=2e-12,rtol=0)
                selected=mass[current.allowed]
                include[g,s]=selected.sum()
                if not 0<include[g,s]<1:
                    raise ValueError('Inclusion probability must be strictly interior.')
                conditioned[g,s]=selected/include[g,s]
                means[g,s]=conditioned[g,s]@current.allowed/current.n
        posterior=[]
        for pattern_name,pattern in patterns.items():
            # Conditional independence of the eight generated candidate observations is
            # an explicit assumption of this synthetic diagnostic, not a MASLD claim.
            log_weights=np.log(q_pre)+np.where(pattern[:,None],np.log(include),np.log1p(-include)).sum(axis=0)
            weights=np.exp(log_weights-log_weights.max())
            q_post=weights/weights.sum()
            posterior.append({'pattern':pattern_name,'inclusion':pattern.astype(int).tolist(),
                              'q_post':q_post.tolist(),'max_absolute_stage_probability_shift':float(np.max(np.abs(q_post-q_pre)))})
        tv=np.abs(conditioned[:,0]-conditioned[:,2]).sum(axis=1)/2
        shift=max(x['max_absolute_stage_probability_shift'] for x in posterior)
        if name=='nuisance_off':
            np.testing.assert_allclose(tv,0.,atol=1e-14,rtol=0)
            if shift>=1e-14:raise ValueError('Conditional-independence control failed.')
        else:
            if not np.max(tv)>1e-5 or not shift>1e-8:
                raise ValueError('Planted nuisance counter-scenario is not distinguishable.')
        result.append({'regime':name,'omega_S':omega_s,'rho_S':rho_s,'inclusion_probabilities':include.tolist(),
                       'conditional_mean_effect_oriented_allele_fraction':means.tolist(),
                       'stage0_vs_stage2_total_variation_per_gene':tv.tolist(),'patterns':posterior})
    report={'status':'exact_synthetic_distribution_checks_passed','kappa':0.,'tau':0.,
            'alpha':alpha.tolist(),'tag_orientation':orientations.tolist(),'q_pre':q_pre.tolist(),
            'candidate_genes':8,'depth':80,'error_probability':.002,'homozygote_phi':.02,
            'genotype_frequencies':[.7,.3,0.],'het_rule':{'min_dp':20,'min_minor_reads':3,'min_minor_pct':15},
            'stage_centering':'fixed generating center1; no test statistic or fit',
            'regimes':result,'stochastic_operations':False,'python':sys.version,
            'numpy':np.__version__,'scipy':scipy.__version__,
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'likelihood_sha256':hashlib.sha256(Path('likelihood_v1.py').read_bytes()).hexdigest(),
            'real_allelic_outcomes_read':False,'fit_performed':False,'P1_size_or_power_evaluated':False,
            'limits':'Synthetic independent candidate genes and fixed known preselection stage law only. '
                     'Shows weak-null count dependence and selection-conditioned stage-law change; '
                     'does not prove P1 is miscalibrated or replace required repeated-null calibration.'}
    (out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':report['status'],
        'max_nuisance_on_total_variation':max(result[1]['stage0_vs_stage2_total_variation_per_gene']),
        'max_stage_probability_shift':max(x['max_absolute_stage_probability_shift'] for x in result[1]['patterns'])}),flush=True)


if __name__=='__main__':main()
