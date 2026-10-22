#!/usr/bin/env python3
"""Synthetic generative-law checks; no B3 allocation, calibration or fits."""
import argparse
import ast
from dataclasses import replace
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
SOURCES=('likelihood_v1.py','repair_p3_bb_gamma_difference_v3.py','calibration_generate_base_v1.py')


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True);ap.add_argument('--inspect',action='store_true');a=ap.parse_args()
    require(sha(Path(__file__))==a.expected_runner_sha and sha(a.pins)==a.expected_pins_sha,'Runner/pins differ')
    pins=json.loads(a.pins.read_text())
    require(set(pins)=={str(HERE/n) for n in SOURCES},'Source family differs')
    for p,h in pins.items():
        require(sha(Path(p))==h,'Source differs');ast.parse(Path(p).read_text())
    if a.inspect:print(json.dumps(dict(source_checks=True,synthetic_only=True,fit_count=0)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
        'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
        'Threads not fixed')
    out=REC/('calibration_generate_base_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        q=archive/p.name;q.write_bytes(p.read_bytes());require(sha(q)==(a.expected_runner_sha if p==Path(__file__) else pins[str(p)]),
            'Archive differs')
    import numpy as np
    import scipy
    from scipy.integrate import quad
    from scipy.special import expit,logsumexp
    from scipy.stats import norm
    original=load('base_generate_original',archive/'likelihood_v1.py')
    gamma=load('base_generate_gamma',archive/'repair_p3_bb_gamma_difference_v3.py')
    base=load('base_generate',archive/'calibration_generate_base_v1.py');start=time.perf_counter()
    rows=[original.Row(0,'synthetic'+str(i),'fixture',20,10,.6,.3,.1,.01,.02,.05,
        1,float(i-1),(0.,),0.,0.,10,2,10) for i in range(3)]
    rows.append(replace(rows[2],gene=1,individual='synthetic3',orientation=-1,background=.2,duplication=-.1))
    model=original.Model(rows,gene_count=5);p=model.initial();p[:2]=[-1.2,.7];p[10]=.08;p[11]=.04
    p[12]=.15;p[13:]=[.05,.03,.04]
    def save(name,value):
        with (out/(name+'.json')).open('x') as f:json.dump(value,f,allow_nan=False,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    save('protocol',dict(source_sha256=pins,seed_zero_tau=2026100701,seed_positive_tau=2026100702,
        zero_tau_replicates=20000,positive_tau_replicates=6000,positive_tau=.2,
        model_G=5,model_Z=1,synthetic_rows=4,independent_unit='simulated dataset; shared stage slope per gene',
        RNG_identity='SeedSequence(seed,spawn_key=(replicate_index,)).spawn(3); independent of job order',
        environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
        Monte_Carlo_check_rule='absolute discrepancy <=6 estimated Monte Carlo SE +1e-8',
        no_fits=True,no_B3_scenarios=True,no_calibration=True,no_human_data=True))
    # Replay checks include latent effects and all pre-selection counts.
    first=base.BaseGenerator(model,p,.2,2026100700,original,gamma)
    second=base.BaseGenerator(model,p,.2,2026100700,original,gamma)
    for _ in range(10):
        arows,akeep,atruth=first.simulate();brows,bkeep,btruth=second.simulate()
        require(arows==brows and akeep==bkeep and np.array_equal(atruth['gene_effects'],btruth['gene_effects'])
            and atruth['genotype_classes']==btruth['genotype_classes'],'Seed replay differs')
    left={i:first.simulate(i) for i in (37,12)};right={i:second.simulate(i) for i in (12,37)}
    for i in (12,37):require(left[i][0]==right[i][0]
        and np.array_equal(left[i][2]['gene_effects'],right[i][2]['gene_effects']),
        'Replicate identity depends on execution order')
    zero=base.BaseGenerator(model,p,0.,2026100701,original,gamma);hist=np.zeros(21,dtype=int);classes=np.zeros(3,dtype=int)
    count=20000;retained=0;class_count_sums=np.zeros(3);negative_classes=np.zeros(3,dtype=int)
    for _ in range(count):
        drawn,kept,truth=zero.simulate();require(np.all(truth['gene_effects']==0),'Zero tau produced effects')
        row=drawn[1];genotype=truth['genotype_classes'][1];classes[genotype]+=1
        class_count_sums[genotype]+=row.a;negative_classes[truth['genotype_classes'][3]]+=1
        keep=truth['included'][1]
        require(keep==(row.a in row.allowed) and tuple(r for r,k in zip(drawn,truth['included']) if k)==kept,
            'Integer inclusion/retained rows differ')
        if keep:hist[row.a]+=1;retained+=1
    row=rows[1];eta,r,_,_=model.predictors(row,p,0.,False);support=row.allowed
    lp=original.mixture_terms(row,np.arange(21),eta,r,False)[0]-math.log(float(row.prior.sum()))
    probabilities=np.exp(lp);require(abs(probabilities.sum()-1)<=1e-12,'Independent unconditional law invalid')
    selection=float(probabilities[support].sum());selection_SE=math.sqrt(selection*(1-selection)/count)
    require(abs(retained/count-selection)<=6*selection_SE+1e-8,'Empirical inclusion disagrees with unconditional law')
    conditioned=probabilities[support]/selection
    literal=np.array([math.exp(original.conditional_terms(replace(row,a=int(k)),eta,r,False)[0]) for k in support])
    require(np.max(np.abs(conditioned-literal))<=1e-12,'Independent conditional normalizer differs')
    frequency=hist[support]/retained;SE=np.sqrt(conditioned*(1-conditioned)/retained)
    require(np.all(np.abs(frequency-conditioned)<=6*SE+1e-8),'Retained-count distribution disagrees')
    weights=row.prior/row.prior.sum();class_SE=np.sqrt(weights*(1-weights)/count)
    require(np.all(np.abs(classes/count-weights)<=6*class_SE+1e-8),'Modeled genotype probabilities disagree')
    class_means=class_count_sums/classes
    expected_means=np.array([20*expit(eta-row.orientation*row.omega),20*row.e,20*(1-row.e)])
    rho=1e-6+(1-2e-6)*expit(r)
    class_variance=20*(expected_means/20)*(1-expected_means/20)*(1+19*np.array([rho,row.phi,row.phi]))
    class_mean_SE=np.sqrt(class_variance/classes)
    require(np.all(np.abs(class_means-expected_means)<=6*class_mean_SE+1e-8),'Genotype-specific count means disagree')
    negative_weights=rows[3].prior/rows[3].prior.sum()
    require(np.all(np.abs(negative_classes/count-negative_weights)<=6*np.sqrt(negative_weights*(1-negative_weights)/count)+1e-8),
        'Negative-orientation genotype probabilities disagree')
    save('zero_tau',dict(replicates=count,retained=retained,retained_fraction=retained/count,
        theoretical_selection=selection,selection_MC_SE=selection_SE,support=support.tolist(),
        observed_frequency=frequency.tolist(),theoretical_conditional=conditioned.tolist(),
        conditional_MC_SE=SE.tolist(),observed_classes=classes.tolist(),theoretical_class_weights=weights.tolist(),
        observed_class_count_means=class_means.tolist(),expected_class_count_means=expected_means.tolist(),
        class_mean_MC_SE=class_mean_SE.tolist(),negative_orientation_classes=negative_classes.tolist(),
        negative_orientation_weights=negative_weights.tolist(),
        independent_conditional_max_difference=float(np.max(np.abs(conditioned-literal)))))
    positive=base.BaseGenerator(model,p,.2,2026100702,original,gamma);N=6000;indicators=[];effect_times_inclusion=[];effects=[]
    covariance_rows=(0,2,3);covariances={j:[] for j in covariance_rows}
    base_means={j:rows[j].n*expit(model.predictors(rows[j],p,0.,False)[0]-rows[j].orientation*rows[j].omega)
        for j in covariance_rows}
    second_effects=[]
    for _ in range(N):
        drawn,kept,truth=positive.simulate();u=float(truth['gene_effects'][0]);keep=truth['included'][2]
        effects.append(u);indicators.append(float(keep));effect_times_inclusion.append(u*keep)
        second_effects.append(float(truth['gene_effects'][1]))
        for j in covariance_rows:
            value=float(truth['gene_effects'][rows[j].gene])*(drawn[j].a-base_means[j])
            covariances[j].append(value if truth['genotype_classes'][j]==0 else 0.)
    row=rows[2];support=row.allowed
    def inclusion(z):
        eta,r,_,_=model.predictors(row,p,.2*z,False)
        return math.exp(float(logsumexp(original.mixture_terms(row,support,eta,r,False)[0]))-math.log(float(row.prior.sum())))
    C,Cerr=quad(lambda z:norm.pdf(z)*inclusion(z),-10.,10.,epsabs=1e-11,epsrel=1e-11)
    uC,uCerr=quad(lambda z:.2*z*norm.pdf(z)*inclusion(z),-10.,10.,epsabs=1e-11,epsrel=1e-11)
    values=(np.asarray(indicators),np.asarray(effect_times_inclusion));expectations=(C,uC);mc=[]
    for values_i,expected,error in zip(values,expectations,(Cerr,uCerr)):
        mean=float(values_i.mean());se=float(values_i.std(ddof=1)/np.sqrt(N))
        require(abs(mean-expected)<=6*se+1e-8+error,'Random-effect-before-selection moment disagrees')
        mc.append(dict(observed=mean,expected=expected,MC_SE=se,quad_error_estimate=error))
    effects=np.asarray(effects);require(abs(float(effects.mean()))<=6*.2/np.sqrt(N)
        and abs(float(effects.var(ddof=1))-.04)<=6*.04*np.sqrt(2/(N-1)),'Pre-selection normal effect moments disagree')
    covariance_checks=[]
    for j in covariance_rows:
        row=rows[j];eta0=model.predictors(row,p,0.,False)[0]-row.orientation*row.omega
        weight=float(row.prior[0]/row.prior.sum())
        expected,error=quad(lambda z:weight*norm.pdf(z)*(.2*z)*
            (row.n*expit(eta0+.2*z*row.stage)-base_means[j]),-10.,10.,epsabs=1e-11,epsrel=1e-11)
        values=np.asarray(covariances[j]);mean=float(values.mean());se=float(values.std(ddof=1)/np.sqrt(N))
        require(abs(expected)>6*se+1e-8,'Covariance fixture cannot distinguish absent random slope')
        require(abs(mean-expected)<=6*se+1e-8+error,'Gene-indexed signed effect/count covariance disagrees')
        covariance_checks.append(dict(row=j,gene=row.gene,stage=row.stage,orientation=row.orientation,
            observed=mean,expected=expected,MC_SE=se,quad_error_estimate=error))
    cross=effects*np.asarray(second_effects);cross_SE=float(cross.std(ddof=1)/np.sqrt(N))
    require(abs(float(cross.mean()))<=6*cross_SE+1e-8,'Independent gene effects disagree')
    # An impossible read-depth predicate yields empty data, never replacement draws.
    empty_model=original.Model([replace(r,depth=21) for r in rows],gene_count=5)
    empty=base.BaseGenerator(empty_model,p,.2,2026100703,original,gamma)
    for _ in range(10):
        drawn,kept,truth=empty.simulate();require(len(drawn)==4 and not kept and not any(truth['included']),
            'Ineligible rows were replaced or retained')
    save('positive_tau',dict(replicates=N,tau=.2,moments=mc,unconditional_effect_mean=float(effects.mean()),
        signed_preselection_het_count_covariances=covariance_checks,
        cross_gene_effect_product_mean=float(cross.mean()),cross_gene_effect_product_MC_SE=cross_SE,
        unconditional_effect_variance=float(effects.var(ddof=1)),expected_effect_variance=.04,
        only_simulation_evaluator_reads_latent_effects=True,empty_datasets_retained=10))
    for path,h in pins.items():require(sha(Path(path))==h,'Input changed during check')
    save('summary',dict(completed=True,seconds=time.perf_counter()-start,replay_datasets=10,
        zero_tau_replicates=count,positive_tau_replicates=N,empty_datasets_retained=10,
        original_count_law_and_selection_checks_passed=True,no_fits=True,no_calibration=True,no_human_data=True))
    print(json.dumps(dict(completed=True,generative_law_checks_passed=True,no_calibration=True)))


if __name__=='__main__':main()
