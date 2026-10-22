"""Synthetic base count generation before inclusion; not a B3 scenario recipe.

The caller must provide a fixed synthetic design, parameters, tau and seed.
Gene-specific stage slopes precede genotype classes and unconditional RNA counts.
Read inclusion then removes rows. No redraw forces an eligible observation.
This module supplies no label-error model, technical plant, allocation or fit.
"""
from dataclasses import replace
import hashlib
import os
from pathlib import Path

ORIGINAL_SHA='41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6'
GAMMA_SHA='546e7ea2f1db21e6fcdd6aed72204d3bbe48c2d818996ab21de5f5ab5215bcf3'


def require(ok,message):
    if not ok:raise ValueError(message)


class BaseGenerator:
    def __init__(self,model,parameters,tau,seed,original,gamma):
        require(os.environ.get('SLURM_JOB_ID','').isdigit()
            and not os.environ.get('SLURM_ARRAY_TASK_ID'),'Individual compute allocation required')
        import numpy as np
        from scipy.special import logsumexp
        for module,digest in ((original,ORIGINAL_SHA),(gamma,GAMMA_SHA)):
            require(hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()==digest,'Count-law source differs')
        require(type(model) is original.Model,'Unmodified original synthetic model required')
        require(model.G<=8 and model.Z<=7 and 0<len(model.rows)<=805,'Synthetic base-generator scope exceeded')
        require(type(seed) is int and seed>=0 and np.isfinite(tau) and tau>=0,'Finite tau and integer seed required')
        require(all(all(type(x) is int and x>=0 for x in (r.n,r.depth,r.minor,r.percent))
            and r.n<=5000 and r.percent<=50 and 0<r.e<1 and 0<=r.phi<1
            and r.orientation in (-1,1) and np.isfinite(r.stage) for r in model.rows),
            'Unsupported count/error/orientation/stage input')
        p=np.asarray(parameters,dtype=float).copy()
        require(p.shape==(model.size,) and np.isfinite(p).all(),'Finite fixed parameters required')
        p.setflags(write=False)
        self.np,self.logsumexp=np,logsumexp
        self.model,self.p,self.tau,self.original,self.gamma=model,p,float(tau),original,gamma
        self.parameter_bytes=p.tobytes();self.rows=tuple(model.rows)
        self.row_snapshots=tuple(tuple(vars(r).items()) for r in self.rows)
        self.seed=seed;self.replicates=0;self.zero_tau_cdfs={}
        self.class_weights=[]
        for row in self.rows:
            weights=np.asarray(row.prior,dtype=float)
            require(weights.shape==(3,) and np.isfinite(weights).all() and np.all(weights>=0)
                and weights.sum()>0,'Nonzero modeled genotype support required')
            # Row.prior omits the uncalled/uncalled state. Condition on modeled
            # support before read inclusion; its common scale cancels in L(a|K).
            weights=weights/weights.sum();weights.setflags(write=False);self.class_weights.append(weights)

    def cdf(self,row,eta,dispersion,genotype):
        np=self.np;k=np.arange(row.n+1)
        if genotype==0:
            lp=self.gamma.stable_bb_terms(k,row.n,eta-row.orientation*row.omega,dispersion,False)[0]
        else:
            mean=row.e if genotype==1 else 1-row.e
            lp=self.original.fixed_bb(k,row.n,mean,row.phi)
        require(lp.shape==(row.n+1,) and np.isfinite(lp).all(),'Unconditional mass invalid')
        logsum=self.logsumexp(lp)
        require(abs(logsum)<=1e-8,'Unconditional mass normalization disagrees')
        cumulative=np.cumsum(np.exp(lp-logsum));cumulative[-1]=1.
        require(np.isfinite(cumulative).all() and np.all(np.diff(cumulative)>=0),'Count CDF invalid')
        cumulative.setflags(write=False);return cumulative

    def simulate(self,replicate_index=None):
        np=self.np
        require(self.p.tobytes()==self.parameter_bytes
            and tuple(tuple(vars(r).items()) for r in self.rows)==self.row_snapshots,'Fixed synthetic design changed')
        index=self.replicates if replicate_index is None else replicate_index
        require(type(index) is int and index>=0,'Nonnegative fixed replicate index required')
        effect_stream,class_stream,count_stream=np.random.SeedSequence(self.seed,spawn_key=(index,)).spawn(3)
        # One shared effect per gene, drawn before either genotype or count.
        effects=np.random.default_rng(effect_stream).normal(0.,self.tau,self.model.G)
        classes_uniform=np.random.default_rng(class_stream).random(len(self.rows))
        count_uniform=np.random.default_rng(count_stream).random(len(self.rows))
        drawn=[];retained=[];classes=[];included=[]
        for ordinal,row in enumerate(self.rows):
            weights=self.class_weights[ordinal];class_cdf=np.cumsum(weights);class_cdf[-1]=1.
            genotype=int(np.searchsorted(class_cdf,classes_uniform[ordinal],side='left'))
            eta,r,_,_=self.model.predictors(row,self.p,float(effects[row.gene]),False)
            key=(ordinal,genotype)
            if self.tau==0 and key in self.zero_tau_cdfs:cumulative=self.zero_tau_cdfs[key]
            else:
                cumulative=self.cdf(row,eta,r,genotype)
                if self.tau==0:self.zero_tau_cdfs[key]=cumulative
            count=int(np.searchsorted(cumulative,count_uniform[ordinal],side='left'))
            require(0<=count<=row.n,'Count escaped integer support')
            generated=replace(row,a=count);lo=max(row.minor,(row.percent*row.n+99)//100)
            keep=row.n>=row.depth and lo<=count<=row.n-lo
            drawn.append(generated);classes.append(genotype);included.append(bool(keep))
            if keep:retained.append(generated)
        self.replicates+=1
        # Fit code receives retained_rows only. Truth remains a separate return
        # value for the synthetic evaluator; do not use it as fit starts/labels.
        truth=dict(gene_effects=effects,genotype_classes=classes,included=included,
            replicate_index=index,simulation_calls=self.replicates,all_template_rows=len(drawn),retained_rows=len(retained))
        return tuple(drawn),tuple(retained),truth
