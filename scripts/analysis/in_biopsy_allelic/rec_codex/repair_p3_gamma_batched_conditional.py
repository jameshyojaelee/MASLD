"""Prepared synthetic conditional mixtures with the checked batched Gamma BB.

The existing conditional assembly and node splitting are inherited. Only its
BB arithmetic is replaced; the older count-prefix tables are not evaluated at
nodes. No model method, optimizer, integration rule or statistical test changes.
This source alone does not admit a conditional, integral, fit or calibration.
"""
import hashlib
import importlib.util
import os
from pathlib import Path
import sys

PINS={
    'repair_p3_batched_conditional.py':'184b3ec43d205290e0c4e2d8df0237a9df770af38116d62762a8b1f9811f3d76',
    'repair_p3_prepared_conditional_v2.py':'289e93b651411042093a2b987717d87cdb5ed156ab7c736e1ac5cf820670fd4b',
    'repair_p3_bb_gamma_batched.py':'5ae40297221af8158de9dad04d32294e163ef2740f344b055d2acad60e9c0b4a',
    'repair_p3_bb_gamma_difference_v3.py':'546e7ea2f1db21e6fcdd6aed72204d3bbe48c2d818996ab21de5f5ab5215bcf3',
    'repair_p3_gamma_parallel_model.py':'06b384d6a6b6185dc5c0eac747f2529a8491c8ac07a7a14167ea876b05cf8fc2',
}


def require(ok,message):
    if not ok:raise RuntimeError(message)


def context_class(directory=None):
    require(os.environ.get('SLURM_JOB_ID','').isdigit()
            and not os.environ.get('SLURM_ARRAY_TASK_ID'),'Individual compute allocation required')
    directory=Path(__file__).parent if directory is None else Path(directory)
    modules={}
    for name,digest in PINS.items():
        source=directory/name
        require(hashlib.sha256(source.read_bytes()).hexdigest()==digest,'Checked conditional source differs: '+name)
        spec=importlib.util.spec_from_file_location('gamma_conditional_'+source.stem,source)
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module
        spec.loader.exec_module(module);modules[name]=module
    import numpy as np
    base=modules['repair_p3_batched_conditional.py']
    gamma=modules['repair_p3_bb_gamma_batched.py']
    scalar=modules['repair_p3_bb_gamma_difference_v3.py']
    snapshots=modules['repair_p3_gamma_parallel_model.py']

    class GammaBatchedConditional(base.BatchedConditional):
        def __init__(self,model,p,rows):
            rows=tuple(rows)
            require(model.G==8 and model.Z==7 and model.size==28 and 0<len(rows)<=256
                    and all(0<=r.n<=5000 for r in rows),'Prepared synthetic scope exceeded')
            self.input_parameters=np.asarray(p)
            require(self.input_parameters.shape==(28,) and self.input_parameters.dtype==np.dtype('float64')
                    and np.isfinite(self.input_parameters).all(),'Finite original parameter vector required')
            self.parameter_bytes=self.input_parameters.tobytes()
            self.row_snapshots=snapshots.row_signatures(rows,np)
            self.input_rows=rows
            super().__init__(model,p,rows)
            self.fixed_parameters=np.asarray(p).copy();self.fixed_parameters.setflags(write=False)
            self.constant_snapshots={key:self.constants(d) for key,d in self.prepared.rows.items()}
            # The inherited constructor makes static tables; all array storage
            # becomes read-only. Node evaluations never use their BB prefixes.
            def freeze(value):
                if isinstance(value,np.ndarray):value.setflags(write=False)
                elif isinstance(value,dict):
                    for v in value.values():freeze(v)
                elif isinstance(value,(list,tuple)):
                    for v in value:freeze(v)
            for d in self.prepared.rows.values():
                for key,value in d.items():
                    if key not in ('row','signature'):freeze(value)
            self.gamma_BB_calls=0
            self.support_calls={'observed':0,'tails':0,'selected':0}

        def constants(self,d):
            return (d['n'],d['lo'],float(d['r']).hex(),float(d['hom']).hex(),
                    tuple((name,d[name].dtype.str,d[name].shape,d[name].tobytes())
                          for name in ('tails','logprior','hom1','hom2')))

        def prefixes(self,eta,d):
            # One node vector suffices; inherited direct-selection branching
            # slices it along the same node axis as all conditional outputs.
            return (eta,)

        def bb(self,k,d,prefixes,need_hessian):
            require(len(prefixes)==1,'Unexpected Gamma node payload')
            self.gamma_BB_calls+=1
            if len(k)==1 and k[0]==d['row'].a:name='observed'
            elif np.array_equal(k,d['tails']):name='tails'
            else:
                require(np.array_equal(k,np.arange(d['lo'],d['n']-d['lo']+1)),
                        'Unknown integer selection support')
                name='selected'
            self.support_calls[name]+=1
            return gamma.bb_terms(k,d['n'],prefixes[0],d['r'],scalar,need_hessian)

        def terms_batch(self,row,etas,linear_r,p,need_hessian=True):
            require(self.input_parameters.tobytes()==self.parameter_bytes
                    and self.fixed_parameters.tobytes()==self.parameter_bytes
                    and self.prepared.p.tobytes()==self.parameter_bytes
                    and np.asarray(p).dtype==np.dtype('float64')
                    and np.asarray(p).shape==(28,) and np.asarray(p).tobytes()==self.parameter_bytes,
                    'Fixed parameters changed')
            d=self.prepared.rows.get(id(row))
            require(d is not None and d['row'] is row
                    and self.constants(d)==self.constant_snapshots[id(row)],'Prepared constants changed')
            ordinal=self.row_ordinals[id(row)]
            require(snapshots.row_signatures((row,),np)[0]==self.row_snapshots[ordinal],
                    'Fixed row changed')
            try:
                return super().terms_batch(row,etas,linear_r,p,need_hessian)
            except ValueError as exc:
                raise RuntimeError('Conditional context input or arithmetic refusal: '+str(exc)) from exc

        terms=terms_batch

        def callback(self,row,etas,linear_r,need_hessian=True):
            return self.terms_batch(row,etas,linear_r,self.fixed_parameters,need_hessian)

    return GammaBatchedConditional
