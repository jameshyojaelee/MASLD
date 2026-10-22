"""Existing count-prefix conditional arithmetic in one shallow G5 fixture.

No conditional formula, integration rule, parameter or statistical test changes.
This small context is separate from the current G8 model and Gamma helpers.
"""
import hashlib
import importlib.util
import os
from pathlib import Path
import sys

BASE_SHA='bf8c1354ea6f4b463370dbe102f080bf18635ebab60106f4adc39e802239948b'


def context_class(directory=None):
    if not os.environ.get('SLURM_JOB_ID','').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID'):
        raise RuntimeError('Individual compute allocation required')
    directory=Path(__file__).parent if directory is None else Path(directory)
    source=directory/'repair_p3_batched_conditional_g5.py'
    if hashlib.sha256(source.read_bytes()).hexdigest()!=BASE_SHA:
        raise RuntimeError('Checked G5 prefix source differs')
    spec=importlib.util.spec_from_file_location('G5_prefix_base',source);base=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=base;spec.loader.exec_module(base)
    import numpy as np

    class PrefixConditional(base.BatchedConditional):
        def __init__(self,model,p,rows):
            rows=tuple(rows)
            if not (model.G==5 and model.Z==1 and model.size==16 and len(rows)==500
                    and all(20<=r.n<=36 for r in rows)):
                raise RuntimeError('Separate 500-row shallow G5 scope required')
            super().__init__(model,p,rows)
            self.fixed_parameters=np.asarray(p).copy();self.fixed_parameters.setflags(write=False)
            self.prefix_bb_calls=0

        def bb(self,k,d,prefixes,need_hessian):
            self.prefix_bb_calls+=1
            return super().bb(k,d,prefixes,need_hessian)

        def callback(self,row,etas,linear_r,need_hessian=True):
            return self.terms_batch(row,etas,linear_r,self.fixed_parameters,need_hessian)

    return PrefixConditional
