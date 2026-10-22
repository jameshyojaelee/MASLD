"""Synthetic native adaptive GH with scoped, checked Gamma BB arithmetic.

The original zero-variance, null-fit, optimizer, mode-search algorithm, bounds,
starts, ties and restart methods remain intact. Positive integrals use the
original requested GH order; there is no automatic order selection or fallback.
This source alone does not admit fitting, calibration or human-data analysis.
"""
import hashlib
import os
from pathlib import Path
import threading
import time
import traceback

ORIGINAL_SHA='41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6'
GAMMA_SHA='546e7ea2f1db21e6fcdd6aed72204d3bbe48c2d818996ab21de5f5ab5215bcf3'
SNAPSHOT_SHA='06b384d6a6b6185dc5c0eac747f2529a8491c8ac07a7a14167ea876b05cf8fc2'
INHERITED=('fit','fit_single','constrained_polish','fit_marginal','fit_marginal_single',
           'bounds','initial','cold_start','informative','predictors','evaluate','sandwich')


class ModelIntegrityError(RuntimeError):
    """Invalid evaluator state, outside the scientific refit failure count."""


def require(ok,message):
    if not ok:raise ModelIntegrityError(message)


def native_gamma_model_class(private,gamma,snapshots,*,event_sink=None):
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    for module,expected in ((private,ORIGINAL_SHA),(gamma,GAMMA_SHA),(snapshots,SNAPSHOT_SHA)):
        require(hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()==expected,
                'Reviewed numerical source differs')
    import numpy as np
    original_bb=private.bb_terms
    gamma_bb=gamma.stable_bb_terms
    require(original_bb.__name__=='bb_terms' and Path(original_bb.__code__.co_filename).resolve()
            ==Path(private.__file__).resolve(),'Original parent BB required')
    require(gamma_bb.__name__=='stable_bb_terms' and Path(gamma_bb.__code__.co_filename).resolve()
            ==Path(gamma.__file__).resolve(),'Checked Gamma function required')
    callables=(private.conditional_terms,private.Model.predictors,private.Model.evaluate,private.Model.marginal)
    busy=False

    class NativeGammaModel(private.Model):
        def marginal(self,p,v,nodes=20,with_hessian=False,gene_terms=None):
            nonlocal busy
            self.last_integration_failure=None
            self.last_evaluated_point=None
            self.last_integration_reports=[]
            self.last_total_error_allowance=None
            require(private.bb_terms is original_bb and callables==(
                private.conditional_terms,private.Model.predictors,private.Model.evaluate,private.Model.marginal),
                'Original parent callables changed')
            if v<=0:
                # Preserve the original negative-v refusal and NaN variance
                # Hessian at the zero boundary, without a Gamma BB patch.
                return super().marginal(p,v,nodes,with_hessian,gene_terms)
            require(gamma.stable_bb_terms is gamma_bb,'Checked Gamma callable changed')
            require(self.G==8 and self.Z==7 and self.size==28 and 0<len(self.rows)<=256
                and all(r.gene in range(8) and 0<=r.n<=5000 for r in self.rows),
                'Prepared small-synthetic scope exceeded')
            require(threading.active_count()==1 and not busy
                and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))>=1,
                'Single-thread nonreentrant evaluator required')
            require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS',
                'MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'Numerical threads must be fixed to one')
            require(np.asarray(p).shape==(28,) and np.asarray(p).dtype==np.dtype('float64')
                and np.isfinite(p).all() and np.isfinite(v) and v>0,'Finite float64 point required')
            require(isinstance(nodes,int) and nodes in (20,40,80,160),
                    'Native quadrature order outside the checked family')
            input_parameters=np.asarray(p);input_bytes=input_parameters.tobytes()
            input_writeable=input_parameters.flags.writeable
            readonly=input_parameters.copy();readonly.setflags(write=False)
            signatures=snapshots.row_signatures(self.rows,np)
            call=getattr(self,'marginal_call_count',0)+1;self.marginal_call_count=call
            self.last_evaluated_point=dict(call=call,parameters=readonly.copy(),variance=float(v),
                native_order=nodes,with_hessian=with_hessian)
            began=time.perf_counter();busy=True
            def unchanged():
                return (private.bb_terms is gamma_bb and gamma.stable_bb_terms is gamma_bb
                    and callables==(private.conditional_terms,private.Model.predictors,private.Model.evaluate,private.Model.marginal)
                    and snapshots.row_signatures(self.rows,np)==signatures
                    and readonly.tobytes()==input_bytes and not readonly.flags.writeable
                    and input_parameters.tobytes()==input_bytes and input_parameters.flags.writeable==input_writeable)
            try:
                private.bb_terms=gamma_bb
                result=super().marginal(readonly,v,nodes,with_hessian,gene_terms)
                require(unchanged(),'Positive native-Gamma evaluator state changed')
                if event_sink:
                    try:
                        event_sink(dict(event='marginal_evaluation',call=call,variance=float(v),
                            native_order=nodes,with_hessian=with_hessian,loglik=float(result[0]),
                            seconds=time.perf_counter()-began,original_mode_search_algorithm=True))
                    except Exception as sink_error:
                        raise ModelIntegrityError('Event callback failed after native-Gamma evaluation') from sink_error
                    require(unchanged(),'Event callback changed native-Gamma evaluator state')
                return result
            except Exception as exc:
                intact=unchanged()
                self.last_integration_failure=dict(call=call,parameters=readonly.copy(),variance=float(v),
                    native_order=nodes,with_hessian=with_hessian,error_type=type(exc).__name__,message=str(exc),
                    fixed_state_unchanged=intact,traceback=traceback.format_exc(),no_successful_return=True)
                if not intact:
                    raise ModelIntegrityError('State changed during native-Gamma evaluation') from exc
                raise
            finally:
                private.bb_terms=original_bb;busy=False

    for name in INHERITED:
        require(getattr(NativeGammaModel,name) is getattr(private.Model,name),'Original inherited method changed: '+name)
    return NativeGammaModel
