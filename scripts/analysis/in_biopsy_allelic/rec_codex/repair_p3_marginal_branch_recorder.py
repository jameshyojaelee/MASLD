#!/usr/bin/env python3
"""Terminal branch receipts only; original fit inputs/returns/exceptions unchanged.

No numerical imports. Actual-model installation requires an admitted synthetic
compute context. The stdlib fake check verifies instrumentation mechanics only.
"""
from contextlib import contextmanager
import functools
import hashlib
import inspect
import json
import math
import numbers
import os
from pathlib import Path
import time
import traceback
from types import MethodType

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
ORIGINAL_SHA='41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6'
NAME='fit_marginal_single'


class BranchRecordingError(RuntimeError):
    """Instrumentation completeness failure; never a numerical fit failure."""


def snapshot(value):
    """Detached JSON snapshot, preserving finite float round-trip values."""
    if value is None or isinstance(value,(str,bool)):return value
    if isinstance(value,numbers.Integral):return int(value)
    if isinstance(value,numbers.Real):
        x=float(value)
        return x if math.isfinite(x) else {'nonfinite_float':x.hex()}
    if isinstance(value,dict):return {str(k):snapshot(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [snapshot(v) for v in value]
    if hasattr(value,'tolist'):return snapshot(value.tolist())
    raise TypeError('Unsupported receipt value: '+type(value).__name__)


class _Recorder:
    def __init__(self,model,expected_function,*,max_calls=16,max_bytes=512*1024,sink=None):
        if NAME in vars(model):raise BranchRecordingError('Existing instance method override; refuse nested installation')
        original=getattr(model,NAME)
        if getattr(original,'__func__',None) is not expected_function:
            raise BranchRecordingError('Original inherited method identity differs')
        signature=inspect.signature(original)
        if list(signature.parameters)!=['start','fixed_v'] or any(p.default is not inspect.Parameter.empty
                for p in signature.parameters.values()):
            raise BranchRecordingError('Expected exact archived bound signature(start,fixed_v)')
        if not 1<=max_calls<=64 or not 1024<=max_bytes<=1024*1024:
            raise BranchRecordingError('Bounded receipt capacity required')
        self.model,self.original,self.expected_function=model,original,expected_function
        self.max_calls,self.max_bytes,self.sink=max_calls,max_bytes,sink
        self.records,self.recording_errors=[],[]
        self.calls,self.omitted,self.saved_bytes=0,0,0
        self.active_operation=None;self.operation_count=0;self.operation_calls=0
        self.restored=False
        recorder=self
        @functools.wraps(expected_function)
        def wrapped(instance,*args,**kwargs):
            # Original call style and objects are forwarded verbatim. Metadata
            # failures are stored separately and cannot become fit failures.
            recorder.calls+=1;recorder.operation_calls+=1
            ordinal=recorder.calls
            record=dict(invocation=ordinal,operation=recorder.active_operation,
                operation_id=recorder.operation_count if recorder.active_operation else None,
                operation_invocation=recorder.operation_calls if recorder.active_operation else None,
                original_function_module=expected_function.__module__,
                original_function_qualname=expected_function.__qualname__)
            try:
                bound=signature.bind(*args,**kwargs)
                start=bound.arguments['start'];fixed=bound.arguments['fixed_v']
                record.update(start=snapshot(start),fixed_v=snapshot(fixed),
                    forwarded_positional_count=len(args),forwarded_keyword_names=list(kwargs))
            except Exception as exc:recorder._problem('input_snapshot',exc,ordinal)
            began=time.perf_counter()
            try:
                result=original(*args,**kwargs)
            except BaseException as exc:
                record.update(status='raised',elapsed_seconds=time.perf_counter()-began,
                    exception_module=type(exc).__module__,exception_type=type(exc).__qualname__,
                    exception_message=str(exc),exception_traceback=''.join(traceback.format_exception(type(exc),exc,exc.__traceback__,limit=8)))
                recorder._emit(record)
                raise  # The original exception object and class are preserved.
            record.update(status='returned',elapsed_seconds=time.perf_counter()-began)
            try:
                record['terminal']={k:snapshot(result[k]) for k in ('parameters','v','loglik','gradient','converged')}
                record['terminal_float_hex']={k:float(result[k]).hex() for k in ('v','loglik')}
                record['terminal_array_metadata']={k:dict(shape=list(getattr(result[k],'shape',())),
                    dtype=str(getattr(result[k],'dtype','Python sequence'))) for k in ('parameters','gradient')}
            except Exception as exc:recorder._problem('terminal_snapshot',exc,ordinal)
            recorder._emit(record)
            return result  # No copy, wrapping or modification of the fit return.
        self.wrapper=MethodType(wrapped,model)
        setattr(model,NAME,self.wrapper)
        if inspect.signature(getattr(model,NAME))!=signature:
            delattr(model,NAME)
            raise BranchRecordingError('Instrumented introspection signature differs')

    def _problem(self,phase,exc,invocation):
        if len(self.recording_errors)<16:
            self.recording_errors.append(dict(phase=phase,invocation=invocation,
                exception_type=type(exc).__name__,message=str(exc)[:1024]))

    def _emit(self,record):
        try:
            body=json.dumps(record,allow_nan=False)
            if len(self.records)>=self.max_calls or self.saved_bytes+len(body.encode())>self.max_bytes:
                self.omitted+=1
                raise BranchRecordingError('Terminal receipt capacity exhausted; original fit remains unchanged')
            # Sink receives a separate JSON copy, so it cannot mutate memory receipts.
            self.records.append(json.loads(body));self.saved_bytes+=len(body.encode())
            if self.sink is not None:self.sink(json.loads(body))
        except Exception as exc:self._problem('receipt_storage_or_sink',exc,record['invocation'])

    @contextmanager
    def operation(self,label):
        """Caller labels one ORIGINAL model.fit_marginal invocation; no fit wrapper.

        Invocation1/2/3 can then be interpreted as warm/cold/boundary restart
        under the guarded source. Direct single-fit calls must not be relabelled.
        """
        if self.active_operation is not None or self.restored:
            raise BranchRecordingError('Nested operation or restored recorder')
        if not isinstance(label,str) or not 1<=len(label)<=128:
            raise BranchRecordingError('Compact operation label required')
        self.active_operation=label;self.operation_count+=1;self.operation_calls=0
        try:yield self
        finally:self.active_operation=None;self.operation_calls=0

    def receipt(self):
        return dict(records=json.loads(json.dumps(self.records)),invocations_seen=self.calls,
            invocations_omitted=self.omitted,recording_errors=list(self.recording_errors),
            complete=(not self.omitted and not self.recording_errors and len(self.records)==self.calls),
            terminal_only=True,optimizer_selection_not_changed=True,
            instance_method_instrumented=True,class_method_unchanged=getattr(type(self.model),NAME) is self.expected_function,
            inherited_instance_identity_unchanged=False,restored=self.restored,
            max_calls=self.max_calls,max_bytes=self.max_bytes,saved_receipt_bytes=self.saved_bytes,
            no_independent_KKT_recomputation=True,no_global_maximum_claim=True)

    def require_complete(self):
        if not self.receipt()['complete']:
            raise BranchRecordingError('Branch evidence incomplete; do not claim a complete branch history')

    def restore(self):
        if self.restored:return
        if vars(self.model).get(NAME) is not self.wrapper:
            raise BranchRecordingError('Concurrent instance-method replacement; refuse overwrite')
        delattr(self.model,NAME)
        if getattr(self.model,NAME).__func__ is not self.expected_function:
            raise BranchRecordingError('Original method did not restore')
        self.restored=True

    def __enter__(self):return self
    def __exit__(self,*exc):self.restore()


def install_recorder(model,expected_function,*,synthetic_admission,max_calls=16,max_bytes=512*1024,sink=None):
    """Compute-only installation on the already admitted exact synthetic fixture.

    Caller must independently guard the fixture/derived-draw and count receipt;
    this helper checks its declared scope, not scientific input ancestry.
    """
    if not os.environ.get('SLURM_JOB_ID','').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID'):
        raise BranchRecordingError('Individual synthetic compute allocation required')
    if (model.G,model.Z,model.size,len(model.rows))!=(8,7,28,173):
        raise BranchRecordingError('Only exact small synthetic G8/Z7/P28/173-row scope admitted')
    if synthetic_admission.get('fixture_root')!=str(FIXTURE) or synthetic_admission.get('participants')!=60:
        raise BranchRecordingError('Already admitted supplied synthetic fixture declaration required')
    digest=synthetic_admission.get('count_vector_sha256','')
    if len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):
        raise BranchRecordingError('Guarded synthetic count receipt required')
    source=Path(inspect.getsourcefile(expected_function)).resolve()
    if not source.is_relative_to(REC.resolve()) or source.parent.name!='executed_sources' or hashlib.sha256(source.read_bytes()).hexdigest()!=ORIGINAL_SHA:
        raise BranchRecordingError('Exact archived original fit method source required')
    recorder=_Recorder(model,expected_function,max_calls=max_calls,max_bytes=max_bytes,sink=sink)
    recorder.synthetic_admission=dict(synthetic_admission)
    return recorder


def fake_check():
    """Stdlib mechanics only. Actual-model numerical verification is not run."""
    returned={'parameters':[1.], 'v':0., 'loglik':-1., 'gradient':[2.,3.], 'converged':True}
    raised=ValueError('original fake exception')
    class Fake:
        def __init__(self):self.calls=[]
        def fit_marginal_single(self,start,fixed_v):
            self.calls.append((start,fixed_v))
            if start is raised:raise raised
            return returned
    model=Fake();start=[4.];fixed=object();original=Fake.fit_marginal_single
    # A non-serializable fixed_v is deliberately used to verify that recording
    # errors cannot change forwarding or returned-object identity.
    with _Recorder(model,original) as recorder:
        result=model.fit_marginal_single(start=start,fixed_v=fixed)
        assert result is returned and model.calls[-1][0] is start and model.calls[-1][1] is fixed
        try:_Recorder(model,original)
        except BranchRecordingError:pass
        else:raise AssertionError('Nested installation accepted')
        try:model.fit_marginal_single(raised,None)
        except ValueError as exc:assert exc is raised
        else:raise AssertionError('Original exception suppressed')
        assert not recorder.receipt()['complete']
    assert model.fit_marginal_single.__func__ is original and NAME not in vars(model)
    with _Recorder(model,original) as recorder:
        with recorder.operation('fake_outer_call'):
            assert model.fit_marginal_single(start,None) is returned
        recorder.require_complete()
        assert recorder.records[0]['terminal']['gradient']==[2.,3.]
    assert model.fit_marginal_single.__func__ is original
    def broken_sink(record):raise RuntimeError('fake sink failure')
    with _Recorder(model,original,sink=broken_sink) as recorder:
        assert model.fit_marginal_single(start,None) is returned
        try:recorder.require_complete()
        except BranchRecordingError:pass
        else:raise AssertionError('Sink failure did not mark evidence incomplete')
    return dict(fake_forwarding_and_return_identity=True,fake_exception_identity=True,
        nested_install_refused=True,restore_passed=True,metadata_failure_preserves_original_call=True,
        sink_failure_preserves_fit_and_marks_incomplete=True,actual_model_compute_verification='not run',
        numerical_imports=False)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fake-check',required=True,action='store_true')
    parser.parse_args()
    print(json.dumps(fake_check(),sort_keys=True))
