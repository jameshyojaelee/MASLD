"""Parallelize only the original three independent null-fit shifts.

The shift body and all post-selection fitting code are compiled from the
hash-pinned original AST. Ordered pool.map preserves the original tie order.
Original arithmetic, bounds, objective, options and polishing remain intact.
"""
import ast
import copy
import hashlib
import multiprocessing
import os
from pathlib import Path
import threading

ORIGINAL_SHA = '41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6'
WORK = None
SHIFT = None


def worker(shift):
    model, start, bounds, objective = WORK
    return SHIFT(model, start, bounds, objective, shift)


def parallel_candidates(model, start, bounds, objective):
    global WORK
    if ((model.G,model.Z,model.size) != (8,7,28) or not 0 < len(model.rows) <= 256
            or not all(0 < r.n <= 5000 for r in model.rows)):
        raise RuntimeError('Small synthetic null-model scope required')
    if (not os.environ.get('SLURM_JOB_ID', '').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID')
            or int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) < 3
            or len(os.sched_getaffinity(0)) < 3 or threading.active_count() != 1):
        raise RuntimeError('Three allocated CPUs and single-thread compute parent required')
    if not all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')):
        raise RuntimeError('One numerical thread per worker required')
    before, row_before = start.tobytes(), repr(model.rows)
    WORK = model, start, bounds, objective
    with multiprocessing.get_context('fork').Pool(3) as pool:
        candidates = pool.map(worker, (0., -.1, .1))
    if start.tobytes() != before or repr(model.rows) != row_before:
        raise RuntimeError('Fixed null optimization inputs changed')
    return candidates


def model_class(base, original):
    global SHIFT
    path = Path(original.__file__)
    if hashlib.sha256(path.read_bytes()).hexdigest() != ORIGINAL_SHA:
        raise RuntimeError('Original null fitter source differs')
    if base.fit_single is not original.Model.fit_single or base.evaluate is not original.Model.evaluate:
        raise RuntimeError('Original fit_single/evaluate required')
    tree = ast.parse(path.read_text())
    model = next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name == 'Model')
    function = next(n for n in model.body if isinstance(n,ast.FunctionDef) and n.name == 'fit_single')
    loops = [n for n in function.body if isinstance(n,ast.For) and isinstance(n.target,ast.Name) and n.target.id == 'shift']
    if len(loops) != 1 or ast.literal_eval(loops[0].iter) != (0., -.1, .1):
        raise RuntimeError('Original shifted candidate loop differs')
    loop = loops[0]
    last = loop.body[-1]
    if not (isinstance(last,ast.Expr) and isinstance(last.value,ast.Call)
            and isinstance(last.value.func,ast.Attribute) and last.value.func.attr == 'append'
            and isinstance(last.value.func.value,ast.Name) and last.value.func.value.id == 'candidates'
            and len(last.value.args) == 1 and isinstance(last.value.args[0],ast.Name) and last.value.args[0].id == 'fit'):
        raise RuntimeError('Original candidate append differs')
    shift_function = ast.FunctionDef(name='isolated_original_shift',
        args=ast.arguments(posonlyargs=[], args=[ast.arg(arg=x) for x in ('self','start','bounds','objective','shift')],
            kwonlyargs=[],kw_defaults=[],defaults=[]),
        body=copy.deepcopy(loop.body[:-1])+[ast.Return(value=ast.Name(id='fit',ctx=ast.Load()))],decorator_list=[])
    altered = copy.deepcopy(function)
    index = function.body.index(loop)
    altered.body[index] = ast.Assign(targets=[ast.Name(id='candidates',ctx=ast.Store())],
        value=ast.Call(func=ast.Name(id='parallel_candidates',ctx=ast.Load()),
                       args=[ast.Name(id=x,ctx=ast.Load()) for x in ('self','start','bounds','objective')],keywords=[]))
    unchanged = all(ast.dump(a,include_attributes=False) == ast.dump(b,include_attributes=False)
                    for i,(a,b) in enumerate(zip(function.body,altered.body)) if i != index)
    if not unchanged or len(function.body) != len(altered.body):
        raise RuntimeError('Code outside candidate loop changed')
    namespace = dict(original.__dict__)
    namespace['parallel_candidates'] = parallel_candidates
    generated = ast.fix_missing_locations(ast.Module(body=[shift_function,altered],type_ignores=[]))
    exec(compile(generated,str(Path(__file__))+'::<original_AST>', 'exec'),namespace)
    SHIFT = namespace['isolated_original_shift']
    result = type('ParallelNullShiftModel',(base,),{'fit_single':namespace['fit_single']})
    result.original_AST_parallel_change = dict(original_source_sha256=ORIGINAL_SHA,
        only_original_shift_loop_replaced=True,post_selection_and_polish_AST_unchanged=True,
        shifts=[0.,-.1,.1],ordered_results=True,workers=3)
    for name in ('fit','fit_marginal','fit_marginal_single','constrained_polish','bounds','initial','cold_start',
                 'informative','predictors','evaluate','sandwich','marginal'):
        if getattr(result,name) is not getattr(base,name):
            raise RuntimeError('Unexpected method override: '+name)
    return result
