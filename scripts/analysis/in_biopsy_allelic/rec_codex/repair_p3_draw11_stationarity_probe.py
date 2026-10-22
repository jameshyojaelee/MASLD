#!/usr/bin/env python3
"""Fixed synthetic null-point Newton probes; no fit or acceptance decision."""
import argparse
from dataclasses import replace
import hashlib
import importlib
import json
import os
from pathlib import Path
import platform
import sys
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
POINT=REC/'repair_p3_gamma_draw_11_22081580'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
SOURCE_DIR=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
PINS={
 'p3_011.json':'1b15b5c46f3d96d770c71cbeb838fed6fcc614964dc0e5379fe9636adf5adacc',
 'draw_counts.json':'9a6ea7d87f798b763e24c046751255948b35b0ce7923ca3df8a413fe39d933c4',
 'protocol.json':'b538f44a0cfccecee08fbc281c69967ee693898aebe85e81127b46715dea1ce4',
 'returned_fit_diagnostic.json':'3d57d38300cc17e7667ed15975e3c378e7db797db2a34aac508b86d53ff3e1e6'}
SOURCES={
 'likelihood_v1.py':'41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6',
 'design_v1.py':'cb3163896f982f39ba94cdbf227a8e7d4c7f85a0a646b26542868e6fa14c9c5e',
 'bfix_v1.py':'324594ce783ba2464921eaee5924065f579060eded6043a48bdf4b7e392bfa56',
 'tests_v1.py':'52cd14d43b649c20347330b64b467dbd81320f10252fe85ce3c538ceef2f092e'}


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def metadata():
    pins={str(POINT/n):h for n,h in PINS.items()}
    pins.update({str(POINT/'executed_sources'/n):h for n,h in SOURCES.items()})
    protocol=json.loads((POINT/'protocol.json').read_text())
    for absolute,digest in protocol['source_and_input_sha256'].items():
        path=Path(absolute).resolve()
        require(path.is_relative_to(REC.resolve()) or path.is_relative_to(FIXTURE.resolve())
            or (path.parent==SOURCE_DIR and path.suffix=='.py'),
            'Input outside synthetic/source scope')
        pins[absolute]=digest
    for absolute,digest in pins.items():require(sha(Path(absolute))==digest,'Pinned input differs: '+absolute)
    require(protocol['draw']==11 and protocol['B']==50,'Exact draw11 required')
    return pins


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args();pins=metadata()
    if args.inspect:
        print(json.dumps({'pins_verified':len(pins),'numeric_imports':False,'no_fitting':True}));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
        'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS',
        'MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'One numerical thread required')
    out=REC/('repair_p3_draw11_stationarity_probe_'+os.environ['SLURM_JOB_ID']);out.mkdir()
    (out/'executed_checker.py').write_bytes(Path(__file__).read_bytes())
    try:
        import numpy as np
        import scipy
        require((platform.python_version(),np.__version__,scipy.__version__)==('3.10.18','2.2.6','1.15.2'),
            'Recorded numerical environment required')
        archive=POINT/'executed_sources';sys.path.insert(0,str(archive))
        loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        for name in SOURCES:
            require(Path(sys.modules[name[:-3]].__file__).resolve()==(archive/name).resolve(),
                'Original source import escaped archive')
        require(Path(loader.FIXTURE).resolve()==FIXTURE.resolve(),'Synthetic fixture differs')
        design,settings=loader.load();template=design.model()
        counts=json.loads((POINT/'draw_counts.json').read_text())
        require(len(template.rows)==len(counts['counts'])==173 and len(design.participants)==60,
            'Synthetic dimensions differ')
        require(counts['full_row_keys']==[[row.individual,settings['gene_order'][row.gene]] for row in template.rows],
            'Synthetic row order differs')
        values=np.asarray(counts['counts'],dtype='<i8')
        require(hashlib.sha256(values.tobytes()).hexdigest()==counts['count_vector_sha256'],
            'Synthetic count vector differs')
        private=sys.modules['likelihood_v1']
        model=private.Model([replace(row,a=int(a)) for row,a in zip(template.rows,values)],gene_count=8)
        saved=json.loads((POINT/'p3_011.json').read_text())['fit']['null']
        p=np.asarray(saved['parameters'],dtype=float);require(p.shape==(28,),'Parameter shape differs')
        bounds=[b if model.informative()[j] else (p[j],p[j]) for j,b in enumerate(model.bounds())]
        def evaluate(point):
            ll,g,h=model.evaluate(point)
            active=np.array([(lo is not None and lo==hi) or
                (lo is not None and point[j]<=lo+1e-10 and g[j]<0) or
                (hi is not None and point[j]>=hi-1e-10 and g[j]>0)
                for j,(lo,hi) in enumerate(bounds)])
            projected=g.copy();projected[active]=0
            require(np.isfinite(ll) and np.isfinite(g).all() and np.isfinite(h).all(),
                'Nonfinite fixed-point calculation')
            return ll,g,h,active,float(np.max(np.abs(projected)))
        ll,g,h,active,maximum=evaluate(p);free=np.flatnonzero(~active)
        require(abs(ll-saved['loglik'])<1e-6,'Saved null likelihood identity differs')
        curvature=-h[np.ix_(free,free)];symmetric=(curvature+curvature.T)/2
        eigenvalues,eigenvectors=np.linalg.eigh(symmetric)
        floor=max(1e-9,float(np.max(np.abs(eigenvalues)))*1e-10)
        regularized=(eigenvectors*np.maximum(eigenvalues,floor))@eigenvectors.T
        proposals={'original_lstsq':np.linalg.lstsq(curvature,g[free],rcond=1e-12)[0],
            'regularized_solve':np.linalg.solve(regularized,g[free])}
        records=[]
        for label,direction in proposals.items():
            step=np.zeros(28);step[free]=direction
            for damping in (1.,.5,.25,.125):
                trial=p+damping*step
                for j,(lo,hi) in enumerate(bounds):trial[j]=np.clip(trial[j],
                    -np.inf if lo is None else lo,np.inf if hi is None else hi)
                tl,tg,th,ta,tm=evaluate(trial)
                records.append({'proposal':label,'damping':damping,'parameters':trial.tolist(),
                    'loglik':tl,'loglik_difference':tl-ll,'gradient':tg.tolist(),
                    'max_projected_gradient':tm,'strict_cutoff_met':tm<1e-6})
                (out/('probe_%02d.json'%len(records))).write_text(
                    json.dumps(records[-1],indent=2,allow_nan=False)+'\n')
        result={'source_and_input_sha256':pins,'checker_sha256':sha(Path(__file__)),
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'count_vector_sha256':counts['count_vector_sha256'],'baseline_loglik':ll,
            'baseline_LL_minus_saved':ll-saved['loglik'],'baseline_max_projected_gradient':maximum,
            'saved_native_converged':saved['converged'],'probes':records,
            'no_proposal_selected':True,'no_fit_result_replaced':True,'no_new_convergence_flag':True,
            'independent_gradient_accuracy_verified':False,'calibration':False}
        for absolute,digest in pins.items():require(sha(Path(absolute))==digest,'Input changed during diagnostic')
        body=json.dumps(result,indent=2,allow_nan=False)+'\n';require(len(body)<1024**2,'Diagnostic cap exceeded')
        (out/'summary.json').write_text(body)
    except BaseException:
        (out/'failure.json').write_text(json.dumps({'traceback':traceback.format_exc(),
            'no_fit_result_replaced':True},indent=2)+'\n');raise


if __name__=='__main__':main()
