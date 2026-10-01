"""Independent Decimal BB/mixture r-score at two frozen synthetic contexts."""
import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
GRID = REC / 'codex_positive_variance_grid_21999349'
ROLE = REC / 'bfix_p3_21998519'
FIX = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
ARCHIVE = ROLE / 'executed_sources'
CAP = 32 << 20
FIXED = {
    GRID / 'protocol.json': '40cafedf06d707d8a6dae53960d006fa10352c863332d9fbb7bc3c40c84cb461',
    GRID / 'summary.json': '94ccdd90a2a116db94c5c0243643c8b9cd31c089de0f55c28db1d97b2a564d3b',
    GRID / 'executed_sources/codex_positive_variance_grid.py': 'dc49d1d807c41e2b81f65ce4c029bd01f278d3a169f99ab6dd44fd21aa462de1',
}
COUNT_SHA = '5f4a1aa455cd09b0d5f31889021ec88b145dcbb50c8e5b600343ea9a3a992691'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def guards():
    result = dict(FIXED)
    for path, digest in result.items():
        require(sha(path) == digest, 'Frozen completed grid differs: '+str(path))
    protocol = json.loads((GRID / 'protocol.json').read_text())
    result.update({Path(p): h for p,h in protocol['source_guards'].items()})
    result.update({GRID / 'executed_sources' / p: h for p,h in protocol['private_source_guards'].items()})
    for path,digest in result.items():
        require(sha(path)==digest,'Frozen source/input differs: '+str(path))
    require(protocol['count_sha256']==COUNT_SHA and protocol['original_source_sha256']==
            '41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6','Wrong archived model/count')
    return result,protocol


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args()
    frozen,protocol=guards()
    plan=dict(contexts=['observed_null','failed_proposal'],gene=2,variance=.5,rows=24,
        digits=[50,70],nodes=[80,160],fd_steps=['0.0003','0.0001','0.00003'],
        fixed_grid=True,independent_decimal_BB_mixture=True,no_fit=True,no_MC=True,
        conversion='Decimal.from_float for every archived numeric input; exact decimal FD steps',
        limits='Double GH nodes/weights; fixed local center; no integration error bound, global mode or calibration claim')
    if args.inspect:
        print(json.dumps(dict(plan=plan,guards={str(p):h for p,h in frozen.items()}),indent=2));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Compute-only allocation required')
    out=REC/('codex_decimal_r_precision_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    own=(Path(__file__),Path(__file__).with_name('run_codex_decimal_r_precision.sbatch'))
    frozen.update({p:sha(p) for p in own})
    for p in own:shutil.copyfile(p,out/p.name)
    started=time.monotonic();phase='imports'
    def save(name,value):
        text=json.dumps(value,indent=2,allow_nan=False)+'\n'
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(text.encode())<=CAP-65536,'Output cap')
        (out/(name+'.json')).write_text(text)
    try:
        # Numerical imports and all data reconstruction occur only under SLURM.
        from decimal import Decimal,localcontext
        import numpy as np
        import scipy
        require((platform.python_version(),np.__version__,scipy.__version__)==
                tuple(protocol['environment'][k] for k in ('python','numpy','scipy')),'Original numeric environment differs')
        require(not any(n in sys.modules for n in ('likelihood_v1','bfix_v1','design_v1','tests_v1')),'Already imported numerical model')
        sys.path.insert(0,str(ARCHIVE))
        loader=importlib.import_module('bfix_v1');likelihood=importlib.import_module('likelihood_v1')
        receipt=json.loads((ROLE/'receipt.json').read_text())
        for name in ('likelihood_v1','bfix_v1','design_v1','tests_v1'):
            p=Path(sys.modules[name].__file__).resolve()
            require(p.parent==ARCHIVE.resolve() and sha(p)==receipt['source_sha256'][name+'.py'],'Archived import differs')
        require(loader.FIXTURE.resolve()==FIX.resolve(),'Fixture escaped approved synthetic source')
        phase='exact_draw_reconstruction'
        design,settings=loader.load();model=design.model()
        null=json.loads((ROLE/'observed_null.json').read_text())
        require(null['converged'] and settings['B']==50 and len(model.rows)==173
                and model.G==8 and model.Z==7 and model.size==28,'Fixture contract differs')
        p0=np.asarray(null['parameters'],float)
        require(p0.tolist()==protocol['parameter_vectors']['observed_null'],'Archived null parameters differ')
        require([(r['individual_id'],r['gene_id']) for r in loader.table('rows.tsv')]==
                [(r.individual,settings['gene_order'][r.gene]) for r in model.rows],'Synthetic row order differs')
        uniforms=np.load(FIX/'draws/p3_uniforms.npy',allow_pickle=False,mmap_mode='r')
        require(uniforms.shape==(50,173,2),'Uniform shape differs')
        drawn=[]
        for row,(u1,u2) in zip(model.rows,uniforms[0]):
            eta,r,_,_=model.predictors(row,p0)
            drawn.append(likelihood.conditional_draw(row,eta,r,u1,u2))
        count_sha=hashlib.sha256(np.asarray([r.a for r in drawn],dtype='<i8').tobytes()).hexdigest()
        require(count_sha==COUNT_SHA,'Exact173 count vector differs')
        rows=[r for r in drawn if r.gene==2]
        require(len(rows)==24 and all(r.a in r.allowed for r in drawn),'Selected synthetic rows differ')
        completed=json.loads((GRID/'summary.json').read_text())
        contexts=[next(c for c in completed['contexts'] if c['parameter']==label and c['v']==.5)
                  for label in plan['contexts']]
        anchors={c['parameter']:next(g for g in c['genes'] if g['gene']==2)['production'][-1] for c in contexts}
        require(all(a['nodes']==80 for a in anchors.values()),'Saved fixed80 anchor missing')
        save('protocol',dict(**plan,count_sha256=count_sha,guards={str(p):h for p,h in frozen.items()},
            saved_centers={k:{j:v[j] for j in ('mode','scale')} for k,v in anchors.items()},
            archived_parameters=protocol['parameter_vectors'],EPS=likelihood.EPS,
            environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform())))
        # Keep every row input; no truth class or external outcome is used.
        save('synthetic_gene_rows',[{k:(list(v) if isinstance(v,tuple) else v) for k,v in r.__dict__.items()} for r in rows])
        def D(x):return Decimal.from_float(float(x))
        def logistic(x):
            if x>=0:return Decimal(1)/(Decimal(1)+(-x).exp())
            e=x.exp();return e/(Decimal(1)+e)
        def pi_value():
            def atan(z):
                term=z;total=z;j=1
                while True:
                    term=-term*z*z;next_total=total+term/Decimal(2*j+1)
                    if next_total==total:return total
                    total=next_total;j+=1
            return 16*atan(Decimal(1)/5)-4*atan(Decimal(1)/239)
        def bb(n,mean,rho,derivative=False,u=None):
            """Full BB probabilities by positive Pochhammer recurrence.

            d(log mass)/dr=(D_a[k]+D_b[n-k]-D_t[n])*d(log t)/dr.
            Fixed homo distributions never call production BB arithmetic.
            """
            one=Decimal(1)
            if rho==0:
                require(not derivative,'Unexpected binomial derivative')
                if mean==0:return [one]+[Decimal(0)]*n,None
                if mean==1:return [Decimal(0)]*n+[one],None
                mass=(one-mean)**n;values=[mass]
                for k in range(n):
                    mass*=Decimal(n-k)/Decimal(k+1)*mean/(one-mean);values.append(mass)
                return values,None
            require(0<mean<1 and 0<rho<1,'Invalid BB probability inputs')
            t=(one-rho)/rho;a=t*mean;b=t*(one-mean)
            mass=one
            for j in range(n):mass*=(b+j)/(t+j)
            values=[mass]
            for k in range(n):
                mass*=Decimal(n-k)/Decimal(k+1)*(a+k)/(b+n-k-1);values.append(mass)
            if not derivative:return values,None
            da=[Decimal(0)];db=[Decimal(0)];dt=Decimal(0)
            for j in range(n):
                da.append(da[-1]+a/(a+j));db.append(db[-1]+b/(b+j));dt+=t/(t+j)
            scores=[(da[k]+db[n-k]-dt)*u for k in range(n+1)]
            return values,scores
        records=[]
        for digits in (50,70):
            with localcontext() as ctx:
                ctx.prec=digits
                one=Decimal(1);epsilon=D(likelihood.EPS);c=one-2*epsilon
                pi=pi_value();fixed_hom={}
                for i,row in enumerate(rows):
                    fixed_hom[i]=(bb(row.n,D(row.e),D(row.phi))[0],bb(row.n,one-D(row.e),D(row.phi))[0])
                for label in plan['contexts']:
                    p=[D(x) for x in protocol['parameter_vectors'][label]];variance=D(.5)
                    mode=D(anchors[label]['mode']);scale=D(anchors[label]['scale'])
                    converted=[]
                    for row in rows:
                        z=[D(x) for x in row.z];stage=D(row.stage)
                        eta=p[2]*(one+p[16]*stage+sum((p[17+j]*x for j,x in enumerate(z)),Decimal(0)))
                        eta-=Decimal(row.orientation)*p[24]*stage
                        eta-=Decimal(row.orientation)*D(row.omega)
                        rr=p[10]+p[25]*stage+p[26]*D(row.background)+p[27]*D(row.duplication)
                        pr,pa,px=D(row.p_ref),D(row.p_alt),D(row.p_x)
                        prior=[2*pr*pa,pr*pr+2*pr*px,pa*pa+2*pa*px]
                        if row.orientation!=1:prior=[prior[0],prior[2],prior[1]]
                        lo=max(row.minor,(row.percent*row.n+99)//100)
                        require(row.n>=row.depth and lo<=row.a<=row.n-lo,'Integer selection differs')
                        converted.append((eta,rr,stage,prior,lo))
                    for nodes in (80,160):
                        phase=f'{label}_digits{digits}_nodes{nodes}'
                        x,w=np.polynomial.hermite.hermgauss(nodes)
                        require(np.isfinite(x).all() and np.isfinite(w).all() and (w>0).all(),'Invalid GH rule')
                        fixed_nodes=[(D(xx),D(ww),mode+Decimal(2).sqrt()*scale*D(xx)) for xx,ww in zip(x,w)]
                        def integral(dr):
                            logterms=[];scores=[];normalization=[]
                            for xx,ww,u in fixed_nodes:
                                logprob=Decimal(0);total_score=Decimal(0)
                                for i,row in enumerate(rows):
                                    eta,rr,stage,prior,lo=converted[i]
                                    m=logistic(eta+u*stage);ell=logistic(rr+dr)
                                    rho=epsilon+c*ell;rp=c*ell*(one-ell)
                                    logt_prime=-rp/(rho*(one-rho))
                                    h,hs=bb(row.n,m,rho,True,logt_prime)
                                    hom1,hom2=fixed_hom[i]
                                    mass=[prior[0]*h[k]+prior[1]*hom1[k]+prior[2]*hom2[k] for k in range(row.n+1)]
                                    numerator=[prior[0]*h[k]*hs[k] for k in range(row.n+1)]
                                    den=sum(mass[lo:row.n-lo+1],Decimal(0))
                                    dden=sum(numerator[lo:row.n-lo+1],Decimal(0))
                                    require(mass[row.a]>0 and den>0,'Nonpositive conditional mixture')
                                    logprob+=mass[row.a].ln()-den.ln()
                                    total_score+=numerator[row.a]/mass[row.a]-dden/den
                                    normalization.append(abs(sum(h,Decimal(0))-one))
                                logterms.append(ww.ln()+xx*xx+logprob-u*u/(2*variance))
                                scores.append(total_score)
                            shift=max(logterms);weights=[(a-shift).exp() for a in logterms]
                            mass=sum(weights,Decimal(0));score=sum((w*s for w,s in zip(weights,scores)),Decimal(0))/mass
                            ll=(Decimal(2).sqrt()*scale).ln()+shift+mass.ln()-(2*pi*variance).ln()/2
                            return ll,score,max(normalization)
                        base_ll,base_score,norm_error=integral(Decimal(0));fds=[]
                        for step in ('0.0003','0.0001','0.00003'):
                            eps=Decimal(step);plus=integral(eps);minus=integral(-eps)
                            fd=(plus[0]-minus[0])/(2*eps)
                            fds.append(dict(step=step,plus_LL=str(plus[0]),minus_LL=str(minus[0]),score_FD=str(fd),
                                FD_minus_analytic=str(fd-base_score),max_BB_normalization_error=str(max(plus[2],minus[2]))))
                        rec=dict(context=label,digits=digits,nodes=nodes,LL=str(base_ll),analytic_r2_score=str(base_score),
                            max_BB_normalization_error=str(norm_error),finite_differences=fds,
                            fixed_double_nodes_weights_SHA256=hashlib.sha256(np.asarray([x,w],dtype='<f8').tobytes()).hexdigest(),
                            saved_grid80_LL=anchors[label]['ll'],saved_grid80_r2_score=anchors[label]['score'][10],
                            LL_minus_saved_float80=str(base_ll-D(anchors[label]['ll'])),
                            score_minus_saved_float80=str(base_score-D(anchors[label]['score'][10])))
                        records.append(rec);save(phase,rec);print(phase,'saved',flush=True)
        phase='comparison_and_post_guards'
        comparisons=[]
        with localcontext() as ctx:
            ctx.prec=70
            for label in plan['contexts']:
                for field in ('LL','analytic_r2_score'):
                    for nodes in (80,160):
                        a=next(r for r in records if r['context']==label and r['nodes']==nodes and r['digits']==50)
                        b=next(r for r in records if r['context']==label and r['nodes']==nodes and r['digits']==70)
                        comparisons.append(dict(context=label,field=field,nodes=nodes,digits70_minus50=str(Decimal(b[field])-Decimal(a[field]))))
                    for digits in (50,70):
                        a=next(r for r in records if r['context']==label and r['nodes']==80 and r['digits']==digits)
                        b=next(r for r in records if r['context']==label and r['nodes']==160 and r['digits']==digits)
                        comparisons.append(dict(context=label,field=field,digits=digits,nodes160_minus80=str(Decimal(b[field])-Decimal(a[field]))))
        for path,digest in frozen.items():require(sha(path)==digest,'Post-run input/source differs: '+str(path))
        require(hashlib.sha256(np.asarray([r.a for r in drawn],dtype='<i8').tobytes()).hexdigest()==COUNT_SHA,'Synthetic counts changed')
        save('summary',dict(status='completed_fixed_precision_diagnostic',records=records,comparisons=comparisons,
            elapsed_seconds=time.monotonic()-started,all_postguards_passed=True,**plan))
    except BaseException:
        (out/'failure.json').write_text(json.dumps(dict(phase=phase,traceback=traceback.format_exc()),indent=2)+'\n');raise


if __name__=='__main__':main()
