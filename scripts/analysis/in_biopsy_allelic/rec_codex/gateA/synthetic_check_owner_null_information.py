#!/usr/bin/env python3
"""Synthetic scientific invariants for new private block assembly; compute only."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import types

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()
    files = [HERE/'owner_null_information.py',Path(__file__),HERE/'run_synthetic_check_owner_null_information.sbatch']
    hashes = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    for p in files[:2]:
        ast.parse(p.read_text())
    if args.inspect:
        print(json.dumps(dict(status='static_metadata_only',sha256=hashes,actual_owner_input_available=False,
              invariants=['row-level archived Jacobian vs new block scatter','orientation and cross covariance','no duplicate inclusion factor','uncapped threshold versus capped Jacobian','complete-F zero-row centering','nuisance basis and units','reduced output whitelist']),indent=2))
        return
    if not os.environ.get('SLURM_JOB_ID','').isdigit():
        raise SystemExit('Compute SLURM required')
    out = REC / ('owner_null_information_synthetic_' + os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    for p in files:
        shutil.copyfile(p,out/p.name)
    try:
        import numpy as np
        import owner_null_information as helper
        from dataclasses import dataclass, replace
        for path, expected in helper.SOURCES.items():
            helper.check(helper.sha(path)==expected,'E_SOURCE')
        # Exact source centering, without importing its loader or likelihood.
        source = helper.HERE.parent/'design_v1.py'
        definitions = [n for n in ast.parse(source.read_text()).body if isinstance(n,ast.ClassDef) and n.name in ('Participant','Design')]
        env = dict(np=np,dataclass=dataclass,replace=replace)
        exec(compile(ast.Module(body=definitions,type_ignores=[]),str(source),'exec'),env)
        # Independent production Jacobian source, extracted without any loader.
        likelihood = helper.HERE.parent/'likelihood_v1.py'
        klass = next(n for n in ast.parse(likelihood.read_text()).body if isinstance(n,ast.ClassDef) and n.name=='Model')
        methods = [n for n in klass.body if isinstance(n,ast.FunctionDef) and n.name in ('__init__','predictors')]
        env_jac = dict(np=np)
        reduced_class = ast.ClassDef(name='Model',bases=[],keywords=[],body=methods,decorator_list=[])
        exec(compile(ast.fix_missing_locations(ast.Module(body=[reduced_class],type_ignores=[])),str(likelihood),'exec'),env_jac)
        rng = np.random.default_rng(20261001)
        raw = rng.normal(size=(20,7))
        raw[:,0] += -8.; raw[:,4] += 5.
        raw[:10,3]=0; raw[10:,3]=1  # cohort constants globally dropped
        raw[0,2]=np.nan; raw[10:,2]=np.nan  # partial indicator in one cohort; all missing in another
        stages = np.tile(np.arange(10)%3,2)
        cohorts = ['toyA']*10+['toyB']*10
        participants = [env['Participant'](f'toy{i:02d}',cohorts[i],int(stages[i]),tuple(raw[i]),1.,0.,0) for i in range(20)]
        design = env['Design'](participants,[],gene_count=0)
        names = ['raw_b','raw_d','raw_ffpe','raw_sex','raw_age10','raw_axis1','raw_axis2']
        masks = ['mask_'+n for n in names]
        indicators = list(design.schema['indicators']); keep=list(design.schema['keep'])
        full_names = names+[names[j]+'_missing' for j in indicators]
        schema = dict(input_mode='raw_complete_F_with_missing_masks',raw_roles=helper.ROLES,raw_count=7,
              raw_columns=names,missing_columns=masks,centering_rule=helper.RULES,indicators=indicators,keep=keep,retained_columns=[full_names[j] for j in keep])
        people=[]
        for i in range(20):
            row=dict(individual_id=f'toy{i:02d}',cohort=cohorts[i],S=int(stages[i]))
            row.update({n:float(raw[i,j]) for j,n in enumerate(names)})
            row.update({n:int(not np.isfinite(raw[i,j])) for j,n in enumerate(masks)})
            people.append(row)
        roster=[{n:r[n] for n in ('individual_id','cohort','S')} for r in people]
        centered = helper._center_private_raw(people,schema)
        z = np.array([[r[n] for n in schema['retained_columns']] for r in centered])
        fill_residual=float(np.max(np.abs(z-design.z)))
        helper.check(fill_residual<1e-12,'E_SYNTH_CENTER_SOURCE')
        gene_rows=[]; weights=[]
        alpha_uncap=[.05,.15,.25,.35]; alpha_cap=[.04,.12,.20,.28]
        for scenario in helper.SCENARIOS:
            for g in range(4):
                gene=f'toy_gene{g}.1'
                gene_rows.append(dict(rho_spec=scenario,gene_id=gene,alpha_gtex_ln_afc=alpha_uncap[g],alpha_plant_capped=alpha_cap[g],expected_included=14.25))
                for i in range(19):  # person19 has no rows but belongs to complete F
                    # W is already inclusion-weighted, with nonzero eta/r cross term.
                    w=.75*np.array([[1+.03*i,.11+g*.01],[.11+g*.01,.7+.02*i]])
                    weights.append(dict(rho_spec=scenario,cohort=cohorts[i],individual_id=f'toy{i:02d}',gene_id=gene,
                         alpha_plant=alpha_cap[g],s_t=1 if (i+g)%2 else -1,S=int(stages[i]),n=40,inclusion_probability=.75,
                         W_eta_eta=w[0,0],W_eta_r=w[0,1],W_r_r=w[1,1]))
        assembled=helper._assemble(people,schema,iter(weights),iter(gene_rows),20,2,76,roster)
        index, cnames,p,_,local,eligibility,_,_=assembled
        row_residual=0.; orientation_residual=0.
        dense={c:np.zeros((8+p+5,8+p+5)) for c in cnames}
        source_rows=[]
        for row in weights:
            if row['rho_spec']!=helper.SCENARIOS[0]:continue
            _,_,v=index[row['individual_id']]
            g=int(row['gene_id'].split('gene')[1].split('.')[0])
            srcrow=types.SimpleNamespace(gene=g,individual=row['individual_id'],stage=v[0],z=tuple(v[3:]),
                    background=v[1],duplication=v[2],orientation=-row['s_t'])
            source_rows.append(srcrow)
        model=env_jac['Model'](source_rows,4)
        parameters=model.initial() if hasattr(model,'initial') else np.zeros(model.size)
        parameters[:4]=alpha_cap
        reference_local={}
        for srcrow,row in zip(source_rows,[r for r in weights if r['rho_spec']==helper.SCENARIOS[0]]):
            _,_,jac,_=model.predictors(srcrow,parameters,need_hessian=False)
            orientation_residual=max(orientation_residual,abs(jac[0,8+1+p]-row['s_t']*srcrow.stage))
            w=np.array([[row['W_eta_eta'],row['W_eta_r']],[row['W_eta_r'],row['W_r_r']]])
            dense[row['cohort']]+=jac.T@w@jac
            coords=[srcrow.gene,4+srcrow.gene,*range(8,model.size)]
            key=(helper.SCENARIOS[0],row['cohort'],row['gene_id'])
            reference_local.setdefault(key,np.zeros((p+7,p+7)))
            j=jac[:,coords]
            reference_local[key]+=j.T@w@j
        for key,block in reference_local.items():row_residual=max(row_residual,float(np.max(np.abs(block-local[key]))))
        helper.check(row_residual<1e-10 and orientation_residual==0,'E_SYNTH_ROW_JAC')
        result=helper._calculate(people,schema,iter(weights),iter(gene_rows),20,2,76,roster)
        direct=helper._projection()(dense,8)
        helper.check(np.isclose(direct['information'],result[0]['I'],rtol=1e-8,atol=1e-8),'E_SYNTH_DENSE_ASSEMBLY')
        helper.check([r['eligible_gene_count'] for r in result[:4]]==[4,3,2,1],'E_SYNTH_UNCAPPED_ELIGIBILITY')
        # Change nuisance units and invertibly rotate the ancestry basis.
        transformed=[dict(r) for r in people]
        for r in transformed:
            r['raw_b']*=100.; r['raw_d']*=.01
            a,b=r['raw_axis1'],r['raw_axis2']
            r['raw_axis1'],r['raw_axis2']=2*a+b,a+3*b
        changed=helper._calculate(transformed,schema,iter(weights),iter(gene_rows),20,2,76,roster)
        basis_error=max(abs(a['I']-b['I']) for a,b in zip(result,changed))
        helper.check(basis_error<1e-6,'E_SYNTH_UNITS_BASIS')
        # Removing a zero-row person changes complete-F centering. It cannot be
        # silently centered on only observed rows; use identical frozen schema.
        subset=helper._center_private_raw(people[:-1],schema)
        zero_row_center_change=max(abs(subset[i]['S_tilde']-centered[i]['S_tilde'])
                                   +abs(subset[i]['b_tilde']-centered[i]['b_tilde']) for i in range(19))
        helper.check(zero_row_center_change>1e-6,'E_SYNTH_ZERO_ROW_ROSTER')
        substituted=[dict(r) for r in people]
        substituted[-1]['individual_id']='toy_replacement_without_rows'
        refused=False
        try:
            helper._assemble(substituted,schema,iter(weights),iter(gene_rows),20,2,76,roster)
        except helper.SafeFailure as error:
            refused=error.args==('E_ROSTER_IDENTITY',)
        helper.check(refused,'E_SYNTH_ZERO_ROW_REPLACEMENT')
        no_cross=[dict(r,W_eta_r=0.) for r in weights]
        diag=helper._calculate(people,schema,iter(no_cross),iter(gene_rows),20,2,76,roster)
        cross_cov_effect=max(abs(a['I']-b['I']) for a,b in zip(result,diag))
        helper.check(cross_cov_effect>1e-8,'E_SYNTH_CROSS_COVARIANCE')
        forbidden={'efficient_score_coefficients','matrix','matrices','local','q','z','gene_ids','individual_ids'}
        helper.check(all(not forbidden & set(r) for r in result),'E_SYNTH_OUTPUT_ALLOWLIST')
        for path,expected in helper.SOURCES.items():
            helper.check(helper.sha(path)==expected,'E_SOURCE_CHANGED')
        for path in files:
            helper.check(helper.sha(path)==hashes[path.name]
                         and helper.sha(out/path.name)==hashes[path.name],'E_PREPARED_SOURCE_CHANGED')
        # Do not save toy matrices or covariates: scalar scientific checks only.
        payload=dict(status='synthetic_invariants_passed',seed=20261001,input_kind='synthetic only',sha256=hashes,
              source_dependency_sha256={p.name:h for p,h in helper.SOURCES.items()},
              environment=dict(python=platform.python_version(),numpy=np.__version__,platform=platform.platform(),
                   threads={n:os.environ.get(n) for n in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')}),
              checks=dict(source_centering_max_abs=fill_residual,independent_row_block_max_abs=row_residual,
                 orientation_residual=orientation_residual,units_basis_information_max_abs=basis_error,
                 zero_row_roster_centering_change=zero_row_center_change,cross_covariance_information_change=cross_cov_effect,
                 zero_row_replacement_refused=True,
                 uncapped_threshold_counts=[4,3,2,1]),owner_private_input_read=False,calibration_or_theta_decision=False)
        (out/'summary.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
        (out/'COMPLETE').write_text('synthetic checks complete\n')
    except Exception:
        (out/'failure.json').write_text(json.dumps(dict(status='synthetic_check_failed',private_values_or_matrices_emitted=False))+'\n')
        raise SystemExit(1) from None


if __name__=='__main__':main()
