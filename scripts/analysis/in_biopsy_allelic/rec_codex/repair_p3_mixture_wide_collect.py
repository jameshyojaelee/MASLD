#!/usr/bin/env python3
"""Stdlib-only reconciliation of the four frozen wider-proposal checks."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

REC=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z')
LAUNCH=REC/'repair_p3_mixture_gh_wide_launch_20261005T154158Z'
POLICY_COLLECTOR=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_mixture_policy_collect.py')
POLICY_COLLECTOR_SHA='39bf973bb579b73222c9ce4cfe7b61cba6837cd3024cb9c18e3df7b3cd9e5991'
KEYS=[('failed0gene2','scale_precision'),('failed0gene2','scale_sqrt'),
      ('positive18','scale_precision'),('positive18','scale_sqrt')]
JOBS=(22041642,22041643,22041644,22041645)
SOURCES={
 'repair_p3_fixed_point_check.py':'c4edba988955b8b922e28a227c4b936061b073df3b2995c60f14785d5565413c',
 'repair_p3_mixture_gh.py':'8a21be3a5d4f7c731b13cdfe5a01165f640be36556fc8227ec6a8c67441828a0',
 'repair_p3_integral.py':'461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3',
 'repair_p3_reference_fast.py':'2eaaf549739f64322d394f0674aadba98f0f2ae6bedda1df56db8e9a8a3949fd',
 'codex_stable_bb_probe.py':'723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488'}
NARROW=[
 (22040801,{'protocol.json':'f5ed0a24122f64662f626c7b68e3c9782caa44e2aa6fe2bb69eacfd12e91ded2',
  'summary.json':'6281b22874dffd55a879164c1ea3f5bfc067e3bfd3fe7f69ead37d4225626b99',
  'independent_baseline.json':'0b97e74d0f9fbdae269048a54ed6b27be9841b884c7ab3219f30ccda27ca6402',
  'candidate_refusal.json':'cc85f5ee6ee37eebed3222ca3c6cc23399953c5ce629048ba1fb8604cbb36bb8'}),
 (22040802,{'protocol.json':'2b5fd8253d9fc945279be24e4e618f5b09c602929b07c2065ef4c8853fdce08b',
  'summary.json':'ea70596bb9e8789d6ae7eca933839d5cbb6f76ae7109b0f6f388c64fba10d60d',
  'independent_baseline.json':'0b97e74d0f9fbdae269048a54ed6b27be9841b884c7ab3219f30ccda27ca6402',
  'candidate_refusal.json':'d878212ebbb225d3ead095e226beca80a15fd433602ea427dab029bb52309d53'}),
 (22040803,{'protocol.json':'151a005b48622ed1218d0b89443742e6cceebf7a9b5b70b7264a1546400d9e8a',
  'summary.json':'7af519663c7267b09f01fcd1c8fab06e28153e6548736f308e2d936900bcc162',
  'independent_baseline.json':'2d4354f212286a42ef1f3f03899675a29c1a971abb408209b043a470ff249d70',
  'candidate_baseline.json':'57634718992e171928ace22310100f32897a73a7a795edd3018f34dd8cd687a7'}),
 (22040804,{'protocol.json':'457d4a0af7f080724cd92bda41ab96d8ad0d68972aa8ce17b81901a3abfd5bbd',
  'summary.json':'e13122fd50d824607f0887eb01b83337f08654fd8162e2a799d409ead59a161e',
  'independent_baseline.json':'c5baa500540131f4b6109ee2186251cfa57bc51af144028f4089bbd94891528f',
  'candidate_refusal.json':'96e4b49f27a880b5c6e677509a12ea2844bd3104d5cffdd62bd0ff59476c52c2'})]


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_policy_collector():
    require(sha(POLICY_COLLECTOR)==POLICY_COLLECTOR_SHA,'Reviewed reconciliation source differs')
    spec=importlib.util.spec_from_file_location('wide_private_policy_reconciliation',POLICY_COLLECTOR)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    # A private module instance reuses the already reviewed stdlib arithmetic.
    # Its own main is not executed and its source file is never changed.
    module.LAUNCH=LAUNCH;module.SOURCES=SOURCES
    return module


def metadata(helper):
    pins,_=helper.read(LAUNCH/'source_hashes.json')
    require(pins==SOURCES,'Frozen wide source manifest differs')
    for name,digest in SOURCES.items():require(sha(LAUNCH/name)==digest,'Frozen wide source differs')
    entries=[json.loads(line) for line in (LAUNCH/'submissions.jsonl').read_text().splitlines()]
    require(tuple(int(e['job_id']) for e in entries)==JOBS,'Exact four wide jobs required')
    for entry,(case,policy) in zip(entries,KEYS):
        require(entry['case']==case and entry['policy']==policy and entry['wide_components'] is True and
                Path(entry['script'])==LAUNCH/f'{case}_{policy}.sbatch','Wide job ownership differs')
    return entries


def wider_geometry(order):
    require(order['wide_components'] is True and order['wide_width_multiplier']==4.,'Wide run option differs')
    geometry=order['proposal_geometry'];components=order['component_contributions']
    require(geometry['wide_components'] is True and geometry['wide_width_multiplier']==4. and
            geometry['original_local_proposals_retained'] is True and
            geometry['components']==[c['proposal'] for c in components], 'Reported wide geometry differs')
    require((len(components)-1)%2==0 and order['prior_policy_nonprior_proposal_count']==len(components)-1,
            'Nonprior wide component family differs')
    residuals=[]
    for i in range(1,len(components),2):
        narrow,wide=components[i]['proposal'],components[i+1]['proposal']
        require(narrow['kind']=='located_maximum' and wide['kind']=='located_maximum_wide' and
                wide['mean_z']==narrow['mean_z'] and wide['base_sd_z']==narrow['sd_z'] and
                wide['width_multiplier']==4. and wide['location']==narrow['location'],
                'A local Gaussian was removed, shifted, or wrongly widened')
        require(all(math.isfinite(v) for v in (narrow['mean_z'],narrow['sd_z'],wide['sd_z'])) and
                narrow['sd_z']>0 and wide['sd_z']>0,'Nonfinite or nonpositive proposal widths')
        residual=wide['sd_z']-4*narrow['sd_z']
        require(residual==0,'Wide Gaussian is not exactly fourfold')
        residuals.append(residual)
    return dict(order=order['order'],z_bound=order['z_bound'],grid_spacing=order['grid_spacing'],
                original_local_proposals=(len(components)-1)//2,added_wide_proposals=(len(components)-1)//2,
                wide_width_relation_max_residual=max(map(abs,residuals),default=0.))


def last_order(record):
    if 'candidate_baseline' in record:return record['candidate_baseline']['integration']['refined']
    orders=[o for o in record['candidate_refusal'].get('integration_reports',[]) or []
            if isinstance(o,dict) and 'component_contributions' in o]
    return orders[-1] if orders else None


def collect(helper,entries):
    # Existing reconciliation retains every wide refusal and compares frozen GL
    # parents. It performs no writes and no numerical-library imports.
    result=helper.collect(entries)
    if not result['ready']:return result
    result['input_hashes'][str(POLICY_COLLECTOR)]=POLICY_COLLECTOR_SHA
    for record,(case,policy),(job,pins) in zip(result['results'],KEYS,NARROW):
        directory=REC/f'repair_p3_fixed_point_{case}_{job}';old={}
        for name,digest in pins.items():
            old[name],actual=helper.read(directory/name)
            require(actual==digest,'Preceding narrow receipt changed')
            result['input_hashes'][str(directory/name)]=digest
        protocol=helper.read(REC/f"repair_p3_fixed_point_{case}_{record['job_id']}"/'protocol.json')[0]
        require(protocol['mixture_wide_components'] is True,
                'Wide protocol option differs')
        for key in helper.FIELDS:require(protocol[key]==old['protocol.json'][key],'Narrow/wide fixed point differs: '+key)
        require(old['protocol.json']['mixture_prior_policy']==policy and old['protocol.json']['mixture_gh'] is True,
                'Preceding narrow policy differs')
        require(old['summary.json']['baseline_only'] is True and old['summary.json']['executed_derivative_axes']==[] and
                old['summary.json']['derivative_checks_completed'] is False,'Preceding narrow summary differs')
        if 'candidate_baseline.json' in old:
            narrow=dict(candidate_baseline=old['candidate_baseline.json'])
            require(old['summary.json']['candidate_integral_returned'] is True,'Narrow accepted baseline mismatch')
            narrow_orders=old['candidate_baseline.json']['integration']['ladder']
        else:
            narrow=dict(candidate_refusal=old['candidate_refusal.json'])
            require(old['summary.json']['candidate_integral_returned'] is False,'Narrow refusal mismatch')
            narrow_orders=[o for o in old['candidate_refusal.json'].get('integration_reports',[]) or []
                           if isinstance(o,dict) and 'component_contributions' in o]
        wide_orders=(record['candidate_baseline']['integration']['ladder'] if 'candidate_baseline' in record else
                     [o for o in record['candidate_refusal'].get('integration_reports',[]) or []
                      if isinstance(o,dict) and 'component_contributions' in o])
        record['wide_geometry_reconciliations']=[wider_geometry(o) for o in wide_orders]
        record['preceding_narrow']=dict(job_id=job,input_hashes=pins,
            candidate_integral_returned=old['summary.json']['candidate_integral_returned'],
            refusal_preserved_at_source='candidate_refusal.json' in old,
            component_reconciliations=[helper.reconcile(o,policy,old['protocol.json']['candidate_settings']) for o in narrow_orders])
        wide_last,narrow_last=last_order(record),last_order(narrow)
        if wide_last and narrow_last:
            record['last_order_wide_minus_narrow_descriptive']=dict(
                wide_order=wide_last['order'],wide_z_bound=wide_last['z_bound'],
                narrow_order=narrow_last['order'],narrow_z_bound=narrow_last['z_bound'],
                logmass_difference=wide_last['logmass']-narrow_last['logmass'],
                gradient_max_difference=helper.max_vector_difference(wide_last['gradient'],narrow_last['gradient']),
                hessian_max_difference=max(helper.max_vector_difference(a,b) for a,b in zip(wide_last['hessian'],narrow_last['hessian'])),
                includes_refused_numerical_approximations=True,independent_accuracy_not_verified=True)
    require(sha(POLICY_COLLECTOR)==POLICY_COLLECTOR_SHA,'Reconciliation source changed')
    result.update(wide_components=True,wide_width_multiplier=4.,
        preceding_narrow_comparisons_descriptive_only=True,all_refusals_remain_preserved=True,
        numerical_policy_selection_not_automatic=True,independent_derivative_accuracy_verified=False,
        full_model_fit_verified=False,statistical_calibration_verified=False)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    parser.add_argument('--expected-collector-sha')
    args=parser.parse_args();helper=load_policy_collector();entries=metadata(helper)
    if args.inspect:
        print(json.dumps(dict(metadata_verified=True,jobs=JOBS,point_policies=KEYS,required_points=4,
            wide_components=True,wide_width_multiplier=4.,stdlib_only=True,result_reads=False,output_written=False)))
        return
    require(args.expected_collector_sha==sha(Path(__file__)),'Reviewed wide collector source differs')
    result=collect(helper,entries)
    if not result['ready']:print(json.dumps(result));return
    body=json.dumps(result,allow_nan=False,indent=2)+'\n'
    require(len(body.encode())<=32*1024**2,'Wide collection output cap exceeded')
    target=LAUNCH/'mixture_wide_collection.json'
    with target.open('x') as handle:handle.write(body)
    print(json.dumps(dict(collection=str(target),completed_points=4,automatic_policy_selection=False)))


if __name__=='__main__':main()
