#!/usr/bin/env python3
"""Collect three pinned synthetic Decimal diagnostics, never accept accuracy."""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
LAUNCH=REC/'repair_p3_decimal_positive_launch_20261005T161736Z'
CASES=((1,22042086,22041792,22),(2,22042087,22041793,24),(3,22042088,22041794,24))
CHECKER_SHA='afa3ef54af7f66ea3c651c85bbc1234b91a78fb8593c3f271573c7a2d15e559f'
HELPER_SHA='6ff2857ce42e408478f9c5dc7cc13ed6c70bcfdd77bd8f2169fff8de9656ac40'
COUNT_SHA='7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad'
FIELDS=('logmass','variance_score','Hvv')
CAP=32*1024**2
PINS = {'repair_p3_decimal_positive_gene1_22042086/point_110_20.json': 'a709d51e8cd9ba8e2eb1cf5f84fcfa48e09964e0c09b1ec5ad1510d6c55e8df3',
 'repair_p3_decimal_positive_gene1_22042086/point_110_40.json': 'e376a6da156ea44fe33aa7f2f04c84b9d8519938cfc74f78b40543e134864ce3',
 'repair_p3_decimal_positive_gene1_22042086/point_110_80.json': 'dfc4035879e9a568ce6253b4b337dc48db14d59662ebb8df65f3712050c69cab',
 'repair_p3_decimal_positive_gene1_22042086/point_80_20.json': '639840c4eaca2f40e817727a7d7e4502eb8b1b14a3854b247feb6f4038c46d0f',
 'repair_p3_decimal_positive_gene1_22042086/point_80_40.json': '23852f4cacd5c026fa3b188d88253c4091a4e586f7ba943788cc7d1a5d7dd0c0',
 'repair_p3_decimal_positive_gene1_22042086/point_80_80.json': 'a49682f7b52eda8922444e172d45d5d30baac88edcc9bed1da66982bc0ccfab9',
 'repair_p3_decimal_positive_gene1_22042086/protocol.json': 'be83fd7434c4f2c97e27c2d1fe1f87982e423a9ed2e7819a79edc593608700e2',
 'repair_p3_decimal_positive_gene1_22042086/refinements.json': '6bcf256cb325034aed31c4edd7bcf149b9dcb630157d171ebaf9ecc71462e25f',
 'repair_p3_decimal_positive_gene1_22042086/rule_110_20.json': '248e32c247eeef3e9fee15e1b1e61c569622c1fe50b466426754e581485735a4',
 'repair_p3_decimal_positive_gene1_22042086/rule_110_40.json': 'bce1fb6e4d6575219b542a9de45ab66931375bb6a7c7b492c40c06b06af04730',
 'repair_p3_decimal_positive_gene1_22042086/rule_110_80.json': '490a3c2ae701f60582735a1732b8c33a80cbe45152f12e96262b2a34c907ac2b',
 'repair_p3_decimal_positive_gene1_22042086/rule_80_20.json': '320f360b1229ec9d78c957e030c3784b65588c49b0d6a2e025ecb16f79efa9a9',
 'repair_p3_decimal_positive_gene1_22042086/rule_80_40.json': '5ede42e13441d1134d71ac775088f30720edb2c8929b25a6bc76d8874fa9e2a0',
 'repair_p3_decimal_positive_gene1_22042086/rule_80_80.json': '26fd89db0e9a223186b1db3ffbb77a9658cce54809ea1512aaad1d5e2ce2ffbf',
 'repair_p3_decimal_positive_gene1_22042086/summary.json': '29c84c74df7a09ab9a65b099f487eba888989986e68a97bb58ce45380c2844a9',
 'repair_p3_decimal_positive_gene2_22042087/point_110_20.json': '8f25ac8e736bc0245b695952cccbc78074d2ca2225996e2e86b8ee1e49937253',
 'repair_p3_decimal_positive_gene2_22042087/point_110_40.json': '5d8b7ce47199d5eb2c119cbf695adbd7661db956639b00be79e2c57311014d62',
 'repair_p3_decimal_positive_gene2_22042087/point_110_80.json': 'd1820d33094e3d5d9350f84f096b14c8c16d528cab37660a6018e697f85e3f8a',
 'repair_p3_decimal_positive_gene2_22042087/point_80_20.json': '3fc07bce1d64d2026588375481f43ce551dcd4480a6db08193b887708a3799ef',
 'repair_p3_decimal_positive_gene2_22042087/point_80_40.json': '3f32593856e2cec2c96bdc7a4225db58fb893df63eda2ea9467fa7fe63afc5e4',
 'repair_p3_decimal_positive_gene2_22042087/point_80_80.json': 'cb535db0c1b007a9a45dd2763f58a9e4ed4373d45ff36090877e6f05bcae27a1',
 'repair_p3_decimal_positive_gene2_22042087/protocol.json': '365907071e59bcae1b2780df16b7fab66dc68a2868f8c24589007d593fc9594d',
 'repair_p3_decimal_positive_gene2_22042087/refinements.json': 'eaf05e52fd1c41c4be88c04eb81f33c01700cde71970894e0f514f034ee6b7ae',
 'repair_p3_decimal_positive_gene2_22042087/rule_110_20.json': '248e32c247eeef3e9fee15e1b1e61c569622c1fe50b466426754e581485735a4',
 'repair_p3_decimal_positive_gene2_22042087/rule_110_40.json': 'bce1fb6e4d6575219b542a9de45ab66931375bb6a7c7b492c40c06b06af04730',
 'repair_p3_decimal_positive_gene2_22042087/rule_110_80.json': '490a3c2ae701f60582735a1732b8c33a80cbe45152f12e96262b2a34c907ac2b',
 'repair_p3_decimal_positive_gene2_22042087/rule_80_20.json': '320f360b1229ec9d78c957e030c3784b65588c49b0d6a2e025ecb16f79efa9a9',
 'repair_p3_decimal_positive_gene2_22042087/rule_80_40.json': '5ede42e13441d1134d71ac775088f30720edb2c8929b25a6bc76d8874fa9e2a0',
 'repair_p3_decimal_positive_gene2_22042087/rule_80_80.json': '26fd89db0e9a223186b1db3ffbb77a9658cce54809ea1512aaad1d5e2ce2ffbf',
 'repair_p3_decimal_positive_gene2_22042087/summary.json': '9ce0a77cb0c48e79a8fb08301d2b4fbea93a41ca9709453a2b1a92aad8ccc6a8',
 'repair_p3_decimal_positive_gene3_22042088/point_110_20.json': '6c0fd56dfd3e9340613a9901c37aaf58335ad5e38e2abb91248d56603a6da34b',
 'repair_p3_decimal_positive_gene3_22042088/point_110_40.json': '87e7f7b2f733a2eab00b3d61688ba0eae8511f4c2340281532623c9026485623',
 'repair_p3_decimal_positive_gene3_22042088/point_110_80.json': 'abddb50ab291ad83f8024f45b0109c7dd7a7a5d125d6c09e00aa7ca6271bba77',
 'repair_p3_decimal_positive_gene3_22042088/point_80_20.json': 'd1f1f935679e4667d0ccba339924e0daf2badb01f60b8e347487c2ea4be9a5dd',
 'repair_p3_decimal_positive_gene3_22042088/point_80_40.json': 'd2ceec066789dc17b9135129403659f793e8121ae0f7e24a4a118046dfddd647',
 'repair_p3_decimal_positive_gene3_22042088/point_80_80.json': 'cd537c6de0f9dd7217a0946a30f4b68e16bdddd0f7b82ee6e814af92823889b6',
 'repair_p3_decimal_positive_gene3_22042088/protocol.json': '3dd3d6594a6fb0f9eb38e95560a558dcd8238d1f89000de870120ea5c1bf0fc4',
 'repair_p3_decimal_positive_gene3_22042088/refinements.json': 'bb765056cbedf082456ea7361acbe67748d06539dd8634c3d35391169ccab0b3',
 'repair_p3_decimal_positive_gene3_22042088/rule_110_20.json': '248e32c247eeef3e9fee15e1b1e61c569622c1fe50b466426754e581485735a4',
 'repair_p3_decimal_positive_gene3_22042088/rule_110_40.json': 'bce1fb6e4d6575219b542a9de45ab66931375bb6a7c7b492c40c06b06af04730',
 'repair_p3_decimal_positive_gene3_22042088/rule_110_80.json': '490a3c2ae701f60582735a1732b8c33a80cbe45152f12e96262b2a34c907ac2b',
 'repair_p3_decimal_positive_gene3_22042088/rule_80_20.json': '320f360b1229ec9d78c957e030c3784b65588c49b0d6a2e025ecb16f79efa9a9',
 'repair_p3_decimal_positive_gene3_22042088/rule_80_40.json': '5ede42e13441d1134d71ac775088f30720edb2c8929b25a6bc76d8874fa9e2a0',
 'repair_p3_decimal_positive_gene3_22042088/rule_80_80.json': '26fd89db0e9a223186b1db3ffbb77a9658cce54809ea1512aaad1d5e2ce2ffbf',
 'repair_p3_decimal_positive_gene3_22042088/summary.json': 'c6fc4f6e5cb2588ffaeff3b637aae48462dde700c8609948e04a4f5c5367b6a6',
 'repair_p3_decimal_positive_launch_20261005T161736Z/source_hashes.json': 'f82b4498e1a81e4dcc8026bbf7a2d2709e1e997f26e3d1520642376dceb22a07',
 'repair_p3_decimal_positive_launch_20261005T161736Z/submissions.jsonl': '631acfc8a0d9193a4212030f1ea739c076ff9b58d99026b0c02e3c172768c1a3',
 'repair_p3_theta4_order_gene1_22041792/order_diagnostic.json': '3703304c193d1a17d0c13ecaf6521db48c62f80dba20458f60533c4a278f5269',
 'repair_p3_theta4_order_gene1_22041792/original_refusal.json': '63581cfe0ffd63d4e6f5445705fd50b99cd99ba9e3669dde54bbd039b403d2f2',
 'repair_p3_theta4_order_gene1_22041792/protocol.json': '922d8065cfc0e32b6b34043d32a1ce41004c61db66a14bc3a1f4de94b74f5b12',
 'repair_p3_theta4_order_gene1_22041792/summary.json': 'd601185914720b949726a474ce5c3d14040564e477eda65baa2490f02b186e23',
 'repair_p3_theta4_order_gene2_22041793/order_diagnostic.json': '71fa9399e8e3487cb955b6b743d52e59cd018a6e88396d1654b205b44c0b9a28',
 'repair_p3_theta4_order_gene2_22041793/original_refusal.json': 'f4a4e82a99d719ecf695e6fceb67ce3ea14aeb3c4f895c73a4d76da6be7842dc',
 'repair_p3_theta4_order_gene2_22041793/protocol.json': '8ad8beaf12666ffe13637769175a5f8ff5b4938fb8f549fd2351601a491500f8',
 'repair_p3_theta4_order_gene2_22041793/summary.json': '4d2066888e4b7033b0d9fc09db72c5684ebc05513483046d6c2a2db4e1a5dc0f',
 'repair_p3_theta4_order_gene3_22041794/order_diagnostic.json': 'aa757da89c9c72efd3bc19128f57999840011dee2f5a743551181b115151f141',
 'repair_p3_theta4_order_gene3_22041794/original_refusal.json': 'd3d2fb0579ef0c6f4bb95c03daf4e1fc9fd15c2ef7d5ad2989475485200b4fac',
 'repair_p3_theta4_order_gene3_22041794/protocol.json': 'e9bc7c74e3693dfab813baa955b250d220c9557778923387f8d51feae6d985f3',
 'repair_p3_theta4_order_gene3_22041794/summary.json': '34079b11b4bebe1fb3742c13816b97fa86e04115d4bb1a535465a24bbf6dfeb6'}


def require(ok,message):
    if not ok: raise ValueError(message)


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    require(path.stat().st_size<=CAP,'Input exceeds32MiB: '+str(path))
    require(sha(path)==PINS[str(path.relative_to(REC))],'Completed receipt changed: '+str(path))
    return json.loads(path.read_text())


def protocol_pins(protocol,directory):
    for name,digest in protocol['source_and_input_sha256'].items():
        path=Path(name)
        require(path.is_relative_to(ROOT),'Pinned input outside project synthetic tree')
        require(sha(path)==digest,'Recorded source/input differs: '+name)
        if path.suffix=='.py':
            require(sha(directory/'executed_sources'/path.name)==digest,'Executed archive differs')


def validate():
    for name,digest in PINS.items(): require(sha(REC/name)==digest,'Receipt/launch pin differs: '+name)
    sources=read(LAUNCH/'source_hashes.json')
    require(sources['repair_p3_decimal_positive_check.py']==CHECKER_SHA
            and sources['repair_p3_decimal_positive.py']==HELPER_SHA,'Frozen Decimal sources differ')
    for name,digest in sources.items(): require(sha(LAUNCH/name)==digest,'Launch source differs')
    submitted=[json.loads(line) for line in (LAUNCH/'submissions.jsonl').read_text().splitlines()]
    require([(x['gene'],int(x['job_id'])) for x in submitted]==[(g,j) for g,j,_,_ in CASES],
            'Exact three submitted gene/jobs differ')
    packets=[]; family_p=None
    for gene,job,gljob,rows in CASES:
        directory=REC/('repair_p3_decimal_positive_gene%d_%d'%(gene,job))
        gl_directory=REC/('repair_p3_theta4_order_gene%d_%d'%(gene,gljob))
        protocol,summary=read(directory/'protocol.json'),read(directory/'summary.json')
        protocol_pins(protocol,directory)
        require(sha(directory/'executed_sources/repair_p3_decimal_positive_check.py')==CHECKER_SHA
                and sha(directory/'executed_sources/repair_p3_decimal_positive.py')==HELPER_SHA,
                'Executed Decimal implementation differs')
        require(protocol['gene']==gene and protocol['theta']=='theta4_perturbed'
            and protocol['fixture_seed']==20260930 and protocol['fixed_v']==.0004
            and protocol['fixed_gene_rows']==rows and protocol['orders']==[20,40,80]
            and protocol['decimal_digits']==[80,110]
            and protocol['full_fixture_count_vector_sha256']==COUNT_SHA
            and protocol['independent_Hvv_accuracy_accepted'] is False,'Decimal identity/scope differs')
        p=protocol['fixed_p']
        parameter_path=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z/params/theta4_perturbed.json'
        parameter=json.loads(parameter_path.read_text())
        expected_p=parameter['alpha']+parameter['r']+[parameter['kappa']]+parameter['delta']+[
            parameter[k] for k in ('omega_S','rho_S','rho_b','rho_d')]
        require(p==expected_p and parameter['v']==.0004,'Exact fixture parameter vector differs')
        require(len(p)==28 and all(isinstance(x,(int,float)) and math.isfinite(x) for x in p),
                'Finite28-coordinate parameter vector required')
        if family_p is None: family_p=p
        require(p==family_p,'Parameters differ across Decimal genes')
        require(summary['completed'] is True and summary['gene']==gene
                and summary['case']=='theta4_perturbed' and summary['probability_evaluations']==280
                and summary['independent_Hvv_accuracy_accepted'] is False
                and summary['full_Hessian_accuracy_verified'] is False,'Summary completion/limits differ')
        records=summary['records']
        require([(r['decimal_digits'],r['order']) for r in records]
                ==[(d,o) for d in (80,110) for o in (20,40,80)],'Six-record family differs')
        for record in records:
            d,o=record['decimal_digits'],record['order']
            require(read(directory/('point_%d_%d.json'%(d,o)))==record,'Saved point differs from summary')
            rule=read(directory/('rule_%d_%d.json'%(d,o)))
            require(rule['order']==o and rule['decimal_digits']==d
                    and len(rule['roots'])==o and len(rule['weights'])==o
                    and record['independent_Hvv_accuracy_accepted'] is False,'Rule/point metadata differs')
            require(all(isinstance(record[field],str) and Decimal(record[field]).is_finite()
                        for field in (*FIELDS,'raw_variance_second_moment')),'Nonfinite Decimal point')
        saved_refinements=read(directory/'refinements.json')
        require(saved_refinements==summary['refinements'] and len(saved_refinements)==7,
                'Saved seven refinements differ')
        gl_protocol=read(gl_directory/'protocol.json'); protocol_pins(gl_protocol,gl_directory)
        for name,digest in protocol['source_and_input_sha256'].items():
            if '/model-a-bfix-20260930T115521Z/' in name:
                require(gl_protocol['source_and_input_sha256'].get(name)==digest,'GL/Decimal fixture pin differs')
        gl=read(gl_directory/'order_diagnostic.json')
        refusal=read(gl_directory/'original_refusal.json')
        gl_summary=read(gl_directory/'summary.json')
        require(gl_protocol['gene']==gene and gl_protocol['theta']=='theta4_perturbed'
                and gl_protocol['fixture_seed']==20260930 and gl_protocol['orders']==[16,32,64,128]
                and gl['p']==p and gl['v']==.0004 and gl['rows']==rows
                and refusal['p']==p and refusal['v']==.0004 and refusal['rows']==rows
                and 'candidate_refusal' in refusal and gl_summary['original_refusal_preserved'] is True,
                'GL point/refusal identity differs')
        require([r['order'] for r in gl['report']['ladder']]==[16,32,64,128],
                'Every GL order required')
        require(gl['report']['no_accuracy_acceptance'] is True
                and gl_summary['no_accuracy_acceptance'] is True,'GL acceptance disposition differs')
        for r in gl['report']['ladder']:
            require(len(r['gradient'])==29 and len(r['hessian'])==29
                    and all(len(row)==29 for row in r['hessian'])
                    and all(math.isfinite(x) for x in (r['logmass'],r['gradient'][-1],r['hessian'][-1][-1])),
                    'Finite GL LL/gv/Hvv required')
        packets.append(dict(gene=gene,job=job,gl_job=gljob,protocol=protocol,summary=summary,
            gl=gl,original_refusal=refusal,saved_refinements=saved_refinements))
    return packets


def collect(packets):
    result=[]
    with localcontext() as context:
        context.prec=110
        for packet in packets:
            records=packet['summary']['records']; refinements=[]
            for d in (80,110):
                series=[r for r in records if r['decimal_digits']==d]
                for a,b in zip(series[:-1],series[1:]):
                    refinements.append(dict(kind='order',digits=d,first_order=a['order'],second_order=b['order'],
                        differences={k:str(Decimal(b[k])-Decimal(a[k])) for k in FIELDS}))
            for order in (20,40,80):
                a,b=[r for r in records if r['order']==order]
                refinements.append(dict(kind='digits',order=order,first_digits=80,second_digits=110,
                    differences={k:str(Decimal(b[k])-Decimal(a[k])) for k in FIELDS}))
            require(refinements==packet['saved_refinements'],'110-digit recomputed refinements differ')
            final=[r for r in records if r['decimal_digits']==110 and r['order']==80][0]
            comparisons=[]
            for gl in packet['gl']['report']['ladder']:
                values={'logmass':gl['logmass'],'variance_score':gl['gradient'][-1],'Hvv':gl['hessian'][-1][-1]}
                comparisons.append(dict(order=gl['order'],GL_values=values,
                    GL_minus_Decimal_110_order80={k:str(Decimal.from_float(float(values[k]))-Decimal(final[k])) for k in FIELDS},
                    comparison_descriptive_only=True))
            result.append(dict(gene=packet['gene'],job=packet['job'],gl_job=packet['gl_job'],
                Decimal_records=records,recomputed_refinements=refinements,
                GL_order_comparisons=comparisons,original_refusal=packet['original_refusal'],
                independent_Hvv_accuracy_accepted=False,full_Hessian_accuracy_verified=False))
    return dict(cases=result,input_sha256=PINS,required_cases=3,completed_cases=3,
        independent_Hvv_accuracy_accepted=False,full_Hessian_accuracy_verified=False,
        no_accuracy_acceptance=True,no_model_fit=True,no_calibration=True,
        uncertainty='Arithmetic/order differences and tail diagnostics do not certify quadrature error.',
        comparison_float_conversion='Exact binary GL values via Decimal.from_float; Decimal reference110digits/order80')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-collector-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args()
    require(sha(Path(__file__).resolve())==args.expected_collector_sha,'Collector source differs')
    packets=validate()
    if args.inspect:
        print(json.dumps(dict(completed_cases=len(packets),input_hashes_checked=True,
            saved_points_checked=True,result_written=False,independent_Hvv_accuracy_accepted=False)))
        return
    body=(json.dumps(collect(packets),indent=2,allow_nan=False)+'\n').encode()
    require(len(body)<=CAP,'Collection exceeds32MiB')
    with (LAUNCH/'decimal_positive_collection.json').open('xb') as handle: handle.write(body)


if __name__=='__main__': main()
