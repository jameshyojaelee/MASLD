#!/usr/bin/env python3
"""Isolated prepared callback check of the immutable composite GL fixed points.

No fits/calibration. Every archived GL ladder and precision refusal is retained.
Only callback preparation, labels and guarded saved-output comparisons are added.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
SNAP = REC/'repair_p3_composite_gl_launch_20261005T132600Z'
SOURCE = SNAP/'repair_p3_fixed_point_check.py'
SOURCE_SHA = '3a7e42b13529f41b0f2bd8eeab17c8f8c5bc66410ba2ccdc734274a59ec29a2c'
KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
GL_SHA = '45b4423811fded58712bd5e85666f9ffc2b6972abbd9bd71964a584ec6d3c59c'
PREPARED = Path(__file__).with_name('repair_p3_prepared_conditional_v2.py')
PREPARED_SHA = '289e93b651411042093a2b987717d87cdb5ed156ab7c736e1ac5cf820670fd4b'
SAVED_GUARDS = {
    "repair_p3_fixed_point_failed0gene2_22039200/protocol.json": "7b5d793be5485a25e5fc7f5b3d97c4ddf30a1a497d697728cd64468ed61f8d40",
    "repair_p3_fixed_point_failed0gene2_22039200/summary.json": "5ad1ea571c71525e6ae1eb18bf625baf1e6779934d64a149e0d16ca5f3137a0b",
    "repair_p3_fixed_point_failed0gene2_22039200/candidate_baseline.json": "edec063451383f9a1386619186c2089d3ef8e6ee34d854e80a2894fdd28e0f9a",
    "repair_p3_fixed_point_positive18_22039201/protocol.json": "ed6c190424fa29501e1377b5b47b12762917e53f46338c3c3f1a408c3f9e4774",
    "repair_p3_fixed_point_positive18_22039201/summary.json": "1daa3275c8d524857c90369b7686fc2234229df765e41f139a4b232c9e323613",
    "repair_p3_fixed_point_positive18_22039201/candidate_baseline.json": "4a576370be2f33ac28bf23052fe283b4a5f614bc30e983bd6a84224da857c722",
    "repair_p3_complex_failed0gene2_22039439/protocol.json": "fe5ea45c4f8c65b474c35b3f1e88d2b7c9036b6e28ab85d1d67abe072c9c49c3",
    "repair_p3_complex_failed0gene2_22039439/summary.json": "af9320034f5bffb24ecf5dedd63094aee676271a5a460b9e3b6318e8138bde7e",
    "repair_p3_complex_failed0gene2_22039439/all_derivatives.json": "d9aeb40099bccde6b42f782e13732c6a1679ca396faab862b3e64a474491948d",
    "repair_p3_complex_failed0gene2_22039440/protocol.json": "fc2858b480761a8160641428e34ef0f10fda5bfe4236dc1c26d3d3b26aee807f",
    "repair_p3_complex_failed0gene2_22039440/summary.json": "0e8ae4b6bda9cb2b848620bdb01760038d4966b5e3aeab36581cb0096e05909c",
    "repair_p3_complex_failed0gene2_22039440/all_derivatives.json": "ab29bc063460e02add39be2c1b1548da62a96cd6180fe5a90c7cabb39631e038",
    "repair_p3_complex_failed0gene2_22039441/protocol.json": "5042f8c5bbae08fb97960b89b3b3ee4a9e984637a25038b7bdd92475f8241eb9",
    "repair_p3_complex_failed0gene2_22039441/summary.json": "438a1d84cda700e2d9ae79d7c75091eb520191a97503d7d3588d7564e97c3b2d",
    "repair_p3_complex_failed0gene2_22039441/all_derivatives.json": "6b6c92abc1cf55be423ace4f8c088b2d66590ce7210466e0ac2677c23455279d",
    "repair_p3_complex_failed0gene2_22039442/protocol.json": "f7e38051cd734af7247ae4d977b8a757ac8afd667f5fe6156e4572f9d64d250e",
    "repair_p3_complex_failed0gene2_22039442/summary.json": "274d42319a9ecf16ea3bb0dd1efc14beee71c5f41e6ace3fb220ec2824c03b67",
    "repair_p3_complex_failed0gene2_22039442/all_derivatives.json": "e7aa8f81cd85a00792affa68fdb6ae07dc66ac4f46aceebc7d793f05db8a49f6",
    "repair_p3_complex_failed0gene2_22039443/protocol.json": "7bdae7a7198aa7516ec351e0213f1c0241ce48b719e59cab1cb188898df44ce3",
    "repair_p3_complex_failed0gene2_22039443/summary.json": "efd376e3f2ffaa983cb42e550d9c4b3c5ecc48f03f8fae7a474f8d40066869d7",
    "repair_p3_complex_failed0gene2_22039443/all_derivatives.json": "df55af4f9a1126e5cfb5e61a87ca0a69e7d6074980dd7d83c7490f2ec45a499c",
    "repair_p3_complex_failed0gene2_22039444/protocol.json": "71dbb90cffbad32b30344bd096e2aa52d943b8a60f4dc2f7c48558adae16cd84",
    "repair_p3_complex_failed0gene2_22039444/summary.json": "f25213d2ed7c5ad88065652882c5c174a82f52d4cab670d0c5c190c4a59b02ab",
    "repair_p3_complex_failed0gene2_22039444/all_derivatives.json": "0a0338f1af8f539aa6cb92ac6b28a1e6154a391d19fc4263c858ccc3bd75a48a",
    "repair_p3_complex_failed0gene2_22039445/protocol.json": "00fbed855181a0cea73dab11b1499377d6469dd1e4f69c6b77a1436a8bf090b7",
    "repair_p3_complex_failed0gene2_22039445/summary.json": "aedc3cf2da9ce8079916ed4305528a1539e0ea6dc8e0d2e9625231df9b622691",
    "repair_p3_complex_failed0gene2_22039445/all_derivatives.json": "c2c9572a00df1ade22afc5736c7f815063cf56394f5f695b733153488fa5d3ee",
    "repair_p3_complex_positive18_22039446/protocol.json": "524d8244d5fbfc9e7c7c3d65257507d2d3f8ffd086cb8abbf2815bd55c537358",
    "repair_p3_complex_positive18_22039446/summary.json": "97e641a8e48446510b9eb58bc5cb521d4d8e3092787ce2ffb3f24e238cae8f46",
    "repair_p3_complex_positive18_22039446/all_derivatives.json": "62da00da72361f0f87c572b8d24ced77732e2432702c2e7d8515311ce9ba73eb",
    "repair_p3_complex_positive18_22039447/protocol.json": "351fa2f255ae718c647ad342986e9c9e6151b5ce7e94b3f239aadb4eb6750351",
    "repair_p3_complex_positive18_22039447/summary.json": "1d802158ad729507753c8681c4b33d7580401ff190514b4294fe46271852153f",
    "repair_p3_complex_positive18_22039447/all_derivatives.json": "ae876570352f438e7313300f90aa610fce60e9312d5047fd93a37d49380d58a9",
    "repair_p3_complex_positive18_22039448/protocol.json": "d04803d6f2956c4e55725b449e7dcfacf0722c984c694bf7d3da0637299e2cc2",
    "repair_p3_complex_positive18_22039448/summary.json": "6b741e13dcea07c9d1160ee73e84358b58f7d85a36e6e2cfb7ce91ffa112785b",
    "repair_p3_complex_positive18_22039448/all_derivatives.json": "5067278e8514501c395ae8a9b6156e201cb20956b4bba328eec5cc1dd291ebda",
    "repair_p3_complex_positive18_22039449/protocol.json": "bfba376f92b16e9a05cd1ab7e9624a53c7ebcd057cf964a602e8ca1f97466650",
    "repair_p3_complex_positive18_22039449/summary.json": "9ac4b23d073bac004c9905a7848f52859e1fa7ca7fef702626e87eb147c79a94",
    "repair_p3_complex_positive18_22039449/all_derivatives.json": "ece1cc3ca482dd681aa50ceab826cf2490a7b5567642208f2cf061731209287d",
    "repair_p3_complex_positive18_22039450/protocol.json": "de31ef38e0c8b21f8372a0440e82a2916f7616d34d2b9c222cfdfefd967a6b3c",
    "repair_p3_complex_positive18_22039450/summary.json": "137cf21d0f306b3639f52eebf6771024c3e1fc6a89dddfd8e98a86644ad0ceb2",
    "repair_p3_complex_positive18_22039450/all_derivatives.json": "e7d990eee4fee2aa690dc4c5a249b4f6342d66e4cf4d30f89386933df724e442",
    "repair_p3_complex_positive18_22039451/protocol.json": "5833a43721d52f28b5b14256e85c5f9375e711ee1952f74403729302ac79b946",
    "repair_p3_complex_positive18_22039451/summary.json": "f3a1926338467b9e1f02bad66ac0ab2ec3805f539eb55a2a3d7fe80606f3aa93",
    "repair_p3_complex_positive18_22039451/all_derivatives.json": "bd65ecbe28e25455d1852894fc1a63d92665c3d281dd27b0d3f86d8e9dca6750",
    "repair_p3_complex_positive18_22039452/protocol.json": "9dd98219a28bf31cbcb716eccd35b5ad3af8be9c58c19d15f9bedb0e2992a41f",
    "repair_p3_complex_positive18_22039452/summary.json": "c98adbb124dcaf68b25544982b64b3e19eb39d8a37a091ad0ee1fed1c538ff81",
    "repair_p3_complex_positive18_22039452/all_derivatives.json": "e191b4fd99090670d4fba2ae91ae0d85233110971e704acbccf5dfe5fabc2b41"
}
BASES = {'failed0gene2':'repair_p3_fixed_point_failed0gene2_22039200',
         'positive18':'repair_p3_fixed_point_positive18_22039201'}
COMPLEX = {'failed0gene2':range(22039439,22039446),
           'positive18':range(22039446,22039453)}
AXES = (2,10,16,17,24,25,28)
CAP = 32*1024**2


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=tuple(BASES), required=True)
    parser.add_argument('--expected-checker-sha', required=True)
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()
    require(sha(Path(__file__)) == args.expected_checker_sha, 'Own checker source changed')
    require(sha(SOURCE) == SOURCE_SHA and sha(PREPARED) == PREPARED_SHA, 'Pinned source changed')
    for relative,digest in SAVED_GUARDS.items():
        require(sha(REC/relative) == digest, 'Guarded completed synthetic output changed: '+relative)
    source = SOURCE.read_text()
    old = """            def integrate_gene(*call_args, **call_kwargs):
                return composite_gl.integrate_gene(*call_args, **call_kwargs,
                                                  expected_kernel_sha=args.expected_kernel_sha,
                                                  kernel_path=archive / KERNEL.name)
"""
    new = """            prepared_module = load_module('prepared_GL_callback_module', archive / 'repair_p3_prepared_conditional_v2.py')
            def integrate_gene(*call_args, **call_kwargs):
                import time as callback_time
                require(len(call_args) == 5, 'Unexpected composite callback argument signature')
                fixed_model, fixed_p, fixed_v, fixed_rows, original_conditional = call_args
                setup_start = callback_time.perf_counter()
                prepared_context = prepared_module.PreparedConditional(fixed_model, fixed_p, fixed_rows)
                setup_seconds = callback_time.perf_counter()-setup_start
                require(prepared_context.original.conditional_terms is original_conditional, 'Preparation conditional source differs')
                def prepared_callback(row, eta, linear_r, need_hessian=True):
                    return prepared_context.terms(row, eta, linear_r, fixed_p, need_hessian)
                integrate_start = callback_time.perf_counter()
                result = composite_gl.integrate_gene(fixed_model, fixed_p, fixed_v, fixed_rows,
                    prepared_callback, **call_kwargs,
                    expected_kernel_sha=args.expected_kernel_sha, kernel_path=archive / KERNEL.name)
                result[3]['prepared_callback'] = dict(setup_seconds=setup_seconds,
                    integration_seconds=callback_time.perf_counter()-integrate_start,
                    context_scope='all ladders at one fixed p only', original_model_globals_unchanged=True)
                require(prepared_context.original.conditional_terms is original_conditional, 'Original conditional changed')
                return result
"""
    require(source.count(old) == 1, 'Expected original explicit callback wrapper differs')
    source = source.replace(old,new)
    label = 'repair_p3_fixed_point_'
    require(source.count(label) == 1, 'Expected original output prefix differs')
    source = source.replace(label,'repair_p3_prepared_gl_')
    namespace = {'__file__':str(SOURCE),'__name__':'isolated_prepared_GL_reconstruction'}
    exec(compile(ast.parse(source),str(SOURCE),'exec'),namespace)
    namespace['GUARDS'][PREPARED] = PREPARED_SHA
    namespace['GUARDS'][Path(__file__).resolve()] = args.expected_checker_sha
    # Completed outputs are small synthetic result JSONs; no private molecular inputs.
    namespace['GUARDS'].update({REC/p:digest for p,digest in SAVED_GUARDS.items()})
    sys.argv = [str(SOURCE),'--case',args.case,'--baseline-only','--fast-reference',
                '--composite-gl-sha',GL_SHA,'--expected-helper-sha',SOURCE_SHA,
                '--expected-kernel-sha',KERNEL_SHA]
    if args.inspect:
        sys.argv.append('--inspect')
    elif not os.environ.get('SLURM_JOB_ID','').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID'):
        raise ValueError('Individual compute allocation required')
    namespace['main']()
    if args.inspect:
        return
    output = REC/('repair_p3_prepared_gl_'+args.case+'_'+os.environ['SLURM_JOB_ID'])
    def save(name, record):
        body = json.dumps(record,indent=2,allow_nan=False)+'\n'
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body.encode()) <= CAP-65536,
                'Prepared GL comparison output cap exhausted')
        with (output/(name+'.json')).open('x') as handle:
            handle.write(body)
            handle.flush()
            os.fsync(handle.fileno())
    phase = 'saved_candidate_comparison'
    try:
        import numpy as np
        completed = json.loads((output/'summary.json').read_text())
        require(completed['candidate_integral_returned'], 'Prepared GL candidate refused; original refusal preserved')
        actual = json.loads((output/'candidate_baseline.json').read_text())
        protocol = json.loads((output/'protocol.json').read_text())
        old_directory = REC/BASES[args.case]
        prior = json.loads((old_directory/'candidate_baseline.json').read_text())
        prior_protocol = json.loads((old_directory/'protocol.json').read_text())
        for key in ('case','draw','fixed_gene','fixed_gene_rows','fixed_p','fixed_v','count_vector_sha256','environment','candidate_settings'):
            require(protocol[key] == prior_protocol[key], 'Saved fixed-point identity differs: '+key)
        require(prior['integration']['refinement_differences_not_certified_error_bounds'] and
                actual['integration']['refinement_differences_not_certified_error_bounds'], 'Expected composite GL records missing')
        require(len(prior['integration']['ladder']) == len(actual['integration']['ladder']) == 9,
                'Original ladder count differs')
        for left,right in zip(prior['integration']['ladder'],actual['integration']['ladder']):
            for key in ('z_bound','grid_spacing','panels','order','node_evaluations'):
                require(left[key] == right[key], 'Original GL ladder differs: '+key)
        require(prior['integration']['likelihood_evaluations'] == actual['integration']['likelihood_evaluations'],
                'Original GL evaluation count differs')
        deltas = dict(logmass=abs(actual['logmass']-prior['logmass']),
            gradient=float(np.max(np.abs(np.asarray(actual['gradient'])-prior['gradient']))),
            hessian=float(np.max(np.abs(np.asarray(actual['hessian'])-prior['hessian']))))
        require(deltas['logmass'] <= 1e-9 and deltas['gradient'] <= 1e-8 and deltas['hessian'] <= 1e-5,
                'Saved composite comparison exceeds unchanged GL budgets')
        phase = 'existing_complex_reference_comparison'
        comparisons = []
        seen_axes = []
        for job in COMPLEX[args.case]:
            matches = [p for p in SAVED_GUARDS if p.endswith('/all_derivatives.json') and
                       p.split('/')[0].endswith('_'+str(job))]
            require(len(matches) == 1, 'Exact complex reference family missing')
            d = (REC/matches[0]).parent
            cp = json.loads((d/'protocol.json').read_text())
            cs = json.loads((d/'summary.json').read_text())
            for key in ('case','draw','fixed_gene','fixed_gene_rows','fixed_p','fixed_v','count_vector_sha256','environment'):
                require(cp[key] == protocol[key], 'Complex point identity differs: '+key)
            require(cs['saved_record_count'] == 6 and cs['successful_record_count'] == 6 and cs['refusal_count'] == 0,
                    'Existing complex records are incomplete')
            records = json.loads((d/'all_derivatives.json').read_text())
            require(len(records) == 6 and len(cp['requested_axes']) == 1, 'Unexpected existing complex shape')
            seen_axes += cp['requested_axes']
            for record in records:
                require(record['completed'] and record['axis'] == cp['requested_axes'][0], 'Complex record refused/mismatched')
                axis = record['axis']
                reference = record['result']
                difference = actual['gradient'][axis]-reference['gradient']
                comparisons.append(dict(source=str(d/'all_derivatives.json'),axis=axis,
                    step=record['step'],bound=record['bound'],candidate_score=actual['gradient'][axis],
                    complex_score=reference['gradient'],difference=difference,
                    reference_estimated_uncertainty=reference['estimated_combined_gradient_uncertainty'],
                    candidate_integration_allowance=actual['integration']['refined']['gradient_error_allowance'][axis],
                    candidate_refinement_difference=actual['integration']['gradient_refinement_difference'][axis],
                    reference_uncertified=True,complex_arithmetic_and_truncation_not_bounded=True,
                    derivative_accuracy_acceptance_not_inferred=True))
        require(tuple(seen_axes) == AXES and len(comparisons) == 42, 'Complete complex axes/records missing')
        for path,digest in namespace['GUARDS'].items():
            require(sha(path) == digest, 'Guarded source/saved output changed after comparisons')
        save('prepared_gl_saved_comparisons',dict(completed=True,case=args.case,
            old_candidate_differences=deltas,all_LLH_exact_equal=all(d == 0 for d in deltas.values()),
            unchanged_GL_ladder_count=len(actual['integration']['ladder']),
            unchanged_GL_evaluations=actual['integration']['likelihood_evaluations'],
            complex_comparisons=comparisons,source_guards=SAVED_GUARDS,
            no_model_fits=True,no_calibration=True,
            independent_complex_reference_not_newly_computed=True,
            no_adoption_or_frozen20_node_agreement_claim=True))
    except Exception:
        save('prepared_gl_comparison_failure',dict(completed=False,phase=phase,
            traceback=traceback.format_exc(),all_original_and_partial_records_preserved=True))
        raise


if __name__ == '__main__':
    main()
