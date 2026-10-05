"""Hash-guarded reuse of completed synthetic fixed-point baseline calculations."""
import hashlib
import json
import os
from pathlib import Path
import re

REC = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z')
FILES = ('protocol.json','candidate_baseline.json','independent_baseline.json','summary.json')
FIELDS = ('case','draw','fixed_gene','fixed_gene_rows','fixed_p','fixed_v',
          'count_vector_sha256','environment','candidate_settings','unchanged_outer_methods',
          'finite_difference_steps','required_combined_derivative_family_axes',
          'reference_arithmetic','composite_gl','original_records_preserved',
          'approximation_revision_not_frozen20_node_agreement')
PINS = ('repair_p3_integral.py','repair_p3_composite_gl.py','repair_p3_reference_fast.py',
        'codex_stable_bb_probe.py','likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py')
SETTINGS = dict(logmass_tolerance=1e-9,gradient_tolerance=1e-8,hessian_tolerance=1e-5,
                quad_relative=2e-12,quad_absolute=2e-12,initial_z_bound=8.,
                maximum_z_bound=96.,grid_spacing=.25,maximum_evaluations=100000,
                maximum_intervals=6000)


def load_reused_baseline(receipt_dir, expected_sha256, current_protocol):
    if not os.environ.get('SLURM_JOB_ID','').isdigit():
        raise ValueError('Compute allocation required for numerical receipt reuse')
    import numpy as np

    def require(ok,message):
        if not ok:
            raise ValueError(message)

    case = current_protocol.get('case')
    require(case in ('positive18','failed0gene2'),'Unsupported positive synthetic case')
    directory = Path(receipt_dir).resolve(strict=True)
    require(directory.parent == REC.resolve() and
            re.fullmatch('repair_p3_fixed_point_'+case+r'_\d+',directory.name) is not None,
            'Receipt must be direct synthetic case child of REC')
    require(set(expected_sha256) == set(FILES),'Exactly four explicit receipt hashes required')
    documents = {}
    for name in FILES:
        path = directory/name
        require(not path.is_symlink() and path.resolve(strict=True).parent == directory,
                'Receipt file escaped its directory')
        body = path.read_bytes()
        require(hashlib.sha256(body).hexdigest() == expected_sha256[name],'Receipt bytes differ')
        documents[name] = json.loads(body)
    prior = documents['protocol.json']
    summary = documents['summary.json']
    for key in FIELDS:
        require(key in prior and key in current_protocol and prior[key] == current_protocol[key],
                'Current fixed point differs: '+key)
    require(prior['candidate_settings'] == SETTINGS and prior['fixed_v'] > 0 and
            len(prior['fixed_p']) == 28 and prior['fixed_gene'] == 2 and prior['fixed_gene_rows'] == 24,
            'Reviewed positive fixed-point dimensions/settings differ')
    require(prior['composite_gl'] is True and
            prior['reference_arithmetic'] == 'independent cumulative integer rising products' and
            prior.get('high_precision_boundary') is False and prior.get('predictor_points_only') is False and
            prior.get('no_model_fit') is True and prior.get('no_calibration') is True and
            prior.get('baseline_only') is True and all(prior['unchanged_outer_methods'].values()),
            'Parent is not the unchanged composite GL baseline calculation')

    def source_pins(protocol):
        result = {}
        for mapping in ('source_guards','own_source_hashes'):
            for path,digest in protocol[mapping].items():
                name = Path(path).name
                if name in PINS:
                    require(name not in result or result[name] == digest,'Conflicting source basename pin')
                    result[name] = digest
        require(set(result) == set(PINS),'Required mathematical source pins missing')
        return result

    require(source_pins(prior) == source_pins(current_protocol),'Mathematical source pins differ')
    require(summary.get('fixed_point_check_completed') is True and
            summary.get('candidate_integral_returned') is True and
            summary.get('no_model_fits') is True and summary.get('no_statistical_calibration') is True and
            summary.get('baseline_only') is True and summary.get('derivative_checks_completed') is False and
            summary.get('executed_derivative_axes') == [] and summary.get('case') == case and
            summary.get('draw') == prior['draw'] and summary.get('positive_variance_gene_scope') == 2 and
            summary.get('finite_difference_steps') == prior['finite_difference_steps'] and
            summary.get('required_combined_derivative_family_axes') == prior['required_combined_derivative_family_axes'],
            'Parent summary does not establish completed baseline-only computation')
    candidate = documents['candidate_baseline.json']
    def finite_numbers(value):
        if isinstance(value,dict):
            return all(finite_numbers(v) for v in value.values())
        if isinstance(value,list):
            return all(finite_numbers(v) for v in value)
        if isinstance(value,(int,float)):
            return bool(np.isfinite(value))
        return True
    require(finite_numbers(candidate) and finite_numbers(documents['independent_baseline.json']) and
            np.isfinite(np.asarray(prior['fixed_p'],dtype=float)).all(),
            'Nonfinite numerical receipt values')
    gradient,hessian = np.asarray(candidate['gradient'],dtype=float),np.asarray(candidate['hessian'],dtype=float)
    require(gradient.shape == (29,) and hessian.shape == (29,29) and
            np.isfinite(gradient).all() and np.isfinite(hessian).all() and
            np.isfinite(candidate['logmass']),'Nonfinite or wrong-sized baseline')
    integration = candidate['integration']
    require(integration.get('same_probability_estimand') is True and
            integration.get('fixed_point_only') is True and
            integration.get('refinement_differences_not_certified_error_bounds') is True and
            integration.get('method') == 'composite Gauss-Legendre on every finite-grid interval',
            'Integration calculation identity differs')
    independent = documents['independent_baseline.json']
    for key in ('first','tighter_expanded'):
        require(np.isfinite(independent[key]['logmass']) and
                np.isfinite(independent[key]['estimated_quad_plus_mass_tail_log_error']) and
                independent[key]['estimated_quad_plus_mass_tail_log_error'] >= 0,
                'Independent scalar baseline invalid')

    def arrays(value):
        if isinstance(value,dict):
            return {k:arrays(v) for k,v in value.items()}
        if isinstance(value,list):
            if value:
                try:
                    array = np.asarray(value,dtype=float)
                    if np.isfinite(array).all():
                        return array
                except (TypeError,ValueError):
                    pass
            return [arrays(v) for v in value]
        return value

    # This is reuse of numerical values, never evidence of derivative accuracy.
    return (float(candidate['logmass']),gradient,hessian,arrays(integration)),independent
