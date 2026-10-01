#!/usr/bin/env python3
"""Owner-private null information assembly; emit reduced scalars only.

No private input is bundled or presumed available. The owner supplies a hashed
complete-F raw design/missing masks and its frozen schema. Matrices and projection
coefficients never enter the JSON output or error messages.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
SOURCES = {
    HERE / 'information_weights_v16.py': '92297f26ae83cd7712b61dd4cc4e6e95d51c842ad6528cab7206bbba54bd5359',
    HERE / 'optimistic_information_v16.py': 'b42b9bf09b0b2b5bf594fbe17b0828694714349727a3cfd1c572d14292e04888',
    HERE / 'information_projection.py': '5cd94a34d5bc3d8bd33d5addb493527374477e973774cda337f7d504f4297153',
    HERE.parent / 'design_v1.py': 'cb3163896f982f39ba94cdbf227a8e7d4c7f85a0a646b26542868e6fa14c9c5e',
    HERE.parent / 'likelihood_v1.py': '41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6',
    ROOT / 'docs/technical/agent_exchange/2026-09-29_model_A_likelihood_spec_v1.md': 'd9b38061a542cce371c09b943a97a43733cd685aeebb770f389e70451d6b3f7a',
}
ROLES = ['b', 'd', 'ffpe', 'sex', 'age10', 'axis_EUR_AFR', 'axis_EUR_EAS']
COHORT_COUNTS = {'GSE130970':54,'GSE135251':160,'GSE162694':83,'GSE213621':208,'GSE240729':53}
RULES = 'complete_F_cohort_mean_fill_all_missing_zero_partial_indicators_center_no_SD_global_zero_drop'
SCENARIOS = ('geuvadis_rho', 'rho_0.02')
THETAS = (0., .1, .2, .3)
FROZEN_TABLE_SHA256 = {
    'weights_path':'027aa5d702825ee49a81595bf53a796e696526404b1763d4d3176e39e7229218',
    'eligibility_path':'043ded253ab93f7aa6594beba360f136b4d8ec7300aa895abd3edef86adfdaf5',
}


class SafeFailure(Exception):
    """Fixed error code only; no private values, paths or tracebacks."""


def check(ok, code):
    if not ok:
        raise SafeFailure(code)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def _np():
    import numpy
    return numpy


def _projection():
    path = HERE / 'information_projection.py'
    check(sha(path) == SOURCES[path], 'E_SOURCE')
    spec = importlib.util.spec_from_file_location('_owner_private_projection', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.project_information


def _rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        check(reader.fieldnames and len(set(reader.fieldnames)) == len(reader.fieldnames), 'E_TABLE')
        for row in reader:
            check(None not in row and all(v is not None for v in row.values()), 'E_TABLE')
            yield row


def _center_private_raw(people, schema):
    np = _np()
    check(schema['input_mode'] == 'raw_complete_F_with_missing_masks'
          and schema['raw_roles'] == ROLES and schema['raw_count'] == 7
          and schema['centering_rule'] == RULES, 'E_SCHEMA')
    columns, masks = schema['raw_columns'], schema['missing_columns']
    check(len(columns) == len(masks) == 7 and len(set(columns+masks)) == 14
          and all(isinstance(n,str) and n for n in columns+masks)
          and not set(columns+masks) & {'individual_id','cohort','S','S_tilde','b_tilde','d_tilde'}, 'E_SCHEMA')
    raw, missing, stages, cohorts = [], [], [], []
    for row in people:
        check(set(row) == {'individual_id','cohort','S',*columns,*masks}, 'E_RAW_FIELDS')
        mask = [int(row[n]) for n in masks]
        check(all(v in (0,1) for v in mask), 'E_MISSING_MASK')
        values = [float('nan') if m else float(row[n]) for n,m in zip(columns,mask)]
        check(all(m or np.isfinite(v) for m,v in zip(mask,values)), 'E_RAW_FINITE')
        raw.append(values); missing.append(mask); stages.append(int(row['S'])); cohorts.append(row['cohort'])
    raw, missing = np.asarray(raw,float), np.asarray(missing,bool)
    cohort = np.asarray(cohorts)
    check(raw.shape == (len(people),7), 'E_RAW_SHAPE')
    indicators = [j for j in range(7) if any(missing[cohort==c,j].any() and not missing[cohort==c,j].all() for c in set(cohorts))]
    check(indicators == schema['indicators'], 'E_INDICATOR_SCHEMA')
    filled = raw.copy()
    for c in sorted(set(cohorts)):
        mask = cohort == c
        for j in range(7):
            observed = mask & ~missing[:,j]
            filled[mask & missing[:,j],j] = filled[observed,j].mean() if observed.any() else 0.
    full = np.column_stack((filled,missing[:,indicators].astype(float)))
    stage = np.asarray(stages,float)
    check(np.isin(stage,(0,1,2)).all(), 'E_STAGE')
    for c in sorted(set(cohorts)):
        mask = cohort == c
        constant = np.ptp(full[mask],axis=0) == 0
        full[mask] -= full[mask].mean(axis=0)
        full[np.ix_(mask,constant)] = 0.
        stage[mask] -= stage[mask].mean()
    keep = np.flatnonzero(np.any(full != 0,axis=0)).tolist()
    check(keep == schema['keep'], 'E_KEEP_SCHEMA')
    names = columns + [columns[j]+'_missing' for j in indicators]
    check([names[j] for j in keep] == schema['retained_columns'], 'E_RETAINED_SCHEMA')
    centered = []
    for i,row in enumerate(people):
        result = dict(individual_id=row['individual_id'],cohort=row['cohort'],S=stages[i],
                      S_tilde=stage[i],b_tilde=full[i,0],d_tilde=full[i,1])
        result.update({names[j]:full[i,j] for j in keep})
        centered.append(result)
    return centered


def _design(people, schema, expected_people=558, expected_cohorts=5):
    np = _np()
    check(schema['raw_roles'] == ROLES and schema['raw_count'] == 7
          and schema['centering_rule'] == RULES, 'E_SCHEMA')
    raw = schema['raw_columns']
    indicators, keep = schema['indicators'], schema['keep']
    check(len(raw) == len(set(raw)) == 7 and all(isinstance(n, str) and n for n in raw), 'E_SCHEMA')
    check(indicators == sorted(set(indicators)) and all(isinstance(i, int) and 0 <= i < 7 for i in indicators), 'E_SCHEMA')
    full = raw + [raw[i] + '_missing' for i in indicators]
    check(keep == sorted(set(keep)) and all(isinstance(i, int) and 0 <= i < len(full) for i in keep), 'E_SCHEMA')
    columns = schema['retained_columns']
    check(columns == [full[i] for i in keep] and len(set(full)) == len(full), 'E_SCHEMA')
    check(len(people) == expected_people and len({r['individual_id'] for r in people}) == expected_people, 'E_ROSTER')
    cohorts = sorted({r['cohort'] for r in people})
    check(len(cohorts) == expected_cohorts, 'E_ROSTER')
    if expected_people == 558:
        check({c:sum(r['cohort']==c for r in people) for c in cohorts} == COHORT_COUNTS, 'E_COMPLETE_F_COHORT_COUNTS')
    index = {}
    arrays = []
    for r in people:
        check(set(r) == {'individual_id','cohort','S','S_tilde','b_tilde','d_tilde', *columns}, 'E_DESIGN_FIELDS')
        check(int(r['S']) in (0, 1, 2), 'E_STAGE')
        values = np.array([float(r[n]) for n in ['S_tilde','b_tilde','d_tilde', *columns]])
        check(np.isfinite(values).all(), 'E_DESIGN_FINITE')
        index[r['individual_id']] = (r['cohort'], int(r['S']), values)
        arrays.append(values)
    all_values = np.stack(arrays)
    largest_centering = 0.
    for c in cohorts:
        mask = np.array([r['cohort'] == c for r in people])
        values = all_values[mask]
        scale = np.maximum(1., np.max(np.abs(values), axis=0))
        residual = float(np.max(np.abs(values.mean(0)) / scale))
        largest_centering = max(largest_centering, residual)
        check(residual <= 1e-10, 'E_CENTERING')
        stage = np.array([int(r['S']) for r in people if r['cohort'] == c], float)
        check(np.allclose(values[:, 0], stage-stage.mean(), rtol=0, atol=1e-10), 'E_STAGE_CENTERING')
        # Constant-within-cohort coordinates must be exactly zero.
        for j in range(values.shape[1]):
            if np.ptp(values[:, j]) == 0:
                check(values[0, j] == 0, 'E_CONSTANT_CENTERING')
    check(not columns or np.all(np.any(all_values[:, 3:] != 0, axis=0)), 'E_GLOBAL_ZERO_DROP')
    for j in (0, 1):
        if j in keep:
            check(np.array_equal(all_values[:, 1+j], all_values[:, 3+keep.index(j)]), 'E_DISPERSION_Z')
        else:
            check(np.all(all_values[:, 1+j] == 0), 'E_DROPPED_DISPERSION_Z')
    return index, cohorts, len(columns), largest_centering


def _assemble(people, schema, weight_rows, gene_rows, expected_people=558, expected_cohorts=5, expected_rows_per_scenario=95983, authoritative_roster=None):
    np = _np()
    check(authoritative_roster is not None and len(authoritative_roster) == expected_people
          and all(set(r)=={'individual_id','cohort','S'} for r in authoritative_roster), 'E_AUTHORITATIVE_ROSTER')
    reference = {(r['individual_id'],r['cohort'],int(r['S'])) for r in authoritative_roster}
    check(len(reference)==expected_people
          and {(r['individual_id'],r['cohort'],int(r['S'])) for r in people} == reference, 'E_ROSTER_IDENTITY')
    people = _center_private_raw(people, schema)
    index, cohorts, p, center_error = _design(people, schema, expected_people, expected_cohorts)
    eligibility = {}
    for row in gene_rows:
        key = (row['rho_spec'], row['gene_id'])
        check(key not in eligibility and key[0] in SCENARIOS, 'E_ELIGIBILITY_JOIN')
        values = tuple(float(row[n]) for n in ('alpha_gtex_ln_afc','alpha_plant_capped','expected_included'))
        check(np.isfinite(values).all() and values[2] >= 0, 'E_ELIGIBILITY_VALUE')
        eligibility[key] = values
    check({k[0] for k in eligibility} == set(SCENARIOS), 'E_SCENARIOS')
    local = {}
    seen = set()
    counts = {s: 0 for s in SCENARIOS}
    inclusion = {k: 0. for k in eligibility}
    min_w_eigenvalue = 0.
    for row in weight_rows:
        check(set(row) == {'rho_spec','cohort','individual_id','gene_id','alpha_plant','s_t','S','n','inclusion_probability','W_eta_eta','W_eta_r','W_r_r'}, 'E_WEIGHT_SCHEMA')
        scenario, gene, person = row['rho_spec'], row['gene_id'], row['individual_id']
        key = (scenario, person, gene)
        check(key not in seen and person in index and (scenario,gene) in eligibility, 'E_WEIGHT_JOIN')
        seen.add(key)
        c, stage, values = index[person]
        check(row['cohort'] == c and int(row['S']) == stage and int(row['s_t']) in (-1,1) and int(row['n']) > 0, 'E_WEIGHT_JOIN')
        alpha, probability = float(row['alpha_plant']), float(row['inclusion_probability'])
        check(np.isfinite(alpha) and np.isfinite(probability) and 0 <= probability <= 1
              and np.isclose(alpha, eligibility[(scenario,gene)][1], rtol=0, atol=1e-12), 'E_WEIGHT_PLANT')
        w = np.array([[float(row['W_eta_eta']),float(row['W_eta_r'])],
                      [float(row['W_eta_r']),float(row['W_r_r'])]])
        check(np.isfinite(w).all(), 'E_WEIGHT_PSD')
        mineig = float(np.linalg.eigvalsh(w)[0])
        min_w_eigenvalue = min(min_w_eigenvalue, mineig)
        check(mineig >= -1e-10*max(1., float(np.max(np.abs(w)))), 'E_WEIGHT_PSD')
        # W already includes inclusion_probability; do not multiply again.
        jac = np.zeros((2, p+7))
        jac[0,0], jac[1,1] = 1., 1.
        jac[0,2] = alpha*values[0]
        jac[0,3:3+p] = alpha*values[3:]
        jac[0,3+p] = int(row['s_t'])*values[0]
        jac[1,-3:] = values[:3]
        blockkey = (scenario,c,gene)
        if blockkey not in local:
            local[blockkey] = np.zeros((p+7,p+7))
        local[blockkey] += jac.T @ w @ jac
        counts[scenario] += 1
        inclusion[(scenario,gene)] += probability
    check(all(n == expected_rows_per_scenario for n in counts.values()), 'E_WEIGHT_CENSUS')
    for key, values in eligibility.items():
        check(np.isclose(inclusion[key], values[2], rtol=1e-9, atol=1e-8), 'E_INCLUSION_SUM')
    return index, cohorts, p, center_error, local, eligibility, counts, min_w_eigenvalue


def _two_stage(local, scenario, genes, cohorts, p, projection):
    np = _np()
    shared = np.zeros((p+5,p+5))
    rank = 0
    for gene in genes:
        block = sum((local.get((scenario,c,gene), np.zeros((p+7,p+7))) for c in cohorts), np.zeros((p+7,p+7)))
        intercept = block[:2,:2]
        scale = np.sqrt(np.maximum(np.diag(intercept), 0))
        active = scale > 0
        inv = np.zeros((2,2))
        if active.any():
            a = intercept[np.ix_(active,active)] / np.outer(scale[active],scale[active])
            ev, vec = np.linalg.eigh(a)
            retained = ev > 1e-10*max(1., float(ev[-1]))
            rank += int(retained.sum())
            normalized_inverse = (vec[:,retained]/ev[retained]) @ vec[:,retained].T
            inv[np.ix_(active,active)] = normalized_inverse / np.outer(scale[active],scale[active])
        cross = block[:2,2:]
        check(np.allclose(intercept @ inv @ cross, cross, rtol=1e-8, atol=1e-8*max(1., float(np.max(np.abs(block))))), 'E_SINGULAR_GENE_COUPLING')
        shared += block[2:,2:] - cross.T @ inv @ cross
    result = projection({'pooled_private':shared}, 0)
    return result['information'], rank


def _calculate(people, schema, weight_rows, gene_rows, expected_people=558, expected_cohorts=5, expected_rows_per_scenario=95983, authoritative_roster=None):
    np = _np()
    values = _assemble(people, schema, weight_rows, gene_rows, expected_people, expected_cohorts, expected_rows_per_scenario, authoritative_roster)
    index, cohorts, p, center_error, local, eligibility, counts, min_w_eigenvalue = values
    projection = _projection()
    results = []
    for scenario in SCENARIOS:
        for theta in THETAS:
            genes = sorted(g for (s,g), a in eligibility.items() if s == scenario and a[2] >= 10 and abs(a[0]) >= theta)
            check(bool(genes), 'E_EMPTY_ELIGIBLE')
            gcount = len(genes)
            size = 2*gcount+p+5
            matrices = {c: np.zeros((size,size)) for c in cohorts}
            for j, gene in enumerate(genes):
                coordinates = [j, gcount+j, *range(2*gcount,size)]
                for c in cohorts:
                    block = local.get((scenario,c,gene))
                    if block is not None:
                        matrices[c][np.ix_(coordinates,coordinates)] += block
            pooled = sum(matrices.values(), np.zeros((size,size)))
            projected = projection(matrices, 2*gcount)
            q = projected['efficient_score_coefficients']
            info = projected['information']
            independent, gene_rank = _two_stage(local, scenario, genes, cohorts, p, projection)
            agreement = abs(info-independent)
            check(agreement <= 1e-6+1e-8*max(abs(info),abs(independent)), 'E_TWO_STAGE_AGREEMENT')
            residual = pooled @ q
            nuisance = np.delete(residual, 2*gcount)
            scale = np.sqrt(np.maximum(np.diag(pooled),0))
            active = scale > 0
            min_scaled_eigenvalue = 0.
            if active.any():
                denominator = np.outer(scale[active],scale[active])
                min_scaled_eigenvalue = min(float(np.linalg.eigvalsh(m[np.ix_(active,active)] / denominator)[0]) for m in matrices.values())
            results.append(dict(rho_spec=scenario, theta=theta, I=info,
                V_c=projected['cohort_information'], signed_R_c=projected['cohort_local_response'],
                identifiable=projected['identifiable'], nuisance_rank=projected['nuisance_rank'], gene_intercept_rank=gene_rank,
                eligible_gene_count=gcount, expected_inclusion=sum(eligibility[(scenario,g)][2] for g in genes),
                complete_F_count=len(index), cohort_counts={c:sum(r[0]==c for r in index.values()) for c in cohorts},
                retained_z_count=p, weight_row_counts=counts,
                checks=dict(complete_F_centering_residual=center_error, minimum_row_weight_eigenvalue=min_w_eigenvalue,
                    minimum_normalized_cohort_eigenvalue=min_scaled_eigenvalue,
                    additional_convention_sensitive_two_stage_information_residual=agreement, V_sum_residual=abs(sum(projected['cohort_information'].values())-info),
                    R_sum_residual=abs(sum(projected['cohort_local_response'].values())-info),
                    pooled_nuisance_orthogonality_residual=float(np.max(np.abs(nuisance),initial=0)),
                    pooled_target_information_residual=abs(float(residual[2*gcount])-info)),
                eligible_gene_order_sha256=hashlib.sha256(json.dumps(genes,separators=(',',':')).encode()).hexdigest()))
    return results


def owner_summary(manifest):
    """Owner-only callable: return exactly eight reduced summaries, never q/M."""
    try:
        check(manifest['contract'] == 'owner_private_null_information_v1' and manifest['development_only'] is True, 'E_MANIFEST')
        check(manifest['plant'] == dict(kappa=0,delta=0,omega_S=0,rho_S=0,rho_b=0,rho_d=0,tau=0), 'E_NULL_PLANT')
        check(manifest['weight_definition'] == 'W_ab = P(included) * E[score_a score_b | included], a,b = eta,r'
              and manifest['fixed_weight_inputs_attested'] is True, 'E_WEIGHT_DEFINITION')
        root = Path(manifest['owner_private_root']).resolve()
        private = [Path(manifest[n]).resolve() for n in ('design_path','schema_path')]
        check(all(p.is_relative_to(root) and not p.is_relative_to(ROOT.resolve()) for p in private), 'E_PRIVATE_BOUNDARY')
        paths = {n:Path(manifest[n]) for n in ('design_path','schema_path','weights_path','eligibility_path','authoritative_roster_path')}
        check(paths['authoritative_roster_path'].resolve() != paths['design_path'].resolve()
              and manifest['roster_source'] == 'owner-frozen authoritative complete development F v1.4', 'E_AUTHORITATIVE_ROSTER')
        expected = manifest['input_sha256']
        check(set(expected) == set(paths), 'E_INPUT_HASH_SCHEMA')
        check(all(expected[n] == h for n,h in FROZEN_TABLE_SHA256.items()), 'E_FROZEN_NULL_TABLE')
        for n,p in paths.items():
            check(sha(p) == expected[n], 'E_INPUT_HASH')
        for p,h in SOURCES.items():
            check(sha(p) == h, 'E_SOURCE')
        schema = json.loads(paths['schema_path'].read_text())
        result = _calculate(list(_rows(paths['design_path'])), schema,
                            _rows(paths['weights_path']), _rows(paths['eligibility_path']),
                            authoritative_roster=list(_rows(paths['authoritative_roster_path'])))
        for n,p in paths.items():
            check(sha(p) == expected[n], 'E_INPUT_CHANGED')
        for p,h in SOURCES.items():
            check(sha(p) == h, 'E_SOURCE_CHANGED')
        return dict(contract='owner_private_null_information_v1', results=result,
                    input_sha256=expected, schema_sha256=expected['schema_path'],
                    source_sha256={p.name:h for p,h in SOURCES.items()},
                    helper_sha256=sha(Path(__file__)),
                    limits='Owner-reviewed reduced null planning summaries only; not formal privacy, lambda/MDE/P2-width/calibration or publication approval.')
    except SafeFailure:
        raise
    except Exception:
        raise SafeFailure('E_PRIVATE_VALIDATION') from None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true')
    parser.add_argument('--owner-manifest')
    parser.add_argument('--output')
    args = parser.parse_args()
    if args.inspect:
        check(not args.owner_manifest and not args.output, 'E_INSPECT_SCOPE')
        check(all(sha(p)==h for p,h in SOURCES.items()), 'E_SOURCE')
        print(json.dumps(dict(contract='owner_private_null_information_v1', required_owner_inputs=['private_complete_raw558_design_and_missing_masks','frozen_schema_with_physical_names','hashed_authoritative_complete_F_roster','public_weights','public_eligibility','input_hash_manifest'],
              schema_raw_roles=ROLES, centering_rule=RULES, emitted='eight reduced scalar scenario/theta summaries only', actual_private_inputs_available=False,
              source_sha256={p.name:h for p,h in SOURCES.items()}, helper_sha256=sha(Path(__file__))), indent=2))
        return
    try:
        check(os.environ.get('SLURM_JOB_ID','').isdigit() and args.owner_manifest and args.output, 'E_OWNER_COMPUTE')
        result = owner_summary(json.loads(Path(args.owner_manifest).read_text()))
        # Explicit allowlisted scalar result only. Never dump caller manifest.
        with Path(args.output).open('x') as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write('\n')
    except Exception:
        # Suppress private exception values, paths, frames and NumPy assertions.
        print(json.dumps(dict(status='owner_private_validation_failed')),file=__import__('sys').stderr)
        raise SystemExit(1) from None


if __name__ == '__main__':
    main()
