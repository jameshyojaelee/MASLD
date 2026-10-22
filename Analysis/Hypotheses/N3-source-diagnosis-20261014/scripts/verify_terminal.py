"""Focused metadata, native-input and numerical checks; no models or resampling."""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def bh(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.] * len(values)
    upper = 1.
    for rank in range(len(values), 0, -1):
        index = order[rank - 1]
        upper = min(upper, values[index] * len(values) / rank)
        result[index] = upper
    return result


checks = {}
fits = {}
for tag, sub, eligible in [('raw', 'results', 45), ('corrected', 'results_corrected', 40)]:
    base = HERE / tag
    out = base / sub
    summary = json.loads((out / 'summary.json').read_text())
    targets = read(out / 'target_classes.tsv')
    tested = [r for r in targets if r['testable'].lower() == 'true']
    assert len(targets) == 51 and len(tested) == eligible
    assert summary['n_donors'] == 70 and summary['n_targets'] == 51
    assert summary['permutation']['B'] == 1000 and summary['permutation']['seed'] == 20260923
    permutations = read(out / 'permutation_summary.tsv')
    assert len(permutations) == 1000 and {int(r['perm']) for r in permutations} == set(range(1, 1001))
    nominal = sum(float(r['frac_nominal']) for r in permutations) / 1000
    assert abs(nominal - summary['permutation']['mean_fraction_nominal_p_lt_0.05']) < 1e-6
    assert summary['permutation']['calibration_pass'] == (.03 <= nominal <= .07)
    empirical = {r['symbol']: float(r['emp_p']) for r in read(out / 'permutation_target_p.tsv')}
    for source, field in [('P.Value', 'q_bh'), ('emp_p', 'q_emp')]:
        q = bh([empirical[r['symbol']] if source == 'emp_p' else float(r[source]) for r in tested])
        assert max(abs(v - float(r[field])) for v, r in zip(q, tested)) < 1e-12
    primary = 'q_bh' if summary['permutation']['calibration_pass'] else 'q_emp'
    assert all(abs(float(r['q_primary']) - float(r[primary])) < 1e-12 for r in tested)
    assert summary['positive_control']['pass'] is True
    assert summary['positive_control']['genes_in_model'] == ['GNMT', 'MAT1A', 'CYP2C19']
    original, changed = read(base / 'donor_meta_original.tsv'), read(base / 'donor_meta.tsv')
    differences = [(a['donor'], k) for a, b in zip(original, changed) for k in a if a[k] != b[k]]
    assert set(differences) == {(d, k) for d in ['SRR14602057', 'SRR14602058'] for k in ['disease_stage_coarse', 'stage_ord']}
    roster = read(out / 'donors_in_model.tsv')
    assert len(roster) == 70
    counts = {stage: sum(r['stage'] == stage for r in roster) for stage in ['Healthy', 'Steatosis', 'Steatohepatitis', 'Cirrhosis']}
    assert list(counts.values()) == [14, 11, 34, 11]
    fits[tag] = dict(n_donors=70, stage_counts=counts, n_target_roster=51, n_testable=eligible,
                     ordinary_BH_hits=sum(float(r['q_bh']) < .05 for r in tested),
                     primary_hits=sum(float(r['q_primary']) < .05 for r in tested),
                     calibration_pass=summary['permutation']['calibration_pass'],
                     nominal_false_positive_fraction=nominal,
                     primary_basis=primary, control=summary['positive_control'])
    checks['terminal_' + tag + '_fit'] = True
for r in read(HERE / 'native_input_sha256.tsv'):
    assert hashlib.sha256(Path(r['input']).read_bytes()).hexdigest() == r['sha256']
for r in read(HERE / 'per_library_input_sha256.tsv'):
    assert hashlib.sha256(Path(r['input']).read_bytes()).hexdigest() == r['sha256']
qc = read(HERE / 'corrected/library_ambient_qc.tsv')
assert len(qc) == 208 and all(r['status'] == 'corrected' and r['profile_source'] == 'empty_droplets' for r in qc)
donors = read(HERE / 'corrected/donor_ambient.tsv')
assert len(donors) == 70 and all(float(r['n_passthrough']) == 0 and float(r['n_no_empty_profile']) == 0 for r in donors)
ambient = json.loads((HERE / 'corrected/ambient/summary.json').read_text())
for tag in ['raw', 'corrected']:
    assert max(ambient['refit_reproduction'][tag]) < 1e-8
historical = (HERE / 'historical_nine_target_diagnostic_roster.txt').read_text().splitlines()
assert len(historical) == 9 and set(historical) == {'LEPR', 'CNR1', 'FGF1', 'PDE4A', 'SLC5A1', 'AR', 'NR1H4', 'NR3C2', 'VDR'}
checks.update(terminal_corrected_fit=True, terminal_ambient_summary=True,
              source_only_metadata_changes=True, native_input_hashes_unchanged=True,
              fixed_target_roster_51=True, tested_families_raw_45_corrected_40=True,
              permutations_raw_1000_corrected_1000=True, seed_20260923_preserved=True,
              independent_BH_and_empirical_fallback_checks=True,
              fixed_historical_diagnostic_roster_9=True, controls_preserved_and_pass=True,
              correction_coverage_verified=True)
result = dict(**checks, fits=fits, correction_coverage=dict(libraries=208, donors=70,
              corrected_libraries=208, passthrough_donors=0),
              limitations=['Native empirical p values use saved target permutation summaries; full target permutation t arrays are not saved, so no independent redraw is claimed.',
                           'No powered negative, clinical target effect or therapeutic efficacy is inferred.'])
(HERE / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
