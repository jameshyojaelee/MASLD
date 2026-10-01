"""Independent C10 recount using approved public and cohort-summary tables only."""
import csv
import hashlib
import json
import math
import os
import platform
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
OBS = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z'
COL = OBS / 'coloc_direction/full-20260930T120450Z'
QC = OBS / 'qc/r1-20260930T121200Z'
inputs = {}


def read(path):
    assert 'rna_genotypes_restricted' not in str(path)
    inputs[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def number(row, key):
    try:
        return float(row[key])
    except (ValueError, KeyError):
        return math.nan


def true(value):
    return value.lower() == 'true'


def groups(rows, key='release'):
    result = defaultdict(list)
    for row in rows:
        result[row[key]].append(row)
    return result


trait = read(COL / 'direction_by_trait.tsv')
pair = read(COL / 'direction_signal_pairs.tsv')
tags = read(COL / 'tag_ld_gene_summary.tsv')
old = read(OBS / 'coloc_direction/full-20260929T140107Z/direction_by_trait.tsv')
old_by_key = {(r['release'], r['ensembl'], r['trait']): int(r['direction_sign']) for r in old}
coloc = {}
for release, rows in groups(trait).items():
    directional = [r for r in rows if int(r['direction_sign'])]
    by_gene = groups(rows, 'ensembl')
    gene_patterns = Counter()
    for gene_rows in by_gene.values():
        signs = {int(r['direction_sign']) for r in gene_rows} - {0}
        gene_patterns['both' if signs == {-1, 1} else 'increase_only' if signs == {1}
                      else 'decrease_only' if signs == {-1} else 'none'] += 1
    pairs = [r for r in pair if r['release'] == release]
    tag_rows = [r for r in tags if r['release'] == release]
    coloc[release] = {
        'trait_rows': len(rows), 'genes': len(by_gene),
        'signs': dict(Counter(r['direction_sign'] for r in rows)),
        'duplicate_gene_trait_keys': len(rows) - len({(r['ensembl'], r['trait']) for r in rows}),
        'gene_direction_patterns': dict(gene_patterns),
        'directional_both_marginals_lt_1e3': sum(number(r, 'gwas_p') < 1e-3 and
                                                number(r, 'eqtl_p_hit2') < 1e-3 for r in directional),
        'directional_both_marginals_lt_1e5': sum(number(r, 'gwas_p') < 1e-5 and
                                                number(r, 'eqtl_p_hit2') < 1e-5 for r in directional),
        'directional_eqtl_fit_missing': sum(not true(r['eqtl_susie_fit_found']) for r in directional),
        'directional_eqtl_conditional_sign_not_true': sum(not true(r['eqtl_conditional_sign_equals_marginal'])
                                                         for r in directional),
        'signal_pairs': len(pairs),
        'ineligible_under_corrected_rule': sum(r['eligible_corrected_rule'] == 'FALSE' for r in pairs),
        'primary_tag_gene_flags': dict(Counter(r['gene_tag_r2_ge_0.8__t3prime_v2_kept'] for r in tag_rows)),
        'primary_tag_genes_with_tags': sum(number(r, 'gene_n_tags__t3prime_v2_kept') > 0 for r in tag_rows),
        'v3_v5_sign_flips': [r['gene'] + ':' + r['trait'] for r in rows
                            if int(r['direction_sign']) * old_by_key.get((release, r['ensembl'], r['trait']), 0) < 0],
        'n4bp2l2_rows': [r for r in rows if r['gene'] == 'N4BP2L2'],
    }

nonpal = [r for r in pair if len(r['gwas_allele1']) == len(r['gwas_allele2']) == 1
          and {r['gwas_allele1'], r['gwas_allele2']} not in ({'A', 'T'}, {'C', 'G'})
          and r['gwas_allele1_follows_file_convention'] != '']
indels = [r for r in pair if max(len(r['gwas_allele1']), len(r['gwas_allele2'])) > 1
          and r['gwas_allele1_follows_file_convention'] != '']
coloc['allele_convention'] = {
    'nonpal_snv_pairs': len(nonpal),
    'nonpal_snv_agree': sum(true(r['gwas_allele1_follows_file_convention']) for r in nonpal),
    'indel_pairs': len(indels),
    'indel_disagree': sum(not true(r['gwas_allele1_follows_file_convention']) for r in indels),
}

cal = read(QC / 'cohort_chrx_calibration.tsv')
model = read(QC / 'cohort_model_inputs.tsv')
t2b = read(ROOT / 'Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-t2b-20260929T220019Z/gateA_cohort_inputs.tsv')
t2b_by_key = {(r['cohort'], r['set']): number(r, 'omega_c') for r in t2b}
qc = {'calibration': {}, 'cohort_model_inputs': model}
for filt, rows in groups(cal, 'site_filter').items():
    testable = [r for r in rows if true(r['testable'])]
    ratios = [number(r, 'ratio_e_model') for r in testable]
    qc['calibration'][filt] = {
        'testable_cells': len(testable), 'ratio_min': min(ratios), 'ratio_max': max(ratios),
        'above_2': sum(x > 2 for x in ratios),
        'max_abs_ratio_recount_difference': max(abs(number(r, 'ratio_e_model') -
            number(r, 'observed_het') / number(r, 'predicted_het_e_model')) for r in testable),
        'zero_event_zero_upper_ci': [r['cohort'] + ':' + r['set'] for r in testable
                                   if number(r, 'observed_het') == 0 and number(r, 'ratio_boot_hi95') == 0],
    }
qc['omega_max_abs_diff_t2b'] = max(abs(number(r, 'omega_c') -
    t2b_by_key[r['cohort'], r['rule_set']]) for r in model)
qc['e_own_cohorts'] = sorted({r['cohort'] for r in model if true(r['e_i_per_library_allowed'])})
qc['ffpe_own_cohorts'] = sorted({r['cohort'] for r in model if true(r['z_ffpe_per_library_allowed'])})
qc['c_own_cohorts'] = sorted({r['cohort'] for r in model if true(r['c_i_per_library_allowed'])})
for r in model:
    for prefix in ('unit_development', 'unit_all'):
        stat, lo, hi = (number(r, 'g1_ffpe_' + prefix + '_' + k) for k in ('stat', 'lo95', 'hi95'))
        n = number(r, 'g1_ffpe_' + prefix + '_sites')
        if r['g1_ffpe_' + prefix + '_verdict'] in ('NA', 'not_computed'):
            assert not math.isfinite(n) or n == 0, (r['cohort'], prefix, 'missing verdict with positive n')
            continue
        verdict = 'undecidable' if n < 10 or not math.isfinite(lo) or (hi-lo)/2 > .02 else (
            'pass' if abs(stat) <= .02 else 'fail')
        assert verdict == r['g1_ffpe_' + prefix + '_verdict'], (r['cohort'], prefix, verdict)
assert qc['omega_max_abs_diff_t2b'] < 1e-12
assert coloc['adopted_2026-08-17']['signs'] == {'-1': 220, '0': 191, '1': 217}
assert coloc['candidate_r2_20260924T053244Z']['signs'] == {'-1': 213, '0': 164, '1': 212}
result = {'coloc': coloc, 'qc': qc, 'inputs_sha256': inputs, 'python': platform.python_version(),
          'scope': 'Public COLOC directions and cohort QC summaries only; no allele outcomes.'}
out = Path(os.environ['CODEX_REC_OUTPUT']) / ('round2_summaries_' + os.environ['SLURM_JOB_ID'] + '.json')
out.write_text(json.dumps(result, indent=2) + '\n')
print('Independent C10 summary recount passed:', out, flush=True)
