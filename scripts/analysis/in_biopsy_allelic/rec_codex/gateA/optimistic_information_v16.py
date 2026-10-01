"""Upper information bound without the two unavailable ancestry nuisance axes.

Uses approved development total depths, priors, non-ancestry z and expected
row weights. This cannot select theta or prove final design power.
"""
import csv
import gzip
import hashlib
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import scipy
from scipy.stats import norm
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from design_v1 import Design, Participant
from information_projection import project_information


ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
BASE = ROOT/'Analysis/MASLD_Model_Benchmark/executions'
REC = BASE/'codex-rec-20260929T142434Z'
PERSONS = BASE/'model-a-gateA-design-20260930T121812Z/design_individuals.tsv'
WEIGHTS = REC/'information_weights_v16_21986233/weights_for_full_z.tsv.gz'
GENES = REC/'expected_counts_v16_21986201/expected_gene_rows.tsv'
COLUMNS = ('z_b_log_e','z_d_dup','z_ffpe','z_sex_female','z_age10')


def table(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path,'rt') as fh:
        yield from csv.DictReader(fh,delimiter='\t')


people = []
for row in table(PERSONS):
    if row['set'] != 'development':
        continue
    raw = tuple(float(row[c]) if row[c] not in ('NA','') else float('nan') for c in COLUMNS)
    people.append(Participant(row['individual_id'], row['cohort'], int(row['S']), raw,
                              math.exp(raw[0]), 0., 0))
assert len(people) == 558
assert len({p.individual for p in people}) == 558
design = Design(people, [], gene_count=0)
person_index = {p.individual:j for j,p in enumerate(design.participants)}
z_count = design.z.shape[1]
local_size = z_count+7  # alpha_g, r_g, then all shared coordinates.
local = defaultdict(lambda:np.zeros((local_size,local_size)))
written = 0
for row in table(WEIGHTS):
    j = person_index[row['individual_id']]
    person = design.participants[j]
    assert person.cohort == row['cohort'] and person.stage == int(row['S'])
    stage = design.stage_centered[j]
    alpha = float(row['alpha_plant'])
    jac = np.zeros((2,local_size))
    jac[0,0] = 1
    jac[0,2] = alpha*stage
    jac[0,3:3+z_count] = alpha*design.z[j]
    jac[0,3+z_count] = int(row['s_t'])*stage
    jac[1,1] = 1
    jac[1,-3:] = [stage,design.background[j],design.duplication[j]]
    w = np.array([[float(row['W_eta_eta']),float(row['W_eta_r'])],
                  [float(row['W_eta_r']),float(row['W_r_r'])]])
    local[(row['rho_spec'],person.cohort,row['gene_id'])] += jac.T@w@jac
    written += 1
assert written == 191966
gene_rows = list(table(GENES))
cohorts = sorted({p.cohort for p in people})
out = Path(os.environ['CODEX_REC_OUTPUT'])/('optimistic_information_v16_'+os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)
results = []
for rho_spec in sorted({row['rho_spec'] for row in gene_rows}):
    for theta in (0.,.1,.2,.3):
        genes = sorted(row['gene_id'] for row in gene_rows if row['rho_spec'] == rho_spec
                       and float(row['expected_included']) >= 10
                       and abs(float(row['alpha_gtex_ln_afc'])) >= theta)
        g_count = len(genes)
        size = 2*g_count+z_count+5
        matrices = {c:np.zeros((size,size)) for c in cohorts}
        for g,gene in enumerate(genes):
            coordinates = [g,g_count+g,*range(2*g_count,size)]
            for cohort in cohorts:
                block = local.get((rho_spec,cohort,gene))
                if block is not None:
                    matrices[cohort][np.ix_(coordinates,coordinates)] += block
        projected = project_information(matrices, 2*g_count)
        assert projected['identifiable']
        information = projected['information']
        # Independent two-stage Schur derivation: first remove each gene's
        # alpha/r intercepts, then remove shared nuisances in the small matrix.
        shared = np.zeros((local_size-2,local_size-2))
        for gene in genes:
            gene_total = sum((local.get((rho_spec,c,gene),np.zeros((local_size,local_size)))
                              for c in cohorts),np.zeros((local_size,local_size)))
            intercept = gene_total[:2,:2]
            cross = gene_total[:2,2:]
            shared += gene_total[2:,2:] - cross.T@np.linalg.solve(intercept,cross)
        schur = project_information({'all':shared},0)
        np.testing.assert_allclose(information,schur['information'],atol=1e-6,rtol=1e-9)
        se = 1/math.sqrt(information)
        result = {'rho_spec':rho_spec,'theta':theta,'eligible_genes':g_count,
                  'non_ancestry_covariates':z_count,'information_upper_bound':information,
                  'standard_error_lower_bound':se,
                  'normal_80pct_detection_lower_forecast_nominal_0.05':float((norm.ppf(.975)+norm.ppf(.8))*se),
                  'normal_80pct_detection_lower_forecast_alpha_0.05_over_3':float((norm.ppf(1-.05/6)+norm.ppf(.8))*se),
                  'cohort_variance_information':projected['cohort_information'],
                  'cohort_local_slope_response':projected['cohort_local_response'],
                  'nuisance_rank':projected['nuisance_rank']}
        result['independent_two_stage_information'] = schur['information']
        results.append(result)
        (out/(rho_spec+'_theta_'+str(theta)+'.json')).write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result),flush=True)
summary = {'biological_n':558,'results':results,
           'centering':'Complete development F per cohort, including participants without an included row.',
           'plant':'kappa=delta=omega_S=dispersion_slopes=tau=0; capped signed GTEx alpha, rho as named.',
           'missing_nuisance_columns':['axis_EUR_AFR','axis_EUR_EAS'],
           'numpy':np.__version__,'scipy':scipy.__version__,
           'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (PERSONS,WEIGHTS,GENES)},
           'limits':'Upper bound on expected information in the specified planning model. Adding ancestry nuisance columns cannot improve information at this same plant. Normal detection forecasts are asymptotic approximations, not CRT power, P2 width, adjusted-effect power or a theta-selection result.'}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
