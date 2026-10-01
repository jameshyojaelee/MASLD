"""Development expected inclusion from approved total depths and public GTEx effects.

No called-heterozygote counts or allele fractions are inputs. These are planning
expectations conditional on the exported cohort-mixture genotype priors.
"""
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import statistics
import sys
from collections import Counter,defaultdict
from functools import lru_cache
from pathlib import Path
import numpy as np
import scipy
from scipy.special import expit,logsumexp

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from likelihood_v1 import bb_logmass,hom_log_inclusion

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
BASE=ROOT/'Analysis/MASLD_Model_Benchmark/executions'
DESIGN=BASE/'model-a-gateA-design-20260930T121812Z'
RULES=BASE/'in-biopsy-allelic-t2b-20260930T123844Z'
QC=BASE/'in-biopsy-allelic-obs-20260926T155651Z/qc/r1c-20260930T123747Z'
GTEX=ROOT/'data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz'


def table(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt') as fh:return list(csv.DictReader(fh,delimiter='\t'))


people={r['individual_id']:r for r in table(DESIGN/'design_individuals.tsv') if r['set']=='development'}
assert len(people)==558 and len({r['cluster_id'] for r in people.values()})==558
cohorts=sorted({r['cohort'] for r in people.values()})
rules={r['cohort']:tuple(int(r[k]) for k in ('min_dp','min_minor_reads','min_minor_pct'))
       for r in table(RULES/'cohort_rules.tsv') if r['set']=='development' and r['cohort'] in cohorts}
inputs={r['cohort']:r for r in table(RULES/'gateA_cohort_inputs.tsv')
        if r['set']=='development' and r['cohort'] in cohorts}
cal={r['cohort']:r for r in table(QC/'cohort_model_inputs.tsv')
     if r['rule_set']=='development' and r['cohort'] in cohorts}
leads={r['gene_id']:r['lead_variant_id'] for r in table(RULES/'tag_table/tag_table.tsv')
       if r['kept_primary']=='True'}
afc={r['gene_id']:r for r in table(GTEX) if r['gene_id'] in leads}
assert set(afc)==set(leads)
assert all(afc[g]['variant_id']==leads[g] for g in leads)
alpha={g:float(afc[g]['log2_aFC'])*math.log(2) for g in leads}
assert all(math.isfinite(v) for v in alpha.values())
plants={g:float(np.clip(v,-3,3)) for g,v in alpha.items()}
rho_specs={'rho_0.02':.02,'geuvadis_rho':float(inputs[cohorts[0]]['rho_het_geuvadis'])}


@lru_cache(maxsize=65536)
def het_call(n,minor,percent,eta,rho):
    lo=max(minor,(percent*n+99)//100)
    if n-lo<lo:return 0.
    mean=expit(eta);concentration=(1-rho)/rho
    return float(np.exp(logsumexp(bb_logmass(np.arange(lo,n-lo+1),n,
                 mean*concentration,(1-mean)*concentration))))


# A symmetric integer filter has the same detection probability after allele reversal.
for n in (10,20,100,1000):
    np.testing.assert_allclose(het_call(n,1,10,.7,.02),het_call(n,1,10,-.7,.02),atol=1e-12)

expected={name:defaultdict(lambda:np.zeros(3)) for name in rho_specs}
cohort_totals={name:defaultdict(lambda:np.zeros(3)) for name in rho_specs}
candidate=Counter();depth_pass=Counter();max_depth=0;subnormalized=0;seen=set()
rounded_prior_rows=0;prior_excess_upper_bound=0.
with gzip.open(DESIGN/'design_rows.tsv.gz','rt') as fh:
    for row in csv.DictReader(fh,delimiter='\t'):
        if row['set']!='development':continue
        person=people[row['individual_id']];c=row['cohort'];g=row['gene_id'];n=int(row['n'])
        assert c==person['cohort'] and g in alpha
        key=(row['individual_id'],g);assert key not in seen;seen.add(key)
        candidate[c]+=1;max_depth=max(max_depth,n)
        d,k,F=rules[c]
        if n<d:continue
        depth_pass[c]+=1
        prior=np.asarray([float(row[v]) for v in ('P_H','P_R','P_A')])
        assert np.isfinite(prior).all() and (prior>=0).all() and prior.sum()<=1+1e-6
        if prior.sum()>1:
            # Independent recount bounded the excess at 9.40e-7. Common
            # rescaling preserves conditional class odds and restores a
            # coherent unconditional expectation at these rounded inputs.
            prior_excess_upper_bound+=prior.sum()-1
            rounded_prior_rows+=1
            prior=prior/prior.sum()
        subnormalized+=int(prior.sum()<1-1e-10)
        s=int(row['s_t']);assert s in (-1,1)
        omega=float(inputs[c]['omega_c']);phi=float(inputs[c]['phi'])
        # v1.4 b exports log(e_model), before v1.6 conservative inflation.
        e=max(math.exp(float(person['z_b_log_e']))*float(cal[c]['e_multiplier_applied']),1e-4)
        assert 0<e<.5
        lo=max(k,(F*n+99)//100)
        hom=math.exp(hom_log_inclusion(n,lo,e,phi)) if lo<=n-lo else 0.
        for name,rho in rho_specs.items():
            true=prior[0]*het_call(n,k,F,plants[g]+s*omega,rho)
            false=(prior[1]+prior[2])*hom
            contribution=np.array([true+false,true,false])
            expected[name][g]+=contribution;cohort_totals[name][c]+=contribution

summaries={}
for name,counts in expected.items():
    thresholds=[]
    for theta in (0.,.1,.2,.3):
        selected=[g for g,v in counts.items() if v[0]>=10 and abs(alpha[g])>=theta]
        thresholds.append({'theta':theta,'eligible_genes':len(selected),
            'expected_rows_eligible':float(sum(counts[g][0] for g in selected)),
            'expected_true_het_rows_eligible':float(sum(counts[g][1] for g in selected)),
            'median_abs_alpha_eligible':statistics.median(abs(plants[g]) for g in selected),
            'planning_tau_max':.125*statistics.median(abs(plants[g]) for g in selected),
            'minimum_distance_expected_rows_to_10':min(abs(v[0]-10) for g,v in counts.items()
                                                        if abs(alpha[g])>=theta)})
    summaries[name]={'rho':rho_specs[name],'genes_with_positive_expectation':int(sum(v[0]>0 for v in counts.values())),
        'total_expected_included':float(sum(v[0] for v in counts.values())),
        'total_expected_true_het':float(sum(v[1] for v in counts.values())),
        'total_expected_hom_false_het':float(sum(v[2] for v in counts.values())),
        'cohort_expected_included_true_false':{c:v.tolist() for c,v in cohort_totals[name].items()},
        'thresholds':thresholds}
paths=[DESIGN/'design_individuals.tsv',DESIGN/'design_rows.tsv.gz',RULES/'cohort_rules.tsv',
       RULES/'gateA_cohort_inputs.tsv',RULES/'tag_table/tag_table.tsv',QC/'cohort_model_inputs.tsv',GTEX]
out=Path(os.environ['CODEX_REC_OUTPUT'])/('expected_counts_v16_'+os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)
with (out/'expected_gene_rows.tsv').open('w') as fh:
    writer=csv.writer(fh,delimiter='\t');writer.writerow(['rho_spec','gene_id','alpha_gtex_ln_afc',
        'alpha_plant_capped','expected_included','expected_true_het','expected_hom_false_het'])
    for name,counts in expected.items():
        for g in sorted(counts):writer.writerow([name,g,alpha[g],plants[g],*counts[g]])
result={'spec':'v1.6','n':558,'participant_counts':dict(Counter(r['cohort'] for r in people.values())),
    'candidate_rows':dict(candidate),'depth_pass_rows':dict(depth_pass),'max_depth':max_depth,
    'rules':rules,'public_primary_tag_genes':len(leads),'depth_passing_rows_subnormalized_prior':subnormalized,
    'depth_passing_rows_tiny_prior_excess_rescaled':rounded_prior_rows,
    'total_expected_rows_change_upper_bound_from_prior_rescaling':prior_excess_upper_bound,
    'summaries':summaries,'inputs_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
    'limits':'Planning genotype priors; capped GTEx effects assumed at chosen tags; rho sensitivity. '
             'No observed allele counts. No Fisher information, MDE, theta choice or Gate A verdict. '
             'Merged priors below one retain excluded X/X mass. Sums above one by <=1e-6 are '
             'commonly rescaled to one; class odds are preserved and the count-change bound is recorded.'}
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'output':str(out),'summaries':summaries}),flush=True)
