"""Load only the shared synthetic B-FIX fixture, independent of its generator."""
import csv
import json
from pathlib import Path
import numpy as np
from design_v1 import Design,Participant
from likelihood_v1 import Row


ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
RAW_COLUMNS=('z_b_log_e','z_d_dup','z_ffpe','z_sex_female','z_age10','z_c')


def table(name):
    with (FIXTURE/name).open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def load():
    settings=json.loads((FIXTURE/'globals.json').read_text())
    genes={g:i for i,g in enumerate(settings['gene_order'])}
    cohorts={r['cohort']:r for r in table('cohorts.tsv')}
    people=[];clusters={}
    for r in table('individuals.tsv'):
        number=lambda key:float(r[key]) if r[key] not in ('NA','') else float('nan')
        count=int(round(np.expm1(number('stage_law_log1p_candidate_rows'))))
        people.append(Participant(r['individual_id'],r['cohort'],int(r['S']),
                                  tuple(number(c) for c in RAW_COLUMNS),number('e_i'),
                                  float(np.expm1(number('stage_law_log1p_median_n'))),count))
        clusters[r['individual_id']]=r['cluster_id']
    rows=[]
    for r in table('rows.tsv'):
        c=cohorts[r['cohort']];orientation=int(r['orientation'])
        assert int(r['s_t'])==-orientation
        assert int(r['a'])==int(r['alt_reads'] if orientation==1 else r['ref_reads'])
        rows.append(Row(genes[r['gene_id']],r['individual_id'],r['cohort'],int(r['n']),int(r['a']),
                        float(r['p_ref']),float(r['p_alt']),float(r['p_x']),float(r['e_i']),
                        settings['phi'],float(c['omega_c']),orientation,0.,(),0.,0.,
                        int(c['min_dp']),int(c['min_minor_reads']),int(c['min_minor_pct'])))
    design=Design(people,rows,clusters)
    names=list(RAW_COLUMNS)+[RAW_COLUMNS[j]+'_missing' for j in design.schema['indicators']]
    retained=[names[j] for j in design.schema['keep']]
    assert retained==settings['design_columns'],(retained,settings['design_columns'])
    assert len(rows)==173 and len(people)==60
    return design,settings


def parameters(path):
    v=json.loads(path.read_text())
    p=np.r_[v['alpha'],v['r'],v['kappa'],v['delta'],v['omega_S'],v['rho_S'],v['rho_b'],v['rho_d']]
    return p,float(v['v'])
