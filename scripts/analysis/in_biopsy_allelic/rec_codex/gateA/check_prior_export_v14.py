"""Audit approved planning genotype probabilities without allelic outcomes."""
import csv
import gzip
import json
import os
from collections import Counter
from pathlib import Path

BASE=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions')
path=BASE/'model-a-gateA-design-20260930T121812Z/design_rows.tsv.gz'
counts=Counter();bad=Counter();examples={};maximum=0.;total=0
with gzip.open(path,'rt') as fh:
    for r in csv.DictReader(fh,delimiter='\t'):
        if r['set']!='development':continue
        values=[float(r[k]) for k in ('P_H','P_R','P_A')]
        weight=sum(values);counts[r['cohort']]+=1;total+=1;maximum=max(maximum,weight)
        if weight>1+1e-10 or min(values)<0:
            bad[r['cohort']]+=1
            key=(r['cohort'],r['gene_id'],r['chrom'],r['pos'])
            if key not in examples:
                examples[key]={'cohort':r['cohort'],'gene':r['gene_id'],'chrom':r['chrom'],
                    'pos':r['pos'],'P_H':values[0],'P_R':values[1],'P_A':values[2],'sum':weight}
result={'path':str(path),'development_rows':total,'rows_by_cohort':dict(counts),
    'invalid_rows_by_cohort':dict(bad),'invalid_rows':sum(bad.values()),'maximum_prior_sum':maximum,
    'distinct_cohort_gene_tag_invalid':len(examples),'invalid_examples':list(examples.values()),
    'scope':'Approved public-prior and total-depth design export; no outcomes or sealed analysis.'}
out=Path(os.environ['CODEX_REC_OUTPUT'])/('prior_export_v14_'+os.environ['SLURM_JOB_ID']+'.json')
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='invalid_examples'},indent=2),flush=True)
print('First invalid examples:',json.dumps(list(examples.values())[:5]),flush=True)
