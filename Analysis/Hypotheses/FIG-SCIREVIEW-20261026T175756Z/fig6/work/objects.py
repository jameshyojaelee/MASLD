import json, collections
F='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/molecular_object_sidecar.jsonl'
c=collections.Counter(); cr=collections.Counter(); keys=None
for l in open(F):
    r=json.loads(l); keys=keys or list(r.keys())
    c[r['candidate_object']]+=1
    if r.get('symbol') in ('GLP1R','HKDC1','MTARC1','THRB'):
        print({k:r[k] for k in r if k in ('symbol','candidate_object','genetics_call','rna_call','primary_evidence_class','recommended_experiment_rule_id','legacy_next_experiment_rule_id','gene_biotype')})
print(c); print([k for k in keys if 'rule' in k or 'route' in k or 'symbol' in k])
