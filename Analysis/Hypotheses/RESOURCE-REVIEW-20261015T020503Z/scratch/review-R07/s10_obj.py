import json, collections
F='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/molecular_object_sidecar.jsonl'
rc=collections.Counter(); oc=collections.Counter(); cross=collections.Counter(); n=0; ids=set(); keys=None
sel={}
for line in open(F):
    d=json.loads(line); n+=1
    if keys is None: keys=list(d.keys())
    rid=d.get('recommended_experiment_rule_id') or d.get('route') 
    rc[rid]+=1; oc[d.get('candidate_object')]+=1; cross[(d.get('primary_evidence_class'),d.get('candidate_object'),rid)]+=1
    ids.add(d.get('ensembl_id'))
    if d.get('symbol') in ('PNPLA3','MTARC1','TM6SF2','HKDC1','THRB','GLP1R','GNMT'): sel[d.get('symbol')]={k:d.get(k) for k in ('candidate_object','open_mechanistic_question','recommended_experiment_rule_id','genetics_call','rna_call','primary_evidence_class')}
print(n,len(ids)); print(keys); print(rc); print(oc)
for k,v in sorted(cross.items(), key=lambda x:-x[1]): print(k,v)
for k,v in sel.items(): print(k,v)
