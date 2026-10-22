# Count r8 molecular-object sidecar records by object type and tissue/spatial fields.
import json, collections
F='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/molecular_object_sidecar.jsonl'
n=0; obj=collections.Counter(); bio=collections.Counter(); keys=None
spatial=collections.Counter(); reg=0
for line in open(F):
    r=json.loads(line); n+=1
    if keys is None: keys=sorted(r)
    obj[r.get('candidate_object')]+=1
    bio[r.get('gene_biotype')]+=1
    aa=json.loads(r['assay_applicability']) if r.get('assay_applicability') else {}
    spatial[aa.get('spatial')]+=1
    s=json.dumps(r).lower()
    if 'regulatory_dna' in s or 'regulatory dna' in s or 'enhancer_object' in s: reg+=1
print('records',n)
print('candidate_object',dict(obj))
print('gene_biotype top',bio.most_common(6))
print('assay_applicability.spatial',dict(spatial))
print('records mentioning regulatory_dna/regulatory dna/enhancer_object anywhere',reg)
print('keys',keys)
