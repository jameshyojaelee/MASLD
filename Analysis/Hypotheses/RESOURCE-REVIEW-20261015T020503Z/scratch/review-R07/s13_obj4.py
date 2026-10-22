import json, collections
F='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/molecular_object_sidecar.jsonl'
C=collections.Counter(); ex={}
for line in open(F):
    d=json.loads(line); a=json.loads(d['assay_applicability'])
    C[(d['rna_call'],a['bulk_rna'],a['proteomics'],d['gene_biotype'])]+=1
    if d['symbol'] in ('GLP1R','HKDC1','THRB'): ex[d['symbol']]=a
for k,v in sorted(C.items()): print(k,v)
print(ex)
