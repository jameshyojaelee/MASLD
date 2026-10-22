import json, collections
F='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/molecular_object_sidecar.jsonl'
c=collections.defaultdict(collections.Counter)
for line in open(F):
    r=json.loads(line)
    aa=json.loads(r['assay_applicability'])
    for k in ['bulk_rna','chromatin','proteomics','spatial']: c['applic_'+k][aa.get(k)]+=1
    c['noncoding_dna_context'][str(r['noncoding_dna_context'])[:60]]+=1
    c['credible_set_id_empty'][r['credible_set_id']=='']+=1
    c['cs_link'][r['credible_set_to_gene_link_status']]+=1
for k,v in c.items(): print(k,dict(v.most_common(6)))
