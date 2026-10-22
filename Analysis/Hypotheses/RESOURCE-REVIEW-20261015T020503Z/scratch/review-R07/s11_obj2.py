import json
F='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/molecular_object_sidecar.jsonl'
for line in open(F):
    d=json.loads(line)
    if d['candidate_object']=='untestable_nomination' and d['primary_evidence_class']!='untested':
        print(d['symbol'],d['ensembl_id_base'],d['gene_biotype'],d['primary_evidence_class'],d['genetics_call'],d['rna_call'],d['mapping_status'],d['counting_annotation_disposition'],d['open_mechanistic_question'])
