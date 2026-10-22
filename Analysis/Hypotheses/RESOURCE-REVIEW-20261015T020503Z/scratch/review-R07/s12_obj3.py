import json, collections
F='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/molecular_object_sidecar.jsonl'
C=collections.defaultdict(collections.Counter)
for line in open(F):
    d=json.loads(line)
    for k in ['assay_applicability','noncoding_dna_context','credible_set_to_gene_link_status','target_link_basis','strand_audit_status','lncrna_genomic_class','mapping_status','gene_biotype','input_state']:
        C[k][str(d[k])]+=1
    C['credible_set_id_nonempty'][bool(d['credible_set_id'])]+=1
    C['coding_pip_mass_nonnull'][d['credible_set_coding_pip_mass'] is not None]+=1
for k,v in C.items(): print(k, dict(v.most_common(8)))
