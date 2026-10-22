"""Count molecular-object types and routes (Fig 6A box, 6C object-route column)."""
import json, collections, pandas as pd
M="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/molecular_object_sidecar.jsonl"
rows=[json.loads(l) for l in open(M)]
d=pd.DataFrame(rows)
print(len(d))
print(d.candidate_object.value_counts())
print(d.next_experiment_rule_id.value_counts())
print(pd.crosstab(d.candidate_object,d.next_experiment_rule_id))
print(pd.crosstab(d.primary_evidence_class,d.legacy_next_experiment_rule_id))
print(d.credible_set_to_gene_link_status.value_counts())
print(d.noncoding_dna_context.value_counts())
aa=d.assay_applicability.map(json.loads).apply(pd.Series)
for c in aa: print(c, aa[c].value_counts().to_dict())
w=d[d.symbol.isin(["GLP1R","HKDC1","MTARC1","THRB","PNPLA3","TM6SF2"])][["symbol","candidate_object","primary_evidence_class","legacy_next_experiment_rule_id","next_experiment_rule_id","genetics_call","rna_call","gene_biotype"]]
print(w.to_string())
out="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig6/recompute/"
pd.crosstab(d.candidate_object,d.next_experiment_rule_id,margins=True).to_csv(out+"fig6a_object_type_by_route.tsv",sep="\t")
pd.crosstab(d.primary_evidence_class,d.legacy_next_experiment_rule_id,margins=True).to_csv(out+"fig6d_class_by_gene_route.tsv",sep="\t")
w.to_csv(out+"fig6c_object_routes.tsv",sep="\t",index=False)
