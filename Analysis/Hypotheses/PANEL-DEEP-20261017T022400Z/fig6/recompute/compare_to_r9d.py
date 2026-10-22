"""After recomputation: compare per-gene recomputed calls with the r9d class columns."""
import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig6/recompute/"
c=pd.read_csv(R+"catalog_calls_recomputed.tsv.gz",sep="\t")
d=pd.read_parquet("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.parquet",
  columns=["ensembl_id","r8_genetics_call_state","r8_rna_call_state","primary_evidence_class","untested_missing_domain","gen_coloc_best_susie_pp4","gen_coloc_best_abf_pp4","rna_log2fc","rna_bh_q","next_experiment_rule_id","coding_cs_state","coding_cs_max_pip","symbol","gen_untestable_reason_r9"])
x=c.merge(d,on="ensembl_id",suffixes=("","_r9d"))
print("rows",len(x))
print("gen call disagree", (x.gen_call!=x.r8_genetics_call_state).sum())
print("rna call disagree", (x.rna_call!=x.r8_rna_call_state).sum())
print("class disagree", (x.class_recomputed!=x.primary_evidence_class).sum())
print(pd.crosstab(x.untested_domain.fillna(""),x.untested_missing_domain.fillna("")))
import numpy as np
print("susie max absdiff", np.nanmax(np.abs(x.susie_best-x.gen_coloc_best_susie_pp4)))
print("abf max absdiff", np.nanmax(np.abs(x.abf_best-x.gen_coloc_best_abf_pp4)))
print("log2fc max absdiff", np.nanmax(np.abs(x.rna_log2fc-x.rna_log2fc_r9d)))
print("q max absdiff", np.nanmax(np.abs(x.rna_q-x.rna_bh_q)))
print(pd.crosstab(x.primary_evidence_class,x.next_experiment_rule_id))
# coding fine-mapped genes by matrix cell
cod=x[x.coding_cs_state=="fine_mapped_coding_pip_ge_0.5"]
print("coding pip>=0.5 genes",len(cod)); print(pd.crosstab(cod.gen_call,cod.rna_call))
print(cod.sort_values("coding_cs_max_pip",ascending=False)[["symbol","coding_cs_max_pip","gen_call","rna_call","class_recomputed"]].to_string())
print(x.gen_untestable_reason_r9.value_counts(dropna=False))
