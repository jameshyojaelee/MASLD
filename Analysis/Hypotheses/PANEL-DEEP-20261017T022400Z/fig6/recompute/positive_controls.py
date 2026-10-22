"""Positive controls: coding fine-mapped variants for PNPLA3, TM6SF2, MTARC1 (and GCKR) from the H87 member map,
plus their eQTL-colocalization and RNA calls from the recomputation."""
import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
m=pd.read_csv(R+"Analysis/Hypotheses/H87-20261015T130324Z/run-20261015T131921Z/h87_coding_member_gene_map.tsv",sep="\t",low_memory=False)
c=pd.read_csv(R+"Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig6/recompute/catalog_calls_recomputed.tsv.gz",sep="\t")
out=[]
for g in ["PNPLA3","TM6SF2","MTARC1","GCKR","HSD17B13","MBOAT7","APOE"]:
    x=m[(m.symbol==g)&(m.canonical_PA==True)]
    r=c[c.symbol==g]
    row={"gene":g}
    if len(x):
        t=x.sort_values("normalized_susie_pip",ascending=False).iloc[0]
        row.update(max_pip=round(t.normalized_susie_pip,3),top_variant_hg19=t.variant_id,pos_grch38=t.pos_grch38,top_study=t.study,
                   consequence=t.gene_T_terms,n_studies_any=x.study.nunique(),n_studies_pip_ge_0_5=x[x.normalized_susie_pip>=0.5].study.nunique())
    if len(r):
        r=r.iloc[0]
        row.update(eqtl_coloc_call=r.gen_call,best_susie=r.susie_best,best_abf=round(r.abf_best,3) if pd.notna(r.abf_best) else None,abf_study=r.abf_study,
                   rna_call=r.rna_call,rna_log2fc=round(r.rna_log2fc,3) if pd.notna(r.rna_log2fc) else None,rna_q=r.rna_q,catalog_class=r.class_recomputed)
    out.append(row)
o=pd.DataFrame(out); print(o.to_string())
o.to_csv(R+"Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig6/recompute/positive_controls.tsv",sep="\t",index=False)
