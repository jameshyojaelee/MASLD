import json, sys
import numpy as np, pandas as pd
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CAT=ROOT+"/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8"
COLOC=ROOT+"/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
BULK=ROOT+"/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
OUT=sys.argv[1]; R={}
b=lambda s: s.astype(str).str.split(".").str[0]
e=pd.read_parquet(CAT+"/passport_evidence_long.parquet"); g=pd.read_parquet(CAT+"/passport_gene_index.parquet")
gen=e[e.evidence_domain=="genetics"].set_index("ensembl_id")
R["gen_direction"]={str(k):int(v) for k,v in gen.direction.value_counts(dropna=False).items()}
R["rna_direction"]={str(k):int(v) for k,v in e[e.evidence_domain=="transcriptomics"].direction.value_counts(dropna=False).items()}
tier=pd.read_csv(ROOT+"/GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t").set_index("study_name").tier
a=pd.read_csv(COLOC+"/susie_coloc_all_gwas.csv",usecols=["gwas_name","ensembl","PP.H3.susie","PP.H4.susie","method","n_cs_pairs_eligible"])
a["ens"]=b(a.ensembl); a=a[a.gwas_name.map(tier).isin([1,2])]
s=a.dropna(subset=["PP.H4.susie"])
gg=s.groupby("ens").agg(h3=("PP.H3.susie","max"),h4=("PP.H4.susie","max"))
R["allg_t12_genes_with_susie"]=len(gg)
hit=gg[(gg.h3>0.8)&(gg.h4<=0.5)]
R["allg_t12_h3gt08_h4le05"]=len(hit)
R["allg_t12_hits_gen_call"]={str(k):int(v) for k,v in pd.Series([gen.call_state.get(x,"absent") for x in hit.index]).value_counts().items()}
sp=pd.read_csv(COLOC+"/susie_coloc_signal_pairs.tsv.gz",sep="\t"); sp["ens"]=b(sp.ensembl)
el=sp[(sp.eligible.astype(str).str.upper()=="TRUE")&sp.gwas.map(tier).isin([1,2])]
R["signal_pair_genes_t12_eligible"]=int(el.ens.nunique())
R["signal_pair_genes_t12_any"]=int(sp[sp.gwas.map(tier).isin([1,2])].ens.nunique())
R["allg_hits_not_in_signalpairs"]=sorted(set(hit.index)-set(el.ens))[:50]
R["n_allg_hits_not_in_signalpairs"]=len(set(hit.index)-set(el.ens))
R["n_eligible_col_summary"]={str(k):int(v) for k,v in s.n_cs_pairs_eligible.fillna(-1).clip(upper=3).value_counts().items()}
# union check
c=pd.read_csv(BULK+"/coefficient_table.tsv.gz",sep="\t",usecols=["gene"]); t12=pd.read_csv(COLOC+"/gene_level_coloc_tier12.csv",usecols=["ensembl"])
U=set(b(c.gene))|set(b(t12.ensembl)); C=set(g.ensembl_id)
R["union"]=len(U); R["catalog_minus_union"]=len(C-U); R["union_minus_catalog"]=len(U-C)
open(OUT+"/recheck_extra.json","w").write(json.dumps(R,indent=1))
print("extra done")
