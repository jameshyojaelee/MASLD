#!/bin/bash
# gate-C-03: can GPT serve as a direction positive control? Extract GPT liver eQTL rows (Broadaway chr8, hg19)
# and GPT COLOC rows from COLOC r2. Read-only on repo files; outputs in this job dir.
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3931
export TMPDIR=$OUT/tmp; mkdir -p $TMPDIR
awk -F'\t' 'NR==1 || $14=="ENSG00000167701" || $13=="GPT"' $ROOT/data/broadaway_eqtl/chr8_marginal_summary_results.tsv > $OUT/gpt_broadaway_chr8.tsv
head -1 $ROOT/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv > $OUT/gpt_coloc_r2.csv
grep -F ',ENSG00000167701,' $ROOT/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv >> $OUT/gpt_coloc_r2.csv || true
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python - <<'PY'
import pandas as pd, numpy as np
O='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3931/'
b=pd.read_csv(O+'gpt_broadaway_chr8.tsv',sep='\t')
print('GPT Broadaway rows',len(b),'genes',b.ENSG.unique().tolist())
b=b.sort_values('PVAL')
print('rows p<1e-3',int((b.PVAL<1e-3).sum()),'p<5e-8',int((b.PVAL<5e-8).sum()))
print('top 5 eQTL variants'); print(b.head(5)[['Variant','POS','NEA','EA','EAF','Beta','SE','PVAL','N']].to_string(index=False))
h=b[b.POS==145733794]
print('row at ALT lead 8:145733794 (H10 GPT lead, ALT-raising allele C, other CG)'); print(h[['Variant','NEA','EA','EAF','Beta','SE','PVAL','N']].to_string(index=False))
c=pd.read_csv(O+'gpt_coloc_r2.csv')
print('COLOC r2 GPT rows',len(c))
if len(c):
    cols=[x for x in c.columns if x.lower() in ('gwas','pp.h4','pp.h3','pp.h4.abf','pp.h3.abf','susie_pp4','susie_pp.h4','method','nsnps')]
    print('columns',list(c.columns)[:40])
    print(c[['gwas_name','PP.H3.abf','PP.H4.abf','PP.H3.susie','PP.H4.susie','n_cs_pairs','method','top_snp','ancestry']].sort_values('PP.H3.abf',ascending=False).head(25).to_string(index=False))
PY
