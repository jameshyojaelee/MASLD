# Debug check (not part of the analysis): why row-level TSS mapping differs between gene_id keying (gate, 0.99721)
# and (gene_id, top_snp chromosome) keying (H74). Reads gene records only; no anchor, no nearest gene.
import re, pandas as pd
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
tier=pd.read_csv(R+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t'); t12=set(tier.study_name[tier.tier.isin([1,2])])
d=pd.read_csv(R+'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv',usecols=['gwas_name','gene','ensembl','chr','PP.H4.susie','top_snp'])
s=d[d.gwas_name.isin(t12)&d['PP.H4.susie'].notna()].copy(); s['snp_chr']=s.top_snp.str.split(':').str[0]
gid={};gc=set()
for l in open(R+'data/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf'):
    if l[0]=='#' or '\tgene\t' not in l: continue
    f=l.split('\t')
    if f[2]!='gene': continue
    full=re.search('gene_id "([^"]+)"',f[8]).group(1); g=re.search('gene_id "([^".]+)',f[8]).group(1); c=f[0].replace('chr','')
    gid.setdefault(g,[]).append((c, full.endswith('_PAR_Y')))
    if not full.endswith('_PAR_Y'): gc.add((g,c))
s['by_gid']=s.ensembl.isin(gid); s['by_gc']=[(e,c) in gc for e,c in zip(s.ensembl,s.snp_chr)]
print('rows',len(s),'mapped by gid',s.by_gid.mean(),'by (gid,snp_chr)',s.by_gc.mean())
x=s[s.by_gid&~s.by_gc]
print(x[['gwas_name','gene','ensembl','chr','top_snp']].to_string()); print({e:gid[e] for e in x.ensembl.unique()})
