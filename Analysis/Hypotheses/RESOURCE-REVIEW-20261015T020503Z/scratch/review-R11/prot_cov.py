import pandas as pd
C='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
d=pd.read_parquet(C)
print(d.primary_evidence_class.value_counts().to_dict())
p=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv',sep='\t',keep_default_na=False,na_values=[''])
print('protein rows',len(p), 'unique genes', p.gene.nunique(), 'padj<0.05', (p.padj<0.05).sum())
ps=set(p.gene); sig=set(p.gene[p.padj<0.05])
d['prot_measured']=d.symbol.isin(ps); d['prot_sig']=d.symbol.isin(sig)
t=d.groupby('primary_evidence_class').agg(n=('symbol','size'),measured=('prot_measured','sum'),sig=('prot_sig','sum'))
t['pct_measured']=(100*t.measured/t.n).round(1)
print(t)
print('symbol collisions among measured:', d[d.prot_measured].symbol.duplicated().sum())
