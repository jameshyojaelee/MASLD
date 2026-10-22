# H86 audit: recompute Catalog, DE-table, deposit and sample counts; direct-join reproduction (known before writing).
import pandas as pd, numpy as np
C='figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
d=pd.read_parquet(C+'passport_gene_index.parquet')
print(d.primary_evidence_class.value_counts().to_dict(), len(d), d.symbol.nunique(), d.ensembl_id.nunique(), int(d.symbol_collision.sum()))
p=pd.read_csv('Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv',sep='\t',keep_default_na=False,na_values=[''])
print('DE rows',len(p),'unique',p.gene.nunique(),'sep',int(p.gene.str.contains('[;,|]').sum()),'sig',int((p.padj<0.05).sum()))
ps=set(p.gene); sig=set(p.gene[p.padj<0.05]); cs=set(d.symbol)
d['q']=d.symbol.isin(ps); d['s']=d.symbol.isin(sig)
print(d.groupby('primary_evidence_class')[['q','s']].sum().to_dict())
print('sig not in Catalog', sorted(sig-cs))
q=pd.read_csv(C+'passport_identity_quarantine.tsv',sep='\t',keep_default_na=False)
qs=set(q.source_symbol)-{''}
ov=ps&qs&cs
print('direct protein symbols also quarantine source_symbol', len(ov), q[q.source_symbol.isin(ov)].reason_code.value_counts().to_dict())
print('APOB in overlap', 'APOB' in ov)
print('NC-2 ceiling pct', round(100*(409+90)/13551,2))
