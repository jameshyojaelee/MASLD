import pandas as pd
W='/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig6/work/'
d=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv')
d['ens']=d.gene.str.split('.').str[0]
m=pd.read_csv(W+'catalog_gene_calls_with_coloc.tsv',sep='\t')
u=m[m.rna_call=='untestable']
hit=u[u.symbol.isin(set(d.symbol.dropna()))]
print('RNA untestable',len(u),' symbol present in 774 fit under another ID',len(hit))
x=hit.merge(d[['symbol','ens','logFC','padj']],on='symbol')
x['sup']=(x.padj<0.05)&(x.logFC.abs()>0.5)
print(x[['symbol','ensembl_id','ens','logFC','padj','sup','call_state','primary_evidence_class']].to_string())
