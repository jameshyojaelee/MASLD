import pandas as pd
W='/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig6/work/'
d=pd.read_parquet('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.parquet')
print(d.shape)
cols=[c for c in d.columns if 'coding' in c or 'class' in c or 'call' in c or 'object' in c or 'route' in c or 'rule' in c]
print(cols)
m=pd.read_csv(W+'catalog_gene_calls_with_coloc.tsv',sep='\t')
x=m.merge(d,on='ensembl_id',suffixes=('','_r9d'))
print('merged',len(x))
pc=[c for c in d.columns if c=='primary_evidence_class']
print('class agree', (x.primary_evidence_class==x['primary_evidence_class_r9d']).mean() if 'primary_evidence_class_r9d' in x else 'n/a')
cs=x.coding_cs_state.value_counts(); print(cs)
pos=x[x.coding_cs_state.astype(str).str.contains('ge_0.5')]
print('coding PIP>=0.5 genes',len(pos))
print(pd.crosstab([pos.call_state,pos.rna_call],'n'))
print(pos.groupby(['call_state','rna_call']).symbol.apply(lambda s:';'.join(sorted(s))).to_string())
for s in ['PNPLA3','TM6SF2','MTARC1','GCKR','HSD17B13']:
    r=x[x.symbol==s]
    print(s, r[['primary_evidence_class','call_state','rna_call','coding_cs_state','coding_cs_max_pip','coding_cs_n_studies','coloc_best_susie_pp4','coloc_best_abf_pp4']].to_dict('records'))
if 'molecular_object' in d.columns or any('object' in c for c in d.columns):
    oc=[c for c in d.columns if 'object' in c]; print(oc)
    for c in oc[:4]:
        print(d[c].value_counts().head(10))
