"""Per-lineage donor n by dataset x stage at >=50 lineage cells (gate-L03-03)."""
import pandas as pd
R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
D = R + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L03/'
r = pd.read_csv(R + 'Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/'
                'metadata/current_source_diagnosis_donor_roster.tsv', sep='\t')
c = pd.read_csv(D + 'c3_donor_lineage_counts.tsv', sep='\t')
print('lens count table datasets:', c.dataset.value_counts().to_dict())
print('example donors:', c.donor.head(3).tolist(), c[c.dataset.isin(['GSE174748','GSE189600'])].donor.tolist())
c['donor'] = c['donor'].str.replace(r'^(GSE174748|GSE189600)_', '', regex=True)
r = r[~r.exclude_stage_analysis.astype(str).str.upper().eq('TRUE')]
m = r.merge(c, on=['donor', 'dataset'], how='left', validate='one_to_one')
order = ['Healthy', 'Steatosis', 'Steatohepatitis', 'Cirrhosis']
for L in ['Hepatocytes', 'Endothelial_cells', 'Fibroblasts', 'Cholangiocytes']:
    k = m[m['n_' + L] >= 50]
    t = pd.crosstab(k.dataset, k.disease_stage_coarse).reindex(columns=order, fill_value=0)
    t['total'] = t.sum(axis=1)
    # donors in datasets with >=2 stage levels (inform the within-dataset slope)
    inf = t[(t[order] > 0).sum(axis=1) >= 2]['total'].sum()
    print('\n==', L, 'n =', len(k), '| donors in datasets with >=2 stages =', inf)
    print(t.to_string())
    print('median lineage cells', k['n_' + L].median(), 'min', k['n_' + L].min())
