"""Stage-eligible donors per lineage at >=50 cells (corrected 2026-10-14 roster)."""
import pandas as pd
R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
D = R + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L03/'
r = pd.read_csv(R + 'Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/'
                'metadata/current_source_diagnosis_donor_roster.tsv', sep='\t')
c = pd.read_csv(D + 'c3_donor_lineage_counts.tsv', sep='\t')
c['donor'] = c['donor'].str.replace(r'^(GSE174748|GSE189600)_', '', regex=True)
r = r[~r.exclude_stage_analysis.astype(str).str.upper().eq('TRUE')]
m = r.merge(c, on=['donor', 'dataset'], how='left', validate='one_to_one')
print('stage-eligible donors', len(r), 'matched', int(m['n_Hepatocytes'].notna().sum()))
order = ['Healthy', 'Steatosis', 'Steatohepatitis', 'Cirrhosis']
for L in ['Hepatocytes', 'Endothelial_cells', 'Fibroblasts', 'Cholangiocytes', 'Macrophages', 'T_cells']:
    k = m[m['n_' + L] >= 50]
    print(L, len(k), dict(zip(order, k.disease_stage_coarse.value_counts().reindex(order).fillna(0).astype(int))),
          dict(k.dataset.value_counts()))
