"""Independent per-donor lineage cell counts from run-level metadata + donor_pairing,
and presence of each donor's runs in the March pseudobulk CSVs (gate-L03-03)."""
import pandas as pd
R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
md = pd.read_csv(R + 'Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/'
                 'metadata/donor_metadata_stage_corrected.tsv', sep='\t')
ros = pd.read_csv(R + 'Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/'
                  'metadata/current_source_diagnosis_donor_roster.tsv', sep='\t')
ros = ros[~ros.exclude_stage_analysis.astype(str).str.upper().eq('TRUE')]
m = {}
for ds in ['GSE244832', 'GSE202379', 'GSE185477']:
    p = pd.read_csv(R + f'data/{ds}/metadata/donor_pairing.csv', dtype=str)
    for _, row in p.iterrows():
        for s in row.rna_srrs.split(';'):
            m[s] = f'{ds}_{row.donor_id}'
md['donor'] = [m.get(s, s) for s in md['sample']]
print('roster donor examples', ros.donor.head(3).tolist(), ros[ros.dataset.isin(['GSE174748','GSE189600'])].donor.tolist())
md = md[md.dataset.isin(ros.dataset.unique())]
Ls = ['Hepatocytes', 'Endothelial_cells', 'Fibroblasts', 'Cholangiocytes']
agg = md.groupby(['donor', 'dataset'])[[f'n_{L}' for L in Ls] + ['n_cells']].sum().reset_index()
x = ros.merge(agg, on=['donor', 'dataset'], how='left', validate='one_to_one')
print('roster donors', len(x), 'matched to run metadata', x.n_cells.notna().sum())
order = ['Healthy', 'Steatosis', 'Steatohepatitis', 'Cirrhosis']
lens = pd.read_csv(R + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L03/c3_donor_lineage_counts.tsv', sep='\t')
lens['donor'] = lens['donor'].str.replace(r'^(GSE174748|GSE189600)_', '', regex=True)
y = x.merge(lens, on=['donor', 'dataset'], suffixes=('', '_lens'))
for L in Ls:
    k = x[x[f'n_{L}'] >= 50]
    diff = (y[f'n_{L}'] - y[f'n_{L}_lens']).abs().max()
    print(L, 'n>=50:', len(k), k.disease_stage_coarse.value_counts().reindex(order).fillna(0).astype(int).tolist(),
          '| max abs diff vs lens table:', diff)
# presence in pseudobulk CSVs: does every run of each eligible donor appear?
inv = {}
for s, d in m.items():
    inv.setdefault(d, []).append(s)
for L in Ls:
    cols = set(pd.read_csv(R + f'Analysis/SingleCell/results_gpu_v2/pseudobulk/{L}_pseudobulk.csv', nrows=0, index_col=0).columns)
    k = x[x[f'n_{L}'] >= 50]
    miss_donor, miss_runs, tot_runs = 0, 0, 0
    for d in k.donor:
        runs = inv.get(d, [d])
        tot_runs += len(runs)
        mr = [r for r in runs if r not in cols]
        miss_runs += len(mr)
        if len(mr) == len(runs):
            miss_donor += 1
    # runs missing from CSV that have lineage cells
    rr = md[md.donor.isin(k.donor)]
    absent_cells = rr[~rr['sample'].isin(cols)][f'n_{L}'].sum()
    print(L, 'CSV cols', len(cols), '| eligible donors fully absent', miss_donor, '| runs absent', miss_runs, 'of', tot_runs,
          '| lineage cells in absent runs', int(absent_cells), 'of', int(rr[f'n_{L}'].sum()))
