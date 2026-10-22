"""Figure 4G node table on the 367 stage-complete disease units (109 GSE162694, 258 GSE213621).

Takes the installed 4G node overlay (fixed geometry, coordinates, memberships) and replaces
continuum_beta with the disease-only stage/sex-adjusted fixed-projection meta coefficient
from review job 7070 (slurm 22252072), which reproduces the installed 435-unit node betas
to 1.3e-15 on the as_published roster. No fitting here; read-only inputs; no randomness.
"""
import sys
import numpy as np
import pandas as pd

ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT = sys.argv[1].rstrip('/') + '/'
NODES = ROOT + 'figures/candidates/continuum-774-panels-20261011T152106Z/molecular_systems_atlas/overlays/all_node_stage_continuum_overlays.tsv.gz'
FE = ROOT + ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F4efgh/'
             'fig4h_disease_only/run_20261015b/feature_stage_continuum_betas_by_roster.tsv.gz')
SYS4H = ROOT + 'figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/system_stage_continuum_by_roster.tsv'

nodes = pd.read_csv(NODES, sep='\t')
fe = pd.read_csv(FE, sep='\t')
checks = []
out = {}
for roster in ['as_published', 'disease_only']:
    f = fe[fe.roster == roster][['feature_id', 'collection', 'stage', 'cont_fixed_projection']]
    x = nodes.merge(f, on=['feature_id', 'collection'], how='left', validate='one_to_one')
    main = x[x.component_role == 'main_system_component']
    unmatched = main[main.cont_fixed_projection.isna()]
    # unmatched nodes must be the nodes with no installed continuum beta
    checks.append(dict(roster=roster, check='unmatched_main_nodes', value=len(unmatched),
                       detail='installed continuum_beta NA for all unmatched: %s' % bool(unmatched.continuum_beta.isna().all())))
    d = (main.cont_fixed_projection - main.continuum_beta).abs()
    checks.append(dict(roster=roster, check='max_abs_diff_vs_installed_continuum_beta', value=float(d.max()), detail=''))
    checks.append(dict(roster=roster, check='median_continuum_beta_main_nodes', value=float(main.cont_fixed_projection.median()), detail=''))
    v = main.cont_fixed_projection.abs().dropna()
    checks.append(dict(roster=roster, check='q98_abs_continuum_beta (colour limit)', value=float(np.quantile(v, 0.98, method='median_unbiased')), detail='type 8'))
    x = x.copy()
    x['installed_continuum_beta_435'] = x['continuum_beta']
    x['continuum_beta'] = x['cont_fixed_projection']
    x['stage_beta'] = x['stage']
    x = x.drop(columns=['cont_fixed_projection', 'stage'])
    # q values and support flags belong to the 435-unit fits; blank them in the disease-only table
    if roster == 'disease_only':
        for c in ['stage_q', 'stage_supported', 'continuum_q', 'continuum_supported', 'support_class']:
            x[c] = np.nan
    out[roster] = x
    x.to_csv(OUT + f'fig4g_nodes_{roster}.tsv.gz', sep='\t', index=False)

# system-level consistency with candidate 4H (collection-balanced medians)
s = pd.read_csv(SYS4H, sep='\t')
for roster in ['as_published', 'disease_only']:
    x = out[roster]
    x = x[x.component_role == 'main_system_component']
    fam = x.dropna(subset=['continuum_beta']).groupby(['community_id', 'collection']).continuum_beta.median()
    S = fam.groupby('community_id').mean()
    ref = s[s.roster == roster].set_index('community_id').cont_fixed_projection
    checks.append(dict(roster=roster, check='system_beta_max_abs_diff_vs_4H_source', value=float((S - ref.reindex(S.index)).abs().max()), detail=''))
    checks.append(dict(roster=roster, check='S01_system_continuum_beta', value=float(S.loc['system_01']), detail=''))
pd.DataFrame(checks).to_csv(OUT + 'fig4g_node_checks.tsv', sep='\t', index=False)
print(pd.DataFrame(checks).to_string())
