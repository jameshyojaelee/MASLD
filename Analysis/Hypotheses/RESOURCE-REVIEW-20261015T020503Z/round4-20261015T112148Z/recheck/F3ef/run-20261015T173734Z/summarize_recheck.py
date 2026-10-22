import pandas as pd, numpy as np
B = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
O = B + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F3ef/run-20261015T173734Z/'
S = B + 'figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/'
coh = ['GSE130970', 'GSE135251', 'GSE162694']
fit = lambda lab, c: pd.read_csv(f'{O}fits/{lab}__{c}.tsv.gz', sep='\t').set_index('gene')
out = []
def rec(k, v): out.append((k, v)); print(k, v)

c = pd.read_csv(O + 'my_counts.tsv', sep='\t'); print(c.to_string())
print(pd.read_csv(O + 'my_reproduction.tsv', sep='\t').to_string())
rel = pd.read_csv(S + 'stage_extension_all_gene_results.tsv', sep='\t', usecols=['gene_id_versioned', 'gene_id_base', 'gene_name', 'logFC', 'FDR', 'axis', 'contrast'])
f1 = rel[(rel.axis == 'fibrosis') & (rel.contrast == 'F1_vs_F0')].set_index('gene_id_versioned')
f4 = rel[(rel.axis == 'fibrosis') & (rel.contrast == 'F4_vs_F0')].set_index('gene_id_versioned')
calls1 = f1.index[f1.FDR < 0.05]; calls4 = f4.index[f4.FDR < 0.05]
w = fit('without_GSE130970', 'F1_vs_F0')
s14 = w.index[w.q < 0.05]; rec('F1 without130970 calls', len(s14)); rec('of those in 1338', len(set(s14) & set(calls1)))
for cc in coh[1:]:
    x = fit('without_' + cc, 'F1_vs_F0'); rec('F1 without ' + cc + ' retained of 1338', int((x.loc[calls1, 'q'] < 0.05).sum()))
w4 = fit('without_GSE130970', 'F4_vs_F0'); rec('F4 without130970 retained of 7043', int((w4.loc[calls4, 'q'] < 0.05).sum()))
# single-cohort same sign
al = pd.concat({cc: fit('only_' + cc, 'F0_to_F1').logFC for cc in coh}, axis=1)
ss = (np.sign(al.loc[calls1]).eq(np.sign(f1.loc[calls1, 'logFC']), axis=0)).all(axis=1).mean(); rec('F1 same sign single-cohort fits', round(ss, 4))
it = pd.read_csv(O + 'interaction_F1_F4.tsv.gz', sep='\t')
for con, calls, ref in [('F1_vs_F0', calls1, f1), ('F4_vs_F0', calls4, f4)]:
    x = it[it.contrast == con].set_index('gene')
    rec(con + ' interaction q<0.05 all genes', int((x.F_q < 0.05).sum()))
    rec(con + ' interaction q<0.05 among calls', int((x.loc[calls, 'F_q'] < 0.05).sum()))
    sg = np.sign(x.loc[calls, coh]).eq(np.sign(ref.loc[calls, 'logFC']), axis=0)
    rec(con + ' same sign all three', round(sg.all(axis=1).mean(), 4))
    for cc in coh: rec(con + ' same sign ' + cc, round(sg[cc].mean(), 4))
# Fig 3F marks under omission
g = pd.read_csv(B + 'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv')
pos = set(g.loc[g.coloc_best_susie_pp4 > 0.5, 'ensembl'])
tr = ['F0_to_F1', 'F1_to_F2', 'F2_to_F3', 'F3_to_F4']
adj = rel[(rel.axis == 'fibrosis_adjacent') & rel.contrast.isin(tr) & rel.gene_id_base.isin(pos) & (rel.FDR < 0.05)]
rec('marks', len(adj))
keepall = None
for cc in coh:
    ks = {}
    kept_ids = set()
    for t in tr:
        m = adj[adj.contrast == t]
        x = fit('without_' + cc, t).loc[m.gene_id_versioned]
        k = (x.q.values < 0.05) & (np.sign(x.logFC.values) == np.sign(m.logFC.values))
        ks[t] = int(k.sum()); kept_ids |= set(zip(m.gene_name[k], [t] * k.sum()))
        if cc == 'GSE130970' and t == 'F0_to_F1': rec('min q F0F1 marks without130970', round(x.q.min(), 4))
    rec('marks kept without ' + cc, ks)
    keepall = kept_ids if keepall is None else keepall & kept_ids
rec('marks kept under all three omissions', sorted(keepall))
# hero genes
hero = rel[rel.gene_name.isin(['GNMT', 'MAT1A', 'CYP2C19'])][['gene_id_versioned', 'gene_name']].drop_duplicates()
rows = []
for lab in ['all3'] + ['without_' + cc for cc in coh]:
    for t in ['F0_to_F1', 'F1_to_F2', 'F2_vs_F0', 'F3_vs_F0', 'F4_vs_F0']:
        x = fit(lab, t)
        for gid, gn in hero.values:
            rows.append((lab, t, gn, x.loc[gid, 'logFC'], x.loc[gid, 'q']))
h = pd.DataFrame(rows, columns=['fit', 'contrast', 'gene', 'logFC', 'q'])
h.to_csv(O + 'hero_recheck.tsv', sep='\t', index=False)
print(h.pivot_table(index=['gene', 'contrast'], columns='fit', values='q').to_string())
late = h[h.contrast.isin(['F2_vs_F0', 'F3_vs_F0', 'F4_vs_F0'])]
print(late.groupby('fit').agg(maxq=('q', 'max'), maxlfc=('logFC', 'max')))
print(late[late.fit != 'all3'].sort_values('q', ascending=False).head(4).to_string())
pd.DataFrame(out, columns=['check', 'value']).to_csv(O + 'recheck_values.tsv', sep='\t', index=False)
