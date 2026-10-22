# Gate C-02, step 1 (login node, detection only, no spot scores).
# Independently re-derive the primary pair set and check the fixed-effect design.
# Uses no spatial statistic and no outcome.
import h5py, numpy as np, pandas as pd, scipy.sparse as sp

B = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT = B + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-C-02/'
T4 = B + 'RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction/'
mem = pd.read_csv(B + 'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv', sep='\t')
fig4d = pd.read_csv(B + 'figures/main/fig4_singlecell_programs/source_tables/current_candidate/fig4_compact_broad_20261005/fig4d_all117_lineage_content.tsv', sep='\t')
dm = pd.read_csv(B + 'Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv', sep='\t')
ed = pd.read_csv(T4 + 'program_interaction_edges.tsv.gz', sep='\t')

print('edges rows', len(ed), '| programs in T4', len(set(ed.program_a) | set(ed.program_b)))
print('fig4d rows', len(fig4d), '| label_matches_content True', int(fig4d.label_matches_content.sum()))
matched = set(fig4d.program_uid[fig4d.label_matches_content])

ed = ed[ed.support != 'membership_entangled'].copy()
ed['cls'] = np.select([ed.support.str.startswith('coordinated_beyond'), ed.support == 'coordinated_global_component_only',
                       ed.support == 'informative_null'], ['coordinated', 'global_only', 'informative_null'], 'other')
ed['cross'] = ed.cell_type_a != ed.cell_type_b
ed['both_matched'] = ed.program_a.isin(matched) & ed.program_b.isin(matched)
print('T4 lineages', sorted(set(ed.cell_type_a) | set(ed.cell_type_b)))
c = ed[ed.cross & ed.both_matched]
print('cross-lineage, both lineage-matched (before coverage):', c.cls.value_counts().to_dict())

# membership weights; all members in denominator (unmapped kept as their own keys); duplicate symbols summed
ok = mem.mapped_symbol_status.isin(['gencode_v49_unique_symbol_confirmed', 'gencode_v49_unambiguous_ensembl_to_symbol'])
mem['key'] = np.where(ok, mem.mapped_symbol, 'UNMAPPED_' + mem.source_gene.astype(str))
w = mem.groupby(['program_uid', 'key'], as_index=False).original_l1_weight.sum()
allw = {p: dict(zip(g.key, g.original_l1_weight)) for p, g in w.groupby('program_uid')}
want = set(mem.mapped_symbol[ok])


def read_sym(g):
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g['categories'][:]])
        return cats[g['codes'][:]]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in g[:]])


prim = list(dm.array_id[dm.group == 'primary'])
pairs = c[['program_a', 'program_b', 'cell_type_a', 'cell_type_b', 'cls', 'r_full', 'r_disjoint', 'r_partial']].reset_index(drop=True)
n_ok50 = np.zeros(len(pairs), int)
n_ok80 = np.zeros(len(pairs), int)
for arr in prim:
    f = h5py.File(B + f'Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled/{arr}.h5ad', 'r')
    sym = read_sym(f['var']['gene_symbol']); X = f['X']
    M = sp.csr_matrix((X['data'][:], X['indices'][:], X['indptr'][:]), shape=tuple(X.attrs['shape'])); f.close()
    keep = np.isin(sym, list(want)); u, inv = np.unique(sym[keep], return_inverse=True)
    agg = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)), shape=(len(inv), len(u)))
    frac = np.asarray(((M[:, np.where(keep)[0]] @ agg) > 0).mean(axis=0)).ravel()
    det = set(u[frac >= 0.01])
    r50, r80 = [], []
    for a, b in zip(pairs.program_a, pairs.program_b):
        ga, gb = allw[a], allw[b]; sh = set(ga) & set(gb)
        da = {k: v for k, v in ga.items() if k not in sh}; db = {k: v for k, v in gb.items() if k not in sh}
        ca = sum(v for k, v in da.items() if k in det) / sum(da.values())
        cb = sum(v for k, v in db.items() if k in det) / sum(db.values())
        r50.append(ca >= 0.5 and cb >= 0.5); r80.append(ca >= 0.8 and cb >= 0.8)
    n_ok50 += np.array(r50); n_ok80 += np.array(r80)
pairs['n_arrays_cov50'] = n_ok50
pairs['n_arrays_cov80'] = n_ok80
pairs.to_csv(OUT + 'g1_cross_matched_pairs_coverage.tsv', sep='\t', index=False)

P = pairs[pairs.n_arrays_cov50 == 32].copy()
print('PRIMARY SET (cov50 disjoint, all 32):', P.cls.value_counts().to_dict(), '| programs', len(set(P.program_a) | set(P.program_b)))
P80 = pairs[pairs.n_arrays_cov80 == 32]
print('cov80 disjoint all 32:', P80.cls.value_counts().to_dict())
P['lp'] = ['-'.join(sorted([a, b])) for a, b in zip(P.cell_type_a, P.cell_type_b)]
print(pd.crosstab(P.lp, P.cls).to_string())

# sign of the bulk raw correlation within each class
for k in ['coordinated', 'informative_null', 'global_only']:
    s = P[P.cls == k]
    print(f'{k}: n {len(s)} | r_disjoint<0 {int((s.r_disjoint < 0).sum())} | r_full<0 {int((s.r_full < 0).sum())} | '
          f'median r_disjoint {s.r_disjoint.median():.3f} | median r_partial {s.r_partial.median():.3f}')

# fixed-effect design for the primary contrast: strata = (program a x lineage b) and (program b x lineage a)
Q = P[P.cls.isin(['coordinated', 'informative_null'])].reset_index(drop=True)
s1 = Q.program_a + '|' + Q.cell_type_b
s2 = Q.program_b + '|' + Q.cell_type_a
strata = sorted(set(s1) | set(s2))
idx = {s: i for i, s in enumerate(strata)}
D = np.zeros((len(Q), len(strata)))
for r, (x, y) in enumerate(zip(s1, s2)):
    D[r, idx[x]] += 1; D[r, idx[y]] += 1
y = (Q.cls == 'coordinated').to_numpy(float)
lab = pd.DataFrame({'s': list(s1) + list(s2), 'y': list(y) + list(y)})
both = lab.groupby('s').y.agg(lambda v: (v.min() == 0) and (v.max() == 1))
rankD = np.linalg.matrix_rank(D)
H = D @ np.linalg.pinv(D)
yt = y - H @ y
print(f'primary pairs {len(Q)} | strata {len(strata)} | strata with both classes {int(both.sum())} | rank(FE) {rankD} | residual df {len(Q) - rankD - 1}')
print(f'var(coordinated) raw {y.var():.4f} | after FE {yt.var():.4f} | share of label variance left {yt.var() / y.var():.3f}')
print(f'effective n for beta, (y~)^T(y~) = {float(yt @ yt):.1f}  vs no-FE n0*n1/n = {y.sum() * (1 - y).sum() / len(y):.1f}')
# leverage of the most influential pairs on beta: weight of each pair = y~_i / sum(y~^2)
wts = np.abs(yt) / (yt @ yt)
print('pairs carrying 50% of |weight| on beta:', int(np.searchsorted(np.cumsum(np.sort(wts)[::-1]) / wts.sum(), 0.5) + 1))
np.save(OUT + 'g1_fe_design_strata.npy', D)
Q.to_csv(OUT + 'g1_primary_pairs.tsv', sep='\t', index=False)
