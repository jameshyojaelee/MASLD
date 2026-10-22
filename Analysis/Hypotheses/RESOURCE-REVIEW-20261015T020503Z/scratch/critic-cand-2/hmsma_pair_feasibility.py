# Feasibility count for a donor-level spatial co-patterning test of bulk-coordinated program pairs in HMSMA.
# Reads only per-array detection. Computes no program score and no spatial statistic.
# For each array: fraction of each program's original L1 weight on member symbols detected in >=1% of spots
# (Ensembl IDs summed to symbols, as in V4). Then counts testable pairs per class and lineage pair, on full and on
# disjoint membership (shared genes removed from both programs, as in the bulk r_disjoint).
import sys
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

B = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT = B + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/critic-cand-2/'
T4 = B + 'RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction/'
mem = pd.read_csv(B + 'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv', sep='\t')
dmap = pd.read_csv(B + 'Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv', sep='\t')
pairs = pd.read_csv(T4 + 'program_interaction_edges.tsv.gz', sep='\t')

mem = mem[mem.mapped_symbol_status.str.contains('confirmed', na=False)].copy()
tot = mem.groupby('program_uid').original_l1_weight.sum()  # denominator: confirmed members (V4 convention)
prog_genes = {p: dict(zip(g.mapped_symbol, g.original_l1_weight)) for p, g in mem.groupby('program_uid')}
want = set(mem.mapped_symbol)


def read_sym(g):
    if isinstance(g, h5py.Group):
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g['categories'][:]])
        return cats[g['codes'][:]]
    return np.array([s.decode() if isinstance(s, bytes) else s for s in g[:]])


def detected_symbols(path):
    f = h5py.File(path, 'r')
    sym = read_sym(f['var']['gene_symbol'])
    X = f['X']
    shape = tuple(X.attrs['shape']) if 'shape' in X.attrs else tuple(X.attrs['h5sparse_shape'])
    M = sp.csr_matrix((X['data'][:], X['indices'][:], X['indptr'][:]), shape=shape)
    f.close()
    keep = np.isin(sym, list(want))
    sub = M[:, np.where(keep)[0]].tocsc()
    s_keep = sym[keep]
    uniq, inv = np.unique(s_keep, return_inverse=True)
    agg = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)), shape=(len(inv), len(uniq)))
    summed = sub @ agg
    det = np.asarray((summed > 0).mean(axis=0)).ravel()
    return set(uniq[det >= 0.01]), M.shape[0]


def coverage(genes, det):
    return sum(w for g, w in genes.items() if g in det)


pairs = pairs[pairs.support != 'membership_entangled'].copy()
pairs['cls'] = np.where(pairs.support.str.startswith('coordinated_beyond'), 'coordinated',
                        np.where(pairs.support == 'coordinated_global_component_only', 'global_only', 'informative_null'))
pairs['lp'] = [ '-'.join(sorted([a, b])) for a, b in zip(pairs.cell_type_a, pairs.cell_type_b)]
pairs['cross'] = pairs.cell_type_a != pairs.cell_type_b

rows = []
cov_rows = []
test_disjoint = {}
for arr in dmap.array_id:
    det, nspots = detected_symbols(B + f'Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled/{arr}.h5ad')
    cov = {p: coverage(g, det) / tot[p] for p, g in prog_genes.items()}
    for p, c in cov.items():
        cov_rows.append((arr, p, c))
    full_ok = pairs.program_a.map(cov).ge(0.8) & pairs.program_b.map(cov).ge(0.8)
    # disjoint membership: drop shared genes from both, coverage relative to the disjoint member weight
    dj = []
    for a, b in zip(pairs.program_a, pairs.program_b):
        ga, gb = prog_genes[a], prog_genes[b]
        sh = set(ga) & set(gb)
        if sh:
            da = {g: w for g, w in ga.items() if g not in sh}
            db = {g: w for g, w in gb.items() if g not in sh}
            ca = coverage(da, det) / max(sum(da.values()), 1e-12)
            cb = coverage(db, det) / max(sum(db.values()), 1e-12)
            dj.append(ca >= 0.8 and cb >= 0.8)
        else:
            dj.append(cov[a] >= 0.8 and cov[b] >= 0.8)
    dj = np.array(dj)
    test_disjoint[arr] = dj
    for view, ok in [('full', full_ok.values), ('disjoint', dj)]:
        t = pairs[ok]
        r = dict(array_id=arr, view=view, n_spots=nspots, n_detected_program_symbols=len(det),
                 n_programs_testable=int(sum(v >= 0.8 for v in cov.values())),
                 cross_coord=int(((t.cls == 'coordinated') & t.cross).sum()),
                 within_coord=int(((t.cls == 'coordinated') & ~t.cross).sum()),
                 cross_global_only=int(((t.cls == 'global_only') & t.cross).sum()),
                 cross_informative_null=int(((t.cls == 'informative_null') & t.cross).sum()),
                 within_noncoord=int(((t.cls != 'coordinated') & ~t.cross).sum()))
        tc = t[t.cross]
        g = tc.groupby('lp').cls.agg(lambda s: ((s == 'coordinated').sum() > 0) and ((s != 'coordinated').sum() > 0))
        r['lineage_pairs_with_both'] = int(g.sum())
        rows.append(r)
    print(arr, rows[-1], flush=True)

res = pd.DataFrame(rows).merge(dmap[['array_id', 'group']], on='array_id')
res.to_csv(OUT + 'hmsma_pair_feasibility_per_array.tsv', sep='\t', index=False)
pd.DataFrame(cov_rows, columns=['array_id', 'program_uid', 'weight_coverage']).to_csv(OUT + 'hmsma_program_coverage_per_array.tsv', sep='\t', index=False)

# common edge set: pair testable (disjoint view) in every primary array
prim = dmap.array_id[dmap.group == 'primary']
common = np.logical_and.reduce([test_disjoint[a] for a in prim])
pairs['testable_all_primary'] = common
pairs['n_primary_arrays_testable'] = np.sum([test_disjoint[a] for a in prim], axis=0)
pairs[['program_a', 'program_b', 'cell_type_a', 'cell_type_b', 'lp', 'cross', 'cls', 'r_partial', 'n_shared',
       'testable_all_primary', 'n_primary_arrays_testable']].to_csv(OUT + 'hmsma_pair_testability.tsv', sep='\t', index=False)
c = pairs[common]
print('COMMON SET (disjoint view, testable in all 32 primary arrays)')
print(pd.crosstab(c.lp, [c.cross, c.cls]))
print('totals', pd.crosstab(c.cross, c.cls))
for th in [32, 30, 28, 25]:
    m = pairs.n_primary_arrays_testable >= th
    print('>=', th, 'primary arrays: cross coordinated', int((m & pairs.cross & (pairs.cls == 'coordinated')).sum()),
          'cross noncoord', int((m & pairs.cross & (pairs.cls != 'coordinated')).sum()),
          'within coordinated', int((m & ~pairs.cross & (pairs.cls == 'coordinated')).sum()))
