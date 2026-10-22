# Gate C-02 feasibility job 3921.
# Computes, on the 32 primary HMSMA arrays:
#   (1) positive controls PC-A (pericentral x periportal hepatocyte panels) and PC-B (scar mesenchyme x scar macrophage panels),
#       Lee's L of residualized panel scores, per donor;
#   (2) label-free noise structure of pair-level Lee's L for the 662 primary pairs (after the program x partner-lineage
#       fixed effects), and the false-positive rate of the proposed donor sign-flip test when the coordinated / null labels
#       are randomly reassigned to pairs (labels carry no information by construction).
# The real coordinated-vs-null contrast is NOT computed here. Only the class sizes are read.
# Approximations vs the candidate's V4 engine: genes are z-scored within each array (V4 z-scores over the pooled 35 arrays).
import os
import time
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components

t0 = time.time()
B = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
G = B + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-C-02/'
OUT = os.environ.get('OUTDIR', B + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3921/')
rng = np.random.default_rng(3921)

mem = pd.read_csv(B + 'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv', sep='\t')
dm = pd.read_csv(B + 'Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv', sep='\t')
Q = pd.read_csv(G + 'g1_primary_pairs.tsv', sep='\t')
n1 = int((Q.cls == 'coordinated').sum()); n0 = int((Q.cls == 'informative_null').sum())
assert n1 + n0 == len(Q)
Q = Q[['program_a', 'program_b', 'cell_type_a', 'cell_type_b']]  # labels dropped from here on
mk = pd.read_csv(B + 'data/external/liver_celltype_markers_combined.tsv', sep='\t')
PROXY = ['ALB', 'APOA1', 'APOA2', 'APOC3', 'APOH', 'AHSG', 'FGA', 'FGB', 'FGG', 'HP', 'SERPINA1', 'TF', 'TTR', 'HPX', 'RBP4']
PAN = {'PC': 'Hepatocyte_PC (pericentral)', 'PP': 'Hepatocyte_PP (periportal)',
       'SAMmes': 'HSC_activated_SAMmes', 'SAMmac': 'SAMmac_scar_associated'}
panels = {k: [g for g in mk.gene[mk.cell_type == v]] for k, v in PAN.items()}

ok = mem.mapped_symbol_status.isin(['gencode_v49_unique_symbol_confirmed', 'gencode_v49_unambiguous_ensembl_to_symbol'])
mem['key'] = np.where(ok, mem.mapped_symbol, 'UNMAPPED_' + mem.source_gene.astype(str))
w = mem.groupby(['program_uid', 'key'], as_index=False).original_l1_weight.sum()
allw = {p: dict(zip(g.key, g.original_l1_weight)) for p, g in w.groupby('program_uid')}
progs = sorted(set(Q.program_a) | set(Q.program_b))
need = sorted(set(k for p in progs for k in allw[p] if not k.startswith('UNMAPPED_')) | set(PROXY) | set(g for v in panels.values() for g in v))
col = {g: i for i, g in enumerate(need)}


def read_sym(g):
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g['categories'][:]])
        return cats[g['codes'][:]]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in g[:]])


prim = list(dm.array_id[dm.group == 'primary'])[:int(os.environ.get('NARR', '32'))]
data = {}
for arr in prim:
    f = h5py.File(B + f'Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled/{arr}.h5ad', 'r')
    sym = read_sym(f['var']['gene_symbol']); X = f['X']
    M = sp.csr_matrix((X['data'][:], X['indices'][:], X['indptr'][:]), shape=tuple(X.attrs['shape']))
    tot = f['obs']['total_counts'][:].astype(float); ngen = f['obs']['n_genes'][:].astype(float)
    xy = f['obsm']['spatial'][:]; f.close()
    keep = np.where(np.isin(sym, need))[0]
    cidx = np.array([col[s] for s in sym[keep]])
    agg = sp.csr_matrix((np.ones(len(keep)), (np.arange(len(keep)), cidx)), shape=(len(keep), len(need)))
    C = (M[:, keep] @ agg).tocsr().astype(np.float32)
    data[arr] = dict(C=C, tot=tot, ngen=ngen, xy=xy)
print(f'loaded {len(data)} arrays in {time.time() - t0:.0f} s', flush=True)

det = {a: np.asarray((d['C'] > 0).mean(axis=0)).ravel() >= 0.01 for a, d in data.items()}
det_all = np.logical_and.reduce([det[a] for a in prim])
pan_use = {k: [g for g in v if g not in PROXY and det_all[col[g]]] for k, v in panels.items()}
for k, v in pan_use.items():
    print(f'panel {k}: {len(panels[k])} listed, {len(v)} outside proxy and detected >=1% in all 32: {v}', flush=True)


def islands(xy, k=6, mult=2.5, min_spots=8):
    tree = cKDTree(xy); dist, nb = tree.query(xy, k=k + 1)
    med = np.median(dist[:, 1][dist[:, 1] > 0]); n = len(xy)
    rows, cols = [], []
    for i in range(n):
        for dd, j in zip(dist[i, 1:], nb[i, 1:]):
            if dd <= mult * med and i != j:
                rows += [i, j]; cols += [j, i]
    A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n)); A.data[:] = 1.0
    _, comp = connected_components(A, directed=False)
    out = []
    for c in np.unique(comp):
        idx = np.where(comp == c)[0]
        if len(idx) < min_spots:
            continue
        Ai = A[idx][:, idx] + sp.identity(len(idx), format='csr')  # self-inclusive
        rs = np.asarray(Ai.sum(axis=1)).ravel()
        out.append((idx, sp.diags(1.0 / rs) @ Ai))  # row-standardized
    return out


def lees_L(Xa, Xb, isl):
    vals, wts = [], []
    for idx, Wn in isl:
        a = Xa[idx] - Xa[idx].mean(0); b = Xb[idx] - Xb[idx].mean(0)
        num = np.sum((Wn @ a) * (Wn @ b), axis=0)
        den = np.sqrt(np.sum(a * a, 0) * np.sum(b * b, 0))
        vals.append(num / den); wts.append(len(idx))
    return np.average(np.vstack(vals), axis=0, weights=np.array(wts, float))


Lmat = np.zeros((len(prim), len(Q)))
pc_rows = []
gene_zone_rows = []
for di, arr in enumerate(prim):
    d = data[arr]; C = d['C']; n = C.shape[0]
    lib = d['tot']
    N = sp.diags((1e4 / lib).astype(np.float32)) @ C; N = N.tocsr(); N.data = np.log1p(N.data)
    Z = N.toarray(); mu = Z.mean(0); sd = Z.std(0, ddof=1); sd[sd <= 0] = 1.0; Z = (Z - mu) / sd
    dm_ = det[arr]
    prox = Z[:, [col[g] for g in PROXY if dm_[col[g]]]].mean(1)
    D = np.column_stack([np.ones(n), prox, np.log1p(lib), np.log1p(d['ngen'])])
    D[:, 1:] = (D[:, 1:] - D[:, 1:].mean(0)) / D[:, 1:].std(0, ddof=1)
    Dp = np.linalg.pinv(D)
    resid = lambda V: V - D @ (Dp @ V)
    isl = islands(d['xy'])
    # pair scores, disjoint membership, weights renormalized over measured (detected) genes
    Wa = np.zeros((len(need), len(Q)), np.float32); Wb = np.zeros((len(need), len(Q)), np.float32)
    for j, (a, b) in enumerate(zip(Q.program_a, Q.program_b)):
        ga, gb = allw[a], allw[b]; sh = set(ga) & set(gb)
        for gset, Wm in ((ga, Wa), (gb, Wb)):
            ks = [k for k in gset if k not in sh and not k.startswith('UNMAPPED_') and dm_[col[k]]]
            ww = np.array([gset[k] for k in ks]); ww = ww / ww.sum()
            Wm[[col[k] for k in ks], j] = ww
    Sa = resid(Z @ Wa); Sb = resid(Z @ Wb)
    Lmat[di] = lees_L(Sa, Sb, isl)
    # positive controls: equal-weight mean z of panel genes, residualized like the programs
    ps = {k: resid(Z[:, [col[g] for g in v]].mean(1, keepdims=True)) for k, v in pan_use.items()}
    LA = lees_L(ps['PC'], ps['PP'], isl)[0]; rA = np.corrcoef(ps['PC'][:, 0], ps['PP'][:, 0])[0, 1]
    LB = lees_L(ps['SAMmes'], ps['SAMmac'], isl)[0]; rB = np.corrcoef(ps['SAMmes'][:, 0], ps['SAMmac'][:, 0])[0, 1]
    # unresidualized PC-A for comparison
    LA_raw = lees_L(Z[:, [col[g] for g in pan_use['PC']]].mean(1, keepdims=True), Z[:, [col[g] for g in pan_use['PP']]].mean(1, keepdims=True), isl)[0]
    pc_rows.append(dict(array_id=arr, n_spots=n, n_islands=len(isl), PC_A_L=LA, PC_A_r=rA, PC_A_L_unresidualized=LA_raw, PC_B_L=LB, PC_B_r=rB))
    # per-gene zone check: residual z of each zonation gene vs leave-gene-out (PC - PP) residual score
    for k in ('PC', 'PP'):
        for g in pan_use[k]:
            pcs = [x for x in pan_use['PC'] if x != g]; pps = [x for x in pan_use['PP'] if x != g]
            axis = resid((Z[:, [col[x] for x in pcs]].mean(1) - Z[:, [col[x] for x in pps]].mean(1))[:, None])[:, 0]
            gz = resid(Z[:, [col[g]]])[:, 0]
            gene_zone_rows.append(dict(array_id=arr, panel=k, gene=g, r_with_PC_minus_PP=np.corrcoef(gz, axis)[0, 1]))
    print(f'{arr} n={n} islands={len(isl)} PC-A L={LA:.4f} r={rA:.4f} | PC-B L={LB:.4f} | {time.time() - t0:.0f} s', flush=True)

pc = pd.DataFrame(pc_rows); pc.to_csv(OUT + 'positive_controls_per_donor.tsv', sep='\t', index=False)
gz = pd.DataFrame(gene_zone_rows)
gzs = gz.groupby(['panel', 'gene']).r_with_PC_minus_PP.agg(['median', lambda s: int((s > 0).sum())]).reset_index()
gzs.columns = ['panel', 'gene', 'median_r_with_PC_minus_PP', 'n_arrays_r_positive']
gzs.to_csv(OUT + 'zonation_gene_check.tsv', sep='\t', index=False)
print('PC-A: donors with L<0:', int((pc.PC_A_L < 0).sum()), '/32 | median L', round(pc.PC_A_L.median(), 4),
      '| unresidualized L<0:', int((pc.PC_A_L_unresidualized < 0).sum()), '| median', round(pc.PC_A_L_unresidualized.median(), 4))
print('PC-B: donors with L>0:', int((pc.PC_B_L > 0).sum()), '/32 | median L', round(pc.PC_B_L.median(), 4))
print(gzs.to_string())

# ---------------- label-free noise structure of pair-level L after the fixed effects
Dfe = np.load(G + 'g1_fe_design_strata.npy')
H = Dfe @ np.linalg.pinv(Dfe)
Lt = Lmat - Lmat @ H.T  # FE-residualized per donor (rows donors, cols pairs)
m = Lt.mean(0)
e = Lt - m
var_m = m.var(ddof=1); var_e = e.var(ddof=1, axis=0).mean()
sig2_pair = var_m - var_e / Lmat.shape[0]
print(f'raw pair-level L: median over donors and pairs {np.median(Lmat):.4f}, SD across pairs of donor-mean {Lmat.mean(0).std(ddof=1):.4f}')
print(f'FE-residual L: SD across pairs of cohort mean {np.sqrt(var_m):.4f} | mean within-pair SD across donors {np.sqrt(var_e):.4f} | '
      f'pair-specific SD shared by all donors {np.sqrt(max(sig2_pair, 0)):.4f} | share of pair-level variance shared across donors {max(sig2_pair, 0) / (max(sig2_pair, 0) + var_e):.3f}')

# random relabelling: same class sizes, labels carry no information
K = int(os.environ.get('NPERM', '4000')); F = 2000
flips = rng.choice([-1.0, 1.0], size=(F, Lmat.shape[0]))
res = []
base = np.r_[np.ones(n1), np.zeros(n0)]
for k in range(K):
    y = rng.permutation(base); yt = y - H @ y
    beta = Lt @ yt / (yt @ yt)  # per donor
    mb = beta.mean(); sdb = beta.std(ddof=1)
    p_flip = (1 + np.sum(flips @ beta / len(beta) >= mb)) / (F + 1)
    res.append((mb, sdb, int((beta > 0).sum()), p_flip))
R = pd.DataFrame(res, columns=['cohort_mean_beta', 'donor_sd_beta', 'n_donors_positive', 'signflip_p_one_sided'])
R.to_csv(OUT + 'random_label_null.tsv', sep='\t', index=False)
print(f'random labels (K={K}): SD of cohort-mean beta across labelings {R.cohort_mean_beta.std(ddof=1):.5f} | '
      f'median donor SD of beta {R.donor_sd_beta.median():.5f} | median donor-level SE {R.donor_sd_beta.median() / np.sqrt(32):.5f}')
print(f'fraction of random labelings with one-sided sign-flip p<0.05: {(R.signflip_p_one_sided < 0.05).mean():.3f} | '
      f'with >=22/32 donors positive: {(R.n_donors_positive >= 22).mean():.3f} | with 32/32 positive: {(R.n_donors_positive == 32).mean():.3f}')
z = 1.645 + 0.842
print(f'MDE (donor sign-flip framing, 0.44 x donor SD): {0.44 * R.donor_sd_beta.median():.5f} Lee L | '
      f'MDE against the pair-label null (2.487 x SD across labelings): {z * R.cohort_mean_beta.std(ddof=1):.5f} Lee L | '
      f'as a share of the SD of FE-residual pair-level L: {z * R.cohort_mean_beta.std(ddof=1) / np.sqrt(var_m):.3f}')
print(f'done {time.time() - t0:.0f} s')
