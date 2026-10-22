# H63 audit recount. Recomputes every n, hash and count stated in PRESPEC_H63.md.
# Computes NO zonation statistic: no center of mass, no slope, no mRNA-protein correlation,
# and reads only the shape (never the values) of the stored zon_struct/com fields.
import hashlib, json, os, sys
import h5py, numpy as np, pandas as pd

ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
R4 = ROOT + '/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z'
OUT = R4 + '/audit/H63'
ACQ = ROOT + '/Analysis/Spatial/candidates/program-context-v2-candidate-2026-08-07/acquisition/yakubovsky2026'
YK = ROOT + '/Analysis/Spatial/data/public_expansion_2026/yakubovsky2026'
SC = ROOT + '/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/scdvp'
CATP = ROOT + '/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
V = YK + '/v.mat'
SCH = ACQ + '/yakubovsky_hdf5_join_schema.tsv'
PRO = SC + '/scdvp_table3_healthy.tsv'
XLSX = SC + '/42255_2026_1459_MOESM3_ESM.xlsx'

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()

files = [
    YK + '/v.mat', YK + '/zon_struct_all_full.mat',
    ACQ + '/yakubovsky_hdf5_join_schema.tsv', ACQ + '/metadata_Clinical_Data.tsv',
    ACQ + '/yakubovsky_source_gate.tsv', ACQ + '/access_and_license_audit.tsv',
    YK + '/Human-liver/Matlab_scripts/a3_visium_12_patients_zonation_reconstruction_for_github.m',
    YK + '/Human-liver/Matlab_scripts/a2_visium_12_patients_initial_QCs_for_github.m',
    YK + '/Human-liver/Matlab_scripts/multilinear_regression_analysis.m',
    YK + '/Human-liver/Matlab_functions/extract_zonation_from_zone_index_for_github.m',
    YK + '/Human-liver/Matlab_functions/import_visium_data_funcion_for_github.m',
    YK + '/Human-liver/Matlab_functions/median_zone_filter_for_github.m',
    XLSX, PRO, CATP,
    R4 + '/gate/NG02/count.py', R4 + '/gate/NG02/count_summary.json', R4 + '/gate/NG02/marker_slopes.tsv',
    R4 + '/write/H63/count_n.py', R4 + '/write/H63/count_n.out', R4 + '/write/H63/landmarks.json',
    R4 + '/write/H63/count_detect.py', R4 + '/write/H63/count_detect.tsv',
    R4 + '/write/H63/detection_fraction.tsv.gz', R4 + '/write/H63/mde.py', R4 + '/write/H63/mde.out',
]
with open(OUT + '/sha256_recomputed.tsv', 'w') as fo:
    fo.write('bytes\tsha256\tpath\n')
    for p in files:
        fo.write(f'{os.path.getsize(p)}\t{sha(p)}\t{os.path.relpath(p, ROOT)}\n')
print('hashes written', flush=True)

# ---- v.mat roster, zones, landmarks, fields (metadata reads) ----
sch = pd.read_csv(SCH, sep='\t')
f = h5py.File(V, 'r')
def strs(ds): return [''.join(chr(c) for c in f[r][()].ravel()) for r in ds[()].ravel()]
genes0 = None; lm = {}; zones = {}; roots = {}
for d in ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8']:
    s = sch[sch.patient == d].set_index('field').hdf5_path
    root = s['mat'].rsplit('/', 1)[0]; roots[d] = root
    g = strs(f[s['gene_name']])
    if genes0 is None: genes0 = g
    z = f[s['zon_struct/zone_index']][()].ravel(); zones[d] = z
    pc = sorted(g[int(i) - 1] for i in f[root + '/LM_pc_ind'][()].ravel())
    pp = sorted(g[int(i) - 1] for i in f[root + '/LM_pp_ind'][()].ravel())
    lm[d] = (pc, pp)
    nm = ''.join(chr(c) for c in f[s['main_feature']][()].ravel())
    fields = {}
    for fld in ['mat', 'indicator_hep_gene', 'zon_struct/com', 'zon_struct_not_hep_norm/com', 'zon_struct/zone_index']:
        p = root + '/' + fld
        fields[fld] = list(f[p].shape) if p in f else None
    hep = f[root + '/indicator_hep_gene'][()].ravel() if (root + '/indicator_hep_gene') in f else None
    print(d, nm, 'same_order', g == genes0, 'n_genes', len(g), 'n_dup', len(g) - len(set(g)),
          'zone_counts', {int(k): int(v) for k, v in zip(*np.unique(z, return_counts=True))},
          'fields', fields, 'hep_indicator_sum', None if hep is None else int(np.nansum(hep)), flush=True)
lm_ident = all(lm[d] == lm['M4'] for d in lm)
LMj = json.load(open(R4 + '/write/H63/landmarks.json'))
print('LM identical across donors', lm_ident, '; equals landmarks.json',
      lm['M4'][0] == LMj['central_LM'] and lm['M4'][1] == LMj['portal_LM'])
LM = set(lm['M4'][0]) | set(lm['M4'][1])
dups = sorted({x for x in genes0 if genes0.count(x) > 1})
print('duplicated symbols', dups)
hep_series = pd.Series(f[roots['M4'] + '/indicator_hep_gene'][()].ravel(), index=genes0) if (roots['M4'] + '/indicator_hep_gene') in f else None

# ---- detection per donor (zones 1-8) ----
det = {}
for d in ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8']:
    s = sch[sch.patient == d].set_index('field').hdf5_path
    X = f[s['mat']][()]
    z = zones[d]
    if X.shape[0] != len(z): X = X.T
    keep = (z >= 1) & (z <= 8)
    det[d] = pd.Series((X[keep] > 0).mean(0), index=genes0)
    print(d, 'spots_zone1_8', int(keep.sum()), 'excluded', int((~keep).sum()),
          'neg_values', int((X < 0).sum()), 'max_bytes', X.nbytes, flush=True)
    del X
D = pd.DataFrame(det); D = D[~D.index.duplicated(keep=False)]
print('symbols kept', len(D))
old = pd.read_csv(R4 + '/write/H63/detection_fraction.tsv.gz', sep='\t', index_col=0)
cmn = [c for c in old.columns if c in D.columns]
print('detection_fraction.tsv.gz max abs diff', float((old[cmn] - D.loc[old.index, cmn]).abs().max().max()),
      'same index', list(old.index) == list(D.index))

# ---- Catalog and protein ----
cat = pd.read_parquet(CATP)
print('catalog rows', len(cat), 'symbol NA', int(cat.symbol.isna().sum()), 'dup symbols', int(cat.symbol.dropna().duplicated().sum()),
      'dup ensembl', int(cat.ensembl_id.duplicated().sum()))
print('class counts', cat.primary_evidence_class.value_counts().to_dict())
print('catalog symbols among duplicated Visium symbols', sorted(set(dups) & set(cat.symbol)))
pro0 = pd.read_csv(PRO, sep='\t')
print('protein rows', len(pro0), 'dup ids', pro0.protein[pro0.protein.duplicated(keep=False)].tolist(),
      'coef NA', int(pro0.coefficient.isna().sum()), 'zonated', pro0.zonated.value_counts().to_dict(),
      'strong', pro0.strong_zonated.value_counts().to_dict())
pro = pro0[~pro0.protein.duplicated(keep=False)]
print('protein unique kept', len(pro))
pm = pro.set_index('protein')
for g in ['CYP2E1', 'ASS1', 'GLUL', 'SDS', 'CPS1', 'PCK1']:
    print('PC-3/markers', g, None if g not in pm.index else (float(pm.coefficient[g]), float(pm.q_value[g]), pm.zonated[g]))
# xlsx vs tsv
xl = pd.ExcelFile(XLSX)
print('xlsx sheets', xl.sheet_names)
for sh in xl.sheet_names:
    t = xl.parse(sh, header=None, nrows=4)
    print('sheet', sh, 'head', t.iloc[:, :4].astype(str).values.tolist())

psym = set(pro.protein); pzon = set(pro.protein[pro.zonated == 'Yes']); pstr = set(pro.protein[pro.strong_zonated == 'Yes'])
sets = {'primary_M4_M8': ['M4', 'M5', 'M6', 'M7', 'M8'], 'sens_M2_M8': ['M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8']}
rows = []; T1 = {}
for sn, dl in sets.items():
    robust = set(D.index[(D[dl] >= 0.05).all(axis=1)])
    for cl in ['ALL', 'genetically_anchored', 'established_state_associated', 'concordant', 'unresolved', 'untested']:
        ss = set((cat if cl == 'ALL' else cat[cat.primary_evidence_class == cl]).symbol)
        inV = ss & set(D.index)
        rows.append(dict(donor_set=sn, cls=cl, n_catalog=len(ss), in_visium=len(inV), robust=len(inV & robust),
                         robust_nonLM=len((inV & robust) - LM), protein=len(ss & psym),
                         T1=len((inV & robust & psym) - LM), T1_zonated=len((inV & robust & pzon) - LM),
                         T1_strong=len((inV & robust & pstr) - LM), robust_protein_LM=len(inV & robust & psym & LM)))
        if cl == 'ALL': T1[sn] = (inV & robust & psym) - LM
out = pd.DataFrame(rows); out.to_csv(OUT + '/recount_sets.tsv', sep='\t', index=False); print(out.to_string())
print('T1 sets identical primary vs sensitivity', T1['primary_M4_M8'] == T1['sens_M2_M8'])
print('catalog&protein', len(set(cat.symbol) & psym), 'of which in Visium', len(set(cat.symbol) & psym & set(D.index)))
print('LM with protein', sorted(LM & psym), len(LM & psym))
rob = set(D.index[(D[sets['primary_M4_M8']] >= 0.05).all(axis=1)])
print('LM robust primary', len(LM & rob), 'LM robust with protein', len(LM & rob & psym))
if hep_series is not None:
    hs = hep_series[~hep_series.index.duplicated(keep=False)]
    t1 = sorted(T1['primary_M4_M8'])
    print('T1 hepatocyte-gene split', int((hs.reindex(t1) == 1).sum()), 'hep;', int((hs.reindex(t1) != 1).sum()), 'other')
# detection distribution of T1 set (diagnostic 7, not outcome)
t1 = sorted(T1['primary_M4_M8'])
print('T1 detection min per donor', D.loc[t1, sets['primary_M4_M8']].min().round(4).to_dict())
print('done')
