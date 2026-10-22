#!/usr/bin/env python3
"""H63: per-gene central-portal position of Gene Catalog r8 genes in healthy human liver.

mRNA layer: live-donor Visium (Yakubovsky et al. 2026, Zenodo 17735587 v.mat), center of
mass (COM) of the zone-mean profile over the source's zones 1-8, per donor, median over the
five non-steatotic donors M4-M8.
Protein layer: Weiss et al. 2026 Nat Metab Supplementary Table 3 (healthy cohort, N=14),
copied, not refit.
Agreement test T1: Spearman rho between the two layers over the 1,337-gene T1 set.

Implements PRESPEC_H63.md (sha256 b638c538...) exactly: stop rules S1-S7 in order, PC-1..3,
NC-1, NC-2, diagnostics 1-9 before T1, T1, D1-D8, per-gene Catalog columns.

--debug: pre-outcome debugging only. Shuffles the Weiss values among non-landmark proteins
(so every agreement statistic is computed on null data), reduces B, and does not stop on S7.
"""
import argparse
import hashlib
import json
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
H63DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRESPEC_SHA = 'b638c538046cc41ce69c1a6c7a057b4e8f65e3b58466cb5469ea40e019453dd9'
PRESPEC_ORIG = ROOT + '/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H63.md'
R4 = 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z'
ACQ = 'Analysis/Spatial/candidates/program-context-v2-candidate-2026-08-07/acquisition/yakubovsky2026'
YK = 'Analysis/Spatial/data/public_expansion_2026/yakubovsky2026'
SC = 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/scdvp'
CATP = 'figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'

# (role, path relative to ROOT, bytes, sha256) exactly as in the prespec inputs table
INPUTS = [
    ('visium_v_mat', YK + '/v.mat', 3996202856, '5fb89af61409c575a713aeacfd570551fbe109b0695915cd61dc0422f651b5d8'),
    ('hdf5_schema', ACQ + '/yakubovsky_hdf5_join_schema.tsv', 15757, '2ee81ecb95de0fed3c52dba46007de291c1b6f6f0ae2dbc2ff61cd44dd50d454'),
    ('clinical', ACQ + '/metadata_Clinical_Data.tsv', 1221, '2ad5b9048bcbe8de1ecd3327eb761ea85ab82a788a2839a215e3833b3d787e2e'),
    ('lipid_gate', ACQ + '/yakubovsky_source_gate.tsv', 1666, 'dedc0ecb6146f45c8afa674d36841014bc9b9f8d70e84d7d584ee05a64942622'),
    ('license', ACQ + '/access_and_license_audit.tsv', 807, 'e0921513135ad19bf8628f31bd99a63e0d53e446836066284459d9416e162f33'),
    ('src_a3', YK + '/Human-liver/Matlab_scripts/a3_visium_12_patients_zonation_reconstruction_for_github.m', 7113, '0cc9f5311ea6147768c323ac51d56e1ca6a2c471371229737a95677b34b8c6f9'),
    ('src_a2', YK + '/Human-liver/Matlab_scripts/a2_visium_12_patients_initial_QCs_for_github.m', 7064, '34438b074ebf536449b154cfcc7d2af9ad3e83da0d5419342fab3af752982839'),
    ('src_histology', YK + '/Human-liver/Matlab_scripts/multilinear_regression_analysis.m', 9609, '8229fe467f80cb8039bd6bd27a0ae7323ab897e51d78438904c4e4decbd2a5ac'),
    ('src_com', YK + '/Human-liver/Matlab_functions/extract_zonation_from_zone_index_for_github.m', 2993, 'a855624ea710cb0b7d2fb34e6aaf2bc66ba4f91e7f5648c05c4ef6e5a5839d75'),
    ('src_bg', YK + '/Human-liver/Matlab_functions/import_visium_data_funcion_for_github.m', 1649, '3f5246a045af536ff3264c455a98875bad704ad3c73809c445816f81e72f10e7'),
    ('src_medfilt', YK + '/Human-liver/Matlab_functions/median_zone_filter_for_github.m', 1497, 'df102876dba3d13ccb70bf4e5e637d94673625492c70c1fa977428ad8d246c9d'),
    ('src_pooled_profile', YK + '/zon_struct_all_full.mat', 203729828, '7764fadcb706d93302cc3e7ce688fb749d998938dc35b463467642a3f0c633f5'),
    ('weiss_xlsx', SC + '/42255_2026_1459_MOESM3_ESM.xlsx', 638091, '71511d6c79d37bf297c95200737bc968a264d2c32031ca9e61097fa7f8e0efd0'),
    ('weiss_tsv', SC + '/scdvp_table3_healthy.tsv', 759265, 'bf560449dbda615d8039fe3896780372384cd2045ca1e86fd1354f6d63285320'),
    ('catalog_r8', CATP, 468456, 'ffd3c4650c77a899e9ebed7327697398c4b6462cfffcd7f0e635cd0d023db06b'),
    ('gate_count_py', R4 + '/gate/NG02/count.py', 3839, 'd16a5ea86ed92bc83aa5b1ce58f72954b4f6fb85e0da2bf77a318b9e723a0385'),
    ('gate_summary', R4 + '/gate/NG02/count_summary.json', 6064, 'f99555f1310fb8bd2bdd784554b3a370c429a565267d1b16593511f18ab56c6e'),
    ('gate_marker_slopes', R4 + '/gate/NG02/marker_slopes.tsv', 1720, 'acb46a0ec1aa643685df7c6a87157d57934729b4f2d94d08521348215474fc84'),
    ('count_n_py', R4 + '/write/H63/count_n.py', 2985, 'd0f099ae0929441ebc64c68a08be8bc95bab0d736a3c1d39a6c17b933eb58ba9'),
    ('count_n_out', R4 + '/write/H63/count_n.out', 2541, '3ec9aaf274598d12a07dee6cc21f1d8beb2431ee24c22bb43f0f6f18bc5b82cb'),
    ('landmarks', R4 + '/write/H63/landmarks.json', 466, 'eeb17d8ae4c5fe1035c63c16b68573240dd83cc128804f9e4f41cac0276063f8'),
    ('count_detect_py', R4 + '/write/H63/count_detect.py', 3237, '724c9a760e4ee5b4321cbe62ef549e782424766d274ef8a043ce2b1a214f642c'),
    ('count_detect_tsv', R4 + '/write/H63/count_detect.tsv', 998, 'a699f6f8c1e691c8161efb0f3abdd49ac514f68404f232f10f6168414c1a604b'),
    ('detection_fraction', R4 + '/write/H63/detection_fraction.tsv.gz', 779236, '9b5bbbafc268636906f0df9647d741c1c78aed5a38175288cddb6d0af7f6ede4'),
    ('mde_py', R4 + '/write/H63/mde.py', 705, '4c23d82a8052d727d458c44da7ef4dbb4922b5698b81dc5e75f8721ad0a0c0dd'),
    ('mde_out', R4 + '/write/H63/mde.out', 353, '4459c709d133826d1687622f0f4f18ce67fcb57f7e52ef2f8f2f85ede8d42479'),
]
P = {r: ROOT + '/' + p for r, p, _, _ in INPUTS}

ALL_DONORS = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8']
PRIMARY = ['M4', 'M5', 'M6', 'M7', 'M8']
SENS = ['M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8']
EXPECTED_LABEL = {d: f'{d}-steatosis' for d in ['M1', 'M2', 'M3']}
EXPECTED_LABEL.update({d: f'{d}-normal' for d in PRIMARY})
# per-zone spot counts (count_n.out, hashed input); zone 9 only in M8
EXPECTED_ZONES = {
    'M1': {1: 299, 2: 461, 3: 506, 4: 571, 5: 615, 6: 607, 7: 492, 8: 441},
    'M2': {1: 110, 2: 210, 3: 228, 4: 272, 5: 271, 6: 299, 7: 246, 8: 170},
    'M3': {1: 153, 2: 346, 3: 382, 4: 395, 5: 410, 6: 371, 7: 379, 8: 285},
    'M4': {1: 145, 2: 293, 3: 294, 4: 343, 5: 346, 6: 301, 7: 298, 8: 249},
    'M5': {1: 98, 2: 211, 3: 246, 4: 251, 5: 234, 6: 230, 7: 238, 8: 158},
    'M6': {1: 69, 2: 185, 3: 212, 4: 200, 5: 197, 6: 158, 7: 196, 8: 141},
    'M7': {1: 126, 2: 272, 3: 325, 4: 373, 5: 333, 6: 345, 7: 308, 8: 217},
    'M8': {1: 157, 2: 310, 3: 317, 4: 279, 5: 353, 6: 371, 7: 332, 8: 225, 9: 213},
}
# zone 1-8 totals stated in the prespec n section
EXPECTED_Z18 = {'M2': 1806, 'M3': 2721, 'M4': 2269, 'M5': 1666, 'M6': 1358, 'M7': 2299, 'M8': 2344}
CLASSES = ['genetically_anchored', 'established_state_associated', 'concordant', 'unresolved', 'untested']
DETECT_MIN = 0.05
TAU_FLOOR = 0.25
RHO_FLOOR = 0.30
ALPHA = 0.05

# ----------------------------------------------------------------------------------------
LOG_T0 = time.time()
DIAG = []          # rows of diagnostics.tsv
REACHED = {}       # numbers reached, for RESULT_H63.md
STOPS = {}         # stop-rule status


def log(*a):
    print(f'[{time.time() - LOG_T0:7.1f}s]', *a, flush=True)


def diag(n, item, key, value):
    if isinstance(value, float):
        value = repr(value)
    DIAG.append((n, item, key, str(value)))
    log(f'D{n}', item, key, value)


def write_diag(out):
    pd.DataFrame(DIAG, columns=['diagnostic', 'item', 'key', 'value']).to_csv(
        out + '/diagnostics.tsv', sep='\t', index=False)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def spear(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return np.nan, int(ok.sum())
    ra = rankdata(a[ok]); rb = rankdata(b[ok])
    if ra.std() == 0 or rb.std() == 0:
        return np.nan, int(ok.sum())
    return float(np.corrcoef(ra, rb)[0, 1]), int(ok.sum())


def boot_ci(rng, a, b, B):
    n = len(a)
    vals = np.empty(B)
    for i in range(B):
        idx = rng.integers(0, n, n)
        vals[i] = spear(a[idx], b[idx])[0]
    lo, hi = np.nanpercentile(vals, [2.5, 97.5])
    return float(lo), float(hi), vals


def com_centered(sums_by_zone, counts):
    """sums_by_zone: 8 x G zone sums of spot-normalized values; counts: 8 zone spot counts.
    Returns c = COM - 4.5 (NaN when the zone-mean profile sums to 0)."""
    m = sums_by_zone / counts[:, None]
    s = m.sum(0)
    com = (np.arange(1, 9, dtype=float)[:, None] * m).sum(0) / np.where(s > 0, s, np.nan)
    return com - 4.5


def onehot(z):
    M = np.zeros((8, len(z)))
    M[z.astype(int) - 1, np.arange(len(z))] = 1.0
    return M


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return 'NA'
    if isinstance(x, (int, np.integer)):
        return f'{int(x):,}'
    return f'{x:.{nd}f}'


def fmtp(p):
    return f'{p:.4g}' if p >= 1e-3 else f'{p:.3e}'


# ----------------------------------------------------------------------------------------
def write_result(out, verdict, rule, detail, extra_sections, debug):
    lines = []
    lines.append('# RESULT_H63' + (' (DEBUG RUN ON NULL PROTEIN DATA, NOT A RESULT)' if debug else ''))
    lines.append('')
    lines.append(f'Prespecification: `PRESPEC_H63.md` sha256 `{PRESPEC_SHA}`. Run directory: `{out}`.')
    lines.append(f'Written {time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())} by `scripts/run_h63.py`; verdict computed in code by the prespec reading rule.')
    lines.append('')
    lines.append(f'## Verdict: {verdict}')
    lines.append('')
    if rule:
        lines.append(f'Stop rule fired: **{rule}**. {detail}')
    else:
        lines.append(detail)
    lines.append('')
    lines.append('## Stop rules (checked in order)')
    lines.append('')
    lines.append('| Rule | Status | Number |')
    lines.append('|---|---|---|')
    for k in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7']:
        st, num = STOPS.get(k, ('not reached', ''))
        lines.append(f'| {k} | {st} | {num} |')
    lines.append('')
    for title, body in extra_sections:
        lines.append(f'## {title}')
        lines.append('')
        lines.extend(body)
        lines.append('')
    lines.append('## Numbers reached (machine-readable copy in `numbers_reached.json`)')
    lines.append('')
    for k, v in REACHED.items():
        lines.append(f'- {k}: {v}')
    lines.append('')
    with open(out + '/RESULT_H63.md', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    with open(out + '/numbers_reached.json', 'w') as fh:
        json.dump({'verdict': verdict, 'rule': rule, 'stops': STOPS, 'reached': REACHED}, fh, indent=1, default=str)
    write_diag(out)


def stop(out, rule, detail, debug):
    STOPS[rule] = ('FIRED', detail)
    log('STOP', rule, detail)
    write_result(out, f'Untestable ({rule})', rule,
                 detail + ' No column enters the Catalog; nothing later was computed.', [], debug)
    sys.exit(0)


# ----------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--debug', action='store_true')
    a = ap.parse_args()
    out = a.out
    debug = a.debug
    if a.seed != 20260923:
        raise SystemExit('seed must be 20260923 (prespec)')
    B_NC1 = 10 if debug else 200
    B_PERM = 200 if debug else 10000
    B_BOOT = 50 if debug else 2000
    log('H63 run', 'DEBUG' if debug else 'FINAL', 'out', out, 'B_NC1', B_NC1, 'B_PERM', B_PERM, 'B_BOOT', B_BOOT)

    # One SeedSequence, four streams (prespec "Seed")
    children = np.random.SeedSequence(a.seed).spawn(4)
    rng0, rng1, rng2, rng3 = (np.random.default_rng(c) for c in children)

    # ============================ S1: hashes, PC-3 ============================
    pre_copy = H63DIR + '/PRESPEC_H63.md'
    h_copy = sha256(pre_copy); h_orig = sha256(PRESPEC_ORIG)
    diag(1, 'prespec', 'sha256_run_copy', h_copy)
    diag(1, 'prespec', 'sha256_plans_original', h_orig)
    if h_copy != PRESPEC_SHA or h_orig != PRESPEC_SHA:
        STOPS['S1'] = ('FIRED', 'prespec hash mismatch')
        stop(out, 'S1', f'PRESPEC_H63.md sha256 {h_copy} / {h_orig} differs from {PRESPEC_SHA}.', debug)
    bad = []
    rows = []
    for role, rel, nbytes, exp in INPUTS:
        pth = ROOT + '/' + rel
        if not os.path.exists(pth):
            bad.append(f'{rel} missing'); continue
        sz = os.path.getsize(pth); h = sha256(pth)
        rows.append((role, rel, sz, h, nbytes, exp, sz == nbytes and h == exp))
        diag(1, 'sha256', rel, f'{h} bytes={sz} match={sz == nbytes and h == exp}')
        if sz != nbytes or h != exp:
            bad.append(f'{rel}: bytes {sz} (expected {nbytes}), sha256 {h} (expected {exp})')
    pd.DataFrame(rows, columns=['role', 'path', 'bytes', 'sha256', 'bytes_expected', 'sha256_expected', 'match']).to_csv(
        out + '/input_sha256.tsv', sep='\t', index=False)
    if bad:
        stop(out, 'S1', 'Input hash mismatch: ' + '; '.join(bad), debug)

    pro0 = pd.read_csv(P['weiss_tsv'], sep='\t')
    pro0['protein'] = pro0['protein'].astype(str)
    diag(1, 'protein', 'rows', len(pro0))
    diag(1, 'protein', 'distinct_ids', pro0.protein.nunique())
    dup_ids = sorted(set(pro0.protein[pro0.protein.duplicated(keep=False)]))
    diag(1, 'protein', 'duplicated_ids_dropped', ','.join(dup_ids))
    pro = pro0[~pro0.protein.duplicated(keep=False)].set_index('protein')
    diag(1, 'protein', 'ids_kept', len(pro))
    pc3 = {}
    pc3_ok = True
    for g, coef_exp, q_exp in [('CYP2E1', -1.82, '6.0e-34'), ('ASS1', 1.89, '1.4e-74')]:
        c = float(pro.loc[g, 'coefficient']); q = float(pro.loc[g, 'q_value'])
        ok = (round(c, 2) == coef_exp) and (f'{q:.1e}' == q_exp)
        pc3[g] = (c, q, ok); pc3_ok &= ok
        diag(1, 'PC-3', g, f'coef={c:.4f} q={q:.3e} expected coef {coef_exp} q {q_exp} pass={ok}')
    REACHED['PC-3'] = {g: f'coef {v[0]:.3f}, q {v[1]:.2e}, pass {v[2]}' for g, v in pc3.items()}
    if not pc3_ok:
        stop(out, 'S1', 'PC-3 failed: protein table values differ from CYP2E1 -1.82 (q 6.0e-34) / ASS1 +1.89 (q 1.4e-74).', debug)
    STOPS['S1'] = ('clear', f'{len(INPUTS)} input hashes and byte counts match; prespec sha256 matches; PC-3 CYP2E1 {pc3["CYP2E1"][0]:.2f}, ASS1 {pc3["ASS1"][0]:+.2f}')
    write_diag(out)

    # ============================ S2: donor roster, zones, landmarks ============================
    LMj = json.load(open(P['landmarks']))
    LM_C = list(LMj['central_LM']); LM_P = list(LMj['portal_LM'])
    LM = set(LM_C) | set(LM_P)

    if debug:
        # shuffle Weiss values jointly among non-landmark proteins: every agreement number becomes null
        drng = np.random.default_rng(999)
        nonlm = [p for p in pro.index if p not in LM]
        cols = ['coefficient', 'p_value', 'q_value', 'zonated', 'strong_zonated']
        perm = drng.permutation(len(nonlm))
        pro.loc[nonlm, cols] = pro.loc[nonlm, cols].values[perm]
        log('DEBUG: Weiss values shuffled among', len(nonlm), 'non-landmark proteins')

    sch = pd.read_csv(P['hdf5_schema'], sep='\t')
    f = h5py.File(P['visium_v_mat'], 'r')

    def strs(ds):
        return [''.join(chr(c) for c in f[r][()].ravel()) for r in ds[()].ravel()]

    def paths(d):
        s = sch[sch.patient == d].set_index('field').hdf5_path
        return s, s['mat'].rsplit('/', 1)[0]

    genes = None
    zones = {}
    s2_bad = []
    hep = None
    for d in ALL_DONORS:
        s, root = paths(d)
        g = strs(f[s['gene_name']])
        if genes is None:
            genes = g
        if g != genes:
            s2_bad.append(f'{d} gene order differs')
        label = ''.join(chr(c) for c in f[s['main_feature']][()].ravel())
        diag(2, 'donor_label', d, label)
        if label != EXPECTED_LABEL[d]:
            s2_bad.append(f'{d} label {label!r} != {EXPECTED_LABEL[d]!r}')
        z = f[s['zon_struct/zone_index']][()].ravel()
        zones[d] = z
        u, cnt = np.unique(z, return_counts=True)
        zc = {int(k): int(v) for k, v in zip(u, cnt)}
        diag(2, 'zone_counts', d, json.dumps(zc))
        if not np.all(np.isin(z, np.arange(1, 10))):
            s2_bad.append(f'{d} zone values outside 1..9')
        if zc != EXPECTED_ZONES[d]:
            s2_bad.append(f'{d} zone counts {zc} != {EXPECTED_ZONES[d]}')
        z18 = int(((z >= 1) & (z <= 8)).sum())
        if d in EXPECTED_Z18 and z18 != EXPECTED_Z18[d]:
            s2_bad.append(f'{d} zone 1-8 spots {z18} != {EXPECTED_Z18[d]}')
        pc = sorted(g[int(i) - 1] for i in f[root + '/LM_pc_ind'][()].ravel())
        pp = sorted(g[int(i) - 1] for i in f[root + '/LM_pp_ind'][()].ravel())
        same = (pc == LM_C and pp == LM_P)
        diag(2, 'landmarks_equal_landmarks_json', d, same)
        if not same:
            s2_bad.append(f'{d} landmarks differ from landmarks.json')
        h = f[root + '/indicator_hep_gene'][()].ravel()
        if hep is None:
            hep = h
        elif not np.array_equal(h, hep):
            log('note: indicator_hep_gene differs in', d)
    if s2_bad:
        stop(out, 'S2', 'Donor roster check failed: ' + '; '.join(s2_bad), debug)
    STOPS['S2'] = ('clear', 'M4-M8 labelled normal, M1-M3 steatosis; zone counts match (M8 213 zone-9 spots excluded); 40 landmarks identical in all 8 donors and equal landmarks.json')
    REACHED['zone 1-8 spots'] = {d: int(((zones[d] >= 1) & (zones[d] <= 8)).sum()) for d in ALL_DONORS}
    write_diag(out)

    G = len(genes)
    gser = pd.Series(genes)
    dupmask = gser.duplicated(keep=False).values
    dup_syms = sorted(set(gser[dupmask]))
    keep_idx = np.where(~dupmask)[0]
    kept = [genes[i] for i in keep_idx]
    pos = {genes[i]: i for i in keep_idx}         # symbol -> column index (non-duplicated only)
    diag(2, 'genes', 'n_visium_symbols', G)
    diag(2, 'genes', 'duplicated_symbols_dropped', ','.join(dup_syms))
    diag(2, 'genes', 'n_kept', len(kept))

    # robust set (pre-outcome, from the hashed detection file) used to restrict NC-1 matrices
    det_file = pd.read_csv(P['detection_fraction'], sep='\t', index_col=0)
    # the prespec's robust set (4,038) is Catalog r8 symbols in Visium with detection >= 0.05 in M4-M8
    cat_syms = set(pd.read_parquet(P['catalog_r8']).symbol.astype(str))
    robust_file = [g for g in det_file.index[(det_file[PRIMARY] >= DETECT_MIN).all(axis=1)] if g in cat_syms]
    rob_cols = np.array([pos[g] for g in robust_file])

    # ============================ read matrices: COM, detection, slope ============================
    C = {}; DET = {}; SLOPE = {}; STORED = {}; NC2 = {}; XROB = {}; ZROB = {}
    for d in ALL_DONORS:
        s, root = paths(d)
        X = f[s['mat']][()]
        z = zones[d]
        if X.shape[0] != len(z):
            X = X.T
        keep = (z >= 1) & (z <= 8)
        if not keep.all():
            X = X[keep]
        z = z[keep].astype(int)
        DET[d] = (X > 0).mean(0)
        tot = X.sum(1)                                      # spot total over all 36,601 genes
        NC2[d] = spear(tot, z)[0]
        if (tot <= 0).any():
            raise SystemExit(f'{d}: spot with zero total')
        X /= tot[:, None]                                   # x_gs
        oh = onehot(z); cnt = oh.sum(1)
        C[d] = com_centered(oh @ X, cnt)
        zc_ = z - z.mean()
        mean_g = X.mean(0)
        SLOPE[d] = (zc_ @ X) / (zc_ @ zc_) / np.where(mean_g > 0, mean_g, np.nan)
        st = {}
        for fld in ['zon_struct/com', 'zon_struct_not_hep_norm/com']:
            pth = root + '/' + fld
            st[fld] = f[pth][()].ravel() if pth in f else None
        STORED[d] = st
        if d in PRIMARY:
            XROB[d] = np.ascontiguousarray(X[:, rob_cols]); ZROB[d] = z.copy()
        log(d, 'read', X.shape, 'zones 1-8 spots', len(z), 'NC-2 rho(spot total, zone)', round(NC2[d], 4))
        del X
    f.close()

    # ============================ S3: set sizes ============================
    Dd = pd.DataFrame({d: DET[d][keep_idx] for d in ALL_DONORS}, index=kept)
    maxdiff = float((det_file[SENS] - Dd.loc[det_file.index, SENS]).abs().max().max())
    diag(3, 'detection', 'max_abs_diff_vs_detection_fraction_tsv', maxdiff)
    cat = pd.read_parquet(P['catalog_r8'])
    cat['symbol'] = cat['symbol'].astype(str)
    psym = set(pro.index); pzon = set(pro.index[pro.zonated == 'Yes']); pstr = set(pro.index[pro.strong_zonated == 'Yes'])
    keptset = set(kept)
    robust = {}
    for sn, dl in [('primary_M4_M8', PRIMARY), ('sens_M2_M8', SENS)]:
        robust[sn] = set(Dd.index[(Dd[dl] >= DETECT_MIN).all(axis=1)])
    T1sets = {}
    cd = pd.read_csv(P['count_detect_tsv'], sep='\t')
    mism = []
    for sn in ['primary_M4_M8', 'sens_M2_M8']:
        rb = robust[sn]
        for cl in ['ALL'] + CLASSES:
            ss = set((cat if cl == 'ALL' else cat[cat.primary_evidence_class == cl]).symbol)
            inV = ss & keptset
            r = dict(in_visium=len(inV), robust=len(inV & rb), robust_nonLM=len((inV & rb) - LM), protein=len(ss & psym),
                     T1=len((inV & rb & psym) - LM), T1_zonated=len((inV & rb & pzon) - LM),
                     T1_strong=len((inV & rb & pstr) - LM), robust_protein_LM=len(inV & rb & psym & LM))
            diag(3, f'set_sizes_{sn}', cl, json.dumps(r))
            if not debug:
                ref = cd[(cd.donor_set == sn) & (cd.cls == cl)].iloc[0]
                for k_mine, k_ref in [('in_visium', 'in_visium'), ('robust', 'robust'), ('robust_nonLM', 'robust_nonLM'),
                                      ('protein', 'protein'), ('T1', 'robust_protein_nonLM'),
                                      ('T1_zonated', 'robust_protein_zonated_nonLM'),
                                      ('T1_strong', 'robust_protein_strong_nonLM'), ('robust_protein_LM', 'robust_protein_LM')]:
                    if r[k_mine] != int(ref[k_ref]):
                        mism.append(f'{sn}/{cl}/{k_mine} {r[k_mine]} vs {int(ref[k_ref])}')
            if cl == 'ALL':
                T1sets[sn] = (inV & rb & psym) - LM
    rob_cat = robust['primary_M4_M8'] & set(cat.symbol)
    n_rob = len(rob_cat); n_T1 = len(T1sets['primary_M4_M8'])
    diag(3, 'S3', 'robust_primary_all_visium_symbols_incl_non_catalog', len(robust['primary_M4_M8']))
    diag(3, 'S3', 'robust_primary_catalog', n_rob)
    diag(3, 'S3', 'T1_set', n_T1)
    diag(3, 'S3', 'T1_identical_primary_vs_sensitivity', T1sets['primary_M4_M8'] == T1sets['sens_M2_M8'])
    diag(3, 'S3', 'robust_primary_catalog_identical_to_detection_file', rob_cat == set(robust_file))
    diag(3, 'S3', 'count_detect_tsv_mismatches', ';'.join(mism) if mism else 'none')
    REACHED['robust primary genes'] = n_rob
    REACHED['T1 set'] = n_T1
    if n_rob != 4038 or n_T1 != 1337 or rob_cat != set(robust_file):
        stop(out, 'S3', f'Set sizes: robust primary {n_rob} (required 4,038), T1 {n_T1} (required 1,337).', debug)
    STOPS['S3'] = ('clear', f'robust primary {n_rob:,}; T1 set {n_T1:,}')
    write_diag(out)

    rob_list = list(robust_file)          # order of NC-1 columns
    rob_pos = {g: i for i, g in enumerate(rob_list)}

    # ============================ S4: COM reproduction ============================
    s4 = {}
    s4_bad = []
    for d in ALL_DONORS:
        vals = {}
        for fld, arr in STORED[d].items():
            if arr is None or len(arr) != G:
                vals[fld] = (np.nan, 0)
                continue
            vals[fld] = spear((C[d] + 4.5)[rob_cols], arr[rob_cols])
        best = np.nanmax([v[0] for v in vals.values()]) if any(np.isfinite(v[0]) for v in vals.values()) else np.nan
        s4[d] = best
        diag(4, 'COM_reproduction_spearman', d, json.dumps({k: [round(v[0], 5) if np.isfinite(v[0]) else None, v[1]] for k, v in vals.items()}) + f' best={best:.5f}' + ('' if d in PRIMARY else ' (not a stop-rule donor)'))
        if d in PRIMARY:
            if not np.isfinite(best):
                s4_bad.append(f'{d}: neither stored COM field present with length {G}')
            elif best < 0.90:
                s4_bad.append(f'{d}: best Spearman {best:.4f} < 0.90')
    REACHED['S4 COM reproduction (best Spearman per primary donor)'] = {d: round(float(s4[d]), 4) for d in PRIMARY}
    if s4_bad:
        stop(out, 'S4', 'COM reproduction failed: ' + '; '.join(s4_bad), debug)
    STOPS['S4'] = ('clear', 'min over M4-M8 of best Spearman = ' + f'{min(s4[d] for d in PRIMARY):.4f} (rule >= 0.90)')
    write_diag(out)

    # per-gene matrices over kept symbols
    Cdf = pd.DataFrame({d: C[d][keep_idx] for d in ALL_DONORS}, index=kept)
    Sdf = pd.DataFrame({d: SLOPE[d][keep_idx] for d in ALL_DONORS}, index=kept)
    cprim = Cdf[PRIMARY]
    c_g = cprim.median(axis=1, skipna=True)
    k_neg = (cprim < 0).sum(axis=1); k_pos = (cprim > 0).sum(axis=1); n_def = cprim.notna().sum(axis=1)
    c_g_sens = Cdf[SENS].median(axis=1, skipna=True)
    coef = pro['coefficient'].astype(float)

    # ============================ PC-1 ============================
    pc1 = {}
    pc1_ok = True
    for d in PRIMARY:
        mc = float(np.nanmedian(Cdf.loc[LM_C, d])); mp = float(np.nanmedian(Cdf.loc[LM_P, d]))
        gl = float(Cdf.loc['GLUL', d])
        ok = (mc < 0) and (mp > 0) and (gl < 0)
        pc1[d] = (mc, mp, gl, ok); pc1_ok &= ok
        diag(5, 'PC-1', d, f'median c central LM {mc:+.4f}; median c portal LM {mp:+.4f}; GLUL c {gl:+.4f}; pass={ok}')
    REACHED['PC-1'] = {d: f'central LM {v[0]:+.3f}, portal LM {v[1]:+.3f}, GLUL {v[2]:+.3f}, pass {v[3]}' for d, v in pc1.items()}
    if not pc1_ok:
        STOPS['S5'] = ('FIRED', 'PC-1 failed')
        stop(out, 'S5', 'PC-1 (orientation) failed in at least one primary donor.', debug)
    STOPS['S5'] = ('clear', 'PC-1 passed in 5 of 5 primary donors')

    # ============================ PC-2 ============================
    lm22 = sorted(LM & rob_cat & psym)
    agree = 0
    pc2_rows = []
    for g in lm22:
        same = np.sign(c_g[g]) == np.sign(coef[g])
        agree += int(same)
        pc2_rows.append((g, 'central' if g in LM_C else 'portal', float(c_g[g]), float(coef[g]), bool(same)))
        diag(5, 'PC-2 gene', g, f'class={"central" if g in LM_C else "portal"} c_g={c_g[g]:+.4f} weiss_coef={coef[g]:+.4f} same_sign={same}')
    pd.DataFrame(pc2_rows, columns=['symbol', 'landmark_class', 'c_g_M4_M8', 'weiss_coefficient', 'same_sign']).to_csv(
        out + '/pc2_landmarks.tsv', sep='\t', index=False)
    diag(5, 'PC-2', 'count_same_sign', f'{agree} of {len(lm22)} (rule >= 15; predicted 19)')
    REACHED['PC-2'] = f'{agree} of {len(lm22)} landmark genes with same sign (rule >= 15 of 22; predicted 19)'
    if agree < 15:
        stop(out, 'S6', f'PC-2 failed: {agree} of {len(lm22)} landmark genes have matching sign (rule >= 15).', debug)
    STOPS['S6'] = ('clear', f'PC-2 {agree} of {len(lm22)}')
    write_diag(out)

    # ============================ NC-1 and NC-2 ============================
    T1 = sorted(T1sets['primary_M4_M8'])
    t1_rob_idx = np.array([rob_pos[g] for g in T1])
    coef_T1 = coef.loc[T1].values
    rp_T1 = rankdata(coef_T1)
    nR = len(rob_list)
    null_cg = np.empty((B_NC1, nR)); null_5of5 = np.empty(B_NC1); null_rho = np.empty(B_NC1)
    null_side5 = np.empty((B_NC1, nR), dtype=bool)
    for b in range(B_NC1):
        cb = np.empty((len(PRIMARY), nR))
        for j, d in enumerate(PRIMARY):
            z = ZROB[d]
            zp = z[rng1.permutation(len(z))]
            oh = onehot(zp)
            cb[j] = com_centered(oh @ XROB[d], oh.sum(1))
        cg = np.median(cb, axis=0)
        side5 = ((cb < 0).sum(0) == 5) | ((cb > 0).sum(0) == 5)
        null_cg[b] = cg; null_side5[b] = side5
        null_5of5[b] = side5.mean()
        null_rho[b] = spear(cg[t1_rob_idx], coef_T1)[0]
    frac5 = float(null_5of5.mean())
    q99 = float(np.quantile(np.abs(null_cg).ravel(), 0.99))
    frac_rho = float((np.abs(null_rho) >= 0.08).mean())
    frac_rho_lt = float((np.abs(null_rho) < 0.08).mean())
    diag(6, 'NC-1', 'mean_fraction_5of5_one_side', f'{frac5:.4f} (predicted 0.04-0.09; rule 0.03-0.10)')
    diag(6, 'NC-1', 'pooled_q99_abs_c_g', f'{q99:.4f} (predicted < 0.10; auditor estimate 0.10-0.40; not a stop rule)')
    diag(6, 'NC-1', 'fraction_perms_abs_rho_ge_0.08', f'{frac_rho:.4f} (rule <= 0.10; predicted |rho| < 0.08 in >= 95%)')
    diag(6, 'NC-1', 'abs_rho_null_max', f'{np.nanmax(np.abs(null_rho)):.4f}')
    pd.DataFrame({'perm': np.arange(1, B_NC1 + 1), 'fraction_5of5': null_5of5, 'rho_T1_vs_protein': null_rho}).to_csv(
        out + '/nc1_permutations.tsv', sep='\t', index=False)
    for d in ALL_DONORS:
        diag(6, 'NC-2 spearman(spot total, zone)', d, f'{NC2[d]:+.4f}' + (' (primary)' if d in PRIMARY else ''))
    nc2_flag = any(abs(NC2[d]) > 0.3 for d in PRIMARY)
    diag(6, 'NC-2', 'slope_flagged_|rho|>0.3_any_primary', nc2_flag)

    # tau_b per decile of minimum detection (detection_fraction.tsv.gz, fixed before the run)
    mindet_rob = det_file.loc[rob_list, PRIMARY].min(axis=1)
    bins, edges = pd.qcut(mindet_rob, 10, labels=False, retbins=True, duplicates='raise')
    bins = bins.values.astype(int)
    tau = {}
    tau_rows = []
    absnull = np.abs(null_cg)
    for bb in range(10):
        cols = np.where(bins == bb)[0]
        qb = float(np.quantile(absnull[:, cols].ravel(), 0.99))
        tb = max(TAU_FLOOR, qb)
        tau[bb] = tb
        fc = float((null_side5[:, cols] & (absnull[:, cols] >= tb)).mean())
        tau_rows.append((bb + 1, edges[bb], edges[bb + 1], len(cols), qb, tb, fc))
        diag(6, 'tau_b', f'bin{bb + 1}', f'min_det [{edges[bb]:.4f},{edges[bb + 1]:.4f}] n={len(cols)} q99_null={qb:.4f} tau_b={tb:.4f} null_false_call_rate={fc:.5f} (predicted < 0.01)')
    pd.DataFrame(tau_rows, columns=['bin', 'min_det_lo', 'min_det_hi', 'n_genes', 'q99_null_abs_c', 'tau_b', 'nc1_false_call_rate']).to_csv(
        out + '/tau_by_detection_decile.tsv', sep='\t', index=False)

    s7 = (frac5 < 0.03) or (frac5 > 0.10) or (frac_rho > 0.10)
    REACHED['NC-1'] = f'5-of-5 fraction {frac5:.4f} (rule 0.03-0.10); pooled q99 |c_g| {q99:.3f}; perms with |rho| >= 0.08: {frac_rho:.3f} (rule <= 0.10)'
    REACHED['tau_b range'] = f'{min(tau.values()):.3f} to {max(tau.values()):.3f}'
    REACHED['NC-2 (primary donors)'] = {d: round(float(NC2[d]), 3) for d in PRIMARY}
    if s7:
        STOPS['S7'] = ('FIRED', REACHED['NC-1'])
    else:
        STOPS['S7'] = ('clear', REACHED['NC-1'])
    write_diag(out)

    # ============================ per-gene calls ============================
    rob_bin = pd.Series(bins, index=rob_list)
    calls = pd.Series(['low_detection' if g in cat_syms else 'not_catalog_gene' for g in kept], index=kept, dtype=object)
    tau_g = pd.Series(np.nan, index=kept)
    for g in rob_list:
        tb = tau[int(rob_bin[g])]; tau_g[g] = tb
        if k_neg[g] == 5 and c_g[g] <= -tb:
            calls[g] = 'central'
        elif k_pos[g] == 5 and c_g[g] >= tb:
            calls[g] = 'portal'
        else:
            calls[g] = 'not_zonated_detected'
    mindet_all = Dd[PRIMARY].min(axis=1)

    def protein_call(sym):
        if sym not in pro.index:
            return 'not_in_table'
        z = pro.loc[sym, 'zonated']; c = float(pro.loc[sym, 'coefficient'])
        if z == 'Yes':
            return 'central' if c < 0 else ('portal' if c > 0 else 'not_zonated')
        if z == 'No':
            return 'not_zonated'
        return 'not_in_table'

    cat = cat.copy()
    cat['zon_mrna_call_yakubovsky2026'] = [calls[s] if s in keptset else 'not_in_assay' for s in cat.symbol]
    cat['zon_protein_call_weiss2026'] = [protein_call(s) for s in cat.symbol]

    # diagnostic 7: detection of T1 set
    for d in PRIMARY:
        v = Dd.loc[T1, d]
        diag(7, 'T1_detection', d, f'min={v.min():.4f} median={v.median():.4f}')
    # diagnostic 8: donor-pair reproducibility over robust non-LM genes
    rnl = sorted(rob_cat - LM)
    pair_rows = []
    for i in range(len(PRIMARY)):
        for j in range(i + 1, len(PRIMARY)):
            r_, n_ = spear(Cdf.loc[rnl, PRIMARY[i]], Cdf.loc[rnl, PRIMARY[j]])
            pair_rows.append((PRIMARY[i], PRIMARY[j], r_, n_))
            diag(8, 'donor_pair_spearman_c', f'{PRIMARY[i]}-{PRIMARY[j]}', f'{r_:.4f} n={n_}')
    pairs = pd.DataFrame(pair_rows, columns=['donor_a', 'donor_b', 'spearman_c', 'n_genes'])
    REACHED['donor-pair Spearman of c over robust non-LM genes (10 pairs)'] = f'median {pairs.spearman_c.median():.3f}, range {pairs.spearman_c.min():.3f} to {pairs.spearman_c.max():.3f}, n {len(rnl):,}'
    # diagnostic 9: calls by class
    tab_m = pd.crosstab(cat.primary_evidence_class, cat.zon_mrna_call_yakubovsky2026, margins=True)
    tab_p = pd.crosstab(cat.primary_evidence_class, cat.zon_protein_call_weiss2026, margins=True)
    for cl, row in tab_m.iterrows():
        diag(9, 'mRNA_call_by_class', cl, json.dumps({k: int(v) for k, v in row.items()}))
    for cl, row in tab_p.iterrows():
        diag(9, 'protein_call_by_class', cl, json.dumps({k: int(v) for k, v in row.items()}))
    tab_m.to_csv(out + '/calls_mrna_by_class.tsv', sep='\t'); tab_p.to_csv(out + '/calls_protein_by_class.tsv', sep='\t')
    write_diag(out)
    log('diagnostics 1-9 written before T1')

    # ============================ T1 ============================
    descr_rows = []
    t1_result = None
    reading = None
    if s7 and not debug:
        log('S7 fired: T1 not computed; columns built, flagged not calibrated')
    else:
        c_T1 = c_g.loc[T1].values
        rho_obs, n_obs = spear(c_T1, coef_T1)
        rho_sp = float(spearmanr(c_T1, coef_T1).statistic)
        rc = rankdata(c_T1)
        rcz = (rc - rc.mean()) / rc.std(); rpz = (rp_T1 - rp_T1.mean()) / rp_T1.std()
        perm_rho = np.empty(B_PERM)
        for i in range(B_PERM):
            perm_rho[i] = (rcz * rpz[rng0.permutation(n_T1)]).mean()
        p_T1 = float(((np.abs(perm_rho) >= abs(rho_obs) - 1e-12).sum() + 1) / (B_PERM + 1))
        lo, hi, boot_rho = boot_ci(rng2, c_T1, coef_T1, B_BOOT)
        rho_d = {d: spear(Cdf.loc[T1, d].values, coef_T1)[0] for d in PRIMARY}
        all_pos = all(rho_d[d] > 0 for d in PRIMARY)
        if p_T1 < ALPHA and rho_obs < 0:
            reading = 'Opposite'
        elif p_T1 >= ALPHA:
            reading = 'Null'
        elif not all_pos:
            reading = 'Not reproduced across donors'
        elif rho_obs >= RHO_FLOOR:
            reading = 'Pass'
        else:
            reading = 'Weak'
        t1_result = dict(rho=rho_obs, rho_scipy=rho_sp, n=n_obs, p=p_T1, ci=(lo, hi), rho_d=rho_d, all_pos=all_pos,
                         B_perm=B_PERM, B_boot=B_BOOT, perm_abs_max=float(np.abs(perm_rho).max()))
        np.savez(out + '/t1_null_and_bootstrap.npz', perm_rho=perm_rho, boot_rho=boot_rho, genes=np.array(T1))
        log('T1', f'rho={rho_obs:.4f} (scipy {rho_sp:.4f}) n={n_obs} p={p_T1:.3e} CI=[{lo:.4f},{hi:.4f}] rho_d={ {k: round(v, 4) for k, v in rho_d.items()} } reading={reading}')
        descr_rows.append(('T1', 'primary M4-M8 median c vs Weiss coefficient', n_obs, rho_obs, lo, hi, p_T1))
        for d in PRIMARY:
            descr_rows.append(('D1', f'rho_d {d}', n_obs, rho_d[d], np.nan, np.nan, np.nan))
        for d in ['M2', 'M3', 'M1']:
            r_, n_ = spear(Cdf.loc[T1, d].values, coef_T1)
            descr_rows.append(('D1', f'rho_d {d} (outside primary set)', n_, r_, np.nan, np.nan, np.nan))

        # stream 3 order: D2, D3, D4 hep, D4 other, D5 zonated, D5 strong, D6 anchored, established, unresolved, untested
        def d_boot(tag, label, genes_):
            genes_ = sorted(genes_)
            a_ = c_g.loc[genes_].values if not tag.startswith('D2') else c_g_sens.loc[genes_].values
            b_ = coef.loc[genes_].values
            r_, n_ = spear(a_, b_)
            lo_, hi_, _ = boot_ci(rng3, a_, b_, B_BOOT)
            descr_rows.append((tag, label, n_, r_, lo_, hi_, np.nan))
            return r_, lo_, hi_, n_

        d_boot('D2', 'sensitivity set M2-M8 median c, same 1,337 genes', T1)
        d_boot('D3', 'T1 genes plus 22 robust landmark genes', set(T1) | set(lm22))
        hep_s = pd.Series(hep[keep_idx], index=kept)
        t1_hep = [g for g in T1 if hep_s[g] == 1]; t1_oth = [g for g in T1 if hep_s[g] != 1]
        d_boot('D4', 'hepatocyte genes (source indicator_hep_gene = 1)', t1_hep)
        d_boot('D4', 'other genes (indicator_hep_gene != 1)', t1_oth)
        d_boot('D5', 'Weiss zonated = Yes', [g for g in T1 if g in pzon])
        d_boot('D5', 'Weiss strong_zonated = Yes', [g for g in T1 if g in pstr])
        clmap = cat.set_index('symbol').primary_evidence_class
        for cl in ['genetically_anchored', 'established_state_associated', 'unresolved', 'untested']:
            d_boot('D6', f'class {cl}', [g for g in T1 if clmap.get(g) == cl])
        conc = [g for g in T1 if clmap.get(g) == 'concordant']
        descr_rows.append(('D6', 'class concordant (counts only; < 10 T1 genes)', len(conc), np.nan, np.nan, np.nan, np.nan))
        # D7 slope (secondary, no interval)
        slope_med = Sdf[PRIMARY].median(axis=1)
        r_, n_ = spear(slope_med.loc[T1].values, coef_T1)
        descr_rows.append(('D7', 'relative OLS slope (zones 1-8), median M4-M8, vs Weiss coefficient', n_, r_, np.nan, np.nan, np.nan))
        r2_, n2_ = spear(slope_med.loc[T1].values, c_T1)
        descr_rows.append(('D7', 'slope median vs COM median (same genes)', n2_, r2_, np.nan, np.nan, np.nan))

    # ============================ per-gene tables ============================
    is_lm = pd.Series([g in LM for g in kept], index=kept)
    sup = pd.DataFrame(index=kept)
    for d in ALL_DONORS:
        sup[f'c_{d}'] = Cdf[d]
    for d in ALL_DONORS:
        sup[f'detection_{d}'] = Dd[d]
    for d in ALL_DONORS:
        sup[f'rel_slope_{d}'] = Sdf[d]
    for d in ALL_DONORS:
        for fld, nm in [('zon_struct/com', 'stored_com'), ('zon_struct_not_hep_norm/com', 'stored_com_not_hep_norm')]:
            arr = STORED[d][fld]
            sup[f'{nm}_{d}'] = arr[keep_idx] if arr is not None else np.nan
    sup['c_median_M4_M8'] = c_g; sup['k_central_M4_M8'] = k_neg; sup['k_portal_M4_M8'] = k_pos; sup['n_donors_defined_M4_M8'] = n_def
    sup['c_median_M2_M8'] = c_g_sens
    sup['k_central_M2_M8'] = (Cdf[SENS] < 0).sum(axis=1); sup['k_portal_M2_M8'] = (Cdf[SENS] > 0).sum(axis=1)
    sup['rel_slope_median_M4_M8'] = Sdf[PRIMARY].median(axis=1)
    sup['min_detection_M4_M8'] = mindet_all
    sup['robust_M4_M8'] = [g in robust['primary_M4_M8'] for g in kept]
    sup['robust_M2_M8'] = [g in robust['sens_M2_M8'] for g in kept]
    sup['detection_decile_bin'] = pd.Series(rob_bin + 1).reindex(kept)
    sup['tau_b'] = tau_g
    sup['mrna_call'] = calls
    sup['landmark'] = ['central_LM' if g in LM_C else ('portal_LM' if g in LM_P else '') for g in kept]
    sup['indicator_hep_gene'] = hep[keep_idx]
    sup['in_T1_set'] = [g in T1sets['primary_M4_M8'] for g in kept]
    sup['weiss_coefficient'] = coef.reindex(kept)
    sup['weiss_q'] = pro['q_value'].reindex(kept)
    sup['weiss_zonated'] = pro['zonated'].reindex(kept)
    sup['weiss_strong_zonated'] = pro['strong_zonated'].reindex(kept)
    clmap = cat.set_index('symbol').primary_evidence_class
    sup['catalog_r8_class'] = clmap.reindex(kept)
    sup.index.name = 'visium_symbol'
    sup.to_csv(out + '/per_gene_supplementary.tsv.gz', sep='\t', float_format='%.6g')
    pd.DataFrame({'visium_symbol': dup_syms, 'note': 'duplicated in the Visium gene list; dropped entirely',
                  'catalog_symbol': [s in set(cat.symbol) for s in dup_syms]}).to_csv(
        out + '/duplicated_visium_symbols.tsv', sep='\t', index=False)

    # Catalog columns
    col = pd.DataFrame({'ensembl_id': cat.ensembl_id.values, 'symbol': cat.symbol.values})
    inv = cat.symbol.isin(keptset).values
    sym = cat.symbol.values
    col['zon_mrna_call_yakubovsky2026'] = cat['zon_mrna_call_yakubovsky2026'].values
    col['zon_mrna_is_landmark'] = [s in LM for s in sym]
    col['zon_mrna_com_median'] = [c_g[s] if ok else np.nan for s, ok in zip(sym, inv)]
    col['zon_mrna_k_central'] = pd.array([int(k_neg[s]) if ok else pd.NA for s, ok in zip(sym, inv)], dtype='Int64')
    col['zon_mrna_k_portal'] = pd.array([int(k_pos[s]) if ok else pd.NA for s, ok in zip(sym, inv)], dtype='Int64')
    col['zon_mrna_min_detection'] = [mindet_all[s] if ok else np.nan for s, ok in zip(sym, inv)]
    col['zon_mrna_com_median_M2_M8'] = [c_g_sens[s] if ok else np.nan for s, ok in zip(sym, inv)]
    col['zon_protein_call_weiss2026'] = cat['zon_protein_call_weiss2026'].values
    col['zon_protein_coef_weiss2026'] = [float(coef[s]) if s in pro.index else np.nan for s in sym]
    col['zon_protein_q_weiss2026'] = [float(pro.loc[s, 'q_value']) if s in pro.index else np.nan for s in sym]
    col['zon_protein_strong_weiss2026'] = [pro.loc[s, 'strong_zonated'] if s in pro.index else None for s in sym]
    if reading == 'Pass':
        mc_ = col.zon_mrna_call_yakubovsky2026.isin(['central', 'portal']).values
        pc_ = col.zon_protein_call_weiss2026.isin(['central', 'portal']).values
        flag = np.where(mc_ & pc_, np.where(col.zon_mrna_call_yakubovsky2026.values == col.zon_protein_call_weiss2026.values, 'agree', 'disagree'),
                        np.where(mc_ | pc_, 'one_layer_only', 'neither'))
        col['zon_layers_flag'] = flag
    base_note = ('mRNA: values derived from the Yakubovsky et al. 2026 Visium deposit (Zenodo 17735587), five non-steatotic '
                 'live-donor livers M4-M8, healthy liver only; protein: Weiss et al. 2026 Nat Metab Supplementary Table 3 (N=14 healthy)')
    notes = []
    for s in sym:
        n_ = base_note
        if s in LM:
            n_ += '; zone axis defined partly by this gene'
        if (s7 and not debug):
            n_ += '; not calibrated (NC-1 failed its rule)'
        notes.append(n_)
    col['zon_source_note'] = notes
    col.to_parquet(out + '/h63_zonation_columns.parquet', index=False)
    dd = [
        ('ensembl_id', 'Catalog r8 Ensembl gene id (key)'),
        ('symbol', 'Catalog r8 symbol (join key to Visium and Weiss symbols)'),
        ('zon_mrna_call_yakubovsky2026', 'central / portal / not_zonated_detected / low_detection (detection < 0.05 in at least one of M4-M8; value shown, no call) / not_in_assay (symbol absent from Visium or one of 10 duplicated symbols). central: robust, 5 of 5 donors c < 0 and median c <= -tau_b; portal: 5 of 5 c > 0 and median c >= tau_b; tau_b = max(0.25, NC-1 null 99th percentile of |c| in the gene detection decile)'),
        ('zon_mrna_is_landmark', 'TRUE for the 40 source landmark genes that define the zone axis'),
        ('zon_mrna_com_median', 'median over M4-M8 of the centered center of mass c = COM - 4.5 of the zone-mean profile of spot-normalized background-subtracted counts, zones 1 (central) to 8 (portal); range -3.5 to +3.5, 0 flat; NA if not in Visium'),
        ('zon_mrna_k_central', 'number of M4-M8 donors with c < 0'),
        ('zon_mrna_k_portal', 'number of M4-M8 donors with c > 0'),
        ('zon_mrna_min_detection', 'minimum over M4-M8 of the fraction of zone 1-8 spots with count > 0'),
        ('zon_mrna_com_median_M2_M8', 'median c over M2-M8 (sensitivity set; M2 and M3 sections are steatotic)'),
        ('zon_protein_call_weiss2026', 'central (Weiss zonated = Yes and coefficient < 0) / portal (Yes and > 0) / not_zonated (No) / not_in_table'),
        ('zon_protein_coef_weiss2026', 'Weiss 2026 zonation coefficient over 20 spatial bins (positive = portal), source-owned'),
        ('zon_protein_q_weiss2026', 'Weiss 2026 q value, copied, not re-corrected'),
        ('zon_protein_strong_weiss2026', 'Weiss 2026 strong_zonated (Yes/No), copied'),
        ('zon_layers_flag', 'agree / disagree / one_layer_only / neither; present only if T1 reads Pass'),
        ('zon_source_note', 'source credit and per-gene caveats'),
    ]
    pd.DataFrame([d_ for d_ in dd if d_[0] in col.columns], columns=['column', 'definition']).to_csv(
        out + '/h63_zonation_columns_dictionary.tsv', sep='\t', index=False)

    # D5 2x2 and D8 named genes
    sections = []
    if t1_result is not None:
        sub = sup.loc[T1]
        pcall = pd.Series([protein_call(g) for g in T1], index=T1)
        both = sub.mrna_call.isin(['central', 'portal']) & pcall.isin(['central', 'portal'])
        t22 = pd.crosstab(sub.mrna_call[both].rename('mrna_call'), pcall[both].rename('protein_call'))
        t22.to_csv(out + '/d5_2x2_T1.tsv', sep='\t')
        allc = cat[cat.zon_mrna_call_yakubovsky2026.isin(['central', 'portal']) & cat.zon_protein_call_weiss2026.isin(['central', 'portal'])]
        t22all = pd.crosstab(allc.zon_mrna_call_yakubovsky2026, allc.zon_protein_call_weiss2026)
        t22all.to_csv(out + '/d5_2x2_all_catalog_incl_landmarks.tsv', sep='\t')
        d8 = cat[(cat.symbol == 'CPS1') | ((cat.primary_evidence_class == 'genetically_anchored')
                                           & cat.zon_mrna_call_yakubovsky2026.isin(['central', 'portal'])
                                           & cat.zon_protein_call_weiss2026.isin(['central', 'portal']))]
        d8 = d8[['symbol', 'ensembl_id', 'primary_evidence_class', 'zon_mrna_call_yakubovsky2026', 'zon_protein_call_weiss2026']].copy()
        d8['c_median_M4_M8'] = [c_g.get(s, np.nan) for s in d8.symbol]
        d8['k_central'] = [k_neg.get(s, np.nan) for s in d8.symbol]
        d8['k_portal'] = [k_pos.get(s, np.nan) for s in d8.symbol]
        d8['min_detection'] = [mindet_all.get(s, np.nan) for s in d8.symbol]
        d8['weiss_coefficient'] = [coef.get(s, np.nan) for s in d8.symbol]
        d8['weiss_q'] = [pro['q_value'].get(s, np.nan) for s in d8.symbol]
        d8['landmark'] = [s in LM for s in d8.symbol]
        d8.to_csv(out + '/d8_named_genes.tsv', sep='\t', index=False)
        descr = pd.DataFrame(descr_rows, columns=['output', 'set', 'n_genes', 'spearman_rho', 'ci_lo_2.5', 'ci_hi_97.5', 'p_two_sided'])
        descr.to_csv(out + '/descriptive_outputs.tsv', sep='\t', index=False)

        r = t1_result
        rd_txt = ', '.join(f'{d} {r["rho_d"][d]:+.3f}' for d in PRIMARY)
        neg_donors = [d for d in PRIMARY if not r['rho_d'][d] > 0]
        neg_rho_txt = ', '.join('%+.2f' % r['rho_d'][d] for d in neg_donors)
        txt = {
            'Pass': f'"For each Catalog gene we give its position along the central-portal axis in five non-steatotic live-donor livers (Visium; Yakubovsky et al. 2026) and its hepatocyte protein zonation (Weiss et al. 2026). Across {r["n"]:,} genes measured in both, the two agreed in direction (Spearman rho {r["rho"]:.2f}, 95% interval {r["ci"][0]:.2f} to {r["ci"][1]:.2f}; positive in each of five donors). Weiss et al. reported a moderate RNA-protein zonation correlation against single-nucleus RNA (their Supplementary Fig. 9)."',
            'Weak': f'"Across {r["n"]:,} genes measured in both layers, mRNA position in five non-steatotic live-donor livers (Visium; Yakubovsky et al. 2026) and hepatocyte protein zonation (Weiss et al. 2026) agreed weakly (Spearman rho {r["rho"]:.2f}, 95% interval {r["ci"][0]:.2f} to {r["ci"][1]:.2f}; positive in each of five donors), so the Catalog shows the two layers in separate columns. Weiss et al. made the earlier RNA-protein comparison against single-nucleus RNA (their Supplementary Fig. 9)."',
            'Not reproduced across donors': f'"Pooled over five donors the layers agreed in direction (rho {r["rho"]:.2f}, {r["ci"][0]:.2f} to {r["ci"][1]:.2f}), but not in donor(s) {", ".join(neg_donors)} (rho_d {neg_rho_txt}); the Catalog lists them separately."',
            'Null': f'"mRNA position in live-donor Visium and hepatocyte protein zonation did not agree at gene level (rho {r["rho"]:.2f}, 95% interval {r["ci"][0]:.2f} to {r["ci"][1]:.2f}); the Catalog lists them separately."',
            'Opposite': f'"mRNA position and hepatocyte protein zonation ran in opposite directions at gene level (rho {r["rho"]:.2f}, {r["ci"][0]:.2f} to {r["ci"][1]:.2f})"; D4 is reported beside it.',
        }[reading]
        sections.append(('T1 (family of one test, alpha 0.05, no multiplicity correction)', [
            f'- Spearman rho = {r["rho"]:.4f} over n = {r["n"]:,} genes (T1 set; gene is the unit of the correlation).',
            f'- Two-sided gene-permutation p = {fmtp(r["p"])} ({r["B_perm"]:,} permutations, stream 0; floor {1 / (r["B_perm"] + 1):.4g}; largest |rho| in the null {r["perm_abs_max"]:.4f}). Genes covary, so this p is too small [Inference]; the donor condition and the 0.30 floor govern the reading.',
            f'- 95% gene-bootstrap percentile interval {r["ci"][0]:.4f} to {r["ci"][1]:.4f} ({r["B_boot"]:,} resamples, stream 2; too narrow when genes covary [Inference]).',
            f'- Donor-level rho_d (each donor\'s c against the same coefficients): {rd_txt}. Positive in {sum(r["rho_d"][d] > 0 for d in PRIMARY)} of 5.',
            f'- Effect-size floor for a combined flag: rho >= 0.30 -> {"met" if r["rho"] >= RHO_FLOOR else "not met"}.',
            f'- Reading: **{reading}**.',
            '',
            'Proposed sentence (returned as data, not installed):',
            '',
            txt,
        ]))
        sections.append(('Descriptive outputs (outside the family; none is called significant)', [
            '| Output | Set | n genes | Spearman rho | 95% interval |', '|---|---|---|---|---|'] + [
            f'| {o} | {s_} | {n_:,} | {fmt(r_)} | {fmt(lo_)} to {fmt(hi_)} |' if np.isfinite(lo_) else f'| {o} | {s_} | {n_:,} | {fmt(r_)} | |'
            for o, s_, n_, r_, lo_, hi_, _p in descr_rows if o != 'T1'] + [
            '', 'D4 prediction: rho higher in hepatocyte genes [Inference]. Per-class rho (D6) is not compared between classes.',
            '', 'D5 2 x 2 (T1 genes called in both layers; rows mRNA call, columns protein call):', '',
            '```', t22.to_string(), '```', '',
            '2 x 2 over all Catalog genes called in both layers, landmarks included:', '', '```', t22all.to_string(), '```',
            '', 'D8 named genes (CPS1 and every genetically anchored gene with both layers called), as is:', '',
            '```', d8.drop(columns=['ensembl_id']).to_string(index=False, float_format=lambda v: f'{v:.3g}'), '```']))
    else:
        sections.append(('T1', ['Not computed: S7 fired (NC-1 outside its rule). Columns built without a combined flag and flagged "not calibrated"; no T1 sentence.']))

    sections.insert(0, ('Controls (predicted vs observed)', [
        f'- PC-3 protein table read: CYP2E1 {pc3["CYP2E1"][0]:.2f} (q {pc3["CYP2E1"][1]:.1e}), ASS1 {pc3["ASS1"][0]:+.2f} (q {pc3["ASS1"][1]:.1e}); predicted -1.82 (6.0e-34) and +1.89 (1.4e-74). Pass.',
        '- PC-1 orientation (circular by design), per primary donor median c over 20 central LM < 0, over 20 portal LM > 0, GLUL < 0: ' + '; '.join(f'{d} {v[0]:+.2f}/{v[1]:+.2f}/{v[2]:+.2f}' for d, v in pc1.items()) + '. Predicted pass in 5 of 5; observed pass in 5 of 5.',
        f'- PC-2 cross-layer join: {agree} of {len(lm22)} robust landmark genes have the same sign in c_g and the Weiss coefficient. Rule >= 15; predicted 19. Disagreeing genes: ' + ', '.join(r_[0] for r_ in pc2_rows if not r_[4]) + '.',
        f'- NC-1 zone labels permuted within donor ({B_NC1} permutations, stream 1): mean fraction of robust genes with 5 of 5 donors on one side {frac5:.4f} (predicted 0.04 to 0.09; rule 0.03 to 0.10); pooled 99th percentile of |c_g| {q99:.3f} (predicted < 0.10, auditor 0.10 to 0.40; not a rule); permutations with |rho| >= 0.08: {frac_rho:.3f} (rule <= 0.10; predicted <= 0.05). S7 {"FIRED" if s7 else "clear"}.',
        '- tau_b by detection decile: ' + ', '.join(f'bin {t[0]} {t[5]:.3f} (false-call {t[6]:.4f})' for t in tau_rows) + '. Predicted false-call rate < 0.01 in every bin.',
        '- NC-2 Spearman(spot total, zone) per primary donor: ' + ', '.join(f'{d} {NC2[d]:+.3f}' for d in PRIMARY) + f'. Slope statistic flagged: {nc2_flag}.',
        '- S4 COM reproduction (best of the two stored fields, robust genes): ' + ', '.join(f'{d} {s4[d]:.4f}' for d in ALL_DONORS) + ' (rule >= 0.90 in M4-M8; predicted >= 0.99).',
        '- Donor-pair Spearman of c over robust non-landmark genes: ' + REACHED['donor-pair Spearman of c over robust non-LM genes (10 pairs)'] + '.',
    ]))
    sections.append(('mRNA and protein calls (Catalog r8 genes)', [
        '```', tab_m.to_string(), '```', '', '```', tab_p.to_string(), '```',
        '', 'Per-class counts depend on the r8 class membership and must be recomputed when the corrected 50-study COLOC rerun is adopted.']))
    sections.append(('Outputs', [
        '- `h63_zonation_columns.parquet` (28,233 Catalog r8 genes, keyed on ensembl_id and symbol) and `h63_zonation_columns_dictionary.tsv`. Candidate file; not written into `catalog_reviewed_r8/`.',
        '- `per_gene_supplementary.tsv.gz` (36,581 kept Visium symbols, both donor sets, all 8 donors), `duplicated_visium_symbols.tsv`.',
        '- `diagnostics.tsv` (diagnostics 1-9 in order, written before T1), `pc2_landmarks.tsv`, `nc1_permutations.tsv`, `tau_by_detection_decile.tsv`, `descriptive_outputs.tsv`, `d5_2x2_*.tsv`, `d8_named_genes.tsv`, `t1_null_and_bootstrap.npz`, `calls_*_by_class.tsv`, `input_sha256.tsv`.',
        '- Gene order inside every resampled set is sorted by symbol (implementation detail; the prespec does not fix it).',
    ]))
    if t1_result is not None:
        REACHED['T1'] = f'rho {t1_result["rho"]:.4f}, 95% interval {t1_result["ci"][0]:.4f} to {t1_result["ci"][1]:.4f}, p {fmtp(t1_result["p"])}, n {t1_result["n"]}, rho_d ' + ', '.join(f'{d} {t1_result["rho_d"][d]:+.4f}' for d in PRIMARY)
    if s7 and not debug:
        verdict = 'Not calibrated (S7)'
        detail = 'NC-1 fell outside its rule; T1 not read; both columns built without a combined flag, flagged "not calibrated".'
        write_result(out, verdict, 'S7', detail, sections, debug)
    else:
        detail = f'All stop rules S1-S7 clear{" (S7 ignored in debug)" if (debug and s7) else ""}. Reading computed from T1 p, rho and the five donor-level rho_d.'
        write_result(out, reading, None, detail, sections, debug)
    log('done')


if __name__ == '__main__':
    main()
