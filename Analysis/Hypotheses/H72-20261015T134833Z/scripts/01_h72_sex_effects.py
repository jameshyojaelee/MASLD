#!/usr/bin/env python3
"""H72: female-minus-male per-allele effects at Catalog enzyme signals, clumped ALT/AST/GGT leads
and chrX, on the raw (U/L) and irnt (sex-specific SD) scales of the Neale round 2 sex-stratified
UK Biobank GWAS.

Implements PRESPEC_H72.md (sha256 9a9724c8...; copied into the run directory). Gates and stop
rules run in the prespec's order; every diagnostic is written before any ALT, AST or GGT
female-minus-male statistic is computed.

--synthetic-outcome (debug only): before anything outcome-related is computed, the female and
male betas of ALT/AST/GGT at gene top SNPs, clumped leads and chrX are replaced by simulated
values (equal-effect null built from the both-sexes irnt beta). Controls (PC1, PC2, NC1) use
real data in both modes, as the prespec allows. Synthetic runs write RESULT_H72_SYNTHETIC.md.
"""
import argparse
import bisect
import hashlib
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import norm

R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
H72 = R + '/Analysis/Hypotheses/H72-20261015T134833Z'
RR = R + '/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z'
R4 = RR + '/round4-20261015T112148Z'
LENS = RR + '/scratch/lens-L06'
SEX = R + '/GWAS/MR_Data/PanUKBB/sex_stratified'
URATE = R + '/data/public_r4/neale_30880'
CAT = R + '/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
COLOC = R + '/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv'
PRESPEC = H72 + '/PRESPEC_H72.md'
PRESPEC_SHA = '9a9724c871df63907cb66d656aba94a688a9ae74d78e27a3dae9aad52db07ecf'
FIELD = {'ALT': '30620', 'AST': '30650', 'GGT': '30730'}
TRAITS = ['ALT', 'AST', 'GGT']
SCALES = ['raw', 'irnt']
Z95 = 1.959964
LO, HI = 2.0 / 3.0, 1.5          # [0.667, 1.5] equivalence bound
F4_LO, F4_HI = 0.90, 1.11
B_BOOT = 10000
CHI2_MED = 0.4549                # as written in S7
PAR = [(60001, 2699520), (154931044, 155260560)]

# Prespec "Inputs with sha256" table (31 files). Bytes None where the table gives none.
INPUTS = [
    ('GWAS/MR_Data/PanUKBB/sex_stratified/30620_irnt.gwas.imputed_v3.female.varorder.tsv.bgz', 463838050, '5177fadd85b2ec089dcb32d63fcd8864a6d0d2d73586a9a7bbba57a0b408a183'),
    ('GWAS/MR_Data/PanUKBB/sex_stratified/30620_irnt.gwas.imputed_v3.male.varorder.tsv.bgz', 463741745, '1c29981b11b3d6e1695d6951f2071886edbca46cb48b4e55a3138754ff2b1328'),
    ('GWAS/MR_Data/PanUKBB/sex_stratified/30650_irnt.gwas.imputed_v3.female.varorder.tsv.bgz', 463844849, '44229f3477c9abb9e1b1404d8f9ec2e0c75bcaa8bc671b32ab09f5bf553d34da'),
    ('GWAS/MR_Data/PanUKBB/sex_stratified/30650_irnt.gwas.imputed_v3.male.varorder.tsv.bgz', 463792252, '0c2e12c5274d27b6869e3d3c1c9d61dca47e7548e9ce383f6457d7558267e4a3'),
    ('GWAS/MR_Data/PanUKBB/sex_stratified/30730_irnt.gwas.imputed_v3.female.varorder.tsv.bgz', 463885011, 'a5adaeed5f8f7f95a449d70ec255c02647898224217e2393ba617f6ab409051e'),
    ('GWAS/MR_Data/PanUKBB/sex_stratified/30730_irnt.gwas.imputed_v3.male.varorder.tsv.bgz', 463719162, '7d866ef92ea8a90a3aaf446f12e7da7b5b7de22f2e50910326a25f2899ec7d78'),
    ('GWAS/MR_Data/PanUKBB/sex_stratified/30620_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz', 463864606, 'cdc9a30d97e4e885d577e32e82ae6a836f196b82f4b3a9bd777a13a2ba7b67c9'),
    ('GWAS/MR_Data/PanUKBB/sex_stratified/30650_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz', 463869240, 'e969f0519e0d0d32247414325f2bc9f75b9289238df94a5f219e994d248b4903'),
    ('GWAS/MR_Data/PanUKBB/sex_stratified/30730_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz', 463753002, 'f13733f39b59fb260eeae0b09a108f0bd93076296c8b670c3ff8aa3c74efd928'),
    ('data/public_r4/neale_30880/30880_irnt.gwas.imputed_v3.female.varorder.tsv.bgz', 463824447, '3712928f5baa517da4c3912db8c8ed6f3c9c2e1e0f0a148e53e1a44435e042af'),
    ('data/public_r4/neale_30880/30880_irnt.gwas.imputed_v3.male.varorder.tsv.bgz', 463693106, '34e7c537ab6e4542094b81b78d95c00b79e8940cd06f8ce77140b9629aefc557'),
    ('GWAS/finemapping/data/sumstats/UKBB_ALT_reformatted_hg19.tsv', 342196068, '695a0154f6c6ef0a0c99ecc2182e28b4455aa9b118f644c0203624d61f1b8d5d'),
    ('GWAS/finemapping/data/sumstats/UKBB_AST_reformatted_hg19.tsv', 342272732, 'fb6e302018c0ad1fd8cca0e3397176e7a626da72dc022960f625bbb1ea838d82'),
    ('GWAS/finemapping/data/sumstats/UKBB_GGT_reformatted_hg19.tsv', 342548337, 'c5d1c94c403c76eef8cad5fcd7e92ac2f9a0fbe9e57670a414945ff3df408b6f'),
    ('GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv', 335138620, 'acdfa306f1ee45383758bc2522602fc2828c4060080b7ad2d96134fc4a399887'),
    ('figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet', 468456, 'ffd3c4650c77a899e9ebed7327697398c4b6462cfffcd7f0e635cd0d023db06b'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L06/ukbb_alt_clumped500kb.tsv', 11133, '7a28fe37cf3399f5cd634184f62591b55a7018acd3b22510fa23458627d98fdf'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L06/ukbb_ast_clumped500kb.tsv', 13235, '303d1731aceef54fcd31af952215f94d7b0f5cf16ccd5bb0a2f3356912b0d338'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L06/ukbb_ggt_clumped500kb.tsv', 14655, 'f3c84f90f867865f737e639ec4896f30ac986278a91adb3ae9ce416d212349f3'),
    ('GWAS/finemapping/config/gwas_registry.tsv', 7868, '6823e42f5c1373785e2a23e08a4e8657fc3c2363d7d4259f6978e8ecfe8b9b8c'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG11/positive_control.txt', 173, '9833bf81f4ab613b58cff0b8b8d1d16619f61ac7736a53e61c2eb1cc4527e901'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG11/pc.py', 1408, 'ed3029e76ad5c18b2cc4f9fe4a30baef465b46b28244990cf7475e13c2462619'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72/count_n.py', None, '8febb07466c0e288327f404ee0def917cea35557d240b29ed171c9c44ae112b3'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72/count_n.out', None, 'a4f67841a47501e3643ec29b3c5b0fca4b5434f2163613fe5d1e58e94d95644f'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72/best_enzyme_row.tsv', None, '056ac0d8f20a5fb09f910f1f024016459592752cc227c49ced8f99e6e4883881'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72/mde.py', None, 'f4c7c3ff1c348d9f959138221f0aba639547bd5b294f4730a858a4376776403f'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72/mde.out', None, '17469e26238d34b8d0c90d284c90aa23ae079f39d037203dd40acd9db6f8b904'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72/mde_rows.tsv', None, 'a0305a800aa7b2c8b03fe6801de35b4951cdc0e2f9308b4e008841461407cad8'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72/count_snps.py', None, 'c3ad673465be259fedb53e1212ed6d84d32883febc6bccd18ebf0aa16447e692'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72/count_snps.out', None, '6f241657793fc88ca6b0da4ad0a131416033a51e303c91d5c79efa25b25319e8'),
    ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72/raw_http.txt', None, 'c21b583935532b08ea951812cba549504b6ece2a2a4c8d9fa234a388c42d1cd9'),
]
RAW_LEN = {'30620_raw.female': 454804739, '30620_raw.male': 454799513, '30650_raw.female': 454497941,
           '30650_raw.male': 454646054, '30730_raw.female': 451369692, '30730_raw.male': 452520192,
           '30880_raw.female': 453341535, '30880_raw.male': 454151669}
# S4 targets
S4 = dict(genes=426, applicable=303, testable_genes=301, f1_tests=281, leads=974, clusters=236)
F1_BY_TRAIT = {'ALT': 87, 'AST': 83, 'GGT': 111}
LEADS_BY_TRAIT = {'ALT': 279, 'AST': 327, 'GGT': 368}
# PC1 reproduction targets (gate positive_control.txt)
PC1 = dict(lead='4:9922167:C:T', bF=-0.47549, bM=-0.23874, zdiff=-43.40)
# chrX coding table (prespec), irnt
CHRX_TABLE = {'ALT': (0.9996, 0.7666, 1.0841, 0.7071), 'AST': (0.9996, 0.7877, 1.1139, 0.7072),
              'GGT': (0.9996, 0.7778, 1.0999, 0.7072)}

COLS = ['variant', 'minor_AF', 'low_confidence_variant', 'n_complete_samples', 'beta', 'se']
DT = {'variant': str, 'minor_AF': float, 'low_confidence_variant': str, 'n_complete_samples': float,
      'beta': float, 'se': float}
CHUNK = 1_000_000

LOGF = None
DIAG = []
CTRL = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if LOGF is not None:
        LOGF.write(s + '\n'); LOGF.flush()


def diag(step, item, value, expected='', status=''):
    DIAG.append(dict(step=step, item=item, value=value, expected=expected, status=status))
    log(f'[diag {step}] {item}: {value}' + (f' (expected {expected})' if expected != '' else '') + (f' -> {status}' if status else ''))


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return path, os.path.getsize(path), h.hexdigest()


# ------------------------------------------------------------------ statistics
def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    idx = np.where(np.isfinite(p))[0]
    m = len(idx)
    if m == 0:
        return q
    ps = p[idx]
    o = np.argsort(ps, kind='mergesort')
    ranked = ps[o] * m / np.arange(1, m + 1)
    qq = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    out = np.empty(m); out[o] = qq
    q[idx] = out
    return q


def holm(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    idx = np.where(np.isfinite(p))[0]
    m = len(idx)
    if m == 0:
        return q
    ps = p[idx]
    o = np.argsort(ps, kind='mergesort')
    adj = np.maximum.accumulate((m - np.arange(m)) * ps[o])
    adj = np.minimum(adj, 1.0)
    out = np.empty(m); out[o] = adj
    q[idx] = out
    return q


def fieller(bF, sF, bM, sM, z=Z95):
    """95% Fieller interval for bF/bM (independent). Bounded iff |bM/sM| > z."""
    bF, sF, bM, sM = map(lambda x: np.asarray(x, float), (bF, sF, bM, sM))
    a = bM ** 2 - z ** 2 * sM ** 2
    b = -2 * bF * bM
    c = bF ** 2 - z ** 2 * sF ** 2
    bounded = np.abs(bM / sM) > z
    disc = np.where(bounded, b * b - 4 * a * c, np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        r1 = (-b - np.sqrt(disc)) / (2 * a)
        r2 = (-b + np.sqrt(disc)) / (2 * a)
    lo = np.where(bounded, np.minimum(r1, r2), -np.inf)
    hi = np.where(bounded, np.maximum(r1, r2), np.inf)
    return lo, hi, bounded


def direction(bF, bM, k=1.0):
    bF = np.asarray(bF, float); bM = k * np.asarray(bM, float)
    out = np.where(np.sign(bF) != np.sign(bM), 'opposite_sign',
                   np.where(np.abs(bF) > np.abs(bM), 'female_larger',
                            np.where(np.abs(bF) < np.abs(bM), 'male_larger', 'equal')))
    out = np.where(np.isfinite(bF) & np.isfinite(bM), out, '')
    return out


def contrast_stats(bF, sF, bM, sM):
    """z_diff, p, Fieller interval, TOST p_eq and direction for autosomal tests."""
    bF, sF, bM, sM = map(lambda x: np.asarray(x, float), (bF, sF, bM, sM))
    z = (bF - bM) / np.sqrt(sF ** 2 + sM ** 2)
    p = 2 * norm.sf(np.abs(z))
    lo, hi, bounded = fieller(bF, sF, bM, sM)
    s = np.sign(bF / sF ** 2 + bM / sM ** 2)
    oF, oM = s * bF, s * bM
    z_lo = (oF - LO * oM) / np.sqrt(sF ** 2 + LO ** 2 * sM ** 2)
    z_hi = (HI * oM - oF) / np.sqrt(sF ** 2 + HI ** 2 * sM ** 2)
    p_eq = np.maximum(norm.sf(z_lo), norm.sf(z_hi))
    inside = bounded & (lo >= LO) & (hi <= HI)
    return dict(z_diff=z, p=p, fieller_lo=lo, fieller_hi=hi, fieller_bounded=bounded,
                fieller_inside=inside, p_eq=p_eq, direction=direction(bF, bM), orient=s)


def assign_states(df):
    """Prespec rules 3-7 for tested rows. df has q_raw, q_irnt, q_eq_raw, q_eq_irnt, dir_*, inside_*."""
    st = []
    for r in df.itertuples():
        sig_r = np.isfinite(r.q_raw) and r.q_raw < 0.05
        sig_i = np.isfinite(r.q_irnt) and r.q_irnt < 0.05
        if sig_r and sig_i:
            st.append('heterogeneous' if r.direction_raw == r.direction_irnt else 'scale_dependent')
        elif sig_r != sig_i:
            st.append('one_scale_only')
        elif (r.fieller_inside_raw and r.fieller_inside_irnt and np.isfinite(r.q_eq_raw) and np.isfinite(r.q_eq_irnt)
              and r.q_eq_raw < 0.05 and r.q_eq_irnt < 0.05):
            st.append('homogeneous_within_bound')
        else:
            st.append('underpowered')
    return st


# ------------------------------------------------------------------ Neale scans
def read_chunks(path):
    return pd.read_csv(path, sep='\t', compression='gzip', usecols=COLS, dtype=DT, chunksize=CHUNK,
                       na_values=['NaN', 'nan', ''], keep_default_na=False)


def first_n(series):
    s = series.dropna()
    return float(s.iloc[0]) if len(s) else None


def scan_trait(job):
    """Stream female/male irnt and raw files of one trait in lockstep (S3) and keep only what H72 needs."""
    t = job['trait']; need = job['need']; raw_ok = job['raw_ok']; synthetic = job['synthetic']
    f = FIELD[t]
    paths = {('irnt', 'F'): f'{SEX}/{f}_irnt.gwas.imputed_v3.female.varorder.tsv.bgz',
             ('irnt', 'M'): f'{SEX}/{f}_irnt.gwas.imputed_v3.male.varorder.tsv.bgz'}
    if raw_ok:
        paths[('raw', 'F')] = f'{job["downloads"]}/{f}_raw.gwas.imputed_v3.female.varorder.tsv.bgz'
        paths[('raw', 'M')] = f'{job["downloads"]}/{f}_raw.gwas.imputed_v3.male.varorder.tsv.bgz'
    if synthetic:
        paths[('both', 'B')] = f'{SEX}/{f}_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz'
    keys = list(paths)
    readers = {k: read_chunks(paths[k]) for k in keys}
    res = dict(trait=t, s3_ok=True, s3_msg='', rows=0, n={}, keep=[], chrx=[], pool=[], raw_ok=raw_ok)
    sdv = {(sc, sx): [] for sc in SCALES for sx in 'FM'}
    ase = {sc: [] for sc in SCALES}
    off = 0; elig = 0
    t0 = time.time()
    it = zip(*[readers[k] for k in keys])
    for chunks in it:
        ch = dict(zip(keys, chunks))
        ref = ch[keys[0]]
        nrow = len(ref)
        v0 = ref['variant'].values
        for k in keys[1:]:
            if len(ch[k]) != nrow or not np.array_equal(ch[k]['variant'].values, v0):
                bad = np.where(ch[k]['variant'].values[:nrow] != v0[:len(ch[k])])[0]
                i = int(bad[0]) if len(bad) else min(nrow, len(ch[k]))
                res['s3_ok'] = False
                res['s3_msg'] = f'{k} differs from {keys[0]} at row {off + i}'
                return res
        for k in keys:
            if k not in res['n'] or res['n'][k] is None:
                res['n'][k] = first_n(ch[k]['n_complete_samples'])
        sp = ref['variant'].str.split(':', n=3, expand=True)
        chrom = sp[0].values; posS = sp[1].values; refa = sp[2].values; alta = sp[3].values
        pos = posS.astype(np.int64)
        isX = chrom == 'X'
        gidx = off + np.arange(nrow)
        low = {k: (ch[k]['low_confidence_variant'].values == 'true') for k in keys}
        maf = {k: ch[k]['minor_AF'].values for k in keys}
        se = {k: ch[k]['se'].values for k in keys}
        nn = {k: ch[k]['n_complete_samples'].values for k in keys}
        # implied SD (PC3, S8) and autosomal se ratio (S9): every 50th global row, autosomes
        s50 = (gidx % 50 == 0) & ~isX
        with np.errstate(invalid='ignore'):
            for sc in SCALES:
                if (sc, 'F') not in ch:
                    continue
                for sx in 'FM':
                    k = (sc, sx)
                    m = s50 & ~low[k] & (maf[k] >= 0.05) & (maf[k] <= 0.5) & np.isfinite(se[k]) & (se[k] != 0) \
                        & np.isfinite(nn[k]) & (nn[k] != 0)
                    sdv[k].extend((se[k][m] * np.sqrt(2 * maf[k][m] * (1 - maf[k][m]) * nn[k][m])).tolist())
                kF, kM = (sc, 'F'), (sc, 'M')
                m = s50 & ~low[kF] & (maf[kF] >= 0.05) & (maf[kF] <= 0.5) & np.isfinite(se[kF]) & (se[kF] != 0) \
                    & np.isfinite(nn[kF]) & (nn[kF] != 0) & ~low[kM] & np.isfinite(se[kM]) & (se[kM] != 0)
                ase[sc].extend((se[kM][m] / se[kF][m]).tolist())
        # rows to keep: needed positions, chrX, NC1 pool
        key = pd.Series(chrom).str.cat(pd.Series(posS), sep=':').values
        mneed = pd.Series(key).isin(need).values
        snv = (pd.Series(refa).str.len().values == 1) & (pd.Series(alta).str.len().values == 1)
        kF, kM = ('irnt', 'F'), ('irnt', 'M')
        e = ~isX & snv & ~low[kF] & ~low[kM] & (maf[kF] >= 0.05) & (maf[kM] >= 0.05)
        cum = elig + np.cumsum(e)
        mpool = e & (cum % 20 == 0)
        elig = int(cum[-1]) if nrow else elig
        sel = mneed | isX | mpool
        if sel.any():
            d = pd.DataFrame({'variant': v0[sel], 'chr': chrom[sel], 'pos': pos[sel], 'ref': refa[sel], 'alt': alta[sel],
                              'key': key[sel], 'is_need': mneed[sel], 'is_x': isX[sel], 'is_pool': mpool[sel], 'snv': snv[sel]})
            for k in keys:
                tag = {'irnt': 'i', 'raw': 'r', 'both': 'b'}[k[0]] + k[1]
                d['low_' + tag] = low[k][sel]
                d['maf_' + tag] = maf[k][sel]
                d['b_' + tag] = ch[k]['beta'].values[sel]
                d['s_' + tag] = se[k][sel]
            res['keep'].append(d)
        off += nrow
    # every reader must be exhausted at the same row (zip stops at the shortest)
    for k in keys:
        try:
            nxt = next(readers[k])
            res['s3_ok'] = False
            res['s3_msg'] = f'{k} has more rows than the shortest file ({off} rows read; next chunk {len(nxt)} rows)'
            return res
        except StopIteration:
            pass
    res['rows'] = off
    res['elig'] = elig
    d = pd.concat(res['keep'], ignore_index=True)
    res['keep'] = None
    if not raw_ok:
        for tag in ('rF', 'rM'):
            d['low_' + tag] = d['low_i' + tag[1]]; d['maf_' + tag] = d['maf_i' + tag[1]]
            d['b_' + tag] = np.nan; d['s_' + tag] = np.nan
    if synthetic:
        # Equal-effect null at outcome rows (needed positions and chrX); pool rows keep real betas (NC1).
        rng = np.random.default_rng(777 + TRAITS.index(t))
        m = (d.is_need | d.is_x).values
        k_m = np.where(d.is_x.values, 0.5, 1.0)
        bB = d.b_bB.values
        bF = bB + rng.standard_normal(len(d)) * d.s_iF.values
        bM = k_m * bB + rng.standard_normal(len(d)) * d.s_iM.values
        d.loc[m, 'b_iF'] = bF[m]; d.loc[m, 'b_iM'] = bM[m]
        d.loc[m, 'b_rF'] = (bF * d.s_rF.values / d.s_iF.values)[m]
        d.loc[m, 'b_rM'] = (bM * d.s_rM.values / d.s_iM.values)[m]
    res['data'] = d
    med = lambda x: float(np.median(x)) if len(x) else float('nan')
    res['sd'] = {f'{sc}_{sx}': med(sdv[(sc, sx)]) for sc in SCALES for sx in 'FM'}
    res['sd_n'] = {f'{sc}_{sx}': len(sdv[(sc, sx)]) for sc in SCALES for sx in 'FM'}
    res['a_se_ratio'] = {sc: med(ase[sc]) for sc in SCALES}
    res['secs'] = round(time.time() - t0, 1)
    res['n'] = {f'{k[0]}_{k[1]}': v for k, v in res['n'].items()}
    return res


def scan_urate(job):
    """PC1 (irnt) and PC2 (raw) at the SLC2A9 region; files read in lockstep until chr4 ends."""
    raw_ok = job['raw_ok']
    paths = {('irnt', 'F'): f'{URATE}/30880_irnt.gwas.imputed_v3.female.varorder.tsv.bgz',
             ('irnt', 'M'): f'{URATE}/30880_irnt.gwas.imputed_v3.male.varorder.tsv.bgz'}
    if raw_ok:
        paths[('raw', 'F')] = f'{job["downloads"]}/30880_raw.gwas.imputed_v3.female.varorder.tsv.bgz'
        paths[('raw', 'M')] = f'{job["downloads"]}/30880_raw.gwas.imputed_v3.male.varorder.tsv.bgz'
    keys = list(paths)
    readers = {k: read_chunks(paths[k]) for k in keys}
    out = dict(n={}, region=None, identity_ok=True)
    parts = []; seen4 = False
    for chunks in zip(*[readers[k] for k in keys]):
        ch = dict(zip(keys, chunks))
        for k in keys:
            if out['n'].get(k) is None:
                out['n'][k] = first_n(ch[k]['n_complete_samples'])
        ref = ch[keys[0]]
        sp = ref['variant'].str.split(':', n=2, expand=True)
        c = sp[0].values; pos = sp[1].values.astype(np.int64)
        m = (c == '4') & (pos >= 9000000) & (pos <= 10500000)
        if m.any():
            seen4 = True
            d = pd.DataFrame({'variant': ref['variant'].values[m], 'pos': pos[m]})
            for k in keys:
                tag = {'irnt': 'i', 'raw': 'r'}[k[0]] + k[1]
                if len(ch[k]) != len(ref) or not np.array_equal(ch[k]['variant'].values[m], d.variant.values):
                    out['identity_ok'] = False
                d['low_' + tag] = ch[k]['low_confidence_variant'].values[m]
                d['b_' + tag] = ch[k]['beta'].values[m]
                d['s_' + tag] = ch[k]['se'].values[m]
            parts.append(d)
        # the region lies early on chr4; stop once a chunk ends past it
        if seen4 and (c[-1] != '4' or pos[-1] > 10500000):
            break
    out['n'] = {f'{k[0]}_{k[1]}': v for k, v in out['n'].items()}
    out['region'] = pd.concat(parts, ignore_index=True) if parts else None
    return out


def read_source_alleles(job):
    t, need = job
    al = {}
    with open(f'{R}/GWAS/finemapping/data/sumstats/UKBB_{t}_reformatted_hg19.tsv') as fh:
        next(fh)
        for line in fh:
            a = line.split('\t', 4)
            k = a[0] + ':' + a[1]
            if k in need:
                al.setdefault(k, set()).add(frozenset((a[2], a[3])))
    return t, al


# ------------------------------------------------------------------ main
def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--downloads', required=True)
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--synthetic-outcome', action='store_true')
    a = ap.parse_args()
    OUT = a.out
    SYN = a.synthetic_outcome
    os.makedirs(OUT, exist_ok=True)
    LOGF = open(OUT + '/run.log', 'a')
    log(f'H72 start {time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}; out {OUT}; synthetic_outcome {SYN}; seed {a.seed}')
    if a.seed != 20260923:
        log('seed differs from the prespec (20260923); refusing to run'); sys.exit(2)
    stop = []        # global stop reasons (S1, S4, S5)
    result_name = 'RESULT_H72_SYNTHETIC.md' if SYN else 'RESULT_H72.md'

    def finish_stopped(reason_lines, verdict):
        write_tables(OUT)
        with open(OUT + '/' + result_name, 'w') as fh:
            fh.write('# H72 result\n\n')
            if SYN:
                fh.write('SYNTHETIC OUTCOME DEBUG RUN: female and male enzyme betas were simulated. Not a result.\n\n')
            fh.write(f'Verdict (computed by the prespec rule): {verdict}\n\n')
            for l in reason_lines:
                fh.write(f'- {l}\n')
            fh.write('\nLimitations sentence (prespec, untestable branch): sex heterogeneity of the genetic anchors was not tested.\n')
        log('VERDICT', verdict)
        sys.exit(0)

    # ---------------------------------------------------------- 1. S1 hashes
    with Pool(8) as pool:
        hs = pool.map(sha256, [PRESPEC] + [R + '/' + p for p, _, _ in INPUTS])
    ph = hs[0][2]
    diag(1, 'prespec sha256', ph, PRESPEC_SHA, 'pass' if ph == PRESPEC_SHA else 'FAIL S1')
    if ph != PRESPEC_SHA:
        stop.append('S1 prespec hash mismatch')
    for (path, nbytes, h), (rel, b0, h0) in zip(hs[1:], INPUTS):
        ok = h == h0 and (b0 is None or b0 == nbytes)
        diag(1, 'input ' + rel, f'{nbytes} bytes sha256 {h}', f'{b0} {h0}', 'pass' if ok else 'FAIL S1')
        if not ok:
            stop.append(f'S1 input mismatch: {rel}')
    if stop:
        finish_stopped(stop, 'stopped_by_rule (S1)')

    # raw downloads (S2, S6): status written by the sbatch download block
    st_path = OUT + '/download_status.tsv'
    raw_ok_trait = {t: True for t in TRAITS}
    urate_raw_ok = True
    dst = pd.read_csv(st_path, sep='\t', header=None, names=['file', 'bytes', 'ok'], dtype=str)
    shas = {}
    if os.path.exists(OUT + '/download_sha256.txt'):
        for line in open(OUT + '/download_sha256.txt'):
            h, p = line.split()
            shas[os.path.basename(p)] = h
    for r in dst.itertuples():
        code = r.file.split('_')[0]; sx = 'female' if '.female.' in r.file else 'male'
        exp = RAW_LEN[f'{code}_raw.{sx}']
        ok = r.ok == '1' and str(exp) == str(r.bytes)
        diag(1, 'raw download ' + r.file, f'{r.bytes} bytes ok={r.ok} sha256 {shas.get(r.file, "NA")}', f'{exp} bytes, gzip -t pass',
             'pass' if ok else ('FAIL S6' if code == '30880' else 'FAIL S2'))
        if not ok:
            if code == '30880':
                urate_raw_ok = False
            else:
                raw_ok_trait = {t: False for t in TRAITS}   # S2: the raw arm is untestable
    if len(dst) != 8:
        diag(1, 'raw download rows', len(dst), 8, 'FAIL S2')
        raw_ok_trait = {t: False for t in TRAITS}; urate_raw_ok = False

    # ---------------------------------------------------------- gene set, best enzyme row (before scanning)
    cat = pd.read_parquet(CAT)
    anc = cat[cat.primary_evidence_class.isin(['genetically_anchored', 'concordant'])].copy()
    co = pd.read_csv(COLOC, usecols=['gwas_name', 'gene', 'ensembl', 'PP.H4.abf', 'PP.H4.susie', 'method', 'top_snp', 'top_snp_PP'],
                     dtype=str)
    co = co[co.gwas_name.isin(['UKBB_ALT', 'UKBB_AST', 'UKBB_GGT']) & co.ensembl.isin(set(anc.ensembl_id))].copy()
    n_enz_rows = len(co); n_dup = int(co.duplicated(['gwas_name', 'ensembl']).sum())
    co['pp_s'] = pd.to_numeric(co['PP.H4.susie'], errors='coerce')
    co['pp_a'] = pd.to_numeric(co['PP.H4.abf'], errors='coerce')
    co['pp'] = np.where(np.isfinite(co.pp_s), co.pp_s, co.pp_a)
    co['pp_src'] = np.where(np.isfinite(co.pp_s), 'susie', 'abf')
    co = co[np.isfinite(co.pp)].sort_values(['ensembl', 'pp', 'gwas_name'], ascending=[True, False, True])
    best = co.drop_duplicates('ensembl').copy()
    best['trait'] = best.gwas_name.str.replace('UKBB_', '', regex=False)
    best['chr'] = best.top_snp.str.split(':').str[0]
    best['pos'] = best.top_snp.str.split(':').str[1].astype(int)
    best = anc[['ensembl_id', 'symbol', 'primary_evidence_class']].merge(best, left_on='ensembl_id', right_on='ensembl', how='left')
    best['applicable'] = best.pp > 0.5
    n_pp_exact = int((best.pp == 0.5).sum())

    leads = {}
    for t in TRAITS:
        L = pd.read_csv(f'{LENS}/ukbb_{t.lower()}_clumped500kb.tsv', sep='\t', header=None,
                        names=['chr', 'pos', 'a1', 'a2', 'beta_both', 'se_both', 'p_both'], dtype={'chr': str})
        L['trait'] = t
        L['key'] = L.chr + ':' + L.pos.astype(str)
        leads[t] = L
    need = {t: set(best[best.trait == t].top_snp) | set(leads[t].key) for t in TRAITS}

    # source alleles at gene top SNPs
    with Pool(3) as pool:
        src = dict(pool.map(read_source_alleles, [(t, set(best[best.trait == t].top_snp)) for t in TRAITS]))

    # ---------------------------------------------------------- Neale scans
    jobs = [dict(trait=t, need=need[t], raw_ok=raw_ok_trait[t], synthetic=SYN, downloads=a.downloads) for t in TRAITS]
    with Pool(4) as pool:
        ur_async = pool.apply_async(scan_urate, (dict(raw_ok=urate_raw_ok, downloads=a.downloads),))
        scans = {r['trait']: r for r in pool.map(scan_trait, jobs)}
        ur = ur_async.get()
    trait_ok = {t: True for t in TRAITS}

    # ---------------------------------------------------------- 2. S3 identity
    for t in TRAITS:
        s = scans[t]
        diag(2, f'{t} variant identity across female/male irnt' + ('/raw' if s['raw_ok'] else '') + ' files',
             'identical' if s['s3_ok'] else s['s3_msg'], 'identical row by row', 'pass' if s['s3_ok'] else 'FAIL S3 (trait stopped)')
        if s['s3_ok']:
            diag(2, f'{t} rows per file', s['rows'], 13791467, 'pass' if s['rows'] == 13791467 else 'note')
            log(f'{t} scan seconds {s["secs"]}')
        else:
            trait_ok[t] = False
    diag(2, 'urate region variant identity irnt/raw', ur['identity_ok'], True, 'pass' if ur['identity_ok'] else 'FAIL (PC2 not usable)')
    if not ur['identity_ok']:
        urate_raw_ok = False

    # ---------------------------------------------------------- 3. n
    n_tab = {}
    for t in TRAITS:
        if not trait_ok[t]:
            continue
        n = scans[t]['n']
        n_tab[t] = n
        for k, v in n.items():
            diag(3, f'{t} n_complete_samples {k}', v)
        if scans[t]['raw_ok']:
            for sx in 'FM':
                ni, nr = n.get(f'irnt_{sx}'), n.get(f'raw_{sx}')
                ok = ni is not None and nr is not None and abs(nr - ni) / ni <= 0.01
                diag(3, f'{t} raw vs irnt n {sx}', f'raw {nr} irnt {ni}', 'within 1%', 'pass' if ok else 'FAIL (raw arm of trait untestable)')
                if not ok:
                    raw_ok_trait[t] = False
    for k, v in ur['n'].items():
        diag(3, f'urate n_complete_samples {k}', v)

    # ---------------------------------------------------------- 4. recounts (S4) and allele resolution
    data = {t: scans[t]['data'] for t in TRAITS if trait_ok[t]}
    diag(4, 'anchored+concordant genes', len(anc), S4['genes'], 'pass' if len(anc) == S4['genes'] else 'FAIL S4')
    diag(4, 'classes', anc.primary_evidence_class.value_counts().to_dict(), "{'genetically_anchored': 400, 'concordant': 26}")
    diag(4, 'UKBB enzyme rows / duplicates', f'{n_enz_rows} / {n_dup}', '1278 / 0', 'pass' if (n_enz_rows, n_dup) == (1278, 0) else 'note')
    diag(4, 'genes with best PP exactly 0.5', n_pp_exact, 0)
    n_app = int(best.applicable.sum())
    diag(4, 'applicable genes (best PP.H4 > 0.5)', n_app, S4['applicable'], 'pass' if n_app == S4['applicable'] else 'FAIL S4')

    need_idx = {t: {k: x for k, x in data[t][data[t].is_need].groupby('key')} for t in data}

    def resolve(t, key, alleles):
        """Return (status, row) for one site against the Neale rows of trait t."""
        rows = need_idx[t].get(key)
        if rows is None:
            return 'no_neale_match', None
        m = rows[[frozenset((r, al)) == alleles for r, al in zip(rows.ref, rows.alt)]]
        if len(m) == 0:
            return 'no_neale_match', None
        if len(m) > 1:
            return 'neale_multi_match', None
        r = m.iloc[0]
        if r.low_iF or r.low_iM or not (r.maf_iF >= 0.01) or not (r.maf_iM >= 0.01):
            return 'lowconf_or_maf', None
        return 'ok', r

    gene_rows = []
    for g in best.itertuples():
        rec = dict(ensembl_id=g.ensembl_id, symbol=g.symbol, catalog_class=g.primary_evidence_class,
                   best_enzyme_trait=g.trait, pp_h4=g.pp, pp_source=g.pp_src, top_snp=g.top_snp, neale_variant='',
                   status='')
        if not g.applicable:
            rec['status'] = 'not_applicable'
        elif not trait_ok.get(g.trait, False):
            rec['status'] = 'untestable_S3'
        else:
            s_ = src[g.trait].get(g.top_snp)
            if not s_:
                rec['status'] = 'untestable_no_source_row'
            elif len(s_) > 1:
                rec['status'] = 'untestable_source_multi'
            else:
                stt, r = resolve(g.trait, g.top_snp, next(iter(s_)))
                if stt != 'ok':
                    rec['status'] = 'untestable_' + stt
                else:
                    rec['status'] = 'testable'; rec['neale_variant'] = r.variant
        gene_rows.append(rec)
    genes = pd.DataFrame(gene_rows)
    for t in TRAITS:
        gt = genes[genes.best_enzyme_trait == t]
        diag(4, f'{t} allele resolution (applicable genes)', gt[gt.status != 'not_applicable'].status.value_counts().to_dict())
    for r in genes[genes.status.str.startswith('untestable')].itertuples():
        diag(4, 'untestable gene', f'{r.symbol} {r.ensembl_id} {r.best_enzyme_trait} {r.top_snp} {r.status}')
    n_testable = int((genes.status == 'testable').sum())
    diag(4, 'testable genes', n_testable, S4['testable_genes'], 'pass' if n_testable == S4['testable_genes'] else 'FAIL S4')
    tests = genes[genes.status == 'testable'].drop_duplicates(['best_enzyme_trait', 'top_snp'])[['best_enzyme_trait', 'top_snp', 'neale_variant']]
    tests = tests.rename(columns={'best_enzyme_trait': 'trait'}).reset_index(drop=True)
    diag(4, 'F1 unique trait x SNP tests', len(tests), S4['f1_tests'], 'pass' if len(tests) == S4['f1_tests'] else 'FAIL S4')
    diag(4, 'F1 tests by trait', tests.trait.value_counts().to_dict(), F1_BY_TRAIT)

    # clusters: 500 kb single-linkage per trait and chromosome
    tests['chr'] = tests.top_snp.str.split(':').str[0]
    tests['pos'] = tests.top_snp.str.split(':').str[1].astype(int)
    tests = tests.sort_values(['trait', 'chr', 'pos']).reset_index(drop=True)
    cl = []; cid = -1; prev = None
    for r in tests.itertuples():
        if prev is None or (r.trait, r.chr) != (prev.trait, prev.chr) or r.pos - prev.pos > 500000:
            cid += 1
        cl.append(cid); prev = r
    tests['cluster'] = cl
    n_cl = tests.cluster.nunique()
    diag(4, '500 kb clusters of testable F1 SNPs', n_cl, S4['clusters'], 'pass' if n_cl == S4['clusters'] else 'FAIL S4')
    diag(4, '500 kb clusters by trait', tests.groupby('trait').cluster.nunique().to_dict(), "{'ALT': 72, 'AST': 70, 'GGT': 94}")

    lead_rows = []
    for t in TRAITS:
        for r in leads[t].itertuples():
            rec = dict(trait=t, lead=r.key, a1=r.a1, a2=r.a2, neale_variant='', status='')
            if not trait_ok[t]:
                rec['status'] = 'untestable_S3'
            else:
                stt, x = resolve(t, r.key, frozenset((r.a1, r.a2)))
                rec['status'] = 'testable' if stt == 'ok' else 'untestable_' + stt
                if stt == 'ok':
                    rec['neale_variant'] = x.variant
            lead_rows.append(rec)
    lds = pd.DataFrame(lead_rows)
    n_leads = int((lds.status == 'testable').sum())
    for t in TRAITS:
        diag(4, f'{t} lead allele resolution', lds[lds.trait == t].status.value_counts().to_dict())
    diag(4, 'testable leads (F2)', n_leads, S4['leads'], 'pass' if n_leads == S4['leads'] else 'FAIL S4')
    s4_fail = [d_ for d_ in DIAG if d_['status'] == 'FAIL S4']
    if s4_fail:
        finish_stopped([f'S4 recount differs: {d_["item"]} = {d_["value"]} (expected {d_["expected"]})' for d_ in s4_fail],
                       'stopped_by_rule (S4)')
    diag(4, 'G1 coverage (>= 200 testable F1 tests)', len(tests), '>= 200', 'pass' if len(tests) >= 200 else 'FAIL G1')
    if len(tests) < 200:
        finish_stopped(['G1 coverage below 200 F1 tests'], 'stopped_by_rule (G1)')

    # attach Neale values to tests and leads
    def attach(df, tcol='trait'):
        parts = []
        for t in TRAITS:
            x = df[df[tcol] == t]
            if len(x) == 0:
                continue
            d = data[t][data[t].is_need].drop_duplicates('variant').set_index('variant')
            cols = ['maf_iF', 'maf_iM', 'low_rF', 'low_rM', 'maf_rF', 'maf_rM', 'b_iF', 's_iF', 'b_iM', 's_iM', 'b_rF', 's_rF', 'b_rM', 's_rM']
            parts.append(x.join(d[cols], on='neale_variant'))
        return pd.concat(parts).loc[df.index]
    tests = attach(tests)
    lds_t = attach(lds[lds.status == 'testable'].copy())
    for nm, df in [('F1 tests', tests), ('F2 leads', lds_t)]:
        mism = int(((df.low_rF != False) | (df.low_rM != False) | (df.maf_rF != df.maf_iF) | (df.maf_rM != df.maf_iM)).sum())
        diag(4, f'{nm}: raw-file low_conf/minor_AF differing from irnt file', mism, 0, 'pass' if mism == 0 else 'note')

    # ---------------------------------------------------------- 5. implied SD (S8, PC3)
    sd = {}
    relabel_s8 = {}
    for t in TRAITS:
        if not trait_ok[t]:
            continue
        sd[t] = scans[t]['sd']
        for k, v in sd[t].items():
            diag(5, f'{t} implied trait SD {k} (n rows {scans[t]["sd_n"][k]})', round(v, 4) if np.isfinite(v) else v,
                 '0.90-1.10' if k.startswith('irnt') else '')
        bad = [k for k in ('irnt_F', 'irnt_M') if not (0.90 <= sd[t][k] <= 1.10)]
        relabel_s8[t] = bool(bad)
        diag(5, f'{t} S8 irnt SD check', 'outside 0.90-1.10: ' + ','.join(bad) if bad else 'within 0.90-1.10', '',
             'S8 relabel to pooled SD units' if bad else 'pass')
        if raw_ok_trait[t]:
            ratio = sd[t]['raw_F'] / sd[t]['raw_M']
            sd[t]['ratio_raw'] = ratio
            pred = '< 1' if t in ('ALT', 'GGT') else ''
            diag(5, f'{t} PC3 SD_F,raw / SD_M,raw', round(ratio, 4), pred,
                 ('met' if ratio < 1 else 'not met') if pred else 'reported')
            CTRL.append(dict(control='PC3', trait=t, scale='raw', statistic='SD_F/SD_M', value=ratio, prediction=pred or 'none',
                             outcome=('met' if ratio < 1 else 'not met') if pred else 'reported'))
        else:
            sd[t]['ratio_raw'] = float('nan')

    # ---------------------------------------------------------- 6. chrX coding (S9)
    chrx_ok = {t: trait_ok[t] for t in TRAITS}
    for t in TRAITS:
        if not trait_ok[t]:
            continue
        d = data[t][data[t].is_x]
        npar = int(sum(((d.pos >= lo) & (d.pos <= hi)).sum() for lo, hi in PAR))
        diag(6, f'{t} chrX rows / PAR rows', f'{len(d)} / {npar}', '427164 / 0')
        for sc in SCALES:
            if sc == 'raw' and not raw_ok_trait[t]:
                continue
            tg = 'i' if sc == 'irnt' else 'r'
            m = (~d['low_' + tg + 'F']) & (~d['low_' + tg + 'M']) & (d['maf_' + tg + 'F'] >= 0.05)
            x = d[m]
            sF, sM, mM, mF = x['s_' + tg + 'F'], x['s_' + tg + 'M'], x['maf_' + tg + 'M'], x['maf_' + tg + 'F']
            ok_ = np.isfinite(sF) & (sF != 0) & np.isfinite(sM) & (sM != 0) & np.isfinite(mM) & (mM != 0)
            maf_ratio = float(np.median(mM[ok_] / mF[ok_]))
            se_x = float(np.median(sM[ok_] / sF[ok_]))
            se_a = scans[t]['a_se_ratio'][sc]
            rr = se_x / se_a
            ok = 0.95 <= maf_ratio <= 1.05 and 0.68 <= rr <= 0.74
            exp = CHRX_TABLE[t] if sc == 'irnt' else '(raw: no table; same S9 bands)'
            diag(6, f'{t} chrX coding {sc}: n {int(ok_.sum())}; male/female MAF; se_M/se_F chrX; autosomes; chrX/autosome',
                 f'{maf_ratio:.4f}; {se_x:.4f}; {se_a:.4f}; {rr:.4f}', exp, 'pass' if ok else 'FAIL S9 (chrX untestable)')
            if not ok:
                chrx_ok[t] = False

    # ---------------------------------------------------------- 7. PC1, PC2
    reg = ur['region']
    x = reg[(reg.low_iF == 'false') & (reg.low_iM == 'false') & np.isfinite(reg.b_iF) & np.isfinite(reg.b_iM)].copy()
    x['zc'] = (x.b_iF / x.s_iF ** 2 + x.b_iM / x.s_iM ** 2) / np.sqrt(1 / x.s_iF ** 2 + 1 / x.s_iM ** 2)
    i = int(np.argmax(np.abs(x.zc.values)))       # first maximum, as the gate's strict '>' scan
    L1 = x.iloc[i]
    zd1 = (L1.b_iF - L1.b_iM) / math.sqrt(L1.s_iF ** 2 + L1.s_iM ** 2)
    sig4 = lambda v: float(f'{v:.4g}')
    pc1_ok = (L1.variant == PC1['lead'] and sig4(L1.b_iF) == sig4(PC1['bF']) and sig4(L1.b_iM) == sig4(PC1['bM'])
              and abs(zd1 - PC1['zdiff']) <= 0.05)
    diag(7, 'PC1 urate irnt lead; beta_F; beta_M; z_diff', f'{L1.variant}; {L1.b_iF}; {L1.b_iM}; {zd1:.2f}',
         f'{PC1["lead"]}; {PC1["bF"]}; {PC1["bM"]}; {PC1["zdiff"]} +/- 0.05', 'pass' if pc1_ok else 'FAIL S5')
    CTRL.append(dict(control='PC1', trait='urate', scale='irnt', statistic='z_diff', value=zd1,
                     prediction='4:9922167:C:T, z_diff -43.40 +/- 0.05, betas to 4 sig. digits', outcome='pass' if pc1_ok else 'FAIL S5'))
    if not pc1_ok:
        finish_stopped([f'S5: PC1 not reproduced ({L1.variant}, z_diff {zd1:.2f})'], 'stopped_by_rule (S5); untestable, no Catalog column')
    pc2_ok = False
    if urate_raw_ok:
        y = reg[reg.variant == L1.variant].iloc[0]
        zd2 = (y.b_rF - y.b_rM) / math.sqrt(y.s_rF ** 2 + y.s_rM ** 2)
        ratio2 = y.b_rF / y.b_rM
        pc2_ok = abs(zd2) >= 4 and abs(y.b_rF) > abs(y.b_rM)
        pred_met = zd2 < -10 and abs(y.b_rF) > abs(y.b_rM) and 1.3 <= ratio2 <= 2.3
        diag(7, 'PC2 urate raw beta_F; se_F; beta_M; se_M; z_diff; ratio', f'{y.b_rF}; {y.s_rF}; {y.b_rM}; {y.s_rM}; {zd2:.2f}; {ratio2:.3f}',
             'S6 needs |z_diff| >= 4 and |beta_F| > |beta_M|; predicted z_diff < -10, ratio 1.3-2.3',
             ('pass' if pc2_ok else 'FAIL S6') + ('; prediction met' if pred_met else '; prediction not fully met'))
        CTRL.append(dict(control='PC2', trait='urate', scale='raw', statistic='z_diff; beta_F/beta_M', value=f'{zd2:.3f}; {ratio2:.4f}',
                         prediction='z_diff < -10, |beta_F|>|beta_M|, ratio 1.3-2.3', outcome=('pass' if pc2_ok else 'FAIL S6') +
                         ('; prediction met' if pred_met else '; prediction not fully met')))
    else:
        diag(7, 'PC2 urate raw', 'raw urate unavailable', '', 'FAIL S6')
        CTRL.append(dict(control='PC2', trait='urate', scale='raw', statistic='', value='', prediction='', outcome='FAIL S6 (unavailable)'))
    raw_arm_failed_global = not pc2_ok
    if raw_arm_failed_global:
        for t in TRAITS:
            raw_ok_trait[t] = False

    # ---------------------------------------------------------- 8. NC1 (S10, S7)
    ss = np.random.SeedSequence(a.seed).spawn(5)
    scale_ok = {(t, sc): trait_ok[t] and (sc == 'irnt' or raw_ok_trait[t]) for t in TRAITS for sc in SCALES}
    gc = {(t, sc): 1.0 for t in TRAITS for sc in SCALES}
    nc1_rows = []
    for ti, t in enumerate(TRAITS):
        if not trait_ok[t]:
            continue
        rng = np.random.default_rng(ss[ti])
        d = data[t][data[t].is_pool].reset_index(drop=True)
        zc = (d.b_iF / d.s_iF ** 2 + d.b_iM / d.s_iM ** 2) / np.sqrt(1 / d.s_iF ** 2 + 1 / d.s_iM ** 2)
        pc = 2 * norm.sf(np.abs(zc))
        fin_ = np.isfinite(d.b_iF) & np.isfinite(d.s_iF) & np.isfinite(d.b_iM) & np.isfinite(d.s_iM) & (d.s_iF > 0) & (d.s_iM > 0)
        if raw_ok_trait[t]:
            fin_ &= np.isfinite(d.b_rF) & np.isfinite(d.s_rF) & np.isfinite(d.b_rM) & np.isfinite(d.s_rM) & (d.s_rF > 0) & (d.s_rM > 0)
        cand = d[fin_ & (pc > 0.5)].reset_index(drop=True)
        order = rng.permutation(len(cand))
        lead_by_chr = {}
        for r in leads[t].itertuples():
            lead_by_chr.setdefault(r.chr, []).append(r.pos)
        for c in lead_by_chr:
            lead_by_chr[c].sort()
        placed = {}; sel = []

        def near(lst, q):
            j = bisect.bisect_left(lst, q - 999999)
            return j < len(lst) and lst[j] <= q + 999999
        cchr = cand['chr'].values; cpos = cand['pos'].values
        for j in order:
            c, q = cchr[j], int(cpos[j])
            if near(lead_by_chr.get(c, []), q) or near(placed.get(c, []), q):
                continue
            bisect.insort(placed.setdefault(c, []), q)
            sel.append(j)
            if len(sel) == 1000:
                break
        nsel = len(sel)
        diag(8, f'{t} NC1 pool (every 20th eligible SNV) / combined p > 0.5 and finite / placed', f'{len(d)} / {len(cand)} / {nsel}',
             '1000 placed (S10: >= 500)', 'pass' if nsel >= 500 else 'FAIL S10 (trait untestable)')
        if nsel < 500:
            scale_ok[(t, 'irnt')] = scale_ok[(t, 'raw')] = False
            CTRL.append(dict(control='NC1', trait=t, scale='both', statistic='placed', value=nsel, prediction='>= 500', outcome='FAIL S10'))
            continue
        x = cand.iloc[sel]
        for sc in SCALES:
            if not scale_ok[(t, sc)]:
                continue
            tg = 'i' if sc == 'irnt' else 'r'
            z = (x['b_' + tg + 'F'] - x['b_' + tg + 'M']) / np.sqrt(x['s_' + tg + 'F'] ** 2 + x['s_' + tg + 'M'] ** 2)
            lam = float(np.median(z ** 2) / CHI2_MED)
            tail = float(np.mean(np.abs(z) > Z95))
            pred = (0.85 <= lam <= 1.15) and (0.035 <= tail <= 0.065)
            if lam > 1.30:
                act = 'FAIL S7 (trait x scale untestable)'; scale_ok[(t, sc)] = False
            elif lam > 1.15:
                act = 'S7 GC-corrected (se x sqrt(lambda))'; gc[(t, sc)] = lam
            else:
                act = 'pass (no correction)'
            diag(8, f'{t} NC1 {sc}: lambda; tail fraction |z|>1.96 (n {nsel})', f'{lam:.4f}; {tail:.4f}',
                 'lambda 0.85-1.15, tail 0.035-0.065', act + ('; prediction met' if pred else '; prediction not met'))
            CTRL.append(dict(control='NC1', trait=t, scale=sc, statistic='lambda; tail', value=f'{lam:.4f}; {tail:.4f}',
                             prediction='lambda 0.85-1.15; tail 0.035-0.065', outcome=act + ('; prediction met' if pred else '; prediction not met')))
            nc1_rows.append(dict(trait=t, scale=sc, n=nsel, lam=lam, tail=tail))
        x[['variant', 'b_iF', 's_iF', 'b_iM', 's_iM', 'b_rF', 's_rF', 'b_rM', 's_rM']].assign(trait=t).to_csv(
            OUT + f'/nc1_snps_{t}.tsv', sep='\t', index=False)

    for t in TRAITS:
        for sc in SCALES:
            diag(8, f'{t} {sc} arm status after stop rules', 'testable' if scale_ok[(t, sc)] else 'untestable',
                 '', f'GC lambda {gc[(t, sc)]:.4f}' if gc[(t, sc)] != 1.0 else '')
    raw_any = any(scale_ok[(t, 'raw')] for t in TRAITS)

    # GC application helper: se used = se x sqrt(lambda)
    def se_used(df, t_col, sc, sx):
        tg = 'i' if sc == 'irnt' else 'r'
        f_ = df[t_col].map(lambda t: math.sqrt(gc[(t, sc)])).astype(float)
        return df['s_' + tg + sx] * f_

    # ---------------------------------------------------------- 9. z_raw vs z_irnt at leads, within sex
    raw_flag = {t: False for t in TRAITS}
    for t in TRAITS:
        x = lds_t[lds_t.trait == t]
        if not (scale_ok[(t, 'raw')] and scale_ok[(t, 'irnt')]) or len(x) < 3:
            continue
        for sx in 'FM':
            r_ = float(np.corrcoef(x['b_r' + sx] / x['s_r' + sx], x['b_i' + sx] / x['s_i' + sx])[0, 1])
            flag = r_ < 0.8
            raw_flag[t] |= flag
            diag(9, f'{t} within-sex corr(z_raw, z_irnt) at leads, {sx} (n {len(x)})', round(r_, 4), '> 0.9 (flag raw arm if < 0.8)',
                 'FLAG raw arm' if flag else ('pass' if r_ > 0.9 else 'below prediction, no flag'))
    for sx in 'FM':
        x = lds_t[[scale_ok[(t, 'raw')] for t in lds_t.trait]]
        if len(x) > 2:
            r_ = float(np.corrcoef(x['b_r' + sx] / x['s_r' + sx], x['b_i' + sx] / x['s_i' + sx])[0, 1])
            diag(9, f'pooled within-sex corr(z_raw, z_irnt) at leads, {sx} (n {len(x)})', round(r_, 4), '> 0.9')

    # ---------------------------------------------------------- 10. bounded Fieller counts
    for nm, df in [('F1', tests), ('F2', lds_t)]:
        for sc in SCALES:
            tg = 'i' if sc == 'irnt' else 'r'
            ok_rows = df[[scale_ok[(t, sc)] for t in df.trait]]
            sM = se_used(ok_rows, 'trait', sc, 'M')
            nb = int((np.abs(ok_rows['b_' + tg + 'M'] / sM) > Z95).sum())
            diag(10, f'{nm} tests with a bounded Fieller interval, {sc}', f'{nb} of {len(ok_rows)}')
            sD = np.sqrt(se_used(ok_rows, 'trait', sc, 'F') ** 2 + sM ** 2)
            if len(ok_rows):
                zc_ = norm.isf(0.05 / len(ok_rows) / 2)
                diag(10, f'{nm} {sc}: median SE(beta_F - beta_M); BH single-hit z; MDE80 (SEs only)',
                     f'{float(np.median(sD)):.5f}; {zc_:.3f}; {(zc_ + norm.isf(0.2)) * float(np.median(sD)):.4f}',
                     'irnt F1: 0.00526; 3.748; 0.0242 (prespec G2)' if (nm == 'F1' and sc == 'irnt') else '')

    log('All diagnostics written. Controls: PC1 pass; PC2 ' + ('pass' if pc2_ok else 'FAIL') + '. Computing outcome families.')
    write_tables(OUT)

    # ========================================================== OUTCOME
    def family(df, label):
        df = df.copy()
        for sc in SCALES:
            tg = 'i' if sc == 'irnt' else 'r'
            ok_rows = np.array([scale_ok[(t, sc)] for t in df.trait])
            sF = se_used(df, 'trait', sc, 'F'); sM = se_used(df, 'trait', sc, 'M')
            st_ = contrast_stats(df['b_' + tg + 'F'], sF, df['b_' + tg + 'M'], sM)
            for k, v in st_.items():
                v = np.asarray(v)
                if v.dtype.kind == 'f':
                    v = np.where(ok_rows, v, np.nan)
                elif v.dtype.kind == 'b':
                    v = np.where(ok_rows, v, False)
                else:
                    v = np.where(ok_rows, v, '')
                df[f'{k}_{sc}'] = v
            df[f'se_F_used_{sc}'] = np.where(ok_rows, sF, np.nan)
            df[f'se_M_used_{sc}'] = np.where(ok_rows, sM, np.nan)
            df[f'q_{sc}'] = bh(df[f'p_{sc}'].values)
            df[f'q_eq_{sc}'] = bh(df[f'p_eq_{sc}'].values)
            df[f'gc_lambda_{sc}'] = [gc[(t, sc)] if scale_ok[(t, sc)] else np.nan for t in df.trait]
            log(f'{label} {sc}: BH family size {int(np.isfinite(df[f"p_{sc}"]).sum())}; TOST family size {int(np.isfinite(df[f"p_eq_{sc}"]).sum())}')
        both_ok = np.array([scale_ok[(t, 'raw')] and scale_ok[(t, 'irnt')] for t in df.trait])
        st_ = assign_states(df)
        df['state'] = np.where(both_ok, st_, 'untestable')
        df['raw_arm_flag_diag9'] = [raw_flag[t] for t in df.trait]
        return df

    F1 = family(tests, 'F1')
    F2 = family(lds_t, 'F2')

    # gene states
    gcols = ['trait', 'top_snp']
    keep_cols = ['b_iF', 's_iF', 'b_iM', 's_iM', 'b_rF', 's_rF', 'b_rM', 's_rM', 'state', 'raw_arm_flag_diag9']
    for sc in SCALES:
        keep_cols += [f'se_F_used_{sc}', f'se_M_used_{sc}', f'z_diff_{sc}', f'p_{sc}', f'q_{sc}', f'fieller_lo_{sc}', f'fieller_hi_{sc}',
                      f'fieller_bounded_{sc}', f'fieller_inside_{sc}', f'p_eq_{sc}', f'q_eq_{sc}', f'gc_lambda_{sc}', f'direction_{sc}']
    g = genes.merge(F1[gcols + keep_cols].rename(columns={'trait': 'best_enzyme_trait'}), on=['best_enzyme_trait', 'top_snp'], how='left')
    g['state'] = np.where(g.status == 'not_applicable', 'not_applicable',
                          np.where(g.status == 'testable', g.state, 'untestable'))
    if raw_arm_failed_global:
        g['state'] = np.where(g.state == 'not_applicable', g.state, 'untestable')
    g['untestable_reason'] = np.where(g.state == 'untestable', np.where(g.status == 'testable', 'stop rule removed trait or scale', g.status), '')
    for sc in SCALES:
        g[f'gc_corrected_{sc}'] = g[f'gc_lambda_{sc}'].fillna(1.0) > 1.0
    ren = {'b_iF': 'beta_F_irnt', 's_iF': 'se_F_irnt_reported', 'b_iM': 'beta_M_irnt', 's_iM': 'se_M_irnt_reported',
           'b_rF': 'beta_F_raw', 's_rF': 'se_F_raw_reported', 'b_rM': 'beta_M_raw', 's_rM': 'se_M_raw_reported'}
    g = g.rename(columns=ren).drop(columns=['status'])
    g.to_csv(OUT + '/catalog_sex_effect_state.tsv', sep='\t', index=False)
    F2.rename(columns=ren).to_csv(OUT + '/locus_sex_effect.tsv', sep='\t', index=False)
    F1.rename(columns=ren).to_csv(OUT + '/gene_arm_tests.tsv', sep='\t', index=False)

    # ---------------------------------------------------------- F3 chrX
    xrows = []
    for t in TRAITS:
        if not trait_ok[t]:
            continue
        d = data[t][data[t].is_x].copy()
        d = d[d.snv & ~d.low_iF & ~d.low_iM & (d.maf_iF >= 0.01) & (d.maf_iM >= 0.01) & np.isfinite(d.b_iF) & np.isfinite(d.b_iM)
              & (d.s_iF > 0) & (d.s_iM > 0)]
        d = d[[not any(lo <= p_ <= hi for lo, hi in PAR) for p_ in d.pos]]
        for model, k in [('XI', 1.0), ('ESC', 2.0)]:
            zc = (d.b_iF / d.s_iF ** 2 + k * d.b_iM / (k ** 2 * d.s_iM ** 2)) / np.sqrt(1 / d.s_iF ** 2 + 1 / (k ** 2 * d.s_iM ** 2))
            pcomb = 2 * norm.sf(np.abs(zc))
            cand = d.assign(zc=zc, pcomb=pcomb)[pcomb < 5e-8]
            # clumping by combined p; |z| is used for ordering because p underflows to 0 for |z| > 37.5 (same order)
            cand = cand.assign(az=-np.abs(cand.zc)).sort_values(['az', 'pos'])
            kept = []
            for r in cand.itertuples():
                if all(abs(r.pos - q) > 500000 for q in kept):
                    kept.append(r.pos)
            L = cand[cand.pos.isin(kept)].drop_duplicates('pos')
            g4 = len(L) >= 5
            diag('F3', f'{t} chrX lead set L_{model}: candidates p<5e-8 / clumped leads', f'{len(cand)} / {len(L)}', '>= 5 leads (G4)',
                 'pass' if (g4 and chrx_ok[t]) else ('FAIL G4 (untestable)' if not g4 else 'S9 untestable'))
            if not (g4 and chrx_ok[t]):
                continue
            for r in L.itertuples():
                rec = dict(trait=t, model=model, variant=r.variant, pos=r.pos, z_combined=r.zc, p_combined=r.pcomb,
                           maf_F=r.maf_iF, maf_M=r.maf_iM)
                for sc in SCALES:
                    tg = 'i' if sc == 'irnt' else 'r'
                    bF, bM = getattr(r, 'b_' + tg + 'F'), getattr(r, 'b_' + tg + 'M')
                    lam = gc[(t, sc)]
                    sF, sM = getattr(r, 's_' + tg + 'F') * math.sqrt(lam), getattr(r, 's_' + tg + 'M') * math.sqrt(lam)
                    if scale_ok[(t, sc)] and np.isfinite(bF) and np.isfinite(bM) and np.isfinite(sF) and np.isfinite(sM):
                        z = (bF - k * bM) / math.sqrt(sF ** 2 + k ** 2 * sM ** 2)
                        rec.update({f'beta_F_{sc}': bF, f'se_F_{sc}': sF, f'beta_M_{sc}': bM, f'se_M_{sc}': sM,
                                    f'z_{sc}': z, f'p_{sc}': 2 * norm.sf(abs(z)), f'gc_lambda_{sc}': lam,
                                    f'direction_{sc}': str(direction([bF], [bM], k)[0])})
                    else:
                        rec.update({f'z_{sc}': np.nan, f'p_{sc}': np.nan})
                xrows.append(rec)
    X = pd.DataFrame(xrows)
    if len(X):
        for model in ('XI', 'ESC'):
            for sc in SCALES:
                m = X.model == model
                X.loc[m, f'q_{sc}'] = bh(X.loc[m, f'p_{sc}'].values)
                log(f'F3 {model} {sc}: BH family size {int(np.isfinite(X.loc[m, f"p_{sc}"]).sum())}')

        def single(r):
            sr = np.isfinite(r.q_raw) and r.q_raw < 0.05
            si = np.isfinite(r.q_irnt) and r.q_irnt < 0.05
            if sr and si:
                return 'rejected' if r.direction_raw == r.direction_irnt else 'scale_dependent'
            if sr != si:
                return 'one_scale_only'
            return 'not_rejected'
        for c in ('direction_raw', 'direction_irnt'):
            if c not in X:
                X[c] = ''
        X['test_state'] = [single(r) for r in X.itertuples()]
        X['joint_state'] = ''
        for t in X.trait.unique():
            xi = X[(X.trait == t) & (X.model == 'XI')].set_index('variant')
            es = X[(X.trait == t) & (X.model == 'ESC')].set_index('variant')
            for v in set(xi.index) & set(es.index):
                xr_, er_ = xi.loc[v, 'test_state'] == 'rejected', es.loc[v, 'test_state'] == 'rejected'
                js = ('consistent_neither' if (xr_ and er_) else 'consistent_XI_only' if (er_ and not xr_)
                      else 'consistent_ESC_only' if (xr_ and not er_) else 'consistent_both')
                X.loc[(X.trait == t) & (X.variant == v), 'joint_state'] = js
    X.to_csv(OUT + '/chrX_sex_effect.tsv', sep='\t', index=False)

    # ---------------------------------------------------------- F4 set-level ratio
    f4 = []
    rng_lead = np.random.default_rng(ss[3]); rng_gene = np.random.default_rng(ss[4])
    for arm, df, rng in [('locus', F2, rng_lead), ('gene', F1, rng_gene)]:
        for t in TRAITS:
            x = df[df.trait == t].copy()
            if arm == 'gene':
                units = x.cluster.values
                uniq = np.unique(units)
                n_units = len(uniq)
                umap = {u: i for i, u in enumerate(uniq)}
                uidx = np.array([umap[u] for u in units]) if len(x) else np.array([], int)
            else:
                n_units = len(x)
                uidx = np.arange(n_units)
            if n_units == 0:
                continue
            idx = rng.integers(0, n_units, size=(B_BOOT, n_units))      # drawn ALT, AST, GGT in order; shared by both scales
            counts = np.zeros((B_BOOT, n_units))
            for bi in range(B_BOOT):
                counts[bi] = np.bincount(idx[bi], minlength=n_units)
            for sc in SCALES:
                if not scale_ok[(t, sc)]:
                    continue
                tg = 'i' if sc == 'irnt' else 'r'
                bF, bM = x['b_' + tg + 'F'].values, x['b_' + tg + 'M'].values
                sF, sM = x[f'se_F_used_{sc}'].values, x[f'se_M_used_{sc}'].values
                w = 1 / (sF ** 2 + sM ** 2)
                s = np.sign(bF / sF ** 2 + bM / sM ** 2)
                num_u = np.bincount(uidx, weights=w * s * bF, minlength=n_units)
                den_u = np.bincount(uidx, weights=w * s * bM, minlength=n_units)
                Rhat = num_u.sum() / den_u.sum()
                num_b = counts @ num_u; den_b = counts @ den_u
                with np.errstate(invalid='ignore', divide='ignore'):
                    Rb = num_b / den_b
                    lr = np.log(Rb)
                nbad = int((~np.isfinite(lr)).sum())
                lrf = lr[np.isfinite(lr)]
                ci95 = np.exp(np.percentile(lrf, [2.5, 97.5]))
                a6 = 0.05 / 6
                ci99 = np.exp(np.percentile(lrf, [100 * a6 / 2, 100 * (1 - a6 / 2)]))
                Bv = len(lrf)
                p = min(1.0, 2 * min((1 + (lrf <= 0).sum()) / (Bv + 1), (1 + (lrf >= 0).sum()) / (Bv + 1)))
                sdr = sd[t]['ratio_raw'] if sc == 'raw' else 1.0
                f4.append(dict(arm=arm, trait=t, scale=sc, n_snps=len(x), n_units=n_units, R=Rhat, ci95_lo=ci95[0], ci95_hi=ci95[1],
                               ci9917_lo=ci99[0], ci9917_hi=ci99[1], p_boot=p, n_boot_nonfinite=nbad,
                               sd_ratio_raw_PC3=sd[t]['ratio_raw'], gc_corrected=gc[(t, sc)] > 1.0,
                               bound_090_111=(ci99[0] >= F4_LO and ci99[1] <= F4_HI)))
    F4 = pd.DataFrame(f4)
    if len(F4):
        for sc in SCALES:
            m = F4.scale == sc
            F4.loc[m, 'p_holm'] = holm(F4.loc[m, 'p_boot'].values)
    F4.to_csv(OUT + '/set_level_ratio.tsv', sep='\t', index=False)

    # ---------------------------------------------------------- descriptive extras (nominal, never called significant)
    extras = []
    g_t = g[g.state != 'not_applicable']
    extras.append('Gene states (all 426): ' + str(g.state.value_counts().to_dict()))
    for cls in ['genetically_anchored', 'concordant']:
        extras.append(f'States, {cls}: ' + str(g[g.catalog_class == cls].state.value_counts().to_dict()))
    g8 = g[(g.pp_h4 >= 0.8) & (g.untestable_reason == '') & (g.state != 'not_applicable')]
    extras.append(f'PP.H4 >= 0.8 subset (tested genes): {len(g8)} genes (prespec 221), {g8.drop_duplicates(["best_enzyme_trait", "top_snp"]).shape[0]} '
                  f'unique trait x SNP (prespec 207); states ' + str(g8.state.value_counts().to_dict()))
    nom = F1[(F1.p_raw < 0.05) & (F1.p_irnt < 0.05) & (F1.direction_raw == F1.direction_irnt)]
    extras.append(f'F1 tests with nominal p < 0.05 on both scales and the same direction (descriptive): {len(nom)} of {len(F1)}')
    ok_ = np.isfinite(F1.z_diff_raw) & np.isfinite(F1.z_diff_irnt)
    if ok_.sum() > 2:
        extras.append(f'Pearson r of z_diff raw vs irnt over F1 tests (descriptive): {np.corrcoef(F1.z_diff_raw[ok_], F1.z_diff_irnt[ok_])[0, 1]:.4f}')
    extras.append('F2 lead states: ' + str(F2.state.value_counts().to_dict()))
    for t in TRAITS:
        x = F2[F2.trait == t]
        extras.append(f'F2 {t} states: ' + str(x.state.value_counts().to_dict()) + '; direction among heterogeneous: '
                      + str(x[x.state == 'heterogeneous'].direction_irnt.value_counts().to_dict()))
    for e_ in extras:
        log(e_)

    write_tables(OUT)
    write_result(OUT, result_name, SYN, g, F1, F2, X, F4, scale_ok, raw_flag, relabel_s8, n_tab, gc, nc1_rows, pc2_ok, extras, sd, chrx_ok)


def write_tables(OUT):
    pd.DataFrame(DIAG).to_csv(OUT + '/diagnostics.tsv', sep='\t', index=False)
    pd.DataFrame(CTRL).to_csv(OUT + '/controls.tsv', sep='\t', index=False)


def fmt_int(n):
    return f'{int(n):,}'


def write_result(OUT, name, SYN, g, F1, F2, X, F4, scale_ok, raw_flag, relabel_s8, n_tab, gc, nc1_rows, pc2_ok, extras, sd, chrx_ok):
    tested = g[~g.state.isin(['not_applicable', 'untestable'])]
    n_tested = len(tested)
    cnt = g.state.value_counts().to_dict()
    N_h = cnt.get('heterogeneous', 0); N_b = cnt.get('homogeneous_within_bound', 0); N_u = cnt.get('underpowered', 0)
    N_one = cnt.get('one_scale_only', 0); N_sd = cnt.get('scale_dependent', 0)
    n_f1 = {sc: int(np.isfinite(F1[f'p_{sc}']).sum()) for sc in SCALES}
    raw_any = any(scale_ok[(t, 'raw')] for t in TRAITS)
    # MDE at run: median SE(diff) irnt x (BH single-hit z + z_0.8)
    sdiff = np.sqrt(F1.se_F_used_irnt ** 2 + F1.se_M_used_irnt ** 2)
    m = int(np.isfinite(sdiff).sum())
    mde = (norm.isf(0.05 / m / 2) + norm.isf(0.2)) * float(np.nanmedian(sdiff)) if m else float('nan')
    nA = n_tab.get('ALT', {})
    if not raw_any or not pc2_ok:
        verdict = 'untestable (raw arm failed: S2/S6): no Catalog column; irnt results descriptive only'
    elif N_h > 0:
        verdict = f'pass_with_heterogeneous_genes (N_h = {N_h})'
    else:
        verdict = 'null (no heterogeneous gene)'
    L = []
    L.append('# H72 result: female-minus-male per-allele effects at Catalog enzyme signals (UK Biobank, Neale round 2)\n')
    if SYN:
        L.append('**SYNTHETIC OUTCOME DEBUG RUN. Female and male enzyme betas were simulated under equal effects. Not a result.**\n')
    L.append(f'Run directory: `{OUT}`. Prespec: `PRESPEC_H72.md` (sha256 {PRESPEC_SHA}). Seed 20260923.\n')
    L.append(f'## Verdict (computed by the prespec rule)\n\n{verdict}\n')
    L.append('## Primary family F1 (Catalog column)\n')
    L.append(f'- Unit: UK Biobank participants in two disjoint sex-specific GWAS (ALT n women {fmt_int(nA.get("irnt_F", 0))}, '
             f'men {fmt_int(nA.get("irnt_M", 0))}). Tests are SNPs, not replicates.')
    L.append(f'- BH family size per scale: raw {n_f1["raw"]}, irnt {n_f1["irnt"]} trait x SNP tests; TOST family the same.')
    L.append(f'- Gene states (426 anchored or concordant genes): heterogeneous {N_h}; scale_dependent {N_sd}; one_scale_only {N_one}; '
             f'homogeneous_within_bound {N_b}; underpowered {N_u}; untestable {cnt.get("untestable", 0)}; not_applicable {cnt.get("not_applicable", 0)}.')
    L.append(f'- Tests with q < 0.05: raw {int((F1.q_raw < 0.05).sum())}, irnt {int((F1.q_irnt < 0.05).sum())}; '
             f'q_eq < 0.05: raw {int((F1.q_eq_raw < 0.05).sum())}, irnt {int((F1.q_eq_irnt < 0.05).sum())}.')
    L.append(f'- Median per-SNP detectable difference (irnt, 80% power, single BH hit): {mde:.4f} SD per allele.')
    het = g[g.state.isin(['heterogeneous', 'scale_dependent', 'one_scale_only', 'homogeneous_within_bound'])]
    if len(het):
        L.append('\n| gene | class | trait | top SNP | state | dir raw | dir irnt | z raw | q raw | z irnt | q irnt | ratio CI irnt |')
        L.append('|---|---|---|---|---|---|---|---|---|---|---|---|')
        for r in het.sort_values(['state', 'q_irnt']).itertuples():
            L.append(f'| {r.symbol} | {r.catalog_class} | {r.best_enzyme_trait} | {r.top_snp} | {r.state} | {r.direction_raw} | {r.direction_irnt} | '
                     f'{r.z_diff_raw:.2f} | {r.q_raw:.3g} | {r.z_diff_irnt:.2f} | {r.q_irnt:.3g} | {r.fieller_lo_irnt:.2f} to {r.fieller_hi_irnt:.2f} |')
    L.append('\n## Proposed reading text (returned as data; nothing installed)\n')
    if verdict.startswith('untestable'):
        L.append('No column. Limitations: sex heterogeneity of the genetic anchors was not tested.')
    elif N_h > 0:
        L.append(f'"At the shared liver-enzyme signal of {N_h} of {n_tested} genetically anchored or concordant genes, the per-allele '
                 f'effect differed between women and men on both the raw and the standardized scale (BH q < 0.05 over {n_f1["irnt"]} SNP '
                 f'tests per scale; UK Biobank, {fmt_int(nA.get("irnt_F", 0))} women and {fmt_int(nA.get("irnt_M", 0))} men for ALT). '
                 f'{N_b} genes were bounded within a female-to-male ratio of 0.67 to 1.5, and {N_u} could not be bounded."'
                 + (f' {N_one} genes differed on one scale only and {N_sd} were scale-dependent.' if (N_one or N_sd) else ''))
    else:
        L.append(f'"No difference between women and men was detected at the enzyme signal of any of the {n_tested} tested genes on both '
                 f'scales (BH q >= 0.05 on at least one scale; median per-SNP detectable difference {mde:.3f} SD per allele on the '
                 f'standardized scale at 80% power for a single BH hit). {N_b} of {n_tested} were bounded within a female-to-male ratio '
                 f'of 0.67 to 1.5 (BH over TOST p on both scales); the remaining {N_u} could not be bounded."'
                 + (f' ({N_one} genes differed on one scale only (BH q < 0.05) and were not detected on the other; {N_sd} were scale-dependent.)'
                    if (N_one or N_sd) else ''))
    if len(F4) and not verdict.startswith('untestable'):
        for arm, lab in (('locus', 'loci'), ('gene', 'gene-arm SNPs')):
            b_ = F4[(F4.arm == arm) & F4.bound_090_111]
            if len(b_):
                L.append(f'"Averaged over {lab}, the female-to-male ratio was ' + '; '.join(
                    f'{r.R:.3f} (99.17% CI {r.ci9917_lo:.3f} to {r.ci9917_hi:.3f}) for {r.trait} on the {r.scale} scale' for r in b_.itertuples())
                    + ' (in UK Biobank sex-stratified GWAS, as previously analysed by Flynn et al. 2020 and Bernabeu et al. 2021)."')
    old_t = 'Most genetically supported genes derive from liver-enzyme rather than direct MASLD traits.'
    if verdict.startswith('untestable'):
        add = 'Sex heterogeneity of the genetic anchors was not tested.'
    else:
        add = (f'Their enzyme signals were compared between women and men in UK Biobank on raw and standardized scales: '
               f'{N_h} of {n_tested} differed on both scales, {N_b} were bounded within a 1.5-fold sex ratio and {N_u} could not be bounded'
               + (f'; {N_one} differed on one scale only' if N_one else '') + (f' and {N_sd} were scale-dependent' if N_sd else '')
               + ' (Supplementary Table).')
    L.append('\nProposed Limitations edit (returned as data, not installed):\n')
    L.append('- file: `docs/manuscript/draft/fig5_discussion.md`, line 232')
    L.append(f'- old_text: "{old_t}"')
    L.append(f'- new_text: "{old_t} {add}"')
    L.append('\n## Secondary families\n')
    L.append('F2 (974 distance-clumped leads; re-analysis of public UK Biobank data, cf. Flynn 2020 and Bernabeu 2021): ' +
             str(F2.state.value_counts().to_dict()) +
             f'; q < 0.05 raw {int((F2.q_raw < 0.05).sum())}, irnt {int((F2.q_irnt < 0.05).sum())}.')
    if len(F4):
        L.append('\nF4 set-level ratio R = sum(w s beta_F) / sum(w s beta_M); bootstrap 10,000 (locus arm: leads; gene arm: 500 kb clusters); '
                 'Holm over 6 per scale; F4-eq bound uses the 99.17% interval inside [0.90, 1.11].\n')
        L.append('| arm | trait | scale | n SNPs | units | R | 95% CI | 99.17% CI | p boot | p Holm | bound | SD_F/SD_M raw (PC3) |')
        L.append('|---|---|---|---|---|---|---|---|---|---|---|---|')
        for r in F4.itertuples():
            L.append(f'| {r.arm} | {r.trait} | {r.scale} | {r.n_snps} | {r.n_units} | {r.R:.4f} | {r.ci95_lo:.4f}-{r.ci95_hi:.4f} | '
                     f'{r.ci9917_lo:.4f}-{r.ci9917_hi:.4f} | {r.p_boot:.3g} | {r.p_holm:.3g} | {r.bound_090_111} | {r.sd_ratio_raw_PC3:.4f} |')
        L.append('')
        for arm in ('locus', 'gene'):
            for t in TRAITS:
                x = F4[(F4.arm == arm) & (F4.trait == t)].set_index('scale')
                if not {'raw', 'irnt'} <= set(x.index):
                    continue
                ri, rr = x.loc['irnt'], x.loc['raw']
                near1_i = ri.ci9917_lo >= F4_LO and ri.ci9917_hi <= F4_HI
                near1_r = rr.ci9917_lo >= F4_LO and rr.ci9917_hi <= F4_HI
                sdr = rr.sd_ratio_raw_PC3
                near_sd_r = np.isfinite(sdr) and rr.ci9917_lo / sdr >= F4_LO and rr.ci9917_hi / sdr <= F4_HI
                if near1_i and near_sd_r:
                    s_ = 'on average, enzyme-locus effects were proportional to each sex\'s trait SD'
                elif near1_r and not near1_i:
                    s_ = 'on average, per-allele effects in U/L were equal'
                else:
                    s_ = 'no summary sentence (R values and intervals reported only)'
                L.append(f'- F4 {arm} arm, {t}: {s_}.')
    if len(X):
        L.append('\nF3 chrX (exploratory): leads per lead set and state: ' +
                 str(X.groupby(['trait', 'model']).size().to_dict()) + '; test states ' +
                 str(X.groupby(['model', 'test_state']).size().to_dict()) + '; joint states ' +
                 str(X[X.joint_state != ''].drop_duplicates(['trait', 'variant']).groupby(['trait', 'joint_state']).size().to_dict())
                 + ' (joint states counted per variant present in both lead sets)')
    else:
        L.append('\nF3 chrX: no testable lead set (G4 or S9).')
    L.append('\n## Controls and stop rules\n')
    for c in CTRL:
        L.append(f'- {c["control"]} {c["trait"]} {c["scale"]}: {c["statistic"]} = {c["value"]}; prediction {c["prediction"]}; {c["outcome"]}')
    L.append('- Arm status: ' + ', '.join(f'{t} {sc} {"testable" if scale_ok[(t, sc)] else "untestable"}'
                                          + (f' (GC lambda {gc[(t, sc)]:.3f})' if gc[(t, sc)] > 1 else '') for t in TRAITS for sc in SCALES))
    L.append('- Raw-arm flag from diagnostic 9 (z_raw vs z_irnt correlation < 0.8): ' + str(raw_flag))
    L.append('- S8 relabel (irnt implied SD outside 0.90-1.10): ' + str(relabel_s8))
    L.append('- chrX coding (S9) testable: ' + str(chrx_ok))
    L.append('\n## Descriptive (nominal, never called significant)\n')
    for e_ in extras:
        L.append(f'- {e_}')
    L.append('\n## Outputs\n\n`catalog_sex_effect_state.tsv`, `gene_arm_tests.tsv`, `locus_sex_effect.tsv`, `chrX_sex_effect.tsv`, '
             '`set_level_ratio.tsv`, `diagnostics.tsv`, `controls.tsv`, `nc1_snps_*.tsv`, `run.log`, figure panels in `figures/`.')
    with open(OUT + '/' + name, 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    log('VERDICT', verdict)


if __name__ == '__main__':
    main()
