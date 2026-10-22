#!/usr/bin/env python3
"""H74: at tier-1/2 liver-trait loci with a SuSiE-COLOC liver eQTL gene (COLOC r2),
how often is the protein-coding gene with the TSS nearest the anchor not a colocalized gene?

Implements PRESPEC_H74.md in the parent directory (sha256 checked below). Every definition,
gate, stop rule, control, diagnostic, test, family and seed comes from that file; section
names in comments refer to it.

Stages:
  --stage controls : input hashes (S1), recount (G1/S2), anchors (S8), S7, D1, D9,
                     PC-A (S4) and PC-B (S3) on real data, then exit. No E1, E2, E3 outcome.
  --stage all      : everything, in the prespec order, and RESULT_H74.md.
  --debug-null     : after the controls, every anchor is replaced by a uniform random
                     position in its window and gene states are shuffled within locus,
                     so all outcome code runs on null data. Used only to debug outcome code.
"""
import argparse, hashlib, json, math, os, re, sys, time
from functools import reduce
import numpy as np
import pandas as pd
from scipy import stats

R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
H74 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRESPEC_COPY = os.path.join(H74, 'PRESPEC_H74.md')
PRESPEC_ORIG = R + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H74.md'
PRESPEC_SHA = 'bc36945db073965e1b73cffbe45ff81ce2f04d7894b7a52e786139dab5d49a64'
R4 = 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/'
COLOC = 'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv'
GTF = 'data/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf'
TIER = 'GWAS/finemapping/config/gwas_trait_tier.tsv'
PAIRS = 'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz'
GATE_LS = R4 + 'gate/NG13/locus_structure.tsv'
# Inputs with sha256 (prespec table). Bytes for the two group-count files are from the
# sentence under the table (515 and 113 bytes).
INPUTS = [
    ('coloc_r2_all_rows', COLOC, 335138620, 'acdfa306f1ee45383758bc2522602fc2828c4060080b7ad2d96134fc4a399887'),
    ('gtf_hg19', GTF, 1171012259, '4b4388ed697f2b02d0450f6da001a43009ffe3eec46cf417139c6b459cf0fcde'),
    ('gwas_tier_table', TIER, 2296, '865f108de0782a708469b7a4e42943075de79b5ae062dc452f4050034bc15fc7'),
    ('gate_count_py', R4 + 'gate/NG13/count.py', 3910, '76fd30828359ca7491131a423d4482ec80ca1dc08c08215550f26c85bf835b41'),
    ('gate_sbatch', R4 + 'gate/NG13/count.sbatch', 710, '006eaf30227e136083238c5b26fb885fa2201c55fa97d45fe624a2dd2077106b'),
    ('gate_log', R4 + 'gate/NG13/count_22251943.log', 1419, '51c23935250719c795829cbb5acbca0ce82da4a6351ede3d830349d5e0165909'),
    ('gate_locus_table', GATE_LS, 8741, '04c4c8711c15be43dd1eef787aae0c54ea48aad13fbe5d5087daf17f24282ddb'),
    ('group_count_cmd', R4 + 'write/H74/count_groups.sh', 515, '403377144441883239f8556e31735bf5a8ceac777ecf98faea0f9712fa3076eb'),
    ('group_count_out', R4 + 'write/H74/count_groups.out', 113, 'ee8b4f867132cae8354e4a225a01db50a7e9f6b4df1131edd52f93e7746d5b3b'),
    ('signal_pairs', PAIRS, 298033, '3cc5c05959c6fe55c0776a5fd2ac58e8796ab5b95a2801b7e48b3667571ee1c3'),
]
MB = 1_000_000
N_LOCI_EXPECTED = 259
N_INF_EXPECTED = 248
B_PERM = 10_000
NC1_DRAWS = 1_000
UNTESTABLE = 'susie_untestable_insufficient_shared_posterior'
PREC = ['colocalized', 'bystander', 'tested_other', 'untested']
NAMED = ['SORT1', 'PNPLA3', 'TM6SF2', 'GCKR', 'HSD17B13', 'MBOAT7', 'APOE']
PCA_SNP = ('1', 109817590)
Z95 = stats.norm.ppf(0.975)


class Tee:
    def __init__(self, path):
        self.f = open(path, 'w')

    def write(self, s):
        sys.__stdout__.write(s); self.f.write(s)

    def flush(self):
        sys.__stdout__.flush(); self.f.flush()


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 24), b''):
            h.update(blk)
    return h.hexdigest()


def wilson(k, n):
    if n == 0:
        return (float('nan'), float('nan'))
    p = k / n; z = Z95; den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (c - h, c + h)


def clopper(k, n, a=0.05):
    lo = 0.0 if k == 0 else stats.beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - a / 2, k + 1, n - k)
    return (float(lo), float(hi))


def q3(x):
    x = np.asarray(x, float)
    return dict(median=float(np.median(x)), q25=float(np.quantile(x, 0.25)),
                q75=float(np.quantile(x, 0.75)), max=float(np.max(x)), n=int(len(x)))


def parse_snp(s):
    m = re.fullmatch(r'(\d+|X|Y):(\d+)', str(s))
    return (m.group(1), int(m.group(2))) if m else None


# ----------------------------------------------------------------------------- state
class Run:
    def __init__(self, out, seed, stage, debug):
        self.out, self.seed, self.stage, self.debug = out, seed, stage, debug
        self.summary = dict(analysis='H74', prespec_sha256=PRESPEC_SHA, seed=seed, stage=stage,
                            debug_null=debug, started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
        self.debug_rng = np.random.default_rng(np.random.SeedSequence([seed, 777])) if debug else None
        self.debug_on = False  # switched on after the controls in --debug-null mode

    def stop(self, rule, msg):
        print(f'\n*** STOP {rule}: {msg}', flush=True)
        self.summary['status'] = 'stopped_by_rule'
        self.summary['stop_rule'] = rule
        self.summary['stop_message'] = msg
        verdict = (f'Untestable ({rule}): {msg}. Per the prespec Reading, no number from this run '
                   f'enters the paper or the Gene Catalog.')
        self.summary['verdict'] = verdict
        self.write_summary()
        with open(os.path.join(self.out, 'RESULT_H74.md'), 'w') as f:
            f.write(('# DEBUG-NULL RUN, NOT A RESULT\n\n' if self.debug else '') +
                    f'# RESULT_H74\n\nStatus: stopped_by_rule ({rule}).\n\nVerdict: {verdict}\n\n'
                    f'Prespec sha256: {PRESPEC_SHA}. Seed {self.seed}. See summary.json and h74.log.\n')
        sys.exit(0)

    def write_summary(self):
        def conv(o):
            if o is pd.NA: return None
            if isinstance(o, (np.integer,)): return int(o)
            if isinstance(o, (np.floating,)): return float(o)
            if isinstance(o, (np.bool_,)): return bool(o)
            if isinstance(o, np.ndarray): return o.tolist()
            raise TypeError(type(o))
        with open(os.path.join(self.out, 'summary.json'), 'w') as f:
            json.dump(self.summary, f, indent=1, default=conv)


# ----------------------------------------------------------------------------- data
def load_gtf():
    """gene records; gene_id with version removed; _PAR_Y skipped; keyed by (gene_id, chromosome)."""
    rec = []; n_par = 0; n_gene = 0
    with open(R + GTF) as fh:
        for l in fh:
            if l[0] == '#' or '\tgene\t' not in l:
                continue
            f = l.split('\t')
            if f[2] != 'gene':
                continue
            n_gene += 1
            full = re.search('gene_id "([^"]+)"', f[8]).group(1)
            if full.endswith('_PAR_Y'):
                n_par += 1
                continue
            gid = re.search('gene_id "([^".]+)', f[8]).group(1)
            gt = re.search('gene_type "([^"]+)"', f[8]).group(1)
            gn = re.search('gene_name "([^"]+)"', f[8]).group(1)
            st, en = int(f[3]), int(f[4])
            tss = st if f[6] == '+' else en
            rec.append((gid, f[0].replace('chr', ''), tss, st, en, f[6], gt, gn))
    G = pd.DataFrame(rec, columns=['gid', 'chr', 'tss', 'start', 'end', 'strand', 'gtype', 'gname'])
    return G, n_gene, n_par


class Annot:
    def __init__(self, G):
        dup = G.duplicated(['gid', 'chr'], keep=False)
        if dup.any():
            raise SystemExit(f'(gene_id, chr) not unique in GTF: {G[dup].head()}')
        self.G = G
        self.key = {(g, c): i for i, (g, c) in enumerate(zip(G.gid, G.chr))}
        self.by_chr = {}
        for kind, sub in (('pc', G[G.gtype == 'protein_coding']), ('all', G)):
            for c, g in sub.groupby('chr'):
                self.by_chr[(kind, c)] = dict(tss=g.tss.values.astype(np.int64), start=g.start.values.astype(np.int64),
                                              end=g.end.values.astype(np.int64), gid=g.gid.values, name=g.gname.values)

    def rec(self, gid, c):
        i = self.key.get((gid, c))
        return None if i is None else self.G.iloc[i]

    def nearest(self, c, x, kind='pc', body=False):
        a = self.by_chr.get((kind, c))
        if a is None:
            return None
        d = np.abs(a['tss'] - x)
        if body:
            d = np.where((a['start'] <= x) & (x <= a['end']), 0, d)
        m = int(d.min()); idx = np.flatnonzero(d == m)
        return dict(dist=m, gids=list(a['gid'][idx]), names=list(a['name'][idx]))

    def n_within(self, c, x, w, kind='pc'):
        a = self.by_chr.get((kind, c))
        return 0 if a is None else int((np.abs(a['tss'] - x) <= w).sum())


def cluster(pos):
    """1-Mb single-linkage clustering of distinct positive top_snp positions; gate order (chr as string)."""
    locs = []
    for c, g in pos.groupby('snp_chr'):
        ps = sorted(g.snp_pos.unique()); cur = [ps[0]]
        for p in ps[1:]:
            if p - cur[-1] <= MB:
                cur.append(p)
            else:
                locs.append((c, cur[0], cur[-1])); cur = [p]
        locs.append((c, cur[0], cur[-1]))
    return locs


def anchor_row(p):
    return p.sort_values(['PP.H4.susie', 'top_snp_PP', 'gwas_name', 'ensembl'],
                         ascending=[False, False, True, True]).iloc[0]


def make_locus(run, A, PAIR, lid, c, a, b, srows, pos, lo=None, hi=None):
    lo = a - MB if lo is None else lo
    hi = b + MB if hi is None else hi
    w = srows[(srows.snp_chr == c) & (srows.snp_pos >= lo) & (srows.snp_pos <= hi)]
    g = w.groupby('ensembl').agg(h4=('PP.H4.susie', 'max'), h3=('PP.H3.susie', 'max'), gene=('gene', 'first'))
    recs = [A.rec(e, c) for e in g.index]
    g['tss'] = [np.nan if r is None else float(r.tss) for r in recs]
    g['gtype'] = [None if r is None else r.gtype for r in recs]
    g['gname'] = [gn if r is None else r.gname for r, gn in zip(recs, g.gene)]
    p = pos[(pos.snp_chr == c) & (pos.snp_pos >= a) & (pos.snp_pos <= b)]
    ar = anchor_row(p)
    pr = PAIR.get((ar.ensembl, ar.gwas_name))
    hits = []
    for e, gw in zip(p.ensembl, p.gwas_name):
        q = PAIR.get((e, gw))
        hits.append(None if q is None else q['hit1'])
    L = dict(id=lid, chr=c, a=int(a), b=int(b), lo=int(lo), hi=int(hi), genes=g, n_pos_rows=len(p),
             n_pos_snps=int(p.snp_pos.nunique()), anchor_row=ar, tier=None,
             hit1_str=None if pr is None else pr['hit1'], pair_PP_H4=None if pr is None else pr['PP.H4'],
             top_str=ar.top_snp, top_pos=int(ar.snp_pos),
             pos_hit1_strs=sorted(set(h for h in hits if h is not None)), n_pos_rows_no_pair=sum(h is None for h in hits))
    ph = parse_snp(L['hit1_str']) if L['hit1_str'] is not None else None
    L['hit1_parsed'] = ph
    L['anchor'] = ph[1] if (ph is not None and ph[0] == c) else None
    L['anchor_top'] = L['top_pos']
    L['hit1_all'] = [parse_snp(h)[1] for h in L['pos_hit1_strs'] if parse_snp(h) is not None and parse_snp(h)[0] == c]
    if run.debug_on:
        debugify(run, L)
    return L


def debugify(run, L):
    """--debug-null only: shuffle gene states within locus and replace every anchor by a random window position."""
    rng = run.debug_rng
    g = L['genes']
    perm = rng.permutation(len(g))
    g[['h4', 'h3']] = g[['h4', 'h3']].values[perm]
    L['anchor'] = int(rng.integers(L['lo'], L['hi'] + 1))
    L['anchor_top'] = int(rng.integers(L['lo'], L['hi'] + 1))
    L['hit1_all'] = [int(rng.integers(L['lo'], L['hi'] + 1)) for _ in L['hit1_all']]
    L['debug'] = True


def gene_class(L, gid, thr=0.5):
    g = L['genes']
    if gid in g.index:
        h4, h3 = g.at[gid, 'h4'], g.at[gid, 'h3']
        if h4 > thr:
            return 'colocalized'
        if h3 > 0.8:
            return 'bystander'
        return 'tested_other'
    return 'untested'


def nearest_class(L, near, thr=0.5):
    cls = [gene_class(L, gid, thr) for gid in near['gids']]
    return min(cls, key=PREC.index), cls


def untested_subclass(D, DIDX, tier34, L, gid):
    idx = DIDX.get(gid)
    if idx is None:
        return 'iv_no_row'
    r = D.iloc[idx]
    r = r[(r.snp_chr == L['chr']) & (r.snp_pos >= L['lo']) & (r.snp_pos <= L['hi'])]
    if len(r) == 0:
        return 'iv_no_row'
    if (r.gwas_name.isin(tier34) & r['PP.H4.susie'].notna()).any():
        return 'i_tier34_susie'
    if (r.method == UNTESTABLE).any():
        return 'ii_susie_untestable'
    meth = set(r.method)
    if not meth <= {'abf_only', 'abf_fallback'}:
        raise SystemExit(f'untested gene {gid} has unexpected methods {meth} at locus {L["id"]}')
    return 'iii_abf_only'


def flag_untestable12(D, DIDX, t12, L, gids):
    out = False
    for gid in gids:
        idx = DIDX.get(gid)
        if idx is None:
            continue
        r = D.iloc[idx]
        r = r[(r.snp_chr == L['chr']) & (r.snp_pos >= L['lo']) & (r.snp_pos <= L['hi']) &
              r.gwas_name.isin(t12) & (r.method == UNTESTABLE)]
        out = out or len(r) > 0
    return out


# ----------------------------------------------------------------------------- E3 machinery
def e3_setup(L, anchor, thr=0.5):
    g = L['genes']; m = g[g.tss.notna()]
    inC = (m.h4 > thr).values
    n, c = len(m), int(inC.sum())
    if c == 0 or c == n:
        return None
    tss = m.tss.values.astype(np.int64)
    d = np.abs(tss - anchor)
    near = d == d.min()
    return dict(id=L['id'], n=n, c=c, inC=inC, ranks=stats.rankdata(d, method='average'), near=near,
                t=int(near.sum()), tss=tss, dist=d, gids=list(m.index), lo=L['lo'], hi=L['hi'])


def exact_TA(setups, inC_override=None):
    """Exact within-locus relabelling null of T_A: convolution of per-locus score distributions
    (hypergeometric K/t for tied nearest sets), on a grid with denominator lcm(t_i)."""
    Lc = reduce(math.lcm, [S['t'] for S in setups], 1)
    dist = np.array([1.0]); obs = 0
    for j, S in enumerate(setups):
        inC = S['inC'] if inC_override is None else inC_override[j]
        c = int(inC.sum()); step = Lc // S['t']
        pmf = np.zeros(Lc + 1)
        for k in range(0, min(S['t'], c) + 1):
            pmf[k * step] += stats.hypergeom.pmf(k, S['n'], c, S['t'])
        dist = np.convolve(dist, pmf)
        obs += int(inC[S['near']].sum()) * step
    units = np.arange(len(dist)) / Lc
    mean = float((units * dist).sum()); sd = float(math.sqrt(((units - mean) ** 2 * dist).sum()))
    p = float(dist[obs:][::-1].sum())
    return dict(lcm=Lc, obs_units=obs, T_A=obs / Lc, p=p, mean=mean, sd=sd, dist=dist, total_mass=float(dist.sum()))


def perm_stats(setups, B, rng, Lc, inC_override=None):
    """Row 0 = observed labels; rows 1..B = within-locus random relabellings with |C_i| fixed.
    Same arithmetic for every row, so observed and null values are comparable without tolerance."""
    TA = np.zeros(B + 1, dtype=np.int64); TB = np.zeros(B + 1)
    for j, S in enumerate(setups):
        inC = S['inC'] if inC_override is None else inC_override[j]
        n, c = S['n'], int(inC.sum())
        lab = np.zeros((B + 1, n), dtype=bool); lab[0] = inC
        if B > 0:
            idx = rng.random((B, n)).argsort(axis=1)[:, :c]
            sub = lab[1:]
            np.put_along_axis(sub, idx, True, axis=1)
        TA += lab[:, S['near']].sum(axis=1) * (Lc // S['t'])
        TB += (lab * (S['ranks'] - 0.5)).sum(axis=1) / (c * n)
    TB /= len(setups)
    return TA, TB


def nc1(setups, rng, ndraw):
    tot = np.zeros(ndraw)
    for S in setups:
        xs = rng.integers(S['lo'], S['hi'] + 1, size=ndraw)
        d = np.abs(S['tss'][None, :] - xs[:, None])
        mask = d == d.min(axis=1, keepdims=True)
        tot += (mask & S['inC'][None, :]).sum(axis=1) / mask.sum(axis=1)
    return tot / len(setups)


def holm(ps):
    order = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (m - r) * ps[i])); adj[i] = run
    return adj


# ----------------------------------------------------------------------------- E1/E2 helpers
def e1_e2(loci, A, key='anchor', thr=0.5, kind='pc', body=False):
    rows = []
    for L in loci:
        x = L[key]
        if x is None:
            continue
        nr = A.nearest(L['chr'], x, kind=kind, body=body)
        cl, _ = nearest_class(L, nr, thr)
        rows.append(dict(id=L['id'], cls=cl, contra=cl != 'colocalized'))
    n = len(rows); k = sum(r['contra'] for r in rows); kb = sum(r['cls'] == 'bystander' for r in rows)
    ku = sum(r['cls'] == 'untested' for r in rows)
    return dict(n=n, E1_k=k, E1=k / n if n else float('nan'), E1_wilson=wilson(k, n), E2_k=kb,
                E2=kb / n if n else float('nan'), untested_k=ku)


def ta_exact_for(loci, key='anchor', thr=0.5):
    setups = [S for S in (e3_setup(L, L[key], thr) for L in loci if L[key] is not None) if S is not None]
    if not setups:
        return dict(n_inf=0)
    ex = exact_TA(setups)
    TA, TB = perm_stats(setups, 0, None, ex['lcm'])
    return dict(n_inf=len(setups), T_A=ex['T_A'], T_A_null_mean=ex['mean'], T_A_exact_p=ex['p'], T_B=float(TB[0]))


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--stage', choices=['controls', 'all'], default='all')
    ap.add_argument('--debug-null', action='store_true')
    ap.add_argument('--B', type=int, default=B_PERM, help='honoured only with --debug-null')
    ap.add_argument('--nc1-draws', type=int, default=NC1_DRAWS, help='honoured only with --debug-null')
    args = ap.parse_args()
    if not os.path.isdir(args.out):
        raise SystemExit(f'--out {args.out} does not exist')
    if args.seed != 20260923:
        raise SystemExit('prespec fixes seed 20260923')
    B = args.B if args.debug_null else B_PERM
    NDRAW = args.nc1_draws if args.debug_null else NC1_DRAWS
    sys.stdout = Tee(os.path.join(args.out, 'h74.log'))
    run = Run(args.out, args.seed, args.stage, args.debug_null)
    S = run.summary
    S['B'] = B; S['nc1_draws'] = NDRAW
    print(f'H74 start {S["started_utc"]} out={args.out} stage={args.stage} debug_null={args.debug_null} B={B} NC1={NDRAW}')
    print('script sha256', sha256(os.path.abspath(__file__)))

    # ---- S1: prespec and input hashes
    t0 = time.time(); bad = []; rows = []
    for nm, path in (('prespec_copy', PRESPEC_COPY), ('prespec_plans', PRESPEC_ORIG)):
        h = sha256(path); ok = h == PRESPEC_SHA; rows.append((nm, path, os.path.getsize(path), h, PRESPEC_SHA, ok))
        if not ok: bad.append(nm)
    for nm, p, nb, h0 in INPUTS:
        h = sha256(R + p); sz = os.path.getsize(R + p); ok = (h == h0) and (sz == nb)
        rows.append((nm, p, sz, h, h0, ok))
        if not ok: bad.append(nm)
    pd.DataFrame(rows, columns=['role', 'path', 'bytes', 'sha256_observed', 'sha256_expected', 'match']).to_csv(
        os.path.join(args.out, 'input_sha256_check.tsv'), sep='\t', index=False)
    print(f'S1 input hash check ({time.time() - t0:.0f}s):', 'all 12 match' if not bad else f'MISMATCH {bad}')
    S['S1'] = dict(passed=not bad, mismatched=bad)
    if bad:
        run.stop('S1', f'sha256 or byte count differs for {bad}')

    # ---- load
    tier = pd.read_csv(R + TIER, sep='\t')
    t12 = set(tier.study_name[tier.tier.isin([1, 2])]); tier34 = set(tier.study_name[tier.tier.isin([3, 4])])
    tmap = tier.set_index('study_name').tier.to_dict()
    plat = set(tier.study_name[tier.study_name.str.startswith('MVP_Platelet')])
    cols = ['gwas_name', 'gene', 'ensembl', 'chr', 'PP.H3.susie', 'PP.H4.susie', 'PP.H4.abf', 'method', 'top_snp', 'top_snp_PP']
    D = pd.read_csv(R + COLOC, usecols=cols)
    sp_ = D.top_snp.str.split(':', expand=True)
    D['snp_chr'] = sp_[0]; D['snp_pos'] = sp_[1].astype(np.int64)
    DIDX = D.groupby('ensembl').indices
    s = D[D.gwas_name.isin(t12) & D['PP.H4.susie'].notna()].copy()
    print(f'rows {len(D)} studies {D.gwas_name.nunique()} methods {D.method.value_counts().to_dict()}')
    print(f'tier12 names {len(t12)}; tier12 susie rows {len(s)} studies {s.gwas_name.nunique()} genes {s.ensembl.nunique()}')
    print('tier12 susie rows with PP.H3.susie NA:', int(s['PP.H3.susie'].isna().sum()))

    t0 = time.time()
    G, n_gene, n_par = load_gtf()
    A = Annot(G)
    print(f'GTF gene records {n_gene}, _PAR_Y skipped {n_par}, kept {len(G)}, protein_coding {int((G.gtype == "protein_coding").sum())} ({time.time() - t0:.0f}s)')

    sp = pd.read_csv(R + PAIRS, sep='\t')
    sp = sp[sp.eligible.astype(str).str.upper() == 'TRUE']
    sp = sp.sort_values(['PP.H4', 'idx1', 'idx2'], ascending=[False, True, True], kind='mergesort')
    best = sp.drop_duplicates(['ensembl', 'gwas'], keep='first')
    PAIR = {(e, g): dict(hit1=h1, hit2=h2, **{'PP.H4': p4}) for e, g, h1, h2, p4 in
            zip(best.ensembl, best.gwas, best.hit1, best.hit2, best['PP.H4'])}
    print(f'signal pairs eligible {len(sp)}, best per (ensembl, gwas) {len(best)}')

    # ---- G1: recount
    pos = s[s['PP.H4.susie'] > 0.5]
    locs = cluster(pos)
    loci = [make_locus(run, A, PAIR, i, c, a, b, s, pos) for i, (c, a, b) in enumerate(locs)]
    for L in loci:
        L['tier'] = int(tmap[L['anchor_row'].gwas_name])
    n_loci = len(loci)
    g_map_rows = s.apply(lambda r: (r.ensembl, r.snp_chr) in A.key, axis=1).mean()
    inf_ids = []
    for L in loci:
        g = L['genes']; m = g[g.tss.notna()]
        if ((m.h4 > 0.5).sum() >= 1) and ((m.h4 <= 0.5).sum() >= 1):
            inf_ids.append(L['id'])
    n_inf = len(inf_ids)
    gate = pd.read_csv(R + GATE_LS, sep='\t')
    mine = pd.DataFrame([dict(chr=L['chr'], start=L['a'], end=L['b'], n_tested=len(L['genes']),
                              n_tested_pc=int((L['genes'].gtype == 'protein_coding').sum()),
                              n_h4=int((L['genes'].h4 > 0.5).sum()),
                              n_h3=int(((L['genes'].h3 > 0.8) & (L['genes'].h4 <= 0.5)).sum())) for L in loci])
    same = len(mine) == len(gate) and (mine[['chr', 'start', 'end', 'n_tested', 'n_h4', 'n_h3']].astype(str).values ==
                                       gate[['chr', 'start', 'end', 'n_tested', 'n_h4', 'n_h3']].astype(str).values).all()
    same_pc = len(mine) == len(gate) and (mine.n_tested_pc.values == gate.n_tested_pc.values).all()
    unmapped_genes = sorted({e for L in loci for e in L['genes'].index[L['genes'].tss.isna()]})
    D1 = dict(n_loci=n_loci, n_informative=n_inf, per_locus_identical_to_gate=bool(same),
              n_tested_pc_identical_to_gate=bool(same_pc),
              n_tested=q3(mine.n_tested), n_tested_mean=float(mine.n_tested.mean()), gate_n_tested_mean=float(gate.n_tested.mean()),
              n_h4_eq1=int((mine.n_h4 == 1).sum()), n_h4_ge2=int((mine.n_h4 >= 2).sum()),
              all_tested_coloc=int((mine.n_tested == mine.n_h4).sum()), n_h3_eq0=int((mine.n_h3 == 0).sum()),
              multi_snp_clusters=int(sum(L['n_pos_snps'] > 1 for L in loci)),
              tss_mapping_rate_rows=float(g_map_rows), tested_genes_without_tss_on_locus_chr=unmapped_genes,
              par_y_records_skipped=n_par, gtf_gene_records=n_gene,
              n_positive_rows=int(len(pos)))
    S['D1'] = D1
    print('\nD1 recount:', json.dumps(D1))
    S['G1'] = dict(n_loci=n_loci, n_informative=n_inf, passed=(n_loci == N_LOCI_EXPECTED and n_inf == N_INF_EXPECTED))
    if not S['G1']['passed']:
        run.stop('S2', f'G1 failed: {n_loci} loci (need 259), {n_inf} informative (need 248)')
    print('G1 passed: 259 loci, 248 informative')
    S['G2'] = dict(E1_ci_halfwidth_at_0p5=round(Z95 * math.sqrt(0.25 / 259), 4), T_A_MDE_gate=0.063)

    # ---- S8: anchors
    bad8 = [L['id'] for L in loci if L['hit1_str'] is None or L['hit1_parsed'] is None or L['hit1_parsed'][0] != L['chr']]
    S['S8'] = dict(passed=not bad8, bad_loci=bad8)
    print(f'S8 anchors: {n_loci - len(bad8)} of {n_loci} loci have a parseable hit1 on the locus chromosome')
    if bad8:
        run.stop('S8', f'anchor hit1 missing, unparseable or off-chromosome at loci {bad8}')
    eqpp = sum(abs(L['pair_PP_H4'] - L['anchor_row']['PP.H4.susie']) < 1e-9 for L in loci)
    print(f'anchor rows whose best eligible pair PP.H4 equals PP.H4.susie: {eqpp} of {n_loci}')
    tie_loci = 0
    for L in loci:
        p = pos[(pos.snp_chr == L['chr']) & (pos.snp_pos >= L['a']) & (pos.snp_pos <= L['b'])]
        top = p[p['PP.H4.susie'] == p['PP.H4.susie'].max()]
        tie_loci += top.top_snp.nunique() > 1
    print(f'loci with max PP.H4.susie tied across different top SNPs (tie-break used): {tie_loci}')

    # ---- D9 anchor provenance (positions only, no gene)
    dd = np.array([abs(L['anchor'] - L['top_pos']) for L in loci])
    D9 = dict(hit1_equals_top_snp=int(sum(L['hit1_str'] == L['top_str'] for L in loci)),
              dist_hit1_to_top_snp=q3(dd),
              hit1_outside_ab=int(sum(not (L['a'] <= L['anchor'] <= L['b']) for L in loci)),
              hit1_outside_window=int(sum(not (L['lo'] <= L['anchor'] <= L['hi']) for L in loci)),
              anchor_rows_PP_H4_abf_le_0p5=int(sum(L['anchor_row']['PP.H4.abf'] <= 0.5 for L in loci)),
              anchor_pair_PPH4_equals_susie=int(eqpp), loci_max_PPH4_tied_across_top_snps=int(tie_loci),
              window_lo_below_1=int(sum(L['lo'] < 1 for L in loci)))
    S['D9'] = D9
    print('D9 anchor provenance:', json.dumps(D9))

    # ---- S7
    no_pc = [L['id'] for L in loci if A.n_within(L['chr'], L['anchor'], 2 * MB) == 0]
    S['S7'] = dict(n_loci_no_pc_within_2Mb=len(no_pc), loci=no_pc, passed=len(no_pc) <= 5)
    print(f'S7: loci with no protein-coding TSS within 2 Mb of the anchor: {len(no_pc)}')
    if len(no_pc) > 5:
        run.stop('S7', f'{len(no_pc)} loci have no protein-coding gene within 2 Mb of the anchor')

    # ---- PC-A (1p13)
    pcA_check = A.nearest(*PCA_SNP)
    print(f'PC-A code check: nearest protein-coding TSS to 1:109817590 = {pcA_check}')
    sort1 = G[(G.gname == 'SORT1') & (G.gtype == 'protein_coding')].iloc[0]
    cand = [L for L in loci if L['chr'] == '1' and L['lo'] <= PCA_SNP[1] <= L['hi'] and sort1.gid in L['genes'].index]
    if len(cand) != 1:
        run.stop('S4', f'1p13 locus not identified uniquely ({len(cand)} candidates)')
    L1 = cand[0]
    nr = A.nearest('1', L1['anchor']); cl1, cls1 = nearest_class(L1, nr)
    C1 = sorted(L1['genes'].gname[L1['genes'].h4 > 0.5])
    pcA = dict(locus=L1['id'], window=[L1['lo'], L1['hi']], anchor=L1['hit1_str'], anchor_top_snp=L1['top_str'],
               nearest_pc=nr['names'], nearest_dist=nr['dist'], nearest_class=cl1, C=C1,
               anchor_is_rs12740374_hg19=L1['hit1_str'] == '1:109817590')
    if pcA['anchor_is_rs12740374_hg19']:
        if nr['names'] != ['PSRC1']:
            S['PC_A'] = pcA
            run.stop('S4', f'1p13 anchor is 1:109817590 but nearest protein-coding gene is {nr["names"]}, not PSRC1')
        pcA['passed'] = cl1 == 'colocalized'
    else:
        pcA['passed'] = cl1 == 'colocalized'
    pcA['mybphl_case'] = 'MYBPHL' in nr['names']
    S['PC_A'] = pcA
    print('PC-A:', json.dumps(pcA, default=str))

    # ---- E3 setups on real anchors (structure only; labels not read for PC-B (i)/(ii))
    setups = [e3_setup(L, L['anchor']) for L in loci if L['id'] in inf_ids]
    assert all(x is not None for x in setups) and len(setups) == N_INF_EXPECTED
    Lc = reduce(math.lcm, [x['t'] for x in setups], 1)
    tied = [(x['id'], x['t']) for x in setups if x['t'] > 1]
    print(f'E3 setups {len(setups)}; loci with tied nearest tested genes: {tied}; lcm {Lc}')
    ss = np.random.SeedSequence(args.seed).spawn(3)
    rng0, rng1, rng2 = (np.random.default_rng(c) for c in ss)

    # ---- PC-B synthetic labels
    ov1 = [x['near'].copy() for x in setups]
    ex1 = exact_TA(setups, ov1); TA1, TB1 = perm_stats(setups, 0, None, ex1['lcm'], ov1)
    ov2 = []
    for x in setups:
        far = np.flatnonzero(x['dist'] == x['dist'].max())
        j = sorted(far, key=lambda k: x['gids'][k])[0]
        v = np.zeros(x['n'], bool); v[j] = True; ov2.append(v)
    ex2 = exact_TA(setups, ov2)
    ov3 = []
    for x in setups:
        v = np.zeros(x['n'], bool); v[rng2.random(x['n']).argsort()[:x['c']]] = True; ov3.append(v)
    ex3 = exact_TA(setups, ov3)
    exr = exact_TA(setups)  # null mean and sd for the real |C_i| (labels of the observed C are not used for mean/sd)
    pcB = dict(i_T_A=ex1['T_A'], i_T_B=float(TB1[0]), i_exact_p=ex1['p'],
               i_pass=(ex1['obs_units'] == len(setups) * ex1['lcm']) and TB1[0] < 0.5 and ex1['p'] < 1e-10,
               ii_T_A=ex2['T_A'], ii_pass=ex2['obs_units'] == 0,
               iii_T_A=ex3['T_A'], iii_null_mean=exr['mean'], iii_null_sd=exr['sd'],
               iii_z=(ex3['T_A'] - exr['mean']) / exr['sd'],
               iii_pass=abs(ex3['T_A'] - exr['mean']) <= 3 * exr['sd'],
               exact_dist_total_mass=exr['total_mass'])
    pcB['passed'] = bool(pcB['i_pass'] and pcB['ii_pass'] and pcB['iii_pass'])
    S['PC_B'] = pcB
    print('PC-B:', json.dumps(pcB, default=float))
    if not pcB['passed']:
        run.stop('S3', 'PC-B synthetic-label check failed: ' + json.dumps(pcB, default=float))
    S['G3'] = dict(PC_A_passed=pcA['passed'], PC_B_passed=pcB['passed'])
    print(f'G3: PC-A {"passed" if pcA["passed"] else "FAILED (not a stop rule unless S4; no Catalog flag)"}; PC-B passed')

    if args.stage == 'controls':
        S['status'] = 'controls_only'
        run.write_summary()
        print('\nstage=controls: stopping before any E1, E2 or E3 outcome.')
        return

    # =================================================================== outcome section
    if args.debug_null:
        run.debug_on = True
        for L in loci:
            debugify(run, L)
        print('\n*** DEBUG-NULL: anchors displaced and gene states shuffled within locus; every number below is null data ***')
        setups = [e3_setup(L, L['anchor']) for L in loci]
        setups = [x for x in setups if x is not None]
        Lc = reduce(math.lcm, [x['t'] for x in setups], 1)
    inf_set = {x['id'] for x in setups}

    # ---- per-locus nearest protein-coding gene (primary anchor)
    LT = []
    for L in loci:
        nr = A.nearest(L['chr'], L['anchor']); cl, cls = nearest_class(L, nr)
        g = L['genes']; Cg = g[g.h4 > 0.5]
        sub = None
        if cl == 'untested':
            subs = [untested_subclass(D, DIDX, tier34, L, gid) for gid, c_ in zip(nr['gids'], cls) if c_ == 'untested']
            sub = sorted(subs)[0]
        flag = flag_untestable12(D, DIDX, t12, L, nr['gids'])
        Cm = Cg[Cg.tss.notna()]
        dC = int(np.abs(Cm.tss.values - L['anchor']).min()) if len(Cm) else None
        C_pc = int((Cg.gtype == 'protein_coding').sum())
        n500 = A.n_within(L['chr'], L['anchor'], 500_000)
        nr0 = A.nearest(L['chr'], L['anchor_top']); cl0, _ = nearest_class(L, nr0)
        L.update(near=nr, near_class=cl, near_classes=cls, untested_sub=sub, untestable_flag=flag, dist_C=dC,
                 C_pc=C_pc, n_pc500=n500, near0=nr0, near0_class=cl0)
        ar = L['anchor_row']
        LT.append(dict(locus=L['id'], chr=L['chr'], start=L['a'], end=L['b'], window_lo=L['lo'], window_hi=L['hi'],
                       anchor_tier=L['tier'], anchor_gwas=ar.gwas_name, anchor_gene=ar.gene, anchor_ensembl=ar.ensembl,
                       anchor_PP_H4_susie=ar['PP.H4.susie'], anchor_PP_H4_abf=ar['PP.H4.abf'],
                       anchor_hit1=L['hit1_str'] if not L.get('debug') else f'debug:{L["anchor"]}',
                       anchor_top_snp=L['top_str'] if not L.get('debug') else f'debug:{L["anchor_top"]}',
                       nearest_pc=';'.join(nr['names']), nearest_pc_ensembl=';'.join(nr['gids']), nearest_pc_dist=nr['dist'],
                       nearest_pc_tie=len(nr['gids']) > 1, nearest_class=cl, nearest_classes=';'.join(cls),
                       untested_subclass=sub, nearest_has_tier12_untestable_row=flag, contradiction=cl != 'colocalized',
                       C=';'.join(sorted(Cg.gname)), n_C=len(Cg), n_C_protein_coding=C_pc, n_tested=len(g),
                       n_tested_mapped=int(g.tss.notna().sum()), dist_nearest_C_tss=dC, n_pc_tss_within_500kb=n500,
                       informative_E3=L['id'] in inf_set,
                       A0_nearest_pc=';'.join(nr0['names']), A0_nearest_dist=nr0['dist'], A0_nearest_class=cl0,
                       A0_contradiction=cl0 != 'colocalized'))
    LT = pd.DataFrame(LT)

    # ---- diagnostics D2-D8 (printed before E1-E3)
    multi = [L for L in loci if len(L['hit1_all']) > 1]
    differ = 0; a1 = []
    for L in loci:
        hs = L['hit1_all'] if L['hit1_all'] else [L['anchor']]
        sets = []; contra = []
        for x in hs:
            nr = A.nearest(L['chr'], x); cl, _ = nearest_class(L, nr)
            sets.append(tuple(sorted(nr['gids']))); contra.append(cl != 'colocalized')
        if len(hs) > 1 and len(set(sets)) > 1:
            differ += 1
        a1.append(float(np.mean(contra)))
    S['D2'] = dict(loci_gt1_distinct_hit1=len(multi), of_which_nearest_pc_differs=differ,
                   loci_with_positive_rows_lacking_pair=int(sum(L['n_pos_rows_no_pair'] > 0 for L in loci)))
    S['D3'] = dict(dist_to_nearest_pc_tss=q3(LT.nearest_pc_dist), dist_to_nearest_C_tss=q3(LT.dist_nearest_C_tss.dropna()))
    S['D4'] = dict(loci_C_without_protein_coding=int((LT.n_C_protein_coding == 0).sum()))
    vc = LT.untested_subclass.value_counts().to_dict()
    S['D5'] = dict(untested_subclasses={k: int(vc.get(k, 0)) for k in ['i_tier34_susie', 'ii_susie_untestable', 'iii_abf_only', 'iv_no_row']},
                   nearest_flagged_tier12_untestable=int(LT.nearest_has_tier12_untestable_row.sum()))
    d6 = {}
    for t_, g_ in LT.groupby('anchor_tier'):
        k = int(g_.contradiction.sum()); n = len(g_)
        d6[f'tier{t_}'] = dict(n=n, E1_k=k, E1=k / n, wilson=wilson(k, n))
    S['D6'] = d6
    q1, q2 = np.quantile(LT.n_pc_tss_within_500kb, [1 / 3, 2 / 3])
    bot = LT[LT.n_pc_tss_within_500kb <= q1]; top = LT[LT.n_pc_tss_within_500kb >= q2]
    S['D8'] = dict(n_pc500=q3(LT.n_pc_tss_within_500kb), cut_low=float(q1), cut_high=float(q2),
                   bottom=dict(n=len(bot), E1_k=int(bot.contradiction.sum()), E1=float(bot.contradiction.mean())),
                   top=dict(n=len(top), E1_k=int(top.contradiction.sum()), E1=float(top.contradiction.mean())))
    # D7 named loci
    NL = []
    for nm in NAMED:
        gr = G[(G.gname == nm) & (G.gtype == 'protein_coding') & G.chr.str.fullmatch(r'\d+')]
        if len(gr) == 0:
            NL.append(dict(named_gene=nm, status='not_in_gtf')); continue
        gr = gr.iloc[0]
        hit = [L for L in loci if L['chr'] == gr.chr and L['lo'] <= gr.tss <= L['hi']]
        if not hit:
            NL.append(dict(named_gene=nm, status='absent', gene_chr=gr.chr, gene_tss=int(gr.tss))); continue
        for L in hit:
            NL.append(dict(named_gene=nm, status='present', gene_chr=gr.chr, gene_tss=int(gr.tss), locus=L['id'],
                           anchor_hit1=LT.loc[L['id'], 'anchor_hit1'], nearest_pc=';'.join(L['near']['names']),
                           nearest_dist=L['near']['dist'], nearest_class=L['near_class'],
                           named_gene_class=gene_class(L, gr.gid), C=LT.loc[L['id'], 'C'],
                           prediction=('agreement (PC-A)' if nm == 'SORT1' else 'none')))
    NL = pd.DataFrame(NL)
    if 'locus' in NL:
        NL['locus'] = NL['locus'].astype('Int64')
    S['D7'] = NL.to_dict(orient='records')
    for k in ['D2', 'D3', 'D4', 'D5', 'D6', 'D8']:
        print(f'{k}:', json.dumps(S[k], default=float))
    print('D7 named loci:\n' + NL.to_string(index=False))

    # ---- E1, E2 and the split
    n = len(LT); k1 = int(LT.contradiction.sum()); kb = int((LT.nearest_class == 'bystander').sum())
    ku = int((LT.nearest_class == 'untested').sum()); ko = int((LT.nearest_class == 'tested_other').sum())
    kc = int((LT.nearest_class == 'colocalized').sum())
    E = dict(n_loci=n, E1_k=k1, E1=k1 / n, E1_wilson=wilson(k1, n), E1_clopper_pearson=clopper(k1, n),
             E2_k=kb, E2=kb / n, E2_wilson=wilson(kb, n), E2_of_contradicted=(kb / k1 if k1 else float('nan')),
             E2_of_contradicted_wilson=wilson(kb, k1),
             split=dict(colocalized=kc, bystander=kb, tested_other=ko, untested=ku,
                        untested_subclasses=S['D5']['untested_subclasses']),
             nearest_ties=LT.loc[LT.nearest_pc_tie, ['locus', 'nearest_pc', 'nearest_classes']].to_dict(orient='records'))
    S['E1_E2'] = E
    print('\nE1/E2:', json.dumps(E, default=float))

    # ---- F1: T_A exact + permutation, T_B permutation
    ex = exact_TA(setups)
    TA, TB = perm_stats(setups, B, rng0, ex['lcm'])
    assert TA[0] == ex['obs_units']
    pA_perm = (1 + int((TA[1:] >= TA[0]).sum())) / (B + 1)
    pB = (1 + int((TB[1:] <= TB[0]).sum())) / (B + 1)
    tol = 3 * math.sqrt(ex['p'] * (1 - ex['p']) / B) + 1 / (B + 1)
    s5_ok = abs(ex['p'] - pA_perm) <= tol
    adj = holm(np.array([ex['p'], pB]))
    F1 = dict(n_informative=len(setups), T_A=ex['T_A'], T_A_rate=ex['T_A'] / len(setups), T_A_null_mean=ex['mean'],
              T_A_null_rate=ex['mean'] / len(setups), T_A_null_sd=ex['sd'], T_A_exact_p=ex['p'], T_A_perm_p=pA_perm,
              S5_tolerance=tol, S5_passed=bool(s5_ok), T_B=float(TB[0]), T_B_null_mean=float(TB[1:].mean()),
              T_B_perm_p=pB, holm_adj=dict(T_A=float(adj[0]), T_B=float(adj[1])),
              rejects=bool((adj < 0.05).any()), tie_lcm=ex['lcm'], B=B)
    S['F1'] = F1
    np.savez_compressed(os.path.join(args.out, 'null_T_A_T_B.npz'), T_A_units=TA, T_B=TB, lcm=ex['lcm'],
                        exact_T_A_dist=ex['dist'])
    print('F1:', json.dumps(F1, default=float))
    if not s5_ok:
        run.write_summary()
        run.stop('S5', f'exact T_A p {ex["p"]:.4g} and permutation p {pA_perm:.4g} differ by more than {tol:.4g}')

    # ---- NC-1 displaced anchor
    rates = nc1(setups, rng1, NDRAW)
    obs_rate = ex['T_A'] / len(setups)
    p_nc1 = (1 + int((rates >= obs_rate - 1e-12).sum())) / (NDRAW + 1)
    NC1 = dict(draws=NDRAW, mean_rate=float(rates.mean()), sd_rate=float(rates.std(ddof=1)),
               q025=float(np.quantile(rates, 0.025)), q975=float(np.quantile(rates, 0.975)),
               observed_rate=obs_rate, nominal_p=p_nc1, predicted_range=[0.15, 0.30],
               mean_in_predicted_range=bool(0.15 <= rates.mean() <= 0.30),
               S6_mean_outside_0p10_0p35=bool(not (0.10 <= rates.mean() <= 0.35)), S6_p_gt_0p05=bool(p_nc1 > 0.05))
    NC1['S6_triggered'] = NC1['S6_mean_outside_0p10_0p35'] or NC1['S6_p_gt_0p05']
    S['NC1'] = NC1
    with np.load(os.path.join(args.out, 'null_T_A_T_B.npz')) as z:
        keep = dict(z)
    np.savez_compressed(os.path.join(args.out, 'null_T_A_T_B.npz'), **keep, nc1_rates=rates)
    print('NC-1:', json.dumps(NC1, default=float))

    # ---- sensitivities A0-A7 and platelet (exact nulls only; no random draws)
    SENS = []

    def add(name, d, note=''):
        r = dict(analysis=name, n_loci=d.get('n'), E1_k=d.get('E1_k'), E1=d.get('E1'),
                 E1_wilson_lo=(d.get('E1_wilson') or (None, None))[0], E1_wilson_hi=(d.get('E1_wilson') or (None, None))[1],
                 E2_k=d.get('E2_k'), E2=d.get('E2'), untested_k=d.get('untested_k'), n_informative=d.get('n_inf'),
                 T_A=d.get('T_A'), T_A_null_mean=d.get('T_A_null_mean'), T_A_exact_p_nominal=d.get('T_A_exact_p'),
                 T_B=d.get('T_B'), note=note)
        SENS.append(r)

    prim = e1_e2(loci, A); prim.update(ta_exact_for(loci)); prim['T_B'] = F1['T_B']
    add('primary', prim, 'hit1 anchor; F1 p values are in summary.json')
    a0 = e1_e2(loci, A, key='anchor_top'); a0.update(ta_exact_for(loci, key='anchor_top'))
    add('A0_top_snp_anchor', a0, 'T_B printed without p (no permutation)')
    add('A1_every_positive_hit1', dict(n=len(a1), E1=float(np.mean(a1))), 'E1 = mean over loci of per-locus contradicted fraction')
    a2 = e1_e2(loci, A, kind='all'); add('A2_nearest_any_gene_type', a2)
    a3 = e1_e2(loci, A, thr=0.8); a3.update(ta_exact_for(loci, thr=0.8)); add('A3_coloc_PPH4_gt_0p8', a3)
    a4 = e1_e2(loci, A, body=True); add('A4_gene_body_distance0', a4)
    d4ids = set(LT.locus[LT.n_C_protein_coding == 0])
    a5 = e1_e2([L for L in loci if L['id'] not in d4ids], A); add('A5_excluding_D4_loci', a5)
    t1 = [L for L in loci if L['tier'] == 1]
    a6 = e1_e2(t1, A); a6.update(ta_exact_for(t1)); add('A6_tier1_only', a6)
    # A7 merge overlapping windows (single linkage of adjacent loci)
    groups = []
    for c in sorted({L['chr'] for L in loci}, key=lambda x: x):
        ls = sorted([L for L in loci if L['chr'] == c], key=lambda L: L['a'])
        cur = [ls[0]]
        for L in ls[1:]:
            if L['a'] - MB <= max(x['b'] for x in cur) + MB:
                cur.append(L)
            else:
                groups.append(cur); cur = [L]
        groups.append(cur)
    merged = []
    for j, grp in enumerate(groups):
        if len(grp) == 1:
            merged.append(grp[0]); continue
        a_, b_ = min(x['a'] for x in grp), max(x['b'] for x in grp)
        M = make_locus(run, A, PAIR, f'm{j}', grp[0]['chr'], a_, b_, s, pos)
        merged.append(M)
    print(f'A7: {len(loci)} loci -> {len(merged)} merged loci')
    a7 = e1_e2(merged, A); a7.update(ta_exact_for(merged))
    add('A7_merged_overlapping_windows', a7, f'{len(loci)} loci -> {len(merged)} merged loci; '
        f'{sum(len(g_) > 1 for g_ in groups)} merged groups')
    # platelet description
    spl = D[D.gwas_name.isin(plat) & D['PP.H4.susie'].notna()].copy()
    posp = spl[spl['PP.H4.susie'] > 0.5]
    lp = [make_locus(run, A, PAIR, f'p{i}', c, a, b, spl, posp) for i, (c, a, b) in enumerate(cluster(posp))]
    nmiss = sum(L['anchor'] is None for L in lp)
    pl = e1_e2(lp, A); add('platelet_description_hit1', pl, f'{len(lp)} platelet loci; {nmiss} without an eligible-pair hit1 (excluded here)')
    pl0 = e1_e2(lp, A, key='anchor_top'); add('platelet_description_top_snp', pl0, f'{len(lp)} platelet loci, top_snp anchor')
    SENS = pd.DataFrame(SENS)
    S['sensitivities'] = SENS.to_dict(orient='records')
    print('\nSensitivities:\n' + SENS.drop(columns=['note']).to_string(index=False))

    # ---- verdict by the Reading rule
    lo, hi = E['E1_wilson']; hw = (hi - lo) / 2
    X, Y, Z = E['E1_k'], E['E2_k'], E['split']['untested']
    lines = []
    lines.append(f'At {X} of 259 tier-1/2 liver-trait loci with a colocalized liver eQTL gene (SuSiE-COLOC PP.H4 > 0.5), '
                 f'the protein-coding gene with the TSS nearest the GWAS lead variant of the colocalizing signal was not a '
                 f'colocalized gene (E1 = {E["E1"]:.3f}, 95% CI {lo:.3f} to {hi:.3f}). At {Y} loci the nearest gene had been '
                 f'tested and showed a distinct liver eQTL signal (PP.H3 > 0.8); at {Z} loci it had no multi-signal '
                 f'colocalization test. A colocalized nearest gene counted as agreement even when it is not the gene with '
                 f'functional evidence: at 1p13 the nearest protein-coding gene to the anchor is '
                 f'{"/".join(pcA["nearest_pc"])} ({pcA["nearest_class"]})'
                 + (', not SORT1.' if 'SORT1' not in pcA['nearest_pc'] else '.'))
    calib = None
    if F1['rejects'] and not NC1['S6_triggered']:
        calib = (f'Colocalized genes lay nearer the lead variant than tested genes placed at random in the same window '
                 f'(observed rate {F1["T_A_rate"]:.3f} against null {F1["T_A_null_rate"]:.3f}; Holm-adjusted p '
                 f'T_A {F1["holm_adj"]["T_A"]:.3g}, T_B {F1["holm_adj"]["T_B"]:.3g}; calibration, not a finding).')
        lines.append(calib)
    elif F1['rejects'] and NC1['S6_triggered']:
        lines.append(f'F1 rejects, but S6 is triggered (NC-1 mean rate {NC1["mean_rate"]:.3f}, NC-1 nominal p '
                     f'{NC1["nominal_p"]:.3g}), so the "nearer than chance" sentence is not written.')
    else:
        lines.append(f'F1 did not reject (Holm-adjusted p T_A {F1["holm_adj"]["T_A"]:.3g}, T_B {F1["holm_adj"]["T_B"]:.3g}; '
                     f'gate MDE 6.3 points). D3, D8 and D9 must be checked for a coordinate or build error before any text '
                     f'is written. No sentence claims proximity carries no information.')
    if E['E1'] < 0.15:
        lines.append('At these loci the nearest gene usually is a colocalized gene, so nearest-gene lists are mostly '
                     'consistent with liver eQTL colocalization here, limited to loci with a colocalized gene.')
    a0diff = abs(a0['E1'] - E['E1']) > hw
    if a0diff:
        lines.append(f'Anchor dependence: with the GWAS-side lead of the colocalizing signal pair (hit1) E1 = {E["E1"]:.3f}; '
                     f'with the coloc.abf top SNP (top_snp) E1 = {a0["E1"]:.3f}.')
    catalog = 'no Catalog flag may be written from this run (PC-A failed)' if not pcA['passed'] else \
        'Catalog nearest-gene flag only after user adoption'
    S['reading'] = dict(sentences=lines, A0_differs_more_than_halfwidth=bool(a0diff), E1_ci_halfwidth=hw,
                        E1_below_0p15=bool(E['E1'] < 0.15), catalog=catalog)
    verdict = (f'Completed. E1 = {X}/259 = {E["E1"]:.3f} (Wilson 95% CI {lo:.3f} to {hi:.3f}); E2 = {Y}/259 = {E["E2"]:.3f} '
               f'({Y}/{X} of contradicted loci); untested nearest gene {Z}/259. F1 (Holm over T_A, T_B) '
               f'{"rejects" if F1["rejects"] else "does not reject"}; S6 {"triggered" if NC1["S6_triggered"] else "not triggered"}; '
               f'PC-A {"passed" if pcA["passed"] else "failed"}, PC-B passed. {catalog}.')
    if args.debug_null:
        verdict = 'DEBUG-NULL (not a result): ' + verdict
    S['verdict'] = verdict
    S['status'] = 'completed'
    S['finished_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    run.write_summary()
    LT.to_csv(os.path.join(args.out, 'locus_table.tsv'), sep='\t', index=False)
    SENS.to_csv(os.path.join(args.out, 'sensitivities.tsv'), sep='\t', index=False)
    NL.to_csv(os.path.join(args.out, 'named_loci.tsv'), sep='\t', index=False)
    write_result(args.out, S, E, F1, NC1, pcA, pcB, SENS, NL, lines, verdict, args.debug_null)
    print('\nVERDICT:', verdict)


def fmt(x, d=3):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return 'NA'
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return f'{x:.{d}g}' if abs(x) < 1e-3 else f'{x:.{d}f}'


def write_result(out, S, E, F1, NC1, pcA, pcB, SENS, NL, lines, verdict, debug):
    w = []
    if debug:
        w.append('# DEBUG-NULL RUN, NOT A RESULT (anchors displaced, gene states shuffled)\n')
    w.append('# RESULT_H74: nearest protein-coding gene versus colocalized genes at tier-1/2 liver-trait loci\n')
    w.append(f'Status: completed. Prespec sha256 {PRESPEC_SHA}; seed {S["seed"]}; B = {F1["B"]}; NC-1 draws {NC1["draws"]}. '
             'Candidate only; not adopted. COLOC release susie_coloc_r2_20260924T053244Z.\n')
    w.append(f'## Verdict (computed in code by the prespec Reading rule)\n\n{verdict}\n')
    w.append('## Sentences allowed by the Reading\n')
    for l in lines:
        w.append(f'- {l}')
    w.append('\n## Primary readouts (estimates, no test; unit = locus, n = 259)\n')
    lo, hi = E['E1_wilson']; cl, ch = E['E1_clopper_pearson']
    w.append(f'- E1 (nearest protein-coding gene not colocalized): {E["E1_k"]}/259 = {E["E1"]:.3f}, Wilson 95% CI '
             f'{lo:.3f} to {hi:.3f} (Clopper-Pearson {cl:.3f} to {ch:.3f}).')
    bl, bh = E['E2_wilson']; cl2, ch2 = E['E2_of_contradicted_wilson']
    w.append(f'- E2 (nearest gene is a bystander, tested, PP.H4 <= 0.5, PP.H3 > 0.8): {E["E2_k"]}/259 = {E["E2"]:.3f}, '
             f'Wilson {bl:.3f} to {bh:.3f}; of contradicted loci {E["E2_k"]}/{E["E1_k"]} = {fmt(E["E2_of_contradicted"])}, '
             f'Wilson {fmt(cl2)} to {fmt(ch2)}.')
    sp = E['split']
    w.append(f'- Split of the 259 nearest genes: colocalized {sp["colocalized"]}, bystander {sp["bystander"]}, tested other '
             f'{sp["tested_other"]}, untested {sp["untested"]} (i tier-3/4 SuSiE {sp["untested_subclasses"]["i_tier34_susie"]}, '
             f'ii SuSiE-untestable {sp["untested_subclasses"]["ii_susie_untestable"]}, iii ABF only '
             f'{sp["untested_subclasses"]["iii_abf_only"]}, iv no row {sp["untested_subclasses"]["iv_no_row"]}). '
             f'Nearest genes flagged for a tier-1/2 SuSiE-untestable row: {S["D5"]["nearest_flagged_tier12_untestable"]}. '
             f'Ties at the nearest protein-coding TSS: {len(E["nearest_ties"])}.')
    w.append('\n## Calibration family F1 = {T_A, T_B}, one-sided, Holm over 2, alpha 0.05 (248 informative loci)\n')
    w.append(f'- T_A (loci whose nearest tested gene is colocalized): {fmt(F1["T_A"])} of {F1["n_informative"]} '
             f'(rate {F1["T_A_rate"]:.3f}); exact-null mean {F1["T_A_null_mean"]:.2f} (rate {F1["T_A_null_rate"]:.4f}), '
             f'SD {F1["T_A_null_sd"]:.2f}; exact p = {fmt(F1["T_A_exact_p"])}; permutation p = {fmt(F1["T_A_perm_p"])} '
             f'(S5 tolerance {fmt(F1["S5_tolerance"])}: {"pass" if F1["S5_passed"] else "FAIL"}); Holm p = {fmt(F1["holm_adj"]["T_A"])}.')
    w.append(f'- T_B (mean normalized distance rank of colocalized genes; null 0.5): {F1["T_B"]:.4f}; permutation null mean '
             f'{F1["T_B_null_mean"]:.4f}; p = {fmt(F1["T_B_perm_p"])} (floor 9.999e-05); Holm p = {fmt(F1["holm_adj"]["T_B"])}.')
    w.append(f'- F1 rejects: {F1["rejects"]}. Every other p value in this run is nominal and descriptive.')
    w.append('\n## Controls\n')
    w.append(f'- PC-A (1p13, not blind): anchor {pcA["anchor"]} (top_snp {pcA["anchor_top_snp"]}); nearest protein-coding '
             f'gene {"/".join(pcA["nearest_pc"])} at {pcA["nearest_dist"]} bp, class {pcA["nearest_class"]}; C = '
             f'{", ".join(pcA["C"])}. Prediction agreement: {"passed" if pcA["passed"] else "FAILED"}.')
    w.append(f'- PC-B (synthetic labels): (i) T_A = {fmt(pcB["i_T_A"])}, T_B = {pcB["i_T_B"]:.4f}, exact p = '
             f'{fmt(pcB["i_exact_p"])}; (ii) T_A = {fmt(pcB["ii_T_A"])}; (iii) T_A = {fmt(pcB["iii_T_A"])} vs null mean '
             f'{pcB["iii_null_mean"]:.2f} (z = {pcB["iii_z"]:.2f}). Passed: {pcB["passed"]}.')
    w.append(f'- NC-1 (displaced anchor, {NC1["draws"]} draws): mean rate {NC1["mean_rate"]:.4f} (2.5 to 97.5%: '
             f'{NC1["q025"]:.4f} to {NC1["q975"]:.4f}); predicted 0.15 to 0.30: '
             f'{"inside" if NC1["mean_in_predicted_range"] else "OUTSIDE"}; observed T_A rate {NC1["observed_rate"]:.4f}; '
             f'NC-1 nominal p = {fmt(NC1["nominal_p"])}. S6 triggered: {NC1["S6_triggered"]}.')
    w.append('\n## Gates and stop rules\n')
    w.append(f'- S1 input hashes: pass. G1/S2: {S["G1"]["n_loci"]} loci, {S["G1"]["n_informative"]} informative: pass. '
             f'S8 anchors: pass. S7: {S["S7"]["n_loci_no_pc_within_2Mb"]} loci without a protein-coding TSS within 2 Mb: pass. '
             f'S4: not triggered. S3: pass. S5: {"pass" if F1["S5_passed"] else "FAIL"}. S6: '
             f'{"triggered" if NC1["S6_triggered"] else "not triggered"}.')
    w.append('\n## Diagnostics\n')
    for k in ['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D8', 'D9']:
        w.append(f'- {k}: `{json.dumps(S[k], default=float)}`')
    w.append('\nD7 named loci:\n')
    try:
        w.append(NL.to_markdown(index=False))
    except ImportError:
        w.append('```\n' + NL.to_string(index=False) + '\n```')
    w.append('\n## Sensitivities (descriptive, nominal, outside F1)\n')
    cols = ['analysis', 'n_loci', 'E1_k', 'E1', 'E1_wilson_lo', 'E1_wilson_hi', 'E2_k', 'untested_k', 'n_informative',
            'T_A', 'T_A_null_mean', 'T_A_exact_p_nominal', 'T_B', 'note']
    w.append('| ' + ' | '.join(cols) + ' |')
    w.append('|' + '---|' * len(cols))
    intcols = {'n_loci', 'E1_k', 'E2_k', 'untested_k', 'n_informative'}

    def cell(c, v):
        if v is None or (isinstance(v, (float, np.floating)) and math.isnan(v)):
            return 'NA'
        if c in intcols:
            return str(int(v))
        if isinstance(v, (float, np.floating)):
            return fmt(float(v))
        return str(v)
    for _, r in SENS.iterrows():
        w.append('| ' + ' | '.join(cell(c, r[c]) for c in cols) + ' |')
    w.append('\n## Files\n\nsummary.json, locus_table.tsv, sensitivities.tsv, named_loci.tsv, null_T_A_T_B.npz, '
             'input_sha256_check.tsv, h74.log, pip_freeze.txt, script_sha256.tsv.')
    with open(os.path.join(out, 'RESULT_H74.md'), 'w') as f:
        f.write('\n'.join(w) + '\n')


if __name__ == '__main__':
    main()
