"""Feasibility gate for L08-03: CLPP between GWAS SuSiE credible sets and GTEx v10 liver
QTL credible sets (eQTL, sQTL, apaQTL), plus small-RNA lead-in-CS overlap, with a
within-chromosome circular-shift null of QTL positions.

Usage: python clpp_gate.py <n_shifts> <outdir>
Inputs are read-only. Seed 20260923.
"""
import sys, os, json, time
import numpy as np
import pandas as pd

NSHIFT = int(sys.argv[1]) if len(sys.argv) > 1 else 100
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(OUT, exist_ok=True)
SEED = 20260923
THR = 0.1
R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
GW = R + '/GWAS/finemapping/results/'
S = R + '/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/'
G = S + 'lens-L08/gtex_v10/'
ME = S + 'gate-L08-03/'
COMP = str.maketrans('ACGT', 'TGCA')

def akey(a, b):
    return ['|'.join(sorted((x, y))) for x, y in zip(a, b)]

t0 = time.time()
# ---------------- GWAS credible sets (hg19 -> hg38) ----------------
c = pd.read_csv(GW + 'credible_sets.csv', dtype={'locus': str}, low_memory=False)
c = c[c.susie_cs > 0].copy()
c['cs_uid'] = c.study + '|' + c.locus + '|' + c.susie_cs.astype(int).astype(str)
c = c.sort_values('susie_pip', ascending=False).drop_duplicates(['cs_uid', 'variant_id'])
lift = pd.read_csv(ME + 'sites_hg19_to_hg38.tsv', sep='\t')
lift['chromosome'] = lift.chr19.str.replace('chr', '').astype(int)
c = c.merge(lift[['chromosome', 'pos19', 'chr38', 'pos38', 'strand38']],
            left_on=['chromosome', 'position'], right_on=['chromosome', 'pos19'], how='inner')
a1 = c.allele1.astype(str).str.upper(); a2 = c.allele2.astype(str).str.upper()
minus = c.strand38 == '-'
a1 = a1.where(~minus, a1.str.translate(COMP)); a2 = a2.where(~minus, a2.str.translate(COMP))
c['al'] = akey(a1, a2)
c['chr'] = c.chr38
arch = pd.read_csv(GW + 'candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_pip_architecture.tsv',
                   sep='\t', dtype={'locus': str})
arch = arch[arch.pip_sum_gate_passed.astype(str).str.lower() == 'true']
arch_uid = set(arch.study + '|' + arch.locus + '|' + arch.susie_cs.astype(int).astype(str))
c['in_resource'] = c.cs_uid.isin(arch_uid)
print('GWAS CS rows', len(c), 'CS', c.cs_uid.nunique(), 'resource gated CS', c[c.in_resource].cs_uid.nunique(),
      'of', len(arch_uid), 'studies', c.study.nunique())
cs_ids = pd.Index(sorted(c.cs_uid.unique()))
c['csi'] = cs_ids.get_indexer(c.cs_uid)
cs_res = np.zeros(len(cs_ids), bool); cs_res[c.loc[c.in_resource, 'csi'].unique()] = True
# CS lead position (hg38) for locus clustering
lead = c.sort_values('susie_pip', ascending=False).drop_duplicates('csi').set_index('csi')
cs_chr = lead.reindex(range(len(cs_ids))).chr.values
cs_pos = lead.reindex(range(len(cs_ids))).pos38.values

chroms = pd.read_csv(ME + 'hg38_chrom_sizes.tsv', sep='\t', header=None, names=['chr', 'len'])
chroms = chroms[chroms.chr != 'chrX']
CH = {k: i for i, k in enumerate(chroms.chr)}
CLEN = chroms.len.values.astype(np.int64)
c['ci'] = c.chr.map(CH).astype(int)
SCALE = np.int64(10**9)
c['skey'] = c.ci.values.astype(np.int64) * SCALE + c.pos38.values.astype(np.int64)
c = c.sort_values('skey').reset_index(drop=True)
g_skey = c.skey.values; g_csi = c.csi.values; g_pip = c.susie_pip.values; g_al = c.al.values

# COLOC r2 SuSiE-positive genes (Resource eQTL colocalization)
t12 = pd.read_csv(GW + 'susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv')
coloc_pos = set(t12.loc[t12.susie_state_t12 == 'positive', 'ensembl'])
print('COLOC r2 tier12 SuSiE-positive genes', len(coloc_pos))

def load_qtl(cls):
    a = pd.read_parquet(G + 'Liver.v10.%s.SuSiE_summary.parquet' % cls)
    p = a.variant_id.str.split('_', expand=True)
    a['chr'] = p[0]; a['pos'] = p[1].astype(np.int64); a['al'] = akey(p[2], p[3])
    a = a[a.chr.isin(CH)].copy()
    a['ci'] = a.chr.map(CH).astype(int)
    a['ens'] = a.phenotype_id.str.extract(r'(ENSG\d+)')[0]
    a['qcs'] = a.phenotype_id + '#' + a.cs_id.astype(str)
    return a.reset_index(drop=True)

def clpp_table(q, qpos, use_allele):
    """Return per (GWAS csi, QTL qcs) CLPP using QTL positions qpos (array aligned to q)."""
    key = q.ci.values.astype(np.int64) * SCALE + qpos
    lo = np.searchsorted(g_skey, key, 'left'); hi = np.searchsorted(g_skey, key, 'right')
    n = hi - lo
    m = n > 0
    if not m.any():
        return pd.DataFrame({'csi': np.array([], dtype=np.int64), 'qcs': np.array([], dtype=object),
                             'ens': np.array([], dtype=object), 'clpp': np.array([], dtype=float)})
    qi = np.repeat(np.nonzero(m)[0], n[m])
    gi = np.concatenate([np.arange(l, h) for l, h in zip(lo[m], hi[m])])
    if use_allele:
        ok = q.al.values[qi] == g_al[gi]
        qi, gi = qi[ok], gi[ok]
    d = pd.DataFrame({'csi': g_csi[gi], 'qcs': q.qcs.values[qi], 'ens': q.ens.values[qi],
                      'clpp': g_pip[gi] * q.pip.values[qi]})
    return d.groupby(['csi', 'qcs', 'ens'], as_index=False).clpp.sum()

def stats(t, eq_hit_cs, cls):
    """t: per pair CLPP. Return counts of GWAS CS with max CLPP >= THR."""
    out = {}
    h = t[t.clpp >= THR]
    h = h.assign(csi=h.csi.astype(np.int64))
    for uni, mask in [('all', np.ones(len(cs_ids), bool)), ('resource', cs_res)]:
        hh = h[mask[h.csi.values]]
        out[uni + '_cs'] = hh.csi.nunique()
        out[uni + '_genes'] = hh.ens.nunique()
        if cls != 'eQTLs':
            b = hh[~hh.csi.isin(eq_hit_cs) & ~hh.ens.isin(coloc_pos)]
            out[uni + '_cs_beyond_eqtl'] = b.csi.nunique()
            out[uni + '_genes_beyond_eqtl'] = b.ens.nunique()
            # distinct 1-Mb loci among beyond-eQTL hits
            if len(b):
                u = np.unique(b.csi.values)
                lk = pd.Series(cs_chr[u]).astype(str) + ':' + pd.Series((cs_pos[u] // 1_000_000).astype(int)).astype(str)
                out[uni + '_loci1Mb_beyond_eqtl'] = lk.nunique()
            else:
                out[uni + '_loci1Mb_beyond_eqtl'] = 0
    return out

rng = np.random.default_rng(SEED)
offsets = rng.integers(0, CLEN, size=(NSHIFT, len(CLEN)))  # same offsets reused across classes

def shifted(q, b):
    off = offsets[b][q.ci.values]
    L = CLEN[q.ci.values]
    return (q.pos.values - 1 + off) % L + 1

rng_local = np.random.default_rng(SEED + 1)
def shifted_local(q, b):
    codes, uniq = pd.factorize(q.qcs)
    mag = rng_local.integers(100_000, 2_000_001, size=len(uniq))
    sgn = rng_local.choice([-1, 1], size=len(uniq))
    off = (mag * sgn)[codes]
    L = CLEN[q.ci.values]
    return (q.pos.values - 1 + off) % L + 1

def draw_pos(q, b, mode):
    return shifted_local(q, b) if mode.startswith('local') else shifted(q, b)

res = {}; nulls = []
eq = load_qtl('eQTLs')
for cls in ['eQTLs', 'sQTLs', 'apaQTLs']:
    q = eq if cls == 'eQTLs' else load_qtl(cls)
    obs = clpp_table(q, q.pos.values, True)
    if cls == 'eQTLs':
        eq_obs = obs
        eq_hit_cs = set(obs.loc[obs.clpp >= THR, 'csi'])
        r = obs[obs.ens == 'ENSG00000069667'].copy()
        r['cs_uid'] = cs_ids[r.csi.values]; r['resource'] = cs_res[r.csi.values]
        r.sort_values('clpp', ascending=False).to_csv(OUT + '/rora_eqtl_clpp_by_gwas_cs.tsv', sep='\t', index=False)
        print('RORA eQTL: GWAS CS with any shared variant', len(r), 'CLPP>=0.1', (r.clpp >= THR).sum(),
              'max CLPP', round(r.clpp.max(), 4) if len(r) else None)
        print(r.sort_values('clpp', ascending=False).head(8)[['cs_uid', 'qcs', 'clpp', 'resource']].to_string())
    so = stats(obs, eq_hit_cs, cls)
    o2 = obs.copy(); o2['cs_uid'] = cs_ids[o2.csi.values]; o2['resource'] = cs_res[o2.csi.values]; o2['class'] = cls
    o2.to_csv(OUT + '/observed_pairs_%s.tsv' % cls, sep='\t', index=False)
    res[cls] = {'observed': so, 'n_qtl_cs': int(q.qcs.nunique()), 'n_qtl_variants': int(q.variant_id.nunique())}
    print(cls, 'observed', so)
    for mode in ['pos_allele', 'pos_only', 'local2Mb_pos_only']:
        for b in range(NSHIFT):
            t = clpp_table(q, draw_pos(q, b, mode), mode == 'pos_allele')
            s = stats(t, eq_hit_cs, cls); s.update(cls=cls, mode=mode, draw=b); nulls.append(s)
    print(cls, 'done', round(time.time() - t0, 1), 's')

# small RNA: lead variant (q<=0.05) inside a GWAS CS
eg = pd.read_csv(G + 'GTEx_Analysis_v10_smallRNA_eQTL/Liver.v10.smallRNA.cis_egenes.txt.gz', sep='\t')
eg = eg[(eg.qval <= 0.05) & eg.chr.isin(CH)].copy()
eg['ci'] = eg.chr.map(CH).astype(int); eg['pos'] = eg.variant_pos.astype(np.int64)
eg['al'] = akey(eg.ref.astype(str), eg.alt.astype(str))
eg['pip'] = 1.0; eg['qcs'] = eg.phenotype_id; eg['ens'] = eg.gene_name.astype(str)
eg = eg.reset_index(drop=True)
obs = clpp_table(eg, eg.pos.values, True)
def srna_stats(t):
    return {'all_cs': t.csi.nunique(), 'resource_cs': t[cs_res[t.csi.values]].csi.nunique() if len(t) else 0,
            'all_srna': t.qcs.nunique()}
so = srna_stats(obs); res['smallRNA_lead_in_cs'] = {'observed': so, 'n_egenes': int(len(eg))}
o2 = obs.copy(); o2['cs_uid'] = cs_ids[o2.csi.values.astype(int)] if len(o2) else []
o2.to_csv(OUT + '/observed_pairs_smallRNA.tsv', sep='\t', index=False)
print('smallRNA observed', so)
for mode in ['pos_allele', 'pos_only', 'local2Mb_pos_only']:
    for b in range(NSHIFT):
        t = clpp_table(eg, draw_pos(eg, b, mode), mode == 'pos_allele')
        s = srna_stats(t) if len(t) else {'all_cs': 0, 'resource_cs': 0, 'all_srna': 0}
        s.update(cls='smallRNA', mode=mode, draw=b); nulls.append(s)

nd = pd.DataFrame(nulls)
nd.to_csv(OUT + '/shift_null_draws.tsv', sep='\t', index=False)
summ = []
for (cls, mode), d in nd.groupby(['cls', 'mode']):
    key = 'eQTLs' if cls == 'eQTLs' else cls
    ob = res[key if key in res else 'smallRNA_lead_in_cs']['observed'] if cls != 'smallRNA' else res['smallRNA_lead_in_cs']['observed']
    for col in [k for k in d.columns if k not in ('cls', 'mode', 'draw')]:
        if col not in ob or d[col].isna().all():
            continue
        x = d[col].values.astype(float)
        p = (1 + (x >= ob[col]).sum()) / (1 + len(x))
        summ.append(dict(cls=cls, mode=mode, stat=col, observed=ob[col], null_mean=x.mean(), null_sd=x.std(ddof=1),
                         null_q95=np.quantile(x, .95), null_max=x.max(), p_upper=p, n_draws=len(x)))
sm = pd.DataFrame(summ)
sm.to_csv(OUT + '/shift_null_summary.tsv', sep='\t', index=False)
print(sm.to_string())
json.dump(res, open(OUT + '/observed.json', 'w'), indent=1, default=int)
print('total seconds', round(time.time() - t0, 1))
