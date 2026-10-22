# H72 prespec counts. Reads NO beta, tstat or pval column of any sex-stratified file.
# Counts: anchored genes by best UKBB-enzyme PP.H4 tier; allele resolution of top SNPs and
# clumped leads against the Neale irnt female/male files; chrX coding and scale diagnostics
# from minor_AF, n and se only.
import csv, gzip, math, statistics as st, sys
from multiprocessing import Pool
import pandas as pd

R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
W = R + '/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72'
L = R + '/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L06'
S = R + '/GWAS/MR_Data/PanUKBB/sex_stratified'
CAT = R + '/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
COLOC = R + '/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv'
PC = {'ALT': '30620', 'AST': '30650', 'GGT': '30730'}
PAR = [(60001, 2699520), (154931044, 155260560)]  # hg19 PAR1, PAR2

def isPAR(p):
    return any(a <= p <= b for a, b in PAR)

def fin(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None

def main():
    out = open(W + '/count_n.out', 'w')
    def p(*a):
        print(*a, file=out, flush=True)
    cat = pd.read_parquet(CAT)
    anc = cat[cat.primary_evidence_class.isin(['genetically_anchored', 'concordant'])]
    p('catalog rows', len(cat), 'anchored+concordant', len(anc), dict(anc.primary_evidence_class.value_counts()))
    ens = set(anc.ensembl_id)
    co = pd.read_csv(COLOC, usecols=['gwas_name', 'gene', 'ensembl', 'chr', 'PP.H4.abf', 'PP.H4.susie', 'method', 'top_snp', 'top_snp_PP'],
                     dtype=str)
    co = co[co.gwas_name.isin(['UKBB_ALT', 'UKBB_AST', 'UKBB_GGT']) & co.ensembl.isin(ens)].copy()
    p('UKBB enzyme rows for anchored genes', len(co), 'genes', co.ensembl.nunique(),
      'duplicate gene x trait rows', int(co.duplicated(['gwas_name', 'ensembl']).sum()))
    co['pp_s'] = co['PP.H4.susie'].map(fin)
    co['pp_a'] = co['PP.H4.abf'].map(fin)
    co['pp'] = co.pp_s.where(co.pp_s.notna(), co.pp_a)
    co['pp_src'] = co.pp_s.notna().map({True: 'susie', False: 'abf'})
    p('method counts', dict(co.method.value_counts()))
    co = co[co.pp.notna()].sort_values(['ensembl', 'pp', 'gwas_name'], ascending=[True, False, True])
    best = co.drop_duplicates('ensembl').copy()
    best['trait'] = best.gwas_name.str.split('_').str[1]
    best['c'] = best.top_snp.str.split(':').str[0]
    best['pos'] = best.top_snp.str.split(':').str[1].astype(int)
    p('genes with a finite UKBB enzyme PP.H4', len(best), '; missing', len(ens) - len(best))
    for lo, hi, lab in [(0.8, 9, '>=0.8'), (0.5, 0.8, '[0.5,0.8)'), (-1, 0.5, '<0.5')]:
        b = best[(best.pp >= lo) & (best.pp < hi)]
        p(f'tier {lab}: genes {len(b)} by trait {dict(b.trait.value_counts())} by pp source {dict(b.pp_src.value_counts())}')
    p('best top_snp chromosomes', sorted(best.c.unique(), key=lambda x: (len(x), x)))
    best[['ensembl', 'gene', 'gwas_name', 'trait', 'pp', 'pp_src', 'method', 'top_snp', 'top_snp_PP']].to_csv(W + '/best_enzyme_row.tsv', sep='\t', index=False)

    # alleles at the top SNP from the COLOC source GWAS (sex-combined UKBB, hg19)
    need = {t: set(zip(best[best.trait == t].c, best[best.trait == t].pos)) for t in PC}
    alle = {}
    for t in PC:
        with open(f'{R}/GWAS/finemapping/data/sumstats/UKBB_{t}_reformatted_hg19.tsv') as f:
            next(f)
            for l in f:
                a = l.split('\t', 4)
                k = (a[0], int(a[1]))
                if k in need[t]:
                    alle.setdefault((t,) + k, []).append(frozenset((a[2], a[3])))
    # clumped leads (no header): chr pos a1 a2 beta se p
    leads = {}
    for t in PC:
        for l in open(f'{L}/ukbb_{t.lower()}_clumped500kb.tsv'):
            a = l.rstrip('\n').split('\t')
            leads.setdefault(t, []).append((a[0], int(a[1]), frozenset((a[2], a[3]))))
    jobs = [(t, sx, {k[1:]: v for k, v in alle.items() if k[0] == t}, leads[t]) for t in PC for sx in ('female', 'male')]
    with Pool(6) as pool:
        res = pool.map(scan, jobs)
    R_ = {(r['t'], r['sx']): r for r in res}
    for t in PC:
        F, M = R_[(t, 'female')], R_[(t, 'male')]
        p(f'--- {t}: rows F {F["rows"]} M {M["rows"]}; n_complete F {F["n"]} M {M["n"]}')
        # gene arm
        bt = best[best.trait == t]
        st_ = {'no_source_row': 0, 'source_multi': 0, 'no_neale_match': 0, 'neale_multi_match': 0, 'lowconf_or_maf': 0, 'usable': 0}
        usable = {'>=0.8': 0, '[0.5,0.8)': 0, '<0.5': 0}
        for r in bt.itertuples():
            k = (r.c, r.pos)
            src = alle.get((t,) + k)
            if not src:
                st_['no_source_row'] += 1; continue
            if len(set(src)) > 1:
                st_['source_multi'] += 1; continue
            m = [v for v in F['site'].get(k, []) if v[0] == src[0]]
            mm = [v for v in M['site'].get(k, []) if v[0] == src[0]]
            if not m or not mm:
                st_['no_neale_match'] += 1; continue
            if len(m) > 1:
                st_['neale_multi_match'] += 1; continue
            if m[0][1] or mm[0][1] or m[0][2] < 0.01 or mm[0][2] < 0.01:
                st_['lowconf_or_maf'] += 1; continue
            st_['usable'] += 1
            tier = '>=0.8' if r.pp >= 0.8 else ('[0.5,0.8)' if r.pp >= 0.5 else '<0.5')
            usable[tier] += 1
        p(f'{t} gene arm (best trait {t}, {len(bt)} genes): {st_}; usable by tier {usable}')
        p(f'{t} locus arm: clumped leads {len(leads[t])}; F status {F["lead_status"]}; M status {M["lead_status"]}')
        p(f'{t} chrX F: rows {F["x_rows"]} nonPAR {F["x_nonpar"]} PAR {F["x_par"]}')
        p(f'{t} chrX diag (nonPAR, lowconf false both, female MAF>=0.05): n {F["x_pairs"]}; median male/female minor_AF {F["x_maf_ratio"]}; '
          f'median se_M/se_F chrX {F["x_se_ratio"]}; median se_M/se_F autosomes {F["a_se_ratio"]}')
        p(f'{t} implied irnt trait SD = median se*sqrt(2*maf*(1-maf)*n), autosomes MAF 0.05-0.5 every 50th row: F {F["sd"]} M {M["sd"]}')
    out.close()

def scan(job):
    t, sx, alle, leads = job
    fn = f'{S}/{PC[t]}_irnt.gwas.imputed_v3.{sx}.varorder.tsv.bgz'
    site = {}
    need = set(alle)
    lead_need = {(c, p) for c, p, _ in leads}
    lead_hit = {}
    rows = 0; n = None; xr = xn = xp = 0; sdv = []
    keepx = sx == 'female'
    # chrX / autosome se and MAF for coding diagnostic, female and male read in lockstep by the female worker
    other = gzip.open(f'{S}/{PC[t]}_irnt.gwas.imputed_v3.male.varorder.tsv.bgz', 'rt') if keepx else None
    if other: next(other)
    xm, xs, asr = [], [], []
    with gzip.open(fn, 'rt') as f:
        next(f)
        for i, l in enumerate(f):
            rows += 1
            a = l.split('\t', 9)  # columns 0..8 only; beta (7) never parsed
            if other is not None:
                b = other.readline().split('\t', 9)
                assert a[0] == b[0], (a[0], b[0])
            v = a[0].split(':')
            c, pos = v[0], int(v[1])
            if n is None and a[4] not in ('', 'NaN'):
                n = a[4]
            low = a[3] == 'true'
            maf = fin(a[2]) or 0.0
            if (c, pos) in need or (c, pos) in lead_need:
                rec = (frozenset((v[2], v[3])), low, maf)
                if (c, pos) in need:
                    site.setdefault((c, pos), []).append(rec)
                if (c, pos) in lead_need:
                    lead_hit.setdefault((c, pos), []).append(rec)
            if c == 'X':
                xr += 1
                if isPAR(pos):
                    xp += 1
                else:
                    xn += 1
                    if other is not None and not low and b[3] != 'true' and maf >= 0.05:
                        sF, sM, mM = fin(a[8]), fin(b[8]), fin(b[2])
                        if sF and sM and mM:
                            xm.append(mM / maf); xs.append(sM / sF)
            elif i % 50 == 0 and not low and 0.05 <= maf <= 0.5:
                s_ = fin(a[8]); nn = fin(a[4])
                if s_ and nn:
                    sdv.append(s_ * math.sqrt(2 * maf * (1 - maf) * nn))
                    if other is not None and b[3] != 'true':
                        sM = fin(b[8])
                        if sM:
                            asr.append(sM / s_)
    ls = {'not_found': 0, 'allele_match': 0, 'multi_match': 0, 'lowconf_or_maf<0.01': 0}
    for c, p_, al in leads:
        m = [r for r in lead_hit.get((c, p_), []) if r[0] == al]
        if not m:
            ls['not_found'] += 1
        elif len(m) > 1:
            ls['multi_match'] += 1
        elif m[0][1] or m[0][2] < 0.01:
            ls['lowconf_or_maf<0.01'] += 1
        else:
            ls['allele_match'] += 1
    med = lambda x: round(st.median(x), 4) if x else None
    return dict(t=t, sx=sx, rows=rows, n=n, site=site, lead_status=ls, x_rows=xr, x_nonpar=xn, x_par=xp,
                x_pairs=len(xm), x_maf_ratio=med(xm), x_se_ratio=med(xs), a_se_ratio=med(asr), sd=med(sdv))

if __name__ == '__main__':
    main()
