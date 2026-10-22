# H72 audit: independent recount of every n, hash and count stated in PRESPEC_H72.md.
# Never reads the female or male beta, tstat or pval column of the ALT, AST or GGT files
# (columns 7, 9, 10). Reads variant, minor_AF, low_confidence_variant, n_complete_samples, se
# from female/male files and beta, se, pval from the both_sexes files. Urate (30880) betas are
# read for PC1 only (already on disk in gate/NG11/positive_control.txt).
import gzip, hashlib, math, os, statistics as st, sys
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy.stats import norm, chi2

SEED = 20260923
R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
RR = R + '/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z'
R4 = RR + '/round4-20261015T112148Z'
A = R4 + '/audit/H72'
W = R4 + '/write/H72'
L = RR + '/scratch/lens-L06'
S = R + '/GWAS/MR_Data/PanUKBB/sex_stratified'
CAT = R + '/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
COLOC = R + '/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv'
FIELD = {'ALT': '30620', 'AST': '30650', 'GGT': '30730'}
PAR = [(60001, 2699520), (154931044, 155260560)]

out = open(A + '/audit_recount.out', 'w')
def p(*a):
    print(*a, file=out, flush=True)

# ---------------------------------------------------------------- 1. hashes and bytes
STATED = [
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

def sha(path):
    h = hashlib.sha256()
    with open(R + '/' + path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return path, os.path.getsize(R + '/' + path), h.hexdigest()

def fin(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None

# ---------------------------------------------------------------- 3. Neale scan per trait
def scan(job):
    t, need, lead_pos = job
    fF, fM, fB = [gzip.open(f'{S}/{FIELD[t]}_irnt.gwas.imputed_v3.{s}.varorder.tsv.bgz', 'rt') for s in ('female', 'male', 'both_sexes')]
    hF, hM, hB = next(fF), next(fM), next(fB)
    assert hF == hM == hB
    rows = 0; nF = nM = nB = None
    xr = xp = 0; xm = []; xs = []; asr = []; sdF = []; sdM = []
    site = {}; pool = []; elig = 0
    xB_gws = []
    for i, (lf, lm, lb) in enumerate(zip(fF, fM, fB)):
        rows += 1
        a = lf.split('\t', 9); b = lm.split('\t', 9); c = lb.rstrip('\n').split('\t')
        if not (a[0] == b[0] == c[0]):
            raise SystemExit(f'{t} variant mismatch at row {i}: {a[0]} {b[0]} {c[0]}')
        if nF is None and a[4] not in ('', 'NaN'): nF = a[4]
        if nM is None and b[4] not in ('', 'NaN'): nM = b[4]
        if nB is None and c[4] not in ('', 'NaN'): nB = c[4]
        v = a[0].split(':'); ch = v[0]; pos = int(v[1])
        lowF = a[3] == 'true'; lowM = b[3] == 'true'
        mafF = fin(a[2]) or 0.0; mafM = fin(b[2]) or 0.0
        k = (ch, pos)
        if k in need or k in lead_pos:
            site.setdefault(k, []).append((frozenset((v[2], v[3])), lowF, lowM, mafF, mafM, fin(a[8]), fin(b[8]), fin(c[7]), fin(c[8])))
        if ch == 'X':
            xr += 1
            if any(lo <= pos <= hi for lo, hi in PAR):
                xp += 1; continue
            if not lowF and not lowM and mafF >= 0.05:
                sF, sM = fin(a[8]), fin(b[8])
                if sF and sM and mafM:
                    xm.append(mafM / mafF); xs.append(sM / sF)
            pB = fin(c[10])
            if pB is not None and pB < 5e-8 and c[3] != 'true' and len(v[2]) == 1 and len(v[3]) == 1 and min(mafF, mafM) >= 0.01:
                xB_gws.append((pos, pB))
            continue
        if i % 50 == 0:
            if not lowF and 0.05 <= mafF <= 0.5:
                s_, nn = fin(a[8]), fin(a[4])
                if s_ and nn:
                    sdF.append(s_ * math.sqrt(2 * mafF * (1 - mafF) * nn))
                    if not lowM:
                        sM = fin(b[8])
                        if sM: asr.append(sM / s_)
            if not lowM and 0.05 <= mafM <= 0.5:
                s_, nn = fin(b[8]), fin(b[4])
                if s_ and nn: sdM.append(s_ * math.sqrt(2 * mafM * (1 - mafM) * nn))
        # NC1 feasibility pool: eligible = SNV, low_conf false both, MAF>=0.05 both; every 20th eligible;
        # proxy selection on both_sexes p > 0.5 (the run uses the female/male IVW p, not read here)
        if not lowF and not lowM and mafF >= 0.05 and mafM >= 0.05 and len(v[2]) == 1 and len(v[3]) == 1:
            elig += 1
            if elig % 20 == 0:
                pB = fin(c[10])
                if pB is not None and pB > 0.5:
                    pool.append((ch, pos))
    # chrX sex-blind proxy leads from both_sexes (indicative only), 500 kb distance clumping by p
    xB_gws.sort(key=lambda x: (x[1], x[0])); xl = []
    for pos, pv in xB_gws:
        if all(abs(pos - q) > 500000 for q in xl): xl.append(pos)
    med = lambda x: round(st.median(x), 4) if x else None
    return dict(t=t, rows=rows, nF=nF, nM=nM, nB=nB, x_rows=xr, x_par=xp, x_pairs=len(xm), x_maf_ratio=med(xm),
                x_se_ratio=med(xs), a_se_ratio=med(asr), sdF=med(sdF), sdM=med(sdM), site=site, pool=pool,
                elig=elig, xB_gws=len(xB_gws), xB_leads=len(xl))

def pc1():
    D = R + '/data/public_r4/neale_30880'
    F = gzip.open(f'{D}/30880_irnt.gwas.imputed_v3.female.varorder.tsv.bgz', 'rt'); M = gzip.open(f'{D}/30880_irnt.gwas.imputed_v3.male.varorder.tsv.bgz', 'rt')
    next(F); next(M); best = None; seen4 = False
    for lf, lm in zip(F, M):
        a = lf.rstrip('\n').split('\t'); b = lm.rstrip('\n').split('\t')
        ch, pos = a[0].split(':')[:2]
        if ch != '4':
            if seen4: break
            continue
        seen4 = True
        pos = int(pos)
        if not (9000000 <= pos <= 10500000): continue
        assert a[0] == b[0]
        if a[3] != 'false' or b[3] != 'false': continue
        bF, sF, bM, sM = fin(a[7]), fin(a[8]), fin(b[7]), fin(b[8])
        if None in (bF, sF, bM, sM): continue
        zc = (bF / sF**2 + bM / sM**2) / math.sqrt(1 / sF**2 + 1 / sM**2)
        if best is None or abs(zc) > abs(best[0]): best = (zc, a[0], bF, sF, bM, sM)
    zc, v, bF, sF, bM, sM = best
    return dict(lead=v, zc=zc, bF=bF, sF=sF, bM=bM, sM=sM, zdiff=(bF - bM) / math.sqrt(sF**2 + sM**2))

def main():
    # 1. hashes
    with Pool(8) as pool:
        hs = pool.map(sha, [x[0] for x in STATED])
    bad = 0
    with open(A + '/sha256_recomputed.tsv', 'w') as f:
        f.write('path\tbytes\tbytes_stated\tsha256\tsha256_stated\tmatch\n')
        for (path, nbytes, h), (_, b0, h0) in zip(hs, STATED):
            ok = (h == h0) and (b0 is None or b0 == nbytes)
            bad += (not ok)
            f.write(f'{path}\t{nbytes}\t{b0}\t{h}\t{h0}\t{ok}\n')
    p('1. hashes: files', len(STATED), 'mismatches', bad)

    # 2. gene arm, independent path
    cat = pd.read_parquet(CAT)
    anc = cat[cat.primary_evidence_class.isin(['genetically_anchored', 'concordant'])]
    p('2. catalog rows', len(cat), 'anchored+concordant', len(anc), anc.primary_evidence_class.value_counts().to_dict(),
      'unique ensembl', anc.ensembl_id.nunique())
    symcol = 'symbol' if 'symbol' in anc.columns else [c for c in anc.columns if 'symbol' in c.lower()][0]
    co = pd.read_csv(COLOC, usecols=['gwas_name', 'gene', 'ensembl', 'PP.H4.abf', 'PP.H4.susie', 'method', 'top_snp'], low_memory=False)
    enz = co[co.gwas_name.isin(['UKBB_ALT', 'UKBB_AST', 'UKBB_GGT'])].copy()
    e = enz[enz.ensembl.isin(set(anc.ensembl_id))].copy()
    p('   enzyme rows joined on ensembl', len(e), 'genes', e.ensembl.nunique(), 'dup gene x trait', int(e.duplicated(['ensembl', 'gwas_name']).sum()))
    es = enz[enz.gene.isin(set(anc[symcol]))]
    p('   (symbol join for comparison with gate: genes', es.gene.nunique(), ')')
    e['pp_s'] = pd.to_numeric(e['PP.H4.susie'], errors='coerce'); e['pp_a'] = pd.to_numeric(e['PP.H4.abf'], errors='coerce')
    e['pp'] = np.where(np.isfinite(e.pp_s), e.pp_s, e.pp_a)
    e['src'] = np.where(np.isfinite(e.pp_s), 'susie', 'abf')
    p('   method x pp source crosstab', pd.crosstab(e.method, e.src).to_dict())
    e = e[np.isfinite(e.pp)]
    e['neg'] = -e.pp
    best = e.sort_values(['ensembl', 'neg', 'gwas_name']).groupby('ensembl', as_index=False).first()
    best['trait'] = best.gwas_name.str.replace('UKBB_', '')
    best['c'] = best.top_snp.str.split(':').str[0]; best['pos'] = best.top_snp.str.split(':').str[1].astype(int)
    p('   genes with best row', len(best), '; best pp exactly 0.5:', int((best.pp == 0.5).sum()),
      '; best pp in [0.49,0.51]:', int(((best.pp >= 0.49) & (best.pp <= 0.51)).sum()))
    for lab, m in [('>=0.8', best.pp >= 0.8), ('[0.5,0.8)', (best.pp >= 0.5) & (best.pp < 0.8)), ('<0.5', best.pp < 0.5)]:
        p(f'   tier {lab}: {int(m.sum())} {best[m].trait.value_counts().to_dict()}')
    p('   best top_snp chromosomes', sorted(best.c.unique(), key=lambda x: (len(x), x)))
    wb = pd.read_csv(W + '/best_enzyme_row.tsv', sep='\t')
    mm = best.merge(wb, on='ensembl', suffixes=('', '_w'))
    p('   agreement with write best_enzyme_row.tsv: rows', len(mm), 'same trait', int((mm.trait == mm.trait_w).sum()), 'same top_snp', int((mm.top_snp == mm.top_snp_w).sum()))
    app = best[best.pp >= 0.5]
    u = app.drop_duplicates(['trait', 'top_snp'])
    u8 = best[best.pp >= 0.8].drop_duplicates(['trait', 'top_snp'])
    def clusters(df, w=500000):
        n = 0
        for (t, ch), d in df.groupby(['trait', 'c']):
            ps = sorted(d.pos); n += 1
            for x, y in zip(ps, ps[1:]):
                if y - x > w: n += 1
        return n
    p('   applicable genes', len(app), 'unique trait x SNP', len(u), u.trait.value_counts().to_dict(), 'clusters', clusters(u))
    p('   PP>=0.8 genes', int((best.pp >= 0.8).sum()), 'unique trait x SNP', len(u8), 'clusters', clusters(u8))

    # COLOC source alleles at needed positions
    need = {t: set(zip(best[best.trait == t].c, best[best.trait == t].pos)) for t in FIELD}
    src = {}
    for t in FIELD:
        with open(f'{R}/GWAS/finemapping/data/sumstats/UKBB_{t}_reformatted_hg19.tsv') as f:
            next(f)
            for l in f:
                a = l.split('\t', 4); k = (a[0], int(a[1]))
                if k in need[t]: src.setdefault((t,) + k, set()).add(frozenset((a[2], a[3])))
    leads = {}
    for t in FIELD:
        leads[t] = []
        for l in open(f'{L}/ukbb_{t.lower()}_clumped500kb.tsv'):
            a = l.rstrip('\n').split('\t'); leads[t].append((a[0], int(a[1]), frozenset((a[2], a[3]))))
        p(f'   leads {t}: {len(leads[t])}; chromosomes {sorted({x[0] for x in leads[t]}, key=lambda x:(len(x),x))[:3]}..; non-SNV {sum(1 for x in leads[t] if any(len(q)!=1 for q in x[2]))}')

    jobs = [(t, need[t], {(c, q) for c, q, _ in leads[t]}) for t in FIELD]
    with Pool(4) as pool:
        ar = pool.apply_async(pc1)
        res = {r['t']: r for r in pool.map(scan, jobs)}
        pc = ar.get()
    p('3. PC1 urate:', pc)

    rows = []
    for t in FIELD:
        r = res[t]
        p(f'--- {t}: rows {r["rows"]}; n F {r["nF"]} M {r["nM"]} both {r["nB"]}; chrX rows {r["x_rows"]} PAR {r["x_par"]}')
        p(f'    chrX diag n {r["x_pairs"]} male/female MAF {r["x_maf_ratio"]} se_M/se_F chrX {r["x_se_ratio"]} autosomes {r["a_se_ratio"]} ratio {r["x_se_ratio"]/r["a_se_ratio"]:.4f}')
        p(f'    implied irnt SD F {r["sdF"]} M {r["sdM"]}')
        p(f'    both_sexes chrX GWS SNVs (MAF>=0.01 both sexes, low_conf false both_sexes) {r["xB_gws"]}; 500 kb clumped {r["xB_leads"]} [proxy only; the run uses the female/male IVW]')
        st_ = dict(no_src=0, src_multi=0, no_match=0, multi=0, lowconf_maf=0, ok=0)
        for g in best[best.trait == t].itertuples():
            k = (g.c, g.pos); s_ = src.get((t,) + k)
            tier = '>=0.8' if g.pp >= 0.8 else ('[0.5,0.8)' if g.pp >= 0.5 else '<0.5')
            if not s_: st_['no_src'] += 1; state = 'no_src'
            elif len(s_) > 1: st_['src_multi'] += 1; state = 'src_multi'
            else:
                al = next(iter(s_)); m = [x for x in r['site'].get(k, []) if x[0] == al]
                if not m: st_['no_match'] += 1; state = 'no_match'
                elif len(m) > 1: st_['multi'] += 1; state = 'multi'
                elif m[0][1] or m[0][2] or m[0][3] < 0.01 or m[0][4] < 0.01: st_['lowconf_maf'] += 1; state = 'lowconf_maf'
                else:
                    st_['ok'] += 1; state = 'ok'
                    x = m[0]; rows.append(dict(kind='gene', trait=t, id=g.ensembl, snp=g.top_snp, pp=g.pp, sF=x[5], sM=x[6], bB=x[7], sB=x[8]))
            if state != 'ok':
                p(f'    gene {g.ensembl} {g.gene} tier {tier} top_snp {g.top_snp} state {state} src {s_} neale {[tuple(sorted(x[0])) for x in r["site"].get(k, [])]}')
        p(f'    gene arm {t}: {st_}')
        ls = dict(not_found=0, multi=0, lowconf_maf=0, ok=0)
        for c, q, al in leads[t]:
            m = [x for x in r['site'].get((c, q), []) if x[0] == al]
            if not m: ls['not_found'] += 1
            elif len(m) > 1: ls['multi'] += 1
            elif m[0][1] or m[0][2] or m[0][3] < 0.01 or m[0][4] < 0.01: ls['lowconf_maf'] += 1
            else:
                ls['ok'] += 1; x = m[0]
                rows.append(dict(kind='lead', trait=t, id=f'{c}:{q}', snp=f'{c}:{q}', pp=None, sF=x[5], sM=x[6], bB=x[7], sB=x[8]))
        p(f'    lead arm {t}: {ls}')
        # NC1 feasibility: greedy thinning, 1 Mb from leads of this trait and from each other
        rng = np.random.default_rng(SEED)
        pl = r['pool']; order = rng.permutation(len(pl))
        lead_by_chr = {}
        for c, q, _ in leads[t]: lead_by_chr.setdefault(c, []).append(q)
        placed = {}
        n_pl = 0
        for j in order:
            c, q = pl[j]
            if any(abs(q - x) < 1000000 for x in lead_by_chr.get(c, [])): continue
            if any(abs(q - x) < 1000000 for x in placed.get(c, [])): continue
            placed.setdefault(c, []).append(q); n_pl += 1
        p(f'    NC1 feasibility: eligible SNVs {r["elig"]}; pool (every 20th, both_sexes p>0.5 proxy) {len(pl)}; greedy 1 Mb placements {n_pl}')

    d = pd.DataFrame(rows)
    g = d[d.kind == 'gene']
    gu = g.drop_duplicates(['trait', 'snp'])
    p('4. testable genes', len(g), 'by tier >=0.8', int((g.pp >= 0.8).sum()), '; unique trait x SNP', len(gu), gu.trait.value_counts().to_dict())
    g8 = g[g.pp >= 0.8].drop_duplicates(['trait', 'snp'])
    p('   PP>=0.8 testable unique trait x SNP', len(g8))
    lead = d[d.kind == 'lead']
    p('   testable leads', len(lead), lead.trait.value_counts().to_dict())

    # 5. power, as in the prespec, then under the audited equivalence rule
    def power(x, label, nsim=2000):
        x = x.copy(); m = len(x)
        x['sD'] = np.sqrt(x.sF**2 + x.sM**2)
        zc = norm.isf(0.05 / m / 2); mde = (zc + norm.isf(0.2)) * x.sD.median()
        nbig = int((x.bB.abs() >= (zc + norm.isf(0.2)) * x.sD).sum())
        # prespec normal approximation (log ratio)
        def pe(r):
            b = abs(r.bB); se = math.sqrt((r.sF / b)**2 + (r.sM / b)**2)
            if b / r.sM < 1.959964: return 0.0
            h = math.log(1.5) - 1.959964 * se
            return max(0.0, 2 * norm.cdf(h / se) - 1) if h > 0 else 0.0
        x['pe'] = x.apply(pe, axis=1)
        p(f'5. {label}: n {m}; median SE(diff) {x.sD.median():.5f}; BH single-hit z {zc:.3f}; MDE80 {mde:.4f}; median |b_both| {x.bB.abs().median():.4f}; n |b_both|>=MDE80 {nbig}; prespec expected bounded {x.pe.sum():.1f}; n P>=0.8 {(x.pe>=0.8).sum()}')
        # simulation under equal effects (bF = bM = bB), oriented positive; one scale only
        rng = np.random.default_rng(SEED)
        b = x.bB.abs().values; sF = x.sF.values; sM = x.sM.values
        z = 1.959964; lo, hi = 2 / 3, 1.5
        n_fieller = []; n_bh = []; n_both = []
        for _ in range(nsim):
            bF = b + rng.standard_normal(m) * sF; bM = b + rng.standard_normal(m) * sM
            # Fieller 95% inside [lo, hi]: bounded and both endpoints inside <=> the two one-sided 2.5% tests reject
            za = (bF - lo * bM) / np.sqrt(sF**2 + lo**2 * sM**2)
            zb = (hi * bM - bF) / np.sqrt(sF**2 + hi**2 * sM**2)
            bounded = np.abs(bM / sM) > z
            fiel = bounded & (za > z) & (zb > z)
            peq = np.maximum(norm.sf(za), norm.sf(zb))
            # BH
            o = np.argsort(peq); ranked = peq[o] * m / np.arange(1, m + 1)
            q = np.minimum.accumulate(ranked[::-1])[::-1]; qq = np.empty(m); qq[o] = np.minimum(q, 1)
            bh = qq < 0.05
            n_fieller.append(fiel.sum()); n_bh.append(bh.sum()); n_both.append((fiel & bh).sum())
        p(f'   simulated (equal effects, one scale, {nsim} draws, seed {SEED}): Fieller-95% bounded calls mean {np.mean(n_fieller):.1f}; '
          f'BH(TOST) q<0.05 calls mean {np.mean(n_bh):.1f}; both rules mean {np.mean(n_both):.1f} (5th-95th pct {np.percentile(n_both,5):.0f}-{np.percentile(n_both,95):.0f})')
        return x
    p('   PP>=0.8 testable genes', int((g.pp >= 0.8).sum()))
    power(gu, 'gene arm (unique trait x SNP)')
    power(lead, 'lead arm')

    # 6. NC1 lambda sampling distribution with 1,000 independent null SNPs
    rng = np.random.default_rng(SEED)
    lam = np.median(rng.chisquare(1, size=(20000, 1000)), axis=1) / chi2.ppf(0.5, 1)
    p(f'6. lambda (n=1000 independent nulls, 20,000 draws): median-chi2 constant {chi2.ppf(0.5,1):.4f}; sd {lam.std():.4f}; '
      f'2.5-97.5 pct {np.percentile(lam,2.5):.3f}-{np.percentile(lam,97.5):.3f}; P(lambda>1.10) {np.mean(lam>1.10):.4f}; '
      f'P(lambda>1.15) {np.mean(lam>1.15):.4f}; P(lambda<0.90) {np.mean(lam<0.90):.4f}; P(outside 0.90-1.10) {np.mean((lam<0.9)|(lam>1.1)):.4f}; '
      f'P(outside 0.85-1.15) {np.mean((lam<0.85)|(lam>1.15)):.4f}')
    pr = 1 - (1 - np.mean(lam > 1.10))**6; pr2 = 1 - (1 - np.mean(lam > 1.15))**6
    p(f'   P(at least one of 6 trait x scale GC-corrected by chance) at 1.10: {pr:.3f}; at 1.15: {pr2:.3f}')
    tail = rng.binomial(1000, 0.05, 200000) / 1000
    p(f'   tail fraction |z|>1.96, n=1000: 2.5-97.5 pct {np.percentile(tail,2.5):.4f}-{np.percentile(tail,97.5):.4f}; P(outside 0.035-0.065) {np.mean((tail<0.035)|(tail>0.065)):.4f}')
    d.to_csv(A + '/audit_power_rows.tsv', sep='\t', index=False)
    out.close()

if __name__ == '__main__':
    main()
