"""NG14 (L06-05) feasibility gate. Counts n, builds the GTEx-defined male-likeness
score, runs the positive control (sex main effect) in GSE130991, computes a
design-based MDE for the sex x age>=51 interaction, and estimates a plausible
interaction effect in open GTEx v8 liver. The GSE130991 interaction itself is NOT
computed (kept blind for the prespec)."""
import gzip, json, sys
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats

np.random.seed(20260923)
R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
OUT = sys.argv[1]
SM = R + '/Analysis/Hypotheses/N6-20261009T135702Z/downloads_untrusted/GSE130991_series_matrix.txt.gz'
GPL = R + '/Analysis/Hypotheses/N6-20261009T135702Z/downloads_untrusted/GPL20265_platform_table.txt'
SB = R + '/data/gtex_sb_eqtl/GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt'
GT = R + '/data/external/allelic_refs/gtex_v8/gene_tpm_2017-06-05_v8_liver.gct.gz'
GP = R + '/data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt'
DEG = R + '/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv'
res = {}

# --- GSE130991 metadata and expression
chars = []
with gzip.open(SM, 'rt') as h:
    for l in h:
        if l.startswith('!Sample_characteristics_ch1'):
            chars.append([x.strip('"\n') for x in l.split('\t')[1:]])
        if l.startswith('!series_matrix_table_begin'):
            break
    hdr = next(h).rstrip('\n').split('\t'); hdr=[x.strip('"') for x in hdr]
recs = [{} for _ in chars[0]]
for row in chars:
    for j, x in enumerate(row):
        if ':' in x:
            k, v = x.split(':', 1); recs[j][k.strip()] = v.strip()
meta = pd.DataFrame(recs, index=hdr[1:])
meta['female'] = (meta['gender'] == 'Female').astype(int)
meta['age'] = meta['age'].astype(float)
meta['old'] = (meta['age'] >= 51).astype(int)
meta['bmi'] = meta['bmi'].astype(float)
meta['hba1c'] = pd.to_numeric(meta['hba1c (%)'], errors='coerce')
meta['statin'] = (meta['statin treatment'] == 'Statin').astype(int)
meta['dmtx'] = (meta['diabetes treatment'] == 'Yes').astype(int)
res['n_by_sex_age51'] = meta.groupby(['gender', 'old']).size().to_dict()
res['n_by_sex_age51'] = {f'{k[0]}_{"ge51" if k[1] else "lt51"}': int(v) for k, v in res['n_by_sex_age51'].items()}
res['n_complete_covariates'] = int(meta[['bmi', 'hba1c', 'statin', 'dmtx']].notna().all(axis=1).sum())
res['statin_by_sex_age'] = {f'{a}_{b}': int(v) for (a, b), v in meta.groupby(['gender', 'old'])['statin'].sum().items()}
res['dmtx_by_sex_age'] = {f'{a}_{b}': int(v) for (a, b), v in meta.groupby(['gender', 'old'])['dmtx'].sum().items()}

# probe -> symbol, chromosome (main-design transcript clusters with a single gene)
gpl = pd.read_csv(GPL, sep='\t', skiprows=1, usecols=['ID', 'seqname', 'gene_assignment', 'category'], dtype=str)
gpl = gpl[gpl['category'].str.startswith('main', na=False)]
def sym(s):
    p = str(s).split(' // ')
    return p[1].strip() if len(p) > 1 and p[1].strip() not in ('---', '') else None
gpl['sym'] = gpl['gene_assignment'].map(sym)
gpl = gpl.dropna(subset=['sym'])
gpl['ID'] = gpl['ID'].astype(str)
res['hta_main_genes'] = int(gpl['sym'].nunique())
_sb = pd.read_csv(SB, sep='\t'); _sb = _sb[_sb['tissue'] == 'Liver']
_g = pd.read_csv(GT, sep='\t', skiprows=2, usecols=['Name', 'Description'])
_e2s = dict(zip(_g['Name'].str.split('.').str[0], _g['Description']))
_need = set(_sb['gene'].str.split('.').str[0].map(_e2s).dropna())
_d = pd.read_csv(DEG, usecols=['symbol', 'padj', 'logFC'])
_need |= set(_d.loc[(_d['padj'] < 0.05) & (_d['logFC'].abs() > 0.5), 'symbol'].dropna())
keep = set(gpl.loc[gpl['sym'].isin(_need), 'ID'])
rows = {}
with gzip.open(SM, 'rt') as h:
    for l in h:
        if l.startswith('!series_matrix_table_begin'):
            next(h); break
    for l in h:
        if l.startswith('!'): break
        f = l.rstrip('\n').split('\t'); pid = f[0].strip('"')
        if pid in keep: rows[pid] = np.array(f[1:], dtype=float)
e = pd.DataFrame(rows, index=hdr[1:]).T
g = gpl.set_index('ID').loc[e.index]
e = e.assign(sym=g['sym'].values, chrom=g['seqname'].values, m=e.mean(axis=1))
e = e.sort_values('m', ascending=False).drop_duplicates('sym')
chrom = e.set_index('sym')['chrom']
E = e.set_index('sym').drop(columns=['chrom', 'm'])
res['hta_genes_loaded'] = int(E.shape[0])

# GTEx liver sex-biased genes (positive effsize = female-biased per this session's XIST check)
gt = pd.read_csv(GT, sep='\t', skiprows=2, usecols=[1, 2]) if False else None
sb = pd.read_csv(SB, sep='\t')
sb = sb[sb['tissue'] == 'Liver'].copy()
sb['ens'] = sb['gene'].str.split('.').str[0]
gtx = pd.read_csv(GT, sep='\t', skiprows=2, index_col=0)
gtx['ens'] = gtx['Name'].str.split('.').str[0]
ens2sym = dict(zip(gtx['ens'], gtx['Description']))
sb['sym'] = sb['ens'].map(ens2sym)
xist = sb[sb['sym'] == 'XIST']['effsize'].tolist()
res['gtex_liver_sbgenes'] = int(len(sb)); res['XIST_liver_effsize'] = xist
sb['chrom'] = sb['sym'].map(chrom)
auto = sb[sb['chrom'].notna() & ~sb['chrom'].isin(['chrX', 'chrY'])]
fem = sorted(set(auto.loc[auto['effsize'] > 0, 'sym'])); mal = sorted(set(auto.loc[auto['effsize'] < 0, 'sym']))
res['autosomal_on_hta_female_biased'] = len(fem); res['autosomal_on_hta_male_biased'] = len(mal)

def score(M, up, dn):
    Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    return Z.reindex(up).mean() - Z.reindex(dn).mean()

meta['male_like'] = score(E, mal, fem).values
X = sm.add_constant(meta[['female', 'old', 'bmi', 'hba1c', 'statin', 'dmtx']].assign(fxo=meta['female'] * meta['old']))
ok = X.notna().all(axis=1)
# positive control: sex main effect (model without interaction, adjusted)
pc = sm.OLS(meta.loc[ok, 'male_like'], X.loc[ok].drop(columns=['fxo'])).fit(cov_type='HC3')
sd_within = np.sqrt(meta.groupby('female')['male_like'].var().mean())
res['pc_female_coef'] = float(pc.params['female']); res['pc_female_p'] = float(pc.pvalues['female'])
res['pc_sex_gap_in_within_sex_SD'] = float(-pc.params['female'] / sd_within)
pc_resid_sd = float(np.sqrt(pc.scale))
res['pc_resid_sd'] = pc_resid_sd
# design-based SE of interaction per unit residual SD (no outcome used)
XtX_inv = np.linalg.inv(X.loc[ok].T.values @ X.loc[ok].values)
se_unit = float(np.sqrt(XtX_inv[list(X.columns).index('fxo'), list(X.columns).index('fxo')]))
z = stats.norm.ppf(1 - 0.025) + stats.norm.ppf(0.8)  # Holm first step alpha 0.05/2, one-sided
z2 = stats.norm.ppf(1 - 0.05 / 2) + stats.norm.ppf(0.8)
res['interaction_se_per_resid_SD'] = se_unit
res['MDE_interaction_resid_SD_onesided_alpha0.025'] = z * se_unit
res['MDE_interaction_resid_SD_twosided_alpha0.05'] = z2 * se_unit
res['MDE_interaction_as_fraction_of_sex_gap'] = float(z * se_unit * pc_resid_sd / (-pc.params['female']))
# within-women only
w = ok & (meta['female'] == 1)
Xw = sm.add_constant(meta.loc[w, ['old', 'bmi', 'hba1c', 'statin', 'dmtx']])
iw = np.linalg.inv(Xw.T.values @ Xw.values)
res['MDE_women_old_vs_young_resid_SD_onesided_0.025'] = float(z * np.sqrt(iw[1, 1]))
res['collinearity_old_statin_corr'] = float(np.corrcoef(meta['old'], meta['statin'])[0, 1])

# secondary outcome coverage: Resource disease score genes on HTA
d = pd.read_csv(DEG)
sig = d[(d['padj'] < 0.05) & (d['logFC'].abs() > 0.5)]
up = set(sig.loc[sig['logFC'] > 0, 'symbol'].dropna()); dn = set(sig.loc[sig['logFC'] < 0, 'symbol'].dropna())
res['disease_up'] = len(up); res['disease_dn'] = len(dn)
res['disease_up_on_hta'] = len(up & set(E.index)); res['disease_dn_on_hta'] = len(dn & set(E.index))
res['disease_up_on_hta_any'] = len(up & set(gpl['sym'])); res['disease_dn_on_hta_any'] = len(dn & set(gpl['sym']))

# --- plausible effect: GTEx v8 liver, women >=50 (age band) shift toward male, beyond men
ph = pd.read_csv(GP, sep='\t')
L = np.log2(gtx.drop(columns=['Name', 'Description', 'ens']).set_index(gtx['Description']).groupby(level=0).mean() + 1)
L = L[(L > 1).mean(axis=1) > 0.5]
subj = ['-'.join(c.split('-')[:2]) for c in L.columns]
pm = ph.set_index('SUBJID').loc[subj]
pm.index = L.columns
pm['female'] = (pm['SEX'] == 2).astype(int)
pm['old'] = pm['AGE'].isin(['50-59', '60-69', '70-79']).astype(int)
pm['male_like'] = score(L, mal, fem).values
res['gtex_n_by_sex_age50'] = {f'{"F" if a else "M"}_{"ge50" if b else "lt50"}': int(v) for (a, b), v in pm.groupby(['female', 'old']).size().items()}
Xg = sm.add_constant(pm[['female', 'old']].assign(fxo=pm['female'] * pm['old'], dthhrdy=pm['DTHHRDY'].fillna(-1)))
fg = sm.OLS(pm['male_like'], Xg).fit(cov_type='HC3')
gap = -fg.params['female']
res['gtex_sex_gap_young'] = float(gap); res['gtex_resid_sd'] = float(np.sqrt(fg.scale))
res['gtex_interaction'] = float(fg.params['fxo']); res['gtex_interaction_se'] = float(fg.bse['fxo'])
res['gtex_interaction_ci95'] = [float(x) for x in fg.conf_int().loc['fxo']]
res['gtex_interaction_resid_SD'] = float(fg.params['fxo'] / np.sqrt(fg.scale))
res['gtex_interaction_frac_gap'] = float(fg.params['fxo'] / gap)
res['gtex_n_genes_expressed_in_score'] = [len(set(mal) & set(L.index)), len(set(fem) & set(L.index))]
json.dump(res, open(OUT + '/gate_counts.json', 'w'), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
