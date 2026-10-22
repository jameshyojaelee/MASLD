# Blind re-check F3cdg (2026-10-15): re-derive R03-01, R03-02 union, R14-07, R16-02, R14-09
# numbers from the release files. Written independently of the producer scripts.
import sys, os, subprocess
import numpy as np, pandas as pd
from scipy.stats import spearmanr

B = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
OUT = sys.argv[1]
R = {}
def put(k, v):
    R[k] = v; print(k, v, flush=True)

# ---------- R03-01: Figure 3C membership by direction ----------
coh = pd.read_csv(f'{B}/figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv', sep='\t',
                  usecols=['gene_id_versioned', 'dataset', 'logFC', 'FDR', 'bh_family_size'])
put('R0301_bh_family_sizes', sorted(coh.bh_family_size.unique().tolist()))
put('R0301_rows_per_dataset', coh.groupby('dataset').size().to_dict())
calls = coh[(coh.FDR < 0.05) & (coh.logFC.abs() > 0.5)].copy()
calls['sgn'] = np.sign(calls.logFC).astype(int)
put('R0301_calls_per_dataset', calls.groupby('dataset').size().to_dict())
put('R0301_memberships', len(calls))
gg = calls.groupby('gene_id_versioned').agg(k=('dataset', 'nunique'), up=('sgn', lambda s: int((s > 0).sum())), dn=('sgn', lambda s: int((s < 0).sum())))
put('R0301_union', len(gg))
put('R0301_k_counts', gg.k.value_counts().sort_index().to_dict())
multi = gg[gg.k >= 2]
mixed = multi[(multi.up > 0) & (multi.dn > 0)]
same = multi[~multi.index.isin(mixed.index)]
put('R0301_multi', len(multi)); put('R0301_mixed', len(mixed)); put('R0301_mixed_pct', round(100 * len(mixed) / len(multi), 2))
put('R0301_mixed_by_k', mixed.k.value_counts().sort_index().to_dict())
put('R0301_same_by_k', same.k.value_counts().sort_index().to_dict())
put('R0301_same_total', len(same))
g126 = set(calls[calls.dataset == 'GSE126848'].gene_id_versioned)
put('R0301_mixed_with_GSE126848', int(mixed.index.isin(g126).sum()))
c4 = calls[calls.dataset != 'GSE126848']
g4 = c4.groupby('gene_id_versioned').agg(k=('dataset', 'nunique'), up=('sgn', lambda s: int((s > 0).sum())), dn=('sgn', lambda s: int((s < 0).sum())))
m4 = g4[g4.k >= 2]
put('R0301_noGSE126848_multi', len(m4)); put('R0301_noGSE126848_mixed', int(((m4.up > 0) & (m4.dn > 0)).sum()))
put('R0301_noGSE126848_mixed_pct', round(100 * ((m4.up > 0) & (m4.dn > 0)).sum() / len(m4), 2))
mx = gg[['up', 'dn']].max(axis=1)
put('R0301_larger_of_up_down_counts', mx.value_counts().sort_index().to_dict())
# per-dataset ribbon counts by tier
def tier(r):
    if r.k == 1: return 'unique'
    if r.up > 0 and r.dn > 0: return 'mixed'
    return f'same{r.k}'
gg['tier'] = gg.apply(tier, axis=1)
rib = calls.merge(gg[['tier']], left_on='gene_id_versioned', right_index=True).groupby(['dataset', 'tier']).size().unstack(fill_value=0)
rib.to_csv(f'{OUT}/r0301_ribbons_recheck.tsv', sep='\t')
put('R0301_ribbon_memberships_by_tier', rib.sum().to_dict())
# pooled direction agreement per cohort
deg = pd.read_csv(f'{B}/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv')
pc = deg[(deg.padj < 0.05) & (deg.logFC.abs() > 0.5)]
put('R0301_pooled_calls', len(pc))
agree = {}
for ds, d in coh.groupby('dataset'):
    m = pc[['gene', 'logFC']].merge(d[['gene_id_versioned', 'logFC']], left_on='gene', right_on='gene_id_versioned', suffixes=('_p', '_c'))
    agree[ds] = (len(m), round(100 * (np.sign(m.logFC_p) == np.sign(m.logFC_c)).mean(), 2))
put('R0301_pooled_direction_agreement_pct', agree)

# ---------- R03-02: Figure 3C union under GSE126848 sensitivities (uses my own R refit table) ----------
gt = pd.read_csv(f'{OUT}/r0302_gse126848_genes_recheck.tsv', sep='\t')
zero_called = set(gt[gt.base_call & gt.zero_all30].gene)
c_nz = calls[~((calls.dataset == 'GSE126848') & calls.gene_id_versioned.isin(zero_called))]
gz = c_nz.groupby('gene_id_versioned').dataset.nunique()
put('R0302_union_without_zero_read_calls', dict(union=len(gz), unique=int((gz == 1).sum()), memberships=len(c_nz)))

# ---------- R14-07: Figure 3G magnitude test ----------
g3 = pd.read_csv(f'{B}/figures/candidates/fig3g-working-layout-774-20261012T011006Z/source_tables/fig3g_raw_treat_interface.tsv', sep='\t')
f = g3[np.isfinite(g3.susie_pp4) & np.isfinite(g3.logFC)]
put('R1407_n', len(f)); put('R1407_states', f.state.value_counts().to_dict())
r_s = spearmanr(f.susie_pp4, f.logFC); r_a = spearmanr(f.susie_pp4, f.logFC.abs())
put('R1407_signed', (r_s.correlation, r_s.pvalue)); put('R1407_abs', (r_a.correlation, r_a.pvalue))
pp = f[f.susie_pp4 > 0.5]
put('R1407_pp4gt05', dict(n=len(pp), up=int((pp.logFC > 0).sum()), down=int((pp.logFC < 0).sum())))
untest = ['ADAM15', 'ASAP3', 'CKB', 'GSTM2', 'LILRB3', 'PSMC5', 'RITA1']
f7 = f[~f.symbol.isin(untest)]
r7 = spearmanr(f7.susie_pp4, f7.logFC.abs())
put('R1407_minus7_untestable', dict(n=len(f7), abs_rho=r7.correlation, abs_p=r7.pvalue, neither=int((f7.state == 'Neither').sum())))

# ---------- R16-02: Figure 3F FDR families ----------
st = pd.read_csv(f'{B}/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_all_gene_results.tsv', sep='\t')
adj = st[st.contrast.isin(['F0_to_F1', 'F1_to_F2', 'F2_to_F3', 'F3_to_F4'])].copy()
put('R1602_rows_per_contrast', adj.groupby('contrast').size().to_dict())
from statsmodels.stats.multitest import multipletests
adj['bh'] = adj.groupby('contrast')['P.Value'].transform(lambda p: multipletests(p, method='fdr_bh')[1])
put('R1602_max_abs_bh_minus_FDR', float((adj.bh - adj.FDR).abs().max()))
t12 = pd.read_csv(f'{B}/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv')
s436 = set(t12[t12.coloc_best_susie_pp4 > 0.5].ensembl)
put('R1602_susie_set', len(s436))
adj['base'] = adj.gene_id_versioned.str.split('.').str[0]
a = adj[adj.base.isin(s436)]
put('R1602_susie_genes_in_stage_table', a.base.nunique())
sel = a[a.FDR < 0.05]
put('R1602_selected_genes', sel.base.nunique())
put('R1602_marks_by_contrast', sel.groupby('contrast').size().to_dict())
sp = a[a.FDR < 0.0125]
put('R1602_split_genes', sp.base.nunique())
put('R1602_split_marks_by_contrast', sp.groupby('contrast').size().to_dict())
put('R1602_split_gene_names', sorted(sp.gene_name.unique().tolist()))
# BH at 0.0125 applied freshly to raw p (equivalent check)
adj['rej0125'] = adj.groupby('contrast')['P.Value'].transform(lambda p: multipletests(p, alpha=0.0125, method='fdr_bh')[0])
put('R1602_split_genes_fresh_BH', a.assign(r=adj.loc[a.index, 'rej0125']).query('r').base.nunique())

# ---------- R14-09: GNMT COLOC rows ----------
master = f'{B}/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv'
hdr = pd.read_csv(master, nrows=0).columns.tolist()
gn = subprocess.run(['grep', ',ENSG00000124713,', master], capture_output=True, text=True).stdout
from io import StringIO
gnmt = pd.read_csv(StringIO(gn), header=None, names=hdr)
gnmt.to_csv(f'{OUT}/r1409_gnmt_all_rows.tsv', sep='\t', index=False)
tier12_studies = None
put('R1409_gnmt_rows', len(gnmt))
put('R1409_gnmt_best_susie', gnmt.sort_values('PP.H4.susie', ascending=False)[['gwas_name', 'PP.H4.susie']].head(3).values.tolist())
put('R1409_gnmt_best_abf', gnmt.sort_values('PP.H4.abf', ascending=False)[['gwas_name', 'PP.H4.abf', 'method']].head(3).values.tolist())
r = gnmt[gnmt.gwas_name == '2021_34841290_NAFLD_EUR'].iloc[0]
put('R1409_nafld_row', dict(abf=r['PP.H4.abf'], susie=r['PP.H4.susie'], method=r['method'], lambda_s_locus=r['lambda_s_locus'], n_cs_pairs=r['n_cs_pairs'], elig=r['eligibility_source'], top_snp=r['top_snp'], top_snp_PP=r['top_snp_PP']))
put('R1409_gnmt_methods', gnmt['method'].value_counts().to_dict())
st_rows = subprocess.run(['grep', '^2021_34841290_NAFLD_EUR,', master], capture_output=True, text=True).stdout
nf = pd.read_csv(StringIO(st_rows), header=None, names=hdr, low_memory=False)
put('R1409_study_rows', len(nf))
put('R1409_study_methods', nf['method'].value_counts(dropna=False).to_dict())
put('R1409_study_susie_finite', int(np.isfinite(nf['PP.H4.susie']).sum()))
put('R1409_study_elig', nf['eligibility_source'].value_counts(dropna=False).to_dict())
put('R1409_study_lambda_nonNA_among_abf_fallback', int(nf.loc[nf['method'] == 'abf_fallback', 'lambda_s_locus'].notna().sum()))
gl = t12[t12.ensembl == 'ENSG00000124713']
put('R1409_tier12_gene_level', gl.to_dict('records'))

pd.Series({k: str(v) for k, v in R.items()}).to_csv(f'{OUT}/recheck_values.tsv', sep='\t', header=False)
