"""F2abce: Fig 2E (R02-06). For each non-European-only multi-signal COLOC gene,
list European SuSiE-positive (and ABF-positive) genes whose lead SNP lies within
1 Mb of the gene's best non-EUR lead SNP. Reads the COLOC r2 release directly."""
import pandas as pd
B = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT = B + 'figures/candidates/review-r3-fixes-20261015T112148Z/F2abce/'
m = pd.read_csv(B + 'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv',
                usecols=['gwas_name', 'gene', 'PP.H4.susie', 'PP.H4.abf', 'top_snp'])
tier = pd.read_csv(B + 'GWAS/finemapping/config/gwas_trait_tier.tsv', sep='\t')
main = set(tier.loc[tier.placement == 'main', 'study_name'])
reg = pd.read_csv(B + 'GWAS/finemapping/config/gwas_registry.tsv', sep='\t')
anc = dict(zip(reg.study_name, reg.ancestry))
m = m[m.gwas_name.isin(main) & m.gene.notna() & (m.gene.astype(str).str.strip() != '') & m.top_snp.notna()].copy()
m['anc'] = m.gwas_name.map(anc)
assert m.anc.notna().all()
m['c'] = m.top_snp.str.split(':').str[0]
m['p'] = m.top_snp.str.split(':').str[1].astype(int)
sig = m[m['PP.H4.susie'] > 0.5]
bins = sig.groupby('gene').anc.apply(set).apply(
    lambda a: 'EUR only' if a == {'EUR'} else ('shared' if 'EUR' in a else 'nonEUR'))
print(bins.value_counts(), len(bins))
mem = pd.read_csv(B + 'figures/main/fig2_genetics/panels/Fig2E_ancestry_unique_coloc_GWS_membership.tsv', sep='\t')
gws = set(mem.loc[mem.non_eur_gws_gate_pass.astype(str) == 'True', 'gene'])
print('GWS genes in installed membership:', len(gws), sorted(gws))
eur = sig[sig.anc == 'EUR']
abf = m[(m.anc == 'EUR') & (m['PP.H4.abf'] > 0.5)]
rows = []
for g in bins[bins == 'nonEUR'].index:
    r = sig[sig.gene == g].sort_values('PP.H4.susie', ascending=False).iloc[0]
    e = eur[(eur.c == r.c) & ((eur.p - r.p).abs() <= 1_000_000)]
    a = abf[(abf.c == r.c) & ((abf.p - r.p).abs() <= 1_000_000)]
    rows.append(dict(gene=g, gws_gate_pass=g in gws, best_gwas=r.gwas_name, best_anc=r.anc, lead=r.top_snp,
                     pp4_susie=round(r['PP.H4.susie'], 3),
                     eur_susie_genes_1Mb=','.join(sorted(e.gene.unique())), n_eur_susie_1Mb=e.gene.nunique(),
                     eur_same_lead=','.join(sorted(e[e.top_snp == r.top_snp].gene.unique())),
                     eur_abf_genes_1Mb=','.join(sorted(a.gene.unique())), n_eur_abf_1Mb=a.gene.nunique()))
t = pd.DataFrame(rows).sort_values(['gws_gate_pass', 'gene'], ascending=[False, True])
t.to_csv(OUT + 'Fig2E_nonEUR_only_locus_check.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
print(t.to_string(index=False))
print('nonEUR-only genes:', len(t))
print('with EUR SuSiE-positive gene within 1 Mb:', (t.n_eur_susie_1Mb > 0).sum())
print('with EUR SuSiE- or ABF-positive gene within 1 Mb:', ((t.n_eur_susie_1Mb > 0) | (t.n_eur_abf_1Mb > 0)).sum())
g9 = t[t.gws_gate_pass]
print('GWS genes:', len(g9), ' without EUR SuSiE-positive neighbour:', list(g9[g9.n_eur_susie_1Mb == 0].gene),
      ' without EUR SuSiE or ABF neighbour:', list(g9[(g9.n_eur_susie_1Mb == 0) & (g9.n_eur_abf_1Mb == 0)].gene))
# 1-Mb single-linkage clusters of the 9 GWS leads
g9 = g9.assign(c=g9.lead.str.split(':').str[0], p=g9.lead.str.split(':').str[1].astype(int)).sort_values(['c', 'p'])
nloc, prev = 0, None
for _, r in g9.iterrows():
    if prev is None or r.c != prev[0] or r.p - prev[1] > 1_000_000:
        nloc += 1
    prev = (r.c, r.p)
print('GWS genes collapse to', nloc, 'loci at 1 Mb')
