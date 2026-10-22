# H74 audit: recount every n stated in PRESPEC_H74 and check structural facts.
# Does NOT compute any anchor-to-gene distance, nearest gene, E1, E2, T_A, T_B or
# the class of any nearest gene. Structure only (loci, tested genes, TSS mapping,
# ties among tested-gene TSS coordinates, window overlap, platelet locus count).
import re, sys, json
import numpy as np, pandas as pd

np.random.seed(20260923)  # no random draws are made; set for the round rule
R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
COLOC = R + 'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv'
GTF = R + 'data/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf'
TIER = R + 'GWAS/finemapping/config/gwas_trait_tier.tsv'

tier = pd.read_csv(TIER, sep='\t')
print('tier table rows', len(tier))
print(tier.groupby(['tier', 'tier_label']).size().to_string())
t12 = set(tier.study_name[tier.tier.isin([1, 2])])
plat = set(tier.study_name[tier.study_name.str.startswith('MVP_Platelet')])
print('tier1/2 study names in table', len(t12), 'platelet', sorted(plat))

cols = ['gwas_name', 'gene', 'ensembl', 'chr', 'PP.H3.susie', 'PP.H4.susie', 'PP.H4.abf',
        'method', 'top_snp', 'top_snp_PP']
d = pd.read_csv(COLOC, usecols=cols)
print('rows', len(d), 'distinct gwas_name', d.gwas_name.nunique())
print('method x PP.H4.susie notna:')
print(pd.crosstab(d.method, d['PP.H4.susie'].notna()).to_string())
print('top_snp NA by method:', d.groupby('method').top_snp.apply(lambda x: x.isna().sum()).to_dict())
print('ensembl has version suffix (contains "."):', d.ensembl.astype(str).str.contains(r'\.').sum())

s = d[d.gwas_name.isin(t12) & d['PP.H4.susie'].notna()].copy()
print('tier12 susie rows', len(s), 'studies', s.gwas_name.nunique(), 'genes', s.ensembl.nunique())
st = tier.set_index('study_name').tier
print('tier12 studies with susie rows by tier:', s.drop_duplicates('gwas_name').gwas_name.map(st).value_counts().to_dict())
print('tier1/2 names without susie rows:', sorted(t12 - set(s.gwas_name)))
print('top_snp chr values (tier12 susie):', sorted(s.top_snp.str.split(':').str[0].unique(), key=lambda x: (len(x), x)))
print('top_snp field counts:', s.top_snp.str.count(':').value_counts().to_dict())
s['pos'] = s.top_snp.str.split(':').str[1].astype(int)
s['chr'] = s.top_snp.str.split(':').str[0]

unt = d[d.gwas_name.isin(t12) & (d.method == 'susie_untestable_insufficient_shared_posterior')]
print('tier12 susie_untestable rows', len(unt), 'genes', unt.ensembl.nunique(),
      'genes of these with any tier12 finite PP.H4.susie row', len(set(unt.ensembl) & set(s.ensembl)))

# GTF gene records
tss = {}; n_gene = 0; n_par = 0; dup = []; pc = 0; idfmt = {}
for l in open(GTF):
    if l[0] == '#':
        continue
    f = l.split('\t')
    if f[2] != 'gene':
        continue
    n_gene += 1
    full = re.search('gene_id "([^"]+)"', f[8]).group(1)
    k = re.sub(r'^ENSG\d+', 'ENSG', full)
    k = re.sub(r'\d+', 'N', k)
    idfmt[k] = idfmt.get(k, 0) + 1
    if full.endswith('_PAR_Y'):
        n_par += 1
        continue
    gid = re.search('gene_id "([^".]+)', f[8]).group(1)
    gt = re.search('gene_type "([^"]+)', f[8]).group(1)
    gn = re.search('gene_name "([^"]+)', f[8]).group(1)
    t = int(f[3]) if f[6] == '+' else int(f[4])
    if gid in tss:
        dup.append((gid, tss[gid][0], f[0]))
    tss[gid] = (f[0].replace('chr', ''), t, gt, gn)
    pc += gt == 'protein_coding'
print('GTF gene records', n_gene, '_PAR_Y skipped', n_par, 'kept unique ids', len(tss),
      'protein_coding kept', sum(v[2] == 'protein_coding' for v in tss.values()))
print('gene_id formats', idfmt)
print('duplicate stripped gene ids after PAR skip', len(dup), dup[:10])
print('GTF chromosomes', sorted(set(v[0] for v in tss.values()), key=lambda x: (len(x), x)))

s['tss'] = s.ensembl.map(lambda g: tss.get(g, (None, None))[1])
s['gtype'] = s.ensembl.map(lambda g: tss.get(g, (None, None, None))[2])
print('tier12 susie rows with TSS', round(s.tss.notna().mean(), 6), 'pc frac', round((s.gtype == 'protein_coding').mean(), 6))
print('tier12 susie genes without TSS', s.loc[s.tss.isna(), 'ensembl'].nunique())


def loci_of(sub):
    pos = sub[sub['PP.H4.susie'] > 0.5]
    out = []
    for c, g in pos.groupby('chr'):
        ps = sorted(g.pos.unique()); cur = [ps[0]]
        for p in ps[1:]:
            if p - cur[-1] <= 1_000_000:
                cur.append(p)
            else:
                out.append((c, cur[0], cur[-1])); cur = [p]
        out.append((c, cur[0], cur[-1]))
    return out


locs = loci_of(s)
print('tier12 positive loci', len(locs))
rows = []; gene_loci = {}
for i, (c, a, b) in enumerate(locs):
    w = s[(s.chr == c) & (s.pos >= a - 1_000_000) & (s.pos <= b + 1_000_000)]
    gs = w.groupby('ensembl').agg(h4=('PP.H4.susie', 'max'), h3=('PP.H3.susie', 'max'), tss=('tss', 'first'))
    for g in gs.index:
        gene_loci.setdefault(g, []).append(i)
    C = gs[gs.h4 > 0.5]; N = gs[gs.h4 <= 0.5]
    m = gs[gs.tss.notna()]
    mC = m[m.h4 > 0.5]; mN = m[m.h4 <= 0.5]
    tie_any = m.tss.duplicated(keep=False).sum()
    npos = w[w['PP.H4.susie'] > 0.5].pos.nunique()
    rows.append(dict(locus=i, chr=c, start=a, end=b, n_tested=len(gs), n_h4=len(C), n_h3=int(((gs.h3 > 0.8) & (gs.h4 <= 0.5)).sum()),
                     n_mapped=len(m), n_mapped_C=len(mC), n_mapped_N=len(mN), n_genes_sharing_tss=int(tie_any),
                     n_pos_snps=npos))
L = pd.DataFrame(rows)
gate = pd.read_csv(R + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG13/locus_structure.tsv', sep='\t')
same = (L[['chr', 'start', 'end', 'n_tested', 'n_h4', 'n_h3']].astype(str).values == gate[['chr', 'start', 'end', 'n_tested', 'n_h4', 'n_h3']].astype(str).values).all()
print('per-locus table identical to gate locus_structure.tsv:', bool(same))
print('n_h4==1', (L.n_h4 == 1).sum(), 'n_h4>=2', (L.n_h4 >= 2).sum(), 'all tested coloc', (L.n_tested == L.n_h4).sum(),
      'n_h3==0', (L.n_h3 == 0).sum(), 'loci with >1 distinct positive top SNP', (L.n_pos_snps > 1).sum())
print('informative (gate rule, n_tested>n_h4)', (L.n_tested > L.n_h4).sum())
print('informative E3 as written (>=1 mapped non-C gene)', (L.n_mapped_N >= 1).sum())
print('informative E3 patched (>=1 mapped C and >=1 mapped non-C)', ((L.n_mapped_N >= 1) & (L.n_mapped_C >= 1)).sum())
print('loci whose C has no mapped TSS', (L.n_mapped_C == 0).sum())
inf = L[(L.n_mapped_N >= 1) & (L.n_mapped_C >= 1)]
p = inf.n_mapped_C / inf.n_mapped
print('patched informative: null mean P', round(p.mean(), 4), 'null mean T_A', round(p.sum(), 2))
print('loci with >=2 tested mapped genes sharing an identical TSS coordinate', (L.n_genes_sharing_tss > 0).sum(),
      'genes involved', int(L.n_genes_sharing_tss.sum()))
# window overlap between adjacent loci on a chromosome
ov = 0
for c, g in L.groupby('chr'):
    g = g.sort_values('start')
    e = (g.end + 1_000_000).values[:-1]; st_ = (g.start - 1_000_000).values[1:]
    ov += int((st_ <= e).sum())
print('adjacent locus pairs with overlapping +-1Mb windows', ov)
multi = {g: v for g, v in gene_loci.items() if len(v) > 1}
print('tested genes counted at >1 locus', len(multi))
L.to_csv(sys.argv[1] + '/audit_locus_structure.tsv', sep='\t', index=False)

# anchor-rule ties (no positions used beyond identity of top_snp)
pos = s[s['PP.H4.susie'] > 0.5].copy()
nt = 0
for i, (c, a, b) in enumerate(locs):
    w = pos[(pos.chr == c) & (pos.pos >= a) & (pos.pos <= b)]
    top = w[w['PP.H4.susie'] == w['PP.H4.susie'].max()]
    if top.top_snp.nunique() > 1:
        nt += 1
print('loci where max PP.H4.susie is tied across different top SNPs', nt)
print('positive tier12 rows with PP.H4.abf <= 0.5 (top_snp comes from ABF):', int((pos['PP.H4.abf'] <= 0.5).sum()), 'of', len(pos))

# platelet structure only
sp = d[d.gwas_name.isin(plat) & d['PP.H4.susie'].notna()].copy()
sp['pos'] = sp.top_snp.str.split(':').str[1].astype(int); sp['chr'] = sp.top_snp.str.split(':').str[0]
print('platelet susie rows', len(sp), 'studies', sp.gwas_name.nunique(), 'positive loci (same 1-Mb rule)', len(loci_of(sp)))
