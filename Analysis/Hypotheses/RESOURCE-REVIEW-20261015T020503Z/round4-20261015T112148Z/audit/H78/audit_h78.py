"""H78 audit: recount every n, key and coverage stated in PRESPEC_H78.md.
No PC/NC control value under the H78 definitions is computed, no 24Q4 common-essential
fraction for any gene set, and no HPA field value for CEG2 or the PC-3 genes.
The gate-style CEG2/NEGv1 LOEUF numbers (34,650-row deciles, symbol join) are recomputed
only to check the provenance of values the prespec already states."""
import re, gzip, hashlib, os, sys
import numpy as np, pandas as pd
np.random.seed(20260923)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500); pd.set_option('display.max_columns', 50)
R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
RR = R + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/'
W = RR + 'write/H78/'
D24 = R + 'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf/source/depmap_24q4/'
CAT = R + 'figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
GN = R + 'data/external/gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv'
HPA = R + 'data/external/hpa/proteinatlas.tsv'
GTF = '/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz'
HART = R + 'data/public_r4/Hart2017_CEG2/'

def hashes(p):
    s, m = hashlib.sha256(), hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            s.update(b); m.update(b)
    return os.path.getsize(p), s.hexdigest(), m.hexdigest()

print('=== A. bytes, sha256, md5')
files = [CAT, GN, HPA, D24 + 'CRISPRGeneEffect.csv', D24 + 'Model.csv', W + 'depmap24q4/CRISPRInferredCommonEssentials.csv',
         W + 'depmap24q4/AchillesCommonEssentialControls.csv', HART + 'CEGv2.txt', HART + 'NEGv1.txt', GTF,
         RR + 'gate/NG17/count_ng17.py', RR + 'gate/NG17/count_22251965.log', W + 'count_h78.py', W + 'count_22252308.log',
         W + 'count_depmap24q4.py', W + 'count24q4_22252330.log', R + 'Analysis/downstream_analysis/essentiality/CRISPRGeneEffect.csv']
rows = []
for p in files:
    b, s, m = hashes(p); rows.append((p.replace(R, ''), b, s, m)); print(p.replace(R, ''), b, s, m)
pd.DataFrame(rows, columns=['path', 'bytes', 'sha256', 'md5']).to_csv(sys.argv[1] + '/sha256_recomputed.tsv', sep='\t', index=False)

print('\n=== B. Catalog')
g = pd.read_parquet(CAT)
print('columns', list(g.columns))
g['ens'] = g.ensembl_id.astype(str).str.split('.').str[0]
print('rows', len(g), 'unique ensembl_id', g.ensembl_id.nunique(), 'unique ens (no version)', g.ens.nunique(),
      'non-ENSG', (~g.ens.str.startswith('ENSG')).sum(), 'with version dot', g.ensembl_id.astype(str).str.contains(r'\.').sum(),
      'PAR_Y', g.ensembl_id.astype(str).str.contains('PAR').sum())
print(pd.crosstab(g.next_experiment_rule_id, g.primary_evidence_class, margins=True).to_string())
print('symbol missing', g.symbol.isna().sum(), 'duplicated symbols', g.symbol.dropna().duplicated().sum())

print('\n=== C. gnomAD v4.1')
gn = pd.read_csv(GN, sep='\t', usecols=['gene', 'gene_id', 'transcript', 'mane_select', 'canonical', 'lof.oe_ci.upper', 'lof.pLI',
                                        'lof.oe_ci.upper_bin_decile', 'constraint_flags', 'syn.oe'], low_memory=False)
print('rows', len(gn), 'mane_select values', gn.mane_select.astype(str).value_counts().to_dict())
m = gn[gn.mane_select.astype(str).str.lower().eq('true')].copy()
ensm = m[m.gene_id.astype(str).str.startswith('ENSG')].copy(); ezm = m[m.gene_id.astype(str).str.fullmatch(r'\d+')].copy()
print('MANE rows', len(m), 'ENSG rows', len(ensm), 'unique ENSG', ensm.gene_id.nunique(), 'Entrez rows', len(ezm), 'unique', ezm.gene_id.nunique(),
      'other', len(m) - len(ensm) - len(ezm))
print('ENSG MANE transcripts startswith ENST', ensm.transcript.astype(str).str.startswith('ENST').sum(),
      '| Entrez MANE transcripts startswith NM_', ezm.transcript.astype(str).str.startswith('NM_').sum())
print('ENSG MANE with LOEUF', ensm['lof.oe_ci.upper'].notna().sum(), 'pLI', ensm['lof.pLI'].notna().sum(),
      'official decile', ensm['lof.oe_ci.upper_bin_decile'].notna().sum(), 'syn.oe', ensm['syn.oe'].notna().sum(),
      'LOEUF present & pLI missing', (ensm['lof.oe_ci.upper'].notna() & ensm['lof.pLI'].isna()).sum())
cf = ensm.constraint_flags
print('constraint_flags NaN', cf.isna().sum(), '| raw value counts (top 8):'); print(cf.astype(str).value_counts().head(8).to_string())
nonempty = cf.notna() & ~cf.astype(str).str.strip().isin(['[]', '', 'nan', 'None'])
print('constraint_flags truly non-empty (ENSG MANE rows)', int(nonempty.sum()), '| among rows with LOEUF', int((nonempty & ensm['lof.oe_ci.upper'].notna()).sum()))
print('ENSG symbols duplicated', ensm.gene.duplicated().sum(), '| Entrez symbols duplicated', ezm.gene.duplicated().sum())
L = ensm.dropna(subset=['lof.oe_ci.upper'])
s2e = ezm.drop_duplicates('gene').set_index('gene')
pair_all = ensm.gene.isin(s2e.index); pair_L = L.gene.isin(s2e.index)
print('ENSG MANE rows with Entrez partner by symbol: all', int(pair_all.sum()), 'of', len(ensm), '| with LOEUF', int(pair_L.sum()), 'of', len(L))
pp = L[pair_L].merge(ezm[['gene', 'gene_id', 'lof.oe_ci.upper', 'lof.pLI']].rename(columns={'gene_id': 'entrez', 'lof.oe_ci.upper': 'loeuf_ez', 'lof.pLI': 'pli_ez'}), on='gene')
d = (pp['lof.oe_ci.upper'] - pp.loeuf_ez).abs()
print('paired rows LOEUF |ENSG-Entrez| : exact equal', int((d == 0).sum()), '<=0.05', int((d <= 0.05).sum()), '>0.2', int((d > 0.2).sum()),
      'Entrez LOEUF NaN', int(pp.loeuf_ez.isna().sum()), 'median', float(d.median()))
pair = ensm[pair_all][['gene_id', 'gene']].merge(ezm[['gene', 'gene_id']].rename(columns={'gene_id': 'entrez'}), on='gene')
pair['entrez'] = pair.entrez.astype(int)
print('pair table rows', len(pair), 'ENSG dup', pair.gene_id.duplicated().sum(), 'Entrez dup', pair.entrez.duplicated().sum())
# recomputed decile (prespec rule) only to check bin sizes; no gene-set fraction
Ls = L.sort_values('gene_id').copy()
Ls['dec'] = pd.qcut(Ls['lof.oe_ci.upper'].rank(method='first'), 10, labels=False)
print('recomputed decile bin sizes', Ls.dec.value_counts().sort_index().tolist(), 'decile0 max LOEUF', float(Ls.loc[Ls.dec == 0, 'lof.oe_ci.upper'].max()),
      '| official decile value counts', ensm['lof.oe_ci.upper_bin_decile'].value_counts().sort_index().to_dict())
print('ties at decile0/1 boundary: LOEUF values shared by dec0 and dec1 rows',
      len(set(Ls.loc[Ls.dec == 0, 'lof.oe_ci.upper']) & set(Ls.loc[Ls.dec == 1, 'lof.oe_ci.upper'])))

print('\n=== D. DepMap 24Q4')
dmh = pd.read_csv(D24 + 'CRISPRGeneEffect.csv', usecols=[0])
cols = pd.read_csv(D24 + 'CRISPRGeneEffect.csv', nrows=0).columns
print('first column name', cols[0], '| lines', len(dmh), 'unique lines', dmh.iloc[:, 0].nunique(), '| genes', len(cols) - 1)
dm = pd.DataFrame({'col': cols[1:]}); dm['sym'] = dm.col.str.split(' \\(').str[0]; dm['entrez'] = dm.col.str.extract(r'\((\d+)\)')[0].astype(int)
print('unique entrez', dm.entrez.nunique(), 'unique sym', dm.sym.nunique())
ce = pd.read_csv(W + 'depmap24q4/CRISPRInferredCommonEssentials.csv'); ce['entrez'] = ce.iloc[:, 0].str.extract(r'\((\d+)\)')[0].astype(int)
ctl = pd.read_csv(W + 'depmap24q4/AchillesCommonEssentialControls.csv'); ctl['entrez'] = ctl.iloc[:, 0].str.extract(r'\((\d+)\)')[0].astype(int)
print('inferred CE', len(ce), 'unique', ce.entrez.nunique(), 'in matrix', ce.entrez.isin(dm.entrez).sum(), '| controls', len(ctl), 'in inferred', ctl.entrez.isin(ce.entrez).sum())
mdl = pd.read_csv(D24 + 'Model.csv', usecols=lambda c: c in ('ModelID', 'OncotreeLineage'))
print('Model.csv rows', len(mdl), '| gene-effect lines in Model', dmh.iloc[:, 0].isin(mdl.ModelID).sum(),
      '| Liver lineage lines in gene-effect', int(dmh.iloc[:, 0].isin(mdl.loc[mdl.OncotreeLineage.astype(str).str.contains('Liver', case=False), 'ModelID']).sum()) if 'OncotreeLineage' in mdl else 'NA')
g2 = g.merge(pair.rename(columns={'gene_id': 'ens'})[['ens', 'entrez']], on='ens', how='left')
g2L = g.merge(pp.rename(columns={'gene_id': 'ens'})[['ens', 'entrez']].assign(entrez=lambda x: x.entrez.astype(int)), on='ens', how='left')
print('Catalog rows after Entrez merge', len(g2), '| Catalog genes with Entrez (all-MANE pairing)', g2.entrez.notna().sum(), '| (LOEUF-only pairing)', g2L.entrez.notna().sum())
print('Catalog Entrez mapping to >1 Catalog ENSG', int(g2.dropna(subset=['entrez']).entrez.duplicated().sum()))
g2['in24'] = g2.entrez.isin(dm.entrez); g2L['in24'] = g2L.entrez.isin(dm.entrez)
t = pd.DataFrame({'n': g2.groupby('next_experiment_rule_id').size(), 'in24_allMANE': g2.groupby('next_experiment_rule_id').in24.sum(),
                  'in24_LOEUFonly': g2L.groupby('next_experiment_rule_id').in24.sum()})
t.loc['ALL'] = [len(g2), g2.in24.sum(), g2L.in24.sum()]; print(t.to_string())
print('inferred CE reachable from Catalog by Entrez', ce.entrez.isin(g2.entrez.dropna().astype(int)).sum())
print('Catalog genes in 24Q4 matrix but no MANE LOEUF', int((g2.in24 & ~g2.ens.isin(L.gene_id)).sum()))

print('\n=== E. HPA')
hcols = pd.read_csv(HPA, sep='\t', nrows=0).columns
bc = [c for c in hcols if 'blood' in c.lower() or 'secretome' in c.lower()]
print('HPA columns matching blood/secretome', bc)
h = pd.read_csv(HPA, sep='\t', usecols=['Gene', 'Ensembl', 'Secretome location', 'Blood concentration - Conc. blood IM [pg/L]', 'Blood concentration - Conc. blood MS [pg/L]'])
print('rows', len(h), 'unique Ensembl', h.Ensembl.nunique(), 'Ensembl NaN', h.Ensembl.isna().sum(), 'non-ENSG', (~h.Ensembl.astype(str).str.startswith('ENSG')).sum())
print('Secretome location value counts:'); print(h['Secretome location'].value_counts(dropna=False).to_string())
sl = h['Secretome location'].dropna()
lev = np.where(sl.str.startswith('Secreted'), 'secreted', np.where(sl.eq('Intracellular and membrane'), 'intracellular_and_membrane',
               np.where(sl.eq('Immunoglobulin genes'), 'immunoglobulin', 'UNMAPPED')))
print('derived levels (non-blank text)', pd.Series(lev).value_counts().to_dict())
for c in ['Blood concentration - Conc. blood IM [pg/L]', 'Blood concentration - Conc. blood MS [pg/L]']:
    v = pd.to_numeric(h[c], errors='coerce'); print(c, 'non-blank', h[c].notna().sum(), 'numeric', v.notna().sum(), 'dtype', h[c].dtype)
gh = g.merge(h.rename(columns={'Ensembl': 'ens'}), on='ens', how='left')
cov = pd.DataFrame({'n': gh.groupby('next_experiment_rule_id').size(),
                    'LOEUF': gh.assign(x=gh.ens.isin(L.gene_id)).groupby('next_experiment_rule_id').x.sum(),
                    'HPA': gh.groupby('next_experiment_rule_id').Gene.count(),
                    'secretome_text': gh.groupby('next_experiment_rule_id')['Secretome location'].count(),
                    'blood_MS': gh.groupby('next_experiment_rule_id')['Blood concentration - Conc. blood MS [pg/L]'].count(),
                    'blood_IM': gh.groupby('next_experiment_rule_id')['Blood concentration - Conc. blood IM [pg/L]'].count()})
cov.loc['ALL'] = cov.sum(); print(cov.to_string())
print('Catalog genes not in HPA', gh.Gene.isna().sum())

print('\n=== F. GENCODE v49')
rows = []
with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.split('\t', 9)
        if p[2] != 'gene': continue
        a = p[8]
        gid = re.search(r'gene_id "([^"]+)"', a).group(1)
        nm = re.search(r'gene_name "([^"]+)"', a); hg = re.search(r'hgnc_id "([^"]+)"', a); gt = re.search(r'gene_type "([^"]+)"', a)
        rows.append((gid.split('.')[0], gid, nm.group(1) if nm else None, hg.group(1) if hg else None, p[0], gt.group(1) if gt else None))
gc = pd.DataFrame(rows, columns=['ens', 'gid_full', 'sym', 'hgnc', 'chr', 'gtype'])
prim = set([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'])
gc['primary'] = gc.chr.isin(prim)
print('gene rows', len(gc), 'unique ENSG (no version)', gc.ens.nunique(), 'with hgnc', gc.hgnc.notna().sum(), 'primary-chrom rows', gc.primary.sum(),
      'PAR_Y rows', gc.gid_full.str.contains('PAR_Y').sum())
print('first non-primary row index', int(np.argmax(~gc.primary.values)), 'last primary row index', int(np.where(gc.primary.values)[0].max()))
gcp = gc[~gc.gid_full.str.endswith('_PAR_Y')].drop_duplicates('ens')
hd = gcp.dropna(subset=['hgnc']).groupby('hgnc').ens.nunique()
print('hgnc_ids with >1 ENSG (all scaffolds)', int((hd > 1).sum()), '| on primary only', int((gcp[gcp.primary].dropna(subset=['hgnc']).groupby('hgnc').ens.nunique() > 1).sum()))
first = gcp.dropna(subset=['hgnc']).drop_duplicates('hgnc')
print('hgnc first-occurrence rows on non-primary scaffold', int((~first.primary).sum()))
cg = g.merge(gcp[['ens', 'sym', 'hgnc', 'chr', 'primary', 'gtype']], on='ens', how='left')
print('Catalog in GENCODE v49', cg.sym.notna().sum(), 'with hgnc', cg.hgnc.notna().sum(), 'on non-primary scaffold', int((cg.primary == False).sum()))

print('\n=== G. Hart CEGv2 / NEGv1 (keys only)')
for nm in ['CEGv2', 'NEGv1']:
    hh = pd.read_csv(HART + f'{nm}.txt', sep='\t'); raw = pd.read_csv(HART + f'{nm}.txt', sep='\t', dtype=str, keep_default_na=False)
    print(nm, 'rows', len(hh), 'ENTREZ dtype', hh.ENTREZ_ID.dtype, 'blank ENTREZ', int((raw.ENTREZ_ID == '').sum()), 'HGNC "NA" or blank', int(raw.HGNC_ID.isin(['NA', '']).sum()),
          'rows with blank GENE', int((raw.GENE == '').sum()))
    print('  rows with blank ENTREZ:', raw[raw.ENTREZ_ID == ''].values.tolist())
    hh['ENTREZ_ID'] = hh.ENTREZ_ID.astype('Int64')
    hm = hh.merge(first[['hgnc', 'ens', 'primary']], left_on='HGNC_ID', right_on='hgnc', how='left')
    hp = hh.merge(gcp[gcp.primary].dropna(subset=['hgnc']).drop_duplicates('hgnc')[['hgnc', 'ens']], left_on='HGNC_ID', right_on='hgnc', how='left')
    print('  HGNC->ENSG first-occurrence', hm.ens.notna().sum(), '(non-primary', int((hm.primary == False).sum()), ') | primary-only', hp.ens.notna().sum(),
          '| first vs primary ENSG disagree', int(((hm.ens != hp.ens) & hm.ens.notna()).sum()))
    print('  with Ensembl MANE LOEUF', hm.ens.isin(L.gene_id).sum(), '| in 24Q4 matrix by Entrez', hh.ENTREZ_ID.isin(dm.entrez).sum(),
          '| in Catalog', hm.ens.isin(g.ens).sum(), '| HGNC in GENCODE but symbol differs', int((hm.ens.notna() & (hm.GENE != hm.merge(gcp[['ens','sym']], on='ens', how='left').sym)).sum()))
    if nm == 'CEGv2':
        print('  CEGv2 Entrez in AchillesCommonEssentialControls (DepMap scaling controls)', int(hh.ENTREZ_ID.isin(ctl.entrez).sum()))

print('\n=== H. PC-3 and spot-check genes: keys and presence only (no values printed)')
pc3 = ['ALB', 'APOA1', 'APOB', 'FGA', 'FGB', 'FGG', 'SERPINA1', 'TF', 'HP', 'AHSG']
spot = ['PNPLA3', 'TM6SF2', 'GCKR', 'HSD17B13', 'MBOAT7', 'APOB', 'ALB', 'RPL3', 'POLR2A', 'MYC']
for lab, S in [('PC-3', pc3), ('spot', spot)]:
    for s in S:
        e = gcp[gcp.sym == s]
        ep = e[e.primary]
        ens1 = ep.ens.iloc[0] if len(ep) else None
        print(lab, s, 'GENCODE ENSG all', len(e), 'primary', len(ep), ens1, 'inCatalog', ens1 in set(g.ens), 'MANE_LOEUF', ens1 in set(L.gene_id),
              'Entrez', int(pair.set_index('gene_id').entrez.get(ens1, -1)), 'in24Q4', int(pair.set_index('gene_id').entrez.get(ens1, -1)) in set(dm.entrez),
              'inHPA', ens1 in set(h.Ensembl))

print('\n=== I. Gate provenance: CEG2/NEGv1 LOEUF with gate deciles (34,650 MANE rows), symbol join (values the prespec already states)')
gg = gn[gn.mane_select.astype(str).str.lower().isin(['true', '1'])].dropna(subset=['lof.oe_ci.upper']).copy()
gg['dec'] = pd.qcut(gg['lof.oe_ci.upper'].rank(method='first'), 10, labels=False)
print('gate MANE rows with LOEUF', len(gg), 'all-row median LOEUF', round(float(gg['lof.oe_ci.upper'].median()), 3), 'syn.oe median', round(float(gg['syn.oe'].median()), 3))
gge = gg[gg.gene_id.astype(str).str.startswith('ENSG')]
for nm in ['CEGv2', 'NEGv1']:
    hh = pd.read_csv(HART + f'{nm}.txt', sep='\t')
    x = gge[gge.gene.isin(hh.GENE)]
    xa = gg[gg.gene.isin(hh.GENE)]
    print(nm, 'symbol-matched ENSG rows', len(x), 'bottom-decile frac', round(float((x.dec == 0).mean()), 3), 'median LOEUF', round(float(x['lof.oe_ci.upper'].median()), 3),
          'syn.oe median', round(float(x['syn.oe'].median()), 3), '| all MANE rows matched', len(xa), 'frac', round(float((xa.dec == 0).mean()), 3))

print('\n=== J. Gate D3 rule on the unrecorded matrix (count only)')
old = pd.read_csv(R + 'Analysis/downstream_analysis/essentiality/CRISPRGeneEffect.csv', index_col=0)
print('old matrix lines', old.shape[0], 'genes', old.shape[1], '| genes with Chronos<-0.5 in >=90% of lines', int(((old < -0.5).mean(0) >= 0.9).sum()))
print('done')
