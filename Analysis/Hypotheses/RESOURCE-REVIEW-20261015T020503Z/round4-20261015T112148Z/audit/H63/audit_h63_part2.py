# H63 audit part 2 (after job 9763 stopped at the xlsx read: openpyxl is only in ~/.local).
# Uses detection_fraction.tsv.gz, which job 9763 reproduced to 1.1e-16. Parses the xlsx with the
# standard library. Computes NO zonation statistic (no COM, slope, or mRNA-protein correlation).
import json, re, zipfile
import h5py, numpy as np, pandas as pd

ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
R4 = ROOT + '/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z'
OUT = R4 + '/audit/H63'
SC = ROOT + '/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/scdvp'
CATP = ROOT + '/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
V = ROOT + '/Analysis/Spatial/data/public_expansion_2026/yakubovsky2026/v.mat'
SCH = ROOT + '/Analysis/Spatial/candidates/program-context-v2-candidate-2026-08-07/acquisition/yakubovsky2026/yakubovsky_hdf5_join_schema.tsv'

# ---- xlsx Table 3 vs TSV ----
z = zipfile.ZipFile(SC + '/42255_2026_1459_MOESM3_ESM.xlsx')
ss = [re.sub(r'<[^>]+>', '', m) for m in re.findall(r'<si>(.*?)</si>', z.read('xl/sharedStrings.xml').decode(), flags=re.S)]
rels = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]+)"', z.read('xl/_rels/workbook.xml.rels').decode()))
rels.update({a: b for b, a in re.findall(r'Target="([^"]+)"[^>]*Id="(rId\d+)"', z.read('xl/_rels/workbook.xml.rels').decode())})
sheets = re.findall(r'<sheet name="([^"]+)"[^>]*r:id="(rId\d+)"', z.read('xl/workbook.xml').decode())
target = [rels[r] for n, r in sheets if n == 'Table 3'][0]
xml = z.read('xl/' + target.lstrip('/').replace('xl/', '')).decode()
def col2i(c):
    n = 0
    for ch in c: n = n * 26 + ord(ch) - 64
    return n - 1
rows = []
for rx in re.findall(r'<row[^>]*>(.*?)</row>', xml, flags=re.S):
    r = {}
    for ref, attrs, body in re.findall(r'<c r="([A-Z]+)\d+"([^>]*)>(.*?)</c>', rx, flags=re.S):
        v = re.search(r'<v>(.*?)</v>', body)
        isv = re.search(r'<is>.*?<t[^>]*>(.*?)</t>', body, flags=re.S)
        if 't="s"' in attrs and v: val = ss[int(v.group(1))]
        elif isv: val = isv.group(1)
        elif v: val = v.group(1)
        else: val = None
        r[col2i(ref)] = val
    rows.append(r)
print('Table 3 xml rows', len(rows))
for i in range(4):
    print('row', i, [rows[i].get(k) for k in sorted(rows[i])][:6], '...', [rows[i].get(k) for k in sorted(rows[i])][-6:])
hdr_i = next(i for i, r in enumerate(rows) if any(str(v).strip() == 'coefficient' for v in r.values()))
hdr = rows[hdr_i]; print('header row', hdr_i, [hdr.get(k) for k in sorted(hdr)])
names = {k: str(v).strip() for k, v in hdr.items() if v is not None}
print('header keys', sorted(hdr)[:3], 'column A in header', 0 in hdr)
names[0] = 'protein'  # column A carries the protein id and has no header text
data = pd.DataFrame([{names.get(k, 'col%d' % k): v for k, v in r.items()} for r in rows[hdr_i + 1:] if r])
print('xlsx data rows', len(data), 'columns', list(data.columns)[:3], list(data.columns)[-6:])
tsv = pd.read_csv(SC + '/scdvp_table3_healthy.tsv', sep='\t')
pcol = [c for c in data.columns if str(c).lower().startswith('protein') or str(c).lower() in ('gene', 'genes', 'gene names')]
print('protein column candidates', pcol)
x = data.rename(columns={pcol[0]: 'protein'}) if pcol else data
x = x[x.protein.notna()]
print('xlsx proteins', len(x), 'tsv proteins', len(tsv), 'same id order', list(x.protein.astype(str)) == list(tsv.protein.astype(str)))
for c in ['coefficient', 'p_value', 'q_value']:
    if c in x.columns:
        a = pd.to_numeric(x[c], errors='coerce').to_numpy(); b = tsv[c].to_numpy()
        print(c, 'max abs diff', float(np.nanmax(np.abs(a - b))) if len(a) == len(b) else 'length differs', 'NA xlsx', int(np.isnan(a).sum()))
for c in ['zonated', 'strong_zonated']:
    if c in x.columns:
        print(c, 'identical', list(x[c].astype(str)) == list(tsv[c].astype(str)), x[c].value_counts().to_dict())

# ---- set sizes by class from the reproduced detection fractions ----
LMj = json.load(open(R4 + '/write/H63/landmarks.json')); LM = set(LMj['central_LM']) | set(LMj['portal_LM'])
D = pd.read_csv(R4 + '/write/H63/detection_fraction.tsv.gz', sep='\t', index_col=0)
cat = pd.read_parquet(CATP)
pro = tsv[~tsv.protein.duplicated(keep=False)]
psym = set(pro.protein); pzon = set(pro.protein[pro.zonated == 'Yes']); pstr = set(pro.protein[pro.strong_zonated == 'Yes'])
sets = {'primary_M4_M8': ['M4', 'M5', 'M6', 'M7', 'M8'], 'sens_M2_M8': ['M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8']}
out = []; T1 = {}
for sn, dl in sets.items():
    robust = set(D.index[(D[dl] >= 0.05).all(axis=1)])
    for cl in ['ALL', 'genetically_anchored', 'established_state_associated', 'concordant', 'unresolved', 'untested']:
        s = set((cat if cl == 'ALL' else cat[cat.primary_evidence_class == cl]).symbol)
        inV = s & set(D.index)
        out.append(dict(donor_set=sn, cls=cl, n_catalog=len(s), in_visium=len(inV), robust=len(inV & robust),
                        robust_nonLM=len((inV & robust) - LM), protein=len(s & psym), T1=len((inV & robust & psym) - LM),
                        T1_zonated=len((inV & robust & pzon) - LM), T1_strong=len((inV & robust & pstr) - LM),
                        robust_protein_LM=len(inV & robust & psym & LM)))
        if cl == 'ALL': T1[sn] = (inV & robust & psym) - LM
out = pd.DataFrame(out); out.to_csv(OUT + '/recount_sets.tsv', sep='\t', index=False); print(out.to_string())
print('T1 sets identical primary vs sensitivity', T1['primary_M4_M8'] == T1['sens_M2_M8'])
print('catalog&protein', len(set(cat.symbol) & psym), 'of which in Visium', len(set(cat.symbol) & psym & set(D.index)))
rob = set(D.index[(D[sets['primary_M4_M8']] >= 0.05).all(axis=1)])
print('LM with protein', len(LM & psym), 'LM robust primary', len(LM & rob), 'LM robust with protein', len(LM & rob & psym))
print('LM robust with protein list', sorted(LM & rob & psym))
t1 = sorted(T1['primary_M4_M8'])
print('T1 detection min per donor', D.loc[t1, sets['primary_M4_M8']].min().round(4).to_dict(),
      'median', D.loc[t1, sets['primary_M4_M8']].median().round(4).to_dict())
# hepatocyte indicator split of T1 (for D4 sizes)
sch = pd.read_csv(SCH, sep='\t'); f = h5py.File(V, 'r')
s = sch[sch.patient == 'M4'].set_index('field').hdf5_path; root = s['mat'].rsplit('/', 1)[0]
genes = [''.join(chr(c) for c in f[r][()].ravel()) for r in f[s['gene_name']][()].ravel()]
hep = pd.Series(f[root + '/indicator_hep_gene'][()].ravel(), index=genes)
same = all(np.array_equal(f[sch[sch.patient == d].set_index('field').hdf5_path['mat'].rsplit('/', 1)[0] + '/indicator_hep_gene'][()].ravel(), hep.values) for d in ['M1','M2','M3','M5','M6','M7','M8'])
hep = hep[~hep.index.duplicated(keep=False)]
print('indicator_hep_gene identical across donors', same, 'values', hep.value_counts().to_dict())
print('T1 hepatocyte split', int((hep.reindex(t1) == 1).sum()), 'hep;', int((hep.reindex(t1) != 1).sum()), 'other')
print('done')
