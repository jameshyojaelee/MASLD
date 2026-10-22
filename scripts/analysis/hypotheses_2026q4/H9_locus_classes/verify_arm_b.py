# Independent re-derivation of arm B program and cell-type counts and Fisher P (pure stdlib).
import csv, sys, math, collections
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
HOT = ROOT + '/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot'
HM = ROOT + '/Analysis/SingleCell/results_gpu_v2/hotspot_modules'
res_dir = sys.argv[1]
D = ['ABO', 'ADH1B', 'APOE', 'GPAM', 'GUSB', 'MTARC1', 'MBOAT7', 'TMC4', 'PNPLA3', 'TECTB', 'TM6SF2']
C = ['C5orf67', 'CEBPA', 'COBLL1', 'GCKR', 'ATPAF2', 'GID4', 'MYO15A', 'HFE', 'IRS1', 'JAZF1', 'MTTP', 'PNPLA2', 'PPARG', 'TRIB1', 'TSC22D2']
cts = ['hepatocytes', 'fibroblasts', 'cholangiocytes', 'macrophages', 'tcells']
univ = {}
for ct in cts:
    with open(f'{HM}/{ct}/autocorr.tsv') as fh:
        r = csv.DictReader(fh, delimiter='\t'); univ[ct] = {x['gene'] for x in r}
mem = list(csv.DictReader(open(HOT + '/program_membership_v2.tsv'), delimiter='\t'))
prog = collections.defaultdict(lambda: dict(src=set(), keys=set()))
ctkeys = collections.defaultdict(set); ctsrc = collections.defaultdict(set); ctuniv = {ct: set(univ[ct]) for ct in cts}
for m in mem:
    pid = m['cell_type'] + '_' + m['module']
    prog[pid]['src'].add(m['source_gene'])
    ks = {m['source_gene'], m['mapped_symbol'], m['canonical_gene']} - {'', 'NA'}
    prog[pid]['keys'] |= ks; ctkeys[m['cell_type']] |= ks; ctsrc[m['cell_type']].add(m['source_gene'])
    ctuniv[m['cell_type']] |= {m['mapped_symbol'], m['canonical_gene']} - {'', 'NA'}
def pval(X, K, N, n):
    if X <= 0: return 1.0
    tot = math.comb(N, n)
    return sum(math.comb(K, i) * math.comb(N - K, n - i) for i in range(X, min(K, n) + 1)) / tot
rows = {}
for name, genes in (('discordant', D), ('concordant', C)):
    for pid, d in prog.items():
        ct = pid.rsplit('_', 1)[0]
        n = sum(1 for g in genes if g in ctuniv[ct]); X = sum(1 for g in genes if g in d['keys'])
        rows[(name, 'program', pid)] = (len(d['src']), len(univ[ct]), n, X, pval(X, len(d['src']), len(univ[ct]), n))
    for ct in cts:
        n = sum(1 for g in genes if g in ctuniv[ct]); X = sum(1 for g in genes if g in ctkeys[ct])
        K = len(ctsrc[ct]); N = len(univ[ct])
        rows[(name, 'celltype', ct)] = (K, N, n, X, pval(X, K, N, n))
got = {}
with open(res_dir + '/program_tests.tsv') as fh:
    for x in csv.DictReader(fh, delimiter='\t'):
        if x['arm'] == 'B': got[(x['class'], x['level'], x['id'])] = (int(x['K']), int(x['N']), int(float(x['n_class_genes_in_universe'])), int(float(x['X'])), float(x['fisher_p']))
bad = 0
for k, v in rows.items():
    g = got.get(k)
    if g is None or g[:4] != v[:4] or abs(g[4] - v[4]) > 1e-9 * max(1, v[4]):
        bad += 1; print('MISMATCH', k, v, g)
print('compared', len(rows), 'mismatches', bad)
top = sorted(rows.items(), key=lambda kv: kv[1][4])[:8]
for k, v in top: print(k, v)
