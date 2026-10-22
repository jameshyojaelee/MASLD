"""NG03 (L07-03) feasibility counts. Read-only on inputs; writes only to this dir.
1. Calibrate per-site CUT&RUN depth from the GEO region count matrix against H25 bait-site census depth.
2. Predict pooled heterozygous depth at all credible-set SNVs in H3K27ac regions (1kGP AF).
3. Count testable variants and credible sets with >= 2 testable variants (within-CS PIP contrast).
4. Count chr22 caQTL positive-control variants predicted testable.
5. Simulate power of the primary within-CS PIP test.
"""
import csv, gzip, subprocess, collections, json, sys
import numpy as np
from scipy import stats

np.random.seed(20260923)
R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
H25 = R + '/Analysis/Hypotheses/H25-20261009T020111Z'
OUT = R + '/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG03'
BCF = '/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools'
KGP = R + '/data/external/allelic_refs/kgp_highcov/1kGP_high_coverage_Illumina.{}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
res = {}

# --- region counts
mat = {}
with gzip.open(R + '/Analysis/MASLD_Model_Benchmark/executions/gse267145-matrix-audit-21064333/raw/GSE267119_H3K27ac.txt.gz', 'rt') as f:
    hdr = [h.strip('"') for h in f.readline().split()]
    for line in f:
        p = line.split()
        mat[p[0].strip('"')] = np.array(p[1:], dtype=float)
res['matrix_regions'] = len(mat); res['matrix_samples'] = len(hdr)
def width(k):
    c, se = k.split(':'); s, e = se.split('-'); return int(e) - int(s) + 1
def chrom_iv(k):
    c, se = k.split(':'); s, e = se.split('-'); return c, int(s), int(e)

# --- calibration on H25 bait sites
tg = {(r['chrom'], r['pos']): r for r in csv.DictReader(open(H25 + '/stage1_targets/targets.tsv'), delimiter='\t')}
cal = []
for r in csv.DictReader(open(H25 + '/stage1_census_20261009T050357Z/site_census.tsv'), delimiter='\t'):
    reg = tg[(r['chrom'], r['pos'])]['region']
    if reg in mat:
        cal.append((float(r['mean_n']), mat[reg].mean() / width(reg)))
cal = np.array(cal)
K = np.median(cal[:, 0] / cal[:, 1])
rho = stats.spearmanr(cal[:, 0], cal[:, 1])[0]
pred = cal[:, 1] * K
res['calibration'] = dict(n_sites=len(cal), K_bp=float(K), spearman=float(rho),
                          median_abs_log2_err=float(np.median(np.abs(np.log2(pred / cal[:, 0])))),
                          observed_median_mean_n=float(np.median(cal[:, 0])))

# --- regions as interval index
regs = collections.defaultdict(list)
for k in mat:
    c, s, e = chrom_iv(k); regs[c].append((s, e, k))
for c in regs: regs[c].sort()
starts = {c: np.array([x[0] for x in regs[c]]) for c in regs}
def find_region(c, pos):
    if c not in regs: return None
    i = np.searchsorted(starts[c], pos, side='right') - 1
    if i >= 0 and regs[c][i][0] <= pos <= regs[c][i][1]: return regs[c][i][2]
    return None

# --- credible sets: hg19 id -> hg38 pos from H25 lifted bed
h38 = {}
for line in open(H25 + '/inputs/cs_all_snv_hg38.sorted.bed'):
    c, s, e, name = line.rstrip('\n').split('\t')[:4]
    h38[name.split('|')[0]] = (c, int(e))
csmem = collections.defaultdict(dict)  # cs_key -> {pos: pip}
nreliable = 0
for x in csv.DictReader(open(R + '/GWAS/finemapping/results/credible_sets.csv')):
    if x['susie_cs'] in ('', 'NA', '0') or x['susie_reliable'] != 'TRUE': continue
    if len(x['allele1']) != 1 or len(x['allele2']) != 1: continue
    vid = x['variant_id']
    if vid not in h38: continue
    nreliable += 1
    key = (x['study'], x['locus'], x['susie_cs'])
    pos = h38[vid]
    csmem[key][pos] = max(csmem[key].get(pos, 0), float(x['susie_pip']))
positions = {}
for key, m in csmem.items():
    for pos in m:
        reg = find_region(pos[0], pos[1])
        if reg: positions[pos] = reg
res['cs_snv_positions_in_regions'] = len(positions)
res['credible_sets_reliable'] = len(csmem)

# --- 1kGP AF
def af_lookup(posset, tag):
    af = {}
    for c in sorted({p[0] for p in posset}):
        bed = OUT + f'/tmp/{tag}.{c}.regions.tsv'
        with open(bed, 'w') as f:
            for cc, p in sorted(posset):
                if cc == c: f.write(f'{cc}\t{p}\n')
        out = subprocess.run([BCF, 'query', '-R', bed, '-i', 'TYPE="snp"', '-f', '%CHROM\t%POS\t%AF\n', KGP.format(c)],
                             capture_output=True, text=True).stdout
        for l in out.splitlines():
            cc, pp, a = l.split('\t')
            try: a = float(a.split(',')[0])
            except ValueError: continue
            af[(cc, int(pp))] = max(af.get((cc, int(pp)), 0), min(a, 1 - a))
    return af

def bed_regions_file(posset, tag):
    pass

af = af_lookup(set(positions), 'cs')
res['cs_positions_with_1kgp_af'] = len(af)

def pooled(pos, reg, maf, mindepth=0):
    d = mat[reg] / width(reg) * K  # predicted per-participant depth
    h = 2 * maf * (1 - maf)
    use = d >= mindepth
    return float((d[use] * h).sum()), float(h * use.sum())

rows = []
for pos, reg in positions.items():
    m = af.get(pos)
    if m is None: continue
    pd, nh = pooled(pos, reg, m)
    pd5, nh5 = pooled(pos, reg, m, 5)
    rows.append((pos, reg, m, pd, nh, pd5, nh5))
def testable(r, five=False):
    pd, nh = (r[5], r[6]) if five else (r[3], r[4])
    return pd >= 100 and nh >= 5
T = {r[0] for r in rows if testable(r)}
T5 = {r[0] for r in rows if testable(r, True)}
res['testable_positions'] = len(T); res['testable_positions_depth5callable'] = len(T5)
res['median_pooled_het_depth_all'] = float(np.median([r[3] for r in rows]))

# --- credible-set structure among testable
def cs_struct(Tset):
    n2 = 0; npairs_diffpip = 0; nvar = 0; samereg = 0; pips = []
    for key, m in csmem.items():
        t = [p for p in m if p in Tset]
        if len(t) >= 2:
            pp = [m[p] for p in t]
            if max(pp) - min(pp) >= 0.1:
                n2 += 1; nvar += len(t); pips += pp
                if len({positions[p] for p in t}) < len(t): samereg += 1
    return dict(cs_ge2_testable_piprange_ge01=n2, variants_in_them=nvar, cs_with_two_in_same_region=samereg, n_highpip_ge01=int(sum(p >= 0.1 for p in pips)))
res['cs_struct'] = cs_struct(T); res['cs_struct_depth5'] = cs_struct(T5)
res['testable_pip_ge0.1'] = int(sum(1 for p in T if max(m.get(p, 0) for m in csmem.values() if p in m) >= 0.1))

# --- chr22 caQTL positive control
cq = {}
with gzip.open(R + '/data/external/currin_2025_caqtl/caQTL_hg38_chr22_with-varBeta-and-MAF_FastQTL-WASP-normal-nominal-updatedConsensusPeaks-25ATAC-PCs-minMAC10-results_1mbFromPeakCenters.txt.gz', 'rt') as f:
    next(f)
    for l in f:
        v, pk, dist, p, b, vb, r2, maf = l.rstrip('\n').split('\t')
        if float(p) < 1e-6 and abs(float(dist)) < 500:
            c, pos, a1, a2 = v.split(':')
            if len(a1) == 1 and len(a2) == 1:
                reg = find_region(c, int(pos))
                if reg: cq[(c, int(pos))] = (reg, float(maf))
cq_t = [k for k, (reg, maf) in cq.items() if pooled(k, reg, maf)[0] >= 100 and pooled(k, reg, maf)[1] >= 5]
chr22_cs = sum(1 for p in positions if p[0] == 'chr22')
res['caqtl_chr22_in_regions'] = len(cq); res['caqtl_chr22_testable'] = len(cq_t)
res['chr22_share_of_cs_positions'] = chr22_cs

json.dump(res, open(OUT + '/gate_counts_dedup_rerun.json', 'w'), indent=1, default=str)
print(json.dumps(res, indent=1, default=str), flush=True)
# --- power simulation of primary test (within-CS weighted regression of |logit ratio| on PIP, permutation within CS)
PD = {r[0]: r[3] for r in rows}
sets=[]
for key,m in csmem.items():
    t=sorted(p for p in m if p in T)
    if len(t)>=2:
        pp=[m[p] for p in t]
        if max(pp)-min(pp)>=0.1: sets.append((key,tuple(t)))
uniq=set(s for _,s in sets)
loci=[]
for s in sorted(uniq):
    hit=[L for L in loci if any(a[0]==b[0] and abs(a[1]-b[1])<500000 for a in s for b in L)]
    new=set(s)
    for L in hit: new|=L; loci.remove(L)
    loci.append(new)
print('CS',len(sets),'distinct variant sets',len(uniq),'independent loci (500kb)',len(loci))
print('distinct testable variants in primary', len(set(v for s in uniq for v in s)))
# reduced power sim on distinct variant sets; tagging: non-causal variant in the same region as a causal one carries its ratio with prob tag
pipof={}
for key,m in csmem.items():
    for p,v in m.items(): pipof.setdefault(tuple(sorted(m)),{})
CSu={}
for key,s in sets:
    m=csmem[key]; CSu.setdefault(s,[]).append([m[p] for p in s])
CS2=[(np.array(np.mean(v,axis=0)), np.array([PD[p] for p in s]), [positions[p] for p in s]) for s,v in CSu.items()]
rng=np.random.default_rng(20260923)
def st(CSd):
    num=den=0
    for pp,y,w in CSd:
        pc=pp-np.average(pp,weights=w); yc=y-np.average(y,weights=w); num+=(w*pc*yc).sum(); den+=(w*pc*pc).sum()
    return num/den
def sim(q,ratio,tag,od=0.01,nperm=100):
    CSd=[]
    for pp,dd,rg in CS2:
        c=rng.random(len(pp))<pp*q; r=np.where(c,ratio,0.5)
        for i in range(len(pp)):
            if not c[i] and any(c[j] and rg[j]==rg[i] for j in range(len(pp))) and rng.random()<tag: r[i]=ratio
        n=np.maximum(dd.round().astype(int),1); x=rng.binomial(n,rng.beta(r*(1-od)/od,(1-r)*(1-od)/od))
        f=(x+.5)/(n+1); CSd.append((pp,np.abs(np.log(f/(1-f))),n/(1+(n-1)*od)))
    s0=st(CSd); nl=[st([(rng.permutation(pp),y,w) for pp,y,w in CSd]) for _ in range(nperm)]
    return (1+sum(v>=s0 for v in nl))/(nperm+1)
for ratio,q,tag in [(0.65,0.25,0),(0.65,0.5,0),(0.65,1.0,0),(0.65,0.5,0.8),(0.65,1.0,0.8),(0.7,0.25,0),(0.7,0.5,0.8)]:
    ps=[sim(q,ratio,tag) for _ in range(40)]
    print('distinct sets',len(CS2),'ratio',ratio,'q',q,'tag',tag,'power',np.mean(np.array(ps)<0.05),flush=True)
