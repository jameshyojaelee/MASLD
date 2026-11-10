"""SS job S0d. Pileup site list for the SteatoSITE BAMs (bcftools mpileup -R), written before any SteatoSITE read,
genotype or count is used. Reads only reference panels (1kGP high-coverage phased panel, GENCODE v49, GTEx v8 liver,
dbSNP-derived identity panel) and frozen inputs; no outcome.
Sources (role column of s1b_targets_annotated.tsv):
  dig1_62   the 10 DIG1 job 62 sites (frozen/dig1/62_targets_annotated.tsv), unchanged, all base changes kept.
  panx1_71  the 7 DIG1 job 71 PANX1 sites (frozen/dig1/71_targets.tsv); alleles as in job 71b.
  tnfsf10   TNFSF10 cis tags for F-SS7 section 45. Lead = TNFSF10 GWAS lead (frozen/dig1/loci.tsv, hg19 3:172294500 G/A;
            rs79287178, GRCh38 chr3:172576710 via dbSNP b157). EUR founders of the 1kGP high-coverage panel (no parents
            listed), phased haplotypes, biallelic SNVs within 500 kb. A site qualifies if it lies in an exon of any
            GENCODE v49 TNFSF10 transcript (UTRs included), EUR MAF >= 0.05, haplotype r2 >= 0.6 with the lead, and
            is not C>T or G>A. All qualifying sites are piled up; the chosen tag is the one with highest r2 (ties:
            higher MAF, then lower position). If none qualifies, the section 45 cis test is untestable.
  identity  the A1 identity panel (in-biopsy-allelic-a1-full-20260924T025858Z/targets, 9,976 SNVs) minus C>T and G>A.
  s3h       F-SS8 universe: biallelic SNVs (1kGP high-coverage, INFO/AF_EUR 0.05-0.95) in the union of GENCODE v49
            exons of the 426 genetically_anchored or concordant r8 genes (frozen catalog), not C>T or G>A, not in an
            exon of any other GENCODE v49 gene. Whether a site is expressed is decided later from read depth
            (PRESPEC_SS_addendum_S3h.md).
  s3h_null  mapping-bias null pool: protein-coding genes outside the 426 with GTEx v8 liver median TPM >= 1 and
            liver eGene pval_beta > 0.2 (job 72a rule); same SNV, exon and base-change rules; 1,500 genes drawn with
            seed 20261025, one SNV per gene drawn with the same generator.
Output: <out>/s1b_targets.tsv (chrom, pos, ref, alt; GRCh38; chr1-22, X; sorted; one row per position) and
<out>/s1b_targets_annotated.tsv, <out>/tnfsf10_ld.tsv, <out>/summary.txt.
Usage: python -I s0d_targets.py <new_out_dir>
"""
import gzip
import random
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

assert len(sys.argv) == 2
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=False)
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
SS = ROOT / "Analysis/Hypotheses/SS-20261025T135237Z"
FZ = SS / "frozen"
BCF = "/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools"
K = ROOT / "data/external/allelic_refs/kgp_highcov"
def kvcf(c):
    """1kGP high-coverage panel file for chromosome c (the chrX file carries a .v2 suffix)."""
    return str(K / f"1kGP_high_coverage_Illumina.{c}.filtered.SNV_INDEL_SV_phased_panel{'.v2' if c == 'chrX' else ''}.vcf.gz")
GTF = "/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/gencode.v49.primary_chr.gtf"
FA = "/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/genome.primary_chr.fa"
IDENT = ROOT / "Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/targets/targets_annotated.tsv"
DEAM = {("C", "T"), ("G", "A")}
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
LEAD = ("chr3", 172576710, "G", "A")  # rs79287178; GWAS ALT-raising allele G (loci.tsv alt_raising_allele)

# verify frozen inputs
chk = subprocess.run(["sha256sum", "-c", "--quiet", "SHA256SUMS"], cwd=FZ, capture_output=True, text=True)
assert chk.returncode == 0, chk.stdout + chk.stderr

# GENCODE v49 exons
ex = []
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        r = line.split("\t", 8)
        if r[2] != "exon":
            continue
        g = re.search(r'gene_id "([^".]+)', r[8]).group(1)
        nm = re.search(r'gene_name "([^"]+)"', r[8]).group(1)
        gt = re.search(r'gene_type "([^"]+)"', r[8]).group(1)
        ex.append((r[0], int(r[3]), int(r[4]), g, nm, gt))
EX = pd.DataFrame(ex, columns=["chrom", "start", "end", "gene", "symbol", "gene_type"]).drop_duplicates()


def merged(df):
    rows = []
    for (c, g), d in df.sort_values(["chrom", "start"]).groupby(["chrom", "gene"], sort=False):
        s0 = e0 = None
        for s, e in zip(d.start, d.end):
            if s0 is None or s > e0 + 1:
                if s0 is not None:
                    rows.append((c, s0, e0, g))
                s0, e0 = s, e
            else:
                e0 = max(e0, e)
        rows.append((c, s0, e0, g))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])


def kgp_snvs(regions, tag):
    """Biallelic SNVs with AF_EUR in 0.05-0.95 inside regions (chrom, start, end, gene); one row per site x gene."""
    out = []
    for c in sorted(set(regions.chrom) & set(CHROMS), key=CHROMS.index):
        if not Path(kvcf(c)).exists():
            print("no 1kGP file for", c); continue
        rf = OUT / f"_{tag}_{c}.regions"
        regions[regions.chrom == c][["chrom", "start", "end"]].sort_values("start").to_csv(rf, sep="\t", header=False, index=False)
        q = subprocess.run(f"{BCF} view -v snps -m2 -M2 -R {rf} {kvcf(c)} 2>/dev/null | "
                           f"{BCF} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%AF_EUR\\n'", shell=True, capture_output=True,
                           text=True, check=True).stdout
        rf.unlink()
        for line in q.splitlines():
            ch, p, ref, alt, af = line.split("\t")
            try:
                af = float(af)
            except ValueError:
                continue
            if len(ref) == 1 and len(alt) == 1 and 0.05 <= af <= 0.95:
                out.append((ch, int(p), ref, alt, af))
    S = pd.DataFrame(out, columns=["chrom", "pos", "ref", "alt", "af_eur"]).drop_duplicates(["chrom", "pos"])
    # assign genes (site inside a merged exon block)
    rows = []
    for c, d in regions.groupby("chrom"):
        sp = np.sort(S[S.chrom == c].pos.values)
        for st, en, g in zip(d.start, d.end, d.gene):
            rows += [(c, int(p), g) for p in sp[np.searchsorted(sp, st, "left"):np.searchsorted(sp, en, "right")]]
    A = pd.DataFrame(rows, columns=["chrom", "pos", "gene"]).drop_duplicates()
    return S.merge(A, on=["chrom", "pos"])


EXM = merged(EX)
# exon owners of every position: count genes whose exons contain the site (for the overlap rule)
def n_exon_genes(sites):
    n = np.zeros(len(sites), int)
    for c in sites.chrom.unique():
        e = EXM[EXM.chrom == c].sort_values("start"); st = e.start.values; en = e.end.values
        idx = np.where(sites.chrom.values == c)[0]
        for i in idx:
            p = sites.pos.values[i]; k = np.searchsorted(st, p, side="right")
            lo = max(0, k - 400)
            n[i] = int(((st[lo:k] <= p) & (en[lo:k] >= p)).sum())
    return n


rows = []
# 1. DIG1 job 62
T62 = pd.read_csv(FZ / "dig1/62_targets_annotated.tsv", sep="\t", header=None,
                  names=["chrom", "start", "end", "rsid", "gene", "ref", "alt", "note"])
for r in T62.itertuples():
    rows.append(dict(chrom=r.chrom, pos=int(r.start), ref=r.ref, alt=r.alt, role="dig1_62", gene=r.gene, id=r.rsid, af_eur=np.nan, note=r.note))
# 2. PANX1 job 71 (alleles from job 71b EXP)
P71 = dict(zip([94128854, 94132981, 94137172, 94144210, 94150836, 94181312, 94181824],
               [("C", "T"), ("T", "G"), ("C", "T"), ("A", "G"), ("C", "G"), ("G", "C"), ("T", "G")]))
T71 = pd.read_csv(FZ / "dig1/71_targets.tsv", sep="\t", header=None, names=["chrom", "start", "end"])
assert set(T71.start) == set(P71)
for r in T71.itertuples():
    ref, alt = P71[int(r.start)]
    rows.append(dict(chrom=r.chrom, pos=int(r.start), ref=ref, alt=alt, role="panx1_71", gene="PANX1",
                     id=f"{r.chrom}:{r.start}", af_eur=np.nan, note="3UTR tag (F-SS6 T3)" if r.start in (94181312, 94181824) else "job 71 site"))

# 3. TNFSF10 LD
lo = pd.read_csv(FZ / "dig1/loci.tsv", sep="\t")
tl = lo[lo.genes == "TNFSF10"].iloc[0]
assert (int(tl.chr), int(tl.pos), tl.a1, tl.a2, tl.alt_raising_allele) == (3, 172294500, "G", "A", "G"), tl
ped = pd.read_csv(K / "20130606_g1k_3202_samples_ped_population.txt", sep=r"\s+", dtype=str)
eur = ped[(ped.Superpopulation == "EUR") & (ped.FatherID == "0") & (ped.MotherID == "0")].SampleID.tolist()
(OUT / "_eur.txt").write_text("\n".join(eur) + "\n")
reg = f"chr3:{LEAD[1] - 500000}-{LEAD[1] + 500000}"
q = subprocess.run(f"{BCF} view -r {reg} -v snps -m2 -M2 -S {OUT / '_eur.txt'} --force-samples {kvcf('chr3')} 2>/dev/null | "
                   f"{BCF} query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n'", shell=True, capture_output=True, text=True, check=True).stdout
(OUT / "_eur.txt").unlink()
pos, refs, alts, H = [], [], [], []
for line in q.splitlines():
    c = line.split("\t")
    g = c[3:]
    h = np.array([[int(x[0]), int(x[2])] if len(x) == 3 and x[1] == "|" and x[0] != "." and x[2] != "." else [-1, -1] for x in g]).ravel()
    if (h < 0).any():
        continue
    f = h.mean()
    if 0.01 <= f <= 0.99:
        pos.append(int(c[0])); refs.append(c[1]); alts.append(c[2]); H.append(h.astype(np.int8))
H = np.array(H); pos = np.array(pos)
li = np.where((pos == LEAD[1]) & (np.array(refs) == LEAD[2]) & (np.array(alts) == LEAD[3]))[0]
assert len(li) == 1, "lead not found in the 1kGP EUR founders"
x = H[li[0]].astype(float); Y = H.astype(float)
xc = x - x.mean(); Yc = Y - Y.mean(1, keepdims=True)
r = (Yc @ xc) / np.sqrt((Yc ** 2).sum(1) * (xc ** 2).sum())
L = pd.DataFrame(dict(chrom="chr3", pos=pos, ref=refs, alt=alts, alt_freq_eur=Y.mean(1), r=r, r2=r ** 2))
L["maf_eur"] = np.minimum(L.alt_freq_eur, 1 - L.alt_freq_eur)
tx = EX[EX.symbol == "TNFSF10"]
L["tnfsf10_exon"] = [bool(((tx.start <= p) & (tx.end >= p)).any()) for p in L.pos]
L["deamination_class"] = [(a, b) in DEAM for a, b in zip(L.ref, L.alt)]
L["qualifies"] = L.tnfsf10_exon & (L.maf_eur >= 0.05) & (L.r2 >= 0.6) & ~L.deamination_class
# r > 0: the tag ALT allele sits on the lead ALT (A) haplotype, i.e. opposite to the ALT-raising, expression-raising G
L["tag_allele_on_lead_G"] = np.where(L.r > 0, L.ref, L.alt)
L.sort_values("r2", ascending=False).to_csv(OUT / "tnfsf10_ld.tsv", sep="\t", index=False)
Q = L[L.qualifies].sort_values(["r2", "maf_eur", "pos"], ascending=[False, False, True])
n_eur = H.shape[1] // 2
for k, t in enumerate(Q.itertuples()):
    rows.append(dict(chrom="chr3", pos=int(t.pos), ref=t.ref, alt=t.alt, role="tnfsf10", gene="TNFSF10", id=f"chr3:{t.pos}",
                     af_eur=t.alt_freq_eur, note=f"r2 {t.r2:.3f} r {t.r:+.3f} with rs79287178; allele on lead G haplotype "
                     f"{t.tag_allele_on_lead_G}; {'CHOSEN TAG' if k == 0 else 'qualifying, not chosen'}"))
tn_ex = L[L.tnfsf10_exon & (L.maf_eur >= 0.05)].sort_values("r2", ascending=False)

# 4. identity panel minus C>T / G>A
I = pd.read_csv(IDENT, sep="\t")
n_id_all = len(I)
I = I[[(a, b) not in DEAM for a, b in zip(I.ref, I.alt)]]
for t in I.itertuples():
    rows.append(dict(chrom=t.chrom, pos=int(t.pos), ref=t.ref, alt=t.alt, role="identity", gene="", id="", af_eur=np.nan, note=f"maf {t.maf}"))

# 5. S3h universe
C = pd.read_parquet(FZ / "catalog/passport_gene_index.parquet")
C = C[C.primary_evidence_class.isin(["genetically_anchored", "concordant"])]
assert len(C) == 426
cat_ids = set(C.ensembl_id.astype(str))
reg_cat = EXM[EXM.gene.isin(cat_ids)]
S3 = kgp_snvs(reg_cat, "cat")
S3["n_exon_genes"] = n_exon_genes(S3)
S3 = S3[[(a, b) not in DEAM for a, b in zip(S3.ref, S3.alt)] & (S3.n_exon_genes == 1)]
sym = dict(zip(C.ensembl_id.astype(str), C.symbol.astype(str)))
for t in S3.itertuples():
    rows.append(dict(chrom=t.chrom, pos=int(t.pos), ref=t.ref, alt=t.alt, role="s3h", gene=sym[t.gene], id=t.gene, af_eur=t.af_eur, note=""))

# 6. null pool
tpm = {}
with gzip.open(ROOT / "data/external/gtex_v8_gene_median_tpm/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz", "rt") as f:
    next(f); next(f); h = next(f).rstrip("\n").split("\t"); j = h.index("Liver")
    for l in f:
        rr = l.rstrip("\n").split("\t"); tpm[rr[0].split(".")[0]] = float(rr[j])
ok = set()
with gzip.open(ROOT / "data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz", "rt") as f:
    h = next(f).rstrip("\n").split("\t"); qb = h.index("pval_beta")
    for l in f:
        rr = l.rstrip("\n").split("\t"); g = rr[0].split(".")[0]
        if float(rr[qb]) > 0.2 and tpm.get(g, 0) >= 1:
            ok.add(g)
pc = set(EX[EX.gene_type == "protein_coding"].gene)
pool_genes = sorted((ok & pc) - cat_ids)
reg_null = EXM[EXM.gene.isin(pool_genes) & EXM.chrom.isin(CHROMS)]
N = kgp_snvs(reg_null, "null")
N["n_exon_genes"] = n_exon_genes(N)
N = N[[(a, b) not in DEAM for a, b in zip(N.ref, N.alt)] & (N.n_exon_genes == 1)]
rnd = random.Random(20261025)
gs = sorted(N.gene.unique()); rnd.shuffle(gs)
gsym = dict(zip(EX.gene, EX.symbol))
for g in gs[:1500]:
    d = N[N.gene == g].sort_values("pos")
    t = d.iloc[rnd.randrange(len(d))]
    rows.append(dict(chrom=t.chrom, pos=int(t.pos), ref=t.ref, alt=t.alt, role="s3h_null", gene=gsym[g], id=g, af_eur=t.af_eur,
                     note=f"GTEx liver TPM {tpm[g]:.1f}"))

A = pd.DataFrame(rows)
A = A[A.chrom.isin(CHROMS)]
A["o"] = A.chrom.map(CHROMS.index)
A = A.sort_values(["o", "pos", "role"])
# reference check against the FASTA (every site)
bed = OUT / "_ref.regions"
A[["chrom", "pos"]].drop_duplicates().assign(e=lambda d: d.pos).to_csv(bed, sep="\t", header=False, index=False)
fa = subprocess.run(f"/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools faidx -r <(awk '{{print $1\":\"$2\"-\"$3}}' {bed}) {FA}",
                    shell=True, executable="/bin/bash", capture_output=True, text=True, check=True).stdout
bed.unlink()
refb = {}
cur = None
for line in fa.splitlines():
    if line.startswith(">"):
        cur = line[1:]
    else:
        refb[cur] = line.strip().upper()
A["fasta_ref"] = [refb[f"{c}:{p}-{p}"] for c, p in zip(A.chrom, A.pos)]
A["ref_ok"] = [fr == r[0].upper() for fr, r in zip(A.fasta_ref, A.ref)]
assert A.ref_ok.all(), A[~A.ref_ok].to_string()
A.drop(columns="o").to_csv(OUT / "s1b_targets_annotated.tsv", sep="\t", index=False)
U = A.drop_duplicates(["chrom", "pos"], keep="first")
U[["chrom", "pos", "ref", "alt"]].to_csv(OUT / "s1b_targets.tsv", sep="\t", header=False, index=False)
summ = [f"identity panel: {n_id_all} sites, {(A.role == 'identity').sum()} after removing C>T/G>A",
        f"S3h catalog genes: 426; with >= 1 site: {S3.gene.nunique()}; sites: {len(S3)}",
        f"null pool genes (non-eGene, TPM >= 1, protein-coding, not catalog): {len(pool_genes)}; with SNV: {len(gs)}; drawn: {min(1500, len(gs))}",
        f"TNFSF10: EUR founders {n_eur}; lead MAF {min(x.mean(), 1 - x.mean()):.4f}; exonic SNVs with MAF >= 0.05: {len(tn_ex)}; qualifying tags: {len(Q)}",
        "TNFSF10 exonic common SNVs (best r2 first):", tn_ex.head(15).to_string(index=False),
        "rows by role:", A.role.value_counts().to_string(), f"unique positions written: {len(U)}"]
(OUT / "summary.txt").write_text("\n".join(summ) + "\n")
print("\n".join(summ))
with open(OUT / "pip_freeze.txt", "w") as fh:
    import importlib.metadata as im
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
