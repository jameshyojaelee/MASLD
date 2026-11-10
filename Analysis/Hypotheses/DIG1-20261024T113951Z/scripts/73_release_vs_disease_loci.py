# DIG1 job 73 (post hoc, after job 72 / section 44). Do liver-enzyme loci split into loci that raise ALT without raising
# liver disease ("release loci") and loci that raise both ("disease loci")?
# Read before this script: at PANX1, COLOC PP.H4 0.97 with MVP ALT EUR, but PP.H2 0.45-0.9 for every MASLD, NASH,
# cirrhosis, PDFF and GGT GWAS (eQTL present, no disease signal). Allele convention checked at PNPLA3 rs738409 and
# TM6SF2 rs58542926 in MVP, FinnGen, UKBB and Pazoki files: beta refers to allele1 in all of them (G and T raise ALT,
# NAFLD, PDFF). No disease z at any ALT locus was read before this header was written.
# Plan (fixed here):
# 1. Universe: genes with max(PP.H4.abf, PP.H4.susie) > 0.8 for MVP_ALT_EUR or MVP_AST_EUR (COLOC r2). Lead SNP per
#    gene = top_snp of its ALT row (AST row if no ALT row passes). Lead SNPs within 500 kb merge into one locus; the
#    locus SNP is the one with largest |z| in MVP ALT EUR. Unit = locus (genes in a locus are not independent).
#    Benchmarks always added: PNPLA3 rs738409 22:44324727, TM6SF2 rs58542926 19:19379549, GCKR rs1260326 2:27730940.
# 2. z = beta/se for allele1, aligned to MVP ALT EUR allele1 (flip on swapped alleles, drop on mismatch), then
#    oriented to the ALT-raising allele. Traits: MVP ALT/AST/NAFLD/Cirrhosis/ChronLiver/Albumin/Platelet EUR, UKBB
#    ALT/AST/GGT, FinnGen NAFLD/NASH, Pazoki PDFF.
# 3. Classes on MVP (same participants, so z ratios compare): k = mean z_NAFLD/z_ALT of PNPLA3 and TM6SF2.
#    disease: z_NAFLD >= 3.29. opposite: z_NAFLD <= -3.29.
#    release: z_ALT >= 5.45, |z_NAFLD| < 1.96 and k*z_ALT - z_NAFLD >= 3. Otherwise intermediate.
# 4. Locus annotations: any colocalized gene in CCL20_STATE_UP (job 69 state.sets), hepatocyte-made (job 68a median
#    hepatocyte share >= 0.85 in all 3 cohorts), bulk disease up (release deg_results padj < 0.05, logFC > 0).
#    Tests: release vs disease loci, Fisher exact for each annotation; family = these 3 tests, Holm.
#    Descriptive: PDFF, cirrhosis, GGT, AST, albumin, platelet z by class.
# Unit: GWAS participants (summary statistics); loci as units for the class comparison. Cross-sectional.
# Usage: python -I 73_release_vs_disease_loci.py <new_out_dir>
import sys, os, gzip
import numpy as np, pandas as pd
from scipy.stats import fisher_exact

out = sys.argv[1]
assert not os.path.exists(out)
os.makedirs(out)
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN = f"{ROOT}/Analysis/Hypotheses/DIG1-20261024T113951Z"
SS = f"{ROOT}/GWAS/finemapping/data/sumstats"
COLOC = f"{ROOT}/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
DEG = f"{ROOT}/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
TRAITS = ["MVP_ALT_EUR", "MVP_AST_EUR", "MVP_NAFLD_EUR", "MVP_Cirrhosis_EUR", "MVP_ChronLiver_EUR", "MVP_Albumin_EUR",
          "MVP_Platelet_EUR", "UKBB_ALT", "UKBB_AST", "UKBB_GGT", "FinnGen_NAFLD", "FinnGen_NASH", "Pazoki_PDFF"]
BENCH = {"PNPLA3_rs738409": (22, 44324727), "TM6SF2_rs58542926": (19, 19379549), "GCKR_rs1260326": (2, 27730940)}

# 1. universe
C = pd.read_csv(COLOC, usecols=["gwas_name", "gene", "ensembl", "PP.H4.abf", "PP.H4.susie", "top_snp"])
C = C[C.gwas_name.isin(["MVP_ALT_EUR", "MVP_AST_EUR"])].copy()
C["h4"] = C[["PP.H4.abf", "PP.H4.susie"]].max(axis=1, skipna=True)
C = C[C.h4 > 0.8].copy()
C["pref"] = (C.gwas_name != "MVP_ALT_EUR").astype(int)
C = C.sort_values(["gene", "pref", "h4"], ascending=[True, True, False]).drop_duplicates("gene")
sp = C.top_snp.astype(str).str.replace("chr", "").str.split(":", expand=True)
C["chr"] = pd.to_numeric(sp[0], errors="coerce"); C["pos"] = pd.to_numeric(sp[1].str.split("_").str[0], errors="coerce")
print("coloc genes:", len(C), "unparsed top_snp:", C.chr.isna().sum())
C = C.dropna(subset=["chr", "pos"]).astype({"chr": int, "pos": int})
want = set(zip(C.chr, C.pos)) | set(BENCH.values())

# 2. extract z
def extract(t):
    keep = []
    for ch in pd.read_csv(f"{SS}/{t}_reformatted_hg19.tsv", sep="\t", chunksize=2_000_000,
                          usecols=["chromosome", "position", "allele1", "allele2", "beta", "se"],
                          dtype={"chromosome": str}):
        ch["chromosome"] = pd.to_numeric(ch.chromosome.str.replace("chr", ""), errors="coerce")
        ch = ch.dropna(subset=["chromosome"])
        m = [(int(a), int(b)) in want for a, b in zip(ch.chromosome, ch.position)]
        keep.append(ch[m])
    d = pd.concat(keep).rename(columns={"chromosome": "chr", "position": "pos"})
    d["chr"] = d.chr.astype(int); d["z"] = d.beta / d.se
    d["allele1"] = d.allele1.str.upper(); d["allele2"] = d.allele2.str.upper()
    return d.drop_duplicates(["chr", "pos", "allele1", "allele2"])

Z = {t: extract(t) for t in TRAITS}
ref = Z["MVP_ALT_EUR"].drop_duplicates(["chr", "pos"], keep=False)[["chr", "pos", "allele1", "allele2", "z"]]
ref = ref.rename(columns={"allele1": "a1", "allele2": "a2", "z": "z_MVP_ALT_EUR"})
W = ref.copy()
for t in TRAITS[1:]:
    d = Z[t].merge(ref[["chr", "pos", "a1", "a2"]], on=["chr", "pos"])
    same = (d.allele1 == d.a1) & (d.allele2 == d.a2); swap = (d.allele1 == d.a2) & (d.allele2 == d.a1)
    d = d[same | swap].copy(); d.loc[swap[same | swap].values, "z"] *= -1
    d = d.drop_duplicates(["chr", "pos"], keep=False)
    W = W.merge(d[["chr", "pos", "z"]].rename(columns={"z": f"z_{t}"}), on=["chr", "pos"], how="left")
sgn = np.sign(W.z_MVP_ALT_EUR)
for c in [c for c in W.columns if c.startswith("z_")]:
    W[c] = W[c] * sgn
W["alt_raising_allele"] = np.where(sgn > 0, W.a1, W.a2)

# loci
G = C.merge(W, on=["chr", "pos"], how="left")
print("coloc genes with MVP ALT z at lead:", G.z_MVP_ALT_EUR.notna().sum(), "of", len(G))
G = G.dropna(subset=["z_MVP_ALT_EUR"]).sort_values(["chr", "pos"]).reset_index(drop=True)
lid, prev, k = [], (None, None), -1
for ch, ps in zip(G.chr, G.pos):
    if ch != prev[0] or ps - prev[1] > 500_000: k += 1
    lid.append(k); prev = (ch, ps)
G["locus"] = lid
hs = pd.read_csv(f"{RUN}/out_hep_share/gene_hep_share.tsv.gz", sep="\t")
hs = hs.groupby("gene").agg(mn=("median_hep_share", "min"), k=("cohort", "nunique"))
HEP = set(hs[(hs.mn >= 0.85) & (hs.k == 3)].index)
STATE = set(open(f"{RUN}/out_state_magma/state.sets").readline().split()[1:])
deg = pd.read_csv(DEG); deg["ens"] = deg.gene.str.split(".").str[0]
UP = set(deg[(deg.padj < 0.05) & (deg.logFC > 0)].ens)
G["state_up"] = G.ensembl.isin(STATE); G["hep_made"] = G.gene.isin(HEP); G["bulk_up"] = G.ensembl.isin(UP)
G.to_csv(f"{out}/coloc_genes_lead_z.tsv", sep="\t", index=False)
lead = G.loc[G.groupby("locus").z_MVP_ALT_EUR.idxmax()]
agg = G.groupby("locus").agg(genes=("gene", lambda x: ",".join(sorted(set(x)))), n_genes=("gene", "nunique"),
                             state_up=("state_up", "any"), hep_made=("hep_made", "any"), bulk_up=("bulk_up", "any"))
zc = [c for c in W.columns if c.startswith("z_")]
L = lead[["locus", "chr", "pos", "a1", "a2", "alt_raising_allele"] + zc].merge(agg, on="locus")
B = W[[(a, b) in set(BENCH.values()) for a, b in zip(W.chr, W.pos)]].copy()
B["genes"] = [next(k for k, v in BENCH.items() if v == (a, b)) for a, b in zip(B.chr, B.pos)]
print("benchmarks:\n", B[["genes"] + zc].round(1).to_string(index=False))
# 3. classes
kk = float(np.mean([B.loc[B.genes == g, "z_MVP_NAFLD_EUR"].iloc[0] / B.loc[B.genes == g, "z_MVP_ALT_EUR"].iloc[0]
                    for g in ["PNPLA3_rs738409", "TM6SF2_rs58542926"]]))
za, zn = L.z_MVP_ALT_EUR, L.z_MVP_NAFLD_EUR
L["cls"] = np.select([zn >= 3.29, zn <= -3.29, (za >= 5.45) & (zn.abs() < 1.96) & (kk * za - zn >= 3)],
                     ["disease", "opposite", "release"], "intermediate")
L.loc[zn.isna(), "cls"] = "no_disease_z"
L["ratio_nafld_alt"] = zn / za
L.to_csv(f"{out}/loci.tsv", sep="\t", index=False)
print(f"k (benchmark z_NAFLD/z_ALT) = {kk:.3f}; loci = {len(L)}")
print(L.cls.value_counts().to_string())
# 4. tests
rows = []
for a in ["state_up", "hep_made", "bulk_up"]:
    r = L[L.cls == "release"][a]; d = L[L.cls == "disease"][a]
    tab = [[int(r.sum()), int((~r).sum())], [int(d.sum()), int((~d).sum())]]
    o, p = fisher_exact(tab)
    rows.append(dict(annotation=a, release_yes=tab[0][0], release_n=len(r), disease_yes=tab[1][0], disease_n=len(d),
                     odds_ratio=o, p=p))
T = pd.DataFrame(rows).sort_values("p")
m = len(T); T["p_holm"] = np.minimum(1, np.maximum.accumulate((m - np.arange(m)) * T.p.values))
T.to_csv(f"{out}/class_annotation_tests.tsv", sep="\t", index=False)
print(T.to_string(index=False))
D = L.groupby("cls")[zc].median().round(2).T
D.to_csv(f"{out}/class_median_z.tsv", sep="\t")
print(D.to_string())
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
show = ["genes", "z_MVP_ALT_EUR", "z_MVP_AST_EUR", "z_MVP_NAFLD_EUR", "z_MVP_Cirrhosis_EUR", "z_Pazoki_PDFF", "z_UKBB_GGT",
        "state_up", "hep_made", "bulk_up"]
for c in ["release", "opposite", "disease"]:
    print(f"\n== {c}"); print(L[L.cls == c].sort_values("z_MVP_ALT_EUR", ascending=False)[show].round(1).to_string(index=False))
os.system(f"/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/pip freeze > {out}/pip_freeze.txt")
