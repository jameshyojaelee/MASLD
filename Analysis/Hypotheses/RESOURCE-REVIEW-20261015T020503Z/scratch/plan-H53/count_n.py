# H53 prespec: count units from the r2 COLOC release. Reads no caQTL data and no outcome.
import pandas as pd, numpy as np
B = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
R2 = B + "/results/susie_coloc_r2_20260924T053244Z"
t = pd.read_csv(R2 + "/gene_level_coloc_tier12.csv")
pos = set(t.loc[t.susie_state_t12 == "positive", "ensembl"])
print("primary SuSiE genes (susie_state_t12 == positive):", len(pos))
s = pd.read_csv(R2 + "/susie_coloc_signal_pairs.tsv.gz", sep="\t")
tier = pd.read_csv(B + "/config/gwas_trait_tier.tsv", sep="\t")[["study_name", "tier"]]
reg = pd.read_csv(B + "/config/gwas_registry.tsv", sep="\t")[["study_name", "ancestry"]]
s = s.merge(tier, left_on="gwas", right_on="study_name", how="left").merge(
    reg, left_on="gwas", right_on="study_name", how="left", suffixes=("", "_r"))
k = s[(s["PP.H4"] > 0.5) & (s.eligible == True)]
kp = k[k.ensembl.isin(pos)]
print("eligible PP.H4>0.5 pairs on primary genes:", len(kp), "genes", kp.ensembl.nunique())
L = kp[kp.tier.isin([1, 2]) & (kp.ancestry == "EUR")]
print("LIVER ARM tier1/2 EUR: pairs", len(L), "genes", L.ensembl.nunique(),
      "distinct hit1", L.hit1.nunique(), "distinct (study,hit1)", L[["gwas", "hit1"]].drop_duplicates().shape[0])
print("primary genes with no tier1/2 EUR pair:", len(pos - set(L.ensembl)))
# platelet comparator: same selection (eligible, PP.H4>0.5), any gene (not restricted to liver-primary genes)
P = k[k.gwas == "MVP_Platelet_EUR"]
print("PLATELET EUR eligible PP.H4>0.5: pairs", len(P), "genes", P.ensembl.nunique(), "distinct hit1", P.hit1.nunique())

def split(h):
    c, p = h.split(":"); return c, int(p)

def cluster(leads, gap=5e5):
    d = pd.DataFrame([split(h) for h in sorted(set(leads))], columns=["chr", "pos"])
    d["chrn"] = d.chr.astype(int); d = d.sort_values(["chrn", "pos"])
    cid, c, prev = [], -1, (None, None)
    for ch, p in zip(d.chr, d.pos):
        if ch != prev[0] or p - prev[1] > gap: c += 1
        cid.append(c); prev = (ch, p)
    d["locus"] = cid; d["hit1"] = d.chr + ":" + d.pos.astype(str); return d

dl = cluster(L.hit1)
print("liver loci (500 kb):", dl.locus.nunique(), "| 250 kb:", cluster(L.hit1, 2.5e5).locus.nunique(),
      "| 1 Mb:", cluster(L.hit1, 1e6).locus.nunique())
lv = dl.groupby("chr").pos.apply(np.array).to_dict()
def far(h):
    c, p = split(h); a = lv.get(c); return a is None or np.all(np.abs(a - p) > 5e5)
pl = sorted(set(P.hit1)); po = [h for h in pl if far(h)]
dp = cluster(po)
print("platelet leads", len(pl), "| >500 kb from every liver lead:", len(po), "| platelet-only loci (500 kb):", dp.locus.nunique())
# per-study distinct leads in the liver arm
print("liver arm distinct leads by study:")
print(L.groupby("gwas").hit1.nunique().sort_values(ascending=False).to_string())
# confirmatory set without chr1 and chr19 (pilot chromosomes)
print("liver loci excluding chr1, chr19:", dl[~dl.chr.isin(["1", "19"])].locus.nunique(),
      "| platelet-only loci excluding chr1, chr19:", dp[~dp.chr.isin(["1", "19"])].locus.nunique())
# same-biobank sensitivity
m = cluster(L.loc[L.gwas == "MVP_ALT_EUR", "hit1"])
print("MVP_ALT_EUR loci (500 kb):", m.locus.nunique(), "| excl chr1/19:", m[~m.chr.isin(["1", "19"])].locus.nunique())
# studies per liver locus
x = L[["gwas", "hit1"]].drop_duplicates().merge(dl[["hit1", "locus"]], on="hit1")
ns = x.groupby("locus").gwas.nunique()
print("studies per liver locus: mean %.2f, loci with >=2 studies %d, max %d" % (ns.mean(), (ns >= 2).sum(), ns.max()))
# tier split of the liver arm
print("liver arm pairs by tier:", L.tier.value_counts().to_dict())
# positive-control loci in the arm (hg19 windows around the published caQTL leads)
for name, ch, lo, hi in [("caPeak21014 chr1 RALGPS2/ANGPTL1 region", "1", 177.9e6, 179.1e6),
                         ("caPeak156441 chr19 region", "19", 41.3e6, 42.4e6)]:
    sub = L[L.hit1.str.startswith(ch + ":")]
    sub = sub[sub.hit1.map(lambda h: lo <= split(h)[1] <= hi)]
    print(name, "liver-arm rows:", sub[["gene", "gwas", "hit1", "idx1", "idx2", "PP.H4"]].to_string(index=False) if len(sub) else "none")
