# Audit H53: symmetric arm counts and gene-filter check. Reads only r2 tables and configs (no caQTL outcome).
import pandas as pd, numpy as np
B = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
R2 = B + "/results/susie_coloc_r2_20260924T053244Z"
t = pd.read_csv(R2 + "/gene_level_coloc_tier12.csv")
pos = set(t.loc[t.susie_state_t12 == "positive", "ensembl"])
s = pd.read_csv(R2 + "/susie_coloc_signal_pairs.tsv.gz", sep="\t")
tier = pd.read_csv(B + "/config/gwas_trait_tier.tsv", sep="\t")[["study_name", "tier"]]
reg = pd.read_csv(B + "/config/gwas_registry.tsv", sep="\t")[["study_name", "ancestry"]]
s = s.merge(tier, left_on="gwas", right_on="study_name", how="left").merge(
    reg, left_on="gwas", right_on="study_name", how="left", suffixes=("", "_r"))
k = s[(s["PP.H4"] > 0.5) & (s.eligible == True)]
Lall = k[k.tier.isin([1, 2]) & (k.ancestry == "EUR")]
L = Lall[Lall.ensembl.isin(pos)]
print("tier1/2 EUR pairs without primary-gene filter:", len(Lall), "| with:", len(L))
P = k[k.gwas == "MVP_Platelet_EUR"]
def split(h): c, p = h.split(":"); return c, int(p)
def cluster(leads, gap=5e5):
    d = pd.DataFrame([split(h) for h in sorted(set(leads))], columns=["chr", "pos"])
    d["chrn"] = d.chr.astype(int); d = d.sort_values(["chrn", "pos"])
    cid, c, prev = [], -1, (None, None)
    for ch, p in zip(d.chr, d.pos):
        if ch != prev[0] or p - prev[1] > gap: c += 1
        cid.append(c); prev = (ch, p)
    d["locus"] = cid; d["hit1"] = d.chr + ":" + d.pos.astype(str); return d
def near(h, dic):
    c, p = split(h); a = dic.get(c); return a is not None and np.any(np.abs(a - p) <= 5e5)
dl = cluster(L.hit1)
pv = pd.DataFrame([split(h) for h in sorted(set(P.hit1))], columns=["chr", "pos"]).groupby("chr").pos.apply(np.array).to_dict()
lv = dl.groupby("chr").pos.apply(np.array).to_dict()
dl["has_plt"] = [near(h, pv) for h in dl.hit1]
lp = dl.groupby("locus").agg(chr=("chr", "first"), has_plt=("has_plt", "any"))
conf = lp[~lp.chr.isin(["1", "19"])]
print("liver loci all chr:", len(lp), "| with platelet lead within 500 kb:", int(lp.has_plt.sum()), "| liver-only:", int((~lp.has_plt).sum()))
print("liver confirmatory (no chr1/19):", len(conf), "| shared:", int(conf.has_plt.sum()), "| liver-only:", int((~conf.has_plt).sum()))
po = [h for h in sorted(set(P.hit1)) if not near(h, lv)]
dp = cluster(po)
print("platelet leads within 500 kb of a liver lead:", len(set(P.hit1)) - len(po), "| platelet-only loci:", dp.locus.nunique(), "| excl chr1/19:", dp[~dp.chr.isin(["1", "19"])].locus.nunique())
m = cluster(L.loc[L.gwas == "MVP_ALT_EUR", "hit1"])
m["has_plt"] = [near(h, pv) for h in m.hit1]
mm = m.groupby("locus").agg(chr=("chr", "first"), has_plt=("has_plt", "any"))
mc = mm[~mm.chr.isin(["1", "19"])]
print("MVP_ALT_EUR confirmatory loci:", len(mc), "| without platelet lead within 500 kb:", int((~mc.has_plt).sum()))
print("chr5 blind-control row:", L[L.hit1 == "5:166834153"][["gene", "gwas", "hit1", "idx1", "idx2", "PP.H4"]].to_string(index=False, header=False))
