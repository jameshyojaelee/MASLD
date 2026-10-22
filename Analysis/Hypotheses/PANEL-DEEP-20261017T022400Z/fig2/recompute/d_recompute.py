"""Independent recomputation of Figure 2D (main, ancestry and direct-subclass panels).

Upstream (analysis outputs, not panel source tables):
  GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/
    credible_set_pip_architecture.tsv   (one row per credible-set instance, class masses)
    credible_set_members_annotated.tsv  (member positions, normalized PIP, consequence)
Steps: PIP-sum gate -> restore sentinel coordinates from member positions ->
single-linkage 1-Mb clustering of sentinels pooled over both trait classes ->
cluster-balanced class fractions -> own cluster bootstrap (seed 20261017).
"""
import os
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN = f"{ROOT}/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out_D")
os.makedirs(OUT, exist_ok=True)

arch = pd.read_csv(f"{RUN}/credible_set_pip_architecture.tsv", sep="\t")
mem = pd.read_csv(f"{RUN}/credible_set_members_annotated.tsv", sep="\t", usecols=[
    "credible_set_uid", "chromosome", "position", "normalized_susie_pip", "consequence_category"])
print("instances", len(arch), "gated", int((arch.pip_sum_gate_passed.astype(str).str.lower() == "true").sum()))
arch["gate"] = arch.pip_sum_gate_passed.astype(str).str.lower() == "true"
a = arch[arch.gate].copy()

# member-level check of class masses
mm = mem.groupby(["credible_set_uid", "consequence_category"]).normalized_susie_pip.sum().unstack(fill_value=0)
chk = a.set_index("credible_set_uid")[["protein_altering_pip_mass", "canonical_splice_pip_mass",
                                       "synonymous_or_utr_pip_mass", "other_noncoding_pip_mass"]]
mm = mm.reindex(chk.index).fillna(0)
for c_arch, c_mem in [("protein_altering_pip_mass", "protein_altering"), ("canonical_splice_pip_mass", "canonical_splice"),
                      ("synonymous_or_utr_pip_mass", "synonymous_or_utr")]:
    if c_mem in mm:
        print("member-vs-arch max abs diff", c_mem, float((mm[c_mem] - chk[c_arch]).abs().max()))
print("member consequence categories", sorted(mem.consequence_category.dropna().unique()))

# restore sentinel coordinate: locus = "chr.pos" stored as decimal, trailing zeros may be lost
span = mem.groupby("credible_set_uid").position.agg(["min", "max"])
a = a.join(span, on="credible_set_uid")
a["chr"] = a.locus.astype(str).str.split(".").str[0]
a["enc"] = a.locus.astype(str).str.split(".").str[1]
def restore(row):
    base = int(row.enc)
    best, bestd = base, None
    for k in range(0, 4):
        p = base * 10 ** k
        d = 0 if row["min"] <= p <= row["max"] else min(abs(p - row["min"]), abs(p - row["max"]))
        if bestd is None or d < bestd:
            best, bestd = p, d
    return pd.Series({"pos": best, "dist_to_members": bestd, "k": int(np.log10(best / base)) if base else 0})
r = a.apply(restore, axis=1)
a = pd.concat([a, r], axis=1)
print("restored (k>0):", int((a.k > 0).sum()), "max distance sentinel-to-member span:", int(a.dist_to_members.max()))

a["cls"] = np.where(a.tier == 1, "direct", "enzyme")
a["pa_spl"] = a.protein_altering_pip_mass + a.canonical_splice_pip_mass

def cluster(df, w=1_000_000):
    lab = pd.Series(index=df.index, dtype=object)
    for c, g in df.groupby("chr"):
        g = g.sort_values("pos")
        cid, last = 0, None
        for i, p in zip(g.index, g.pos):
            if last is not None and p - last > w:
                cid += 1
            lab[i] = f"{c}_{cid}"
            last = p
    return lab

a["cluster"] = cluster(a)
cl = a.groupby("cluster").cls.agg(lambda s: set(s))
n_dir = sum("direct" in s for s in cl); n_enz = sum("enzyme" in s for s in cl)
n_sh = sum(len(s) == 2 for s in cl)
print(f"clusters total {len(cl)} direct {n_dir} enzyme {n_enz} shared {n_sh}")
print("instances by class", a.cls.value_counts().to_dict())

MASS = ["protein_altering_pip_mass", "canonical_splice_pip_mass", "synonymous_or_utr_pip_mass", "other_noncoding_pip_mass"]
def frac(df, col):
    """cluster-balanced fraction: mean over clusters of the mean CS mass in that cluster."""
    return df.groupby("cluster")[col].mean().mean()

rows = []
for c in ["direct", "enzyme"]:
    d = a[a.cls == c]
    rows.append({"class": c, "n_sets": len(d), "n_clusters": d.cluster.nunique(),
                 **{m: 100 * frac(d, m) for m in MASS}, "pa_spl": 100 * frac(d, "pa_spl"),
                 "cs_weighted_pa_spl": 100 * d.pa_spl.mean()})
main = pd.DataFrame(rows)
print(main.round(3).to_string())
main.to_csv(f"{OUT}/D_main_fractions.tsv", sep="\t", index=False)

# ancestry and subclass strata
def stratum(df):
    out = {}
    for c in ["direct", "enzyme"]:
        d = df[df.cls == c]
        out[c] = (100 * frac(d, "pa_spl") if len(d) else np.nan, len(d), d.cluster.nunique())
    return out
res = []
for name, mask in [("pooled", a.index == a.index), ("EUR", a.ancestry == "EUR"), ("AMR", a.ancestry == "AMR"),
                   ("nonEUR", a.ancestry != "EUR")]:
    s = stratum(a[mask])
    res.append({"stratum": name, "direct_pct": s["direct"][0], "direct_sets": s["direct"][1], "direct_clusters": s["direct"][2],
                "enzyme_pct": s["enzyme"][0], "enzyme_sets": s["enzyme"][1], "enzyme_clusters": s["enzyme"][2],
                "diff_pp": s["direct"][0] - s["enzyme"][0]})
eur = res[1]; amr = res[2]
res.append({"stratum": "EUR+AMR fixed 16/26,10/26",
            "direct_pct": 16 / 26 * eur["direct_pct"] + 10 / 26 * amr["direct_pct"],
            "enzyme_pct": 16 / 26 * eur["enzyme_pct"] + 10 / 26 * amr["enzyme_pct"]})
res[-1]["diff_pp"] = res[-1]["direct_pct"] - res[-1]["enzyme_pct"]
enz_all = a[a.cls == "enzyme"]
for name, traits in [("PDFF", ["PDFF"]), ("NAFLD+NASH", ["NAFLD", "NASH"])]:
    d = a[(a.cls == "direct") & a.trait.isin(traits)]
    res.append({"stratum": f"subclass {name}", "direct_pct": 100 * frac(d, "pa_spl"), "direct_sets": len(d),
                "direct_clusters": d.cluster.nunique(), "enzyme_pct": 100 * frac(enz_all, "pa_spl"),
                "enzyme_sets": len(enz_all), "enzyme_clusters": enz_all.cluster.nunique()})
    res[-1]["diff_pp"] = res[-1]["direct_pct"] - res[-1]["enzyme_pct"]
    print(name, "studies", sorted(d.study.unique()), "ancestries", d.ancestry.value_counts().to_dict())
# EUR enzyme comparator for PDFF (reasonable alternative)
d = a[(a.cls == "direct") & (a.trait == "PDFF")]
e = a[(a.cls == "enzyme") & (a.ancestry == "EUR")]
res.append({"stratum": "alt: PDFF vs EUR enzyme", "direct_pct": 100 * frac(d, "pa_spl"), "enzyme_pct": 100 * frac(e, "pa_spl"),
            "enzyme_sets": len(e), "enzyme_clusters": e.cluster.nunique()})
res[-1]["diff_pp"] = res[-1]["direct_pct"] - res[-1]["enzyme_pct"]
strata = pd.DataFrame(res)
print(strata.round(3).to_string())
strata.to_csv(f"{OUT}/D_strata.tsv", sep="\t", index=False)

# Own physical-cluster bootstrap (not the producer's stream): resample clusters with replacement;
# a sampled shared cluster carries both classes.
rng = np.random.default_rng(20261017)
def boot(df, B=10000):
    cm = df.groupby(["cluster", "cls"]).pa_spl.mean().unstack()
    clusters = cm.index.values
    vals = cm.values  # columns direct, enzyme
    cols = list(cm.columns)
    idx = rng.integers(0, len(clusters), size=(B, len(clusters)))
    out = np.empty(B)
    for b in range(B):
        v = vals[idx[b]]
        dmean = np.nanmean(v[:, cols.index("direct")]) if "direct" in cols else np.nan
        emean = np.nanmean(v[:, cols.index("enzyme")])
        out[b] = 100 * (dmean - emean)
    return np.nanpercentile(out, [2.5, 97.5]), np.nanpercentile(out, [100 / 120, 100 - 100 / 120]), int(np.isnan(out).sum())
bt = []
for name, df in [("pooled", a), ("EUR", a[a.ancestry == "EUR"]), ("AMR", a[a.ancestry == "AMR"]), ("nonEUR", a[a.ancestry != "EUR"]),
                 ("PDFF_vs_all_enzyme", a[((a.cls == "direct") & (a.trait == "PDFF")) | (a.cls == "enzyme")]),
                 ("NAFLDNASH_vs_all_enzyme", a[((a.cls == "direct") & a.trait.isin(["NAFLD", "NASH"])) | (a.cls == "enzyme")])]:
    ci95, ci98, nnan = boot(df)
    bt.append({"stratum": name, "lo95": ci95[0], "hi95": ci95[1], "lo98.333": ci98[0], "hi98.333": ci98[1], "nan_draws": nnan})
bt = pd.DataFrame(bt)
print(bt.round(2).to_string())
bt.to_csv(f"{OUT}/D_own_bootstrap.tsv", sep="\t", index=False)
a.to_csv(f"{OUT}/D_gated_instances_with_clusters.tsv.gz", sep="\t", index=False)
