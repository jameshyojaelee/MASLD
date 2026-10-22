"""Independent re-derivation of Fig 2D coding plus splice fractions (pooled,
per ancestry, subclass, concentration). Reads release files only.
Seed 20260923."""
import zlib
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = ROOT + "/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2d/out"
R = ROOT + "/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816"
D, E = "tier1_direct_masld_pdff", "tier2_liver_enzyme"
GAP = 1_000_000
NB = 10000

d = pd.read_csv(R + "/multiplicity/corrected_locus_coordinates.tsv", sep="\t")
assert len(d) == 6415 and d.pip_sum_gate_passed.all()
d["num"] = d.protein_altering_pip_mass + d.canonical_splice_pip_mass
d["den"] = d[["protein_altering_pip_mass", "canonical_splice_pip_mass",
              "synonymous_or_utr_pip_mass", "other_noncoding_pip_mass"]].sum(axis=1)
d["chrom"] = d.chrom.astype(str)


def clusters(frame):
    """1-Mb single linkage on sentinel positions, per chromosome."""
    f = frame.sort_values(["chrom", "pos"])
    ids = pd.Series(index=f.index, dtype=object)
    for c, b in f.groupby("chrom", sort=False):
        new = np.concatenate([[True], np.diff(b.pos.to_numpy()) > GAP])
        ids.loc[b.index] = [f"{c}:{i}" for i in np.cumsum(new)]
    return ids.reindex(frame.index)


def frac(frame, cl):
    g = frame.assign(cl=cl).groupby(["trait_scope", "cl"])[["num", "den"]].mean()
    out = {}
    for s in (D, E):
        if s in g.index.get_level_values(0):
            gs = g.loc[s]
            out[s] = (gs.num.sum() / gs.den.sum(), len(gs))
        else:
            out[s] = (np.nan, 0)
    return out, g


log = []
def p(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    log.append(s)

# pooled
d["pcl"] = clusters(d)
o, gpool = frac(d, d.pcl)
p("POOLED f_direct", o[D], "f_enzyme", o[E], "diff", o[D][0] - o[E][0], "ratio", o[D][0] / o[E][0])
present = d.groupby("pcl").trait_scope.agg(lambda x: set(x))
p("physical", len(present), "shared", sum(len(x) == 2 for x in present))

# pooled bootstrap (my own design: resample physical clusters)
def boot_pooled(frame, cl, fit_id, nb=NB):
    g = frame.assign(cl=cl).groupby(["cl", "trait_scope"])[["num", "den"]].mean().unstack(fill_value=0.0)
    na, da = g[("num", D)].to_numpy(), g[("den", D)].to_numpy()
    ne, de = g[("num", E)].to_numpy(), g[("den", E)].to_numpy()
    rng = np.random.default_rng([20260923, zlib.crc32(fit_id.encode())])
    W = np.zeros((nb, len(g)))
    idx = rng.integers(0, len(g), size=(nb, len(g)))
    for i in range(nb):
        W[i] = np.bincount(idx[i], minlength=len(g))
    with np.errstate(divide="ignore", invalid="ignore"):
        fa = W @ na / (W @ da)
        fb = W @ ne / (W @ de)
    ok = np.isfinite(fa) & np.isfinite(fb)
    diff = (fa - fb)[ok]
    return diff, fa[ok], fb[ok]

diff, fa, fb = boot_pooled(d, d.pcl, "check:pooled")
p("POOLED my boot diff 95%", np.percentile(diff, [2.5, 97.5]), "ratio 95%",
  np.nanpercentile(np.where(fb > 0, fa / fb, np.nan), [2.5, 97.5]))

# per-stratum point estimates (clusters formed within stratum)
rows = []
strata = {"EUR": ["EUR"], "AMR": ["AMR"], "nonEUR": ["AFR", "AMR", "EAS", "SAS"],
          "AFR": ["AFR"], "EAS": ["EAS"]}
scl = {}
for name, anc in strata.items():
    sub = d[d.ancestry.isin(anc)]
    cl = clusters(sub)
    scl[name] = cl
    o, _ = frac(sub, cl)
    pr = sub.assign(cl=cl).groupby("cl").trait_scope.agg(lambda x: set(x))
    p(f"STRATUM {name}: f_dir {o[D][0]:.6f} ({o[D][1]} cl) f_enz {o[E][0]:.6f} ({o[E][1]} cl)"
      f" diff {o[D][0]-o[E][0]:.6f} ratio {o[D][0]/o[E][0] if o[E][0] else np.nan:.6f} phys {len(pr)} shared {sum(len(x)==2 for x in pr)}")
    rows.append(dict(stratum=name, f_direct=o[D][0], f_enzyme=o[E][0], diff=o[D][0] - o[E][0],
                     n_cl_direct=o[D][1], n_cl_enzyme=o[E][1]))

# shared bootstrap over the pooled clusters holding EUR or AMR rows (prespec design)
ea = d[d.ancestry.isin(["EUR", "AMR"])]
pool_ids = sorted(ea.pcl.unique(), key=lambda s: (int(s.split(":")[0]) if s.split(":")[0].isdigit() else 99, int(s.split(":")[1])))
p("pooled clusters holding EUR/AMR:", len(pool_ids))
pos = {c: i for i, c in enumerate(pool_ids)}
# per pooled cluster: summed within-stratum cluster-mean num/den for (stratum, class)
arms = {}
for s in ("EUR", "AMR"):
    sub = d[d.ancestry == s].copy()
    sub["scl"] = scl[s].loc[sub.index]
    g = sub.groupby(["trait_scope", "scl"]).agg(num=("num", "mean"), den=("den", "mean"), pcl=("pcl", "first"), npcl=("pcl", "nunique"))
    assert (g.npcl == 1).all(), "within-stratum cluster spans >1 pooled cluster"
    for x in (D, E):
        gx = g.loc[x]
        vn = np.zeros(len(pool_ids)); vd = np.zeros(len(pool_ids))
        for _, r in gx.iterrows():
            vn[pos[r.pcl]] += r.num; vd[pos[r.pcl]] += r.den
        arms[(s, x)] = (vn, vd)
rng = np.random.default_rng([20260923, zlib.crc32(b"H58:C:boot")])
idx = rng.integers(0, len(pool_ids), size=(NB, len(pool_ids)))
W = np.stack([np.bincount(r, minlength=len(pool_ids)) for r in idx]).astype(float)
F = {}
for k, (vn, vd) in arms.items():
    with np.errstate(divide="ignore", invalid="ignore"):
        F[k] = W @ vn / (W @ vd)
DE = F[("EUR", D)] - F[("EUR", E)]
DA = F[("AMR", D)] - F[("AMR", E)]
wE = 16 / 26
C1 = wE * DE + (1 - wE) * DA
for nm, v in (("C1 EUR+AMR", C1), ("C2 EUR", DE), ("C3 AMR", DA)):
    ok = np.isfinite(v)
    p(f"{nm}: finite {ok.sum()} 95% {np.percentile(v[ok],[2.5,97.5])} 98.333% {np.percentile(v[ok],[100/120,100-100/120])} sd {v[ok].std(ddof=1):.5f}")
e = {r["stratum"]: r for r in rows}
p("C1 point", wE * e["EUR"]["diff"] + (1 - wE) * e["AMR"]["diff"],
  "fD", wE * e["EUR"]["f_direct"] + (1 - wE) * e["AMR"]["f_direct"],
  "fE", wE * e["EUR"]["f_enzyme"] + (1 - wE) * e["AMR"]["f_enzyme"])
# check weights: direct den share
for s in ("EUR", "AMR"):
    p("direct den sum", s, arms[(s, D)][1].sum())

# non-EUR pooled bootstrap (own design)
sub = d[d.ancestry != "EUR"]
diff, _, _ = boot_pooled(sub, scl["nonEUR"], "check:nonEUR")
p("nonEUR my boot diff 95%", np.percentile(diff, [2.5, 97.5]))

# subclass split
d["sub"] = np.where(d.trait_scope == E, "ENZ", np.where(d.trait == "PDFF", "PDFF", "NAFLD_NASH"))
for sc in ("PDFF", "NAFLD_NASH"):
    sub = d[d["sub"].isin([sc, "ENZ"])]
    o1, _ = frac(sub, sub.pcl)  # pooled cluster ids
    cl2 = clusters(sub)
    o2, _ = frac(sub, cl2)       # re-clustered on subset
    diff, fa, fb = boot_pooled(sub, sub.pcl, f"check:sub:{sc}")
    p(f"SUBCLASS {sc}: n_cs {int((sub['sub']==sc).sum())} pooledIDs f {o1[D][0]:.6f} ({o1[D][1]} cl) vs {o1[E][0]:.6f} ({o1[E][1]} cl) diff {o1[D][0]-o1[E][0]:.6f};"
      f" reclustered f {o2[D][0]:.6f} ({o2[D][1]}) vs {o2[E][0]:.6f} ({o2[E][1]}); my boot diff95 {np.percentile(diff,[2.5,97.5])}")

# concentration of direct sum
gd = d[d.trait_scope == D].groupby("pcl").agg(m=("num", "mean"), dm=("den", "mean"), n=("num", "size"),
                                              pos=("pos", "min"), chrom=("chrom", "first"),
                                              studies=("study", lambda x: ";".join(sorted(set(x)))),
                                              has_eur=("ancestry", lambda x: "EUR" in set(x)))
gd = gd.sort_values("m", ascending=False)
gd["share"] = gd.m / gd.m.sum()
gd["cum"] = gd.share.cumsum()
p(gd[["chrom", "pos", "n", "m", "share", "cum", "has_eur"]].to_string())
p("direct clusters with a EUR set:", int(gd.has_eur.sum()), "of", len(gd))
p("top5 cum share", gd.cum.iloc[4])
# PDFF share of the direct sum: each cluster mean splits by member CS
dd = d[d.trait_scope == D].copy()
dd["ncl"] = dd.groupby("pcl").num.transform("size")
dd["contrib"] = dd.num / dd.ncl
p("PDFF share of summed direct cluster means", dd[dd.trait == "PDFF"].contrib.sum() / dd.contrib.sum())
# removal of top-k direct clusters
enz = d[d.trait_scope == E].groupby("pcl")[["num", "den"]].mean()
fE_all = enz.num.sum() / enz.den.sum()
for k in range(1, 6):
    top = gd.index[:k]
    rest = gd.drop(top)
    fD = rest.m.sum() / rest.dm.sum()
    enz2 = enz.drop([c for c in top if c in enz.index])
    fE2 = enz2.num.sum() / enz2.den.sum()
    p(f"drop top{k}: direct-rows-only ratio {fD/fE_all:.4f}; whole-physical-cluster ratio {fD/fE2:.4f}")

# EUR enzyme study counts and non-EUR enzyme sets
p("enzyme non-EUR sets", int(((d.trait_scope == E) & (d.ancestry != "EUR")).sum()))
pd.DataFrame(rows).to_csv(OUT + "/coding_strata_check.tsv", sep="\t", index=False)
gd.to_csv(OUT + "/direct_concentration_check.tsv", sep="\t")
open(OUT + "/02_coding.log", "w").write("\n".join(log) + "\n")
