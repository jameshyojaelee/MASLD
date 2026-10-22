"""Independent context fractions: adopted flags and position-repaired flags.
Credible-set-weighted (two definitions), pooled cluster-weighted differences,
and AMR differences with 99% intervals from a shared bootstrap over pooled
context clusters holding EUR or AMR rows. Seed 20260923."""
import zlib
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = ROOT + "/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2d/out"
D, E = "tier1_direct_masld_pdff", "tier2_liver_enzyme"
GAP = 1_000_000
NB = 10000
co = pd.read_csv(ROOT + "/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv", sep="\t")
m = pd.read_csv(ROOT + "/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv", sep="\t", low_memory=False)
m = m[m.credible_set_uid.isin(co.credible_set_uid)].copy()
s = pd.read_csv(OUT + "/sites_reannotated_check.tsv", sep="\t")
s = s.rename(columns={"chrom": "chromosome", "pos": "position"})
m = m.merge(s, on=["chromosome", "position"], how="left", validate="m:1")
print("gated member rows", len(m))

def tf(x):
    return x.astype(str).str.lower() == "true"

cons = m.consequence.astype(str)
utr = cons.str.contains(r"(?:^|,)[35]_prime_UTR_variant(?:$|,)", regex=True)
syn = cons.str.contains(r"(?:^|,)synonymous_variant(?:$|,)", regex=True)
m["nc"] = (m.consequence_category == "other_noncoding") | ((m.consequence_category == "synonymous_or_utr") & utr & ~syn)
pip = m.normalized_susie_pip

# unresolved share by ancestry (rows, all gated rows; and noncoding rows; and PIP)
m["unres_adopt"] = (m.promoter_context_status != "resolved") | (m.lineage_accessibility_status != "resolved")
for lab, sub in (("all gated rows", m), ("noncoding rows", m[m.nc])):
    t = sub.groupby("ancestry").agg(rows=("unres_adopt", "size"), unres=("unres_adopt", "sum"))
    t["share"] = t.unres / t.rows
    w = sub.assign(w=sub.normalized_susie_pip * sub.unres_adopt).groupby("ancestry").agg(pip=("normalized_susie_pip", "sum"), upip=("w", "sum"))
    t["pip_share"] = w.upip / w.pip
    print("UNRESOLVED", lab); print(t.to_string())
print("sites with 1:1 map", s.pos38.notna().sum(), "of", len(s))

# masses
def masses(flags):
    pr_res, ln_res, pr, ab, ac = flags
    out = pd.DataFrame({
        "credible_set_uid": m.credible_set_uid, "trait_scope": m.trait_scope, "ancestry": m.ancestry,
        "noncoding": pip * m.nc,
        "promoter": pip * m.nc * (pr_res & pr),
        "abc": pip * m.nc * ab,
        "accessible": pip * m.nc * (ln_res & ac),
        "none": pip * m.nc * (pr_res & ln_res & ~pr & ~ab & ~ac),
        "unresolved": pip * m.nc * (~pr_res | ~ln_res),
    })
    g = out.groupby(["credible_set_uid", "trait_scope", "ancestry"], as_index=False).sum()
    return g[g.noncoding > 0]

adopt = (m.promoter_context_status == "resolved", m.lineage_accessibility_status == "resolved",
         tf(m.promoter_proximal), tf(m.abc_enhancer_overlap), tf(m.lineage_accessible))
mapped = m.pos38.notna()
rep = (mapped, mapped, mapped & tf(m.promoter), tf(m.abc), mapped & tf(m.accessible))
CT = ["promoter", "abc", "accessible", "none", "unresolved"]
res = {}
for lab, fl in (("adopted", adopt), ("repaired", rep)):
    g = masses(fl)
    res[lab] = g
    print(f"== {lab}: context CS direct {int((g.trait_scope==D).sum())} enzyme {int((g.trait_scope==E).sum())}")
    for x in (D, E):
        gx = g[g.trait_scope == x]
        eq = {c: (gx[c] / gx.noncoding).mean() for c in CT}
        mw = {c: gx[c].sum() / gx.noncoding.sum() for c in CT}
        print(f"  {x[:5]} CS-equal-weight " + " ".join(f"{c}={eq[c]*100:.2f}%" for c in CT))
        print(f"  {x[:5]} mass-weighted    " + " ".join(f"{c}={mw[c]*100:.2f}%" for c in CT))

# flag agreement on rows resolved before repair
rr = m.promoter_context_status == "resolved"
for a, b in (("promoter_proximal", "promoter"), ("lineage_accessible", "accessible"), ("abc_enhancer_overlap", "abc")):
    print("agreement on adopted-resolved rows", a, (tf(m.loc[rr, a]) == tf(m.loc[rr, b])).mean(), "n", int(rr.sum()))

# cluster estimator
pos = co.set_index("credible_set_uid")[["chrom", "pos"]]
def clusters(frame):
    f = frame.join(pos, on="credible_set_uid").sort_values(["chrom", "pos"])
    ids = pd.Series(index=f.index, dtype=object)
    for c, b in f.groupby("chrom", sort=False):
        new = np.concatenate([[True], np.diff(b.pos.to_numpy()) > GAP])
        ids.loc[b.index] = [f"{c}:{i}" for i in np.cumsum(new)]
    return ids.reindex(frame.index)

def point(frame, cl, col):
    g = frame.assign(cl=cl).groupby(["trait_scope", "cl"])[[col, "noncoding"]].mean()
    fa = g.loc[D][col].sum() / g.loc[D].noncoding.sum()
    fb = g.loc[E][col].sum() / g.loc[E].noncoding.sum()
    return fa, fb, len(g.loc[D]), len(g.loc[E])

for lab in ("adopted", "repaired"):
    g = res[lab].copy()
    g["pcl"] = clusters(g)
    for col in CT:
        fa, fb, na, nb = point(g, g.pcl, col)
        # bootstrap physical clusters
        w = g.groupby(["pcl", "trait_scope"])[[col, "noncoding"]].mean().unstack(fill_value=0.0)
        rng = np.random.default_rng([20260923, zlib.crc32(f"check:ctx:{lab}:{col}".encode())])
        idx = rng.integers(0, len(w), size=(NB, len(w)))
        W = np.stack([np.bincount(r, minlength=len(w)) for r in idx]).astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            da = W @ w[(col, D)].to_numpy() / (W @ w[("noncoding", D)].to_numpy())
            db = W @ w[(col, E)].to_numpy() / (W @ w[("noncoding", E)].to_numpy())
        dd = (da - db)[np.isfinite(da) & np.isfinite(db)]
        print(f"POOLED {lab} {col}: {fa:.6f} vs {fb:.6f} diff {fa-fb:.6f} clusters {na}/{nb} my95 {np.percentile(dd,[2.5,97.5])}")

# AMR family X with shared bootstrap over pooled context clusters holding EUR/AMR rows
g = res["repaired"].copy()
g["pcl"] = clusters(g)
ea = g[g.ancestry.isin(["EUR", "AMR"])]
pids = sorted(ea.pcl.unique(), key=lambda z: (int(z.split(":")[0]) if z.split(":")[0].isdigit() else 99, int(z.split(":")[1])))
print("pooled context clusters with EUR/AMR rows", len(pids))
ix = {c: i for i, c in enumerate(pids)}
amr = g[g.ancestry == "AMR"].copy()
amr["scl"] = clusters(amr)
rng = np.random.default_rng([20260923, zlib.crc32(b"H58:X:boot")])
idx = rng.integers(0, len(pids), size=(NB, len(pids)))
W = np.stack([np.bincount(r, minlength=len(pids)) for r in idx]).astype(float)
for col in CT:
    fa, fb, na, nb = point(amr, amr.scl, col)
    vec = {}
    for x in (D, E):
        t = amr[amr.trait_scope == x].groupby("scl").agg(n=(col, "mean"), d=("noncoding", "mean"), p=("pcl", "first"), k=("pcl", "nunique"))
        assert (t.k == 1).all()
        vn = np.zeros(len(pids)); vd = np.zeros(len(pids))
        for _, r in t.iterrows():
            vn[ix[r.p]] += r.n; vd[ix[r.p]] += r.d
        vec[x] = (vn, vd)
    with np.errstate(divide="ignore", invalid="ignore"):
        da = W @ vec[D][0] / (W @ vec[D][1]); db = W @ vec[E][0] / (W @ vec[E][1])
    dd = (da - db)[np.isfinite(da) & np.isfinite(db)]
    print(f"AMR {col}: {fa:.6f} vs {fb:.6f} diff {fa-fb:.6f} clusters {na}/{nb} 95% {np.percentile(dd,[2.5,97.5])} 99% {np.percentile(dd,[0.5,99.5])}")
