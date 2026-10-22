# Feasibility gate for candidate C-01 (signed TWAS z vs disease-state DE t).
# Blinded: the observed primary statistic (TWAS z paired with its own gene's DE t) is never computed.
import glob, json
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr, norm

B = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
OT = B + "RNA-seq/results/causal_inference/otters_broadaway/"
REL = B + "RNA-seq/results/manuscript_release/candidates/"
SEED = 20261015
NDRAW = 2000
rng = np.random.default_rng(SEED)

# ---------- A. DE table ----------
de = pd.read_csv(REL + "bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
print("A DE rows", len(de), "cols", list(de.columns))
de["canon"] = (de.padj < 0.05) & (de.logFC.abs() > 0.5)
print("A canonical", int(de.canon.sum()), "up", int((de.canon & (de.logFC > 0)).sum()),
      "down", int((de.canon & (de.logFC < 0)).sum()), "NaN t", int(de.t.isna().sum()))
de = de[~de.gene.str.contains("PAR_Y")].copy()
de["sid"] = de.gene.str.replace(r"\.\d+$", "", regex=True)
print("A non-PAR_Y rows", len(de), "unique sid", de.sid.nunique())

# ---------- B. OTTERS P0.05 ----------
def load(t, model="P0.05"):
    fs = sorted(glob.glob(OT + f"{t}/chr*/{model}.txt"))
    tw = pd.concat([pd.read_csv(f, sep="\t") for f in fs])
    n0 = len(tw)
    tw = tw.dropna(subset=["FUSION_Z"])
    tw["mid"] = (tw.GeneStart + tw.GeneEnd) / 2
    return tw, n0, len(fs)

tw = {}
for t in ["ukbb_alt", "ghodsian", "ukbb_ast"]:
    d, n0, nf = load(t)
    print(f"B {t}: files {nf} rows {n0} nonNaN {len(d)} dupTarget {int(d.TargetID.duplicated().sum())} "
          f"mean z^2 {np.mean(d.FUSION_Z**2):.3f} |z|>3 {int((d.FUSION_Z.abs()>3).sum())} |z|>2 {int((d.FUSION_Z.abs()>2).sum())}")
    tw[t] = d.drop_duplicates("TargetID")

# orientation anchors (P0.05 model)
for g, nm in [("ENSG00000167701", "GPT"), ("ENSG00000120053", "GOT1"), ("ENSG00000125166", "GOT2")]:
    vals = []
    for t in ["ukbb_alt", "ukbb_ast", "ghodsian"]:
        s = tw[t].set_index("TargetID").FUSION_Z
        vals.append(f"{t} {s.get(g, np.nan):.2f}")
    print("B anchor", nm, g, "; ".join(vals))

# ---------- C. join and chromosome sizes ----------
J = {}
for t in ["ukbb_alt", "ghodsian"]:
    m = tw[t].merge(de[["sid", "t", "AveExpr", "canon", "logFC"]], left_on="TargetID", right_on="sid")
    m = m.sort_values(["CHROM", "mid"]).reset_index(drop=True)
    sizes = m.groupby("CHROM").size()
    print(f"C {t}: joined {len(m)} canonical {int(m.canon.sum())} (up {int((m.canon&(m.logFC>0)).sum())}, down {int((m.canon&(m.logFC<0)).sum())}) "
          f"|z|>3 {int((m.FUSION_Z.abs()>3).sum())} |z|>3&canon {int(((m.FUSION_Z.abs()>3)&m.canon).sum())} "
          f"chroms {len(sizes)} min chr n {sizes.min()} (chr {sizes.idxmin()}) max {sizes.max()}")
    print(f"C {t}: rho(z, AveExpr) {spearmanr(m.FUSION_Z, m.AveExpr)[0]:.4f}; rho(t, AveExpr) {spearmanr(m.t, m.AveExpr)[0]:.4f}")
    J[t] = m

# ---------- D. COLOC genes ----------
col = pd.read_csv(B + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
print("D coloc rows", len(col), "has cols", [c for c in ["ensembl", "coloc_best_susie_pp4"] if c in col.columns])
gen = set(col.loc[col.coloc_best_susie_pp4 > 0.5, "ensembl"].astype(str).str.replace(r"\.\d+$", "", regex=True))
print("D genes SuSiE PP.H4>0.5:", len(gen))
for t in J:
    print(f"D {t}: colocalized genes in join {int(J[t].TargetID.isin(gen).sum())}")

# ---------- E. LOCO tables ----------
for g in ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"]:
    f = glob.glob(REL + f"bulk-loco-{g}-20261009T141221Z/deg_results.csv")
    if not f:
        print("E LOCO", g, "MISSING"); continue
    l = pd.read_csv(f[0])
    l["sid"] = l.gene.astype(str).str.replace(r"\.\d+$", "", regex=True)
    n = J["ukbb_alt"].TargetID.isin(set(l.sid)).sum()
    print(f"E LOCO drop {g}: rows {len(l)} has t {'t' in l.columns} has AveExpr {'AveExpr' in l.columns} genes shared with ALT join {int(n)}")

# ---------- F. blinded null for the PRIMARY stratified statistic ----------
def zrank(x):
    r = rankdata(x); return (r - r.mean()) / r.std()

def strat_stat(xvals, y_r_by_dec, dec_idx):
    # weighted mean of per-decile Spearman; weights = n per decile
    tot = 0.0; ntot = 0
    for k, idx in dec_idx.items():
        xr = zrank(xvals[idx])
        tot += np.mean(xr * y_r_by_dec[k]) * len(idx); ntot += len(idx)
    return tot / ntot

def prep(m):
    dec = pd.qcut(m.AveExpr.rank(method="first"), 10, labels=False).values
    dec_idx = {k: np.where(dec == k)[0] for k in range(10)}
    y_r = {k: zrank(m.t.values[idx]) for k, idx in dec_idx.items()}
    bounds = m.groupby("CHROM").indices
    keys = sorted(bounds)
    return dec_idx, y_r, bounds, keys

def shifted(x, bounds, keys):
    out = x.copy()
    for k in keys:
        idx = bounds[k]; n = len(idx); mn = min(500, n // 4)
        out[idx] = np.roll(x[idx], int(rng.integers(mn, n - mn + 1)))
    return out

res = {}
for t, m in J.items():
    dec_idx, y_r, bounds, keys = prep(m)
    x = m.FUSION_Z.values
    null = np.array([strat_stat(shifted(x, bounds, keys), y_r, dec_idx) for _ in range(NDRAW)])
    sd = null.std(); q = np.percentile(null, [2.5, 97.5])
    mde = (norm.isf(0.025 / 2) + norm.isf(0.2)) * sd
    mde05 = (norm.isf(0.05 / 2) + norm.isf(0.2)) * sd
    print(f"F {t}: stratified null n={len(m)} mean {null.mean():.5f} SD {sd:.5f} 2.5/97.5 {q[0]:.4f}/{q[1]:.4f} "
          f"indepSD {1/np.sqrt(len(m)):.5f} MDE(alpha .025 2-sided, 80%) {mde:.4f} MDE(alpha .05) {mde05:.4f}")
    res[t] = dict(n=len(m), null_sd=float(sd), mde_holm1=float(mde), mde_05=float(mde05))
    # G. jackknife SE on 20 null-shifted pairings (blinded), 200 contiguous blocks in genomic order
    jk = []
    blocks = np.array_split(np.arange(len(m)), 200)
    for _ in range(20):
        xs = shifted(x, bounds, keys)
        dfx = pd.DataFrame({"x": xs, "t": m.t.values, "A": m.AveExpr.values})
        vals = []
        for b in blocks:
            keep = np.ones(len(m), bool); keep[b] = False
            d2 = dfx[keep]
            dec = pd.qcut(d2.A.rank(method="first"), 10, labels=False).values
            s = 0.0
            for k in range(10):
                ii = dec == k
                s += np.mean(zrank(d2.x.values[ii]) * zrank(d2.t.values[ii])) * ii.sum()
            vals.append(s / len(d2))
        vals = np.array(vals)
        jk.append(np.sqrt((len(vals) - 1) / len(vals) * np.sum((vals - vals.mean()) ** 2)))
    print(f"G {t}: jackknife SE over 20 null pairings median {np.median(jk):.5f} range {min(jk):.5f}-{max(jk):.5f}; 95% CI half-width ~ {1.96*np.median(jk):.4f}")
    res[t]["jk_se_median"] = float(np.median(jk))

# ---------- H. positive control with the EXACT P0.05 z ----------
c = pd.read_csv(OT + "concordance/otters_vs_spredixcan_concordance.csv")
pc = c[["ensembl_id", "spx_z", "otters_z"]].merge(tw["ghodsian"][["TargetID", "CHROM", "mid", "FUSION_Z"]],
                                                  left_on="ensembl_id", right_on="TargetID")
pc = pc[np.isfinite(pc.spx_z)].sort_values(["CHROM", "mid"]).reset_index(drop=True)
print(f"H concordance file otters_z equals P0.05 z in {np.mean(np.isclose(pc.otters_z, pc.FUSION_Z)):.3f} of {len(pc)} genes")
rho_best = spearmanr(pc.otters_z, pc.spx_z)[0]
rho_p05 = spearmanr(pc.FUSION_Z, pc.spx_z)[0]
bnd = pc.groupby("CHROM").indices; ks = sorted(bnd)
xr = zrank(pc.FUSION_Z.values); yr = zrank(pc.spx_z.values)
nullpc = np.array([np.mean(zrank(shifted(pc.FUSION_Z.values, bnd, ks)) * yr) for _ in range(NDRAW)])
p = (1 + np.sum(np.abs(nullpc) >= abs(np.mean(xr * yr)))) / (NDRAW + 1)
print(f"H OTTERS-vs-SPrediXcan ghodsian n {len(pc)}: rho(best_z) {rho_best:.4f} rho(P0.05 z) {rho_p05:.4f} "
      f"sign agreement P0.05 {np.mean(np.sign(pc.FUSION_Z)==np.sign(pc.spx_z)):.3f} null SD {nullpc.std():.4f} null max|.| {np.abs(nullpc).max():.4f} p {p:.5f}")

# cross-trait ALT vs NAFLD with P0.05 z
ct = tw["ukbb_alt"][["TargetID", "CHROM", "mid", "FUSION_Z"]].merge(
    tw["ghodsian"][["TargetID", "FUSION_Z"]].rename(columns={"FUSION_Z": "z2"}), on="TargetID")
print(f"H cross-trait ALT vs NAFLD P0.05 z: n {len(ct)} spearman {spearmanr(ct.FUSION_Z, ct.z2)[0]:.4f}")

# ---------- I. power under sparse alternatives (synthetic pairing, blinded) ----------
# Random pairing of the real z and real (t, AveExpr, canon) marginals; then impose sign agreement
# with probability q among a gene subset S. Reports the expected stratified statistic.
def sim(m, subset_rule, q, reps=200):
    out = []
    x = m.FUSION_Z.values
    for _ in range(reps):
        perm = rng.permutation(len(m))
        t = m.t.values[perm].copy(); A = m.AveExpr.values[perm]; can = m.canon.values[perm]
        S = subset_rule(x, can)
        agree = rng.random(S.sum()) < q
        tt = np.abs(t[S]) * np.where(agree, np.sign(x[S]), -np.sign(x[S]))
        t[S] = tt
        dec = pd.qcut(pd.Series(A).rank(method="first"), 10, labels=False).values
        s = 0.0
        for k in range(10):
            ii = dec == k
            s += np.mean(zrank(x[ii]) * zrank(t[ii])) * ii.sum()
        out.append(s / len(m))
    return np.mean(out), int(S.sum())

for t, m in J.items():
    for lab, rule in [("|z|>3", lambda x, c: np.abs(x) > 3), ("|z|>2", lambda x, c: np.abs(x) > 2),
                      ("|z|>3 & canonical", lambda x, c: (np.abs(x) > 3) & c)]:
        for q in [0.6, 0.69, 0.8, 1.0]:
            e, nS = sim(m, rule, q, reps=50)
            print(f"I {t}: S={lab} (|S|~{nS}) sign agreement {q}: expected stratified rho {e:.4f}  vs MDE {res[t]['mde_holm1']:.4f}")

json.dump(res, open(B + "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-C-01/gate_c01.json", "w"), indent=1)
