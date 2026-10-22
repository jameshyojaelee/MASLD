"""Feasibility gate for L09-01 (equivalence negatives in bulk). Read-only on repo files.
Independent re-implementation; does not import the lens script."""
import sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
COH = f"{ROOT}/figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv"
POOL = f"{ROOT}/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
LOCO = f"{ROOT}/RNA-seq/results/manuscript_release/candidates/bulk-loco-{{c}}-20261009T141221Z/deg_results.csv"
STAGE = f"{ROOT}/figures/candidates/fig3e-s3-774-nas0ref-20261012T003338Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv"
NEG = f"{ROOT}/data/published_degs/PRJNA512027/replicated_16_genes.csv"
HK = sys.argv[1] if len(sys.argv) > 1 else None
COHORTS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"]
rng = np.random.default_rng(20260923)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def re_meta(Y, V, trunc=True):
    """Y, V: genes x k arrays. DerSimonian-Laird tau2, HK SE (truncated at 1 if trunc)."""
    k = Y.shape[1]
    W = 1 / V
    fe = (W * Y).sum(1) / W.sum(1)
    Q = (W * (Y - fe[:, None]) ** 2).sum(1)
    c = W.sum(1) - (W ** 2).sum(1) / W.sum(1)
    tau2 = np.maximum(0, (Q - (k - 1)) / c)
    Ws = 1 / (V + tau2[:, None])
    mu = (Ws * Y).sum(1) / Ws.sum(1)
    qhk = (Ws * (Y - mu[:, None]) ** 2).sum(1) / (k - 1)
    if trunc:
        qhk = np.maximum(qhk, 1.0)
    se_nom = np.sqrt(1 / Ws.sum(1))
    se = np.sqrt(qhk) * se_nom
    I2 = np.where(Q > 0, np.maximum(0, (Q - (k - 1)) / Q), 0)
    return mu, se, se_nom, tau2, Q, I2, fe, np.sqrt(1 / W.sum(1))


def tost_p(mu, se, m, df):
    if df is None:
        return np.maximum(stats.norm.sf((mu + m) / se), stats.norm.cdf((mu - m) / se))
    return np.maximum(stats.t.sf((mu + m) / se, df), stats.t.cdf((mu - m) / se, df))


d = pd.read_csv(COH, sep="\t")
print("## A. cohort file")
print("rows", len(d), "genes", d.gene_id_versioned.nunique(), "datasets", sorted(d.dataset.unique()))
print("contrast", d.contrast.unique().tolist())
print(d.groupby("dataset")[["n_control", "n_disease", "n_participants"]].first().to_string())
print("covariates by dataset", d.groupby("dataset").covariates.first().to_dict())
print("t==0 rows", int((d.t == 0).sum()), "NaN logFC", int(d.logFC.isna().sum()), "NaN t", int(d.t.isna().sum()))
d["SE"] = d.logFC / d.t
print("SE nonfinite", int((~np.isfinite(d.SE)).sum()), "SE<=0", int((d.SE <= 0).sum()))
print("per-cohort SE median:", d.groupby("dataset").SE.median().round(4).to_dict())
Y = d.pivot(index="gene_id_versioned", columns="dataset", values="logFC")[COHORTS]
V = d.pivot(index="gene_id_versioned", columns="dataset", values="SE")[COHORTS] ** 2
print("genes with all 5 cohorts", int(Y.notna().all(1).sum()))
genes = Y.index
sym = d.drop_duplicates("gene_id_versioned").set_index("gene_id_versioned").gene_name.reindex(genes)

mu, se, se_nom, tau2, Q, I2, fe, se_fe = re_meta(Y.values, V.values, trunc=True)
_, se_un, *_ = re_meta(Y.values, V.values, trunc=False)
print("\n## B. reproduce pre-look and alternatives (BH over", len(genes), "genes)")
res = {}
for m in (0.25, 0.5):
    p_tr = tost_p(mu, se, m, 4)
    p_un = tost_p(mu, se_un, m, 4)
    p_fe = tost_p(fe, se_fe, m, None)
    res[m] = p_tr
    print(f"margin {m}: RE-DL truncHK t4 {int((bh(p_tr) < 0.05).sum())}; RE-DL untruncHK t4 {int((bh(p_un) < 0.05).sum())}; FE inverse-variance z {int((bh(p_fe) < 0.05).sum())}")
p05 = res[0.5]
q05 = bh(p05)
call = q05 < 0.05
R = int(call.sum())
alpha_eff = 0.05 * R / len(genes)
tcrit = stats.t.ppf(1 - alpha_eff, 4)
print("max p among calls", p05[call].max(), "BH-effective one-sided alpha 0.05*R/m =", round(alpha_eff, 5), "t4 crit", round(tcrit, 3), "-> SE_HK must be <", round(0.5 / tcrit, 4))
print("pooled RE SE_HK quantiles 5/25/50/75/95:", np.quantile(se, [.05, .25, .5, .75, .95]).round(4).tolist())
print("fraction genes truncation active (qHK<1):", round(float((se == se_nom).mean()), 3))
print("tau2>0 fraction", round(float((tau2 > 0).mean()), 3), "tau median among tau2>0", round(float(np.sqrt(np.median(tau2[tau2 > 0]))), 4))
print("I2 median among calls", round(float(np.median(I2[call])), 3), "; among all", round(float(np.median(I2)), 3))
print("|mu| among calls quantiles 50/95/max:", np.quantile(np.abs(mu[call]), [.5, .95, 1]).round(4).tolist())

# pooled shipped fit (condition 2) and canonical
P = pd.read_csv(POOL).set_index("gene")
P = P.reindex(genes)
print("\n## C. pooled shipped fit")
print("pooled genes matched", int(P.logFC.notna().sum()), "of", len(genes))
canon = ((P.padj < 0.05) & (P.logFC.abs() > 0.5)).values
print("canonical calls", int(canon.sum()))
ci90 = (P.logFC - 1.645 * P.SE > -0.5) & (P.logFC + 1.645 * P.SE < 0.5)
c2 = ci90.values
print("cond1 (RE TOST BH)", R, "; cond1&2 (pooled 90% CI in +-0.5)", int((call & c2).sum()))

# LOCO condition (3) as stated: TOST nominal 0.05 in each LOCO pooled fit (z approx, df ~ hundreds)
print("\n## D. leave-one-cohort-out")
loco_ok = np.ones(len(genes), bool)
for c in COHORTS:
    L = pd.read_csv(LOCO.format(c=c)).set_index("gene").reindex(genes)
    miss = L.logFC.isna().values
    pl = tost_p(L.logFC.values, L.SE.values, 0.5, None)
    ok = (pl < 0.05) & ~miss
    print(c, "LOCO pooled fit: genes missing", int(miss.sum()), "; calls missing from LOCO", int((call & miss).sum()), "; calls failing nominal TOST", int((call & ~miss & ~ok).sum()))
    loco_ok &= ok
print("cond1&2&3 (stated rule)", int((call & c2 & loco_ok).sum()))
# LOCO with the same RE estimator (k=4, t3)
re4_ok = np.ones(len(genes), bool)
for i, c in enumerate(COHORTS):
    keep = [j for j in range(5) if j != i]
    m4, s4, *_ = re_meta(Y.values[:, keep], V.values[:, keep], trunc=True)
    p4 = tost_p(m4, s4, 0.5, 3)
    ok = p4 < 0.05
    print(c, "removed: RE k=4 t3 nominal TOST passes among calls", int((call & ok).sum()), "of", R)
    re4_ok &= ok
print("cond1&2 & RE-LOCO(k=4,t3) in all five", int((call & c2 & re4_ok).sum()))
final = call & c2 & loco_ok & ~canon
print("final (stated rule, canonical -> indeterminate)", int(final.sum()), "; calls that are canonical", int((call & canon).sum()))

# stage conflict
print("\n## E. conflict with stage contrasts (pooled 774 stage fit)")
S = pd.read_csv(STAGE, sep="\t", usecols=["gene_id_versioned", "axis", "contrast", "logFC", "FDR"])
S = S[S.axis.isin(["fibrosis", "NAS"])]
for con in ["F3_vs_F0", "F4_vs_F0", "NAS5_8_vs_NAS0"]:
    s = S[S.contrast == con].set_index("gene_id_versioned").reindex(genes)
    hit = ((s.FDR < 0.05) & (s.logFC.abs() > 0.5)).values
    hit1 = ((s.FDR < 0.05) & (s.logFC.abs() > 1.0)).values
    print(con, "BH<0.05 & |logFC|>0.5 among final", int((final & hit).sum()), "of", int(final.sum()), "; with |logFC|>1:", int((final & hit1).sum()))

# negative control
print("\n## F. negative control (16 Moylan genes)")
neg = pd.read_csv(NEG).Gene_ID.tolist()
for g in neg:
    idx = np.where(sym.values == g)[0]
    if len(idx) == 0:
        print(g, "absent")
        continue
    i = idx[0]
    print(f"{g}\tRE mu {mu[i]:.3f}\tSE_HK {se[i]:.3f}\tTOST p {p05[i]:.3g}\tq {q05[i]:.3g}\tcall {bool(call[i])}\tfinal {bool(final[i])}\tpooled logFC {P.logFC.values[i]:.3f}\tcohort logFC {np.round(Y.values[i], 2).tolist()}")

# power / MDE by simulation, using each gene's estimated V and tau2 as truth
print("\n## G. power of TOST at margin 0.5 (t4 crit at BH-effective alpha), truth mu in {0, 0.25}")
sub = rng.choice(len(genes), 3000, replace=False)
B = 400
for true_mu in (0.0, 0.25, 0.4):
    pw = []
    for i in sub:
        v = V.values[i]
        t2 = tau2[i]
        yy = true_mu + rng.normal(0, np.sqrt(t2), (B, 5)) + rng.normal(0, np.sqrt(v), (B, 5))
        mm, ss, *_ = re_meta(yy, np.tile(v, (B, 1)), trunc=True)
        pw.append(np.mean((mm - tcrit * ss > -0.5) & (mm + tcrit * ss < 0.5)))
    pw = np.array(pw)
    print(f"true mu {true_mu}: median power {np.median(pw):.3f}; fraction genes with power>=0.8 {np.mean(pw >= 0.8):.3f}; expected calls ~{np.mean(pw)*len(genes):.0f}")
# analytic SE threshold for 80% power at mu=0 ignoring SE noise: m/(tcrit+z0.9)
print("SE_nominal needed for ~80% power at mu=0 (ignoring HK noise):", round(0.5 / (tcrit + stats.norm.ppf(0.9)), 4), "; genes with se_nom below it:", int((se_nom < 0.5 / (tcrit + stats.norm.ppf(0.9))).sum()))

out = pd.DataFrame({"gene": genes, "symbol": sym.values, "mu": mu, "se_hk": se, "se_nom": se_nom, "tau2": tau2, "I2": I2,
                    "p_tost05": p05, "q_tost05": q05, "cond2": c2, "loco_pooled_ok": loco_ok, "loco_re4_ok": re4_ok,
                    "canonical": canon, "final": final, "AveExpr": P.AveExpr.values})
out.to_csv("gate_per_gene.tsv.gz", sep="\t", index=False)

if HK:
    print("\n## H. positive control feasibility (housekeeping set)")
    h = pd.read_csv(HK, sep="\t")
    base = pd.Series(genes).str.split(".").str[0].values
    in_u = np.isin(base, h.ensembl_gene.dropna().unique())
    in_s = np.isin(sym.values, h.gene_symbol.unique())
    print("HK set unique symbols", h.gene_symbol.nunique(), "unique ensembl", h.ensembl_gene.nunique())
    print("in universe by ensembl", int(in_u.sum()), "by symbol", int(in_s.sum()))
    hk = in_u | in_s
    print("HK among canonical", int((hk & canon).sum()), "expected by chance", round(canon.sum() * hk.sum() / len(genes), 1))
    dec = pd.qcut(P.AveExpr.values, 10, labels=False)
    tab = pd.DataFrame({"dec": dec, "hk": hk, "final": final}).groupby("dec").agg(n=("hk", "size"), hk=("hk", "sum"), final=("final", "sum"), hk_final=("final", lambda x: 0))
    tab["hk_final"] = pd.DataFrame({"dec": dec, "x": hk & final}).groupby("dec").x.sum()
    print(tab.to_string())
    # Mantel-Haenszel common OR across deciles with both HK and non-HK genes present
    num = den = 0.0
    for k_, g in pd.DataFrame({"dec": dec, "hk": hk, "f": final}).groupby("dec"):
        a = ((g.hk) & (g.f)).sum(); b = ((g.hk) & (~g.f)).sum(); c = ((~g.hk) & (g.f)).sum(); dd = ((~g.hk) & (~g.f)).sum(); n = len(g)
        num += a * dd / n; den += b * c / n
    print("MH OR (descriptive, gate pre-look):", round(num / den, 3) if den > 0 else "undefined")
    # simpler non-equivalence comparator: 'not BH significant' (absence of call) as the naive negative
    notsig = (P.padj >= 0.05).values
    num = den = 0.0
    for k_, g in pd.DataFrame({"dec": dec, "hk": hk, "f": notsig}).groupby("dec"):
        a = ((g.hk) & (g.f)).sum(); b = ((g.hk) & (~g.f)).sum(); c = ((~g.hk) & (g.f)).sum(); dd = ((~g.hk) & (~g.f)).sum(); n = len(g)
        num += a * dd / n; den += b * c / n
    print("MH OR for the naive rule 'pooled padj>=0.05' (comparator):", round(num / den, 3) if den > 0 else "undefined")

print("\n## I. heterogeneity among final calls (claim scope: mean vs each cohort)")
F = d.pivot(index="gene_id_versioned", columns="dataset", values="FDR")[COHORTS].reindex(genes).values
Yv = Y.values
any_cohort = ((F < 0.05) & (np.abs(Yv) > 0.5)).any(1)
print("final calls with >=1 cohort BH<0.05 and |logFC|>0.5:", int((final & any_cohort).sum()), "of", int(final.sum()))
for j, c in enumerate(COHORTS):
    print(" ", c, int((final & (F[:, j] < 0.05) & (np.abs(Yv[:, j]) > 0.5)).sum()))
# 95% prediction interval (Higgins 2009 form, t_{k-2}=t3)
pi_half = stats.t.ppf(0.975, 3) * np.sqrt(tau2 + se ** 2)
pi_in = (mu - pi_half > -0.5) & (mu + pi_half < 0.5)
print("final calls whose 95% prediction interval lies inside +-0.5:", int((final & pi_in).sum()), "of", int(final.sum()))
print("final calls with I2>=0.5:", int((final & (I2 >= 0.5)).sum()))
