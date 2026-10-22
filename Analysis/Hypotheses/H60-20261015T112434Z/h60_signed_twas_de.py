# H60: signed agreement between liver TWAS z (OTTERS P0.05, UKBB ALT) and disease-vs-control DE t (774-unit refit).
# Implements PRESPEC_H60.md (sha256 420dbb48...). Phase A (diagnostics, controls, stop rules) is written and hashed
# before Phase B computes any observed z-with-t statistic.
# Usage: python h60_signed_twas_de.py <OUT>       (env H60_DRYRUN=1: t vectors permuted, small B; for debugging only)
import glob, hashlib, json, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import rankdata, skew, norm

RUN = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H60-20261015T112434Z"
PRESPEC_SHA = "420dbb486398504d2c9bcbdb59f4d4790784fa3edafa79c453942e1ae17bc8fb"
B_ = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
OT = B_ + "RNA-seq/results/causal_inference/otters_broadaway/"
REL = B_ + "RNA-seq/results/manuscript_release/candidates/"
PRIM = REL + "bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
SEED = 20260923
OUT = sys.argv[1]
DRY = os.environ.get("H60_DRYRUN", "0") == "1"
NSHIFT = 200 if DRY else 10000
NPC5 = 4 if DRY else 100
NPC5_SHIFT = 200 if DRY else 1000
NCAL = 10 if DRY else 200
NBOOT = 100 if DRY else 2000
NBLK = 200
T0 = time.time()


def log(*a):
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)


def sha_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


if sha_file(os.path.join(RUN, "PRESPEC_H60.md")) != PRESPEC_SHA:
    sys.exit("PRESPEC hash mismatch")
chk = pd.read_csv(os.path.join(OUT, "input_hashes_checked.tsv"), sep="\t")
if (chk.status == "MISMATCH").any():
    sys.exit("input hash mismatch recorded by h60_verify_inputs.py")
log("prespec and input hashes verified; DRYRUN =", DRY)

# ------------------------------------------------------------------ core statistics
def zrank(x, axis=-1):
    r = rankdata(x, axis=axis)
    m = r.mean(axis=axis, keepdims=True)
    s = r.std(axis=axis, keepdims=True)
    return (r - m) / s


def decile_groups(A, k=10):
    """A in gene order; pd.qcut of rank(method='first'). Returns list of index arrays."""
    if k == 1:
        return [np.arange(len(A))]
    d = pd.qcut(pd.Series(np.asarray(A)).rank(method="first"), k, labels=False).values
    return [np.where(d == j)[0] for j in range(k)]


def strat(x, y, groups):
    n = sum(len(g) for g in groups)
    return sum(np.mean(zrank(x[g]) * zrank(y[g])) * len(g) for g in groups) / n


def chrom_bounds(chrom):
    chrom = np.asarray(chrom)
    keys = list(dict.fromkeys(chrom))  # order of appearance (genome order)
    st = [np.where(chrom == k)[0] for k in keys]
    for s in st:
        assert np.all(np.diff(s) == 1), "chromosome not contiguous"
    return st


def shift_amounts(bounds, nd, rng):
    n = np.array([len(s) for s in bounds])
    m = np.minimum(500, n // 4)
    return rng.integers(m, n - m + 1, size=(nd, len(bounds)))


def shifted_index(bounds, amounts):
    """amounts (b x nchr) -> index matrix (b x N) so that x[idx] is the shifted vector."""
    b = amounts.shape[0]
    N = sum(len(s) for s in bounds)
    idx = np.empty((b, N), dtype=np.int64)
    for c, s in enumerate(bounds):
        n = len(s)
        j = np.arange(n)[None, :]
        idx[:, s] = s[0] + (j - amounts[:, [c]]) % n  # np.roll by amount
    return idx


def null_strat(x, y, groups, bounds, nd, rng, chunk=500):
    """Shift x within chromosomes; groups and y stay with the gene. Returns null array."""
    amts = shift_amounts(bounds, nd, rng)
    yr = [zrank(y[g]) for g in groups]
    N = len(x)
    out = np.empty(nd)
    for a in range(0, nd, chunk):
        X = x[shifted_index(bounds, amts[a:a + chunk])]
        tot = np.zeros(X.shape[0])
        for g, yy in zip(groups, yr):
            tot += (zrank(X[:, g], axis=1) @ yy)  # = mean(prod) * n_k
        out[a:a + chunk] = tot / N
    return out


def jackknife(x, y, A, k=10, nblk=NBLK, return_vals=False):
    N = len(x)
    blocks = np.array_split(np.arange(N), nblk)
    vals = []
    for b in blocks:
        keep = np.ones(N, bool); keep[b] = False
        vals.append(strat(x[keep], y[keep], decile_groups(A[keep], k)))
    vals = np.array(vals)
    se = np.sqrt((nblk - 1) / nblk * np.sum((vals - vals.mean()) ** 2))
    return (se, vals) if return_vals else se


def pval(null, obs):
    return (1 + np.sum(np.abs(null) >= abs(obs))) / (len(null) + 1)


def null_diag(null, n, label):
    d = dict(null=label, n=n, mean=null.mean(), sd=null.std(), q025=np.percentile(null, 2.5),
             q975=np.percentile(null, 97.5), max_abs=np.abs(null).max(), skewness=skew(null),
             inv_sqrt_n=1 / np.sqrt(n), ndraw=len(null))
    d["mean_offset_gt_0.002"] = abs(d["mean"]) > 0.002
    return d


def wilson(k, n, z=1.959964):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


# ------------------------------------------------------------------ data
def load_otters(t, model="P0.05"):
    fs = [OT + f"{t}/chr{c}/{model}.txt" for c in range(1, 23)]
    d = pd.concat([pd.read_csv(f, sep="\t") for f in fs]).dropna(subset=["FUSION_Z"])
    d["mid"] = (d.GeneStart + d.GeneEnd) / 2
    assert not d.TargetID.duplicated().any()
    return d


def load_de(path):
    d = pd.read_csv(path)
    d["sid"] = d.gene.astype(str).str.replace(r"\.\d+$", "", regex=True)
    assert not d.sid.duplicated().any()
    return d


de = load_de(PRIM)
de["canon"] = (de.padj < 0.05) & (de.logFC.abs() > 0.5)
r1 = pd.read_csv(os.path.join(OUT, "refit_R1_hepadj_toptable.tsv.gz"), sep="\t")
r1["sid"] = r1.gene.str.replace(r"\.\d+$", "", regex=True)
de = de.merge(r1[["sid", "t"]].rename(columns={"t": "t_adj"}), on="sid", how="left")
assert de.t_adj.notna().all()
loco = {}
for g in ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621", "NONE"]:
    loco[g] = load_de(REL + f"bulk-loco-{g}-20261009T141221Z/deg_results.csv")

if DRY:  # destroy the real pairing for debugging
    prng = np.random.default_rng(1)
    p = prng.permutation(len(de))
    de["t"] = de.t.values[p]; de["t_adj"] = de.t_adj.values[p]
    for g in loco:
        loco[g]["t"] = loco[g].t.values[prng.permutation(len(loco[g]))]
    log("DRYRUN: t vectors permuted")

tw = {t: load_otters(t) for t in ["ukbb_alt", "ghodsian", "ukbb_ast", "ukbb_ggt"]}
tw_p001 = load_otters("ukbb_alt", "P0.001")


def join(twd, ded, cols=("t", "t_adj", "AveExpr", "canon", "logFC")):
    cols = [c for c in cols if c in ded.columns]
    m = twd.merge(ded[["sid"] + cols], left_on="TargetID", right_on="sid")
    return m.sort_values(["CHROM", "mid", "TargetID"]).reset_index(drop=True)


J = {t: join(tw[t], de) for t in tw}
alt = J["ukbb_alt"]
N = len(alt)

# COLOC r2 genes and 1 Mb neighbourhood (count_sets.py recipe)
col = pd.read_csv(B_ + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
gen = set(col.loc[col.coloc_best_susie_pp4 > 0.5, "ensembl"].astype(str).str.replace(r"\.\d+$", "", regex=True))


def coloc_keep(m):
    ing = m.TargetID.isin(gen).values
    keep = np.ones(len(m), bool)
    for ch, g in m[ing].groupby("CHROM"):
        idx = np.where((m.CHROM == ch).values)[0]
        d = np.min(np.abs(m.loc[idx, "mid"].values[:, None] - g.mid.values[None, :]), axis=1)
        keep[idx[d <= 1e6]] = False
    return keep


alt_keep = coloc_keep(alt)

diag = []  # (section, item, value)


def D(sec, item, val):
    diag.append((sec, item, val))


# ------------------------------------------------------------------ Phase A: pre-observation diagnostics
log("Phase A")
for _, r in chk.iterrows():
    D("input_hash", r["item"], f"{r['sha256']} {r['status']}")
for t, m in J.items():
    sz = m.groupby("CHROM").size()
    D("join", f"{t}_n", len(m)); D("join", f"{t}_n_chrom", len(sz))
    D("join", f"{t}_min_chr_n", int(sz.min())); D("join", f"{t}_max_chr_n", int(sz.max()))
    D("join", f"{t}_per_chrom", ";".join(f"{k}:{v}" for k, v in sz.items()))
    D("join", f"{t}_canonical", int(m.canon.sum())); D("join", f"{t}_absz_gt3", int((m.FUSION_Z.abs() > 3).sum()))
    D("join", f"{t}_decile_sizes", ";".join(str(len(g)) for g in decile_groups(m.AveExpr.values)))
D("join", "ukbb_alt_coloc_in_join", int(alt.TargetID.isin(gen).sum()))
D("join", "ukbb_alt_coloc_excluded_n", int(alt_keep.sum()))
D("join", "ukbb_alt_coloc_excluded_min_chr_n", int(alt[alt_keep].groupby("CHROM").size().min()))
D("join", "ukbb_alt_coloc_complement_n", int((~alt_keep).sum()))
D("join", "ukbb_alt_absz3_and_canonical", int(((alt.FUSION_Z.abs() > 3) & alt.canon).sum()))
mhc = ((alt.CHROM == 6) & (alt.mid >= 28477797) & (alt.mid <= 33448354)).values
D("join", "ukbb_alt_MHC_n", int(mhc.sum()))


def sp(a, b):
    return float(np.mean(zrank(np.asarray(a, float)) * zrank(np.asarray(b, float))))


for t in ["ukbb_alt", "ghodsian"]:
    D("rho_aveexpr", f"rho(z,AveExpr)_{t}", sp(J[t].FUSION_Z, J[t].AveExpr))
D("rho_aveexpr", "rho(t_primary,AveExpr)_ALTjoin", sp(alt.t, alt.AveExpr))
D("rho_aveexpr", "rho(t_adj,AveExpr)_ALTjoin", sp(alt.t_adj, alt.AveExpr))
for g in ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"]:
    lj = join(tw["ukbb_alt"], loco[g], cols=("t", "AveExpr"))
    D("rho_aveexpr", f"rho(t_LOCO_{g},AveExpr)_ALTjoin", sp(lj.t, lj.AveExpr))

# composition refit diagnostics
rc = pd.read_csv(os.path.join(OUT, "refit_checks.tsv"), sep="\t").set_index("item")["value"]
for k, v in rc.items():
    D("composition_refit", k, v)
D("composition_refit", "spearman(t_adj,t_primary)_all_DE_genes", sp(de.t_adj, de.t))
r0_pass = str(rc["R0_pass"]).upper() == "TRUE"
ln = loco["NONE"].set_index("sid").t.reindex(de.sid)
D("loco_none", "max_abs_t_diff_vs_primary", float(np.nanmax(np.abs(ln.values - de.t.values))) if not DRY else "DRYRUN")

# --- primary / qualifier / COLOC-excluded nulls (no observed statistic)
z = alt.FUSION_Z.values
A = alt.AveExpr.values
t_p = alt.t.values
t_a = alt.t_adj.values
G = decile_groups(A)
bnd = chrom_bounds(alt.CHROM.values)
nulls = {}
nulls["primary"] = null_strat(z, t_p, G, bnd, NSHIFT, np.random.default_rng(SEED))
nulls["qualifier"] = null_strat(z, t_a, G, bnd, NSHIFT, np.random.default_rng(SEED))
ak = alt[alt_keep].reset_index(drop=True)
Gk = decile_groups(ak.AveExpr.values); bk = chrom_bounds(ak.CHROM.values)
nulls["coloc_excluded"] = null_strat(ak.FUSION_Z.values, ak.t.values, Gk, bk, NSHIFT, np.random.default_rng(SEED))
nulls["coloc_excluded_adj"] = null_strat(ak.FUSION_Z.values, ak.t_adj.values, Gk, bk, NSHIFT, np.random.default_rng(SEED))
ndiag = [null_diag(nulls["primary"], N, "primary"), null_diag(nulls["qualifier"], N, "qualifier"),
         null_diag(nulls["coloc_excluded"], len(ak), "coloc_excluded"),
         null_diag(nulls["coloc_excluded_adj"], len(ak), "coloc_excluded_adj")]
log("primary/qualifier/coloc nulls done")


# jackknife SE on 20 null pairings
def jk_null(m, nrep=20):
    x = m.FUSION_Z.values; y = m.t.values; Am = m.AveExpr.values
    b = chrom_bounds(m.CHROM.values)
    idx = shifted_index(b, shift_amounts(b, nrep, np.random.default_rng(SEED)))
    return np.array([jackknife(x[i], y, Am) for i in idx])


jkn_p = jk_null(alt); jkn_k = jk_null(ak)
D("jk_null", "primary_median_SE_20_null_pairings", float(np.median(jkn_p)))
D("jk_null", "primary_range", f"{jkn_p.min():.5f}-{jkn_p.max():.5f}")
D("jk_null", "coloc_excluded_median_SE_20_null_pairings", float(np.median(jkn_k)))
D("jk_null", "coloc_excluded_range", f"{jkn_k.min():.5f}-{jkn_k.max():.5f}")
log("null jackknife done")


# calibration table
def impose(rng, m, q, subset):
    x = m.FUSION_Z.values
    perm = rng.permutation(len(m))
    t = m.t.values[perm].copy(); Ap = m.AveExpr.values[perm]; can = m.canon.values[perm]
    S = (np.abs(x) > 3) if subset == "a" else ((np.abs(x) > 3) & can)
    agree = rng.random(S.sum()) < q
    t[S] = np.abs(t[S]) * np.where(agree, np.sign(x[S]), -np.sign(x[S]))
    return x, t, Ap, int(S.sum())


cal = []
rng = np.random.default_rng(SEED)
for subset in ["a", "b"]:
    for q in [0.55, 0.60, 0.65, 0.69, 0.80, 1.0]:
        vals, ns = [], []
        for _ in range(NCAL):
            x, t, Ap, nS = impose(rng, alt, q, subset)
            vals.append(strat(x, t, decile_groups(Ap))); ns.append(nS)
        cal.append(dict(subset="ALT |z|>3" if subset == "a" else "ALT |z|>3 and canonical", q=q,
                        expected_rho=np.mean(vals), sd_rho=np.std(vals), mean_nS=np.mean(ns), reps=NCAL))
cal = pd.DataFrame(cal)
log("calibration done")

# ------------------------------------------------------------------ controls
ctrl = []


def unstrat_test(x, y, chrom, rng):
    G1 = [np.arange(len(x))]
    b = chrom_bounds(chrom)
    obs = strat(x, y, G1)
    nl = null_strat(x, y, G1, b, NSHIFT, rng)
    return obs, pval(nl, obs), nl


# PC1
conc = pd.read_csv(OT + "concordance/otters_vs_spredixcan_concordance.csv")
conc = conc[conc.gwas == "ghodsian"]
pc = conc[["ensembl_id", "spx_z"]].merge(tw["ghodsian"][["TargetID", "CHROM", "mid", "FUSION_Z"]],
                                         left_on="ensembl_id", right_on="TargetID")
pc = pc[np.isfinite(pc.spx_z)].sort_values(["CHROM", "mid", "TargetID"]).reset_index(drop=True)
o, p, nl = unstrat_test(pc.FUSION_Z.values, pc.spx_z.values, pc.CHROM.values, np.random.default_rng(SEED))
sa = float(np.mean(np.sign(pc.FUSION_Z) == np.sign(pc.spx_z)))
ctrl.append(dict(control="PC1", stat=o, p=p, n=len(pc), extra=f"sign agreement {sa:.3f}; min chr n {pc.groupby('CHROM').size().min()}",
                 predicted="rho 0.4377, sign agr 0.665, p at floor", rule="rho>0 and p<=0.001", passed=(o > 0) and (p <= 0.001)))
ndiag.append(null_diag(nl, len(pc), "PC1"))
# PC2
pc2 = tw["ukbb_alt"][["TargetID", "CHROM", "mid", "FUSION_Z"]].merge(
    tw["ghodsian"][["TargetID", "FUSION_Z"]].rename(columns={"FUSION_Z": "z2"}), on="TargetID")
pc2 = pc2.sort_values(["CHROM", "mid", "TargetID"]).reset_index(drop=True)
o, p, nl = unstrat_test(pc2.FUSION_Z.values, pc2.z2.values, pc2.CHROM.values, np.random.default_rng(SEED))
ctrl.append(dict(control="PC2", stat=o, p=p, n=len(pc2), extra="-", predicted="rho 0.1402, p at floor",
                 rule="rho>0 and p<=0.001", passed=(o > 0) and (p <= 0.001)))
ndiag.append(null_diag(nl, len(pc2), "PC2"))
# PC3
zs = {t: tw[t].set_index("TargetID").FUSION_Z for t in tw}
gpt = float(zs["ukbb_alt"].get("ENSG00000167701", np.nan))
ctrl.append(dict(control="PC3", stat=gpt, p=np.nan, n=1,
                 extra=f"GOT2 AST z {zs['ukbb_ast'].get('ENSG00000125166', np.nan):.2f} (printed, not a rule); GOT1 AST z {zs['ukbb_ast'].get('ENSG00000120053', np.nan):.2f} (printed, not a rule)",
                 predicted="+8.97", rule="GPT ALT z > 3", passed=gpt > 3))
# PC4
ts = de.set_index("sid").t
c1, c2 = float(ts["ENSG00000108821"]), float(ts["ENSG00000164692"])
ctrl.append(dict(control="PC4", stat=min(c1, c2), p=np.nan, n=2, extra=f"COL1A1 t {c1:.3f}; COL1A2 t {c2:.3f}",
                 predicted="both > 0", rule="both t > 0", passed=(c1 > 0) and (c2 > 0)))
# NC1, NC2 (ALT and NAFLD joins)
for t in ["ukbb_alt", "ghodsian"]:
    m = J[t]
    for nc, y, pred in [("NC1", m.AveExpr.values, {"ukbb_alt": "0.0140 (p~0.10)", "ghodsian": "-0.0103 (p~0.22)"}),
                        ("NC2", np.log(m.n_snps.values), {"ukbb_alt": "-0.0053 (p~0.56)", "ghodsian": "-0.0029 (p~0.72)"})]:
        o, p, nl = unstrat_test(m.FUSION_Z.values, y, m.CHROM.values, np.random.default_rng(SEED))
        ctrl.append(dict(control=f"{nc}_{t}", stat=o, p=p, n=len(m), extra="-", predicted=pred[t], rule="p>=0.01", passed=p >= 0.01))
        ndiag.append(null_diag(nl, len(m), f"{nc}_{t}"))
log("PC1-4, NC done")

# PC5: pipeline recovery through the full direction rule
rng = np.random.default_rng(SEED)
rec = []
for r in range(NPC5):
    x, t, Ap, nS = impose(rng, alt, 0.60, "a")
    Gp = decile_groups(Ap)
    o = strat(x, t, Gp)
    nl = null_strat(x, t, Gp, bnd, NPC5_SHIFT, rng)
    p = pval(nl, o)
    se = jackknife(x, t, Ap)
    rec.append(dict(rep=r, rho=o, p=p, se=se, lo=o - 1.96 * se, hi=o + 1.96 * se,
                    recovered=(p < 0.05) and (o - 1.96 * se > 0) and (o > 0)))
rec = pd.DataFrame(rec)
recovery = rec.recovered.mean()
ctrl.append(dict(control="PC5", stat=recovery, p=np.nan, n=NPC5,
                 extra=f"mean synthetic rho {rec.rho.mean():.4f}; median p {rec.p.median():.4f}; median SE {rec.se.median():.4f}",
                 predicted=">=0.90 [Estimate]", rule="recovery >= 0.80", passed=recovery >= 0.80))
log("PC5 done", recovery)
ctrl = pd.DataFrame(ctrl)

# ------------------------------------------------------------------ stop rules
c = ctrl.set_index("control")
nd = {d["null"]: d for d in ndiag}
altsz = alt.groupby("CHROM").size()
stops = [
    dict(rule="S1", value=f"ALT join {N}; min chr n {altsz.min()}", triggered=(N < 10000) or (altsz.min() < 50), consequence="untestable"),
    dict(rule="S2", value=f"PC1 pass {c.loc['PC1','passed']}; PC2 pass {c.loc['PC2','passed']}",
         triggered=not (c.loc["PC1", "passed"] and c.loc["PC2", "passed"]), consequence="untestable"),
    dict(rule="S3", value=f"GPT ALT z {gpt:.3f}", triggered=not c.loc["PC3", "passed"], consequence="untestable"),
    dict(rule="S4", value=f"COL1A1 t {c1:.3f}; COL1A2 t {c2:.3f}", triggered=not c.loc["PC4", "passed"], consequence="untestable"),
    dict(rule="S5", value="; ".join(f"{k} p {c.loc[k,'p']:.4f}" for k in c.index if k.startswith("NC")),
         triggered=bool((c.loc[[k for k in c.index if k.startswith("NC")], "p"] < 0.01).any()), consequence="untestable"),
    dict(rule="S6", value=f"PC5 recovery {recovery:.3f}", triggered=recovery < 0.80, consequence="untestable"),
    dict(rule="S7", value=f"primary null SD {nd['primary']['sd']:.5f}; median null-pairing jackknife SE {np.median(jkn_p):.5f}",
         triggered=(nd["primary"]["sd"] > 0.0107) or (np.median(jkn_p) > 0.0153), consequence="reading 3 removed"),
    dict(rule="S8", value=f"R0 pass {r0_pass}; max|t diff| {rc['R0_max_abs_t_diff']}", triggered=not r0_pass,
         consequence="no direction may be stated"),
]
stops = pd.DataFrame(stops)

# ------------------------------------------------------------------ write Phase A and close it
pre_files = {}
def wr(df, name):
    p = os.path.join(OUT, name); df.to_csv(p, sep="\t", index=False); pre_files[name] = p

wr(pd.DataFrame(diag, columns=["section", "item", "value"]), "diagnostics.tsv")
wr(pd.DataFrame(ndiag), "diagnostics_nulls.tsv")
wr(cal, "diagnostics_calibration.tsv")
wr(ctrl, "controls.tsv")
wr(rec, "controls_PC5_replicates.tsv")
wr(stops, "stop_rules.tsv")
np.savez_compressed(os.path.join(OUT, "nulls_pre_observation.npz"), **nulls)
pre_files["nulls_pre_observation.npz"] = os.path.join(OUT, "nulls_pre_observation.npz")
with open(os.path.join(OUT, "diagnostics_complete.txt"), "w") as f:
    f.write(f"Phase A closed {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} before any observed z-with-t statistic. DRYRUN={DRY}\n")
    for k, p in pre_files.items():
        f.write(f"{sha_file(p)}  {k}\n")
log("Phase A closed")
print(ctrl.to_string()); print(stops.to_string())

untestable = stops[(stops.rule.isin(["S1", "S2", "S3", "S4", "S5", "S6"])) & stops.triggered]


def write_result(lines):
    with open(os.path.join(OUT, "RESULT_H60.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


hdr = ["# RESULT H60: signed agreement of ALT-TWAS z with disease-vs-control DE t", "",
       f"Run dir: {OUT}. Prespec sha256 {PRESPEC_SHA}. Seed {SEED}. Shift draws {NSHIFT}. DRYRUN={DRY}.",
       "Inputs: OTTERS P0.05 (Broadaway liver eQTL weights; run of 2026-03-24, not part of any adopted release); "
       "DE t from bulk-corrected-controls-20261008T162308Z (774 genetic-identity clusters, not verified donors); "
       "unit of the statistic is the gene; null keeps genomic order (within-chromosome circular shift).", ""]
if len(untestable) and not DRY:  # DRYRUN continues past stop rules for debugging only
    write_result(hdr + ["## Verdict: UNTESTABLE (stop rule)", ""] +
                 [f"- {r.rule}: {r.value}" for r in untestable.itertuples()] +
                 ["", "No primary statistic was computed."])
    log("STOPPED by rule", list(untestable.rule))
    sys.exit(0)

# ------------------------------------------------------------------ Phase B: observed statistics
log("Phase B")
S7 = bool(stops.set_index("rule").loc["S7", "triggered"])
S8 = bool(stops.set_index("rule").loc["S8", "triggered"])


def full(x, y, Am, null, label):
    G_ = decile_groups(Am)
    o = strat(x, y, G_)
    se, vals = jackknife(x, y, Am, return_vals=True)
    p = pval(null, o)
    d = dict(set=label, n=len(x), rho=o, se=se, lo=o - 1.96 * se, hi=o + 1.96 * se, p_shift=p,
             null_sd=null.std(), null_mean=null.mean(), se_ratio=se / null.std(),
             jk_sign_flip_blocks=int(np.sum(np.sign(vals) != np.sign(o))))
    if abs(null.mean()) > 0.002:
        d["p_tail"] = min(1.0, 2 * min((1 + np.sum(null >= o)) / (len(null) + 1), (1 + np.sum(null <= o)) / (len(null) + 1)))
    return d


prim = full(z, t_p, A, nulls["primary"], "primary ALT (t primary)")
qual = full(z, t_a, A, nulls["qualifier"], "qualifier ALT (t hepatocyte-adjusted)")
log("primary", prim); log("qualifier", qual)

excl0 = lambda d: (d["lo"] > 0) or (d["hi"] < 0)
inside = lambda d: (d["lo"] >= -0.03) and (d["hi"] <= 0.03)
if prim["p_shift"] < 0.05 and excl0(prim):
    reading = 1 if prim["rho"] > 0 else 2
elif prim["p_shift"] >= 0.05 and inside(prim) and not S7:
    reading = 3
else:
    reading = 4
qual_pass = (not S8) and np.sign(qual["rho"]) == np.sign(prim["rho"]) and qual["p_shift"] < 0.05 and excl0(qual)

fixed_seq = None
coloc = full(ak.FUSION_Z.values, ak.t.values, ak.AveExpr.values, nulls["coloc_excluded"], "D2 COLOC-excluded (t primary)")
if reading in (1, 2) and qual_pass:
    fixed_seq = bool(coloc["p_shift"] < 0.05 and excl0(coloc) and np.sign(coloc["rho"]) == np.sign(prim["rho"]))

qstar = None
if reading == 3:
    ca = cal[cal.subset == "ALT |z|>3"].sort_values("q")
    hit = ca[ca.expected_rho > prim["hi"]]
    qstar = float(hit.q.iloc[0]) if len(hit) else None

if reading in (1, 2):
    if qual_pass:
        verdict = ("risk-concordant signed agreement" if reading == 1 else "opposed direction") + " (survives hepatocyte-fraction adjustment)"
    elif S8:
        verdict = ("risk-concordant" if reading == 1 else "opposed") + " unadjusted direction, qualifier unavailable (S8); no direction may be stated"
    else:
        verdict = "unadjusted agreement that does not survive hepatocyte-fraction adjustment; no direction stated"
elif reading == 3:
    verdict = "bounded genome-wide rank agreement inside +/-0.03" + ("" if inside(qual) else " (without composition adjustment)")
else:
    verdict = "inconclusive at this resolution"
log("READING", reading, verdict)

# ------------------------------------------------------------------ descriptives
desc = [prim, qual, coloc]


def desc_set(m, label, ycol="t", k=10, rng_seed=SEED):
    m = m.reset_index(drop=True)
    x = m.FUSION_Z.values; y = m[ycol].values; Am = m.AveExpr.values
    nl = null_strat(x, y, decile_groups(Am, k), chrom_bounds(m.CHROM.values), NSHIFT, np.random.default_rng(rng_seed))
    d = full(x, y, Am, nl, label) if k == 10 else None
    if k == 1:
        o = strat(x, y, [np.arange(len(x))])
        se = jackknife(x, y, Am, k=1)
        d = dict(set=label, n=len(x), rho=o, se=se, lo=o - 1.96 * se, hi=o + 1.96 * se, p_shift=pval(nl, o), null_sd=nl.std(), null_mean=nl.mean())
    ndiag_post.append(null_diag(nl, len(x), label))
    return d


ndiag_post = []
desc.append(desc_set(alt, "D1 unstratified ALT", k=1))
desc.append(desc_set(alt[~alt_keep], "D2 COLOC complement (t primary)"))
for g in ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"]:
    lj = join(tw["ukbb_alt"], loco[g], cols=("t", "AveExpr"))
    lj = lj[lj.TargetID.isin(set(alt.TargetID))]  # genes shared with the ALT join (prespec D3)
    desc.append(desc_set(lj, f"D3 LOCO drop {g}"))
loo = []
for ch in sorted(alt.CHROM.unique()):
    mm = alt[alt.CHROM != ch]
    loo.append(dict(chrom=ch, n=len(mm), rho=strat(mm.FUSION_Z.values, mm.t.values, decile_groups(mm.AveExpr.values))))
loo = pd.DataFrame(loo)
desc.append(desc_set(J["ghodsian"], "D7 NAFLD (t primary)"))
desc.append(desc_set(J["ghodsian"], "D7 NAFLD (t adjusted)", ycol="t_adj"))
desc.append(desc_set(J["ukbb_ast"], "D8 AST"))
desc.append(desc_set(J["ukbb_ggt"], "D8 GGT"))
desc.append(desc_set(join(tw_p001, de), "D9 ALT OTTERS P0.001"))
desc.append(desc_set(alt[~mhc], "D10 ALT without MHC"))
desc.append(desc_set(ak, "D11 COLOC-excluded (t adjusted)", ycol="t_adj"))
desc.append(desc_set(alt[~alt_keep], "D11 COLOC complement (t adjusted)", ycol="t_adj"))
log("descriptives D1-D4, D7-D11 done")


# D5 / D6: strong-signal subset
def d5(ycol, label):
    S = np.abs(z) > 3
    sub = alt[S].reset_index(drop=True)
    x = sub.FUSION_Z.values; y = sub[ycol].values; Am = sub.AveExpr.values
    o = strat(x, y, decile_groups(Am))
    agr = float(np.mean(np.sign(x) == np.sign(y)))
    exp_ = float(np.mean(x > 0) * np.mean(y > 0) + np.mean(x < 0) * np.mean(y < 0))
    newblk = (sub.CHROM.diff().fillna(1) != 0) | (sub.mid.diff().fillna(1e9) > 1e6)
    blk = newblk.cumsum().values
    bl = [np.where(blk == b)[0] for b in np.unique(blk)]
    rng = np.random.default_rng(SEED)
    bo, ba = [], []
    for _ in range(NBOOT):
        ii = np.concatenate([bl[j] for j in rng.integers(0, len(bl), len(bl))])
        bo.append(strat(x[ii], y[ii], decile_groups(Am[ii])))
        xx, yy = x[ii], y[ii]
        ba.append(np.mean(np.sign(xx) == np.sign(yy)))
    # nominal shift p, regenerating the |z|>3 selection on every draw
    rng = np.random.default_rng(SEED)
    amts = shift_amounts(bnd, NSHIFT, rng)
    yfull = alt[ycol].values
    nl_r, nl_a = [], []
    for ii in (row for a in range(0, NSHIFT, 500) for row in shifted_index(bnd, amts[a:a + 500])):
        zs_ = z[ii]; s_ = np.abs(zs_) > 3
        xs_, ys_ = zs_[s_], yfull[s_]
        nl_r.append(strat(xs_, ys_, decile_groups(A[s_])))
        nl_a.append(np.mean(np.sign(xs_) == np.sign(ys_)) - (np.mean(xs_ > 0) * np.mean(ys_ > 0) + np.mean(xs_ < 0) * np.mean(ys_ < 0)))
    nl_r, nl_a = np.array(nl_r), np.array(nl_a)
    ndiag_post.append(null_diag(nl_r, int(S.sum()), label + " rho_strat"))
    return dict(set=label, n=int(S.sum()), n_blocks=len(bl), rho=o, lo=np.percentile(bo, 2.5), hi=np.percentile(bo, 97.5),
                p_shift=pval(nl_r, o), sign_agreement=agr, sign_agreement_expected_indep=exp_,
                sign_agreement_boot_lo=np.percentile(ba, 2.5), sign_agreement_boot_hi=np.percentile(ba, 97.5),
                p_shift_sign_agreement_minus_expected=pval(nl_a, agr - exp_))


d5p = d5("t", "D5 ALT |z|>3 (t primary)")
d5a = d5("t_adj", "D11 ALT |z|>3 (t adjusted)")
s6 = alt[(np.abs(z) > 3) & alt.canon.values]
k6 = int(np.sum(np.sign(s6.FUSION_Z) == np.sign(s6.t)))
w6 = wilson(k6, len(s6))
d6 = dict(set="D6 |z|>3 and canonical", n=len(s6), sign_agree=k6, frac=k6 / len(s6), wilson_lo=w6[0], wilson_hi=w6[1])
log("D5/D6 done")

# D12 figure source table
def wd_rank(y):
    out = np.empty(N)
    for g in G:
        out[g] = zrank(y[g])
    return out


bins = pd.qcut(pd.Series(z), 20, labels=False).values
yr_p, yr_a = wd_rank(t_p), wd_rank(t_a)
d12 = pd.DataFrame([dict(bin=b + 1, n=int((bins == b).sum()), z_min=z[bins == b].min(), z_max=z[bins == b].max(),
                         z_mean=z[bins == b].mean(), mean_std_rank_t_primary=yr_p[bins == b].mean(),
                         se_primary=yr_p[bins == b].std(ddof=1) / np.sqrt((bins == b).sum()),
                         mean_std_rank_t_adjusted=yr_a[bins == b].mean(),
                         se_adjusted=yr_a[bins == b].std(ddof=1) / np.sqrt((bins == b).sum())) for b in range(20)])

# ------------------------------------------------------------------ write Phase B
desc_df = pd.DataFrame(desc)
desc_df.to_csv(os.path.join(OUT, "estimates.tsv"), sep="\t", index=False)
loo.to_csv(os.path.join(OUT, "D4_leave_one_chromosome_out.tsv"), sep="\t", index=False)
pd.DataFrame([d5p, d5a]).to_csv(os.path.join(OUT, "D5_D11_absz3_subset.tsv"), sep="\t", index=False)
pd.DataFrame([d6]).to_csv(os.path.join(OUT, "D6_absz3_canonical_sign.tsv"), sep="\t", index=False)
d12.to_csv(os.path.join(OUT, "D12_figure_source_binned.tsv"), sep="\t", index=False)
pd.DataFrame(ndiag_post).to_csv(os.path.join(OUT, "diagnostics_nulls_post.tsv"), sep="\t", index=False)
post = pd.DataFrame([dict(item="SE_ratio_primary", value=prim["se_ratio"], flag=prim["se_ratio"] > 1.25),
                     dict(item="SE_ratio_qualifier", value=qual["se_ratio"], flag=qual["se_ratio"] > 1.25),
                     dict(item="jk_sign_flip_blocks_primary", value=prim["jk_sign_flip_blocks"], flag=""),
                     dict(item="jk_sign_flip_blocks_qualifier", value=qual["jk_sign_flip_blocks"], flag="")])
post.to_csv(os.path.join(OUT, "diagnostics_post_observation.tsv"), sep="\t", index=False)

f4 = lambda v: f"{v:.4f}"
loco_rows = desc_df[desc_df.set.str.startswith("D3")]
L = hdr + [f"## Verdict: reading {reading}: {verdict}", ""]
L += ["## Primary (family of one: ALT rho_strat, two-sided alpha 0.05; no correction needed)", "",
      f"- n genes {prim['n']}; rho_strat {f4(prim['rho'])}; jackknife 95% CI {f4(prim['lo'])} to {f4(prim['hi'])} (SE {prim['se']:.5f}, 200 blocks); "
      f"shift p {prim['p_shift']:.4g} ({NSHIFT} within-chromosome circular shifts; floor {1/(NSHIFT+1):.1e}).",
      f"- Shift-null SD {prim['null_sd']:.5f}; SE ratio (jackknife SE / null SD) {prim['se_ratio']:.3f}"
      + (" > 1.25: the shift null understates the variance; the direction decision rests on the jackknife CI." if prim["se_ratio"] > 1.25 else "."),
      f"- Jackknife blocks whose removal changes the sign of rho: {prim['jk_sign_flip_blocks']} of 200.", "",
      "## Composition qualifier (hepatocyte fraction added to the design)", "",
      f"- R0 reproduction: pass {r0_pass} (max abs t difference {rc['R0_max_abs_t_diff']}). VIF hepatocyte fraction {float(rc['R1_VIF_hepatocyte_fraction']):.3f}; "
      f"R1 canonical calls {rc['R1_canonical_calls_padj005_abslogFC05']}.",
      f"- Adjusted rho_strat {f4(qual['rho'])}; CI {f4(qual['lo'])} to {f4(qual['hi'])}; shift p {qual['p_shift']:.4g}; SE ratio {qual['se_ratio']:.3f}. "
      f"Qualifier passes (same sign, p<0.05, CI excludes 0): {qual_pass}. Adjusted CI inside +/-0.03: {inside(qual)}.", ""]
L += ["## Fixed-sequence step (COLOC-excluded set)", "",
      (f"- Tested: pass {fixed_seq}." if fixed_seq is not None else "- Not tested (primary did not pass with the qualifier); agreement beyond colocalized loci is not shown.")
      + f" COLOC-excluded rho {f4(coloc['rho'])}, CI {f4(coloc['lo'])} to {f4(coloc['hi'])}, shift p {coloc['p_shift']:.4g} (n {coloc['n']}).", ""]
if reading == 3:
    L += [f"- q* (smallest imposed sign agreement among ALT |z|>3 whose expected rho exceeds the upper 95% limit {f4(prim['hi'])}): {qstar}.",
          "- The bound does not exclude full sign agreement among the 103 |z|>3 canonical genes (expected rho about 0.0225).", ""]
L += ["## Stop rules", ""] + [f"- {r.rule}: {r.value}; triggered {r.triggered} ({r.consequence})" for r in stops.itertuples()] + [""]
L += ["## Controls (predicted before running)", ""] + [f"- {r.control}: stat {r.stat:.4g}; p {r.p:.4g}; n {r.n}; {r.extra}; predicted {r.predicted}; rule {r.rule}; pass {r.passed}" for r in ctrl.itertuples()]
L += ["- No control tests the real TWAS-to-DE pairing itself.", ""]
L += ["## Descriptive estimates (outside the family; every p here is nominal)", ""]
L += [f"- {r.set}: n {r.n}; rho {f4(r.rho)}; CI {f4(r.lo)} to {f4(r.hi)}; nominal shift p {r.p_shift:.4g}" for r in desc_df.iloc[3:].itertuples()]
L += [f"- D3 LOCO range of rho: {f4(loco_rows.rho.min())} to {f4(loco_rows.rho.max())}; fits with sign opposite to primary: "
      + (", ".join(loco_rows.loc[np.sign(loco_rows.rho) != np.sign(prim['rho']), 'set']) or "none"),
      f"- D4 leave-one-chromosome-out range: {f4(loo.rho.min())} to {f4(loo.rho.max())}",
      f"- D5 ALT |z|>3 (n {d5p['n']}, {d5p['n_blocks']} locus blocks): rho_strat {f4(d5p['rho'])}, block-bootstrap 95% CI {f4(d5p['lo'])} to {f4(d5p['hi'])}, nominal shift p {d5p['p_shift']:.4g}; "
      f"sign agreement {d5p['sign_agreement']:.4f} (expected under independence {d5p['sign_agreement_expected_indep']:.4f}; bootstrap CI {d5p['sign_agreement_boot_lo']:.4f} to {d5p['sign_agreement_boot_hi']:.4f}; nominal shift p {d5p['p_shift_sign_agreement_minus_expected']:.4g})",
      f"- D11 adjusted |z|>3: rho_strat {f4(d5a['rho'])}, CI {f4(d5a['lo'])} to {f4(d5a['hi'])}, nominal p {d5a['p_shift']:.4g}; sign agreement {d5a['sign_agreement']:.4f} vs {d5a['sign_agreement_expected_indep']:.4f}",
      f"- D6 |z|>3 and canonical: {d6['sign_agree']}/{d6['n']} same sign ({d6['frac']:.3f}; Wilson 95% {d6['wilson_lo']:.3f} to {d6['wilson_hi']:.3f})",
      "- D7 NAFLD estimates are labelled inconclusive at this signal level if inside +/-0.03; they are never a bound.", "",
      "Files: estimates.tsv, controls.tsv, stop_rules.tsv, diagnostics*.tsv, D4/D5/D6/D12 tables. Phase A hashes in diagnostics_complete.txt."]
write_result(L)
log("RESULT_H60.md written")
