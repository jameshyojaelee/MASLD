"""N2a analysis (PRESPEC_N2a.md). Runs only for a set whose gate passed.
Usage: python -I 06_analysis.py RUN_ROOT GATE_DIR OUT_DIR SET     (SET = primary | secondary; OUT_DIR must not exist)
Unit: the woman. Family: T0, T1, T2 with Holm. Seed 20260923 for every draw."""
import sys, os, csv, json, collections
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from n2a_lib import load_counts, site_record

root, gate_dir, outdir, SET = sys.argv[1:5]
SEED, REPS, REPS_B = 20260923, 1000, 200
gate = json.load(open(os.path.join(gate_dir, "gate_summary.json")))
if not gate[SET]["gate_pass"]:
    sys.exit("gate not passed for " + SET)
if os.path.exists(outdir):
    sys.exit("refusing to overwrite " + outdir)
os.makedirs(outdir)
rng = np.random.default_rng(SEED)

cls_of = {}
with open(os.path.join(root, "inputs", "sites.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        cls_of[(r["chrom"], int(r["pos"]))] = r["cls"]
pairs = [r for r in csv.DictReader(open(os.path.join(root, "inputs", "n2a_pairs.tsv")), delimiter="\t") if r["include"] == "True" and r["set"] == SET]
pairs.sort(key=lambda r: (r["cohort"], r["person_id"]))
gate_k = {r["person_id"]: int(r["k_x_informative"]) for r in csv.DictReader(open(os.path.join(gate_dir, "gate_per_woman.tsv")), delimiter="\t")}
pairs = [r for r in pairs if gate_k[r["person_id"]] >= 5]


def raw_counts(a, b):
    """(ref1, depth1, ref2, depth2) over all reads, REF allele as the counted allele."""
    return a[0] + a[2], a[0] + a[1] + a[2] + a[3], b[0] + b[2], b[0] + b[1] + b[2] + b[3]


# columns: bin, x1, n1, x2, n2, ref1, depth1, ref2, depth2, is_inactive_class
W = []
for pr in pairs:
    c1 = load_counts(os.path.join(root, "counts", pr["lib1"] + ".counts.tsv.gz"))
    c2 = load_counts(os.path.join(root, "counts", pr["lib2"] + ".counts.tsv.gz"))
    d = dict(pr=pr, X=[], Xs=[], A=[])
    for key in sorted(c1.keys() & c2.keys()):
        rec, recs = site_record(c1[key], c2[key]), site_record(c2[key], c1[key])
        if rec is None:
            continue
        rc = raw_counts(c1[key], c2[key])
        if cls_of[key].startswith("X_"):
            d["X"].append(rec + rc + (1 if cls_of[key] == "X_inactive" else 0,))
            d["Xs"].append(recs + (0, 0, 0, 0, 0))
        else:
            d["A"].append(rec + rc + (0,))
    d["X"], d["Xs"], d["A"] = (np.array(d[k], dtype=float).reshape(-1, 10) for k in ("X", "Xs", "A"))
    W.append(d)


def delta(arr):
    n1, n2 = arr[:, 2].sum(), arr[:, 4].sum()
    return arr[:, 3].sum() / n2 - arr[:, 1].sum() / n1 if n1 > 0 and n2 > 0 else np.nan


# ---- beta-binomial machinery (mean p, intraclass correlation rho) ----
def bb_ll(x, n, p, rho):
    s = 1.0 / rho - 1.0
    a, b = p * s, (1 - p) * s
    return (gammaln(n + 1) - gammaln(x + 1) - gammaln(n - x + 1) + gammaln(x + a) + gammaln(n - x + b) - gammaln(n + a + b)
            + gammaln(a + b) - gammaln(a) - gammaln(b))


def fit_rho(groups, mixture=False):
    """groups: list of (x, n). One p per group (or a mixture of p and 1-p over sites with p in [0.5, 1]); shared rho. ML."""
    def group_ll(x, n, rho):
        if mixture:
            f = lambda th: -np.sum(np.logaddexp(bb_ll(x, n, th, rho), bb_ll(x, n, 1 - th, rho)) - np.log(2))
            r = minimize_scalar(f, bounds=(0.5, 0.999), method="bounded")
        else:
            f = lambda p: -np.sum(bb_ll(x, n, p, rho))
            r = minimize_scalar(f, bounds=(1e-4, 1 - 1e-4), method="bounded")
        return -r.fun, r.x

    def total(lr):
        rho = 10 ** lr
        return -sum(group_ll(x, n, rho)[0] for x, n in groups)
    r = minimize_scalar(total, bounds=(-5, -0.5), method="bounded")
    rho = 10 ** r.x
    return rho, [group_ll(x, n, rho)[1] for x, n in groups]


def var_ratio(n, p, rho):
    p = min(max(p, 0.01), 0.99)
    return np.sum(n * p * (1 - p) * (1 + (n - 1) * rho)) / np.sum(n) ** 2


# ---- per-woman Delta, overdispersion per cohort ----
cohorts = sorted({d["pr"]["cohort"] for d in W})
rho_c = {}
for c in cohorts:
    groups = []
    for d in W:
        if d["pr"]["cohort"] == c and len(d["X"]):
            a = d["X"]; a = a[a[:, 2] > 0]
            groups += [(a[:, 1], a[:, 2]), (a[:, 3], a[:, 4])]
    rho_c[c] = fit_rho(groups)[0]
res = []
for d in W:
    a = d["X"]; c = d["pr"]["cohort"]; a_ok = a[a[:, 2] > 0]
    x1, n1, x2, n2 = a_ok[:, 1], a_ok[:, 2], a_ok[:, 3], a_ok[:, 4]
    R1, R2 = x1.sum() / n1.sum(), x2.sum() / n2.sum()
    d["delta"], d["R1"], d["R2"] = R2 - R1, R1, R2
    d["se"] = float(np.sqrt(var_ratio(n1, R1, rho_c[c]) + var_ratio(n2, R2, rho_c[c])))
    d["k"], d["k_valid"] = len(a), len(a_ok)
    s = d["Xs"]; s = s[s[:, 2] > 0]
    d["delta_swapped"] = delta(s)
    ina = a[(a[:, 9] == 1) & (a[:, 2] > 0)]
    d["delta_inactive"] = delta(ina) if len(ina) >= 5 else np.nan
    d["k_inactive"] = len(ina)
    d["abs_stage"] = float(d["pr"]["abs_stage_change"])


# ---- autosomal technical-noise null: matched draws of her autosomal informative sites ----
def draw_idx(r, xbins, A, reps):
    cols = []
    short = 0
    for b, m in xbins.items():
        pool = np.where(A[:, 0] == b)[0]
        if len(pool) == 0:
            pool = np.arange(len(A)); short += 1
        if len(pool) >= m:
            cols.append(pool[np.argsort(r.random((reps, len(pool))), axis=1)[:, :m]])
        else:
            cols.append(pool[r.integers(0, len(pool), size=(reps, m))]); short += 1
    return np.concatenate(cols, axis=1), short


def delta_from_idx(A, idx):
    with np.errstate(divide="ignore", invalid="ignore"):
        return A[idx, 3].sum(1) / A[idx, 4].sum(1) - A[idx, 1].sum(1) / A[idx, 2].sum(1)


for d in W:
    xb = collections.Counter(d["X"][:, 0].astype(int)); xb = {float(k): v for k, v in xb.items()}
    d["idx"], d["short"] = draw_idx(rng, xb, d["A"], REPS)
    d["null"] = delta_from_idx(d["A"], d["idx"])
    d["p_null"] = (1 + np.sum(np.abs(d["null"]) >= abs(d["delta"]))) / (REPS + 1)

tests = []
Dl = np.array([d["delta"] for d in W]); SE = np.array([d["se"] for d in W]); n_w = len(W)
# T0: DerSimonian-Laird mean of signed delta against 0
w = 1 / SE ** 2
mu_f = np.sum(w * Dl) / w.sum(); Q = np.sum(w * (Dl - mu_f) ** 2)
tau2 = max(0.0, (Q - (n_w - 1)) / (w.sum() - np.sum(w ** 2) / w.sum()))
ws = 1 / (SE ** 2 + tau2); mu = np.sum(ws * Dl) / ws.sum(); se_mu = 1 / np.sqrt(ws.sum())
p0 = 2 * stats.norm.sf(abs(mu / se_mu))
wil = stats.wilcoxon(Dl)
tests.append(dict(test="T0 signed delta vs 0 (DL random effects)", n_women=n_w, estimate=mu, se=se_mu, p_raw=p0, note="tau2=%.5f; Wilcoxon p=%.4g" % (tau2, wil.pvalue)))
# T1: mean |delta| against the autosomal null
obs = np.mean(np.abs(Dl))
null_means = np.nanmean([np.abs(d["null"]) for d in W], axis=0)
p1 = (1 + np.sum(null_means >= obs)) / (REPS + 1)
tests.append(dict(test="T1 mean |delta| vs autosomal null", n_women=n_w, estimate=obs, se=float(np.std(null_means, ddof=1)), p_raw=p1,
                  note="null mean %.4f; null 95th pct %.4f" % (np.mean(null_means), np.percentile(null_means, 95))))
# T2: slope of |delta| on |stage change| with cohort fixed effect; permutation within cohort
y = np.abs(Dl); x = np.array([d["abs_stage"] for d in W]); coh = np.array([d["pr"]["cohort"] for d in W])


def slope(y, x, coh):
    cols = [np.ones(len(y)), x] + [(coh == c).astype(float) for c in cohorts[1:] if (coh == c).any()]
    return np.linalg.lstsq(np.column_stack(cols), y, rcond=None)[0][1]


s_obs = slope(y, x, coh)
perm = np.empty(REPS)
for r in range(REPS):
    xp = x.copy()
    for c in cohorts:
        ix = np.where(coh == c)[0]
        xp[ix] = x[rng.permutation(ix)]
    perm[r] = slope(y, xp, coh)
p2 = (1 + np.sum(perm >= s_obs)) / (REPS + 1)
rho_s = stats.spearmanr(x, y)
tests.append(dict(test="T2 slope of |delta| on |stage change|, cohort fixed effect", n_women=n_w, estimate=s_obs, se=float(np.std(perm, ddof=1)), p_raw=p2,
                  note="one-sided slope>0; Spearman rho=%.3f p=%.3g; stage-change values: %s" % (rho_s.statistic, rho_s.pvalue, {int(k): v for k, v in collections.Counter(x.astype(int)).items()})))
order = np.argsort([t["p_raw"] for t in tests]); m = len(tests); run = 0
for rank, i in enumerate(order):
    run = max(run, (m - rank) * tests[i]["p_raw"]); tests[i]["p_holm"] = min(1.0, run)

# ---- sensitivity (outside the family) ----
sens = []
sw = np.array([d["delta_swapped"] for d in W])
sens.append(dict(analysis="swapped direction (phase by biopsy 2)", n_women=int(np.isfinite(sw).sum()), mean_abs_delta=float(np.nanmean(np.abs(sw))),
                 slope_abs_delta_on_abs_stage=float(slope(np.abs(sw)[np.isfinite(sw)], x[np.isfinite(sw)], coh[np.isfinite(sw)]))))
ia = np.array([d["delta_inactive"] for d in W]); ok = np.isfinite(ia)
sens.append(dict(analysis="genes classed inactive only (>=5 sites)", n_women=int(ok.sum()), mean_abs_delta=float(np.mean(np.abs(ia[ok]))) if ok.any() else np.nan,
                 slope_abs_delta_on_abs_stage=float(slope(np.abs(ia[ok]), x[ok], coh[ok])) if ok.sum() > 3 else np.nan))
for c in cohorts:
    k = coh == c
    sens.append(dict(analysis="cohort " + c, n_women=int(k.sum()), mean_abs_delta=float(np.mean(y[k])),
                     slope_abs_delta_on_abs_stage=float(np.polyfit(x[k], y[k], 1)[0]) if k.sum() > 3 and np.ptp(x[k]) > 0 else np.nan))
loo = [slope(np.delete(y, i), np.delete(x, i), np.delete(coh, i)) for i in range(n_w)]
sens.append(dict(analysis="leave-one-woman-out slope range", n_women=n_w - 1, mean_abs_delta=np.nan, slope_abs_delta_on_abs_stage="%.4f to %.4f" % (min(loo), max(loo))))

# magnitude-only mixture: per biopsy, per-site phase integrated out, skew in [0.5, 1]
rho_mix_c = {}
for c in cohorts:
    g = []
    for d in W:
        if d["pr"]["cohort"] == c:
            g += [(d["X"][:, 5], d["X"][:, 6]), (d["X"][:, 7], d["X"][:, 8])]
    rho_mix_c[c] = fit_rho(g, mixture=True)[0]


def theta_hat(xr, n, rho):
    f = lambda th: -np.sum(np.logaddexp(bb_ll(xr, n, th, rho), bb_ll(xr, n, 1 - th, rho)))
    return minimize_scalar(f, bounds=(0.5, 0.999), method="bounded").x


for d in W:
    c = d["pr"]["cohort"]; rho = rho_mix_c[c]
    d["m1"], d["m2"] = (theta_hat(d["X"][:, 5], d["X"][:, 6], rho) - 0.5, theta_hat(d["X"][:, 7], d["X"][:, 8], rho) - 0.5)
    d["dm"] = d["m2"] - d["m1"]
    nulls = []
    for r in range(REPS_B):
        ix = d["idx"][r]
        nulls.append(theta_hat(d["A"][ix, 7], d["A"][ix, 8], rho) - theta_hat(d["A"][ix, 5], d["A"][ix, 6], rho))
    d["dm_null"] = np.array(nulls)
    d["p_dm"] = (1 + np.sum(np.abs(d["dm_null"]) >= abs(d["dm"]))) / (REPS_B + 1)
dm = np.array([d["dm"] for d in W])
nm = np.mean([np.abs(d["dm_null"]) for d in W], axis=0)
sens.append(dict(analysis="magnitude-only mixture: mean |m2 - m1| vs autosomal null (200 draws)", n_women=n_w, mean_abs_delta=float(np.mean(np.abs(dm))),
                 slope_abs_delta_on_abs_stage="p=%.4g; slope %.4f" % ((1 + np.sum(nm >= np.mean(np.abs(dm)))) / (REPS_B + 1), slope(np.abs(dm), x, coh))))

# ---- write ----
def wr(path, rows, cols):
    with open(os.path.join(outdir, path), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(("%.6g" % r[c]) if isinstance(r.get(c), (float, np.floating)) else str(r.get(c, "")) for c in cols) + "\n")


wr("tests.tsv", tests, ["test", "n_women", "estimate", "se", "p_raw", "p_holm", "note"])
wr("sensitivity.tsv", sens, ["analysis", "n_women", "mean_abs_delta", "slope_abs_delta_on_abs_stage"])
rows = []
for d in W:
    pr = d["pr"]
    rows.append(dict(cohort=pr["cohort"], person_id=pr["person_id"], k_x_sites=d["k"], k_valid=d["k_valid"], R1_crossfit=d["R1"], R2=d["R2"], delta=d["delta"],
                     se_bb=d["se"], p_vs_own_autosomal_null=d["p_null"], auto_null_sd=float(np.std(d["null"], ddof=1)), delta_swapped=d["delta_swapped"],
                     k_inactive=d["k_inactive"], delta_inactive_only=d["delta_inactive"], stage1=pr["stage1"], stage2=pr["stage2"], abs_stage_change=d["abs_stage"],
                     m1=d["m1"], m2=d["m2"], delta_m=d["dm"], p_delta_m_vs_autosomal=d["p_dm"], auto_pool_shortfall=d["short"]))
wr("per_woman.tsv", rows, list(rows[0].keys()))
json.dump(dict(set=SET, n_women=n_w, seed=SEED, reps=REPS, reps_mixture=REPS_B, rho_bb_by_cohort=rho_c, rho_mixture_by_cohort=rho_mix_c,
               note="T0 two-sided; T1, T2 one-sided; Holm over T0, T1, T2"), open(os.path.join(outdir, "summary.json"), "w"), indent=1)
for t in tests:
    print("%-62s n=%d est=%.4f p=%.4g holm=%.4g | %s" % (t["test"], t["n_women"], t["estimate"], t["p_raw"], t["p_holm"], t["note"]))
for s in sens:
    print(s)
