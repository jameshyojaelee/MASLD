"""Test drivers for Model A-final (spec v1 section 5 with v1.2, v1.7, v1.8 and v1.8.4): P1 conditional
randomization with supplied S* draws, P2 cluster bootstrap with supplied indices and lambda
draws, P3 conditional parametric bootstrap with supplied uniforms. Draws come from files so that
two implementations can be compared draw by draw.

v1.8.4 (SPEC_v184_ADDENDUM.md, H-18): every section-3 refit runs K = 3 starts (warm = the reported observed
MLE, cold = section 7 start of the refit data, probe j = 1 with the item 7 entropy) and keeps the start chosen
by fitting.select_k (eps_tie,3 = 1e-6). P3 section 4 is the item 6.1 pair from the draw's kept section-3 theta
(v = 0 and v = tau_max^2), kept by the v1.8 table at K = 2 with eps_tie4 on reference-integrator l; the LR uses
l3* (item 6.3) and item 6.4. A failed refit is an exceedance. P2 tau-profile fits keep v1.7 rule 5.

Every function also returns, under the key "_fits", the Fit objects behind each draw as (role, Fit, extra)
so that the caller can export every start's end point.
"""
import numpy as np
from scipy.stats import norm

import fitting as ft
import model_a as ma
from betabinom import bb_logpmf

TIE = 1e-6
BAND = 1e-3
LAMBDA_MIN = 0.05
EPS3 = 1e-6            # v1.8 section 2 item 3, section-3 fits


def mc_se(p, B):
    return float(np.sqrt(p * (1.0 - p) / B))


def kappa_z(theta, D):
    se, cond = ft.sandwich_se(theta, D)
    return float(theta[D.L.kappa] / se), se, cond


def probe_entropies_fixed(t, b):
    """B-FIX and real data: probe start j = 1 of draw b of test t: [[20260923, 13, t, b, 1]] (v1.8.4 item 7)."""
    return [ft.probe_entropy_fixed_refit(t, b, 1)]


def refit3(D, warm, probe_entropies):
    """v1.8.4 item 2: section-3 refit from warm, cold and the probe start(s); kept by fitting.select_k.
    Returns (kept Fit or None, rule, kept index or None, list of the K Fits)."""
    return ft.refit_section3(D, warm, probe_entropies)


def _start_rec(fits, prefix=""):
    """v1.8.4 item 10 record per start: l - l_warm, converged, attempts, seconds (and kappa, status)."""
    rec = {}
    l0 = fits[0].loglik
    for name, f in zip(ft.START_NAMES, fits):
        rec[f"{prefix}loglik_{name}"] = float(f.loglik)
        rec[f"{prefix}dl_{name}_vs_warm"] = float(f.loglik - l0)
        rec[f"{prefix}converged_{name}"] = bool(f.converged)
        rec[f"{prefix}attempts_{name}"] = int(f.attempts)
        rec[f"{prefix}seconds_{name}"] = float(f.seconds)
        rec[f"{prefix}status_{name}"] = str(f.status)
    return rec


def _fit_items(fits, base, D):
    return [(f"{base}/{n}", f, dict(kind="section3")) for n, f in zip(ft.START_NAMES, fits)]


def _p(ind):
    return (1 + int(np.sum(ind))) / (len(ind) + 1)


# --- P1 ----------------------------------------------------------------------------------------

def run_p1(prob, fit3, Sstar, probe_entropies=None):
    """z = kappa_hat / sandwich SE; refit section 3 for every S* draw (K = 3); exceedance when
    |z*| >= |z| - 1e-6 or the refit fails (fitting.select_k).
    probe_entropies(b) gives the probe entropies of draw b (default: B-FIX/real data, t = 1)."""
    pe = probe_entropies or (lambda b: probe_entropies_fixed(1, b))
    D = prob.data()
    z, se, cond = kappa_z(fit3.theta, D)
    L = prob.L
    draws, items = [], []
    for b, S in enumerate(Sstar):
        Db = prob.data(S=S)
        ent = pe(b)
        f, rule, k, fits = refit3(Db, fit3.theta, ent)
        items += [(b, *it) for it in _fit_items(fits, "refit3", Db)]
        rec = dict(draw=b, rule=rule, failed=f is None, probe_entropy=repr(ent))
        rec.update(_start_rec(fits))
        if f is None:
            rec.update(kept=None, kept_index=None, kappa=None, se=None, z=None, exceed=1, in_band=False)
        else:
            zb, seb, _ = kappa_z(f.theta, Db)
            rec.update(kept=ft.START_NAMES[k], kept_index=k, kappa=float(f.theta[L.kappa]), se=seb,
                       z=zb, exceed=int(abs(zb) >= abs(z) - TIE),
                       in_band=bool(abs(abs(zb) - abs(z)) <= BAND))
        draws.append(rec)
    B = len(draws)
    p = _p([d["exceed"] for d in draws])
    return dict(kappa_hat=float(fit3.theta[L.kappa]), se=se, z=z, info_condition_number=cond, B=B,
                p=p, mc_se=mc_se(p, B), failed_draws=[d["draw"] for d in draws if d["failed"]],
                cold_kept_draws=[d["draw"] for d in draws if d["kept"] == "cold"],
                probe_kept_draws=[d["draw"] for d in draws if d["kept"] == "probe1"],
                band_draws=[d["draw"] for d in draws if d["in_band"]],
                exceed_indicators=[d["exceed"] for d in draws], draws=draws, _fits=items)


# --- P2 ----------------------------------------------------------------------------------------

def _kappa_components(kap, conv, lam, Delta):
    """p_up and p_lo with invalid lambda draws and failed refits counted as exceedances."""
    invalid = ~np.isfinite(lam) | (lam <= LAMBDA_MIN)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = kap / lam
    bad = invalid | ~conv
    up = bad | (ratio >= Delta)
    lo = bad | (ratio <= -Delta)
    B = len(kap)
    return dict(p_up=_p(up), p_lo=_p(lo), p_up_mc_se=mc_se(_p(up), B),
                p_lo_mc_se=mc_se(_p(lo), B),
                ratio_interval90=[float(x) for x in np.quantile(ratio[~bad], [0.05, 0.95])],
                ratio_draws=[None if b else float(r) for r, b in zip(ratio, bad)],
                up_indicators=up.astype(int).tolist(), lo_indicators=lo.astype(int).tolist())


def run_p2(prob, fit3, fit4, mult, lambdas, tau_max, Delta, v_scale, probe_entropies=None):
    """kappa component per source from resample refits and lambda draws; tau component from
    the section-4 profile signed root at tau0 = tau_max * lambda^(0.05) (v1.2 item 1).
    l_p(tau_hat) and l_p(tau0) use the reference integrator for every gene (v1.8 section 2
    item 3); the l_p(tau0) profile is refit warm (observed section-4 MLE) and cold (section 7
    start) and kept by the v1.8 rule; a failed profile gives p_tau = 1."""
    pe = probe_entropies or (lambda b: probe_entropies_fixed(2, b))
    L = prob.L
    recs, items = [], []
    for b, m in enumerate(mult):
        Db = prob.data(mult=m)
        ent = pe(b)
        f, rule, k, fits = refit3(Db, fit3.theta, ent)
        items += [(b, *it) for it in _fit_items(fits, "refit3", Db)]
        rec = dict(draw=b, rule=rule, failed=f is None, probe_entropy=repr(ent),
                   kept=None if f is None else ft.START_NAMES[k], kept_index=k,
                   kappa=np.nan if f is None else float(f.theta[L.kappa]))
        rec.update(_start_rec(fits))
        for n, ff in zip(ft.START_NAMES, fits):
            rec[f"kappa_{n}"] = float(ff.theta[L.kappa])
        recs.append(rec)
    B = len(recs)
    v_hat = float(fit4.theta[L.v])
    tau_hat = float(np.sqrt(v_hat))
    ll_hat = max(fit4.loglik_ref, fit3.loglik) if v_hat > 0 else fit3.loglik
    D = prob.data()
    kap = np.array([r["kappa"] for r in recs])
    conv = np.array([not r["failed"] for r in recs])
    out = dict(B=B, draws=recs, failed_draws=[r["draw"] for r in recs if r["failed"]],
               cold_kept_draws=[r["draw"] for r in recs if r["kept"] == "cold"],
               probe_kept_draws=[r["draw"] for r in recs if r["kept"] == "probe1"], _fits=items,
               tau_hat=tau_hat, loglik_at_tau_hat=ll_hat, tau_hat_at_v_upper=bool(fit4.v_at_upper),
               sources={})
    for s, lam in lambdas.items():
        lam = np.asarray(lam, float)
        invalid = ~np.isfinite(lam) | (lam <= LAMBDA_MIN)
        lam05 = float(np.quantile(np.where(invalid, 0.0, lam), 0.05))
        rec = dict(invalid_lambda_draws=np.flatnonzero(invalid).tolist(), lambda_q05=lam05)
        rec.update(_kappa_components(kap, conv, lam, Delta))
        if lam05 <= LAMBDA_MIN:
            rec.update(tau0=None, p_tau=1.0, note="lambda_q05 <= 0.05")
        else:
            tau0 = tau_max * lam05
            fpw = ft.fit_section4(D, fit4.theta, tau0 ** 2, v_scale, fix_v=True)
            fpc = ft.fit_section4(D, ft.start_values(D), tau0 ** 2, v_scale, fix_v=True)
            fp, rule = ft.select(fpw, fpc, ft.eps_tie4(D), use_ref=True)
            items.append(("", f"profile_tau0.{s}/warm", fpw, dict(kind="section4_profile")))
            items.append(("", f"profile_tau0.{s}/cold", fpc, dict(kind="section4_profile")))
            rec.update(tau0=tau0, profile_rule=rule,
                       loglik_profile_tau0_warm_ref=fpw.loglik_ref,
                       loglik_profile_tau0_cold_ref=fpc.loglik_ref,
                       loglik_profile_tau0_warm_gh=fpw.loglik,
                       profile_converged_warm=fpw.converged, profile_converged_cold=fpc.converged,
                       profile_status_warm=fpw.status, profile_status_cold=fpc.status)
            if fp is None:
                rec.update(loglik_diff=None, signed_root=None, p_tau=1.0,
                           note="profile refit failed: p_tau = 1")
            else:
                diff = ll_hat - fp.loglik_ref
                r = float(np.sign(tau0 - tau_hat) * np.sqrt(2.0 * max(diff, 0.0)))
                rec.update(loglik_diff=diff, signed_root=r, p_tau=float(norm.sf(r)),
                           profile_kept="warm" if fp is fpw else "cold")
        rec["p_source"] = max(rec["p_up"], rec["p_lo"], rec["p_tau"])
        out["sources"][s] = rec
    out["p"] = max(r["p_source"] for r in out["sources"].values())
    return out


# --- P3 ----------------------------------------------------------------------------------------

def _inverse_cdf(cdf, u):
    """Smallest index k with cdf[k] >= u (the last index if rounding leaves u above cdf[-1])."""
    return min(int(np.searchsorted(cdf, u, side="left")), len(cdf) - 1)


def simulate_counts(D, theta, U):
    """Section 5.3: per row, class from U[:, 0] over (H, R, A) with weights P(G) C_G(n), then a
    from U[:, 1] over K_c(n) in increasing a under the truncated f_G."""
    R = D.rows
    x1, x2, _, _ = ma.predictors(theta, D)
    logC_H = ma.row_eval(x1, x2, R, 0).logC_H
    plus = R.orientation == 1
    logC_R = np.where(plus, R.logC_O0, R.logC_O1)
    logC_A = np.where(plus, R.logC_O1, R.logC_O0)
    logw = np.stack([R.logPH + logC_H, R.logPR + logC_R, R.logPA + logC_A], axis=1)
    mu_H = ma.expit(x1)
    rho_H = ma.EPS_RHO + (1.0 - 2.0 * ma.EPS_RHO) * ma.expit(x2)
    a_new = np.empty(len(R.n), dtype=np.int64)
    cls = np.empty(len(R.n), dtype=np.int64)
    for i in range(len(R.n)):
        pw = np.exp(logw[i] - logw[i].max())
        c = _inverse_cdf(np.cumsum(pw) / pw.sum(), U[i, 0])
        cls[i] = c
        a = np.arange(R.lo[i], R.hi[i] + 1)
        if c == 0:
            lf = bb_logpmf(a, R.n[i], mu_H[i], rho_H[i])
        else:
            # R: tag-REF homozygote; with orientation +1 a counts tag ALT (errors), else tag REF
            mu = R.e[i] if (c == 1) == plus[i] else 1.0 - R.e[i]
            lf = bb_logpmf(a, R.n[i], mu, R.phi)
        pm = np.exp(lf - lf.max())
        a_new[i] = a[_inverse_cdf(np.cumsum(pm) / pm.sum(), U[i, 1])]
    return a_new, cls


def section4_pair_lr(D, theta3_kept, l3, v_max):
    """v1.8.4 items 6.1-6.4 for one data set: the pair from the kept section-3 theta, the kept member,
    l3* and LR. Returns a dict (kept4 None when the pair table fails)."""
    kept, f0, f1, rule4 = ft.fit_section4_two_starts(D, theta3_kept, v_max)
    ok4 = kept is f0 or kept is f1
    l3s, contrib = ft.ell3_star(l3, [f0, f1], D)
    out = dict(kept4=kept if ok4 else None, member0=f0, member1=f1, rule4=rule4, l3star=l3s,
               l3star_contrib=contrib, ok4=ok4,
               kept4_index=None if not ok4 else (0 if kept is f0 else 1))
    out["LR"] = ft.lr_v184(kept, l3s, D) if ok4 else None
    out["LR_gh"] = ft.lr_v184_gh(kept, l3s, D) if ok4 else None
    out["missed_v0"] = bool(l3s - l3 > ft.MISSED_V0_MARGIN)
    return out


def lr_fit(D, theta3_warm, v_max, probe_entropies):
    """P3 draw: section-3 refit at K = 3 (item 2), then the item 6.1 pair from the kept section-3 theta,
    l3* (6.3) and LR (6.4)."""
    f3, rule3, k3, fits3 = refit3(D, theta3_warm, probe_entropies)
    out = dict(f3=f3, rule3=rule3, k3=k3, fits3=fits3, ok=False, LR=None, LR_gh=None, pair=None)
    if f3 is None:
        return out
    pr = section4_pair_lr(D, f3.theta, f3.loglik, v_max)
    out.update(pair=pr, ok=pr["ok4"], LR=pr["LR"], LR_gh=pr["LR_gh"])
    return out


def run_p3(prob, fit3, LR_obs, U, v_max, probe_entropies=None):
    """LR* per draw with every included row (i, g, tag, n) held fixed; draws generated at the reported observed
    section-3 theta (fit3); p = (1 + #{LR* >= LR - 1e-6}) / (B + 1); a failed draw counts as an exceedance."""
    pe = probe_entropies or (lambda b: probe_entropies_fixed(3, b))
    D = prob.data()
    L = prob.L
    draws, items, counts = [], [], []
    for b in range(U.shape[0]):
        a_b, cls = simulate_counts(D, fit3.theta, U[b])
        counts.append((a_b.copy(), cls.copy()))
        Db = prob.data(rows=prob.rows.with_counts(a_b))
        ent = pe(b)
        r = lr_fit(Db, fit3.theta, v_max, ent)
        f3, pr = r["f3"], r["pair"]
        items += [(b, *it) for it in _fit_items(r["fits3"], "refit3", Db)]
        rec = dict(draw=b, LR=r["LR"], LR_gh=r["LR_gh"], failed=not r["ok"], probe_entropy=repr(ent),
                   rule3=r["rule3"], kept3=None if f3 is None else ft.START_NAMES[r["k3"]],
                   kept3_index=r["k3"], loglik3=None if f3 is None else f3.loglik,
                   kappa3=None if f3 is None else float(f3.theta[L.kappa]),
                   n_class_H=int(np.sum(cls == 0)), n_class_R=int(np.sum(cls == 1)),
                   n_class_A=int(np.sum(cls == 2)), sum_a=int(a_b.sum()))
        rec.update(_start_rec(r["fits3"], "s3_"))
        if pr is not None:
            f4 = pr["kept4"]
            for i, (nm, m) in enumerate((("v0", pr["member0"]), ("vmax", pr["member1"]))):
                l0 = pr["l3star_contrib"][i][1]
                items.append((b, f"pair/{nm}", m, dict(kind="section4", l3_at_v0=l0)))
                rec[f"m{i}_v_hat"] = float(m.theta[L.v])
                rec[f"m{i}_loglik_ref"] = float(m.loglik_ref)
                rec[f"m{i}_loglik_gh"] = float(m.loglik)
                rec[f"m{i}_converged"] = bool(m.converged)
                rec[f"m{i}_attempts"] = int(m.attempts)
                rec[f"m{i}_seconds"] = float(m.seconds)
                rec[f"m{i}_status"] = str(m.status)
                rec[f"m{i}_l3_at_v0"] = l0
            rec.update(rule4=pr["rule4"], kept4=None if f4 is None else ("v0", "vmax")[pr["kept4_index"]],
                       loglik4_ref=None if f4 is None else f4.loglik_ref,
                       loglik4_gh=None if f4 is None else f4.loglik,
                       v_hat=None if f4 is None else float(f4.theta[L.v]),
                       v_at_upper=bool(f4 is not None and f4.v_at_upper),
                       n_reference_genes4=None if f4 is None else f4.n_reference_genes,
                       repolished4=bool(f4 is not None and f4.repolished),
                       l3star=pr["l3star"], l3star_minus_l3=pr["l3star"] - f3.loglik,
                       missed_v0=pr["missed_v0"])
        rec["exceed"] = int((not r["ok"]) or r["LR"] >= LR_obs - TIE)
        rec["in_band"] = bool(r["ok"] and abs(r["LR"] - LR_obs) <= BAND)
        draws.append(rec)
    B = len(draws)
    p = _p([d["exceed"] for d in draws])
    return dict(LR_obs=LR_obs, B=B, p=p, mc_se=mc_se(p, B),
                failed_draws=[d["draw"] for d in draws if d["failed"]],
                band_draws=[d["draw"] for d in draws if d["in_band"]],
                v_upper_draws=[d["draw"] for d in draws if d.get("v_at_upper")],
                LR_on_v_upper=[d["draw"] for d in draws if d.get("v_at_upper") and d["LR"]],
                missed_v0_draws=[d["draw"] for d in draws if d.get("missed_v0")],
                cold_kept_draws=[d["draw"] for d in draws if d["kept3"] == "cold"],
                probe_kept_draws=[d["draw"] for d in draws if d["kept3"] == "probe1"],
                exceed_indicators=[d["exceed"] for d in draws], draws=draws, _fits=items, _counts=counts)
