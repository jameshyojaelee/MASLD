"""Test drivers for Model A-final (spec v1 section 5 with v1.2, v1.7 and v1.8): P1 conditional
randomization with supplied S* draws, P2 cluster bootstrap with supplied indices and lambda
draws, P3 conditional parametric bootstrap with supplied uniforms. Draws come from files so that
two implementations can be compared draw by draw.

Every refit is run from the warm start (the observed MLE) and the cold start (section 7), and
the kept fit follows the v1.8 section 2 item 3 table (fitting.select): eps_tie = 1e-6 for
section-3 fits and max(1e-6, 1e-8 G_fit) for section-4 fits, section-4 end points compared with
the reference integrator for every gene. A refit that the table marks failed is an exceedance.
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


def refit3(D, warm):
    """Section-3 refit from the observed MLE (warm) and the section 7 start (cold)."""
    fw, fc = ft.fit_section3(D, start=warm), ft.fit_section3(D)
    kept, rule = ft.select(fw, fc, EPS3)
    return kept, rule, fw, fc


def refit4(D, warm_theta4, v_max):
    """Section-4 refit: the section 5.3 two v starts from the observed section-4 MLE (warm) and
    from the section 7 start (cold); kept fit by the v1.8 rule on reference-integrator l."""
    w4 = ft.fit_section4_two_starts(D, warm_theta4, v_max)[0]
    c4 = ft.fit_section4_two_starts(D, ft.start_values(D), v_max)[0]
    kept, rule = ft.select(w4, c4, ft.eps_tie4(D), use_ref=True)
    return kept, rule, w4, c4


def _p(ind):
    return (1 + int(np.sum(ind))) / (len(ind) + 1)


# --- P1 ----------------------------------------------------------------------------------------

def run_p1(prob, fit3, Sstar):
    """z = kappa_hat / sandwich SE; refit section 3 for every S* draw; exceedance when
    |z*| >= |z| - 1e-6 or the refit fails (v1.8 start-selection table)."""
    D = prob.data()
    z, se, cond = kappa_z(fit3.theta, D)
    L = prob.L
    draws = []
    for b, S in enumerate(Sstar):
        Db = prob.data(S=S)
        f, rule, fw, fc = refit3(Db, fit3.theta)
        rec = dict(draw=b, rule=rule, loglik_warm=fw.loglik, loglik_cold=fc.loglik,
                   converged_warm=fw.converged, converged_cold=fc.converged,
                   kappa_warm=float(fw.theta[L.kappa]), kappa_cold=float(fc.theta[L.kappa]),
                   status_warm=fw.status, status_cold=fc.status, failed=f is None,
                   proj_max_grad_warm=fw.proj_max_grad, proj_max_grad_cold=fc.proj_max_grad)
        if f is None:
            rec.update(kept=None, kappa=None, se=None, z=None, exceed=1, in_band=False)
        else:
            zb, seb, _ = kappa_z(f.theta, Db)
            rec.update(kept="warm" if f is fw else "cold", kappa=float(f.theta[L.kappa]), se=seb,
                       z=zb, exceed=int(abs(zb) >= abs(z) - TIE),
                       in_band=bool(abs(abs(zb) - abs(z)) <= BAND))
        draws.append(rec)
    B = len(draws)
    p = _p([d["exceed"] for d in draws])
    return dict(kappa_hat=float(fit3.theta[L.kappa]), se=se, z=z, info_condition_number=cond, B=B,
                p=p, mc_se=mc_se(p, B), failed_draws=[d["draw"] for d in draws if d["failed"]],
                cold_kept_draws=[d["draw"] for d in draws if d["kept"] == "cold"],
                band_draws=[d["draw"] for d in draws if d["in_band"]], draws=draws)


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


def run_p2(prob, fit3, fit4, mult, lambdas, tau_max, Delta, v_scale):
    """kappa component per source from resample refits and lambda draws; tau component from
    the section-4 profile signed root at tau0 = tau_max * lambda^(0.05) (v1.2 item 1).
    l_p(tau_hat) and l_p(tau0) use the reference integrator for every gene (v1.8 section 2
    item 3); the l_p(tau0) profile is refit warm (observed section-4 MLE) and cold (section 7
    start) and kept by the v1.8 rule; a failed profile gives p_tau = 1."""
    L = prob.L
    recs = []
    for b, m in enumerate(mult):
        f, rule, fw, fc = refit3(prob.data(mult=m), fit3.theta)
        recs.append(dict(draw=b, rule=rule, kappa_warm=float(fw.theta[L.kappa]),
                         loglik_warm=fw.loglik, converged_warm=fw.converged,
                         kappa_cold=float(fc.theta[L.kappa]), loglik_cold=fc.loglik,
                         converged_cold=fc.converged, failed=f is None,
                         kept=None if f is None else ("warm" if f is fw else "cold"),
                         kappa=np.nan if f is None else float(f.theta[L.kappa])))
    B = len(recs)
    v_hat = float(fit4.theta[L.v])
    tau_hat = float(np.sqrt(v_hat))
    ll_hat = max(fit4.loglik_ref, fit3.loglik) if v_hat > 0 else fit3.loglik
    D = prob.data()
    kap = np.array([r["kappa"] for r in recs])
    conv = np.array([not r["failed"] for r in recs])
    out = dict(B=B, draws=recs, failed_draws=[r["draw"] for r in recs if r["failed"]],
               cold_kept_draws=[r["draw"] for r in recs if r["kept"] == "cold"],
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


def lr_fit(D, theta3_warm, theta4_warm, v_max):
    """Section-3 and section-4 refits with warm and cold starts (v1.8 section 2 item 3) and
    LR = 2 (max(l4_ref, l3) - l3)."""
    f3, rule3, f3w, f3c = refit3(D, theta3_warm)
    f4, rule4, f4w, f4c = refit4(D, theta4_warm, v_max)
    ok = f3 is not None and f4 is not None
    LR = ft.lr_statistic(f3, f4) if ok else None
    LR_gh = ft.lr_statistic_gh(f3, f4) if ok else None
    return dict(f3=f3, f4=f4, rule3=rule3, rule4=rule4, f3w=f3w, f3c=f3c, f4w=f4w, f4c=f4c,
                ok=ok, LR=LR, LR_gh=LR_gh)


def run_p3(prob, fit3, fit4, LR_obs, U, v_max):
    """LR* per draw with every included row (i, g, tag, n) held fixed; p = (1 + #{LR* >= LR -
    1e-6}) / (B + 1); a failed refit counts as an exceedance."""
    D = prob.data()
    L = prob.L
    draws = []
    for b in range(U.shape[0]):
        a_b, cls = simulate_counts(D, fit3.theta, U[b])
        Db = prob.data(rows=prob.rows.with_counts(a_b))
        r = lr_fit(Db, fit3.theta, fit4.theta, v_max)
        f3, f4 = r["f3"], r["f4"]
        rec = dict(draw=b, LR=r["LR"], LR_gh=r["LR_gh"], failed=not r["ok"],
                   rule3=r["rule3"], rule4=r["rule4"],
                   loglik3=None if f3 is None else f3.loglik,
                   loglik4_ref=None if f4 is None else f4.loglik_ref,
                   loglik4_gh=None if f4 is None else f4.loglik,
                   v_hat=None if f4 is None else float(f4.theta[L.v]),
                   v_at_upper=bool(f4 is not None and f4.v_at_upper),
                   kappa3=None if f3 is None else float(f3.theta[L.kappa]),
                   kept3=None if f3 is None else ("warm" if f3 is r["f3w"] else "cold"),
                   kept4=None if f4 is None else ("warm" if f4 is r["f4w"] else "cold"),
                   loglik3_warm=r["f3w"].loglik, loglik3_cold=r["f3c"].loglik,
                   loglik4_warm_ref=r["f4w"].loglik_ref, loglik4_cold_ref=r["f4c"].loglik_ref,
                   status4_warm=r["f4w"].status, status4_cold=r["f4c"].status,
                   n_reference_genes4=None if f4 is None else f4.n_reference_genes,
                   repolished4=bool(f4 is not None and f4.repolished),
                   n_class_H=int(np.sum(cls == 0)), n_class_R=int(np.sum(cls == 1)),
                   n_class_A=int(np.sum(cls == 2)), sum_a=int(a_b.sum()))
        rec["exceed"] = int((not r["ok"]) or r["LR"] >= LR_obs - TIE)
        rec["in_band"] = bool(r["ok"] and abs(r["LR"] - LR_obs) <= BAND)
        draws.append(rec)
    B = len(draws)
    p = _p([d["exceed"] for d in draws])
    return dict(LR_obs=LR_obs, B=B, p=p, mc_se=mc_se(p, B),
                failed_draws=[d["draw"] for d in draws if d["failed"]],
                band_draws=[d["draw"] for d in draws if d["in_band"]],
                v_upper_draws=[d["draw"] for d in draws if d["v_at_upper"]],
                LR_on_v_upper=[d["draw"] for d in draws if d["v_at_upper"] and d["LR"]],
                draws=draws)
