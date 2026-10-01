"""Test drivers for Model A-final (spec v1 section 5 with v1.2): P1 conditional randomization
with supplied S* draws, P2 cluster bootstrap with supplied indices and lambda draws, P3
conditional parametric bootstrap with supplied uniforms. Draws come from files so that two
implementations can be compared draw by draw.

Every refit is done twice: from the observed MLE (spec section 7, "warm", the primary
result) and from the section 7 cold start. The "best" policy keeps the converged fit with the
higher log-likelihood. The section-3 likelihood can have several local maxima in a resample
(dispersion intercepts and slopes on their bounds), so the two policies can differ.
"""
import numpy as np
from scipy.stats import norm

import fitting as ft
import model_a as ma
from betabinom import bb_logpmf

TIE = 1e-6
BAND = 1e-3
LAMBDA_MIN = 0.05
BETTER = 1e-6          # log-likelihood gain that counts as a different optimum


def mc_se(p, B):
    return float(np.sqrt(p * (1.0 - p) / B))


def kappa_z(theta, D):
    se, cond = ft.sandwich_se(theta, D)
    return float(theta[D.L.kappa] / se), se, cond


def refits3(D, warm):
    """Section-3 refit warm-started from the observed MLE and from the section 7 start."""
    return ft.fit_section3(D, start=warm), ft.fit_section3(D)


def better(a, b):
    """Converged before failed, then the higher log-likelihood."""
    if a.converged != b.converged:
        return a if a.converged else b
    return b if b.loglik > a.loglik + BETTER else a


def _p(ind):
    return (1 + int(np.sum(ind))) / (len(ind) + 1)


# --- P1 ----------------------------------------------------------------------------------------

def run_p1(prob, fit3, Sstar):
    """z = kappa_hat / sandwich SE; refit section 3 for every S* draw; exceedance when
    |z*| >= |z| - 1e-6 or the refit fails."""
    D = prob.data()
    z, se, cond = kappa_z(fit3.theta, D)
    L = prob.L
    draws = []
    for b, S in enumerate(Sstar):
        Db = prob.data(S=S)
        fw, fc = refits3(Db, fit3.theta)
        rec = dict(draw=b, loglik_warm=fw.loglik, loglik_cold=fc.loglik,
                   kappa_cold=float(fc.theta[L.kappa]),
                   cold_better=bool(fc.converged and fc.loglik > fw.loglik + BETTER))
        for pol, f in (("warm", fw), ("best", better(fw, fc))):
            rec[f"converged_{pol}"] = f.converged
            rec[f"kappa_{pol}"] = float(f.theta[L.kappa])
            if f.converged:
                zb, seb, _ = kappa_z(f.theta, Db)
                rec.update({f"se_{pol}": seb, f"z_{pol}": zb,
                            f"exceed_{pol}": int(abs(zb) >= abs(z) - TIE),
                            f"in_band_{pol}": bool(abs(abs(zb) - abs(z)) <= BAND)})
            else:
                rec.update({f"se_{pol}": None, f"z_{pol}": None, f"exceed_{pol}": 1,
                            f"in_band_{pol}": False})
        rec["proj_max_grad_warm"] = fw.proj_max_grad
        draws.append(rec)
    B = len(draws)
    out = dict(kappa_hat=float(fit3.theta[L.kappa]), se=se, z=z, info_condition_number=cond, B=B,
               cold_better_draws=[d["draw"] for d in draws if d["cold_better"]], draws=draws)
    for pol in ("warm", "best"):
        p = _p([d[f"exceed_{pol}"] for d in draws])
        out[f"p_{pol}"] = p
        out[f"mc_se_{pol}"] = mc_se(p, B)
        out[f"failed_draws_{pol}"] = [d["draw"] for d in draws if not d[f"converged_{pol}"]]
        out[f"band_draws_{pol}"] = [d["draw"] for d in draws if d[f"in_band_{pol}"]]
    return out


# --- P2 ----------------------------------------------------------------------------------------

def _kappa_components(kap, conv, lam, Delta):
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
    the section-4 profile signed root at tau0 = tau_max * lambda^(0.05) (v1.2 item 1)."""
    L = prob.L
    recs = []
    for b, m in enumerate(mult):
        fw, fc = refits3(prob.data(mult=m), fit3.theta)
        fb = better(fw, fc)
        recs.append(dict(draw=b, kappa_warm=float(fw.theta[L.kappa]), loglik_warm=fw.loglik,
                         converged_warm=fw.converged, kappa_cold=float(fc.theta[L.kappa]),
                         loglik_cold=fc.loglik, converged_cold=fc.converged,
                         kappa_best=float(fb.theta[L.kappa]), converged_best=fb.converged,
                         cold_better=bool(fc.converged and fc.loglik > fw.loglik + BETTER)))
    B = len(recs)
    tau_hat = float(np.sqrt(fit4.theta[L.v]))
    ll_hat = max(fit4.loglik, fit3.loglik)
    D = prob.data()
    out = dict(B=B, draws=recs, cold_better_draws=[r["draw"] for r in recs if r["cold_better"]],
               tau_hat=tau_hat, loglik_at_tau_hat=ll_hat, sources={})
    for s, lam in lambdas.items():
        lam = np.asarray(lam, float)
        invalid = ~np.isfinite(lam) | (lam <= LAMBDA_MIN)
        lam05 = float(np.quantile(np.where(invalid, 0.0, lam), 0.05))
        rec = dict(invalid_lambda_draws=np.flatnonzero(invalid).tolist(), lambda_q05=lam05)
        for pol in ("warm", "best"):
            kap = np.array([r[f"kappa_{pol}"] for r in recs])
            conv = np.array([r[f"converged_{pol}"] for r in recs])
            rec[pol] = _kappa_components(kap, conv, lam, Delta)
        if lam05 <= LAMBDA_MIN:
            rec.update(tau0=None, p_tau_warm=1.0, p_tau_best=1.0, note="lambda_q05 <= 0.05")
        else:
            tau0 = tau_max * lam05
            fpw = ft.fit_section4(D, fit4.theta, tau0 ** 2, v_scale, fix_v=True)
            fpc = ft.fit_section4(D, ft.start_values(D), tau0 ** 2, v_scale, fix_v=True)
            rec.update(tau0=tau0, loglik_profile_tau0_warm=fpw.loglik,
                       loglik_profile_tau0_cold=fpc.loglik, profile_converged_warm=fpw.converged,
                       profile_converged_cold=fpc.converged,
                       profile_proj_max_grad_warm=fpw.proj_max_grad)
            for pol, llp in (("warm", fpw.loglik), ("best", better(fpw, fpc).loglik)):
                diff = ll_hat - llp
                r = float(np.sign(tau0 - tau_hat) * np.sqrt(2.0 * max(diff, 0.0)))
                rec.update({f"loglik_diff_{pol}": diff, f"signed_root_{pol}": r,
                            f"p_tau_{pol}": float(norm.sf(r))})
        for pol in ("warm", "best"):
            rec[f"p_source_{pol}"] = max(rec[pol]["p_up"], rec[pol]["p_lo"], rec[f"p_tau_{pol}"])
        out["sources"][s] = rec
    for pol in ("warm", "best"):
        out[f"p_{pol}"] = max(r[f"p_source_{pol}"] for r in out["sources"].values())
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


def lr_fit(D, theta3_start, v_max):
    f3 = ft.fit_section3(D, start=theta3_start)
    best, f40, f41 = ft.fit_section4_two_starts(D, f3.theta, v_max)
    return f3, best, f40, f41, ft.lr_statistic(f3, best)


def run_p3(prob, fit3, LR_obs, U, v_max):
    """LR* per draw with every included row (i, g, tag, n) held fixed; p = (1 + #{LR* >= LR -
    1e-6}) / (B + 1); a failed refit counts as an exceedance."""
    D = prob.data()
    L = prob.L
    draws = []
    for b in range(U.shape[0]):
        a_b, cls = simulate_counts(D, fit3.theta, U[b])
        Db = prob.data(rows=prob.rows.with_counts(a_b))
        f3, best, f40, f41, LR = lr_fit(Db, fit3.theta, v_max)
        ok = f3.converged and best.converged
        rec = dict(draw=b, LR=LR, loglik3=f3.loglik, loglik4=best.loglik,
                   v_hat=float(best.theta[L.v]), kappa3=float(f3.theta[L.kappa]),
                   converged=ok, start_kept="v0" if best is f40 else "vmax",
                   n_class_H=int(np.sum(cls == 0)), n_class_R=int(np.sum(cls == 1)),
                   n_class_A=int(np.sum(cls == 2)), sum_a=int(a_b.sum()),
                   exceed=int((not ok) or LR >= LR_obs - TIE),
                   in_band=bool(abs(LR - LR_obs) <= BAND))
        # best policy: also refit section 3 from the section 7 start; if that optimum is
        # higher, add section-4 fits started from it and recompute LR over all fits
        fc = ft.fit_section3(Db)
        rec.update(loglik3_cold=fc.loglik, kappa3_cold=float(fc.theta[L.kappa]),
                   cold_better=bool(fc.converged and fc.loglik > f3.loglik + BETTER))
        if rec["cold_better"]:
            b4c, _, _ = ft.fit_section4_two_starts(Db, fc.theta, v_max)
            b4 = better(best, b4c)
            LRb = ft.lr_statistic(fc, b4)
            okb = fc.converged and b4.converged
            rec.update(LR_best=LRb, v_hat_best=float(b4.theta[L.v]), converged_best=okb,
                       exceed_best=int((not okb) or LRb >= LR_obs - TIE))
        else:
            rec.update(LR_best=LR, v_hat_best=rec["v_hat"], converged_best=ok,
                       exceed_best=rec["exceed"])
        draws.append(rec)
    B = len(draws)
    p = _p([d["exceed"] for d in draws])
    pb = _p([d["exceed_best"] for d in draws])
    return dict(LR_obs=LR_obs, B=B, p_warm=p, mc_se_warm=mc_se(p, B), p_best=pb,
                mc_se_best=mc_se(pb, B),
                failed_draws=[d["draw"] for d in draws if not d["converged"]],
                band_draws=[d["draw"] for d in draws if d["in_band"]],
                cold_better_draws=[d["draw"] for d in draws if d["cold_better"]], draws=draws)
