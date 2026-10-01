"""Model A-final likelihood (spec v1 with v1.1 and v1.2): oriented truncated genotype mixture,
section-3 likelihood with analytic gradient and Hessian, section-4 adaptive Gauss-Hermite
marginal in v = tau^2, and the boundary score at v = 0.

Coordinates (spec section 7): theta = [alpha_1..alpha_G, r_1..r_G, kappa, delta (design
column order), omega_S, gamma_fib, gamma_mye (kappa_A only), rho_S, rho_b, rho_d, v].
Section 3 ignores v. Units are natural logits.
"""
import numpy as np
import scipy.sparse as sp
from numpy.polynomial.hermite import hermgauss
from scipy.special import expit, logsumexp

from betabinom import (bb_logpmf, digamma_diff, lchoose, lgamma_diff, log_mass_over_K,
                       seg_logsumexp, seg_ranges, trigamma_diff)
from het_rule import truncation_bounds

EPS_RHO = 1e-6          # smooth rho floor (section 3.2)
E_FLOOR = 1e-4          # homozygote error floor (section 2.3)
SQRT2 = np.sqrt(2.0)
LOG_HALF = np.log(0.5)


# ---------------------------------------------------------------------------------------------
# Rows: fixed per-row inputs


class Rows:
    """Included rows (individual x gene, chosen tag): counts, K_c(n), merged priors (section 2.1)
    mapped to the oriented classes (section 2.2), and the fixed homozygote terms."""

    def __init__(self, n, a, rule_d, rule_k, rule_F, orientation, p_ref, p_alt, p_x, e, phi):
        self.n = np.asarray(n, dtype=np.int64)
        self.orientation = np.asarray(orientation, dtype=np.int64)
        self.s_t = -self.orientation            # +1 when the lead-ALT allele is the tag REF
        self.rule = (np.asarray(rule_d), np.asarray(rule_k), np.asarray(rule_F))
        lo, hi, nonempty = truncation_bounds(self.n, *self.rule)
        if not nonempty.all():
            raise ValueError("a row has an empty truncation set")
        self.lo, self.hi = lo, hi
        p_ref, p_alt, p_x = (np.asarray(v, float) for v in (p_ref, p_alt, p_x))
        PH = 2.0 * p_ref * p_alt
        PR = p_ref ** 2 + 2.0 * p_ref * p_x
        PA = p_alt ** 2 + 2.0 * p_alt * p_x
        plus = self.orientation == 1
        self.logPH = np.log(PH)
        self.logPO0 = np.log(np.where(plus, PR, PA))   # homozygous for the lead-REF-side allele
        self.logPO1 = np.log(np.where(plus, PA, PR))   # homozygous for the lead-ALT-side allele
        self.logPR, self.logPA = np.log(PR), np.log(PA)
        self.e = np.maximum(np.asarray(e, float), E_FLOOR)
        self.phi = float(phi)
        self.logC_O0 = log_mass_over_K(self.n, lo, hi, self.e, self.phi)
        self.logC_O1 = log_mass_over_K(self.n, lo, hi, 1.0 - self.e, self.phi)
        self.logcD = np.logaddexp(self.logPO0 + self.logC_O0, self.logPO1 + self.logC_O1)
        self._build_tails()
        self.set_counts(a)

    def _build_tails(self):
        nr = len(self.n)
        rows = np.arange(nr)
        s_lo, a_lo = seg_ranges(np.zeros(nr, np.int64), self.lo)            # a < lo
        s_hi, a_hi = seg_ranges(self.hi + 1, self.n - self.hi)              # a > n - lo
        self.tail_row = np.concatenate([rows[s_lo], rows[s_hi]])
        self.tail_a = np.concatenate([a_lo, a_hi])
        self.tail_lchoose = lchoose(self.n[self.tail_row], self.tail_a)

    def set_counts(self, a):
        a = np.asarray(a, dtype=np.int64)
        if ((a < self.lo) | (a > self.hi)).any():
            raise ValueError("an oriented count lies outside K_c(n)")
        self.a = a
        lf0 = bb_logpmf(a, self.n, self.e, self.phi)
        lf1 = bb_logpmf(a, self.n, 1.0 - self.e, self.phi)
        self.logcN = np.logaddexp(self.logPO0 + lf0, self.logPO1 + lf1)
        self.a_lchoose = lchoose(self.n, a)

    def take(self, idx):
        """Rows restricted to idx (fixed terms sliced, not recomputed)."""
        out = object.__new__(Rows)
        for name in ("n", "orientation", "s_t", "lo", "hi", "logPH", "logPO0", "logPO1", "logPR",
                     "logPA", "e", "logC_O0", "logC_O1", "logcD", "a", "logcN", "a_lchoose"):
            setattr(out, name, getattr(self, name)[idx])
        out.rule = tuple(np.broadcast_to(r, self.n.shape)[idx] for r in self.rule)
        out.phi = self.phi
        out._build_tails()
        return out

    def with_counts(self, a):
        out = self.take(np.arange(len(self.n)))
        out.set_counts(a)
        return out


# ---------------------------------------------------------------------------------------------
# Row likelihood and its derivatives in (x1, x2) = (logit p_ig, dispersion linear predictor)


class _RowState:
    """Per-row BB quantities at (x1, x2); `terms` evaluates log BB(a) and its x-derivatives."""

    def __init__(self, x1, x2, n, order):
        self.n = n
        self.mu, self.omu = expit(x1), expit(-x1)
        sig = expit(x2)
        rho = EPS_RHO + (1.0 - 2.0 * EPS_RHO) * sig
        t = (1.0 - rho) / rho
        self.t = t
        self.al, self.be = self.mu * t, self.omu * t
        self.lgc = -lgamma_diff(t, n)
        if order >= 1:
            self.Dt = digamma_diff(t, n)
            self.m1 = self.mu * self.omu
            rp = (1.0 - 2.0 * EPS_RHO) * sig * (1.0 - sig)      # d rho / d x2
            self.tp = -rp / rho ** 2                            # d t / d x2
        if order >= 2:
            self.Tt = trigamma_diff(t, n)
            rpp = rp * (1.0 - 2.0 * sig)
            self.tpp = -rpp / rho ** 2 + 2.0 * rp ** 2 / rho ** 3

    def terms(self, j, a, lch, order):
        al, be, n = self.al[j], self.be[j], self.n[j]
        lf = lch + lgamma_diff(al, a) + lgamma_diff(be, n - a) + self.lgc[j]
        if order == 0:
            return (lf,)
        mu, omu, t, m1, tp = self.mu[j], self.omu[j], self.t[j], self.m1[j], self.tp[j]
        Da = digamma_diff(al, a)
        Db = digamma_diff(be, n - a)
        Dab = Da - Db
        Dmix = mu * Da + omu * Db - self.Dt[j]
        d1 = m1 * t * Dab
        d2 = tp * Dmix
        if order == 1:
            return lf, d1, d2
        Ta = trigamma_diff(al, a)
        Tb = trigamma_diff(be, n - a)
        d11 = m1 * (omu - mu) * t * Dab + (m1 * t) ** 2 * (Ta + Tb)
        d12 = m1 * tp * (Dab + t * (mu * Ta - omu * Tb))
        d22 = self.tpp[j] * Dmix + tp ** 2 * (mu ** 2 * Ta + omu ** 2 * Tb - self.Tt[j])
        return lf, d1, d2, d11, d12, d22


class RowEval:
    """Result holder: logL per row; (order >= 1) g1, g2 = d log L / d(x1, x2); (order 2) g11,
    g12, g22."""


def log_C_H(st, R, order):
    """log C_H(n) = log sum over K of BB(a | n, p, rho) and its x-derivatives (section 7).

    T = mass of the two excluded tails by log-sum-exp; log C_H = log1p(-T) when T < 0.5,
    else a direct log-sum-exp over K for that row.
    """
    nr = len(R.n)
    tl = st.terms(R.tail_row, R.tail_a, R.tail_lchoose, order)
    logT = seg_logsumexp(tl[0], R.tail_row, nr)
    direct = logT >= LOG_HALF
    C = -np.expm1(logT)                        # 1 - T
    with np.errstate(divide="ignore", invalid="ignore"):
        logC = np.log1p(-np.exp(logT))
    out = [logC]
    if order >= 1:
        f = np.exp(tl[0])
        C = np.where(direct, 1.0, C)           # rows with T >= 0.5 are filled below
        c1 = -np.bincount(R.tail_row, f * tl[1], nr) / C
        c2 = -np.bincount(R.tail_row, f * tl[2], nr) / C
        out += [c1, c2]
        if order >= 2:
            out.append(-np.bincount(R.tail_row, f * (tl[3] + tl[1] ** 2), nr) / C - c1 ** 2)
            out.append(-np.bincount(R.tail_row, f * (tl[4] + tl[1] * tl[2]), nr) / C - c1 * c2)
            out.append(-np.bincount(R.tail_row, f * (tl[5] + tl[2] ** 2), nr) / C - c2 ** 2)
    if direct.any():
        rd = np.flatnonzero(direct)
        seg, a = seg_ranges(R.lo[rd], R.hi[rd] - R.lo[rd] + 1)
        j = rd[seg]
        kt = st.terms(j, a, lchoose(R.n[j], a), order)
        lC = seg_logsumexp(kt[0], seg, len(rd))
        out[0][rd] = lC
        if order >= 1:
            pi = np.exp(kt[0] - lC[seg])
            e1 = np.bincount(seg, pi * kt[1], len(rd))
            e2 = np.bincount(seg, pi * kt[2], len(rd))
            out[1][rd], out[2][rd] = e1, e2
            if order >= 2:
                out[3][rd] = np.bincount(seg, pi * (kt[3] + kt[1] ** 2), len(rd)) - e1 ** 2
                out[4][rd] = np.bincount(seg, pi * (kt[4] + kt[1] * kt[2]), len(rd)) - e1 * e2
                out[5][rd] = np.bincount(seg, pi * (kt[5] + kt[2] ** 2), len(rd)) - e2 ** 2
    return out


def row_eval(x1, x2, R, order):
    """Row contribution L_ig = sum_G P(G) f_G(a|n) / sum_G P(G) C_G(n) (section 2.3).

    Only the H class depends on (x1, x2); the homozygote numerator and normalizers are fixed.
    """
    st = _RowState(x1, x2, R.n, order)
    nr = len(R.n)
    num = st.terms(np.arange(nr), R.a, R.a_lchoose, order)
    ch = log_C_H(st, R, order)
    logN = np.logaddexp(R.logPH + num[0], R.logcN)
    logD = np.logaddexp(R.logPH + ch[0], R.logcD)
    out = RowEval()
    out.logL = logN - logD
    out.logC_H = ch[0]
    if order >= 1:
        wN, vN = np.exp(R.logPH + num[0] - logN), np.exp(R.logcN - logN)      # vN = 1 - wN
        wD, vD = np.exp(R.logPH + ch[0] - logD), np.exp(R.logcD - logD)
        out.g1 = wN * num[1] - wD * ch[1]
        out.g2 = wN * num[2] - wD * ch[2]
        if order >= 2:
            out.g11 = wN * num[3] + wN * vN * num[1] ** 2 - wD * ch[3] - wD * vD * ch[1] ** 2
            out.g12 = (wN * num[4] + wN * vN * num[1] * num[2]
                       - wD * ch[4] - wD * vD * ch[1] * ch[2])
            out.g22 = wN * num[5] + wN * vN * num[2] ** 2 - wD * ch[5] - wD * vD * ch[2] ** 2
    return out


# ---------------------------------------------------------------------------------------------
# Parameter layout and design


class Layout:
    """Section 7 coordinate order."""

    def __init__(self, gene_names, design_columns, kappa_A=False):
        G, p = len(gene_names), len(design_columns)
        self.G, self.p, self.kappa_A = G, p, kappa_A
        self.gene_names, self.design_columns = list(gene_names), list(design_columns)
        self.alpha = np.arange(0, G)
        self.r = np.arange(G, 2 * G)
        self.kappa = 2 * G
        self.delta = np.arange(2 * G + 1, 2 * G + 1 + p)
        i = 2 * G + 1 + p
        self.omega_S = i
        i += 1
        if kappa_A:
            self.gamma_fib, self.gamma_mye = i, i + 1
            i += 2
        self.rho_S, self.rho_b, self.rho_d, self.v = i, i + 1, i + 2, i + 3
        self.size = i + 4
        names = ([f"alpha[{g}]" for g in gene_names] + [f"r[{g}]" for g in gene_names]
                 + ["kappa"] + [f"delta[{c}]" for c in design_columns] + ["omega_S"]
                 + (["gamma_fib", "gamma_mye"] if kappa_A else [])
                 + ["rho_S", "rho_b", "rho_d", "v"])
        self.names = names
        lb = np.full(self.size, -np.inf)
        ub = np.full(self.size, np.inf)
        lb[self.alpha], ub[self.alpha] = -6.0, 6.0
        lb[self.r], ub[self.r] = -15.0, 5.0
        lb[self.kappa], ub[self.kappa] = -0.9, 0.9
        for k in (self.rho_S, self.rho_b, self.rho_d, self.omega_S):
            lb[k], ub[k] = -3.0, 3.0
        lb[self.v] = 0.0
        self.lb, self.ub = lb, ub

    def from_dict(self, d):
        th = np.zeros(self.size)
        th[self.alpha] = d["alpha"]
        th[self.r] = d["r"]
        th[self.kappa] = d["kappa"]
        th[self.delta] = d["delta"]
        th[self.omega_S] = d["omega_S"]
        if self.kappa_A:
            th[self.gamma_fib], th[self.gamma_mye] = d.get("gamma_fib", 0.0), d.get("gamma_mye", 0.0)
        th[self.rho_S], th[self.rho_b], th[self.rho_d] = d["rho_S"], d["rho_b"], d["rho_d"]
        th[self.v] = d.get("v", 0.0)
        return th

    def to_dict(self, th):
        return dict(zip(self.names, (float(x) for x in th)))


def design_spec(cohort, zraw, base_names):
    """Fixed z columns (section 3.3; v1.2 item 5): retained base columns in order, then one
    missing indicator per base column that is partly missing within some cohort, in base order.
    A column that is 0 in every cohort after centring is dropped."""
    spec = [("value", j, base_names[j]) for j in range(len(base_names))]
    for j in range(len(base_names)):
        miss = np.isnan(zraw[:, j])
        if any(0 < miss[cohort == c].sum() < (cohort == c).sum() for c in np.unique(cohort)):
            spec.append(("missing", j, base_names[j] + "_missing"))
    Z = centre_z(cohort, zraw, np.ones(len(cohort)), spec)
    return [s for s, col in zip(spec, Z.T) if np.any(col != 0.0)]


def centre_z(cohort, zraw, mult, spec):
    """z columns centred within cohort with multiplicities; missing values set to the cohort
    mean of the observed values before centring; all-missing -> 0 in that cohort."""
    Z = np.zeros((len(cohort), len(spec)))
    for c in np.unique(cohort):
        sel = (cohort == c) & (mult > 0)
        if not sel.any():
            continue
        m = mult[sel]
        for k, (kind, j, _) in enumerate(spec):
            x = zraw[sel, j]
            miss = np.isnan(x)
            if kind == "missing":
                ind = miss.astype(float)
                Z[sel, k] = ind - np.sum(m * ind) / m.sum()
            elif (~miss).any():
                fill = np.sum(m[~miss] * x[~miss]) / m[~miss].sum()
                xf = np.where(miss, fill, x)
                Z[sel, k] = xf - np.sum(m * xf) / m.sum()
    return Z


def centre_stage(cohort, S, mult):
    St = np.zeros(len(S))
    for c in np.unique(cohort):
        sel = (cohort == c) & (mult > 0)
        if sel.any():
            St[sel] = S[sel] - np.sum(mult[sel] * S[sel]) / mult[sel].sum()
    return St


class Problem:
    """Raw inputs of one fit set: individuals (cohort, S, raw z, cluster), included rows,
    genes and fixed cohort inputs. `data()` builds the centred Data for a fit."""

    def __init__(self, ind_ids, ind_cohort, S, zraw, base_names, cluster, rows, row_ind,
                 row_gene, gene_names, omega_c, w_fib=None, w_mye=None, kappa_A=False,
                 b_name="z_b_log_e", d_name="z_d_dup"):
        self.ind_ids = list(ind_ids)
        self.cohort = np.asarray(ind_cohort)
        self.S = np.asarray(S, float)
        self.zraw = np.asarray(zraw, float)
        self.base_names = list(base_names)
        self.cluster = np.asarray(cluster)
        self.rows = rows
        self.row_ind = np.asarray(row_ind, dtype=np.int64)
        self.row_gene = np.asarray(row_gene, dtype=np.int64)
        self.omega_row = np.asarray(omega_c, float)          # omega_c of the row's cohort
        self.spec = design_spec(self.cohort, self.zraw, self.base_names)
        self.design_columns = [s[2] for s in self.spec]
        self.b_col = self.design_columns.index(b_name)
        self.d_col = self.design_columns.index(d_name)
        G = len(gene_names)
        self.w_fib = np.zeros(G) if w_fib is None else np.asarray(w_fib, float)
        self.w_mye = np.zeros(G) if w_mye is None else np.asarray(w_mye, float)
        self.L = Layout(gene_names, self.design_columns, kappa_A)

    def data(self, S=None, mult=None, rows=None):
        S = self.S if S is None else np.asarray(S, float)
        mult = np.ones(len(self.ind_ids)) if mult is None else np.asarray(mult, float)
        rows = self.rows if rows is None else rows
        St = centre_stage(self.cohort, S, mult)
        Z = centre_z(self.cohort, self.zraw, mult, self.spec)
        w_row = mult[self.row_ind]
        keep = np.flatnonzero(w_row > 0)
        return Data(self, rows.take(keep) if len(keep) < len(w_row) else rows, keep, St, Z, mult)


class Data:
    """Centred per-row design, row weights (multiplicities) and the free-parameter masks."""

    def __init__(self, prob, rows, keep, St_ind, Z_ind, mult):
        L = prob.L
        self.L, self.rows = L, rows
        ind = prob.row_ind[keep]
        self.ind = ind
        self.gene = prob.row_gene[keep]
        self.w = mult[ind]
        self.St = St_ind[ind]
        self.Z = Z_ind[ind]
        self.b = self.Z[:, prob.b_col]
        self.d = self.Z[:, prob.d_col]
        self.s_t = rows.s_t.astype(float)
        self.omega = prob.omega_row[keep]
        self.wfib = prob.w_fib[self.gene]
        self.wmye = prob.w_mye[self.gene]
        self.nr = len(ind)
        self.cluster_codes, self.row_cluster = np.unique(prob.cluster[ind], return_inverse=True)
        self.gene_has_rows = np.bincount(self.gene, minlength=L.G) > 0
        free = np.ones(L.size, bool)
        free[L.alpha[~self.gene_has_rows]] = False
        free[L.r[~self.gene_has_rows]] = False
        zero_col = ~np.any(self.Z != 0.0, axis=0)
        free[L.delta[zero_col]] = False
        free[L.v] = False
        self.free3 = free
        self.free4 = free.copy()
        self.free4[L.v] = True
        # J2 does not depend on theta
        r = np.arange(self.nr)
        rows_ = np.concatenate([r, r, r, r])
        cols = np.concatenate([L.r[self.gene], np.full(self.nr, L.rho_S), np.full(self.nr, L.rho_b),
                               np.full(self.nr, L.rho_d)])
        vals = np.concatenate([np.ones(self.nr), self.St, self.b, self.d])
        self.J2 = sp.csr_matrix((vals, (rows_, cols)), shape=(self.nr, L.size))
        self.Agg = sp.csr_matrix((np.ones(self.nr), (self.gene, r)), shape=(L.G, self.nr))


def predictors(theta, D):
    """x1 = logit p_ig (section 3 form, u = 0) and x2 = dispersion linear predictor."""
    L = D.L
    al = theta[L.alpha][D.gene]
    q = 1.0 + theta[L.kappa] * D.St + D.Z @ theta[L.delta]
    if L.kappa_A:
        q = q + (theta[L.gamma_fib] * D.wfib + theta[L.gamma_mye] * D.wmye) * D.St
    x1 = al * q + D.s_t * (D.omega + theta[L.omega_S] * D.St)
    x2 = theta[L.r][D.gene] + theta[L.rho_S] * D.St + theta[L.rho_b] * D.b + theta[L.rho_d] * D.d
    return x1, x2, al, q


def jacobian1(D, al, q):
    """d x1 / d theta per row (sparse, rows x P)."""
    L, nr = D.L, D.nr
    r = np.arange(nr)
    rows_ = [r, r, np.repeat(r, L.p), r]
    cols = [L.alpha[D.gene], np.full(nr, L.kappa), np.tile(L.delta, nr), np.full(nr, L.omega_S)]
    vals = [q, al * D.St, (al[:, None] * D.Z).ravel(), D.s_t * D.St]
    if L.kappa_A:
        rows_ += [r, r]
        cols += [np.full(nr, L.gamma_fib), np.full(nr, L.gamma_mye)]
        vals += [al * D.wfib * D.St, al * D.wmye * D.St]
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows_), np.concatenate(cols))),
                         shape=(nr, L.size))


def _second_order_x1(D, G1):
    """sum over rows of G1 * d^2 x1 / d theta^2 (alpha_g crossed with kappa, delta, gamma)."""
    L = D.L
    H = np.zeros((L.size, L.size))
    shared = [(L.kappa, D.St)] + [(L.delta[k], D.Z[:, k]) for k in range(L.p)]
    if L.kappa_A:
        shared += [(L.gamma_fib, D.wfib * D.St), (L.gamma_mye, D.wmye * D.St)]
    for col, x in shared:
        s = np.bincount(D.gene, G1 * x, L.G)
        H[L.alpha, col] += s
        H[col, L.alpha] += s
    return H


def _assemble(D, J1, G1, G2, G11=None, G12=None, G22=None):
    grad = J1.T @ G1 + D.J2.T @ G2
    if G11 is None:
        return grad, None
    J2 = D.J2
    H = (J1.T @ sp.diags(G11) @ J1 + J1.T @ sp.diags(G12) @ J2 + J2.T @ sp.diags(G12) @ J1
         + J2.T @ sp.diags(G22) @ J2).toarray()
    H += _second_order_x1(D, G1)
    return grad, H


def loglik3(theta, D, order=1):
    """Section-3 log-likelihood sum_i w_i log L_ig; order 1 adds the gradient, order 2 the
    Hessian (both in the full section 7 coordinates; the v entries are 0)."""
    x1, x2, al, q = predictors(theta, D)
    re = row_eval(x1, x2, D.rows, order)
    ll = float(np.sum(D.w * re.logL))
    if order == 0:
        return ll
    J1 = jacobian1(D, al, q)
    w = D.w
    if order == 1:
        grad, _ = _assemble(D, J1, w * re.g1, w * re.g2)
        return ll, grad
    grad, H = _assemble(D, J1, w * re.g1, w * re.g2, w * re.g11, w * re.g12, w * re.g22)
    return ll, grad, H


def row_scores(theta, D):
    """Per-row score vectors (sparse rows x P) at theta, section 3."""
    x1, x2, al, q = predictors(theta, D)
    re = row_eval(x1, x2, D.rows, 1)
    J1 = jacobian1(D, al, q)
    return sp.diags(D.w * re.g1) @ J1 + sp.diags(D.w * re.g2) @ D.J2


def boundary_score(theta, D):
    """d ell / d v at v = 0: 1/2 sum_g [(sum_i d_u log L_ig)^2 + sum_i d2_u log L_ig] at u = 0."""
    x1, x2, _, _ = predictors(theta, D)
    re = row_eval(x1, x2, D.rows, 2)
    G = D.L.G
    s1 = np.bincount(D.gene, D.w * re.g1 * D.St, G)
    s2 = np.bincount(D.gene, D.w * re.g11 * D.St ** 2, G)
    return 0.5 * float(np.sum(s1 ** 2 + s2))


def gene_loglik3(theta, D):
    """Per-gene section-3 log-likelihood (used for the section 7 starting values)."""
    x1, x2, _, _ = predictors(theta, D)
    re = row_eval(x1, x2, D.rows, 0)
    return np.bincount(D.gene, D.w * re.logL, D.L.G)


# ---------------------------------------------------------------------------------------------
# Section 4: adaptive Gauss-Hermite in v = tau^2


def gene_modes(x1b, x2, D, v, u0=None, tol=1e-10, max_iter=200):
    """u_hat_g = mode of h(u) = sum_i w_i log L_ig(u) - u^2/(2v) by Newton to |du| < tol, and
    sigma_hat_g = (-h''(u_hat_g))^(-1/2)."""
    G, gene, St, w = D.L.G, D.gene, D.St, D.w
    u = np.zeros(G) if u0 is None else np.array(u0, float)
    cap = 10.0 * np.sqrt(v)
    for _ in range(max_iter):
        re = row_eval(x1b + u[gene] * St, x2, D.rows, 2)
        hp = np.bincount(gene, w * re.g1 * St, G) - u / v
        hpp = np.bincount(gene, w * re.g11 * St ** 2, G) - 1.0 / v
        hpp = np.where(hpp < 0, hpp, -1.0 / v)
        step = np.clip(-hp / hpp, -cap, cap)
        u = u + step
        if np.max(np.abs(step)) < tol:
            break
    else:
        raise RuntimeError("mode search for u_g did not converge")
    re = row_eval(x1b + u[gene] * St, x2, D.rows, 2)
    hpp = np.bincount(gene, w * re.g11 * St ** 2, G) - 1.0 / v
    if np.any(hpp >= 0):
        raise RuntimeError("h''(u_hat) >= 0")
    return u, 1.0 / np.sqrt(-hpp)


def loglik4(theta, D, order=1, nodes=20, detail=False):
    """Section-4 log-likelihood sum_g log int N(u; 0, v) prod_i L_ig(u)^{w_i} du by adaptive
    Gauss-Hermite with `nodes` nodes. v = 0 is the separate branch (section 3, with d/dv the
    boundary score). Gradients and Hessian hold (u_hat, sigma_hat) fixed."""
    L = D.L
    v = float(theta[L.v])
    if v == 0.0:
        res = loglik3(theta, D, order)
        if order == 0:
            return res
        res = list(res)
        res[1] = res[1].copy()
        res[1][L.v] = boundary_score(theta, D)
        if order == 2:
            res[2] = res[2].copy()
            res[2][L.v, :] = np.nan
            res[2][:, L.v] = np.nan
        return tuple(res)
    x1b, x2, al, q = predictors(theta, D)
    G, gene, St, w = L.G, D.gene, D.St, D.w
    u_hat, sig = gene_modes(x1b, x2, D, v)
    xk, wk = hermgauss(nodes)
    logW = np.log(wk) + xk ** 2
    U = u_hat[None, :] + SQRT2 * sig[None, :] * xk[:, None]        # nodes x G
    Hk = np.empty((nodes, G))
    ev = []
    for k in range(nodes):
        re = row_eval(x1b + U[k][gene] * St, x2, D.rows, order)
        Hk[k] = np.bincount(gene, w * re.logL, G) - U[k] ** 2 / (2.0 * v)
        ev.append(re)
    A = logW[:, None] + Hk
    lse = logsumexp(A, axis=0)
    has = D.gene_has_rows
    ell_g = lse + np.log(SQRT2 * sig) - 0.5 * np.log(2.0 * np.pi * v)
    ll = float(np.sum(ell_g[has]))
    if order == 0:
        return (ll, dict(ell_g=ell_g, u_hat=u_hat, sigma=sig)) if detail else ll
    pi = np.exp(A - lse[None, :])                                  # nodes x G
    P = pi[:, gene]                                                # nodes x rows
    E1 = sum(P[k] * ev[k].g1 for k in range(nodes))                # posterior means per row
    E2 = sum(P[k] * ev[k].g2 for k in range(nodes))
    J1 = jacobian1(D, al, q)
    q_v = U ** 2 / (2.0 * v ** 2)                                  # d h_k / d v
    Eq_v = np.sum(pi * q_v, axis=0)
    dv = float(np.sum((Eq_v - 1.0 / (2.0 * v))[has]))
    if order == 1:
        grad, _ = _assemble(D, J1, w * E1, w * E2)
        grad[L.v] = dv
        return ll, grad
    E11 = sum(P[k] * ev[k].g11 for k in range(nodes))
    E12 = sum(P[k] * ev[k].g12 for k in range(nodes))
    E22 = sum(P[k] * ev[k].g22 for k in range(nodes))
    grad, H = _assemble(D, J1, w * E1, w * E2, w * E11, w * E12, w * E22)
    grad[L.v] = dv
    # between-node covariance of d h_k / d(theta, v) within each gene
    gidx = np.arange(G)
    blocks = []
    for k in range(nodes):
        sq = np.sqrt(P[k])
        Mk = D.Agg @ (sp.diags(w * sq * (ev[k].g1 - E1)) @ J1
                      + sp.diags(w * sq * (ev[k].g2 - E2)) @ D.J2)
        vcol = np.sqrt(pi[k]) * (q_v[k] - Eq_v) * has
        Mk = Mk + sp.csr_matrix((vcol, (gidx, np.full(G, L.v))), shape=(G, L.size))
        blocks.append(Mk)
    M = sp.vstack(blocks).tocsr()
    H += (M.T @ M).toarray()
    H[L.v, L.v] += float(np.sum((np.sum(pi * (-U ** 2 / v ** 3), axis=0)
                                 + 1.0 / (2.0 * v ** 2))[has]))
    return ll, grad, H


def gh_gene_integral_quad(theta, D, g, v, width=40.0):
    """Reference ell_g by scipy.integrate.quad (epsabs 1e-12) for the accuracy test."""
    from scipy.integrate import quad
    x1b, x2, _, _ = predictors(theta, D)
    sel = D.gene == g
    R = D.rows.take(np.flatnonzero(sel))
    x1g, x2g, Sg, wg = x1b[sel], x2[sel], D.St[sel], D.w[sel]

    def h(u):
        return float(np.sum(wg * row_eval(x1g + u * Sg, x2g, R, 0).logL)) - u ** 2 / (2.0 * v)

    u_hat, sig = gene_modes(x1b, x2, D, v)
    c, s = u_hat[g], sig[g]
    h0 = h(c)
    val, err = quad(lambda u: np.exp(h(u) - h0), c - width * s, c + width * s,
                    epsabs=1e-12, epsrel=1e-13, limit=500, points=[c])
    return h0 + np.log(val) - 0.5 * np.log(2.0 * np.pi * v), err
