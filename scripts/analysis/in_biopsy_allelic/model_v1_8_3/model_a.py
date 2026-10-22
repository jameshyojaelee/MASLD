"""Model A-final likelihood (spec v1 to v1.8): oriented truncated genotype mixture,
section-3 likelihood with analytic gradient and Hessian, section-4 adaptive Gauss-Hermite
marginal in v = tau^2 with the v1.8 section 1 mode search, global check and reference
integrator, v in [0, 1], and the boundary score at v = 0.

Coordinates (spec section 7): theta = [alpha_1..alpha_G, r_1..r_G, kappa, delta (design
column order), omega_S, gamma_fib, gamma_mye (kappa_A only), rho_S, rho_b, rho_d, v].
Section 3 ignores v. Units are natural logits.
"""
import numpy as np
import scipy.sparse as sp
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
from scipy.special import expit, logsumexp

from betabinom import (bb_logpmf, digamma_diff, digamma_diff_rows, lchoose, lgamma_diff,
                       lgamma_diff_rows, log_mass_over_K, row_halves, seg_logsumexp, seg_ranges,
                       trigamma_diff, trigamma_diff_rows)
from het_rule import truncation_bounds

EPS_RHO = 1e-6          # smooth rho floor (section 3.2)
E_FLOOR = 1e-4          # homozygote error floor (section 2.3)
SQRT2 = np.sqrt(2.0)
LOG_HALF = np.log(0.5)

# v1.8 section 1 (mode search, global check, reference integrator, v bound)
V_UPPER = 1.0           # item 10: v in [0, 1]
GOLDEN = 0.3819660112501051
MODE_DU = 1e-10         # item 3 stop: |du| < 1e-10 ...
MODE_DH = 1e-8          # ... and |h'(c)| < 1e-8 max(1, |h(c)|)
MODE_MAX_ITER = 500     # implementation guard; a search that hits it is treated as not converged
SAME_MODE = 1e-8        # item 4 (i): modes more than 1e-8 max(1, |u|) apart
GRID_MIN, GRID_CAP = 129, 8193
FLAT_H2 = -1e-12        # items 4 (iii) and 5
SECOND_MODE_GAP = 25.0  # item 5
H_TIE = 1e-10           # item 4: ties in h
H_ROUND = 1e-12         # rounding allowance in the item 3 comparison (implementation, item 9)
REF_TAIL = 37.0         # item 6: W
REF_AGREE = 1e-10       # item 6: doubling stops when two successive values agree within 1e-10
REF_MAX_PANELS = 8192
GH_REF_TOL = 1e-8       # item 7
CHUNK_ROWS = 200_000    # rows per replicated block in per-gene evaluations (memory only)
GL_X, GL_W = leggauss(20)
GH, REF = 0, 1          # integrator codes


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
        # v1.8.1 (P-E): x-only halves of the gamma-function differences, once per row
        self.hal, self.hbe = row_halves(self.al, order), row_halves(self.be, order)
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
        n = self.n[j]
        nma = n - a
        lf = (lch + lgamma_diff_rows(self.al, self.hal, j, a)
              + lgamma_diff_rows(self.be, self.hbe, j, nma) + self.lgc[j])
        if order == 0:
            return (lf,)
        mu, omu, t, m1, tp = self.mu[j], self.omu[j], self.t[j], self.m1[j], self.tp[j]
        Da = digamma_diff_rows(self.al, self.hal, j, a)
        Db = digamma_diff_rows(self.be, self.hbe, j, nma)
        Dab = Da - Db
        Dmix = mu * Da + omu * Db - self.Dt[j]
        d1 = m1 * t * Dab
        d2 = tp * Dmix
        if order == 1:
            return lf, d1, d2
        Ta = trigamma_diff_rows(self.al, self.hal, j, a)
        Tb = trigamma_diff_rows(self.be, self.hbe, j, nma)
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
        lb[self.v], ub[self.v] = 0.0, V_UPPER        # v1.8 section 1 item 10
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
        self.n_genes_with_rows = int(self.gene_has_rows.sum())
        self._gidx = None
        self._grows = {}

    def gene_idx(self, g):
        """Row indices of gene g (built once per Data)."""
        if self._gidx is None:
            order = np.argsort(self.gene, kind="stable")
            bounds = np.searchsorted(self.gene[order], np.arange(self.L.G + 1))
            self._gidx = [order[bounds[k]:bounds[k + 1]] for k in range(self.L.G)]
        return self._gidx[g]

    def gene_rows(self, g):
        if g not in self._grows:
            self._grows[g] = self.rows.take(self.gene_idx(g))
        return self._grows[g]


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
# Section 4: v = tau^2. Adaptive Gauss-Hermite 20 (spec section 4) with the v1.8 section 1 mode
# search (item 3), global check (item 4), integrator rule (items 5 and 7) and composite
# Gauss-Legendre reference integrator (item 6).


class Sec4State:
    """Carried between the evaluations of one fit: the previous mode per gene (item 3), whether
    the gene was already evaluated in this fit (item 4 (iv)), the integrator chosen at the last
    evaluation (item 8) and event counts."""

    def __init__(self, G):
        self.u_prev = np.full(G, np.nan)
        self.seen = np.zeros(G, bool)
        self.last_choice = None
        self.events = dict(evaluations=0, global_checks=0, grid_cap=0, second_mode=0,
                           flat_mode=0, starts_disagree=0, nonconcave_iterate=0,
                           start_not_converged=0, reference_checks=0, gh_ref_mismatch=0,
                           reference_used=0, max_reference_panels=0, max_grid_points=0)


class _Ctx:
    """h_g(u) = sum_i w_i log L_ig(u) - u^2/(2v) at fixed (theta, v): all genes at one u each,
    or one gene at many u (rows replicated node-major in blocks of about CHUNK_ROWS rows)."""

    def __init__(self, D, x1b, x2, v):
        self.D, self.x1b, self.x2, self.v = D, x1b, x2, v

    def h_all(self, u, order):
        D, G = self.D, self.D.L.G
        re = row_eval(self.x1b + u[D.gene] * D.St, self.x2, D.rows, order)
        h = np.bincount(D.gene, D.w * re.logL, G) - u ** 2 / (2.0 * self.v)
        if order < 2:
            return h, None, None
        h1 = np.bincount(D.gene, D.w * re.g1 * D.St, G) - u / self.v
        h2 = np.bincount(D.gene, D.w * re.g11 * D.St ** 2, G) - 1.0 / self.v
        return h, h1, h2

    def h_active(self, u, order, genes_active):
        """v1.8.1 (P-A): h_all evaluated on the rows of the active genes only. Entries of the
        active genes equal h_all's (same rows, same order within each gene); entries of the
        other genes are not evaluated (h = -u^2/(2v), h' = -u/v, h'' = -1/v) and must not be read."""
        genes_active = np.asarray(genes_active, bool)
        if genes_active.all():
            return self.h_all(u, order)
        D, G = self.D, self.D.L.G
        sel = np.flatnonzero(genes_active[D.gene])
        gs, Ss, ws = D.gene[sel], D.St[sel], D.w[sel]
        re = row_eval(self.x1b[sel] + u[gs] * Ss, self.x2[sel], D.rows.take(sel), order)
        h = np.bincount(gs, ws * re.logL, G) - u ** 2 / (2.0 * self.v)
        if order < 2:
            return h, None, None
        h1 = np.bincount(gs, ws * re.g1 * Ss, G) - u / self.v
        h2 = np.bincount(gs, ws * re.g11 * Ss ** 2, G) - 1.0 / self.v
        return h, h1, h2

    def gene_blocks(self, g, u, order):
        D = self.D
        idx = D.gene_idx(g)
        Rg = D.gene_rows(g)
        m = len(idx)
        kc = max(1, CHUNK_ROWS // max(m, 1))
        x1g, x2g, Sg = self.x1b[idx], self.x2[idx], D.St[idx]
        for s in range(0, len(u), kc):
            uu = u[s:s + kc]
            k = len(uu)
            Rr = Rg if k == 1 else Rg.take(np.tile(np.arange(m), k))
            x1 = np.tile(x1g, k) + np.repeat(uu, m) * np.tile(Sg, k)
            yield slice(s, s + k), uu, row_eval(x1, np.tile(x2g, k), Rr, order), m

    def h_gene(self, g, u, order):
        u = np.atleast_1d(np.asarray(u, float))
        idx = self.D.gene_idx(g)
        wg, Sg = self.D.w[idx], self.D.St[idx]
        h, h1, h2 = np.empty(len(u)), np.empty(len(u)), np.empty(len(u))
        for sl, uu, re, m in self.gene_blocks(g, u, order):
            k = len(uu)
            h[sl] = (wg * re.logL.reshape(k, m)).sum(1) - uu ** 2 / (2.0 * self.v)
            if order >= 2:
                h1[sl] = (wg * Sg * re.g1.reshape(k, m)).sum(1) - uu / self.v
                h2[sl] = (wg * Sg ** 2 * re.g11.reshape(k, m)).sum(1) - 1.0 / self.v
        return (h, None, None) if order < 2 else (h, h1, h2)


def _local_search(f, a, c, b, ha, hc, hb, h1c, h2c, active=None, subset=False):
    """v1.8 section 1 item 3, run for several independent searches at once. Keeps a < c < b with
    h(c) >= max(h(a), h(b)); a Newton step from c when h''(c) < 0 and the step lands inside
    (a, b), otherwise a golden-section step into the larger subinterval; the point with the
    larger h becomes c. Stops when |du| < 1e-10 and |h'(c)| < 1e-8 max(1, |h(c)|).
    f(x) returns (h, h', h'') at one point per search; with subset=True f(x, act) is called and
    only the entries of the searches still active are read (v1.8.1, P-A). Returns c, h(c),
    h''(c), converged, and whether any iterate had h'' >= 0."""
    a, c, b, ha, hc, hb, h1c, h2c = (np.array(x, float) for x in (a, c, b, ha, hc, hb, h1c, h2c))
    K = len(c)
    done = np.zeros(K, bool) if active is None else ~np.asarray(active, bool)
    ok = np.zeros(K, bool)
    nonconcave = (h2c >= 0) & ~done
    for _ in range(MODE_MAX_ITER):
        act = ~done
        if not act.any():
            break
        with np.errstate(divide="ignore", invalid="ignore"):
            xn = c - h1c / h2c
        newton = (h2c < 0) & (xn > a) & (xn < b)
        xg = np.where((b - c) > (c - a), c + GOLDEN * (b - c), c - GOLDEN * (c - a))
        x = np.where(act, np.where(newton, xn, xg), c)
        hx, h1x, h2x = f(x, act) if subset else f(x)
        nonconcave |= act & (h2x >= 0)
        du = np.abs(x - c)
        # rounding allowance (v1.8 item 9): h carries ~1e-12 relative rounding, so a Newton
        # step that lowers |h'| is accepted when h drops by no more than H_ROUND max(1, |h|)
        up = act & ((hx > hc) | (newton & (np.abs(h1x) < np.abs(h1c))
                                 & (hx >= hc - H_ROUND * np.maximum(1.0, np.abs(hc)))))
        down = act & ~up
        left = x < c
        a_new = np.where(up & ~left, c, np.where(down & left, x, a))
        b_new = np.where(up & left, c, np.where(down & ~left, x, b))
        ha = np.where(up & ~left, hc, np.where(down & left, hx, ha))
        hb = np.where(up & left, hc, np.where(down & ~left, hx, hb))
        a, b = a_new, b_new
        c = np.where(up, x, c)
        hc, h1c, h2c = [np.where(up, new, old) for new, old in ((hx, hc), (h1x, h1c), (h2x, h2c))]
        conv = act & (du < MODE_DU) & (np.abs(h1c) < MODE_DH * np.maximum(1.0, np.abs(hc)))
        stall = act & ~conv & ((b - a) <= 4.0 * np.finfo(float).eps * np.maximum(1.0, np.abs(c)))
        ok |= conv
        done |= conv | stall
    return c, hc, h2c, ok, nonconcave


def _distinct(cands):
    """Local maxima (u, h, h'') merged when within 1e-8 max(1, |u|) of a higher one."""
    out = []
    for cnd in sorted(cands, key=lambda t: -t[1]):
        if all(abs(cnd[0] - o[0]) > SAME_MODE * max(1.0, abs(o[0])) for o in out):
            out.append(cnd)
    return out


def _pick(cands):
    """Item 4: the local maximum with the largest h; ties within 1e-10 in h go to the smaller
    |u|, and a tie in |u| to the negative u."""
    hmax = max(t[1] for t in cands)
    return min((t for t in cands if t[1] >= hmax - H_TIE), key=lambda t: (abs(t[0]), t[0]))


def _sigma(h2):
    return (-h2) ** -0.5 if h2 < 0 else np.nan


def _global_check(ctx, g, U, cands, state):
    """Item 4: scan h on [-U, U] with N = max(129, ceil(4U/sigma_hat) + 1) points (cap 8,193)
    and run the local search from every grid point whose h exceeds both neighbours."""
    best = _pick(cands) if cands else None
    sig = _sigma(best[2]) if best is not None else np.nan
    cap = not np.isfinite(sig)
    N = GRID_CAP if cap else max(GRID_MIN, int(np.ceil(4.0 * U / sig)) + 1)
    if N > GRID_CAP:
        N, cap = GRID_CAP, True
    state.events["max_grid_points"] = max(state.events["max_grid_points"], N)
    grid = np.linspace(-U, U, N)
    hg = ctx.h_gene(g, grid, 0)[0]
    k = np.flatnonzero((hg[1:-1] > hg[:-2]) & (hg[1:-1] > hg[2:])) + 1
    if len(k):
        hc, h1c, h2c = ctx.h_gene(g, grid[k], 2)
        uu, hh, hh2, ok, _ = _local_search(lambda x: ctx.h_gene(g, x, 2), grid[k - 1], grid[k],
                                           grid[k + 1], hg[k - 1], hc, hg[k + 1], h1c, h2c)
        cands = cands + [(float(x), float(y), float(z)) for x, y, z, o in zip(uu, hh, hh2, ok) if o]
    return _distinct(cands), cap


def _modes(ctx, state, check_all):
    """Modes for every gene (item 3, two starts), the global check where item 4 triggers, and
    the item 5 integrator choice (GH or REF; item 7 is applied in loglik4)."""
    D, v = ctx.D, ctx.v
    G = D.L.G
    has = D.gene_has_rows
    zero = np.zeros(G)
    h0, h10, h20 = ctx.h_all(zero, 2)
    U = np.sqrt(2.0 * v * np.maximum(-h0, 0.0))
    hUp, hUm = ctx.h_all(U, 0)[0], ctx.h_all(-U, 0)[0]
    f = lambda x, act: ctx.h_active(x, 2, act)  # noqa: E731  (v1.8.1, P-A)
    u1, hu1, h2u1, ok1, nc1 = _local_search(f, -U, zero, U, hUm, h0, hUp, h10, h20, active=has,
                                            subset=True)
    prev = state.u_prev
    use2 = has & np.isfinite(prev) & (prev > -U) & (prev < U)
    ok2 = np.zeros(G, bool)
    nc2 = np.zeros(G, bool)
    if use2.any():
        p0 = np.where(use2, prev, 0.0)
        hp, h1p, h2p = ctx.h_all(p0, 2)
        use2 &= hp >= np.maximum(hUm, hUp)
        if use2.any():
            u2, hu2, h2u2, ok2, nc2 = _local_search(f, -U, p0, U, hUm, hp, hUp, h1p, h2p,
                                                    active=use2, subset=True)
            ok2 &= use2
    ev = state.events
    out = dict(U=U, u_hat=np.zeros(G), h_hat=np.zeros(G), h2_hat=np.full(G, -1.0 / v),
               sigma=np.full(G, np.sqrt(v)), choice=np.full(G, GH), checked=np.zeros(G, bool),
               modes=[[] for _ in range(G)])
    for g in np.flatnonzero(has):
        cands = [(float(u1[g]), float(hu1[g]), float(h2u1[g]))] if ok1[g] else []
        if ok2[g]:
            cands.append((float(u2[g]), float(hu2[g]), float(h2u2[g])))
        disagree = bool(ok1[g] and ok2[g]
                        and abs(u1[g] - u2[g]) > SAME_MODE * max(1.0, abs(u1[g])))
        nonconc = bool(nc1[g] or (use2[g] and nc2[g]))
        cands = _distinct(cands)
        best = _pick(cands) if cands else None
        flat = best is None or best[2] > FLAT_H2
        trig = check_all or not state.seen[g] or disagree or nonconc or flat
        ev["starts_disagree"] += disagree
        ev["nonconcave_iterate"] += nonconc
        ev["start_not_converged"] += best is None
        if trig:
            ev["global_checks"] += 1
            cands, cap = _global_check(ctx, g, U[g], cands, state)
            if not cands:
                raise RuntimeError(f"no converged local maximum for gene {g}")
            best = _pick(cands)
            second = any(t is not best and t[1] >= best[1] - SECOND_MODE_GAP for t in cands)
            flat = best[2] > FLAT_H2
            ev["grid_cap"] += cap
            ev["second_mode"] += second
            ev["flat_mode"] += flat
            if cap or second or flat:
                out["choice"][g] = REF
            out["checked"][g] = True
        out["u_hat"][g], out["h_hat"][g], out["h2_hat"][g] = best
        out["sigma"][g] = _sigma(best[2])
        out["modes"][g] = cands
    return out


def _composite_gl(edges, n):
    nodes, wts = [], []
    for e0, e1, k in zip(edges[:-1], edges[1:], n):
        e = np.linspace(e0, e1, k + 1)
        mid, half = 0.5 * (e[:-1] + e[1:]), 0.5 * (e[1:] - e[:-1])
        nodes.append((mid[:, None] + half[:, None] * GL_X[None, :]).ravel())
        wts.append((half[:, None] * GL_W[None, :]).ravel())
    return np.concatenate(nodes), np.log(np.concatenate(wts))


def _reference(ctx, g, md, state):
    """Item 6: composite Gauss-Legendre, 20 nodes per panel, on [-W, W] with panel edges at +-W
    and at every local mode found; panels no wider than the narrowest mode's sigma; the panel
    count doubles until two successive values agree within 1e-10 (cap 8,192 panels)."""
    v = ctx.v
    cands = md["modes"][g]
    sigs = [s for s in (_sigma(t[2]) for t in cands) if np.isfinite(s)]
    sig_min = min(sigs) if sigs else np.sqrt(v)          # no defined sigma: prior SD
    s_hat = md["sigma"][g] if np.isfinite(md["sigma"][g]) else sig_min
    W = np.sqrt(2.0 * v * (-md["h_hat"][g] + REF_TAIL + max(0.0, np.log(np.sqrt(v) / s_hat))))
    inner = [t[0] for t in cands if -W < t[0] < W]
    edges = np.unique(np.concatenate([[-W, W], inner]))
    n = np.maximum(1, np.ceil(np.diff(edges) / sig_min)).astype(np.int64)
    prev = None
    while n.sum() <= REF_MAX_PANELS:
        nodes, logwt = _composite_gl(edges, n)
        hv = ctx.h_gene(g, nodes, 0)[0]
        ell = float(logsumexp(logwt + hv)) - 0.5 * np.log(2.0 * np.pi * v)
        state.events["max_reference_panels"] = max(state.events["max_reference_panels"],
                                                   int(n.sum()))
        if prev is not None and abs(ell - prev) <= REF_AGREE:
            return dict(ell=ell, nodes=nodes, logwt=logwt, h=hv, panels=int(n.sum()))
        prev = ell
        n = 2 * n
    raise RuntimeError(f"reference integrator: gene {g} did not agree within "
                       f"{REF_MAX_PANELS} panels")


PI_SKIP = 1e-40     # nodes with posterior weight below this add nothing to a double-precision sum


def _gene_node_terms(ctx, g, r, order, J1):
    """Fixed-node derivatives for gene g on the reference nodes: per-row posterior means of the
    row derivatives, E[d h/d v], and (order 2) the between-node covariance of d h_k / d(theta, v)
    on the gene's own columns. Same rule as GH20 with (u_hat, sigma_hat) held fixed (item 6)."""
    D, v, L = ctx.D, ctx.v, ctx.D.L
    idx = D.gene_idx(g)
    wg = D.w[idx]
    m = len(idx)
    A = r["logwt"] + r["h"]
    pi = np.exp(A - logsumexp(A))
    keep = pi >= PI_SKIP
    nodes, pi = r["nodes"][keep], pi[keep]
    qv = nodes ** 2 / (2.0 * v ** 2)
    t = dict(E1=np.zeros(m), E2=np.zeros(m), Eq_v=float(pi @ qv))
    if order >= 2:
        J1g, J2g = J1[idx], D.J2[idx]
        cols = np.union1d(np.union1d(J1g.indices, J2g.indices), [L.v]).astype(np.int64)
        pos_v = int(np.searchsorted(cols, L.v))
        J1d, J2d = J1g[:, cols].toarray(), J2g[:, cols].toarray()
        nc = len(cols)
        t.update(E11=np.zeros(m), E12=np.zeros(m), E22=np.zeros(m), cols=cols,
                 hvv=float(pi @ (-nodes ** 2 / v ** 3)) + 1.0 / (2.0 * v ** 2))
        S1, S2 = np.zeros(nc), np.zeros((nc, nc))
    for sl, uu, re, _ in ctx.gene_blocks(g, nodes, order):
        k = len(uu)
        p = pi[sl]
        g1, g2 = re.g1.reshape(k, m), re.g2.reshape(k, m)
        t["E1"] += p @ g1
        t["E2"] += p @ g2
        if order >= 2:
            t["E11"] += p @ re.g11.reshape(k, m)
            t["E12"] += p @ re.g12.reshape(k, m)
            t["E22"] += p @ re.g22.reshape(k, m)
            Dk = (g1 * wg) @ J1d + (g2 * wg) @ J2d
            Dk[:, pos_v] += qv[sl]
            S1 += p @ Dk
            S2 += (Dk * p[:, None]).T @ Dk
    if order >= 2:
        t["cov"] = S2 - np.outer(S1, S1)
    return t


def loglik4(theta, D, order=1, nodes=20, detail=False, state=None, final=False, force=None,
            all_reference=False):
    """Section-4 log-likelihood sum_g log int N(u; 0, v) prod_i L_ig(u)^{w_i} du.

    GH `nodes` (20) centred at the v1.8 mode with sigma_hat = (-h''(u_hat))^(-1/2), or the
    reference integrator where v1.8 section 1 items 4, 5 and 7 require it. `state` carries the
    previous modes within one fit (fresh state: every gene is at its first evaluation, so the
    global check runs for all). `final` marks a fit's final converged point (item 4 (v)).
    `force` (array of GH/REF per gene) holds the integrators fixed (item 8 re-polish);
    `all_reference` uses the reference integrator for every gene (v1.8 section 2 item 3).
    Gradients and Hessian hold the nodes fixed. v = 0 is the separate branch (section 3, with
    d/dv the boundary score)."""
    L = D.L
    v = float(theta[L.v])
    if v == 0.0:
        if state is not None:
            state.last_choice = None
        res = loglik3(theta, D, order)
        if order == 0:
            return (res, dict(ell_g=gene_loglik3(theta, D), integrator=None)) if detail else res
        res = list(res)
        res[1] = res[1].copy()
        res[1][L.v] = boundary_score(theta, D)
        if order == 2:
            res[2] = res[2].copy()
            res[2][L.v, :] = np.nan
            res[2][:, L.v] = np.nan
        return tuple(res)
    if state is None:
        state = Sec4State(L.G)
    ev_count = state.events
    ev_count["evaluations"] += 1
    x1b, x2, al, q = predictors(theta, D)
    ctx = _Ctx(D, x1b, x2, v)
    G, gene, St, w = L.G, D.gene, D.St, D.w
    has = D.gene_has_rows
    md = _modes(ctx, state, final or all_reference)
    choice = md["choice"].copy() if force is None else np.asarray(force).copy()
    if all_reference:
        choice[:] = REF
    choice[~has] = GH
    u_hat = md["u_hat"]
    sig_gh = np.where(np.isfinite(md["sigma"]), md["sigma"], np.sqrt(v))
    run_gh = not all_reference
    ell_gh = np.full(G, np.nan)
    if run_gh:
        xk, wk = hermgauss(nodes)
        logW = np.log(wk) + xk ** 2
        U = u_hat[None, :] + SQRT2 * sig_gh[None, :] * xk[:, None]        # nodes x G
        Hk = np.empty((nodes, G))
        evs = []
        for k in range(nodes):
            re = row_eval(x1b + U[k][gene] * St, x2, D.rows, order)
            Hk[k] = np.bincount(gene, w * re.logL, G) - U[k] ** 2 / (2.0 * v)
            evs.append(re)
        A = logW[:, None] + Hk
        lse = logsumexp(A, axis=0)
        ell_gh = lse + np.log(SQRT2 * sig_gh) - 0.5 * np.log(2.0 * np.pi * v)
    ell_g = ell_gh.copy()
    refs = {}
    if force is None and not all_reference:
        for g in np.flatnonzero(md["checked"] & (choice == GH) & has):
            r = _reference(ctx, g, md, state)
            ev_count["reference_checks"] += 1
            if abs(r["ell"] - ell_gh[g]) > GH_REF_TOL:
                choice[g] = REF
                ev_count["gh_ref_mismatch"] += 1
                refs[g] = r
    for g in np.flatnonzero((choice == REF) & has):
        if g not in refs:
            refs[g] = _reference(ctx, g, md, state)
        ell_g[g] = refs[g]["ell"]
    ev_count["reference_used"] += len(refs)
    state.u_prev = u_hat.copy()
    state.seen |= has
    state.last_choice = choice.copy()
    ll = float(np.sum(ell_g[has]))
    info = dict(ell_g=ell_g, ell_gh=ell_gh, u_hat=u_hat, sigma=md["sigma"], integrator=choice,
                checked=md["checked"], n_modes=np.array([len(c) for c in md["modes"]]),
                reference_panels={int(g): r["panels"] for g, r in refs.items()})
    if order == 0:
        return (ll, info) if detail else ll
    nr = D.nr
    E1, E2 = np.zeros(nr), np.zeros(nr)
    if order >= 2:
        E11, E12, E22 = np.zeros(nr), np.zeros(nr), np.zeros(nr)
    Eq_v, hvv = np.zeros(G), np.zeros(G)
    if run_gh:
        pi = np.exp(A - lse[None, :])                              # nodes x G
        P = pi[:, gene]                                            # nodes x rows
        E1 = sum(P[k] * evs[k].g1 for k in range(nodes))
        E2 = sum(P[k] * evs[k].g2 for k in range(nodes))
        q_v = U ** 2 / (2.0 * v ** 2)                              # d h_k / d v
        Eq_v = np.sum(pi * q_v, axis=0)
        if order >= 2:
            E11 = sum(P[k] * evs[k].g11 for k in range(nodes))
            E12 = sum(P[k] * evs[k].g12 for k in range(nodes))
            E22 = sum(P[k] * evs[k].g22 for k in range(nodes))
            hvv = np.sum(pi * (-U ** 2 / v ** 3), axis=0) + 1.0 / (2.0 * v ** 2)
    J1 = jacobian1(D, al, q)
    terms = {}
    for g, r in refs.items():
        t = _gene_node_terms(ctx, g, r, order, J1)
        idx = D.gene_idx(g)
        E1[idx], E2[idx] = t["E1"], t["E2"]
        Eq_v[g] = t["Eq_v"]
        if order >= 2:
            E11[idx], E12[idx], E22[idx] = t["E11"], t["E12"], t["E22"]
            hvv[g] = t["hvv"]
        terms[g] = t
    dv = float(np.sum((Eq_v - 1.0 / (2.0 * v))[has]))
    if order == 1:
        grad, _ = _assemble(D, J1, w * E1, w * E2)
        grad[L.v] = dv
        return (ll, grad, info) if detail else (ll, grad)
    grad, H = _assemble(D, J1, w * E1, w * E2, w * E11, w * E12, w * E22)
    grad[L.v] = dv
    if run_gh:
        # between-node covariance of d h_k / d(theta, v) within each GH gene
        gh_gene = ((choice == GH) & has).astype(float)
        gh_row = gh_gene[gene]
        gidx = np.arange(G)
        blocks = []
        for k in range(nodes):
            sq = np.sqrt(P[k]) * gh_row
            Mk = D.Agg @ (sp.diags(w * sq * (evs[k].g1 - E1)) @ J1
                          + sp.diags(w * sq * (evs[k].g2 - E2)) @ D.J2)
            vcol = np.sqrt(pi[k]) * (q_v[k] - Eq_v) * gh_gene
            Mk = Mk + sp.csr_matrix((vcol, (gidx, np.full(G, L.v))), shape=(G, L.size))
            blocks.append(Mk)
        M = sp.vstack(blocks).tocsr()
        H += (M.T @ M).toarray()
    for g, t in terms.items():
        H[np.ix_(t["cols"], t["cols"])] += t["cov"]
    H[L.v, L.v] += float(np.sum(hvv[has]))
    return (ll, grad, H, info) if detail else (ll, grad, H)


def gh_gene_integral_quad(theta, D, g, v, width=40.0):
    """Reference ell_g by scipy.integrate.quad (epsabs 1e-12) for the section 4 accuracy test,
    centred at the v1.8 mode (global check run for every gene)."""
    from scipy.integrate import quad
    th = np.array(theta, float)
    th[D.L.v] = v
    _, det = loglik4(th, D, 0, detail=True, final=True)
    x1b, x2, _, _ = predictors(th, D)
    sel = D.gene == g
    R = D.rows.take(np.flatnonzero(sel))
    x1g, x2g, Sg, wg = x1b[sel], x2[sel], D.St[sel], D.w[sel]

    def h(u):
        return float(np.sum(wg * row_eval(x1g + u * Sg, x2g, R, 0).logL)) - u ** 2 / (2.0 * v)

    c = det["u_hat"][g]
    s = det["sigma"][g] if np.isfinite(det["sigma"][g]) else np.sqrt(v)
    h0 = h(c)
    val, err = quad(lambda u: np.exp(h(u) - h0), c - width * s, c + width * s,
                    epsabs=1e-12, epsrel=1e-13, limit=500, points=[c])
    return h0 + np.log(val) - 0.5 * np.log(2.0 * np.pi * v), err
