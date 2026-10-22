"""v1.8.3-gpu: torch float64 evaluators for model_a.loglik3, model_a.loglik4 and model_a.boundary_score.

Selected by model_a's backend switch (MODEL_A_BACKEND=torch). Built from the throwaway port in
A6_gateB/gpu-spike-s4-20261016T121756Z/scripts/torch_sec4.py, with these changes:
  - deterministic mode: torch.use_deterministic_algorithms(True), CUBLAS_WORKSPACE_CONFIG=:4096:8 (set before
    CUDA starts), TF32 off;
  - the fixed-node derivative terms of genes on the reference integrator (model_a._gene_node_terms) run on the
    device too (the spike fell back to the CPU code);
  - the item 3 / item 4 local searches call model_a._local_search itself (subset=True), not a copy;
  - one device model per Data object, cached on the object, no reference cycle;
  - kappa_A = True is refused (no Gate B or B-FIX problem uses it).

What runs where. On the device: every row likelihood term and its x-derivatives (beta-binomial terms, the
truncation tails and direct sums of log C_H, the mixture over genotype classes), all per-gene sums h_g(u), h', h'',
the GH nodes, the gradient and Hessian assembly, the reference-integrator panels and the fixed-node terms. On the
host (numpy, the CPU code's own control flow): the item 3 bracketing/Newton/golden decisions per search, the
item 4 grid peak picking, _distinct/_pick, the item 5/7 integrator decisions, the item 6 panel doubling decision,
the log-sum-exp over GH nodes and panels (scipy, as in the CPU code), and the Sec4State bookkeeping.

Arithmetic order differs from the CPU code in the reductions (index_add_ / matmul / sum versus np.bincount /
scipy.sparse / np.sum), so results agree to rounding, not bitwise. Within one device type, deterministic mode gives
the same bits run to run.

Batching across datasets is not implemented: a ragged batch over several datasets would change the grouping of
the segment sums of each dataset, so per-dataset arithmetic would depend on the batch.
"""
import os

os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

import time  # noqa: E402
import weakref  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402
from numpy.polynomial.hermite import hermgauss  # noqa: E402
from scipy.special import logsumexp  # noqa: E402

import model_a as ma  # noqa: E402

torch.use_deterministic_algorithms(True)
torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False

DEVICE = os.environ.get("MODEL_A_TORCH_DEVICE", "cuda")
SERIES_X = 30.0
EPS_RHO = ma.EPS_RHO
LOG_HALF = float(np.log(0.5))
SQRT2 = float(np.sqrt(2.0))
DT = torch.float64

# Seconds and calls per evaluator (host wall time of the call, device synchronized at the end because every call
# returns host arrays). Read by drivers that report per-unit device time.
STATS = dict(calls={}, seconds={}, dev_models_built=0, dev_model_build_s=0.0)


_DEV_ERRORS = tuple(e for e in (getattr(torch, "OutOfMemoryError", None), getattr(torch, "AcceleratorError", None))
                    if e is not None)


class DevFatal(BaseException):
    """Device failure (out of memory, CUDA error). A BaseException so that the fit code's except clauses do not
    count it as a failed fit; the run stops instead."""


def call(name, D, theta, *args, **kw):
    key = name if name == "boundary_score" else f"{name}.o{args[0] if args else kw.get('order', 1)}"
    t0 = time.perf_counter()
    try:
        dev = device_model(D)
        out = getattr(dev, name)(np.asarray(theta, float), *args, **kw)
    except _DEV_ERRORS as err:
        raise DevFatal(repr(err)) from err
    STATS["calls"][key] = STATS["calls"].get(key, 0) + 1
    STATS["seconds"][key] = STATS["seconds"].get(key, 0.0) + time.perf_counter() - t0
    return out


def stats_snapshot():
    return dict(calls=dict(STATS["calls"]), seconds=dict(STATS["seconds"]),
                dev_models_built=STATS["dev_models_built"], dev_model_build_s=STATS["dev_model_build_s"],
                total_eval_s=float(sum(STATS["seconds"].values())))


def device_model(D):
    dev = getattr(D, "_torch_dev", None)
    if dev is None or dev.data_ref() is not D:
        t0 = time.perf_counter()
        dev = Dev(D, DEVICE)
        D._torch_dev = dev
        STATS["dev_models_built"] += 1
        STATS["dev_model_build_s"] += time.perf_counter() - t0
    return dev


# ------------------------------------------------------------------------------- special functions (betabinom.py)
def _lser(y):
    y2 = y * y
    return (1 / 12 - (1 / 360 - (1 / 1260 - (1 / 1680 - (1 / 1188 - (691 / 360360
            - 1 / (156 * y2)) / y2) / y2) / y2) / y2) / y2) / y


def _dser(y):
    y2 = y * y
    return -(1 / 12 - (1 / 120 - (1 / 252 - (1 / 240 - (1 / 132 - (691 / 32760
             - 1 / (12 * y2)) / y2) / y2) / y2) / y2) / y2) / y2


def _tser(y):
    y2 = y * y
    return -(1 / 30 - (1 / 42 - (1 / 30 - (5 / 66 - (691 / 2730 - 7 / (6 * y2)) / y2) / y2)
             / y2) / y2) / (y2 * y2 * y)


def lgamma_diff(x, m):
    """betabinom.lgamma_diff: exact difference below x = 30, Stirling series from 30 up."""
    big = x >= SERIES_X
    small = torch.lgamma(x + m) - torch.lgamma(x)
    xb = torch.where(big, x, torch.full_like(x, SERIES_X))
    z = xb + m
    ser = (xb - 0.5) * torch.log1p(m / xb) + m * torch.log(z) - m + (_lser(z) - _lser(xb))
    return torch.where(big, ser, small)


def digamma_diff(x, m):
    big = x >= SERIES_X
    small = torch.digamma(x + m) - torch.digamma(x)
    xb = torch.where(big, x, torch.full_like(x, SERIES_X))
    z = xb + m
    ser = torch.log1p(m / xb) + m / (2 * xb * z) + (_dser(z) - _dser(xb))
    return torch.where(big, ser, small)


def trigamma_diff(x, m):
    big = x >= SERIES_X
    two = torch.full_like(x, 2.0)
    small = torch.special.zeta(two, x + m) - torch.special.zeta(two, x)
    xb = torch.where(big, x, torch.full_like(x, SERIES_X))
    z = xb + m
    ser = (-m / (xb * z) - m * (xb + z) / (2 * xb ** 2 * z ** 2)
           - m * (xb ** 2 + xb * z + z ** 2) / (6 * xb ** 3 * z ** 3)
           + (_tser(z) - _tser(xb)))
    return torch.where(big, ser, small)


def lchoose(n, a):
    return torch.lgamma(n + 1.0) - torch.lgamma(a + 1.0) - torch.lgamma(n - a + 1.0)


def segsum(x, seg, nseg):
    return torch.zeros(nseg, dtype=x.dtype, device=x.device).index_add_(0, seg, x)


def seg_lse(x, seg, nseg):
    m = torch.full((nseg,), -float("inf"), dtype=x.dtype, device=x.device)
    m = m.scatter_reduce(0, seg, x, reduce="amax", include_self=True)
    mf = torch.where(torch.isfinite(m), m, torch.zeros_like(m))
    s = torch.zeros(nseg, dtype=x.dtype, device=x.device).index_add_(0, seg, torch.exp(x - mf[seg]))
    return torch.log(s) + mf


def ragged(starts, counts, dev):
    """Flattened ranges: segment id, absolute index and offset for starts[j] .. starts[j] + counts[j] - 1."""
    total = int(counts.sum())
    seg = torch.repeat_interleave(torch.arange(len(counts), device=dev), counts, output_size=total)
    first = torch.cumsum(counts, 0) - counts
    offs = torch.arange(total, device=dev) - first[seg]
    return seg, starts[seg] + offs, offs


def chunks_by_cost(cost, budget):
    """Contiguous [s, e) pieces of about `budget` cumulative cost (a piece may exceed it by one item)."""
    cost = np.asarray(cost, np.int64)
    if len(cost) == 0:
        return []
    cs = np.cumsum(cost)
    pid = (cs - cost) // max(int(budget), 1)
    cut = np.flatnonzero(np.diff(pid)) + 1
    b = np.concatenate([[0], cut, [len(cost)]])
    return list(zip(b[:-1].tolist(), b[1:].tolist()))


def logaddexp(a, b):
    m = torch.maximum(a, b)
    out = m + torch.log1p(torch.exp(-torch.abs(a - b)))
    return torch.where(torch.isinf(m) & (m < 0), m, out)


class _RS:
    """model_a._RowState on the device (the shape halves are recomputed per entry, not cached per row)."""

    def __init__(self, x1, x2, n, order):
        self.n = n
        self.mu, self.omu = torch.sigmoid(x1), torch.sigmoid(-x1)
        sig = torch.sigmoid(x2)
        rho = EPS_RHO + (1.0 - 2.0 * EPS_RHO) * sig
        t = (1.0 - rho) / rho
        self.t = t
        self.al, self.be = self.mu * t, self.omu * t
        self.lgc = -lgamma_diff(t, n)
        if order >= 1:
            self.Dt = digamma_diff(t, n)
            self.m1 = self.mu * self.omu
            rp = (1.0 - 2.0 * EPS_RHO) * sig * (1.0 - sig)
            self.tp = -rp / rho ** 2
        if order >= 2:
            self.Tt = trigamma_diff(t, n)
            rpp = rp * (1.0 - 2.0 * sig)
            self.tpp = -rpp / rho ** 2 + 2.0 * rp ** 2 / rho ** 3

    def terms(self, j, a, lch, order):
        g = (lambda v: v) if j is None else (lambda v: v[j])  # noqa: E731
        al, be, n = g(self.al), g(self.be), g(self.n)
        lf = lch + lgamma_diff(al, a) + lgamma_diff(be, n - a) + g(self.lgc)
        if order == 0:
            return (lf,)
        mu, omu, t, m1, tp = g(self.mu), g(self.omu), g(self.t), g(self.m1), g(self.tp)
        Da = digamma_diff(al, a)
        Db = digamma_diff(be, n - a)
        Dab = Da - Db
        Dmix = mu * Da + omu * Db - g(self.Dt)
        d1 = m1 * t * Dab
        d2 = tp * Dmix
        if order == 1:
            return lf, d1, d2
        Ta = trigamma_diff(al, a)
        Tb = trigamma_diff(be, n - a)
        d11 = m1 * (omu - mu) * t * Dab + (m1 * t) ** 2 * (Ta + Tb)
        d12 = m1 * tp * (Dab + t * (mu * Ta - omu * Tb))
        d22 = g(self.tpp) * Dmix + tp ** 2 * (mu ** 2 * Ta + omu ** 2 * Tb - g(self.Tt))
        return lf, d1, d2, d11, d12, d22


# ------------------------------------------------------------------------------- device model of one Data object
class Dev:
    def __init__(self, D, device, budget0=12_000_000, budget12=4_000_000, direct_budget=6_000_000):
        L = D.L
        if L.kappa_A:
            raise NotImplementedError("backend_torch covers kappa_A = False only")
        self.data_ref = weakref.ref(D)
        self.L, self.dev, self.dt = L, torch.device(device), DT
        self.budget = {0: budget0, 1: budget12, 2: budget12}
        self.direct_budget = direct_budget
        R = D.rows
        dev, dt = self.dev, self.dt
        f = lambda v: torch.as_tensor(np.asarray(v, float), dtype=dt, device=dev)   # noqa: E731
        i = lambda v: torch.as_tensor(np.asarray(v, np.int64), device=dev)          # noqa: E731
        self.f, self.i = f, i
        nr, G = D.nr, L.G
        self.nr, self.G, self.P = nr, G, L.size
        self.n, self.lo, self.hi = f(R.n), i(R.lo), i(R.hi)
        self.klen = self.hi - self.lo + 1
        self.a, self.a_lch = f(R.a), f(R.a_lchoose)
        self.logPH, self.logcN, self.logcD = f(R.logPH), f(R.logcN), f(R.logcD)
        order = np.argsort(R.tail_row, kind="stable")
        self.t_a, self.t_lch = f(R.tail_a[order]), f(R.tail_lchoose[order])
        tcount = np.bincount(R.tail_row, minlength=nr)
        self.tcount = i(tcount)
        self.tstart = i(np.concatenate([[0], np.cumsum(tcount)[:-1]]))
        gorder = np.argsort(D.gene, kind="stable")             # the row order of Data.gene_idx
        gcount = np.bincount(D.gene, minlength=G)
        self.gstart_np = np.concatenate([[0], np.cumsum(gcount)[:-1]]).astype(np.int64)
        self.gcount_np = gcount.astype(np.int64)
        self.g_rows = i(gorder)
        self.gcount = i(gcount)
        self.gstart = i(self.gstart_np)
        self.gcost = gcount + np.bincount(D.gene, weights=tcount, minlength=G).astype(np.int64)
        self.rowcost = 1 + tcount
        self.gene = i(D.gene)
        self.St, self.w, self.Z, self.b, self.d = f(D.St), f(D.w), f(D.Z), f(D.b), f(D.d)
        self.s_t, self.omega = f(D.s_t), f(D.omega)
        self.ix_alpha, self.ix_r, self.ix_delta = i(L.alpha), i(L.r), i(L.delta)
        self.arange_nr = torch.arange(nr, device=dev)
        J2 = torch.zeros((nr, self.P), dtype=dt, device=dev)
        J2[self.arange_nr, self.ix_r[self.gene]] = 1.0
        J2[:, L.rho_S], J2[:, L.rho_b], J2[:, L.rho_d] = self.St, self.b, self.d
        self.J2d = J2
        self.has = D.gene_has_rows.copy()
        self.theta = None

    # ---------------------------------------------------------------- predictors (model_a.predictors)
    def set_theta(self, theta):
        L = self.L
        self.theta = theta
        th = self.f(theta)
        self.v = float(theta[L.v])
        al = th[self.ix_alpha][self.gene]
        q = 1.0 + th[L.kappa] * self.St + self.Z @ th[self.ix_delta]
        self.al, self.q = al, q
        self.x1b = al * q + self.s_t * (self.omega + th[L.omega_S] * self.St)
        self.x2 = th[self.ix_r][self.gene] + th[L.rho_S] * self.St + th[L.rho_b] * self.b + th[L.rho_d] * self.d

    # ---------------------------------------------------------------- row likelihood (model_a.row_eval, log_C_H)
    def row_eval(self, ridx, x1, order):
        N = int(ridx.numel())
        n = self.n[ridx]
        st = _RS(x1, self.x2[ridx], n, order)
        num = st.terms(None, self.a[ridx], self.a_lch[ridx], order)
        seg, tix, _ = ragged(self.tstart[ridx], self.tcount[ridx], self.dev)
        tl = st.terms(seg, self.t_a[tix], self.t_lch[tix], order)
        logT = seg_lse(tl[0], seg, N)
        direct = logT >= LOG_HALF
        C = -torch.expm1(logT)
        logC = torch.log1p(-torch.exp(logT))
        ch = [logC]
        if order >= 1:
            fT = torch.exp(tl[0])
            Cc = torch.where(direct, torch.ones_like(C), C)
            c1 = -segsum(fT * tl[1], seg, N) / Cc
            c2 = -segsum(fT * tl[2], seg, N) / Cc
            ch += [c1, c2]
            if order >= 2:
                ch.append(-segsum(fT * (tl[3] + tl[1] ** 2), seg, N) / Cc - c1 ** 2)
                ch.append(-segsum(fT * (tl[4] + tl[1] * tl[2]), seg, N) / Cc - c1 * c2)
                ch.append(-segsum(fT * (tl[5] + tl[2] ** 2), seg, N) / Cc - c2 ** 2)
            del fT
        del tl, seg, tix
        rd = torch.nonzero(direct).reshape(-1)
        if rd.numel():
            kl_all = self.klen[ridx][rd]
            lo_all = self.lo[ridx][rd]
            for s, e in chunks_by_cost(kl_all.cpu().numpy(), self.direct_budget):
                rds, kl, lo = rd[s:e], kl_all[s:e], lo_all[s:e]
                seg2, _, offs = ragged(torch.zeros_like(kl), kl, self.dev)
                j = rds[seg2]
                a = (lo[seg2] + offs).to(self.dt)
                kt = st.terms(j, a, lchoose(n[j], a), order)
                m = e - s
                lC = seg_lse(kt[0], seg2, m)
                ch[0][rds] = lC
                if order >= 1:
                    pi = torch.exp(kt[0] - lC[seg2])
                    e1 = segsum(pi * kt[1], seg2, m)
                    e2 = segsum(pi * kt[2], seg2, m)
                    ch[1][rds], ch[2][rds] = e1, e2
                    if order >= 2:
                        ch[3][rds] = segsum(pi * (kt[3] + kt[1] ** 2), seg2, m) - e1 ** 2
                        ch[4][rds] = segsum(pi * (kt[4] + kt[1] * kt[2]), seg2, m) - e1 * e2
                        ch[5][rds] = segsum(pi * (kt[5] + kt[2] ** 2), seg2, m) - e2 ** 2
        logPH, logcN, logcD = self.logPH[ridx], self.logcN[ridx], self.logcD[ridx]
        logN = logaddexp(logPH + num[0], logcN)
        logD = logaddexp(logPH + ch[0], logcD)
        out = dict(logL=logN - logD, logC_H=ch[0])
        if order >= 1:
            wN, vN = torch.exp(logPH + num[0] - logN), torch.exp(logcN - logN)
            wD, vD = torch.exp(logPH + ch[0] - logD), torch.exp(logcD - logD)
            out["g1"] = wN * num[1] - wD * ch[1]
            out["g2"] = wN * num[2] - wD * ch[2]
            if order >= 2:
                out["g11"] = wN * num[3] + wN * vN * num[1] ** 2 - wD * ch[3] - wD * vD * ch[1] ** 2
                out["g12"] = (wN * num[4] + wN * vN * num[1] * num[2]
                              - wD * ch[4] - wD * vD * ch[1] * ch[2])
                out["g22"] = wN * num[5] + wN * vN * num[2] ** 2 - wD * ch[5] - wD * vD * ch[2] ** 2
        return out

    def row_eval_all(self, x1, order):
        """All rows once, in pieces by tail budget (memory only); outputs concatenated in row order."""
        pieces = chunks_by_cost(self.rowcost, self.budget[order])
        outs = [self.row_eval(self.arange_nr[s:e], x1[s:e], order) for s, e in pieces]
        if len(outs) == 1:
            return outs[0]
        return {k: torch.cat([o[k] for o in outs]) for k in outs[0]}

    # ---------------------------------------------------------------- h_g(u) for many (gene, u)
    def h_multi(self, gk, uk, order):
        """h_g(u) = sum_i w_i log L_ig(u) - u^2/(2v) for queries (gk[k], uk[k]); order 2 adds h', h''."""
        gk = np.asarray(gk, np.int64)
        uk = np.asarray(uk, float)
        K = len(gk)
        h = np.empty(K)
        h1 = np.empty(K) if order >= 2 else None
        h2 = np.empty(K) if order >= 2 else None
        v = self.v
        for s, e in chunks_by_cost(self.gcost[gk], self.budget[order]):
            g_t, u_t = self.i(gk[s:e]), self.f(uk[s:e])
            seg, pos, _ = ragged(self.gstart[g_t], self.gcount[g_t], self.dev)
            ridx = self.g_rows[pos]
            Sg = self.St[ridx]
            re = self.row_eval(ridx, self.x1b[ridx] + u_t[seg] * Sg, order)
            wr = self.w[ridx]
            m = e - s
            h[s:e] = (segsum(wr * re["logL"], seg, m) - u_t ** 2 / (2.0 * v)).cpu().numpy()
            if order >= 2:
                h1[s:e] = (segsum(wr * Sg * re["g1"], seg, m) - u_t / v).cpu().numpy()
                h2[s:e] = (segsum(wr * Sg ** 2 * re["g11"], seg, m) - 1.0 / v).cpu().numpy()
        return h, h1, h2

    def _f_sub(self, genes):
        """f(x, act) for model_a._local_search(subset=True): one search per entry of `genes`."""
        genes = np.asarray(genes, np.int64)

        def f(x, act):
            hx, h1x, h2x = (np.full(len(x), np.nan) for _ in range(3))
            idx = np.flatnonzero(act)
            if len(idx):
                hx[idx], h1x[idx], h2x[idx] = self.h_multi(genes[idx], x[idx], 2)
            return hx, h1x, h2x
        return f

    # ---------------------------------------------------------------- assembly (model_a._assemble)
    def assemble1(self, G1, G2):
        L, G = self.L, self.G
        g = torch.zeros(self.P, dtype=self.dt, device=self.dev)
        g[self.ix_alpha] = segsum(G1 * self.q, self.gene, G)
        g[L.kappa] = (G1 * self.al * self.St).sum()
        g[self.ix_delta] = (G1 * self.al) @ self.Z
        g[L.omega_S] = (G1 * self.s_t * self.St).sum()
        g[self.ix_r] = segsum(G2, self.gene, G)
        g[L.rho_S] = (G2 * self.St).sum()
        g[L.rho_b] = (G2 * self.b).sum()
        g[L.rho_d] = (G2 * self.d).sum()
        return g

    def J1rows(self, ridx=None):
        """Dense d x1 / d theta (model_a.jacobian1) for rows ridx (all rows if None)."""
        L = self.L
        if ridx is None:
            ridx = self.arange_nr
        m = int(ridx.numel())
        ar = torch.arange(m, device=self.dev)
        al, St = self.al[ridx], self.St[ridx]
        J1 = torch.zeros((m, self.P), dtype=self.dt, device=self.dev)
        J1[ar, self.ix_alpha[self.gene[ridx]]] = self.q[ridx]
        J1[:, L.kappa] = al * St
        J1[:, self.ix_delta] = al[:, None] * self.Z[ridx]
        J1[:, L.omega_S] = self.s_t[ridx] * St
        return J1

    def assemble2(self, G1, G2, G11, G12, G22):
        L, G = self.L, self.G
        grad = self.assemble1(G1, G2)
        J1, J2 = self.J1rows(), self.J2d
        H = J1.T @ (G11[:, None] * J1)
        X = J1.T @ (G12[:, None] * J2)
        H += X + X.T
        del X
        H += J2.T @ (G22[:, None] * J2)
        del J1
        shared = [(L.kappa, self.St)] + [(int(L.delta[k]), self.Z[:, k]) for k in range(L.p)]
        for col, x in shared:
            s = segsum(G1 * x, self.gene, G)
            H[self.ix_alpha, col] += s
            H[col, self.ix_alpha] += s
        return grad, H

    # ---------------------------------------------------------------- section 3
    def loglik3(self, theta, order=1):
        self.set_theta(theta)
        re = self.row_eval_all(self.x1b, order)
        w = self.w
        ll = float((w * re["logL"]).sum())
        if order == 0:
            return ll
        if order == 1:
            return ll, self.assemble1(w * re["g1"], w * re["g2"]).cpu().numpy()
        grad, H = self.assemble2(w * re["g1"], w * re["g2"], w * re["g11"], w * re["g12"], w * re["g22"])
        return ll, grad.cpu().numpy(), H.cpu().numpy()

    def _boundary_from(self, re):
        G = self.G
        s1 = segsum(self.w * re["g1"] * self.St, self.gene, G)
        s2 = segsum(self.w * re["g11"] * self.St ** 2, self.gene, G)
        return 0.5 * float((s1 ** 2 + s2).sum())

    def boundary_score(self, theta):
        self.set_theta(theta)
        return self._boundary_from(self.row_eval_all(self.x1b, 2))

    # ---------------------------------------------------------------- item 4 global check, batched over genes
    def global_check_batch(self, genes, U, cands_map, state):
        specs = []
        for g in genes:
            cands = cands_map[g]
            best = ma._pick(cands) if cands else None
            sig = ma._sigma(best[2]) if best is not None else np.nan
            cap = not np.isfinite(sig)
            N = ma.GRID_CAP if cap else max(ma.GRID_MIN, int(np.ceil(4.0 * U[g] / sig)) + 1)
            if N > ma.GRID_CAP:
                N, cap = ma.GRID_CAP, True
            state.events["max_grid_points"] = max(state.events["max_grid_points"], N)
            specs.append((g, np.linspace(-U[g], U[g], N), cap))
        if not specs:
            return {}
        qg = np.concatenate([np.full(len(gr), g, np.int64) for g, gr, _ in specs])
        qu = np.concatenate([gr for _, gr, _ in specs])
        hall = self.h_multi(qg, qu, 0)[0]
        off = 0
        sg, sa, sc, sb, sha, shb = [], [], [], [], [], []
        for g, grid, cap in specs:
            hg = hall[off:off + len(grid)]
            off += len(grid)
            k = np.flatnonzero((hg[1:-1] > hg[:-2]) & (hg[1:-1] > hg[2:])) + 1
            for kk in k:
                sg.append(g)
                sa.append(grid[kk - 1])
                sc.append(grid[kk])
                sb.append(grid[kk + 1])
                sha.append(hg[kk - 1])
                shb.append(hg[kk + 1])
        out = {g: (list(cands_map[g]), cap) for g, _, cap in specs}
        if sg:
            sg = np.asarray(sg, np.int64)
            hc, h1c, h2c = self.h_multi(sg, np.asarray(sc), 2)
            uu, hh, hh2, ok, _ = ma._local_search(self._f_sub(sg), sa, sc, sb, sha, hc, shb, h1c, h2c,
                                                  subset=True)
            for x, y, z, o, g in zip(uu, hh, hh2, ok, sg):
                if o:
                    out[int(g)][0].append((float(x), float(y), float(z)))
        return {g: (ma._distinct(c), cap) for g, (c, cap) in out.items()}

    # ---------------------------------------------------------------- model_a._modes
    def modes(self, state, check_all):
        v, G = self.v, self.G
        has = self.has
        allg = np.arange(G)
        zero = np.zeros(G)
        h0, h10, h20 = self.h_multi(allg, zero, 2)
        U = np.sqrt(2.0 * v * np.maximum(-h0, 0.0))
        hUp, hUm = self.h_multi(allg, U, 0)[0], self.h_multi(allg, -U, 0)[0]
        f = self._f_sub(allg)
        u1, hu1, h2u1, ok1, nc1 = ma._local_search(f, -U, zero, U, hUm, h0, hUp, h10, h20, active=has,
                                                   subset=True)
        prev = state.u_prev
        use2 = has & np.isfinite(prev) & (prev > -U) & (prev < U)
        ok2 = np.zeros(G, bool)
        nc2 = np.zeros(G, bool)
        if use2.any():
            p0 = np.where(use2, prev, 0.0)
            hp, h1p, h2p = self.h_multi(allg, p0, 2)
            use2 &= hp >= np.maximum(hUm, hUp)
            if use2.any():
                u2, hu2, h2u2, ok2, nc2 = ma._local_search(f, -U, p0, U, hUm, hp, hUp, h1p, h2p, active=use2,
                                                           subset=True)
                ok2 &= use2
        ev = state.events
        out = dict(U=U, u_hat=np.zeros(G), h_hat=np.zeros(G), h2_hat=np.full(G, -1.0 / v),
                   sigma=np.full(G, np.sqrt(v)), choice=np.full(G, ma.GH), checked=np.zeros(G, bool),
                   modes=[[] for _ in range(G)])
        trig_genes, cmap = [], {}
        for g in np.flatnonzero(has):
            cands = [(float(u1[g]), float(hu1[g]), float(h2u1[g]))] if ok1[g] else []
            if ok2[g]:
                cands.append((float(u2[g]), float(hu2[g]), float(h2u2[g])))
            disagree = bool(ok1[g] and ok2[g]
                            and abs(u1[g] - u2[g]) > ma.SAME_MODE * max(1.0, abs(u1[g])))
            nonconc = bool(nc1[g] or (use2[g] and nc2[g]))
            cands = ma._distinct(cands)
            best = ma._pick(cands) if cands else None
            flat = best is None or best[2] > ma.FLAT_H2
            trig = check_all or not state.seen[g] or disagree or nonconc or flat
            ev["starts_disagree"] += disagree
            ev["nonconcave_iterate"] += nonconc
            ev["start_not_converged"] += best is None
            if trig:
                ev["global_checks"] += 1
                trig_genes.append(int(g))
                cmap[int(g)] = cands
            else:
                out["u_hat"][g], out["h_hat"][g], out["h2_hat"][g] = best
                out["sigma"][g] = ma._sigma(best[2])
                out["modes"][g] = cands
        res = self.global_check_batch(trig_genes, U, cmap, state)
        for g in trig_genes:
            cands, cap = res[g]
            if not cands:
                raise RuntimeError(f"no converged local maximum for gene {g}")
            best = ma._pick(cands)
            second = any(t is not best and t[1] >= best[1] - ma.SECOND_MODE_GAP for t in cands)
            flat = best[2] > ma.FLAT_H2
            ev["grid_cap"] += cap
            ev["second_mode"] += second
            ev["flat_mode"] += flat
            if cap or second or flat:
                out["choice"][g] = ma.REF
            out["checked"][g] = True
            out["u_hat"][g], out["h_hat"][g], out["h2_hat"][g] = best
            out["sigma"][g] = ma._sigma(best[2])
            out["modes"][g] = cands
        return out

    # ---------------------------------------------------------------- item 6 reference integrator, batched
    def reference_batch(self, genes, md, state):
        v = self.v
        work = {}
        for g in genes:
            cands = md["modes"][g]
            sigs = [s for s in (ma._sigma(t[2]) for t in cands) if np.isfinite(s)]
            sig_min = min(sigs) if sigs else np.sqrt(v)
            s_hat = md["sigma"][g] if np.isfinite(md["sigma"][g]) else sig_min
            W = np.sqrt(2.0 * v * (-md["h_hat"][g] + ma.REF_TAIL + max(0.0, np.log(np.sqrt(v) / s_hat))))
            inner = [t[0] for t in cands if -W < t[0] < W]
            edges = np.unique(np.concatenate([[-W, W], inner]))
            n = np.maximum(1, np.ceil(np.diff(edges) / sig_min)).astype(np.int64)
            work[g] = dict(edges=edges, n=n, prev=None)
        done = {}
        active = list(genes)
        while active:
            nodes_l, lw_l = [], []
            for g in active:
                wk = work[g]
                if wk["n"].sum() > ma.REF_MAX_PANELS:
                    raise RuntimeError(f"reference integrator: gene {g} did not agree within "
                                       f"{ma.REF_MAX_PANELS} panels")
                nodes, logwt = ma._composite_gl(wk["edges"], wk["n"])
                nodes_l.append(nodes)
                lw_l.append(logwt)
            qg = np.concatenate([np.full(len(nd), g, np.int64) for g, nd in zip(active, nodes_l)])
            hv_all = self.h_multi(qg, np.concatenate(nodes_l), 0)[0]
            off, nxt = 0, []
            for g, nodes, logwt in zip(active, nodes_l, lw_l):
                wk = work[g]
                hv = hv_all[off:off + len(nodes)]
                off += len(nodes)
                ell = float(logsumexp(logwt + hv)) - 0.5 * np.log(2.0 * np.pi * v)
                state.events["max_reference_panels"] = max(state.events["max_reference_panels"],
                                                           int(wk["n"].sum()))
                if wk["prev"] is not None and abs(ell - wk["prev"]) <= ma.REF_AGREE:
                    done[g] = dict(ell=ell, nodes=nodes, logwt=logwt, h=hv, panels=int(wk["n"].sum()))
                else:
                    wk["prev"] = ell
                    wk["n"] = 2 * wk["n"]
                    nxt.append(g)
            active = nxt
        return done

    # ---------------------------------------------------------------- model_a._gene_node_terms on the device
    def gene_node_terms(self, g, r, order):
        """Fixed-node derivatives for gene g on its reference nodes (same rule as model_a._gene_node_terms):
        posterior-weighted row derivatives, E[dh/dv] and, at order 2, the between-node covariance of
        dh_k/d(theta, v) on the gene's own columns. Returns device tensors for the row terms."""
        L, v = self.L, self.v
        s0, m = int(self.gstart_np[g]), int(self.gcount_np[g])
        ridx = self.g_rows[s0:s0 + m]
        wg, Sg = self.w[ridx], self.St[ridx]
        A = r["logwt"] + r["h"]
        pi = np.exp(A - logsumexp(A))
        keep = pi >= ma.PI_SKIP
        nodes, pi = r["nodes"][keep], pi[keep]
        qv = nodes ** 2 / (2.0 * v ** 2)
        z = lambda: torch.zeros(m, dtype=self.dt, device=self.dev)  # noqa: E731
        t = dict(E1=z(), E2=z(), Eq_v=float(pi @ qv))
        if order >= 2:
            cols = np.unique(np.concatenate([[L.alpha[g], L.kappa], L.delta, [L.omega_S, L.r[g], L.rho_S,
                                                                               L.rho_b, L.rho_d, L.v]]))
            cols = cols.astype(np.int64)
            pos_v = int(np.searchsorted(cols, L.v))
            cols_t = self.i(cols)
            J1d = self.J1rows(ridx)[:, cols_t]
            J2d = self.J2d[ridx][:, cols_t]
            nc = len(cols)
            t.update(E11=z(), E12=z(), E22=z(), cols=cols,
                     hvv=float(pi @ (-nodes ** 2 / v ** 3)) + 1.0 / (2.0 * v ** 2))
            S1 = torch.zeros(nc, dtype=self.dt, device=self.dev)
            S2 = torch.zeros((nc, nc), dtype=self.dt, device=self.dev)
        kc = max(1, ma.CHUNK_ROWS // max(m, 1))
        kc = max(1, min(kc, self.budget[order] // max(int(self.gcost[g]), 1)))
        x1g = self.x1b[ridx]
        for s in range(0, len(nodes), kc):
            uu = self.f(nodes[s:s + kc])
            k = int(uu.numel())
            p = self.f(pi[s:s + kc])
            x1 = (x1g[None, :] + uu[:, None] * Sg[None, :]).reshape(-1)
            re = self.row_eval(ridx.repeat(k), x1, order)
            g1, g2 = re["g1"].reshape(k, m), re["g2"].reshape(k, m)
            t["E1"] += p @ g1
            t["E2"] += p @ g2
            if order >= 2:
                t["E11"] += p @ re["g11"].reshape(k, m)
                t["E12"] += p @ re["g12"].reshape(k, m)
                t["E22"] += p @ re["g22"].reshape(k, m)
                Dk = (g1 * wg) @ J1d + (g2 * wg) @ J2d
                Dk[:, pos_v] += self.f(qv[s:s + kc])
                S1 += p @ Dk
                S2 += (Dk * p[:, None]).T @ Dk
        if order >= 2:
            t["cov"] = (S2 - torch.outer(S1, S1)).cpu().numpy()
        t["ridx"] = ridx
        return t

    # ---------------------------------------------------------------- section 4 (model_a.loglik4)
    def loglik4(self, theta, order=1, nodes=20, detail=False, state=None, final=False, force=None,
                all_reference=False):
        L = self.L
        v = float(theta[L.v])
        G = self.G
        if v == 0.0:
            if state is not None:
                state.last_choice = None
            self.set_theta(theta)
            if order == 0:
                re = self.row_eval_all(self.x1b, 0)
                res = float((self.w * re["logL"]).sum())
                if detail:
                    # model_a.gene_loglik3 is the CPU code; the per-gene split is evaluated here on the device
                    ell_g = segsum(self.w * re["logL"], self.gene, G).cpu().numpy()
                    return res, dict(ell_g=ell_g, integrator=None)
                return res
            re = self.row_eval_all(self.x1b, 2)
            w = self.w
            ll = float((w * re["logL"]).sum())
            if order == 1:
                grad = self.assemble1(w * re["g1"], w * re["g2"]).cpu().numpy()
                grad[L.v] = self._boundary_from(re)
                return ll, grad
            grad, H = self.assemble2(w * re["g1"], w * re["g2"], w * re["g11"], w * re["g12"], w * re["g22"])
            grad, H = grad.cpu().numpy(), H.cpu().numpy()
            grad[L.v] = self._boundary_from(re)
            H[L.v, :] = np.nan
            H[:, L.v] = np.nan
            return ll, grad, H
        if state is None:
            state = ma.Sec4State(L.G)
        ev_count = state.events
        ev_count["evaluations"] += 1
        self.set_theta(theta)
        has = self.has
        md = self.modes(state, final or all_reference)
        choice = md["choice"].copy() if force is None else np.asarray(force).copy()
        if all_reference:
            choice[:] = ma.REF
        choice[~has] = ma.GH
        u_hat = md["u_hat"]
        sig_gh = np.where(np.isfinite(md["sigma"]), md["sigma"], np.sqrt(v))
        run_gh = not all_reference
        ell_gh = np.full(G, np.nan)
        gk = None
        if run_gh:
            xk, wk = hermgauss(nodes)
            logW = np.log(wk) + xk ** 2
            U = u_hat[None, :] + SQRT2 * sig_gh[None, :] * xk[:, None]       # nodes x G
            gk = self._gh_nodes(self.f(U), order)
            Hk = gk["Hk"].cpu().numpy()
            A = logW[:, None] + Hk
            lse = logsumexp(A, axis=0)
            ell_gh = lse + np.log(SQRT2 * sig_gh) - 0.5 * np.log(2.0 * np.pi * v)
        ell_g = ell_gh.copy()
        refs = {}
        if force is None and not all_reference:
            chk = [int(g) for g in np.flatnonzero(md["checked"] & (choice == ma.GH) & has)]
            rr = self.reference_batch(chk, md, state)
            for g in chk:
                ev_count["reference_checks"] += 1
                if abs(rr[g]["ell"] - ell_gh[g]) > ma.GH_REF_TOL:
                    choice[g] = ma.REF
                    ev_count["gh_ref_mismatch"] += 1
                    refs[g] = rr[g]
        need = [int(g) for g in np.flatnonzero((choice == ma.REF) & has) if int(g) not in refs]
        if need:
            refs.update(self.reference_batch(need, md, state))
        for g in np.flatnonzero((choice == ma.REF) & has):
            ell_g[g] = refs[int(g)]["ell"]
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
        nr = self.nr
        z = lambda: torch.zeros(nr, dtype=self.dt, device=self.dev)   # noqa: E731
        E1, E2 = z(), z()
        if order >= 2:
            E11, E12, E22 = z(), z(), z()
        Eq_v, hvv = np.zeros(G), np.zeros(G)
        if run_gh:
            pi = np.exp(A - lse[None, :])                               # nodes x G
            P = self.f(pi)[:, self.gene]                                # nodes x rows
            E1 = (P * gk["g1"]).sum(0)
            E2 = (P * gk["g2"]).sum(0)
            q_v = U ** 2 / (2.0 * v ** 2)
            Eq_v = np.sum(pi * q_v, axis=0)
            if order >= 2:
                E11 = (P * gk["g11"]).sum(0)
                E12 = (P * gk["g12"]).sum(0)
                E22 = (P * gk["g22"]).sum(0)
                hvv = np.sum(pi * (-U ** 2 / v ** 3), axis=0) + 1.0 / (2.0 * v ** 2)
        terms = {}
        for g, r in refs.items():
            t = self.gene_node_terms(g, r, order)
            ridx = t["ridx"]
            E1[ridx], E2[ridx] = t["E1"], t["E2"]
            Eq_v[g] = t["Eq_v"]
            if order >= 2:
                E11[ridx], E12[ridx], E22[ridx] = t["E11"], t["E12"], t["E22"]
                hvv[g] = t["hvv"]
            terms[g] = t
        dv = float(np.sum((Eq_v - 1.0 / (2.0 * v))[has]))
        w = self.w
        if order == 1:
            grad = self.assemble1(w * E1, w * E2).cpu().numpy()
            grad[L.v] = dv
            return (ll, grad, info) if detail else (ll, grad)
        grad, H = self.assemble2(w * E1, w * E2, w * E11, w * E12, w * E22)
        if run_gh:
            gh_gene = ((choice == ma.GH) & has).astype(float)
            gh_row = self.f(gh_gene)[self.gene]
            blocks = []
            for k in range(nodes):
                sq = torch.sqrt(P[k]) * gh_row
                Mk = self._gene_rows_to_cols(w * sq * (gk["g1"][k] - E1), w * sq * (gk["g2"][k] - E2))
                Mk[:, L.v] += self.f(np.sqrt(pi[k]) * (q_v[k] - Eq_v) * gh_gene)
                blocks.append(Mk)
            M = torch.cat(blocks, 0)
            H += M.T @ M
            del M, blocks
        H = H.cpu().numpy()
        grad = grad.cpu().numpy()
        grad[L.v] = dv
        for g, t in terms.items():
            H[np.ix_(t["cols"], t["cols"])] += t["cov"]
        H[L.v, L.v] += float(np.sum(hvv[has]))
        return (ll, grad, H, info) if detail else (ll, grad, H)

    def _gene_rows_to_cols(self, a_, b_):
        """Agg @ (diag(a) J1 + diag(b) J2) as a dense G x P matrix."""
        L, G = self.L, self.G
        gidx = torch.arange(G, device=self.dev)
        M = torch.zeros((G, self.P), dtype=self.dt, device=self.dev)
        M[gidx, self.ix_alpha] = segsum(a_ * self.q, self.gene, G)
        M[:, L.kappa] = segsum(a_ * self.al * self.St, self.gene, G)
        M[:, self.ix_delta] = torch.zeros((G, L.p), dtype=self.dt, device=self.dev).index_add_(
            0, self.gene, (a_ * self.al)[:, None] * self.Z)
        M[:, L.omega_S] = segsum(a_ * self.s_t * self.St, self.gene, G)
        M[gidx, self.ix_r] = segsum(b_, self.gene, G)
        M[:, L.rho_S] = segsum(b_ * self.St, self.gene, G)
        M[:, L.rho_b] = segsum(b_ * self.b, self.gene, G)
        M[:, L.rho_d] = segsum(b_ * self.d, self.gene, G)
        return M

    def _gh_nodes(self, U_t, order):
        """Row terms at the GH nodes for all rows: Hk (nodes x G) and per-node row derivatives."""
        K = U_t.shape[0]
        G, nr = self.G, self.nr
        per_node = int(self.rowcost.sum())
        kc = max(1, self.budget[order] // max(per_node, 1))
        Hk = torch.empty((K, G), dtype=self.dt, device=self.dev)
        keys = ([] if order == 0 else ["g1", "g2"]) + (["g11", "g12", "g22"] if order >= 2 else [])
        out = {k: torch.empty((K, nr), dtype=self.dt, device=self.dev) for k in keys}
        Urow = U_t[:, self.gene]                                       # nodes x rows
        for s in range(0, K, kc):
            e = min(K, s + kc)
            m = e - s
            ridx = self.arange_nr.repeat(m)
            x1 = (self.x1b[None, :] + Urow[s:e] * self.St[None, :]).reshape(-1)
            re = self.row_eval(ridx, x1, order)
            ll = self.w[None, :] * re["logL"].reshape(m, nr)
            for kk in range(m):
                Hk[s + kk] = segsum(ll[kk], self.gene, G) - U_t[s + kk] ** 2 / (2.0 * self.v)
            for k in keys:
                out[k][s:e] = re[k].reshape(m, nr)
        out["Hk"] = Hk
        return out


def environment():
    """Device, driver and library versions for run headers."""
    import ctypes
    import platform
    out = dict(host=platform.node(), torch=torch.__version__, torch_cuda=torch.version.cuda,
               cudnn=torch.backends.cudnn.version() if torch.backends.cudnn.is_available() else None,
               deterministic=torch.are_deterministic_algorithms_enabled(),
               cublas_workspace_config=os.environ.get("CUBLAS_WORKSPACE_CONFIG"), device=DEVICE)
    try:
        libc = ctypes.CDLL("libc.so.6")
        libc.gnu_get_libc_version.restype = ctypes.c_char_p
        out["glibc"] = libc.gnu_get_libc_version().decode()
    except Exception as err:  # noqa: BLE001
        out["glibc"] = repr(err)
    if DEVICE.startswith("cuda") and torch.cuda.is_available():
        p = torch.cuda.get_device_properties(0)
        out["gpu"] = dict(name=p.name, cc=f"{p.major}.{p.minor}", mem_gb=p.total_memory / 1e9,
                          sms=p.multi_processor_count)
        try:
            import subprocess
            out["nvidia_driver"] = subprocess.run(
                ["nvidia-smi", "--query-gpu=driver_version", "--format=csv,noheader"], capture_output=True,
                text=True, timeout=30).stdout.strip()
        except Exception as err:  # noqa: BLE001
            out["nvidia_driver"] = repr(err)
    return out
