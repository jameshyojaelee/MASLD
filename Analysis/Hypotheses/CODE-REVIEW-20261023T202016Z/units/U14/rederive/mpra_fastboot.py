"""Vectorised paired block bootstrap of the Fisher-mean Spearman score, for mpra_7b.py and mpra_ladder.py.

A block-bootstrap resample is the element set with multiplicities m_e (an element enters once per draw
of its block). Spearman on the expanded vector uses average ranks; with multiplicities, a group of tied
values with total multiplicity c, preceded by s copies of smaller values, has average rank s + (c + 1) / 2.
The correlation is the multiplicity-weighted Pearson correlation of those ranks, which equals scipy's
rankdata-then-Pearson on the expanded vector (checked against it in check_against_rankdata).
"""
import numpy as np
from scipy.stats import rankdata


def draw_block_multiplicities(rng, n_blocks, B):
    """(B, n_blocks) counts; one integers() call per draw, as in the per-draw loop."""
    counts = np.zeros((B, n_blocks))
    for b in range(B):
        counts[b] = np.bincount(rng.integers(0, n_blocks, n_blocks), minlength=n_blocks)
    return counts


def _ranks(vals, M):
    uniq, inv = np.unique(vals, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    starts = np.r_[0, np.flatnonzero(np.diff(inv[order])) + 1]
    Cg = np.add.reduceat(M[:, order], starts, axis=1)
    S = np.cumsum(Cg, axis=1) - Cg
    return (S + (Cg + 1) / 2)[:, inv]


def fisher_scores(Y, P, M):
    """Y: (2, n) outcome per context; P: (k, 2, n) predictions; M: (B, n) element multiplicities.
    Returns (B, k) tanh(mean over contexts of arctanh(Spearman))."""
    N = M.sum(1, keepdims=True)
    mid = (N + 1) / 2
    z = np.zeros((M.shape[0], P.shape[0]))
    for j in range(Y.shape[0]):
        ry = _ranks(Y[j], M) - mid
        vy = (M * ry * ry).sum(1)
        for i in range(P.shape[0]):
            rp = _ranks(P[i, j], M) - mid
            rho = (M * ry * rp).sum(1) / np.sqrt(vy * (M * rp * rp).sum(1))
            z[:, i] += np.arctanh(rho)
    return np.tanh(z / Y.shape[0])


def check_against_rankdata(Y, P, M, elem_rows):
    """Recompute the first rows of M with rankdata on expanded index vectors; return max abs difference."""
    fast = fisher_scores(Y, P, M)
    worst = 0.0
    for b in range(M.shape[0]):
        idx = np.repeat(np.arange(M.shape[1]), M[b].astype(int))
        zs = np.zeros(P.shape[0])
        for j in range(Y.shape[0]):
            ry = rankdata(Y[j, idx]); ry -= ry.mean()
            rp = rankdata(P[:, j, idx], axis=1); rp -= rp.mean(1, keepdims=True)
            zs += np.arctanh((rp @ ry) / np.sqrt((rp ** 2).sum(1) * (ry ** 2).sum()))
        worst = max(worst, float(np.abs(np.tanh(zs / Y.shape[0]) - fast[b]).max()))
    return worst
