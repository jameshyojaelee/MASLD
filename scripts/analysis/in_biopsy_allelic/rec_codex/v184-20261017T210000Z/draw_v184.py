"""REC v1.8.4: P3 conditional draw with the genotype class order the fixture specifies.

Defect in likelihood_v1.conditional_draw (sha256 41ebc253...): it takes the inverse CDF of
column-0 uniforms over the class index order (0, 1, 2). Row.prior (likelihood_v1.py:82) stores
index 1 as hom tag-REF and index 2 as hom tag-ALT for orientation +1, but swaps them
(prior[[0, 2, 1]]) for orientation -1, so for those rows the draw order was (H, A, R).
B-FIX globals.json p3_uniform_use: "column 0 picks the genotype class by inverse CDF over
(H, R, A) in that order with weights P(G)C_G(n); column 1 picks a by inverse CDF over K_c(n) in
increasing a". The generator labels class 0 H, 1 R (hom tag-REF), 2 A (hom tag-ALT) for every
orientation. This version walks the classes in genotype order (H, R, A) for every row.
Class weights, components and the count step are unchanged from likelihood_v1.
"""
from dataclasses import replace

import numpy as np
from scipy.special import logsumexp

from likelihood_v1 import bb_terms, fixed_bb

# Index positions of (H, R, A) in Row.prior / the component stack, by orientation.
GENOTYPE_ORDER = {1: (0, 1, 2), -1: (0, 2, 1)}


def class_weights(row, eta, linear_r):
    k = row.allowed
    lp = bb_terms(k, row.n, eta - row.orientation * row.omega, linear_r)[0]
    parts = np.stack((lp, fixed_bb(k, row.n, row.e, row.phi), fixed_bb(k, row.n, 1 - row.e, row.phi)))
    with np.errstate(divide='ignore'):
        class_logs = np.log(row.prior) + logsumexp(parts, axis=1)
    return k, parts, np.exp(class_logs - logsumexp(class_logs))


def conditional_draw_v184(row, eta, linear_r, uniform_class, uniform_count):
    """Return (drawn row, class index chosen, class index the v1 order would have chosen)."""
    k, parts, class_prob = class_weights(row, eta, linear_r)
    order = GENOTYPE_ORDER[int(row.orientation)]
    position = min(2, int(np.searchsorted(np.cumsum(class_prob[list(order)]), uniform_class, side='left')))
    g = order[position]
    g_v1 = min(2, int(np.searchsorted(np.cumsum(class_prob), uniform_class, side='left')))
    p = np.exp(parts[g] - logsumexp(parts[g]))
    chosen = min(len(k) - 1, int(np.searchsorted(np.cumsum(p), uniform_count, side='left')))
    return replace(row, a=int(k[chosen])), g, g_v1
