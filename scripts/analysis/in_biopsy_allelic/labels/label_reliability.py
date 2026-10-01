#!/usr/bin/env python3
"""Model A, layer B-LABEL: test-retest reliability of recorded fibrosis stage, from labels only (spec v1 5.2).

Inputs
  --gse193066-placement  gse193066_donor_keyed_placement.tsv: sample_id (run) with a "1st biopsy" and a
                         "2nd biopsy" row, source_fibrosis_stage 0-4 (58 repeat-biopsy donors).
  --gse193066-source     kleiner-reliability-20260901T122243Z: results/pairs.tsv (the same 58 pairs, checked)
                         and results/kleiner_reliability.json (within-pair SD 0.6433 [0.5170, 0.7600]).
  --gse193066-series     GSE193066_series_matrix.txt.gz: the 58 pairs are re-derived from the GEO record.
  --gse193066-crosswalk  source participant crosswalk (run_id -> participant_token), checked against the placement.
  --pxd-source           kleiner-reliability-replication-pxd051911-20260901T211358Z: results/pairs.tsv
                         arms A1 (V1->V2, 27 pairs, PRIMARY in its prespec) and A2 (V1->V3 without V2, 30 pairs,
                         design-matched), and results/kleiner_replication.json. Re-derived from --pxd-meta.
  --frozen-crosswalk     seal frozen_crosswalk.tsv, checked against frozen_crosswalk.sha256 in the same directory
                         before it is read. unit_library marks the one library per development or sealed
                         individual; S (0 = F0-1, 1 = F2, 2 = F3-4), C (source control).
  --kinship-pairs        identity/kinship_related_pairs.tsv (same_individual pairs; checked against the crosswalk).
  --metadata             unified_metadata.csv: fibrosis_stage 0-4 (Kleiner cohorts) and condition (GSE213621).
  --design-individuals   Gate A design_individuals.tsv (job 1208): z columns without the axes, cluster_id, S.
  --ancestry             restricted library_ancestry_pcs.tsv: axis_EUR_AFR, axis_EUR_EAS (read, never written).
  --r2-moments           r2_moments_observed.json (job 1290): observed raw and adjusted R2 of S on z per F cohort.
  --p2-resamples         p2_resample_development.tsv.gz (job 1290): the fixed P2 resamples of F (spec 5.2).
  --r2-by-resample       r2_by_resample.tsv.gz (job 1290): adjusted R2 per P2 resample and cohort (checked).
  --codex-lambda         optional: Codex C8 lambda_558.json, compared at the point.
Target (spec v1 1.1, 5.2)
  F: unit_library True, run_role development, C = 0, S observed, cohort in GSE130970, GSE135251, GSE162694,
  GSE213621, GSE240729 (558). F_sealed (237, adds GSE174478) is reported for the sealed run (spec 6), marginal
  lambda only. Each cohort's latent stage distribution is fitted from that set's own labels.
Sources and arms (pair tables; fitted per arm, deposits never pooled)
  GSE193066 (58 pairs), PXD051911_A1 (27 pairs)   the two spec 5.2 sources.
  PXD051911_A2, PXD051911_joint                   sensitivities (joint: one var_e shared by A1 and A2).
  GSE193066_same29, GSE193066_drop11              provisional identity sensitivities.
  GSE213621_pairs25                               25 genotype-same library pairs, S only, symmetric: a selected set.
Specs per arm
  Y5: latent Z_v = T + e_v (+ delta at the repeat visit), T ~ N(mu, var_T), e_v ~ N(0, var_e); recorded stage
      Y = k when k - 0.5 < Z <= k + 0.5 (cutpoints 0.5 ... 3.5 fixed). MLE on the 5 x 5 table.
      Y5_rho_nonneg: the same with var_T >= 0. S3: the same model on the S-collapsed 3 x 3 table.
  Y5_het (sources only): error SD exp(a + b (T - 2)) depends on the latent level; b = 0 is Y5_rho_nonneg.
lambda
  lambda_c = Var(E[S|T]) / Var(S) at cohort c's latent stage distribution with the arm's var_e; the slope
  attenuation for kappa per unit of the true score E[S|T]. lambda_cond_c = (lambda_c - R2_c) / (1 - R2_c), R2_c the
  adjusted within-cohort R2 of S on the spec z columns and ancestry axes (v1.2 items 6, 7), floored at 0.
  Carried spec per source: the smaller lambda of Y5 and S3 at the point (= the larger var_e), fixed across draws.
  Pools: headcount mean of lambda_c and lambda_cond_c (the spec pool with I_c replaced by n_c; B-MODEL re-pools
  with I_c from the per-cohort columns).
Bootstrap and fit checks
  Seed 20260926. 2,000 draws: pairs resampled within arm; at draw b the F target is P2 resample b and R2 is
  recomputed on it; F_sealed individuals resampled within cohort. Parametric-bootstrap goodness of fit per arm
  and spec (1,000 draws, refit each). Latent-shape sensitivity at the sources: grid NPMLE and LP bounds.
Outputs (--out; label counts and aggregate summaries only, no identifiers)
  label_reliability.json, lambda_source_draws.tsv.gz (the P2 consumer file), carried_spec.tsv,
  lambda_by_cohort.tsv, het_sensitivity.tsv, latent_shape_sensitivity.tsv, goodness_of_fit.tsv, pair_tables.tsv,
  target_stage_counts.tsv, bootstrap_draws.npz, bootstrap_draws.tsv.gz.
  --restricted-out: pairs_used.tsv and target_participants.tsv (one row per person). No file is overwritten.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog, minimize
from scipy.special import owens_t
from scipy.stats import chi2, multivariate_normal, norm

SEED = 20260926
N_BOOT = 2000
N_GOF = 1000
GOF_ALPHA = 0.05
SD_GUARD_TOL = 0.02
INVALID_LAMBDA = 0.05          # spec 5.2 and v1.2 item 1: a draw at or below this is invalid
CHUNK = 250
CUTS_Y = np.array([0.5, 1.5, 2.5, 3.5])
CUTS_S = np.array([1.5, 2.5])
Y_VALUES = np.arange(5.0)
S_VALUES = np.arange(3.0)
CUTS = {"Y5": CUTS_Y, "S3": CUTS_S}
VALS = {"Y5": Y_VALUES, "S3": S_VALUES}
LATENT_GRID = np.linspace(-4.0, 8.0, 2401)


def bounds_for(fix_delta=False, rho_nonneg=False):
    return [(-6.0, 10.0), (np.log(0.05), np.log(20.0)), (0.0, 3.8) if rho_nonneg else (-3.8, 3.8),
            (0.0, 0.0) if fix_delta else (-4.0, 4.0)]


PAIR_BOUNDS = bounds_for()
SINGLE_BOUNDS = PAIR_BOUNDS[:2]

F_COHORTS = ["GSE130970", "GSE135251", "GSE162694", "GSE213621", "GSE240729"]    # sorted: the P2 strata order
F_COUNTS = {"GSE130970": 54, "GSE135251": 160, "GSE162694": 83, "GSE213621": 208, "GSE240729": 53}   # spec 1.1
SEALED_COHORTS = F_COHORTS + ["GSE174478"]
SEALED_COUNTS = {"GSE130970": 13, "GSE135251": 42, "GSE162694": 22, "GSE213621": 55, "GSE240729": 13,
                 "GSE174478": 92}
SETS = {"F": ("development", F_COHORTS), "F_sealed": ("sealed", SEALED_COHORTS)}
Z_COLS = ["z_b_log_e", "z_d_dup", "z_ffpe", "z_sex_female", "z_age10", "axis_EUR_AFR", "axis_EUR_EAS", "z_c"]
GSE213621_CONDITION = {"Fibrosis_F0F1": 0, "Fibrosis_F2": 1, "Fibrosis_F3F4": 2}
PXD_ARMS = {"PXD051911_A1": "A1_initial_to_at_surgery",
            "PXD051911_A2": "A2_initial_to_followup_no_surgery_sample"}
KLEINER_LEVELS = {"F0": 0, "F1": 1, "F2": 2, "F3": 3, "F4": 4}

SPECS = {"Y5": dict(scale="Y5", fix_delta=False, rho_nonneg=False),
         "Y5_rho_nonneg": dict(scale="Y5", fix_delta=False, rho_nonneg=True),
         "S3": dict(scale="S3", fix_delta=False, rho_nonneg=False),
         "S3_sym": dict(scale="S3", fix_delta=True, rho_nonneg=False),
         "S3_sym_rho_nonneg": dict(scale="S3", fix_delta=True, rho_nonneg=True)}
KLEINER_SPECS = ["Y5", "Y5_rho_nonneg", "S3"]
PAIRS25_SPECS = ["S3_sym", "S3_sym_rho_nonneg"]
GOF_SPECS = {"Y5", "S3", "S3_sym"}
SOURCES = ["GSE193066", "PXD051911_A1"]
SOURCE_SPECS = ["Y5", "S3"]
ARM_ORDER = ["GSE193066", "PXD051911_A1", "PXD051911_A2", "PXD051911_joint",
             "GSE193066_same29", "GSE193066_drop11", "GSE213621_pairs25"]
ARM_ROLE = {"GSE193066": "source (spec v1 5.2)",
            "PXD051911_A1": "source (spec v1 5.2; PRIMARY arm of the PXD051911 replication prespec)",
            "PXD051911_A2": "sensitivity (DESIGN-MATCHED arm of the PXD051911 replication prespec)",
            "PXD051911_joint": "sensitivity (A1 and A2 with one shared var_e)",
            "GSE193066_same29": "sensitivity, provisional identity (29 pairs called one person)",
            "GSE193066_drop11": "sensitivity, provisional identity (47 pairs; 11 pairs called two people dropped)",
            "GSE213621_pairs25": "sensitivity, selected set of 25 genotype-same library pairs, S only, symmetric"}
PRIMARY_Y5_ARMS = ["GSE193066", "PXD051911_A1", "PXD051911_A2"]   # the recorded arms: SD agreement is a guard
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:8.1f} s] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------- probability model
def bvn_cdf(h, k, rho):
    """P(X <= h, Y <= k) for a standard bivariate normal with correlation rho (Owen's T form)."""
    h, k = np.broadcast_arrays(np.asarray(h, float), np.asarray(k, float))
    h = np.where(h == 0.0, 1e-12, h)
    k = np.where(k == 0.0, 1e-12, k)
    s = np.sqrt(1.0 - rho * rho)
    beta = np.where(h * k > 0, 0.0, 0.5)
    return (0.5 * norm.cdf(h) + 0.5 * norm.cdf(k)
            - owens_t(h, (k - rho * h) / (h * s)) - owens_t(k, (h - rho * k) / (k * s)) - beta)


def cell_probs(cuts, mu1, mu2, sd, rho):
    """Joint probabilities of two ordinal reads cut from a bivariate normal; (K+1) x (K+1)."""
    u, v = (cuts - mu1) / sd, (cuts - mu2) / sd
    m = len(cuts)
    F = np.zeros((m + 2, m + 2))
    F[1:m + 1, 1:m + 1] = bvn_cdf(u[:, None], v[None, :], rho)
    F[m + 1, 1:m + 1] = norm.cdf(v)
    F[1:m + 1, m + 1] = norm.cdf(u)
    F[m + 1, m + 1] = 1.0
    P = F[1:, 1:] - F[:-1, 1:] - F[1:, :-1] + F[:-1, :-1]
    return np.clip(P, 1e-300, None)


def pair_moments(P, vals):
    """Means, variances, covariance and E[(v2 - v1)^2] of a joint table over scored levels."""
    p1, p2 = P.sum(1), P.sum(0)
    m1, m2 = p1 @ vals, p2 @ vals
    v1, v2 = p1 @ vals ** 2 - m1 ** 2, p2 @ vals ** 2 - m2 ** 2
    cov = vals @ P @ vals - m1 * m2
    d2 = float(((vals[None, :] - vals[:, None]) ** 2 * P).sum())
    return dict(m1=m1, m2=m2, v1=v1, v2=v2, cov=cov, d2=d2)


def single_occasion_corr(cuts, vals, mu, sd, rho):
    """Corr of two independent reads at one occasion (no shift): the reliability of the scored label."""
    mo = pair_moments(cell_probs(cuts, mu, mu, sd, rho), vals)
    return mo["cov"] / mo["v1"], mo


def unpack(theta):
    mu, log_sd, z, delta = theta
    return mu, float(np.exp(log_sd)), float(np.tanh(z)), delta


def scale_of(N):
    return "Y5" if np.shape(N)[0] == 5 else "S3"


def pair_negll(theta, N):
    mu, sd, rho, delta = unpack(theta)
    return -(N * np.log(cell_probs(CUTS[scale_of(N)], mu, mu + delta, sd, rho))).sum()


def pair_table(y1, y2, k=5):
    N = np.zeros((k, k))
    np.add.at(N, (np.asarray(y1, int), np.asarray(y2, int)), 1)
    return N


def harmonize(y):
    """Kleiner stage 0-4 to S: 0 = F0-1, 1 = F2, 2 = F3-4."""
    y = np.asarray(y, float)
    return (y >= 2).astype(float) + (y >= 3)


def table_for(y1, y2, input_scale, spec_scale):
    if input_scale == "Y5" and spec_scale == "S3":
        y1, y2 = harmonize(y1), harmonize(y2)
    return pair_table(y1, y2, 5 if spec_scale == "Y5" else 3)


def moment_start(N):
    vals = VALS[scale_of(N)]
    p1, p2 = N.sum(1) / N.sum(), N.sum(0) / N.sum()
    m1, m2 = p1 @ vals, p2 @ vals
    v = 0.5 * (p1 @ vals ** 2 - m1 ** 2 + p2 @ vals ** 2 - m2 ** 2)
    cov = vals @ (N / N.sum()) @ vals - m1 * m2
    r = cov / v if v > 0 else 0.0
    return np.array([m1, np.log(max(np.sqrt(v), 0.3)), np.arctanh(np.clip(r, -0.9, 0.9)), m2 - m1])


def at_bound(x, bounds, tol=1e-6):
    return any((lo != hi) and (abs(xi - lo) < tol or abs(xi - hi) < tol) for xi, (lo, hi) in zip(x, bounds))


def fit_pairs(N, warm=None, fix_delta=False, rho_nonneg=False):
    """MLE of (mu, log sd_Z, atanh rho, delta) on a 5 x 5 or 3 x 3 pair table; best of a moment and a warm start."""
    bounds = bounds_for(fix_delta, rho_nonneg)
    starts = [moment_start(N)] + ([np.asarray(warm, float)] if warm is not None else [])
    best = None
    for s in starts:
        s = np.array([np.clip(v, lo, hi) for v, (lo, hi) in zip(s, bounds)])
        r = minimize(pair_negll, s, args=(N,), method="L-BFGS-B", bounds=bounds)
        if best is None or r.fun < best.fun:
            best = r
    return dict(theta=best.x, negll=float(best.fun), converged=bool(best.success),
                at_bound=at_bound(best.x, bounds))


# ---- joint fit: one var_e shared by two arms of one deposit
def joint_bounds(rho_nonneg=False):
    z = (0.0, 3.8) if rho_nonneg else (-3.8, 3.8)
    return [(-6.0, 10.0), (-6.0, 10.0), z, z, (-4.0, 4.0), (-4.0, 4.0), (np.log(1e-4), np.log(100.0))]


def joint_arm_params(theta):
    """theta = (mu_1, mu_2, atanh rho_1, atanh rho_2, delta_1, delta_2, log var_e); sd_a^2 = var_e / (1 - rho_a)."""
    var_e = float(np.exp(theta[6]))
    out = []
    for a in (0, 1):
        rho = float(np.tanh(theta[2 + a]))
        out.append((float(theta[a]), float(np.sqrt(var_e / (1.0 - rho))), rho, float(theta[4 + a])))
    return out, var_e


def joint_negll(theta, Ns):
    arms, _ = joint_arm_params(theta)
    return sum(-(N * np.log(cell_probs(CUTS[scale_of(N)], mu, mu + d, sd, rho))).sum()
               for N, (mu, sd, rho, d) in zip(Ns, arms))


def joint_start(Ns, var_e):
    s = [moment_start(N) for N in Ns]
    th = [s[0][0], s[1][0], 0.0, 0.0, s[0][3], s[1][3], np.log(max(var_e, 1e-3))]
    for a in (0, 1):
        th[2 + a] = np.arctanh(np.clip(1.0 - var_e / np.exp(2 * s[a][1]), -0.9, 0.9))
    return np.array(th)


def fit_joint(Ns, warm=None, rho_nonneg=False, var_e_starts=()):
    bounds = joint_bounds(rho_nonneg)
    s = [moment_start(N) for N in Ns]
    ve0 = float(np.mean([(1 - np.tanh(x[2])) * np.exp(2 * x[1]) for x in s]))
    starts = [joint_start(Ns, v) for v in (max(ve0, 0.02), *var_e_starts)]
    if warm is not None:
        starts.append(np.asarray(warm, float))
    best = None
    for st in starts:
        st = np.array([np.clip(v, lo, hi) for v, (lo, hi) in zip(st, bounds)])
        r = minimize(joint_negll, st, args=(Ns,), method="L-BFGS-B", bounds=bounds)
        if best is None or r.fun < best.fun:
            best = r
    return dict(theta=best.x, negll=float(best.fun), converged=bool(best.success),
                at_bound=at_bound(best.x, bounds))


# ---- single-read target fits
def fit_single(counts, cuts=CUTS_Y):
    """Single-read ordinal probit with the fixed cutpoints: (mu_Z, sd_Z) of a cohort's stage labels."""
    counts = np.asarray(counts, float)
    vals = np.arange(len(cuts) + 1.0)

    def probs(t):
        return np.diff(np.concatenate([[0.0], norm.cdf((cuts - t[0]) / np.exp(t[1])), [1.0]]))

    def nll(t):
        return -(counts * np.log(np.clip(probs(t), 1e-300, None))).sum()

    p = counts / counts.sum()
    m = p @ vals + (cuts[0] - 0.5)
    sd0 = np.sqrt(max(p @ vals ** 2 - (p @ vals) ** 2, 0.09))
    r = minimize(nll, [m, np.log(sd0)], method="L-BFGS-B", bounds=SINGLE_BOUNDS)
    mu, sd = float(r.x[0]), float(np.exp(r.x[1]))
    fitted_p = probs(r.x)
    return dict(mu_Z=mu, sd_Z=sd, converged=bool(r.success), at_bound=at_bound(r.x, SINGLE_BOUNDS),
                fitted_counts=fitted_p * counts.sum(), fitted_p=fitted_p)


def exact_three_level_probit(counts):
    """Closed form for 3 levels and cutpoints 1.5, 2.5: the probit is exactly identified."""
    p = np.asarray(counts, float) / np.sum(counts)
    z0, z1 = norm.ppf(p[0]), norm.ppf(p[0] + p[1])
    sd = (CUTS_S[1] - CUTS_S[0]) / (z1 - z0)
    return CUTS_S[0] - sd * z0, sd


# ---- heteroscedastic sensitivity: the read error SD depends on the latent level
HET_Z = np.linspace(-7.0, 7.0, 801)
HET_W = norm.pdf(HET_Z) / norm.pdf(HET_Z).sum()
HET_CENTRE = 2.0
HET_BOUNDS = [(-6.0, 10.0), (np.log(0.02), np.log(20.0)), (np.log(0.02), np.log(5.0)), (-3.0, 3.0), (-4.0, 4.0)]


def het_read_probs(cuts, t, a, b):
    """P(read = k | latent level t) when the read error SD is exp(a + b (t - 2)); rows t, columns k."""
    sd = np.exp(a + b * (t - HET_CENTRE))
    cdf = norm.cdf((cuts[None, :] - t[:, None]) / sd[:, None])
    return np.diff(np.hstack([np.zeros((len(t), 1)), cdf, np.ones((len(t), 1))]), axis=1)


def het_cell_probs(cuts, mu, sd_T, a, b, delta):
    """Two reads, conditionally independent given T ~ N(mu, sd_T^2); the repeat read is at T + delta."""
    t = mu + sd_T * HET_Z
    P = (HET_W[:, None] * het_read_probs(cuts, t, a, b)).T @ het_read_probs(cuts, t + delta, a, b)
    return np.clip(P, 1e-300, None)


def het_negll(theta, N):
    mu, lsd, a, b, delta = theta
    return -(N * np.log(het_cell_probs(CUTS[scale_of(N)], mu, np.exp(lsd), a, b, delta))).sum()


def het_from_hom(theta_hom):
    """Y5 (mu, log sd, atanh rho, delta) with rho >= 0 as the b = 0 point of the heteroscedastic model."""
    mu, sd, rho, delta = unpack(theta_hom)
    rho = float(np.clip(rho, 0.02, 0.98))
    return np.array([mu, np.log(np.sqrt(rho) * sd), np.log(np.sqrt(1 - rho) * sd), 0.0, delta])


def fit_het(N, starts):
    best = None
    for st in starts:
        st = np.array([np.clip(v, lo, hi) for v, (lo, hi) in zip(st, HET_BOUNDS)])
        r = minimize(het_negll, st, args=(N,), method="L-BFGS-B", bounds=HET_BOUNDS)
        if best is None or r.fun < best.fun:
            best = r
    return dict(theta=best.x, negll=float(best.fun), converged=bool(best.success),
                at_bound=at_bound(best.x, HET_BOUNDS))


def het_fit_single(counts, cuts, a, b, mu0, sd0):
    """A cohort's latent (mu_T, sd_T) from its single-read labels under the arm's error model (a, b)."""
    counts = np.asarray(counts, float)

    def nll(v):
        p = HET_W @ het_read_probs(cuts, v[0] + np.exp(v[1]) * HET_Z, a, b)
        return -(counts * np.log(np.clip(p, 1e-300, None))).sum()

    s0 = np.sqrt(max(sd0 ** 2 - np.exp(2 * a), 0.05))
    r = minimize(nll, [mu0, np.log(s0)], method="L-BFGS-B", bounds=HET_BOUNDS[:2])
    return float(r.x[0]), float(np.exp(r.x[1])), bool(r.success), at_bound(r.x, HET_BOUNDS[:2])


def het_cohort_lambda(mu, sd_T, a, b):
    ps = het_read_probs(CUTS_S, mu + sd_T * HET_Z, a, b)
    g, g2 = ps @ S_VALUES, ps @ S_VALUES ** 2
    es = float(HET_W @ g)
    var, cov = float(HET_W @ g2 - es ** 2), float(HET_W @ g ** 2 - es ** 2)
    return cov / var, es, var, cov


# ---------------------------------------------------------------- summaries
def lambda_conditional(lam, r2):
    """Attenuation of the S coefficient in a model with nuisance covariates z (label error independent of z):
    (lambda - R2) / (1 - R2), R2 = within-cohort R2 of observed S on z. Spec v1.2 item 6 uses the adjusted R2."""
    lam, r2 = np.asarray(lam, float), np.asarray(r2, float)
    return (lam - r2) / (1.0 - r2)


def cohort_lambda(var_e, mu, sd):
    var_T = sd ** 2 - var_e
    lam, mo = single_occasion_corr(CUTS_S, S_VALUES, mu, sd, max(var_T, 0.0) / sd ** 2)
    return float(lam), mo, bool(var_T <= 0)


def add_pools(row, targets):
    """Headcount pools of lambda and lambda_cond, and the within-cohort variance pool, per target set."""
    for s in SETS:
        mem = [t for t in targets if t["set"] == s]
        if not mem:
            continue
        n = np.array([t["n"] for t in mem], float)
        keys = [t["key"] for t in mem]
        row[f"lambda_hc_{s}"] = float(n @ np.array([row[f"lambda_{k}"] for k in keys]) / n.sum())
        row[f"lambda_within_{s}"] = float(n @ np.array([row[f"CovS_{k}"] for k in keys]) /
                                          (n @ np.array([row[f"VarS_{k}"] for k in keys])))
        if all(f"lambdacond_{k}" in row for k in keys):
            row[f"lambdacond_hc_{s}"] = float(n @ np.array([row[f"lambdacond_{k}"] for k in keys]) / n.sum())
    return row


def target_block(var_e, targets, r2=None):
    """lambda and S moments per target cohort, lambda_cond where an R2 is given, and the pools."""
    r2 = r2 or {}
    row, n_floored = {}, 0
    for t in targets:
        k = t["key"]
        lam, mo, floored = cohort_lambda(var_e, t["mu_Z"], t["sd_Z"])
        row.update({f"lambda_{k}": lam, f"ES_{k}": float(mo["m1"]), f"VarS_{k}": float(mo["v1"]),
                    f"CovS_{k}": float(mo["cov"])})
        if k in r2:
            row[f"R2adj_{k}"] = float(r2[k])
            row[f"lambdacond_{k}"] = float(lambda_conditional(lam, r2[k]))
        n_floored += int(floored)
    row["n_target_cohorts_floored"] = n_floored
    return add_pools(row, targets)


def het_target_block(a, b, targets, r2):
    row, n_bound = {}, 0
    fts = [t for t in targets if t["set"] == "F"]
    for t in fts:
        k = t["key"]
        mu, sd_T, _, bnd = het_fit_single(t["counts"], CUTS[t["scale"]], a, b, t["mu_Z"], t["sd_Z"])
        lam, es, var, cov = het_cohort_lambda(mu, sd_T, a, b)
        row.update({f"lambda_{k}": lam, f"ES_{k}": es, f"VarS_{k}": var, f"CovS_{k}": cov, f"sdT_{k}": sd_T})
        if k in r2:
            row[f"R2adj_{k}"] = float(r2[k])
            row[f"lambdacond_{k}"] = float(lambda_conditional(lam, r2[k]))
        n_bound += int(bnd)
    row["n_target_fits_at_bound"] = n_bound
    return add_pools(row, fts)


def deposit_block(mu, sd, rho, delta, scale):
    mo = pair_moments(cell_probs(CUTS[scale], mu, mu + delta, sd, rho), VALS[scale])
    lam_in, _ = single_occasion_corr(CUTS_S, S_VALUES, mu, sd, rho)
    return dict(implied_within_pair_sd=float(np.sqrt(mo["d2"] / 2.0)), implied_mean_delta=float(mo["m2"] - mo["m1"]),
                lambda_in_deposit=float(lam_in))


def summarise_fit(theta, targets, scale="Y5", r2=None):
    mu, sd, rho, delta = unpack(theta)
    var_e, var_T = (1.0 - rho) * sd ** 2, rho * sd ** 2
    row = dict(mu=float(mu), sd_Z=sd, rho=rho, delta=float(delta), var_T=var_T, var_e=var_e,
               sd_e=float(np.sqrt(var_e)))
    row.update(deposit_block(mu, sd, rho, delta, scale))
    row.update(target_block(var_e, targets, r2))
    return row


def summarise_joint(theta, targets, scale="Y5", r2=None, names=("A1", "A2")):
    arms, var_e = joint_arm_params(theta)
    row = dict(var_e=var_e, sd_e=float(np.sqrt(var_e)), rho=float(min(a[2] for a in arms)))
    for name, (mu, sd, rho, delta) in zip(names, arms):
        row.update({f"mu_{name}": mu, f"sd_Z_{name}": sd, f"rho_{name}": rho, f"delta_{name}": delta,
                    f"var_T_{name}": rho * sd ** 2})
        row.update({f"{k}_{name}": v for k, v in deposit_block(mu, sd, rho, delta, scale).items()})
    row.update(target_block(var_e, targets, r2))
    return row


def summarise_het(theta, targets, r2):
    mu, lsd, a, b, delta = (float(v) for v in theta)
    mo = pair_moments(het_cell_probs(CUTS_Y, mu, np.exp(lsd), a, b, delta), Y_VALUES)
    row = dict(mu=mu, sd_T=float(np.exp(lsd)), a=a, b=b, delta=delta,
               **{f"sd_e_at_{str(t).replace('.', 'p')}": float(np.exp(a + b * (t - HET_CENTRE)))
                  for t in (0.5, 1.5, 2.5, 3.5)},
               implied_within_pair_sd=float(np.sqrt(mo["d2"] / 2.0)))
    row.update(het_target_block(a, b, targets, r2))
    return row


def empirical(y1, y2, input_scale="Y5"):
    y1, y2 = np.asarray(y1, float), np.asarray(y2, float)
    s1, s2 = (harmonize(y1), harmonize(y2)) if input_scale == "Y5" else (y1, y2)
    corr_s = np.corrcoef(s1, s2)[0, 1] if s1.std() > 0 and s2.std() > 0 else np.nan
    return dict(empirical_within_pair_sd=float(np.sqrt(((y2 - y1) ** 2).sum() / (2.0 * len(y1)))),
                empirical_mean_delta=float((y2 - y1).mean()), empirical_corr_S=float(corr_s),
                empirical_within_pair_sd_S=float(np.sqrt(((s2 - s1) ** 2).sum() / (2.0 * len(s1)))),
                empirical_exact_agreement=int((y1 == y2).sum()))


def deviance(N, P):
    E = P * N.sum()
    m = N > 0
    return float(2 * (N[m] * np.log(N[m] / E[m])).sum())


# ---------------------------------------------------------------- covariate R2 (spec 3.3, v1.2 items 5-7)
def design_matrix(Z, S):
    """Centred design for one cohort sample (rows may repeat): an all-missing column is dropped, a partly missing
    column is mean-filled and gets an indicator appended after the base columns in base-column order."""
    cols, miss = [], []
    for j in range(Z.shape[1]):
        v = Z[:, j]
        m = np.isnan(v)
        if m.all():
            continue
        v = np.where(m, np.nanmean(v), v)
        cols.append(v - v.mean())
        if m.any():
            miss.append(m.astype(float))
    cols += [m - m.mean() for m in miss]
    X = np.column_stack(cols) if cols else np.zeros((len(S), 0))
    return X, S - S.mean()


def r2_of(Z, S):
    """(n, rank p, raw R2, adjusted R2 floored at 0) of S on the within-cohort design."""
    X, s = design_matrix(Z, np.asarray(S, float))
    n = len(s)
    p = int(np.linalg.matrix_rank(X)) if X.shape[1] else 0
    sts = float(s @ s)
    if sts == 0 or p == 0:
        return n, p, np.nan, np.nan
    beta, *_ = np.linalg.lstsq(X, s, rcond=None)
    res = s - X @ beta
    raw = float(1 - res @ res / sts)
    adj = max(0.0, 1 - (1 - raw) * (n - 1) / (n - p - 1)) if n - p - 1 > 0 else np.nan
    return n, p, raw, adj


# ---------------------------------------------------------------- latent-shape sensitivity
def latent_kernel(cuts, var_e, grid=LATENT_GRID):
    sd_e = np.sqrt(max(var_e, 1e-10))
    cdf = norm.cdf((cuts[:, None] - grid[None, :]) / sd_e)
    return np.diff(np.vstack([np.zeros(len(grid)), cdf, np.ones(len(grid))]), axis=0)


def latent_shape_sensitivity(counts, cuts, var_e, n_iter=20000):
    """lambda for S at one cohort with a free latent distribution and the arm's var_e.

    NPMLE: the grid mixture that maximizes the likelihood of the observed single-read counts (EM).
    LP bounds: min and max of lambda over all grid mixtures that reproduce the observed proportions exactly;
    Var(S) and E[S] are then fixed at the observed values, so Var(E[S|T]) is linear in the weights.
    """
    counts = np.asarray(counts, float)
    p = counts / counts.sum()
    K = latent_kernel(cuts, var_e)
    g = latent_kernel(CUTS_S, var_e).T @ S_VALUES        # E[S | T = t]
    s_of_level = harmonize(np.arange(5)) if len(cuts) == 4 else S_VALUES
    m_obs = p @ s_of_level
    v_obs = p @ s_of_level ** 2 - m_obs ** 2
    w = np.full(K.shape[1], 1.0 / K.shape[1])
    for _ in range(n_iter):
        f = K @ w
        w = w * (K.T @ (p / np.clip(f, 1e-300, None)))
    f = K @ w
    m_np = w @ g
    ps = latent_kernel(CUTS_S, var_e) @ w
    v_np = ps @ S_VALUES ** 2 - (ps @ S_VALUES) ** 2
    cov_np = w @ g ** 2 - m_np ** 2
    out = dict(npmle_lambda=float(cov_np / v_np), npmle_max_abs_marginal_error=float(np.abs(f - p).max()),
               npmle_cov=float(cov_np), npmle_var=float(v_np), obs_ES=float(m_obs), obs_VarS=float(v_obs))
    # LP bounds at the observed marginal; if no mixture reproduces it at this var_e, at the NPMLE-fitted marginal
    for marginal, b in (("observed", p), ("npmle_fitted", f)):
        m_b = b @ s_of_level
        v_b = b @ s_of_level ** 2 - m_b ** 2
        res = {s: linprog(sg * g ** 2, A_eq=K, b_eq=b, bounds=(0, None), method="highs")
               for s, sg in (("min", 1.0), ("max", -1.0))}
        if all(r.status == 0 for r in res.values()):
            break
    out.update(lp_marginal=marginal if all(r.status == 0 for r in res.values()) else "infeasible",
               lp_var_S=float(v_b))
    for sense, r in res.items():
        ok = r.status == 0
        cov = float(g ** 2 @ r.x - m_b ** 2) if ok else np.nan
        out[f"lp_{sense}_cov"] = cov
        out[f"lp_{sense}_lambda"] = float(cov / v_b) if ok else np.nan
    return out


# ---------------------------------------------------------------- inputs
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path, **kw):
    return pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""], **kw)


def gse193066_pairs(placement_path, source_dir):
    p = read_tsv(placement_path, dtype={"source_fibrosis_stage": "Int64"})
    p["participant"] = p["donor_id"].str.split("::").str[1]
    first = p[p["biopsy"] == "1st biopsy"].set_index("participant")
    second = p[p["biopsy"] == "2nd biopsy"].set_index("participant")
    if first.index.duplicated().any() or second.index.duplicated().any():
        raise SystemExit("GSE193066: a donor has two rows for the same biopsy")
    both = sorted(set(first.index) & set(second.index))
    pairs = pd.DataFrame({"participant": both,
                          "y1": first.loc[both, "source_fibrosis_stage"].astype(int).values,
                          "y2": second.loc[both, "source_fibrosis_stage"].astype(int).values,
                          "run1": first.loc[both, "sample_id"].values, "run2": second.loc[both, "sample_id"].values})
    rec = read_tsv(Path(source_dir) / "results/pairs.tsv")
    rec = rec.rename(columns={"participant_id": "participant", "f1": "y1", "f2": "y2"})
    rec = rec[["participant", "y1", "y2"]].sort_values("participant").reset_index(drop=True)
    same = pairs[["participant", "y1", "y2"]].reset_index(drop=True).equals(rec.astype({"y1": int, "y2": int}))
    j = json.load(open(Path(source_dir) / "results/kleiner_reliability.json"))
    b = j["bootstrap_all_58"]["sigma_within"]
    recorded = dict(point=b["point"], ci_lo=b["ci_lo"], ci_hi=b["ci_hi"],
                    source=str(Path(source_dir) / "results/kleiner_reliability.json") + " bootstrap_all_58.sigma_within")
    return pairs, same, recorded


def gse193066_series_pairs(series_path):
    """(participant, y1, y2) re-derived from the GEO series matrix: title HUnafldNNN[_1|_2], biopsy, fibrosis stage."""
    import gzip
    rows = {}
    with gzip.open(series_path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title") or line.startswith("!Sample_characteristics_ch1"):
                f = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
                if f[0] == "!Sample_title":
                    rows["title"] = f[1:]
                else:
                    key = f[1].split(":")[0].strip()
                    if all(":" in x for x in f[1:]) and all(x.split(":")[0].strip() == key for x in f[1:]):
                        rows[key] = [x.split(":", 1)[1].strip() for x in f[1:]]
    df = pd.DataFrame({"title": rows["title"], "biopsy": rows["biopsy"], "fib": rows["fibrosis stage"]})
    df["participant"] = df["title"].str.split("_").str[0]
    wide = df.pivot(index="participant", columns="biopsy", values="fib").dropna()
    out = pd.DataFrame({"participant": wide.index, "y1": wide["1st biopsy"].astype(int).values,
                        "y2": wide["2nd biopsy"].astype(int).values}).sort_values("participant")
    return out.reset_index(drop=True), int(len(df)), int(df["participant"].nunique())


def pxd_pairs(source_dir, meta_path):
    rec = read_tsv(Path(source_dir) / "results/pairs.tsv")
    j = json.load(open(Path(source_dir) / "results/kleiner_replication.json"))
    df = read_tsv(meta_path, dtype=str)
    bad = sorted(set(df["kleiner_fibrosis_grade"]) - set(KLEINER_LEVELS))
    if bad:
        raise SystemExit(f"PXD051911: unmapped kleiner_fibrosis_grade levels {bad}")
    df["visit"] = df["unique_identifier"].str.split("_").str[0]
    df["fib"] = df["kleiner_fibrosis_grade"].map(KLEINER_LEVELS)
    if df.duplicated(["patient_name", "visit"]).any():
        raise SystemExit("PXD051911: a patient has two rows at one visit")
    wide = df.pivot(index="patient_name", columns="visit", values="fib")
    has = {v: wide[v].notna() if v in wide else pd.Series(False, index=wide.index) for v in ("V1", "V2", "V3")}
    rederived = {"PXD051911_A1": wide[has["V1"] & has["V2"]][["V1", "V2"]],
                 "PXD051911_A2": wide[has["V1"] & has["V3"] & ~has["V2"]][["V1", "V3"]]}
    out = {}
    for arm, key in PXD_ARMS.items():
        r = rec[rec["arm"] == key].rename(columns={"patient": "participant", "f1": "y1", "f2": "y2"})
        r = r[["participant", "y1", "y2"]].sort_values("participant").reset_index(drop=True).astype({"y1": int, "y2": int})
        d = rederived[arm]
        d = pd.DataFrame({"participant": d.index, "y1": d.iloc[:, 0].astype(int).values,
                          "y2": d.iloc[:, 1].astype(int).values}).sort_values("participant").reset_index(drop=True)
        b = j["bootstrap"][key]["sigma_e_stages"]
        recorded = dict(point=b["point"], ci_lo=b["ci_lo"], ci_hi=b["ci_hi"],
                        source=str(Path(source_dir) / "results/kleiner_replication.json") + f" bootstrap.{key}.sigma_e_stages")
        out[arm] = (r, bool(r.equals(d)), recorded)
    return out


def read_frozen_crosswalk(path):
    rec = Path(path).parent / "frozen_crosswalk.sha256"
    expected = rec.read_text().split()[0]
    got = sha256(path)
    if got != expected:
        raise SystemExit(f"frozen crosswalk sha256 {got} != {expected} recorded in {rec}")
    return read_tsv(path, dtype=str), dict(path=str(path), sha256=got, record=str(rec), matches=True)


def target_participants(x, meta):
    """Staged unit individuals of F and F_sealed (spec v1 1.1) from the frozen seal crosswalk.

    unit_library True (one library per development or sealed individual), C = 0, S observed. Kleiner cohorts: Y is
    the 0-4 metadata stage and must reproduce the frozen S. GSE213621: S must equal the coded condition label.
    """
    u = x[x["unit_library"] == "True"].copy()
    if u["individual_id"].isna().any() or u["individual_id"].duplicated().any():
        raise SystemExit("frozen crosswalk: unit_library is not exactly one library per individual")
    if not u["run_role"].isin(["development", "sealed"]).all():
        raise SystemExit("frozen crosswalk: a unit library is neither development nor sealed")
    staged = (u["C"].astype(float) == 0) & u["S"].notna()
    keep = pd.Series(False, index=u.index)
    for s, (role, cohorts) in SETS.items():
        keep |= (u["run_role"] == role) & u["cohort"].isin(cohorts)
    t = u[staged & keep].copy()
    t["set"] = t["run_role"].map({role: s for s, (role, _) in SETS.items()})
    t["S"] = t["S"].astype(float).astype(int)
    stage = meta.set_index("sample_id")["fibrosis_stage"]
    cond = meta.set_index("sample_id")["condition"]
    g = t["cohort"] == "GSE213621"
    t["Y"] = np.nan
    ystr = t.loc[~g, "run"].map(stage)
    if not ystr.isin(["0", "1", "2", "3", "4"]).all():
        raise SystemExit(f"target: staged Kleiner unit libraries without an integer 0-4 stage: "
                         f"{t.loc[~g][~ystr.isin(['0', '1', '2', '3', '4'])]['run'].tolist()[:5]}")
    t.loc[~g, "Y"] = ystr.astype(int)
    ok_k = (harmonize(t.loc[~g, "Y"]) == t.loc[~g, "S"]).all()
    ok_g = (t.loc[g, "run"].map(cond).map(GSE213621_CONDITION) == t.loc[g, "S"]).all()
    if not (ok_k and ok_g):
        raise SystemExit("target: frozen S disagrees with the metadata stage or condition label")
    counts = {s: {c: int(v) for c, v in t[t["set"] == s]["cohort"].value_counts().sort_index().items()} for s in SETS}
    t = t.sort_values(["set", "cohort", "individual_id"]).reset_index(drop=True)
    return t[["individual_id", "set", "cohort", "run", "Y", "S"]], counts


def identity_classes(pairs, x):
    """Provisional identity class of each GSE193066 pair from the frozen crosswalk individual_id."""
    ind = x.set_index("run")["individual_id"]
    i1, i2 = pairs["run1"].map(ind), pairs["run2"].map(ind)
    return np.where(i1.isna() | i2.isna(), "uncheckable", np.where(i1 == i2, "same", "different"))


def gse213621_pairs(x, meta, kinship_path):
    """Two-library individuals of GSE213621 in the frozen crosswalk, S of each library from the condition label."""
    g = x[x["cohort"] == "GSE213621"]
    multi = g[g["individual_id"].notna()].groupby("individual_id")["run"].apply(sorted)
    multi = multi[multi.str.len() >= 2]
    if (multi.str.len() > 2).any():
        raise SystemExit("GSE213621: an individual has more than two libraries")
    cond = meta.set_index("sample_id")["condition"].map(GSE213621_CONDITION)
    pairs = pd.DataFrame({"participant": multi.index, "run1": multi.str[0].values, "run2": multi.str[1].values})
    pairs["y1"], pairs["y2"] = pairs["run1"].map(cond), pairs["run2"].map(cond)
    if pairs[["y1", "y2"]].isna().any().any():
        raise SystemExit("GSE213621: a same-person library has no staged condition label")
    pairs[["y1", "y2"]] = pairs[["y1", "y2"]].astype(int)
    k = read_tsv(kinship_path, dtype=str)
    k = k[k["relation"] == "same_individual"]
    kp = {frozenset(r) for r in zip(k["#IID1"], k["IID2"])}
    ours = {frozenset(r) for r in zip(pairs["run1"], pairs["run2"])}
    kin_in = {p for p in kp if all(r in set(g["run"]) for r in p)}
    return pairs.sort_values("participant").reset_index(drop=True), ours == kin_in


def f_covariates(tp, design_path, ancestry_path):
    """F individuals with the spec z columns: design_individuals.tsv plus the ancestry axes joined by the unit run."""
    d = pd.read_csv(design_path, sep="\t", dtype={"individual_id": str, "cluster_id": str})
    d = d[d["set"] == "development"]
    anc = pd.read_csv(ancestry_path, sep="\t", usecols=["run", "axis_EUR_AFR", "axis_EUR_EAS"], dtype={"run": str})
    f = tp[tp["set"] == "F"]
    m = f.merge(d.drop(columns=["set"]), on="individual_id", how="left", suffixes=("", "_design"))
    m = m.merge(anc, on="run", how="left")
    same = (len(d) == len(f) and set(d["individual_id"]) == set(f["individual_id"])
            and (m["cohort"] == m["cohort_design"]).all() and (m["S"] == m["S_design"]).all())
    return m, d, bool(same)


def p2_resample_index(path, fz, design, n_boot):
    """Row positions into each F cohort frame for P2 resample b = 0..n_boot-1 (clusters expanded within cohort)."""
    d = pd.read_csv(path, sep="\t", dtype={"cluster_id": str, "stratum": str})
    n_file = int(d["draw"].max()) + 1
    d = d[d["draw"] < n_boot].sort_values(["draw", "stratum", "slot"])
    pos = {c: dict(zip(g["individual_id"], range(len(g)))) for c, g in fz.items()}
    members = design.groupby("cluster_id")["individual_id"].apply(list).to_dict()
    home = dict(zip(design["individual_id"], design["cohort"]))
    idx = [dict() for _ in range(n_boot)]
    for (b, c), g in d.groupby(["draw", "stratum"], sort=False):
        idx[b][c] = np.array([pos[c][i] for k in g["cluster_id"] for i in members[k] if home[i] == c], int)
    complete = all(set(r) == set(F_COHORTS) for r in idx)
    sizes = all(len(idx[b][c]) == len(fz[c]) for b in range(n_boot) for c in F_COHORTS) if complete else False
    return idx, dict(draws_in_file=n_file, complete=complete, resample_size_equals_cohort_size=sizes,
                     max_cluster_size=int(max(len(v) for v in members.values())))


# ---------------------------------------------------------------- workers
def interval(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] if v.size else [None, None]


def fit_spec(arm, spec, tables, warm=None, var_e_starts=()):
    sp = SPECS[spec]
    if arm["kind"] == "joint":
        return fit_joint(tables, warm=warm, rho_nonneg=sp["rho_nonneg"], var_e_starts=var_e_starts)
    return fit_pairs(tables, warm=warm, fix_delta=sp["fix_delta"], rho_nonneg=sp["rho_nonneg"])


def summarise_spec(arm, spec, theta, targets, r2):
    if arm["kind"] == "joint":
        return summarise_joint(theta, targets, SPECS[spec]["scale"], r2)
    return summarise_fit(theta, targets, SPECS[spec]["scale"], r2)


def arm_tables(arm, spec, idx=None):
    scale = SPECS[spec]["scale"]
    if arm["kind"] == "joint":
        return [table_for(m["y1"][ii], m["y2"][ii], m["input_scale"], scale)
                for m, ii in zip(arm["members"], idx if idx is not None else [slice(None)] * 2)]
    ii = idx if idx is not None else slice(None)
    return table_for(arm["y1"][ii], arm["y2"][ii], arm["input_scale"], scale)


def arm_empirical(arm, idx=None):
    if arm["kind"] == "joint":
        out = {}
        for m, ii in zip(arm["members"], idx if idx is not None else [slice(None)] * 2):
            out.update({f"{k}_{m['name']}": v for k, v in empirical(m["y1"][ii], m["y2"][ii], m["input_scale"]).items()})
        return out
    ii = idx if idx is not None else slice(None)
    return empirical(arm["y1"][ii], arm["y2"][ii], arm["input_scale"])


def run_bootstrap(task):
    """Draws task['draws'] of one arm and spec; target_draws[b] and r2_draws[b] are P2 resample b."""
    arm, spec, idx, warm = task["arm"], task["spec"], task["idx"], task["warm"]
    rows = []
    for b in task["draws"]:
        ib = [x[b] for x in idx] if arm["kind"] == "joint" else idx[b]
        f = fit_spec(arm, spec, arm_tables(arm, spec, ib), warm=warm)
        rows.append(dict(arm=arm["name"], spec=spec, draw=b, converged=f["converged"], at_bound=f["at_bound"],
                         **summarise_spec(arm, spec, f["theta"], task["target_draws"][b], task["r2_draws"][b]),
                         **arm_empirical(arm, ib)))
    return rows


def run_het_bootstrap(task):
    arm, idx, warm = task["arm"], task["idx"], task["warm"]
    rows = []
    for b in task["draws"]:
        N = arm_tables(arm, "Y5", idx[b])
        f = fit_het(N, [warm, het_from_hom(moment_start(N))])
        rows.append(dict(arm=arm["name"], spec="Y5_het", draw=b, converged=f["converged"], at_bound=f["at_bound"],
                         **summarise_het(f["theta"], task["target_draws"][b], task["r2_draws"][b])))
    return rows


def run_gof(task):
    arm, spec, sims, warm = task["arm"], task["spec"], task["sims"], task["warm"]
    devs, fails = [], 0
    n_sim = len(sims[0]) if arm["kind"] == "joint" else len(sims)
    for b in range(n_sim):
        tabs = [s[b] for s in sims] if arm["kind"] == "joint" else sims[b]
        if spec == "Y5_het":
            f = fit_het(tabs, [warm])
            mu, lsd, a, bb, d = f["theta"]
            devs.append(deviance(tabs, het_cell_probs(CUTS_Y, mu, np.exp(lsd), a, bb, d)))
        else:
            f = fit_spec(arm, spec, tabs, warm=warm)
            devs.append(point_deviance(arm, spec, f["theta"], tabs))
        fails += int(not f["converged"])
    return dict(arm=arm["name"], spec=spec, devs=np.array(devs), n_not_converged=fails)


def point_deviance(arm, spec, theta, tables):
    scale = SPECS[spec]["scale"]
    if arm["kind"] == "joint":
        ps, _ = joint_arm_params(theta)
        return sum(deviance(N, cell_probs(CUTS[scale], mu, mu + d, sd, rho)) for N, (mu, sd, rho, d) in zip(tables, ps))
    mu, sd, rho, d = unpack(theta)
    return deviance(tables, cell_probs(CUTS[scale], mu, mu + d, sd, rho))


def simulate_tables(arm, spec, theta, tables, rng, B):
    if spec == "Y5_het":
        mu, lsd, a, b, d = theta
        P = het_cell_probs(CUTS_Y, mu, np.exp(lsd), a, b, d).ravel()
        return rng.multinomial(int(tables.sum()), P / P.sum(), size=B).reshape(B, 5, 5).astype(float)
    scale = SPECS[spec]["scale"]
    k = len(CUTS[scale]) + 1
    if arm["kind"] == "joint":
        ps, _ = joint_arm_params(theta)
        out = []
        for N, (mu, sd, rho, d) in zip(tables, ps):
            P = cell_probs(CUTS[scale], mu, mu + d, sd, rho).ravel()
            out.append(rng.multinomial(int(N.sum()), P / P.sum(), size=B).reshape(B, k, k).astype(float))
        return out
    mu, sd, rho, d = unpack(theta)
    P = cell_probs(CUTS[scale], mu, mu + d, sd, rho).ravel()
    return rng.multinomial(int(tables.sum()), P / P.sum(), size=B).reshape(B, k, k).astype(float)


def bvn_check(rng, n=60):
    """Compare bvn_cdf with scipy's multivariate normal CDF at random points, including sign changes."""
    h, k = rng.uniform(-4, 4, n), rng.uniform(-4, 4, n)
    rho = rng.uniform(-0.95, 0.95, n)
    ours = np.array([bvn_cdf(a, b, r) for a, b, r in zip(h, k, rho)], float)
    ref = np.array([multivariate_normal.cdf([a, b], mean=[0, 0], cov=[[1, r], [r, 1]], abseps=1e-12, releps=1e-12)
                    for a, b, r in zip(h, k, rho)])
    return float(np.max(np.abs(ours - ref)))


def rnd(d):
    return {k: (round(float(v), 6) if isinstance(v, (float, np.floating)) else v) for k, v in d.items()}


def chunks(n):
    return [list(range(i, min(i + CHUNK, n))) for i in range(0, n, CHUNK)]


# ---------------------------------------------------------------- run
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for flag in ("--gse193066-placement", "--gse193066-source", "--gse193066-series", "--gse193066-crosswalk",
                 "--pxd-source", "--pxd-meta", "--frozen-crosswalk", "--kinship-pairs", "--metadata",
                 "--design-individuals", "--ancestry", "--r2-moments", "--p2-resamples", "--r2-by-resample",
                 "--out", "--restricted-out"):
        ap.add_argument(flag, required=True)
    ap.add_argument("--codex-lambda", default=None)
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    ap.add_argument("--n-gof", type=int, default=N_GOF)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    a = ap.parse_args()

    out, rout = Path(a.out), Path(a.restricted_out)
    names = ["label_reliability.json", "lambda_source_draws.tsv.gz", "carried_spec.tsv", "lambda_by_cohort.tsv",
             "het_sensitivity.tsv", "latent_shape_sensitivity.tsv", "goodness_of_fit.tsv", "bootstrap_draws.npz",
             "bootstrap_draws.tsv.gz", "pair_tables.tsv", "target_stage_counts.tsv"]
    rnames = ["pairs_used.tsv", "target_participants.tsv"]
    clash = [str(out / f) for f in names if (out / f).exists()] + [str(rout / f) for f in rnames if (rout / f).exists()]
    if clash:
        raise SystemExit(f"refusing to overwrite: {clash}")
    out.mkdir(parents=True, exist_ok=True)
    rout.mkdir(parents=True, exist_ok=True)

    guards, checks = [], []

    def guard(name, passed, detail):
        guards.append(dict(guard=name, passed=bool(passed), detail=detail))

    def check(name, passed, detail):
        checks.append(dict(check=name, passed=bool(passed), detail=detail))

    guard("bvn_cdf_matches_scipy", (err := bvn_check(np.random.default_rng(SEED))) < 1e-7,
          f"max |ours - scipy| over 60 random points = {err:.2e} (tolerance 1e-7)")

    # ---- pair inputs
    x, xw = read_frozen_crosswalk(a.frozen_crosswalk)
    guard("frozen_crosswalk_sha256_matches_record", xw["matches"], f"{xw['sha256']} = first field of {xw['record']}")
    meta = pd.read_csv(a.metadata, keep_default_na=False, na_values=[""], dtype=str)
    if meta["sample_id"].duplicated().any():
        raise SystemExit("metadata: duplicated sample_id")

    g_pairs, g_same, g_rec = gse193066_pairs(a.gse193066_placement, a.gse193066_source)
    guard("GSE193066_pairs_equal_recorded_pairs", g_same,
          f"{len(g_pairs)} pairs from the placement table vs the 58 rows of {a.gse193066_source}/results/pairs.tsv")
    sm_pairs, sm_n, sm_part = gse193066_series_pairs(a.gse193066_series)
    guard("GSE193066_pairs_equal_series_matrix",
          sm_pairs.equals(g_pairs[["participant", "y1", "y2"]].reset_index(drop=True)),
          f"{a.gse193066_series} (sha256 {sha256(a.gse193066_series)[:12]}...): {sm_n} samples, {sm_part} participants, "
          f"{len(sm_pairs)} pairs")
    cw = read_tsv(a.gse193066_crosswalk, dtype=str).set_index("run_id")["participant_token"]
    guard("GSE193066_placement_runs_match_source_crosswalk",
          (g_pairs["run1"].map(cw).eq(g_pairs["participant"]) & g_pairs["run2"].map(cw).eq(g_pairs["participant"])).all(),
          f"run1 and run2 of every pair map to its participant_token in {a.gse193066_crosswalk}")
    g_pairs["identity_class"] = identity_classes(g_pairs, x)
    ic = g_pairs["identity_class"].value_counts().to_dict()
    guard("GSE193066_identity_classes_reproduce_STATUS",
          (ic.get("same", 0), ic.get("different", 0), ic.get("uncheckable", 0)) == (29, 11, 18),
          f"frozen-crosswalk individual_id of run1 and run2: same {ic.get('same', 0)}, different "
          f"{ic.get('different', 0)}, uncheckable {ic.get('uncheckable', 0)} (docs/STATUS.md GSE193066 row: 29, 11, 18)")

    arms = {"GSE193066": dict(name="GSE193066", kind="single", input_scale="Y5",
                              y1=g_pairs["y1"].to_numpy(), y2=g_pairs["y2"].to_numpy())}
    recorded = {"GSE193066": g_rec}
    pairs_out = [g_pairs.assign(arm="GSE193066")]
    for arm, (pairs, same, rec) in pxd_pairs(a.pxd_source, a.pxd_meta).items():
        arms[arm] = dict(name=arm, kind="single", input_scale="Y5", y1=pairs["y1"].to_numpy(), y2=pairs["y2"].to_numpy())
        recorded[arm] = rec
        pairs_out.append(pairs.assign(arm=arm))
        guard(f"{arm}_pairs_equal_rederived_from_meta", same,
              f"{len(pairs)} recorded pairs vs the pairing re-derived from {a.pxd_meta}")
    arms["PXD051911_joint"] = dict(name="PXD051911_joint", kind="joint", input_scale="Y5",
                                   members=[dict(arms["PXD051911_A1"], name="A1"), dict(arms["PXD051911_A2"], name="A2")])
    for name, keep in (("GSE193066_same29", ["same"]), ("GSE193066_drop11", ["same", "uncheckable"])):
        sub = g_pairs[g_pairs["identity_class"].isin(keep)]
        arms[name] = dict(name=name, kind="single", input_scale="Y5", y1=sub["y1"].to_numpy(), y2=sub["y2"].to_numpy())
        pairs_out.append(sub.assign(arm=name))
    p25, p25_ok = gse213621_pairs(x, meta, a.kinship_pairs)
    disc = p25[p25["y1"] != p25["y2"]]
    kinds = pd.Series([tuple(sorted(r)) for r in zip(disc["y1"], disc["y2"])]).value_counts().to_dict()
    guard("GSE213621_same_person_pairs_reproduce_STATUS",
          p25_ok and len(p25) == 25 and len(disc) == 11 and kinds == {(1, 2): 6, (0, 1): 4, (0, 2): 1},
          f"{len(p25)} two-library individuals in the frozen crosswalk (equal to the GSE213621 same_individual pairs of "
          f"{a.kinship_pairs}: {p25_ok}); {len(disc)} carry two stage labels: "
          f"{ {f'S{k[0]}-S{k[1]}': v for k, v in kinds.items()} } (docs/STATUS.md GSE213621 row: 25; 11 = 6 F2/F3F4, "
          f"4 F0F1/F2, 1 F0F1/F3F4)")
    arms["GSE213621_pairs25"] = dict(name="GSE213621_pairs25", kind="single", input_scale="S3",
                                     y1=p25["y1"].to_numpy(), y2=p25["y2"].to_numpy())
    pairs_out.append(p25.assign(arm="GSE213621_pairs25", identity_class="same"))
    arms = {k: arms[k] for k in ARM_ORDER}
    arm_specs = {k: (PAIRS25_SPECS if k == "GSE213621_pairs25" else KLEINER_SPECS) for k in arms}

    # ---- targets: F and F_sealed
    tp, tcounts = target_participants(x, meta)
    guard("F_counts_equal_spec_1_1", tcounts["F"] == F_COUNTS,
          f"F staged development unit individuals by cohort {tcounts['F']} (total {sum(tcounts['F'].values())}); "
          f"spec v1 1.1: {F_COUNTS} (558)")
    guard("F_sealed_counts", tcounts["F_sealed"] == SEALED_COUNTS,
          f"F_sealed by cohort {tcounts['F_sealed']} (total {sum(tcounts['F_sealed'].values())}); "
          f"expected {SEALED_COUNTS} (237)")
    fcov, design, same_design = f_covariates(tp, a.design_individuals, a.ancestry)
    guard("F_equals_gateA_design_individuals", same_design,
          f"558 F individual_id, cohort and S equal the development rows of {a.design_individuals}")
    guard("F_axes_present", fcov[["axis_EUR_AFR", "axis_EUR_EAS"]].notna().all().all(),
          f"ancestry axes found for every F unit run in {a.ancestry}")
    fz = {c: fcov[fcov["cohort"] == c].sort_values("individual_id").reset_index(drop=True) for c in F_COHORTS}

    # observed R2, reproduced against the lead's moments (job 1290)
    r2_obs, r2_raw_obs, r2_rows = {}, {}, []
    rec_m = json.load(open(a.r2_moments))
    for c in F_COHORTS:
        n, p, raw, adj = r2_of(fz[c][Z_COLS].to_numpy(float), fz[c]["S"].to_numpy(float))
        r2_obs[f"F_{c}"], r2_raw_obs[f"F_{c}"] = adj, raw
        r2_rows.append(dict(cohort=c, n=n, rank_p=p, r2_raw=raw, r2_adjusted=adj,
                            recorded_raw=rec_m[c]["r2_raw"], recorded_adjusted=rec_m[c]["r2_adjusted"]))
    r2_diff = max(max(abs(r["r2_raw"] - r["recorded_raw"]), abs(r["r2_adjusted"] - r["recorded_adjusted"]))
                  for r in r2_rows)
    guard("R2_reproduces_r2_moments_observed", r2_diff < 1e-10,
          f"max |ours - {a.r2_moments}| over raw and adjusted R2 = {r2_diff:.2e}; adjusted: "
          + ", ".join(f"{r['cohort']} {r['r2_adjusted']:.4f}" for r in r2_rows))
    log("inputs read")

    # ---- target single-read fits and fit tests
    rng_streams = np.random.SeedSequence(SEED).spawn(len(ARM_ORDER) + 4)
    stream = dict(zip(ARM_ORDER + ["targets_sealed", "gof", "target_gof", "het_gof"], rng_streams))
    rng_tgof = np.random.default_rng(stream["target_gof"])
    targets, tc_rows, counts_of = [], [], {}
    for s, (_, cohorts) in SETS.items():
        for c in cohorts:
            tc = tp[(tp["set"] == s) & (tp["cohort"] == c)]
            scale = "S3" if c == "GSE213621" else "Y5"
            cnt = (np.bincount(tc["S"], minlength=3) if scale == "S3" else np.bincount(tc["Y"].astype(int), minlength=5))
            key = f"{s}_{c}"
            counts_of[key] = cnt
            f = fit_single(cnt, CUTS[scale])
            targets.append(dict(key=key, set=s, cohort=c, scale=scale, n=int(cnt.sum()), counts=cnt,
                                mu_Z=f["mu_Z"], sd_Z=f["sd_Z"]))
            exp_c = f["fitted_counts"]
            x2 = float(((cnt - exp_c) ** 2 / exp_c).sum())
            df = len(cnt) - 1 - 2
            row = dict(set=s, cohort=c, label_scale=scale, n=int(cnt.sum()),
                       **{f"Y{k}": (int(cnt[k]) if scale == "Y5" else "") for k in range(5)},
                       **{f"S{v}": int((tc["S"] == v).sum()) for v in range(3)},
                       **{f"fitted_{'Y' if scale == 'Y5' else 'S'}{k}": round(float(exp_c[k]), 3) for k in range(len(cnt))},
                       mu_Z=f["mu_Z"], sd_Z=f["sd_Z"], converged=f["converged"], at_bound=f["at_bound"],
                       gof_df=df, gof_pearson_X2=round(x2, 4))
            if df > 0:
                sims = rng_tgof.multinomial(int(cnt.sum()), f["fitted_p"] / f["fitted_p"].sum(), size=a.n_gof)
                xb = []
                for sm in sims:
                    fs = fit_single(sm, CUTS[scale])
                    xb.append(((sm - fs["fitted_counts"]) ** 2 / fs["fitted_counts"]).sum())
                row.update(gof_p_asymptotic=float(chi2.sf(x2, df)),
                           gof_p_parametric_bootstrap=float((1 + np.sum(np.array(xb) >= x2 - 1e-9)) / (a.n_gof + 1)),
                           gof_note="")
            else:
                mu_x, sd_x = exact_three_level_probit(cnt)
                row.update(gof_p_asymptotic=np.nan, gof_p_parametric_bootstrap=np.nan,
                           gof_note=f"3 levels, 2 parameters: exactly identified, no fit test (closed form mu_Z "
                                    f"{mu_x:.4f}, sd_Z {sd_x:.4f})")
                guard(f"{key}_single_read_fit_equals_closed_form",
                      abs(mu_x - f["mu_Z"]) < 1e-3 and abs(sd_x - f["sd_Z"]) < 1e-3,
                      f"numerical mu_Z {f['mu_Z']:.6f}, sd_Z {f['sd_Z']:.6f}; closed form {mu_x:.6f}, {sd_x:.6f}")
            tc_rows.append(row)
    guard("target_single_read_fits_converged", all(r["converged"] and not r["at_bound"] for r in tc_rows),
          "; ".join(f"{r['set']} {r['cohort']} n={r['n']} mu_Z={r['mu_Z']:.3f} sd_Z={r['sd_Z']:.3f}" for r in tc_rows))
    log("target fits done")

    # ---- target draws: F = P2 resample b (spec 5.2) with R2 recomputed on it; F_sealed resampled here
    p2_idx, p2_info = p2_resample_index(a.p2_resamples, fz, design, a.n_boot)
    guard("P2_resamples_cover_every_draw_and_cohort",
          p2_info["complete"] and p2_info["resample_size_equals_cohort_size"] and p2_info["draws_in_file"] >= a.n_boot,
          f"{a.p2_resamples}: {p2_info['draws_in_file']} draws in file, first {a.n_boot} used; every draw has the 5 F "
          f"strata: {p2_info['complete']}; resample size = cohort size: {p2_info['resample_size_equals_cohort_size']}; "
          f"largest development cluster {p2_info['max_cluster_size']}")
    label_of = {c: (fz[c]["S"] if c == "GSE213621" else fz[c]["Y"]).to_numpy(int) for c in F_COHORTS}
    zmat = {c: fz[c][Z_COLS].to_numpy(float) for c in F_COHORTS}
    smat = {c: fz[c]["S"].to_numpy(float) for c in F_COHORTS}
    base = {t["key"]: t for t in targets}
    rng_ts = np.random.default_rng(stream["targets_sealed"])
    sealed_sims = {t["key"]: rng_ts.multinomial(t["n"], t["counts"] / t["n"], size=a.n_boot)
                   for t in targets if t["set"] == "F_sealed"}
    target_draws, r2_draws, r2_draw_rows = [], [], []
    for b in range(a.n_boot):
        tb, rb = [], {}
        for c in F_COHORTS:
            ii = p2_idx[b][c]
            key = f"F_{c}"
            t = base[key]
            cnt = np.bincount(label_of[c][ii], minlength=len(t["counts"]))
            fs = fit_single(cnt, CUTS[t["scale"]])
            tb.append(dict(t, counts=cnt, mu_Z=fs["mu_Z"], sd_Z=fs["sd_Z"]))
            n, p, raw, adj = r2_of(zmat[c][ii], smat[c][ii])
            rb[key] = adj
            r2_draw_rows.append(dict(draw=b, cohort=c, n=n, rank_p=p, r2_raw=raw, r2_adjusted=adj))
        for key, sims in sealed_sims.items():
            fs = fit_single(sims[b], CUTS[base[key]["scale"]])
            tb.append(dict(base[key], counts=sims[b], mu_Z=fs["mu_Z"], sd_Z=fs["sd_Z"]))
        target_draws.append(tb)
        r2_draws.append(rb)
    r2d = pd.DataFrame(r2_draw_rows)
    rec_r2 = pd.read_csv(a.r2_by_resample, sep="\t")
    cmp = r2d.merge(rec_r2, on=["draw", "cohort"], suffixes=("", "_rec"), how="left")
    d_r2 = float(np.nanmax(np.abs(np.r_[cmp["r2_raw"] - cmp["r2_raw_rec"], cmp["r2_adjusted"] - cmp["r2_adjusted_rec"]])))
    guard("R2_by_resample_reproduces_record", cmp["r2_adjusted_rec"].notna().all() and d_r2 < 1e-10,
          f"per draw and cohort, max |ours - {a.r2_by_resample}| = {d_r2:.2e} over {len(cmp)} rows")
    r2_shift = {c: dict(observed_adjusted=r2_obs[f"F_{c}"],
                        resample_median_adjusted=float(r2d.loc[r2d["cohort"] == c, "r2_adjusted"].median()),
                        resample_p2p5_p97p5=interval(r2d.loc[r2d["cohort"] == c, "r2_adjusted"]))
                for c in F_COHORTS}
    log("target draws done")

    # ---- point fits
    point, fits, point_tables, table_rows = {}, {}, {}, []
    for name, arm in arms.items():
        emp = arm_empirical(arm)
        for spec in arm_specs[name]:
            tabs = arm_tables(arm, spec)
            ve_starts = ()
            if arm["kind"] == "joint":
                ve_starts = tuple(point[(m, spec)]["var_e"] for m in ("PXD051911_A1", "PXD051911_A2"))
            f = fit_spec(arm, spec, tabs, var_e_starts=ve_starts)
            fits[(name, spec)] = f
            point_tables[(name, spec)] = tabs
            point[(name, spec)] = dict(summarise_spec(arm, spec, f["theta"], targets, r2_obs), **emp,
                                       negll=f["negll"], converged=f["converged"], at_bound=f["at_bound"])
            if spec in ("Y5", "S3", "S3_sym") and arm["kind"] == "single":
                for i in range(tabs.shape[0]):
                    table_rows.append(dict(arm=name, scale=SPECS[spec]["scale"], first=i,
                                           **{f"repeat_{j}": int(tabs[i, j]) for j in range(tabs.shape[0])}))

    # SD agreement: a guard for the recorded arms' Y5 fits; a check with a consequence for every other fit
    sd_ok = {}
    for (name, spec), r in point.items():
        scale = SPECS[spec]["scale"]
        for m in (["A1", "A2"] if arms[name]["kind"] == "joint" else [None]):
            sfx = f"_{m}" if m else ""
            imp = r[f"implied_within_pair_sd{sfx}"]
            emp_sd = r[f"empirical_within_pair_sd{sfx}" if scale == "Y5" else f"empirical_within_pair_sd_S{sfx}"]
            label = f"{name}{sfx}_{spec}_implied_vs_empirical_within_pair_sd"
            passed = abs(imp - emp_sd) <= SD_GUARD_TOL
            sd_ok[(name, spec)] = sd_ok.get((name, spec), True) and passed
            detail = f"{'Y' if scale == 'Y5' else 'S'} scale: model-implied {imp:.4f}, empirical {emp_sd:.4f}, |diff| " \
                     f"{abs(imp - emp_sd):.4f} (tolerance {SD_GUARD_TOL})"
            (guard if name in PRIMARY_Y5_ARMS and spec == "Y5" else check)(label, passed, detail)
    for name in ("GSE193066", "PXD051911_A1", "PXD051911_A2"):
        e = point[(name, "Y5")]["empirical_within_pair_sd"]
        guard(f"{name}_loader_reproduces_recorded_within_pair_sd", abs(e - recorded[name]["point"]) < 1e-12,
              f"sqrt(sum d^2 / 2n) = {e:.16f}; recorded {recorded[name]['point']:.16f} ({recorded[name]['source']})")

    # carried spec per source: the larger var_e = the smaller lambda in every cohort (lambda_c falls with var_e)
    carried, carried_rows = {}, []
    for src in SOURCES:
        ve = {s: point[(src, s)]["var_e"] for s in SOURCE_SPECS}
        cs = max(SOURCE_SPECS, key=lambda s: ve[s])
        other = [s for s in SOURCE_SPECS if s != cs][0]
        smaller_everywhere = all(point[(src, cs)][f"lambda_{t['key']}"] <= point[(src, other)][f"lambda_{t['key']}"]
                                 + 1e-12 for t in targets)
        guard(f"{src}_carried_spec_smaller_in_every_target", smaller_everywhere,
              f"{cs} var_e {ve[cs]:.6f} >= {other} var_e {ve[other]:.6f}; lambda under {cs} <= under {other} in every "
              f"target cohort, so also in every pool and after the covariate correction")
        carried[src] = cs
        carried_rows.append(dict(source=src, carried_spec=cs, var_e_carried=ve[cs], other_spec=other,
                                 var_e_other=ve[other],
                                 lambdacond_hc_F_carried=point[(src, cs)]["lambdacond_hc_F"],
                                 lambdacond_hc_F_other=point[(src, other)]["lambdacond_hc_F"],
                                 rule="spec v1 5.2: the smaller of Y5 and S3 at the point estimate, fixed across draws"))

    # shared var_e in PXD051911: likelihood ratio against separate A1 and A2 fits (1 df, asymptotic)
    joint_lrt = {}
    for spec in ("Y5", "S3"):
        stat = 2.0 * (fits[("PXD051911_joint", spec)]["negll"] - fits[("PXD051911_A1", spec)]["negll"]
                      - fits[("PXD051911_A2", spec)]["negll"])
        joint_lrt[spec] = dict(lrt=float(stat), df=1, p_asymptotic=float(chi2.sf(max(stat, 0.0), 1)),
                               var_e_A1=point[("PXD051911_A1", spec)]["var_e"],
                               var_e_A2=point[("PXD051911_A2", spec)]["var_e"],
                               var_e_shared=point[("PXD051911_joint", spec)]["var_e"])

    # no-shift sensitivity at the point (Y5 single-table arms)
    no_shift = {}
    for name, arm in arms.items():
        if arm["kind"] == "single" and arm["input_scale"] == "Y5":
            f = fit_pairs(arm_tables(arm, "Y5"), fix_delta=True)
            r = summarise_fit(f["theta"], targets, r2=r2_obs)
            no_shift[name] = dict(var_e=r["var_e"], lambda_hc_F=r["lambda_hc_F"], lambdacond_hc_F=r["lambdacond_hc_F"],
                                  negll=f["negll"], converged=f["converged"])

    # heteroscedastic sensitivity (sources, Y5)
    het = {}
    for src in SOURCES:
        N = point_tables[(src, "Y5")]
        h0 = het_from_hom(fits[(src, "Y5_rho_nonneg")]["theta"])
        starts = [h0] + [np.r_[h0[:3], bb, h0[4]] for bb in (-0.5, 0.5, 1.0)]
        f = fit_het(N, starts)
        row = summarise_het(f["theta"], targets, r2_obs)
        lrt = 2.0 * (fits[(src, "Y5_rho_nonneg")]["negll"] - f["negll"])
        emp_sd = point[(src, "Y5")]["empirical_within_pair_sd"]
        mu, lsd, a_, b_, d_ = f["theta"]
        het[src] = dict(fit=f, point=dict(row, negll=f["negll"], converged=f["converged"], at_bound=f["at_bound"]),
                        deviance=deviance(N, het_cell_probs(CUTS_Y, mu, np.exp(lsd), a_, b_, d_)),
                        lrt_vs_homoscedastic=dict(lrt=float(lrt), df=1, p_asymptotic=float(chi2.sf(max(lrt, 0.0), 1))),
                        sd_check=dict(implied=row["implied_within_pair_sd"], empirical=emp_sd,
                                      passed=abs(row["implied_within_pair_sd"] - emp_sd) <= SD_GUARD_TOL))
    log("point fits done")

    # ---- bootstrap and goodness of fit (all random draws made here; workers are deterministic)
    idx = {}
    for name in ("GSE193066", "PXD051911_A1", "PXD051911_A2", "GSE193066_same29", "GSE193066_drop11",
                 "GSE213621_pairs25"):
        n = len(arms[name]["y1"])
        rng = np.random.default_rng(stream[name])
        idx[name] = np.stack([rng.integers(0, n, n) for _ in range(a.n_boot)])
    idx["PXD051911_joint"] = [idx["PXD051911_A1"], idx["PXD051911_A2"]]
    def draw_slice(dr):
        return dict(draws=dr, target_draws={b: target_draws[b] for b in dr}, r2_draws={b: r2_draws[b] for b in dr})

    boot_tasks = [dict(arm=arms[k], spec=s, idx=idx[k], warm=fits[(k, s)]["theta"], **draw_slice(dr))
                  for k in arms for s in arm_specs[k] for dr in chunks(a.n_boot)]
    het_tasks = [dict(arm=arms[k], idx=idx[k], warm=het[k]["fit"]["theta"], **draw_slice(dr))
                 for k in SOURCES for dr in chunks(a.n_boot)]
    gof_keys = [(k, s) for k in arms for s in arm_specs[k] if s in GOF_SPECS]
    gof_streams = dict(zip(gof_keys, stream["gof"].spawn(len(gof_keys))))
    gof_tasks = [dict(arm=arms[k], spec=s, warm=fits[(k, s)]["theta"],
                      sims=simulate_tables(arms[k], s, fits[(k, s)]["theta"], point_tables[(k, s)],
                                           np.random.default_rng(gof_streams[(k, s)]), a.n_gof)) for k, s in gof_keys]
    het_streams = dict(zip(SOURCES, stream["het_gof"].spawn(len(SOURCES))))
    gof_tasks += [dict(arm=arms[k], spec="Y5_het", warm=het[k]["fit"]["theta"],
                       sims=simulate_tables(arms[k], "Y5_het", het[k]["fit"]["theta"], point_tables[(k, "Y5")],
                                            np.random.default_rng(het_streams[k]), a.n_gof)) for k in SOURCES]
    with ProcessPoolExecutor(max_workers=max(1, a.workers)) as ex:
        gof_futs = [ex.submit(run_gof, t) for t in gof_tasks]
        boot_futs = [ex.submit(run_bootstrap, t) for t in boot_tasks]
        het_futs = [ex.submit(run_het_bootstrap, t) for t in het_tasks]
        gof_res = [f.result() for f in gof_futs]
        boot_rows = [r for f in boot_futs for r in f.result()]
        het_rows = [r for f in het_futs for r in f.result()]
    draws = pd.DataFrame(boot_rows).sort_values(["arm", "spec", "draw"]).reset_index(drop=True)
    hdraws = pd.DataFrame(het_rows).sort_values(["arm", "draw"]).reset_index(drop=True)
    log("bootstrap and goodness of fit done")

    gof_rows, gof = [], {}
    for g in gof_res:
        k = (g["arm"], g["spec"])
        d_obs = (het[g["arm"]]["deviance"] if g["spec"] == "Y5_het"
                 else point_deviance(arms[g["arm"]], g["spec"], fits[k]["theta"], point_tables[k]))
        p = float((1 + np.sum(g["devs"] >= d_obs - 1e-9)) / (len(g["devs"]) + 1))
        gof[k] = dict(deviance=d_obs, p=p, rejected=p < GOF_ALPHA)
        gof_rows.append(dict(arm=g["arm"], spec=g["spec"], deviance=round(d_obs, 4), n_sim=len(g["devs"]),
                             sim_deviance_median=round(float(np.median(g["devs"])), 4),
                             sim_deviance_95pct=round(float(np.percentile(g["devs"], 95)), 4),
                             p_parametric_bootstrap=p, model_rejected_at_0p05=p < GOF_ALPHA,
                             n_sim_refits_not_converged=g["n_not_converged"]))

    # ---- per arm/spec summaries
    results = {}
    rho_cols = lambda df: [c for c in df.columns if c == "rho" or c.startswith("rho_A")]
    for (name, spec), r in point.items():
        dr = draws[(draws["arm"] == name) & (draws["spec"] == spec)]
        num = [c for c in dr.columns if c not in ("arm", "spec", "draw", "converged", "at_bound")
               and pd.api.types.is_numeric_dtype(dr[c]) and dr[c].notna().any()]
        ci = {c: interval(dr[c]) for c in num}
        rc = rho_cols(dr)
        results.setdefault(name, dict(role=ARM_ROLE[name], n_pairs=(len(arms[name]["y1"]) if arms[name]["kind"] == "single"
                                                                    else [len(m["y1"]) for m in arms[name]["members"]]),
                                      recorded_within_pair_sd=recorded.get(name), specs={}))
        entry = dict(point=rnd(r), ci95=ci,
                     bootstrap=dict(n_draws=int(len(dr)), unit="pair (participant), resampled within arm; F target = "
                                    "P2 resample b with R2 recomputed; F_sealed individuals resampled within cohort",
                                    n_not_converged=int((~dr["converged"]).sum()), n_at_bound=int(dr["at_bound"].sum()),
                                    n_rho_negative=int((dr[rc].min(axis=1) < 0).sum()),
                                    n_var_e_below_0p01=int((dr["var_e"] < 0.01).sum()),
                                    n_draws_any_target_floored=int((dr["n_target_cohorts_floored"] > 0).sum())))
        if (name, spec) in gof:
            entry["goodness_of_fit"] = dict(deviance=round(gof[(name, spec)]["deviance"], 4),
                                            p_parametric_bootstrap=gof[(name, spec)]["p"],
                                            label="model rejected" if gof[(name, spec)]["rejected"] else "not rejected")
        if arms[name]["kind"] == "single":
            N = point_tables[(name, spec)]
            mu, sd, rho, d = unpack(fits[(name, spec)]["theta"])
            E = cell_probs(CUTS[SPECS[spec]["scale"]], mu, mu + d, sd, rho) * N.sum()
            k = N.shape[0]
            dd = np.arange(-(k - 1), k)
            entry["fit_check"] = dict(
                delta_values=dd.tolist(),
                observed_delta_counts=[int(sum(N[i, i + j] for i in range(k) if 0 <= i + j < k)) for j in dd],
                expected_delta_counts=[round(float(sum(E[i, i + j] for i in range(k) if 0 <= i + j < k)), 3) for j in dd],
                observed_first=N.sum(1).astype(int).tolist(), expected_first=[round(v, 3) for v in E.sum(1)],
                observed_repeat=N.sum(0).astype(int).tolist(), expected_repeat=[round(v, 3) for v in E.sum(0)])
        results[name]["specs"][spec] = entry
    for name in no_shift:
        results[name]["no_shift_sensitivity_Y5"] = rnd(no_shift[name])

    # ---- P2 consumer file: sources x {Y5, S3} x draw, F target
    fkeys = [f"F_{c}" for c in F_COHORTS]
    cons_cols = (["var_e", "rho", "converged", "at_bound"] + [f"lambda_{k}" for k in fkeys]
                 + [f"R2adj_{k}" for k in fkeys] + [f"lambdacond_{k}" for k in fkeys]
                 + ["lambda_hc_F", "lambda_within_F", "lambdacond_hc_F"])
    cons = draws[draws["arm"].isin(SOURCES) & draws["spec"].isin(SOURCE_SPECS)][["arm", "spec", "draw"] + cons_cols]
    cons = cons.rename(columns={"arm": "source"})
    cons.insert(2, "carried", cons.apply(lambda r: carried[r["source"]] == r["spec"], axis=1))
    lc = cons["lambdacond_hc_F"]
    cons["invalid_lambdacond_hc_F"] = ~np.isfinite(lc) | (lc <= INVALID_LAMBDA) | ~cons["converged"]

    def source_summary(src, spec):
        c = cons[(cons["source"] == src) & (cons["spec"] == spec)]
        v = c["lambdacond_hc_F"].to_numpy(float)
        inv = c["invalid_lambdacond_hc_F"].to_numpy(bool)
        r = point[(src, spec)]
        return dict(spec=spec, carried=carried[src] == spec, var_e_point=r["var_e"], var_e_ci95=interval(c["var_e"]),
                    lambda_hc_F_point=r["lambda_hc_F"], lambda_hc_F_ci95=interval(c["lambda_hc_F"]),
                    lambdacond_hc_F_point=r["lambdacond_hc_F"], lambdacond_hc_F_ci95=interval(v),
                    lambdacond_hc_F_draw_median=float(np.nanmedian(v)),
                    lambdacond_hc_F_p05_invalid_as_0=float(np.percentile(np.where(inv, 0.0, v), 5)),
                    n_invalid_draws=int(inv.sum()), n_draws=int(len(c)),
                    per_cohort={cc: dict(lambda_point=r[f"lambda_F_{cc}"], lambda_ci95=interval(c[f"lambda_F_{cc}"]),
                                         R2adj_observed=r2_obs[f"F_{cc}"],
                                         lambdacond_point=r[f"lambdacond_F_{cc}"],
                                         lambdacond_ci95=interval(c[f"lambdacond_F_{cc}"]),
                                         lambdacond_rawR2_point=float(lambda_conditional(r[f"lambda_F_{cc}"],
                                                                                         r2_raw_obs[f"F_{cc}"])))
                                for cc in F_COHORTS},
                    gof=("model rejected" if gof[(src, spec)]["rejected"] else "not rejected"),
                    gof_p=gof[(src, spec)]["p"], sd_check=("passed" if sd_ok[(src, spec)] else "failed"))

    src_summary = {src: {spec: source_summary(src, spec) for spec in SOURCE_SPECS} for src in SOURCES}

    # ---- heteroscedastic sensitivity summary and the prespecified rule
    het_tab = []
    for src in SOURCES:
        hd = hdraws[hdraws["arm"] == src]
        hp = het[src]["point"]
        cs = carried[src]
        cval = point[(src, cs)]["lambdacond_hc_F"]
        het[src]["summary"] = dict(
            lambda_hc_F_point=hp["lambda_hc_F"], lambda_hc_F_ci95=interval(hd["lambda_hc_F"]),
            lambdacond_hc_F_point=hp["lambdacond_hc_F"], lambdacond_hc_F_ci95=interval(hd["lambdacond_hc_F"]),
            b_point=hp["b"], b_ci95=interval(hd["b"]), gof_p=gof[(src, "Y5_het")]["p"],
            carried_spec=cs, carried_gof=src_summary[src][cs]["gof"], carried_lambdacond_hc_F=cval,
            het_below_carried=bool(hp["lambdacond_hc_F"] < cval),
            n_draws=int(len(hd)), n_not_converged=int((~hd["converged"]).sum()), n_at_bound=int(hd["at_bound"].sum()))
        for tgt, k in [(t["cohort"], t["key"]) for t in targets if t["set"] == "F"] + [("pool_headcount_F", "hc_F")]:
            het_tab.append(dict(source=src, target=tgt, lambda_point=hp[f"lambda_{k}"],
                                lambda_ci95_lo=interval(hd[f"lambda_{k}"])[0],
                                lambda_ci95_hi=interval(hd[f"lambda_{k}"])[1],
                                lambdacond_point=hp[f"lambdacond_{k}"],
                                lambdacond_ci95_lo=interval(hd[f"lambdacond_{k}"])[0],
                                lambdacond_ci95_hi=interval(hd[f"lambdacond_{k}"])[1],
                                homoscedastic_carried_lambdacond_point=point[(src, cs)][f"lambdacond_{k}"],
                                a=hp["a"], b=hp["b"], sd_e_at_0p5=hp["sd_e_at_0p5"], sd_e_at_1p5=hp["sd_e_at_1p5"],
                                sd_e_at_2p5=hp["sd_e_at_2p5"], sd_e_at_3p5=hp["sd_e_at_3p5"],
                                lrt_b0=het[src]["lrt_vs_homoscedastic"]["lrt"],
                                lrt_b0_p=het[src]["lrt_vs_homoscedastic"]["p_asymptotic"],
                                gof_p=gof[(src, "Y5_het")]["p"], sd_check_passed=het[src]["sd_check"]["passed"]))

    # ---- lambda table (long; every arm, spec and target)
    lam_rows = []
    for (name, spec), r in point.items():
        ci = results[name]["specs"][spec]["ci95"]
        gf = gof.get((name, spec))
        common = dict(arm=name, arm_role=ARM_ROLE[name], spec=spec,
                      carried=bool(name in SOURCES and carried[name] == spec), var_e_point=r["var_e"],
                      var_e_ci95_lo=ci["var_e"][0], var_e_ci95_hi=ci["var_e"][1],
                      gof_p=(gf["p"] if gf else np.nan), model_rejected=(gf["rejected"] if gf else ""),
                      sd_check_passed=sd_ok[(name, spec)])
        for t in targets:
            k = t["key"]
            row = dict(common, set=t["set"], target=t["cohort"], n_target=t["n"], lambda_point=r[f"lambda_{k}"],
                       ci95_lo=ci[f"lambda_{k}"][0], ci95_hi=ci[f"lambda_{k}"][1],
                       E_S=r[f"ES_{k}"], Var_S=r[f"VarS_{k}"], Cov_S1S2=r[f"CovS_{k}"],
                       Var_S_ci95=json.dumps(ci[f"VarS_{k}"]), Cov_S1S2_ci95=json.dumps(ci[f"CovS_{k}"]))
            if t["set"] == "F":
                row.update(R2_adjusted=r2_obs[k], R2_raw=r2_raw_obs[k], lambdacond_point=r[f"lambdacond_{k}"],
                           lambdacond_ci95_lo=ci[f"lambdacond_{k}"][0], lambdacond_ci95_hi=ci[f"lambdacond_{k}"][1],
                           lambdacond_rawR2_point=float(lambda_conditional(r[f"lambda_{k}"], r2_raw_obs[k])))
            lam_rows.append(row)
        for s in SETS:
            for pool, key in (("pool_headcount", f"lambda_hc_{s}"), ("pool_within", f"lambda_within_{s}")):
                row = dict(common, set=s, target=pool, n_target=sum(t["n"] for t in targets if t["set"] == s),
                           lambda_point=r[key], ci95_lo=ci[key][0], ci95_hi=ci[key][1])
                if s == "F" and pool == "pool_headcount":
                    row.update(lambdacond_point=r["lambdacond_hc_F"], lambdacond_ci95_lo=ci["lambdacond_hc_F"][0],
                               lambdacond_ci95_hi=ci["lambdacond_hc_F"][1])
                lam_rows.append(row)
        lam_rows.append(dict(common, set="deposit", target="in_deposit",
                             n_target=results[name]["n_pairs"] if isinstance(results[name]["n_pairs"], int)
                             else sum(results[name]["n_pairs"]),
                             lambda_point=r.get("lambda_in_deposit", np.nan),
                             ci95_lo=ci.get("lambda_in_deposit", [np.nan])[0],
                             ci95_hi=ci.get("lambda_in_deposit", [np.nan, np.nan])[1]))
    lam = pd.DataFrame(lam_rows)

    # ---- rho >= 0 sensitivity summary: how far the 2.5% bound moves
    rho_sens = {}
    for name, specs in arm_specs.items():
        base_s, nonneg = specs[0], specs[1]
        rho_sens[name] = dict(
            statistic="lambdacond_hc_F",
            unconstrained=dict(spec=base_s, point=point[(name, base_s)]["lambdacond_hc_F"],
                               ci95=results[name]["specs"][base_s]["ci95"]["lambdacond_hc_F"],
                               n_rho_negative=results[name]["specs"][base_s]["bootstrap"]["n_rho_negative"]),
            rho_nonneg=dict(spec=nonneg, point=point[(name, nonneg)]["lambdacond_hc_F"],
                            ci95=results[name]["specs"][nonneg]["ci95"]["lambdacond_hc_F"]))

    # ---- latent-shape sensitivity (sources, both specs, F cohorts)
    shape_rows = []
    for src in SOURCES:
        for spec in SOURCE_SPECS:
            r = point[(src, spec)]
            per = {}
            for t in [t for t in targets if t["set"] == "F"]:
                s_ = latent_shape_sensitivity(t["counts"], CUTS[t["scale"]], r["var_e"])
                per[t["key"]] = s_
                shape_rows.append(dict(source=src, spec=spec, var_e=r["var_e"], target=t["cohort"],
                                       lambda_normal=r[f"lambda_{t['key']}"], **s_))
            ft = [t for t in targets if t["set"] == "F"]
            n = np.array([t["n"] for t in ft], float)
            row = dict(source=src, spec=spec, var_e=r["var_e"], target="pool_headcount_F", lambda_normal=r["lambda_hc_F"],
                       npmle_lambda=float(n @ np.array([per[t["key"]]["npmle_lambda"] for t in ft]) / n.sum()))
            for sense in ("min", "max"):
                row[f"lp_{sense}_lambda"] = float(n @ np.array([per[t["key"]][f"lp_{sense}_lambda"] for t in ft]) / n.sum())
            row["lp_marginal"] = ";".join(sorted({per[t["key"]]["lp_marginal"] for t in ft}))
            shape_rows.append(row)
    shape = pd.DataFrame(shape_rows)
    log("latent shape done")

    # ---- Codex C8 comparison at the point
    codex = None
    if a.codex_lambda:
        cj = json.load(open(a.codex_lambda))
        diffs = []
        for src, sd_ in cj["sources"].items():
            if src not in arms:
                continue
            for spec, v in sd_["specifications"].items():
                for c, lam_c in v["per_cohort"].items():
                    diffs.append(dict(source=src, spec=spec, cohort=c, codex=lam_c, ours=point[(src, spec)][f"lambda_F_{c}"],
                                      var_e_codex=v["var_e"], var_e_ours=point[(src, spec)]["var_e"]))
                diffs.append(dict(source=src, spec=spec, cohort="headcount_mean", codex=v["headcount_mean"],
                                  ours=point[(src, spec)]["lambda_hc_F"], var_e_codex=v["var_e"],
                                  var_e_ours=point[(src, spec)]["var_e"]))
        dd = pd.DataFrame(diffs)
        dd["abs_diff"] = (dd["codex"] - dd["ours"]).abs()
        src_max = dd[dd["source"].isin(SOURCES)]["abs_diff"].max()
        codex = dict(path=a.codex_lambda, max_abs_diff_sources=float(src_max),
                     max_abs_diff_all=float(dd["abs_diff"].max()),
                     by_source_spec=dd.groupby(["source", "spec"])["abs_diff"].max().round(8).reset_index()
                     .to_dict("records"),
                     codex_carried={s: cj["sources"][s]["conservative_spec_by_headcount"] for s in SOURCES},
                     codex_ci95={s: cj["sources"][s].get("ci95") for s in SOURCES})
        check("Codex_C8_lambda_558_agrees_at_the_sources", src_max < 1e-4,
              f"per-cohort and headcount lambda for GSE193066 and PXD051911_A1 (Y5, S3) vs {a.codex_lambda}: max |diff| "
              f"{src_max:.2e}; all arms {dd['abs_diff'].max():.2e}; carried specs ours {carried}, Codex "
              f"{codex['codex_carried']}")

    # ---- statements (numbers read from this run)
    gsd = {c: float(np.sqrt(((g_pairs.loc[g_pairs["identity_class"] == c, "y2"] -
                              g_pairs.loc[g_pairs["identity_class"] == c, "y1"]) ** 2).sum() /
                            (2.0 * (g_pairs["identity_class"] == c).sum()))) for c in ("same", "different", "uncheckable")}
    gexact = {c: int((g_pairs.loc[g_pairs["identity_class"] == c, "y1"] ==
                      g_pairs.loc[g_pairs["identity_class"] == c, "y2"]).sum()) for c in ("same", "different", "uncheckable")}
    ident = {}
    for n_ in ("GSE193066", "GSE193066_drop11", "GSE193066_same29"):
        cs = max(SOURCE_SPECS, key=lambda s: point[(n_, s)]["var_e"])
        ident[n_] = dict(spec=cs, lambdacond_hc_F=point[(n_, cs)]["lambdacond_hc_F"],
                         ci95=results[n_]["specs"][cs]["ci95"]["lambdacond_hc_F"])
    pxd_first = {n_: int(point_tables[(n_, "Y5")].sum(1)[2:].sum()) for n_ in ("PXD051911_A1", "PXD051911_A2")}
    rejected = sorted(f"{k[0]} {k[1]} (p {v['p']:.3f})" for k, v in gof.items() if v["rejected"])
    tg_rej = [f"{r['set']} {r['cohort']} (X2 {r['gof_pearson_X2']:.2f}, bootstrap p {r['gof_p_parametric_bootstrap']:.3f})"
              for r in tc_rows if r["gof_df"] > 0 and r["gof_p_parametric_bootstrap"] < GOF_ALPHA]
    ss = {s: src_summary[s][carried[s]] for s in SOURCES}
    statements = [
        "Estimand: lambda corrects kappa_U to a slope per unit of the classical true score E[S|T]. It is not the "
        "slope per step of the error-free stage cut(T).",
        "Target: F, the 558 staged development unit individuals of spec v1 1.1 (GSE130970 54, GSE135251 160, "
        "GSE162694 83, GSE213621 208, GSE240729 53). Each cohort's latent stage distribution is fitted from F's own "
        "labels (0-4 Kleiner stage; GSE213621 the 3-level condition label, exactly identified). F_sealed (237) is "
        "reported for the sealed run only (spec 6), marginal lambda.",
        "Carried spec per source (spec 5.2): the smaller of Y5 and S3 at the point estimate, fixed across draws. "
        "lambda_c falls as var_e rises for any fixed target distribution, so the carried spec is the one with the "
        f"larger var_e and is smaller in every cohort and pool: GSE193066 {carried['GSE193066']}, PXD051911_A1 "
        f"{carried['PXD051911_A1']}.",
        "Covariates: lambda_cond_c = (lambda_c - R2_c) / (1 - R2_c) with R2_c the adjusted within-cohort R2 of S on "
        "the spec z columns and the two ancestry axes, floored at 0 (v1.2 items 5-7), reproduced to 1e-10 from the "
        "lead's r2_moments_observed.json. It assumes the label error is independent of z. The raw-R2 value is a "
        "sensitivity (lambdacond_rawR2_point).",
        "Pools: lambdacond_hc_F is the spec 5.2 pool with I_c replaced by the headcount n_c; B-MODEL re-pools the "
        "per-cohort lambdacond_F_<cohort> columns with the Gate A I_c.",
        "Draws: draw b uses pair resample b for var_e and P2 resample b (p2_resample_development.tsv.gz, "
        "SeedSequence([20260923, 2, b])) for the target fits and R2, so lambda*_b is index-aligned with kappa*_b. "
        "R2 on a with-replacement resample is larger than the observed adjusted R2 (duplicated rows fit better), "
        "which lowers lambda*_b relative to the point; see r2_resample_shift.",
        "Invalid draws (spec 5.2, v1.2 item 1): lambdacond_hc_F non-finite, <= 0.05 or unconverged. They are flagged "
        "in invalid_lambdacond_hc_F and entered as 0 for the 5th percentile.",
        "GSE213621 stage is a 3-level condition label per library. Genotypes place 25 same-person library pairs in it "
        f"and {len(disc)} carry two different labels (docs/STATUS.md GSE213621 row); the Model A unit set drops the "
        "discordant persons. Its lambda uses the transported var_e, the same assumption as every other target. The "
        "GSE213621_pairs25 arm estimates var_e from the 25 pairs "
        f"(point {point[('GSE213621_pairs25', 'S3_sym')]['var_e']:.3f}); it is a selected set, and its discordance "
        "suggests the transported var_e may understate GSE213621 label error.",
        "Each arm is a test-retest across two physically distinct biopsies, so var_e combines reader variability, "
        "biopsy sampling and real change between the biopsies. lambda is a lower bound on the single-occasion "
        "reliability relevant to Model A only if target reads are no noisier than the arm's and the between-biopsy "
        "change is independent of T. Biopsy sampling error applies in part to the RNA-bearing sample as well.",
        "Sources are not pooled. PXD051911 enters through arm A1 (spec 5.2; PRIMARY in the replication prespec, "
        "kleiner-reliability-replication-pxd051911-20260901T211358Z/prespec/PRESPECIFICATION.md:107). A2 and the "
        "joint A1+A2 fit are sensitivities. The joint fit's implied within-pair SD misses A1 or A2 by more than 0.02 "
        f"(Y5 {'passed' if sd_ok[('PXD051911_joint', 'Y5')] else 'failed'}, S3 "
        f"{'passed' if sd_ok[('PXD051911_joint', 'S3')] else 'failed'}), so it is not carriable; shared-var_e LRT Y5 "
        f"{joint_lrt['Y5']['lrt']:.2f} (1 df, p {joint_lrt['Y5']['p_asymptotic']:.3f}).",
        f"PXD051911 has few first reads at F2 or above (A1 {pxd_first['PXD051911_A1']}, A2 "
        f"{pxd_first['PXD051911_A2']}), so its error at the S cutpoints 1.5 and 2.5 rests mostly on the F0/F1 boundary "
        "under a homoscedastic equal-interval latent. The S3 fit uses only those cutpoints; Y5_het lets the error SD "
        "change with the latent level.",
        "GSE193066 identity (provisional calls from the frozen-crosswalk individual_id; reproduces docs/STATUS.md): "
        f"{ic.get('same', 0)} of the 58 pairs are called one person, {ic.get('different', 0)} pairs are called as two "
        f"different people, {ic.get('uncheckable', 0)} cannot be checked. Within-pair SD: called one person "
        f"{gsd['same']:.4f} ({gexact['same']} exact agreements), called two people {gsd['different']:.4f} "
        f"({gexact['different']}), uncheckable {gsd['uncheckable']:.4f} ({gexact['uncheckable']}). lambdacond_hc_F by "
        "identity status (carried spec each): "
        + ", ".join(f"{k} {v['spec']} {v['lambdacond_hc_F']:.3f}" for k, v in ident.items())
        + ". Provisional until the contamination-aware identity rerun; spec Appendix B2 applies if GSE193066 is "
          "withdrawn.",
        "rho < 0 draws (var_T < 0 in the source) are kept in the primary intervals; they raise var_e and lower "
        "lambda, which is conservative for P2. The rho >= 0 fits are the admissible-model sensitivity "
        "(rho_nonneg_sensitivity). Percentile intervals next to the var_e = 0 boundary are approximate.",
        "Goodness of fit: parametric bootstrap of the deviance per arm and spec, refit each draw, rejection at p < "
        f"0.05. Rejected: {', '.join(rejected) if rejected else 'none'}. Rule fixed before this run: when a source's "
        "carried spec is rejected, the Y5_het fit is its required sensitivity. The carried var_e draws stay the P2 "
        "input (the heteroscedastic model has no single var_e); if Y5_het gives a smaller lambdacond_hc_F than the "
        "carried value, that is flagged to B-SPEC as the less favourable value.",
        "Target latent shape: the normal single-read fit is rejected at p < 0.05 (parametric bootstrap, 2 df) in "
        f"{', '.join(tg_rej) if tg_rej else 'no cohort'}. latent_shape_sensitivity.tsv gives lambda under a "
        "nonparametric latent (grid NPMLE) and the range over all latent distributions that reproduce each cohort's "
        "observed labels; per-cohort lambda is not to be published without it.",
        "lambda depends on the stage distribution it is evaluated at; the transported quantity is var_e on the "
        "latent stage scale, not a reliability coefficient.",
        "Carried values at F (headcount pool, lambda_cond): "
        + "; ".join(f"{s} {ss[s]['spec']} {ss[s]['lambdacond_hc_F_point']:.4f} [{ss[s]['lambdacond_hc_F_ci95'][0]:.4f}, "
                    f"{ss[s]['lambdacond_hc_F_ci95'][1]:.4f}], invalid draws {ss[s]['n_invalid_draws']}"
                    for s in SOURCES) + ".",
    ]

    report = dict(
        purpose="Model A spec v1 5.2 lambda: reliability of the harmonized 3-level stage S from label-only repeat "
                "pairs, per source, at the development fit set F; sources never pooled.",
        seed=SEED, n_boot=int(a.n_boot), n_gof=int(a.n_gof),
        model=dict(latent="Z_v = T + e_v (+ delta at the repeat visit), T ~ N(mu, var_T), e_v ~ N(0, var_e)",
                   cutpoints=CUTS_Y.tolist(), S_cutpoints=CUTS_S.tolist(), specs=SPECS,
                   het="Y5_het: reads conditionally independent given T ~ N(mu, sd_T^2); error SD exp(a + b (T - 2)) "
                       "at the level read (T, or T + delta at the repeat); b = 0 is Y5_rho_nonneg",
                   rho="var_T / (var_T + var_e); allowed in (-1, 1) except in the *_rho_nonneg specs",
                   joint="PXD051911_joint: shared var_e; arm-specific mu, sd_Z (rho_a = 1 - var_e / sd_Z,a^2), delta",
                   implied_within_pair_sd="sqrt(E[(Y2 - Y1)^2] / 2) on the fit's scale (Y 0-4 for Y5, S 0-2 for S3)",
                   lambda_definition="Var(E[S|T]) / Var(S) at the cohort's latent distribution with the arm's var_e "
                                     "(var_T,c = sd_Z,c^2 - var_e, floored at 0): kappa per unit of E[S|T]",
                   lambda_conditional="(lambda_c - R2_c) / (1 - R2_c), R2_c adjusted, floored at 0 (v1.2 item 6)",
                   pools=dict(lambda_hc="sum_c n_c lambda_c / sum_c n_c",
                              lambdacond_hc="sum_c n_c lambda_cond_c / sum_c n_c (spec 5.2 with I_c -> n_c)",
                              lambda_within="sum_c n_c Cov_c(S1, S2) / sum_c n_c Var_c(S)")),
        guards=guards, guards_failed=[g["guard"] for g in guards if not g["passed"]],
        checks=checks, checks_failed=[c["check"] for c in checks if not c["passed"]],
        checks_consequence="checks are not exit-status guards; a fit that fails its SD check is labelled "
                           "'SD check failed' and is not carriable",
        carried_spec=carried_rows,
        sources=src_summary,
        heteroscedastic_sensitivity={s: dict(summary=het[s]["summary"], lrt_vs_homoscedastic=het[s]["lrt_vs_homoscedastic"],
                                             sd_check=het[s]["sd_check"], point=rnd(het[s]["point"]))
                                     for s in SOURCES},
        pxd051911_shared_var_e_lrt=joint_lrt,
        goodness_of_fit=dict(rule=f"parametric bootstrap, {a.n_gof} draws from the fitted model at the arm's n, "
                                  "refit each; p = (1 + #{D_sim >= D_obs}) / (draws + 1); rejected if p < 0.05",
                             arms=gof_rows),
        r2=dict(observed=r2_rows, resample_shift=r2_shift),
        arms=results,
        rho_nonneg_sensitivity=rho_sens,
        GSE193066_identity_range=ident,
        codex_c8_comparison=codex,
        targets=dict(rule="Frozen seal crosswalk (sha256 checked): unit_library True, C = 0, S observed; F = run_role "
                          "development in the five F cohorts; F_sealed = run_role sealed in those plus GSE174478. "
                          "Kleiner cohorts: Y = metadata fibrosis_stage, checked against S; GSE213621: S checked "
                          "against the condition label.",
                     counts=tcounts, cohorts=tc_rows, p2_resamples=p2_info),
        statements=statements,
        inputs={k: dict(path=str(p), sha256=sha256(p)) for k, p in {
            "gse193066_placement": Path(a.gse193066_placement),
            "gse193066_pairs": Path(a.gse193066_source) / "results/pairs.tsv",
            "gse193066_json": Path(a.gse193066_source) / "results/kleiner_reliability.json",
            "gse193066_series": Path(a.gse193066_series), "gse193066_crosswalk": Path(a.gse193066_crosswalk),
            "pxd_pairs": Path(a.pxd_source) / "results/pairs.tsv",
            "pxd_json": Path(a.pxd_source) / "results/kleiner_replication.json",
            "pxd_meta": Path(a.pxd_meta), "frozen_crosswalk": Path(a.frozen_crosswalk),
            "kinship_pairs": Path(a.kinship_pairs), "metadata": Path(a.metadata),
            "design_individuals": Path(a.design_individuals), "ancestry": Path(a.ancestry),
            "r2_moments": Path(a.r2_moments), "p2_resamples": Path(a.p2_resamples),
            "r2_by_resample": Path(a.r2_by_resample), "script": Path(__file__).resolve(),
            **({"codex_lambda": Path(a.codex_lambda)} if a.codex_lambda else {})}.items()},
        environment=dict(python=sys.version.split()[0], workers=int(a.workers),
                         pip_freeze=subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True,
                                                   text=True).stdout.splitlines()),
    )

    npz = {}
    for (name, spec), dr in pd.concat([draws, hdraws]).groupby(["arm", "spec"]):
        dr = dr.sort_values("draw")
        for col in dr.columns.drop(["arm", "spec"]):
            if (pd.api.types.is_numeric_dtype(dr[col]) or pd.api.types.is_bool_dtype(dr[col])) and dr[col].notna().any():
                npz[f"{name}__{spec}__{col}"] = dr[col].to_numpy()

    cons.to_csv(out / "lambda_source_draws.tsv.gz", sep="\t", index=False)
    pd.DataFrame(carried_rows).to_csv(out / "carried_spec.tsv", sep="\t", index=False)
    lam.to_csv(out / "lambda_by_cohort.tsv", sep="\t", index=False)
    pd.DataFrame(het_tab).to_csv(out / "het_sensitivity.tsv", sep="\t", index=False)
    shape.to_csv(out / "latent_shape_sensitivity.tsv", sep="\t", index=False)
    pd.DataFrame(gof_rows).to_csv(out / "goodness_of_fit.tsv", sep="\t", index=False)
    pd.DataFrame(table_rows).to_csv(out / "pair_tables.tsv", sep="\t", index=False)
    pd.DataFrame([{k: v for k, v in r.items()} for r in tc_rows]).to_csv(out / "target_stage_counts.tsv", sep="\t",
                                                                         index=False)
    pd.concat([draws, hdraws]).to_csv(out / "bootstrap_draws.tsv.gz", sep="\t", index=False)
    np.savez_compressed(out / "bootstrap_draws.npz", **npz)
    pd.concat(pairs_out)[["arm", "participant", "run1", "run2", "y1", "y2", "identity_class"]].to_csv(
        rout / "pairs_used.tsv", sep="\t", index=False)
    tp.to_csv(rout / "target_participants.tsv", sep="\t", index=False)
    with open(out / "label_reliability.json", "w") as fh:
        json.dump(report, fh, indent=1, default=lambda v: v.item() if hasattr(v, "item") else
                  (v.tolist() if isinstance(v, np.ndarray) else str(v)))
    log("written")

    print(json.dumps(dict(guards_failed=report["guards_failed"], checks_failed=report["checks_failed"],
                          gof={f"{r['arm']}|{r['spec']}": r["p_parametric_bootstrap"] for r in gof_rows},
                          carried=carried,
                          sources={s: {k: v for k, v in ss[s].items() if k != "per_cohort"} for s in SOURCES},
                          het={s: het[s]["summary"] for s in SOURCES},
                          codex=(codex and dict(max_abs_diff_sources=codex["max_abs_diff_sources"]))),
                     indent=1, default=str))
    return 1 if report["guards_failed"] else 0


if __name__ == "__main__":
    sys.exit(main())
