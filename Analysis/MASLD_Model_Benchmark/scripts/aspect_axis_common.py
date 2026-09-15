#!/usr/bin/env python3
"""Shared instrument for the four-arm aspect-axis decomposition.

One instrument, four arms. Every arm is loaded to the same contract -- a
participant axis, an aspect table on the participant axis, and a feature matrix
on the participant axis -- so the per-arm decompositions are comparable and the
arms are never pooled at the feature level.

The four arms are FOUR DIFFERENT INSTRUMENTS (bulk RNA on two scales, liver
protein, single-cell pseudobulk). Nothing here concatenates their features or
their participants. The only pooling is of decision statistics.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import pathlib
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata
from scipy.stats import t as student_t

ROOT = pathlib.Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BENCH = ROOT / "Analysis/MASLD_Model_Benchmark"


# --------------------------------------------------------------------------
# integrity helpers
# --------------------------------------------------------------------------
def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_sha256(obj) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


class JoinError(RuntimeError):
    pass


def require_join(observed: int, expected: int, label: str) -> int:
    """A join that lands on the wrong number of rows raises. A zero join always
    raises. This is the guard that fires when a key format silently changes."""
    if observed == 0:
        raise JoinError(f"{label}: join produced ZERO rows (expected {expected})")
    if observed != expected:
        raise JoinError(f"{label}: join produced {observed} rows, expected {expected}")
    return observed


def read_tsv_strict(path, *, sep="\t", label=None, **kwargs) -> pd.DataFrame:
    """Read a delimited table after measuring its field count per line.

    Ragged headers are endemic in this project: an unnamed pandas index makes the
    header one field shorter than the body, and every positional reader then
    returns the NEIGHBOURING column. This measures NF on the raw bytes before
    pandas is allowed to align anything, and records the finding.
    """
    label = label or str(path)
    opener = gzip.open if str(path).endswith(".gz") else open
    widths = []
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i > 2000:
                break
            widths.append(line.rstrip("\n").count(sep) + 1)
    header_width = widths[0]
    body = set(widths[1:]) if len(widths) > 1 else set()
    ragged = bool(body) and (header_width not in body or len(body) > 1)
    frame = pd.read_csv(path, sep=sep, **kwargs)
    frame.attrs["nf_header"] = header_width
    frame.attrs["nf_body"] = sorted(body)
    frame.attrs["ragged_header"] = ragged
    if ragged and header_width + 1 in body:
        raise JoinError(
            f"{label}: ragged header ({header_width} vs {sorted(body)}) -- this is "
            "an unnamed pandas index; resolve by name, never by position")
    return frame


# --------------------------------------------------------------------------
# tie structure and power, computed from labels only
# --------------------------------------------------------------------------
def tie_ceiling(values) -> float:
    """Largest |Spearman| an ordinal label can reach against ANY untied
    continuous variable, given only its tie structure.

    Spearman with ties is the Pearson correlation of midranks. For an untied x
    whose ordering respects the label's blocks, cov(R_y, R_x) = var(R_y), so the
    bound is sd(midrank_y) / sd(1..n). It depends on the marginal alone.
    """
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    n = v.size
    if n < 3:
        return float("nan")
    ry = rankdata(v)
    sd_y = ry.std()
    sd_x = np.arange(1, n + 1, dtype=float).std()
    if sd_x == 0:
        return float("nan")
    return float(sd_y / sd_x)


def tie_ceiling_from_counts(counts) -> float:
    """Same bound expressed from a marginal level-count vector, so the pinned
    published constants can be reproduced without touching participant data."""
    counts = [int(c) for c in counts if int(c) > 0]
    values = np.concatenate([np.full(c, i, float) for i, c in enumerate(counts)])
    return tie_ceiling(values)


def mde_spearman(n: int, k_covariates: int, alpha: float, power: float = 0.80
                 ) -> float:
    """Smallest |partial Spearman| detectable at the stated alpha and power.

    Fisher z on the partial correlation with df = n - 3 - k. Collinearity among
    the aspects does not enter here: the partial correlation is already the
    standardised residual association, so multicollinearity moves what is
    ATTAINABLE (reported separately as the residual rank-variance fraction), not
    what is detectable.
    """
    df = n - 3 - k_covariates
    if df <= 0:
        return float("nan")
    z = norm.isf(alpha / 2.0) + norm.ppf(power)
    return float(np.tanh(z / np.sqrt(df)))


# --------------------------------------------------------------------------
# association machinery
# --------------------------------------------------------------------------
def rank_columns(X: np.ndarray) -> np.ndarray:
    """Column-wise rank transform, average ties, NaN-free input required."""
    return np.apply_along_axis(rankdata, 0, X).astype(float)


def residualise(Y: np.ndarray, C: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(C, Y, rcond=None)
    return Y - C @ beta


def standardise(Y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    s = Y.std(axis=0)
    good = s > 0
    Z = np.zeros_like(Y)
    Z[:, good] = (Y[:, good] - Y[:, good].mean(axis=0)) / s[good]
    return Z, good


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    v = p[ok]
    n = v.size
    if n == 0:
        return q
    order = np.argsort(v)
    ranked = v[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0.0, 1.0)
    q[ok] = out
    return q


def rho_to_p(rho: np.ndarray, df: int) -> np.ndarray:
    r = np.clip(rho, -0.999999, 0.999999)
    tstat = r * np.sqrt(df / (1.0 - r ** 2))
    return 2.0 * student_t.sf(np.abs(tstat), df)


@dataclass
class DirectionResult:
    exposure: str
    adjusted_for: list
    n: int
    family_size: int
    df: int
    count_bh05: int
    max_abs_rho: float
    residual_rank_variance_fraction: float
    null_draws: int
    null_count_mean: float
    null_count_p95: float
    null_count_max: int
    perm_p: float
    null_maxrho_p95: float
    verdict: str
    secondary_axis_permutation: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    def as_dict(self):
        d = dict(self.__dict__)
        return d


COLLINEARITY_TOLERANCE = 1e-10

# A BH family in the thousands against a permutation null that is a point mass
# at zero makes ANY count of one look extreme: the permutation p is then an
# artefact of the null's discreteness, not evidence. Five is the smallest count
# at which the criterion cannot be carried by a single feature. This is a
# judgment call and it is prespecified rather than discovered.
MIN_DISCOVERY_COUNT = 5


def partial_association(Zfeat_ranked: np.ndarray, x_rank: np.ndarray,
                        C: np.ndarray):
    """Return (rho, a_res_standardised, residual_variance_fraction, ok_mask).

    Zfeat_ranked is the rank-transformed feature matrix. Both it and the
    exposure are residualised on the SAME design C, then correlated. C must
    already carry an intercept column.
    """
    a = residualise(x_rank.reshape(-1, 1), C).ravel()
    x_centred = x_rank - x_rank.mean()
    denom = float((x_centred ** 2).sum())
    frac = float((a ** 2).sum() / denom) if denom > 0 else float("nan")
    # Exact collinearity never lands on exactly zero in floating point. An
    # exposure that the design already determines leaves a residual whose
    # variance is a rounding artefact, and correlating with it manufactures a
    # count out of noise. The tolerance is on the FRACTION of the exposure's
    # rank variance that survives, so it is scale free.
    sa = a.std()
    if not np.isfinite(frac) or frac < COLLINEARITY_TOLERANCE or sa == 0:
        return None, None, frac, None
    a = (a - a.mean()) / sa
    R = residualise(Zfeat_ranked, C)
    Zs, ok = standardise(R)
    rho = np.full(Zfeat_ranked.shape[1], np.nan)
    rho[ok] = (Zs[:, ok] * a[:, None]).sum(axis=0) / len(a)
    return rho, a, frac, ok


def _count_and_max(rho, df, ok):
    p = np.full(rho.shape, np.nan)
    p[ok] = rho_to_p(rho[ok], df)
    q = bh(p)
    return int(np.nansum(q < 0.05)), float(np.nanmax(np.abs(rho[ok]))) if ok.any() else float("nan")


def run_direction(Zfeat_ranked: np.ndarray, x_rank: np.ndarray, C: np.ndarray,
                  *, exposure: str, adjusted_for: list, rng: np.random.Generator,
                  n_perm: int, block: int = 250,
                  x_raw_for_secondary: np.ndarray | None = None,
                  C_secondary_builder=None) -> DirectionResult:
    """Freedman-Lane permutation of the residualised exposure.

    Permuting rows of the Z-residualised FEATURE matrix and permuting the
    Z-residualised EXPOSURE are algebraically the same operation on the
    correlation, so the cheap form is used. This preserves the exposure's
    dependence on the adjustment set, which permuting the RAW exposure does not.
    The raw-exposure null is still computed, as a secondary, because the frozen
    2026-08-27 runs used it and the two must be comparable.
    """
    n = Zfeat_ranked.shape[0]
    k = C.shape[1] - 1
    df = n - 2 - k
    rho, a, frac, ok = partial_association(Zfeat_ranked, x_rank, C)
    if rho is None:
        return DirectionResult(
            exposure=exposure, adjusted_for=list(adjusted_for), n=n,
            family_size=int(Zfeat_ranked.shape[1]), df=df, count_bh05=0,
            max_abs_rho=float("nan"), residual_rank_variance_fraction=frac,
            null_draws=0, null_count_mean=float("nan"),
            null_count_p95=float("nan"), null_count_max=0, perm_p=float("nan"),
            null_maxrho_p95=float("nan"), verdict="NOT_APPLICABLE_COLLINEAR",
            notes=["the exposure is exactly determined by the adjustment set; "
                   "the residual is identically zero and no test exists"])
    count, maxrho = _count_and_max(rho, df, ok)

    R = residualise(Zfeat_ranked, C)
    Zs, okf = standardise(R)
    Zs = Zs[:, okf]
    null_counts = np.empty(n_perm, int)
    null_max = np.empty(n_perm, float)
    done = 0
    while done < n_perm:
        b = min(block, n_perm - done)
        A = np.empty((n, b))
        for j in range(b):
            A[:, j] = a[rng.permutation(n)]
        RH = (Zs.T @ A) / n                              # family x b
        P = rho_to_p(RH, df)
        for j in range(b):
            q = bh(P[:, j])
            null_counts[done + j] = int(np.sum(q < 0.05))
            null_max[done + j] = float(np.max(np.abs(RH[:, j])))
        done += b
    perm_p = float((1 + int((null_counts >= count).sum())) / (1 + n_perm))

    secondary = {}
    if x_raw_for_secondary is not None and C_secondary_builder is not None:
        sc = np.empty(min(n_perm, 1000), int)
        for j in range(sc.size):
            xp = x_raw_for_secondary[rng.permutation(n)]
            r2, a2, _, ok2 = partial_association(
                Zfeat_ranked, rankdata(xp), C_secondary_builder(xp))
            sc[j] = 0 if r2 is None else _count_and_max(r2, df, ok2)[0]
        secondary = {
            "scheme": "permute the RAW exposure then re-residualise (the frozen "
                      "2026-08-27 scheme; it does not preserve the exposure's "
                      "dependence on the adjustment set)",
            "n_draws": int(sc.size),
            "null_count_mean": float(sc.mean()),
            "null_count_p95": float(np.percentile(sc, 95)),
            "perm_p": float((1 + int((sc >= count).sum())) / (1 + sc.size)),
        }

    p95 = float(np.percentile(null_counts, 95))
    unrestricted = "UNIQUE" if (perm_p < 0.05 and count > p95) else None
    if perm_p < 0.05 and count > p95 and count >= MIN_DISCOVERY_COUNT:
        verdict = "UNIQUE"
    elif unrestricted == "UNIQUE":
        verdict = "UNIQUE_BUT_BELOW_MINIMUM_COUNT"
    elif count == 0 or count <= p95:
        verdict = ("EMPTY_UNDERPOWERED"
                   if maxrho < float(np.percentile(null_max, 95))
                   else "EMPTY_NOT_INDEPENDENT")
    else:
        verdict = "NOT_UNIQUE"
    return DirectionResult(
        exposure=exposure, adjusted_for=list(adjusted_for), n=n,
        family_size=int(Zs.shape[1]), df=df, count_bh05=count,
        max_abs_rho=maxrho, residual_rank_variance_fraction=frac,
        null_draws=n_perm, null_count_mean=float(null_counts.mean()),
        null_count_p95=p95, null_count_max=int(null_counts.max()),
        perm_p=perm_p, null_maxrho_p95=float(np.percentile(null_max, 95)),
        verdict=verdict, secondary_axis_permutation=secondary,
        notes=([] if verdict != "UNIQUE_BUT_BELOW_MINIMUM_COUNT" else
               [f"count {count} clears its permutation null but is below the "
                f"prespecified minimum of {MIN_DISCOVERY_COUNT}; the verdict "
                "would rest on fewer than five features"]))


def stouffer(p_values, weights, *, n_perm=None):
    """Weighted Stouffer over per-arm permutation p values.

    A permutation p is bounded away from 0 and 1 by construction: with B draws
    it lies in [1/(B+1), 1]. An unclamped p of exactly 1 sends Fisher z to
    -inf and destroys the combination, so both tails are clamped to the
    resolution the permutation actually has.
    """
    p = np.asarray(p_values, float)
    w = np.asarray(weights, float)
    if n_perm:
        lo = 1.0 / (n_perm + 1.0)
        p = np.clip(p, lo, 1.0 - lo)
    else:
        p = np.clip(p, 1e-6, 1.0 - 1e-6)
    z = norm.isf(p)
    combined = float((w * z).sum() / np.sqrt((w ** 2).sum()))
    return combined, float(norm.sf(combined)), (w * z / np.sqrt((w ** 2).sum())).tolist()
