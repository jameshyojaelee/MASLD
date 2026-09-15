#!/usr/bin/env python3
"""Resolve pretraining exposure in EXPRESSION space, because the metadata cannot.

Foundation-model pretraining corpora are routinely published de-identified. The corpus
this method was developed against ships 581,503 x 20,010 expression values with an `obs`
of exactly two columns -- a row index and a train/test split -- and an empty `uns`. There
is no GSE, GSM, SRR, PRJNA or E-MTAB anywhere in the deposit. So "our evaluation cohort
does not appear in the corpus metadata" is not a finding: there is no namespace for it to
appear in. Today essentially every benchmark of a pretrained encoder simply assumes the
answer to "was my evaluation set in this model's pretraining data?"

What such a deposit DOES carry is the expression matrix itself. The question is therefore
answerable by measurement, and this tool measures it.

THE STATISTIC. Correlate every query sample against every corpus row; per query sample
keep the argmax corpus row and the max r. The maximum ALONE cannot decide anything: the
query's quantification pipeline is not the corpus's, so a same-sample pair is attenuated,
and two samples of the same tissue are highly correlated anyway. In the reference
application the corpus's own nearest-neighbour background had median r 0.9878, sitting
right on top of where the query samples landed.

What decides it is a TWO-LEG criterion, and both legs are load-bearing:

  LEG 1, STRUCTURE. A study copied into a corpus occupies ONE contiguous run of about its
  own size and maps ONE-TO-ONE, because corpora are built by concatenating studies. Two
  order-free quantities capture that: `injectivity`, the fraction of samples matching
  DISTINCT corpus rows -- nearest neighbours to a foreign study collide, several query
  samples picking the same row -- and `window_over_n`, the narrowest window of corpus rows
  holding 90% of the matches divided by the cohort size. A merely similar study scatters
  its matches wherever its own neighbours happen to sit.

  LEG 2, MAGNITUDE. The median max correlation must fall inside a same-sample band that is
  CALIBRATED FROM PROVEN BIJECTIONS, not invented. This is the subtle part. The band is
  measured on cohorts the offset test has already proven are in the corpus; those samples
  are literally the same libraries the query pipeline quantified, so their correlation IS
  the cross-pipeline attenuation for this pipeline pair. No threshold is guessed anywhere,
  and if no bijection is proven the tool reports that the band is NOT CALIBRABLE and calls
  nothing, rather than falling back on a number.

THE BIJECTION TEST. The largest set of samples sharing a single constant offset
(corpus_row minus rank within the cohort). Its null permutes the QUERY ORDER, which
destroys the correspondence while preserving every correlation and every matched corpus
row -- so a large aligned block cannot be manufactured by the cohort merely being similar
to a dense corpus neighbourhood. The reported bar is the EXCEEDANCE COUNT at a fixed draw
count. A null maximum is not a bar: it grows with the number of draws.

WHAT THIS CAN AND CANNOT SETTLE. A high match under both legs is strong evidence of
inclusion. A low match is "no near-duplicate detected" and NOT "not in the corpus" -- the
exposure state stays `unknown`, which is a state the vocabulary already has and which is
never folded into clean. Read MODEL_CARD.md before quoting anything from the output.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from corpus_io import open_corpus

HERE = Path(__file__).resolve().parent
TOOL = "pretraining-exposure-v1"
DEFAULT_SEED = 1103

NEGATIVE_READING = (
    "A cohort not called encoder_seen is reported `unknown`: no near-duplicate was "
    "detected. That is NOT evidence the cohort is absent from the corpus. Differing "
    "quantification pipelines attenuate the correlation, and a corpus carrying no "
    "accessions can never yield a `clean_declared` state."
)


# --------------------------------------------------------------------------------------
# The state vocabulary is IMPORTED, never restated. A convention copied instead of
# imported is a defect that no reading of either copy can detect.
# --------------------------------------------------------------------------------------
def load_exposure_vocabulary() -> tuple[object, str]:
    """Import `exposure.py` -- the per-study exposure states and `disposition()`.

    Searched in order: $EXPOSURE_VOCABULARY, a copy beside this file, the sibling
    encoder-benchmark release. Raises with every path tried rather than defining a
    second, silently divergent copy of the vocabulary.
    """
    tried = []
    cands = []
    env = os.environ.get("EXPOSURE_VOCABULARY")
    if env:
        cands.append(Path(env))
    cands += [HERE / "exposure.py", HERE.parent / "encoder-benchmark-v1" / "exposure.py"]
    for c in cands:
        tried.append(str(c))
        if c.is_file():
            spec = importlib.util.spec_from_file_location("exposure_vocabulary", c)
            mod = importlib.util.module_from_spec(spec)
            # Registered BEFORE exec: @dataclass resolves its own module through
            # sys.modules, and a module loaded by path is not there by default.
            sys.modules["exposure_vocabulary"] = mod
            spec.loader.exec_module(mod)
            for need in ("CLEAN_STATES", "CONFOUNDED_STATES", "disposition", "ExposureLedger"):
                if not hasattr(mod, need):
                    raise ImportError(f"{c}: exposure vocabulary is missing `{need}`")
            return mod, str(c)
    raise ImportError(
        "the exposure state vocabulary was not found. This tool imports it rather than "
        "restating it, so it cannot run without it. Set $EXPOSURE_VOCABULARY to an "
        "exposure.py, or place one beside this file. Tried: " + "; ".join(tried)
    )


# --------------------------------------------------------------------------------------
# Backends. numpy alone is sufficient; torch is used only when a device is available.
# --------------------------------------------------------------------------------------
class _NumpyBackend:
    name = "numpy-cpu"

    def array(self, a):
        return np.ascontiguousarray(a, dtype=np.float32)

    def znorm(self, A):
        """Row-centre and row-normalise, so a dot product IS the Pearson correlation."""
        A = A - A.mean(axis=1, keepdims=True)
        n = np.linalg.norm(A, axis=1, keepdims=True)
        n[n == 0] = 1.0
        return A / n

    def gram(self, Qz, Rz):
        return Qz @ Rz.T

    def full(self, shape, value, dtype="f"):
        return np.full(shape, value, dtype=np.float32 if dtype == "f" else np.int64)

    def merge_topk(self, bv, bi, C, start, k, mask=None):
        if mask is not None:
            C = np.where(mask, np.float32(-2.0), C)
        cols = np.arange(start, start + C.shape[1], dtype=np.int64)
        v = np.concatenate([bv, C], axis=1)
        i = np.concatenate([bi, np.broadcast_to(cols, C.shape)], axis=1)
        part = np.argpartition(-v, k - 1, axis=1)[:, :k]
        r = np.arange(v.shape[0])[:, None]
        vv, ii = v[r, part], i[r, part]
        o = np.argsort(-vv, axis=1, kind="stable")
        return vv[r, o], ii[r, o]

    def numpy(self, a):
        return np.asarray(a)


class _TorchBackend:
    def __init__(self, device):
        import torch

        self.t = torch
        self.dev = device
        self.name = f"torch-{device}"

    def array(self, a):
        return self.t.as_tensor(np.ascontiguousarray(a, dtype=np.float32), device=self.dev)

    def znorm(self, A):
        A = A - A.mean(dim=1, keepdim=True)
        n = self.t.linalg.norm(A, dim=1, keepdim=True)
        n[n == 0] = 1.0
        return A / n

    def gram(self, Qz, Rz):
        return Qz @ Rz.T

    def full(self, shape, value, dtype="f"):
        kw = dict(device=self.dev)
        if dtype == "f":
            return self.t.full(shape, float(value), dtype=self.t.float32, **kw)
        return self.t.full(shape, int(value), dtype=self.t.long, **kw)

    def merge_topk(self, bv, bi, C, start, k, mask=None):
        if mask is not None:
            C = C.masked_fill(mask, -2.0)
        cols = self.t.arange(start, start + C.shape[1], device=self.dev)
        v = self.t.cat([bv, C], dim=1)
        i = self.t.cat([bi, cols.expand(C.shape[0], -1)], dim=1)
        tv, ti = self.t.topk(v, k, dim=1)
        return tv, self.t.gather(i, 1, ti)

    def numpy(self, a):
        return a.detach().cpu().numpy()


def make_backend(device: str):
    """`auto` uses CUDA when torch reports one, otherwise numpy. Never required."""
    if device == "numpy":
        return _NumpyBackend()
    try:
        import torch
    except Exception:
        if device == "cuda":
            raise RuntimeError("--device cuda but torch is not importable")
        return _NumpyBackend()
    if device == "cuda" or (device == "auto" and torch.cuda.is_available()):
        if not torch.cuda.is_available():
            raise RuntimeError("--device cuda but torch reports no CUDA device")
        return _TorchBackend("cuda")
    if device == "torch-cpu":
        return _TorchBackend("cpu")
    return _NumpyBackend()


# --------------------------------------------------------------------------------------
# Structure statistics (LEG 1) and the bijection test.
# --------------------------------------------------------------------------------------
def block_stats(idx: np.ndarray, frac: float = 0.9) -> dict:
    """Order-FREE evidence of identity, which the offset run is not.

    The offset run needs the query row order to match the deposit's. It usually does not,
    so it understates. Two things do not depend on order:

      window       the narrowest window of corpus rows holding `frac` of the matches. A
                   study copied into the corpus occupies one contiguous run, so the width
                   is about the study's size. A merely similar study puts the matches
                   wherever its own neighbours sit, scattered.
      injectivity  the fraction of DISTINCT corpus rows matched. Nearest neighbours to a
                   foreign study collide. Matching n samples to n distinct rows is what a
                   duplicated study looks like.
    """
    s = np.sort(np.asarray(idx, dtype=np.int64))
    n = len(s)
    k = max(2, int(np.ceil(frac * n)))
    k = min(k, n)
    w = int((s[k - 1:] - s[:n - k + 1]).min()) + 1
    return {
        "narrowest_window_holding_90pct": w,
        "window_width_over_n": round(w / n, 4),
        "distinct_fraction": round(len(set(s.tolist())) / n, 4),
        "n_distinct_corpus_rows": int(len(set(s.tolist()))),
        "corpus_row_range": [int(s[0]), int(s[-1])],
        "span": int(s[-1] - s[0] + 1),
    }


def largest_constant_offset_run(idx: np.ndarray) -> tuple[int, int]:
    """Largest number of samples on one constant offset (corpus_row minus rank).

    The cohort's own row order is its rank. Returns (count, offset).
    """
    idx = np.asarray(idx, dtype=np.int64)
    off = idx - np.arange(len(idx), dtype=np.int64)
    u, c = np.unique(off, return_counts=True)
    j = int(np.argmax(c))
    return int(c[j]), int(u[j])


def permutation_null(idx: np.ndarray, n_perm: int, rng: np.random.Generator) -> np.ndarray:
    """Null for the aligned run: permute the QUERY ORDER.

    This destroys the correspondence while preserving every correlation and every matched
    corpus row, so a large aligned block cannot be manufactured by mere similarity to a
    dense corpus neighbourhood. It is the increment's own null, not a label-free MDE.
    """
    idx = np.asarray(idx, dtype=np.int64)
    n = len(idx)
    out = np.empty(n_perm, dtype=np.int64)
    for i in range(n_perm):
        out[i] = largest_constant_offset_run(idx[rng.permutation(n)])[0]
    return out


def analytic_expected_max_run(idx: np.ndarray, max_pairs: int = 40_000_000) -> float:
    """Expected largest constant-offset run under a random query order.

    A DIAGNOSTIC, not a bar. Under a permutation each of the n samples lands on one of the
    achievable offsets with probability p(o), which is exactly the (matched row, position)
    pair count normalised. The number of offsets carrying at least k samples is then
    approximately Poisson with mean E_k = sum_o C(n,k) p(o)^k (1-p(o))^(n-k), and

        E[M] = sum_{k>=1} P(M >= k) ~= sum_{k>=1} (1 - exp(-E_k)).

    Computed in log space so a large n does not overflow the binomial coefficient. Its
    only use is checking that the permutation machinery behaves; the reported bar is the
    exceedance count.
    """
    idx = np.asarray(idx, dtype=np.int64)
    n = len(idx)
    if n < 2:
        return float(n)
    if n * n <= max_pairs:
        offs = idx[:, None] - np.arange(n, dtype=np.int64)[None, :]
        _, cnt = np.unique(offs, return_counts=True)
        p = cnt.astype(np.float64) / cnt.sum()
    else:
        # Uniform fallback over the achievable offset range for a very large cohort.
        m = int(idx.max() - idx.min() + 1) + n - 1
        p = np.full(m, 1.0 / m)
    lg = math.lgamma
    lp, l1p = np.log(p), np.log1p(-np.minimum(p, 1 - 1e-15))
    total = 0.0
    for k in range(1, n + 1):
        logc = lg(n + 1) - lg(k + 1) - lg(n - k + 1)
        ek = float(np.sum(np.exp(logc + k * lp + (n - k) * l1p)))
        total += 1.0 - math.exp(-min(ek, 700.0))
        if k > 1 and ek < 1e-12:
            break
    return total


# --------------------------------------------------------------------------------------
# The scan itself: running argmax over corpus chunks, never the whole corpus.
# --------------------------------------------------------------------------------------
@dataclass
class QueryBlock:
    name: str
    Z: object                        # backend array, row-normalised on the shared genes
    self_rows: np.ndarray | None     # corpus row to exclude per query row, or None


def scan_corpus(reader, blocks: Sequence[QueryBlock], corpus_cols: np.ndarray,
                backend, chunk_rows: int = 8192, topk: int = 5,
                progress_every: int = 20, quiet: bool = False) -> dict:
    """Stream the corpus once, keeping only a running top-k per query row.

    Peak memory is one chunk of the corpus plus the query, whatever the corpus size.
    """
    n_ref = reader.n_rows
    best_v = {b.name: backend.full((b.Z.shape[0], topk), -2.0, "f") for b in blocks}
    best_i = {b.name: backend.full((b.Z.shape[0], topk), -1, "i") for b in blocks}
    self_rows = {}
    for b in blocks:
        if b.self_rows is not None:
            if isinstance(backend, _TorchBackend):
                self_rows[b.name] = backend.t.as_tensor(
                    b.self_rows.astype(np.int64), device=backend.dev)[:, None]
            else:
                self_rows[b.name] = b.self_rows.astype(np.int64)[:, None]

    t0 = time.time()
    n_chunks = 0
    for s, e, R in reader.iter_chunks(chunk_rows):
        Rz = backend.znorm(backend.array(R[:, corpus_cols]))
        for b in blocks:
            C = backend.gram(b.Z, Rz)
            mask = None
            if b.name in self_rows:
                if isinstance(backend, _TorchBackend):
                    rows = backend.t.arange(s, e, device=backend.dev)[None, :]
                else:
                    rows = np.arange(s, e, dtype=np.int64)[None, :]
                mask = self_rows[b.name] == rows
            best_v[b.name], best_i[b.name] = backend.merge_topk(
                best_v[b.name], best_i[b.name], C, s, topk, mask)
        n_chunks += 1
        if not quiet and n_chunks % progress_every == 0:
            el = time.time() - t0
            print(f"  {e}/{n_ref} corpus rows  {el / 60:5.1f} min  "
                  f"eta {el / max(e, 1) * (n_ref - e) / 60:5.1f} min", flush=True)
    return {
        "max": {k: backend.numpy(v).astype(np.float64) for k, v in best_v.items()},
        "idx": {k: backend.numpy(v).astype(np.int64) for k, v in best_i.items()},
        "n_chunks": n_chunks,
        "seconds": round(time.time() - t0, 2),
    }


# --------------------------------------------------------------------------------------
# LEG 2: the same-sample band, calibrated from proven bijections.
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Criterion:
    """Every number a call depends on, in one object that is written to the output."""
    injectivity_min: float = 0.85
    window_over_n_max: float = 3.0
    band_quantile: float = 0.05
    proven_aligned_fraction_min: float = 0.90
    proven_max_exceedances: int = 0
    n_permutations: int = 2000
    topk: int = 5


def calibrate_band(cohorts: dict, matches: dict, crit: Criterion) -> dict | None:
    """The same-sample cross-pipeline band, MEASURED from cohorts the bijection proved.

    Those samples are the same libraries the query pipeline quantified, so their
    correlation IS the attenuation between the two pipelines. Nothing here is assumed.

    Returns None when no cohort is proven -- the caller must then refuse to apply leg 2
    rather than substitute a default. A cohort that calibrates the band cannot be tested
    by it; its evidence is leg 1 plus the permutation test, and `calibrates_band` records
    which cohorts those are.
    """
    proven = [
        c for c, v in cohorts.items()
        if v["n_exceedances"] <= crit.proven_max_exceedances
        and v["aligned_fraction"] >= crit.proven_aligned_fraction_min
        and v["distinct_fraction"] >= crit.injectivity_min
    ]
    if not proven:
        return None
    parts = []
    for c in proven:
        idx, mx = matches[c]
        on = (idx - np.arange(len(idx))) == cohorts[c]["constant_offset"]
        parts.append(mx[on])
    band = np.concatenate(parts)
    return {
        "from_cohorts": sorted(proven),
        "n_samples": int(len(band)),
        "min": float(band.min()),
        f"p{int(crit.band_quantile * 100):02d}": float(np.quantile(band, crit.band_quantile)),
        "median": float(np.median(band)),
        "max": float(band.max()),
        "lower_bound_applied": float(np.quantile(band, crit.band_quantile)),
    }


# --------------------------------------------------------------------------------------
# Block confirmation: mutual argmax and diagonal vs off-diagonal, on the matched window.
# --------------------------------------------------------------------------------------
def block_confirm(reader, corpus_cols: np.ndarray, Qz_np: np.ndarray,
                  matched_rows: np.ndarray) -> dict:
    """Pull the matched corpus rows and score the cohort against them, all pairs.

    The correspondence used is the ASSIGNMENT (each query row's own argmax), not the
    identity permutation, so this is meaningful for a cohort whose deposit order differs
    from the query's. `mutual_argmax` counts query rows i whose best block row is b AND
    for which b's best query row is i -- an argmax over both rows and columns.
    """
    rows = np.unique(matched_rows)
    B = reader.take_rows(rows.tolist())[:, corpus_cols].astype(np.float64)
    B = B - B.mean(axis=1, keepdims=True)
    nb = np.linalg.norm(B, axis=1, keepdims=True)
    nb[nb == 0] = 1.0
    B /= nb
    C = Qz_np.astype(np.float64) @ B.T                      # (n_query, n_block)
    pos = {int(r): k for k, r in enumerate(rows)}
    assign = np.array([pos[int(r)] for r in matched_rows], dtype=np.int64)
    n = C.shape[0]
    diag = C[np.arange(n), assign]
    off = C.copy()
    off[np.arange(n), assign] = -np.inf
    row_ok = int((C.argmax(axis=1) == assign).sum())
    col_best = C.argmax(axis=0)
    mutual = int(sum(1 for i in range(n)
                     if C[i].argmax() == assign[i] and col_best[assign[i]] == i))
    finite_off = off[np.isfinite(off)]
    return {
        "n_query": n,
        "n_block_rows": int(len(rows)),
        "diag_mean": float(diag.mean()),
        "diag_min": float(diag.min()),
        "offdiag_mean": float(finite_off.mean()) if finite_off.size else float("nan"),
        "offdiag_max": float(finite_off.max()) if finite_off.size else float("nan"),
        "row_argmax_on_assignment": row_ok,
        "mutual_argmax": mutual,
        "min_margin_over_best_other": float((diag - off.max(axis=1)).min()),
        "median_margin_over_best_other": float(np.median(diag - off.max(axis=1))),
    }


# --------------------------------------------------------------------------------------
# The audit.
# --------------------------------------------------------------------------------------
def audit(query: np.ndarray, query_genes: Sequence[str], cohorts: Sequence[str],
          corpus_path, *, corpus_dataset: str = "X", corpus_gene_key: str | None = None,
          corpus_genes: Sequence[str] | None = None, assume_aligned: bool = False,
          crit: Criterion = Criterion(), seed: int = DEFAULT_SEED,
          chunk_rows: int = 8192, n_background: int = 1500, device: str = "auto",
          band_lower_override: float | None = None, declared: dict | None = None,
          model_id: str = "unnamed_model", do_block_confirm: bool = True,
          quiet: bool = False) -> dict:
    """Full per-cohort exposure audit. `query` is (n_samples, n_genes)."""
    vocab, vocab_path = load_exposure_vocabulary()
    reader = open_corpus(corpus_path, dataset=corpus_dataset, gene_key=corpus_gene_key,
                         gene_ids=list(corpus_genes) if corpus_genes is not None else None)
    try:
        query = np.asarray(query, dtype=np.float32)
        cohorts = np.asarray([str(c) for c in cohorts])
        if query.shape[0] != len(cohorts):
            raise ValueError(f"{query.shape[0]} query samples but {len(cohorts)} cohort labels")

        # ---- shared gene vocabulary -------------------------------------------------
        if assume_aligned:
            if query.shape[1] != reader.n_cols:
                raise ValueError(
                    f"--assume-aligned but query has {query.shape[1]} genes and the corpus "
                    f"{reader.n_cols}")
            qcols = np.arange(query.shape[1], dtype=np.int64)
            ccols = qcols.copy()
            shared_note = "positional, --assume-aligned"
        else:
            if reader.gene_ids is None:
                raise ValueError(
                    "the corpus carries no gene identifiers. Pass --corpus-genes FILE, or "
                    "--assume-aligned if the two matrices are known to share a column order.")
            cpos = {}
            for i, g in enumerate(reader.gene_ids):
                cpos.setdefault(g, i)
            shared = [g for g in query_genes if g in cpos]
            if len(shared) < 50:
                raise ValueError(
                    f"only {len(shared)} genes shared between the query and the corpus. "
                    f"The two matrices must share a vocabulary.")
            qpos = {g: i for i, g in enumerate(query_genes)}
            qcols = np.array([qpos[g] for g in shared], dtype=np.int64)
            ccols = np.array([cpos[g] for g in shared], dtype=np.int64)
            shared_note = "intersection of gene identifiers"

        Q = query[:, qcols]
        finite = np.isfinite(Q).all(axis=0)
        n_dropped = int((~finite).sum())
        if n_dropped:
            qcols, ccols, Q = qcols[finite], ccols[finite], Q[:, finite]
        if Q.shape[1] < 50:
            raise ValueError(f"only {Q.shape[1]} usable shared genes after dropping non-finite")

        backend = make_backend(device)
        if not quiet:
            print(f"  backend {backend.name}; {Q.shape[0]} query samples over "
                  f"{Q.shape[1]} shared genes; corpus {reader.n_rows} rows", flush=True)

        # ---- background probes: corpus rows scored against the corpus, self excluded --
        rng = np.random.default_rng(seed)
        n_bg = int(min(n_background, reader.n_rows))
        bg_rows = np.sort(rng.choice(reader.n_rows, size=n_bg, replace=False))
        bg = reader.take_rows(bg_rows.tolist())[:, ccols]

        blocks = [
            QueryBlock("__query__", backend.znorm(backend.array(Q)), None),
            QueryBlock("__background__", backend.znorm(backend.array(bg)), bg_rows),
        ]
        scan = scan_corpus(reader, blocks, ccols, backend, chunk_rows=chunk_rows,
                           topk=crit.topk, quiet=quiet)
        qmax = scan["max"]["__query__"][:, 0]
        qidx = scan["idx"]["__query__"][:, 0]
        qgap = scan["max"]["__query__"][:, 0] - scan["max"]["__query__"][:, 1]
        bgmax = scan["max"]["__background__"][:, 0]

        # ---- per cohort --------------------------------------------------------------
        rng_perm = np.random.default_rng(seed)
        per, matches = {}, {}
        for c in sorted(set(cohorts.tolist())):
            m = cohorts == c
            idx, mx, gp = qidx[m], qmax[m], qgap[m]
            n = int(m.sum())
            run, off = largest_constant_offset_run(idx)
            null = permutation_null(idx, crit.n_permutations, rng_perm)
            rec = {
                "n": n,
                "largest_aligned_run": run,
                "aligned_fraction": round(run / n, 4),
                "constant_offset": off,
                "n_permutations": int(crit.n_permutations),
                "n_exceedances": int((null >= run).sum()),
                "exceedance_rate": float((null >= run).mean()),
                "perm_null_mean_run": float(null.mean()),
                "analytic_expected_max_run": round(analytic_expected_max_run(idx), 4),
                "max_correlation_median": float(np.median(mx)),
                "max_correlation_min": float(mx.min()),
                "max_correlation_max": float(mx.max()),
                "gap_max_minus_second_median": float(np.median(gp)),
            }
            rec.update(block_stats(idx))
            per[c] = rec
            matches[c] = (idx, mx)

        # ---- LEG 2: the band ---------------------------------------------------------
        band = calibrate_band(per, matches, crit)
        if band_lower_override is not None:
            band_mode = "operator_supplied"
            band_lower = float(band_lower_override)
            band = (band or {}) | {"lower_bound_applied": band_lower,
                                   "operator_supplied": True}
        elif band is not None:
            band_mode = "from_proven_bijections"
            band_lower = band["lower_bound_applied"]
        else:
            band_mode = "not_calibrable"
            band_lower = None

        # ---- the two-leg call --------------------------------------------------------
        calibrators = set(band["from_cohorts"]) if (band and "from_cohorts" in band) else set()
        state, why = {}, {}
        for c, v in per.items():
            leg1 = (v["distinct_fraction"] >= crit.injectivity_min
                    and v["window_width_over_n"] <= crit.window_over_n_max)
            leg2 = band_lower is not None and v["max_correlation_median"] >= band_lower
            v["leg1_structure_passes"] = bool(leg1)
            v["leg2_magnitude_passes"] = bool(leg2)
            v["leg2_decidable"] = band_lower is not None
            v["calibrates_band"] = c in calibrators
            # `encoder_seen` needs BOTH legs. Anything short of that is `unknown`: a corpus
            # with no accessions cannot yield `clean_declared`, and absence of evidence is
            # not evidence of absence. The vocabulary already has a word for that.
            state[c] = "encoder_seen" if (leg1 and leg2) else "unknown"
            why[c] = (
                f"injectivity {v['distinct_fraction']}, 90% of matches inside a window "
                f"{v['window_width_over_n']}x the cohort size, median r "
                f"{v['max_correlation_median']:.4f}"
                + (f", band lower bound {band_lower:.4f}" if band_lower is not None
                   else ", band NOT CALIBRABLE so leg 2 could not be applied")
                + f", aligned run {v['largest_aligned_run']}/{v['n']} with "
                  f"{v['n_exceedances']}/{v['n_permutations']} exceedances")

        # ---- block confirmation for the called cohorts -------------------------------
        confirm = {}
        if do_block_confirm:
            Qz_np = _NumpyBackend().znorm(np.asarray(Q, dtype=np.float64))
            for c in sorted(state):
                if state[c] != "encoder_seen":
                    continue
                m = cohorts == c
                confirm[c] = block_confirm(reader, ccols, Qz_np[m], qidx[m])

        # ---- ledger, dispositions, declared-vs-measured ------------------------------
        ledger = vocab.ExposureLedger.from_per_study(model_id, state)
        disp = {c: vocab.disposition(s) for c, s in state.items()}
        dvm = None
        if declared:
            dvm = {}
            for c, d in declared.items():
                if c not in state:
                    continue
                dvm[c] = {
                    "declared": d, "measured": state[c],
                    "declared_disposition": vocab.disposition(d),
                    "measured_disposition": disp[c],
                    "agree": vocab.disposition(d) == disp[c],
                }

        out = {
            "tool": TOOL,
            "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "model_id": model_id,
            "exposure_vocabulary_module": vocab_path,
            "seed": int(seed),
            "backend": backend.name,
            "corpus": reader.describe() | {"chunks_read": scan["n_chunks"],
                                           "scan_seconds": scan["seconds"]},
            "genes": {
                "query": int(query.shape[1]), "corpus": int(reader.n_cols),
                "shared_used": int(Q.shape[1]), "dropped_nonfinite": n_dropped,
                "alignment": shared_note,
            },
            "criterion": asdict(crit),
            "corpus_background": {
                "n_probes": n_bg,
                "note": ("corpus rows scored against the corpus with self excluded: the "
                         "same-pipeline, different-sample distribution, which is the right "
                         "reference for how high a maximum gets with no duplicate"),
                "median_nn_correlation": float(np.median(bgmax)),
                "p90_nn_correlation": float(np.quantile(bgmax, 0.90)),
                "p99_nn_correlation": float(np.quantile(bgmax, 0.99)),
                "self_duplication_rate": float((bgmax > 0.99999).mean()),
            },
            "same_sample_band": band,
            "band_calibration": band_mode,
            "cohorts": per,
            "block_confirmation": confirm,
            "exposure_state_measured": state,
            "exposure_reason": why,
            "disposition": disp,
            "ledger": {
                "model_id": ledger.model_id, "clean": list(ledger.clean),
                "confounded": list(ledger.confounded), "unresolved": list(ledger.unresolved),
                "has_clean_study": bool(ledger.has_clean_study),
            },
            "declared_vs_measured": dvm,
            "negative_reading": NEGATIVE_READING,
        }
        out["_matches"] = {c: (matches[c][0], matches[c][1]) for c in matches}
        return out
    finally:
        reader.close()


# --------------------------------------------------------------------------------------
# I/O helpers and CLI.
# --------------------------------------------------------------------------------------
def read_matrix(path, orient: str = "genes") -> tuple[np.ndarray, list[str], list[str]]:
    """Read a query matrix. `genes`: rows are genes, columns are samples (the default,
    matching how expression matrices are usually deposited). `samples`: the transpose."""
    p = Path(path)
    if p.suffix == ".npz":
        z = np.load(str(p), allow_pickle=False)
        X = z["X"]
        genes = [str(g) for g in z["genes"]]
        samples = [str(s) for s in z["samples"]] if "samples" in z else \
            [str(i) for i in range(X.shape[0])]
        return np.asarray(X, dtype=np.float32), genes, samples
    opener = __import__("gzip").open if str(p).endswith(".gz") else open
    with opener(str(p), "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        names, rows, vals = header[1:], [], []
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            rows.append(parts[0])
            vals.append(parts[1:])
    M = np.asarray(vals, dtype=np.float32)
    if orient == "genes":
        return M.T, rows, names            # -> (samples, genes)
    return M, names, rows


def read_cohorts(path, samples: Sequence[str]) -> list[str]:
    """Two-column sample_id -> cohort table, with or without a header."""
    pairs = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            if not line.strip():
                continue
            a, _, b = line.rstrip("\n").partition("\t")
            if i == 0 and b.strip().lower() in {"cohort", "study", "group", "label"}:
                continue
            pairs[a.strip()] = b.strip()
    missing = [s for s in samples if s not in pairs]
    if missing:
        raise KeyError(f"{len(missing)} query samples have no cohort label, e.g. {missing[:5]}")
    return [pairs[s] for s in samples]


class _Enc(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return super().default(o)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="exposure_scan.py",
        description="Resolve pretraining exposure in expression space, per cohort.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="A negative result reads 'no near-duplicate detected', never "
               "'not in the corpus'. See MODEL_CARD.md.")
    p.add_argument("--query", required=True,
                   help="query matrix: TSV (rows=genes, cols=samples) or .npz with X/genes/samples")
    p.add_argument("--query-orient", default="genes", choices=("genes", "samples"))
    p.add_argument("--cohorts", required=True, help="TSV: sample_id <TAB> cohort")
    p.add_argument("--corpus", required=True,
                   help=".h5ad/.h5, .npy, .npz or .tsv(.gz) of corpus rows x genes")
    p.add_argument("--corpus-dataset", default="X")
    p.add_argument("--corpus-gene-key", default=None,
                   help="key under var/ carrying gene identifiers (h5ad)")
    p.add_argument("--corpus-genes", default=None,
                   help="text file, one gene id per corpus column")
    p.add_argument("--assume-aligned", action="store_true",
                   help="the two matrices already share a column order")
    p.add_argument("--out", required=True)
    p.add_argument("--emit-matches", default=None, help="npz of per-sample matched rows")
    p.add_argument("--model-id", default="unnamed_model")
    p.add_argument("--declared", default=None,
                   help="JSON mapping cohort -> declared exposure state, for reconciliation")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument("--chunk-rows", type=int, default=8192)
    p.add_argument("--topk", type=int, default=5)
    p.add_argument("--n-background", type=int, default=1500)
    p.add_argument("--n-permutations", type=int, default=2000)
    p.add_argument("--injectivity-min", type=float, default=0.85)
    p.add_argument("--window-over-n-max", type=float, default=3.0)
    p.add_argument("--band-quantile", type=float, default=0.05)
    p.add_argument("--band-lower", type=float, default=None,
                   help="override the calibrated band lower bound. Use only with a band "
                        "measured on proven same-sample pairs for THIS pipeline pair.")
    p.add_argument("--device", default="auto",
                   choices=("auto", "cuda", "torch-cpu", "numpy"))
    p.add_argument("--no-block-confirm", action="store_true")
    p.add_argument("--quiet", action="store_true")
    return p


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    X, genes, samples = read_matrix(a.query, a.query_orient)
    cohorts = read_cohorts(a.cohorts, samples)
    corpus_genes = None
    if a.corpus_genes:
        corpus_genes = [l.strip() for l in open(a.corpus_genes) if l.strip()]
    declared = json.load(open(a.declared)) if a.declared else None
    crit = Criterion(injectivity_min=a.injectivity_min,
                     window_over_n_max=a.window_over_n_max,
                     band_quantile=a.band_quantile,
                     n_permutations=a.n_permutations,
                     topk=a.topk)
    res = audit(X, genes, cohorts, a.corpus, corpus_dataset=a.corpus_dataset,
                corpus_gene_key=a.corpus_gene_key, corpus_genes=corpus_genes,
                assume_aligned=a.assume_aligned, crit=crit, seed=a.seed,
                chunk_rows=a.chunk_rows, n_background=a.n_background, device=a.device,
                band_lower_override=a.band_lower, declared=declared,
                model_id=a.model_id, do_block_confirm=not a.no_block_confirm,
                quiet=a.quiet)
    m = res.pop("_matches")
    Path(a.out).write_text(json.dumps(res, indent=2, cls=_Enc) + "\n")
    if a.emit_matches:
        np.savez_compressed(a.emit_matches,
                            sample_id=np.asarray(samples), cohort=np.asarray(cohorts),
                            **{f"{c}__corpus_row": m[c][0] for c in m},
                            **{f"{c}__max_r": m[c][1] for c in m})
    if not a.quiet:
        for c in sorted(res["exposure_state_measured"]):
            print(f"  {c:<28} {res['exposure_state_measured'][c]:<14} "
                  f"{res['disposition'][c]:<11} {res['exposure_reason'][c]}")
        print(f"  band calibration: {res['band_calibration']}")
        print(f"  wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
