#!/usr/bin/env python3
"""Donor-held-out cis-accessibility comparison with explicit paired-RNA endpoints.

Aggregates are seeded, disjoint 50-cell partitions within donors. No imputed
expression or embedding is used. Every link fit is confined to training donors.
Dense allocations are limited to a gene's cis window.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import multiprocessing
import os
from pathlib import Path
import socket
import subprocess
import sys
import time
from types import MappingProxyType, ModuleType

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite
# Load SciPy's lazily imported MMIO library before threadpoolctl scans loaded
# libraries; otherwise the very first MatrixMarket call can escape its limit.
try:
    from scipy.io._fast_matrix_market import _fmm_core  # noqa: F401
except ImportError:
    _fmm_core = None  # Older SciPy uses the serial MatrixMarket implementation.
from threadpoolctl import threadpool_limits


def stable_seed(seed, value):
    def canonical(item):
        if isinstance(item, np.generic):
            return item.item()
        if isinstance(item, tuple):
            return tuple(canonical(x) for x in item)
        if isinstance(item, list):
            return [canonical(x) for x in item]
        return item
    value = canonical(value)
    return (int(seed) + int(hashlib.sha256(str(value).encode()).hexdigest()[:8], 16)) % (2**32)


def require_compute():
    if not os.environ.get("SLURM_JOB_ID") or "login" in socket.gethostname():
        raise RuntimeError("Network fitting requires a compute allocation")


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def read_matrix(root, names):
    for name in names:
        path = root / name
        if path.exists():
            with threadpool_limits(limits=1):
                return (sparse.load_npz(path) if path.suffix == ".npz" else sparse.csr_matrix(mmread(path))).tocsr()
    diagnostics = root / "diagnostics.json"
    if diagnostics.exists():
        chunks = json.loads(diagnostics.read_text()).get("chunks", [])
        if chunks:
            raise MemoryError("Chunked atlas input: select bounded donor/cell rows before network execution")
    raise FileNotFoundError(f"No matrix in {root}: {names}")


def cell_keys(cells):
    fields = ["cohort", "donor_id", "cell_id"]
    if not set(fields).issubset(cells):
        raise ValueError(f"cells.tsv requires {fields}")
    keys = cells[fields].agg("|".join, axis=1)
    if keys.duplicated().any():
        raise ValueError("Duplicated cohort/donor/cell identity")
    return keys


def condition_eligibility(cells, settings, technical_pilot=False):
    """Apply the prespecified disease population before any cell aggregation."""
    keep = np.ones(len(cells), dtype=bool)
    if not technical_pilot:
        for cohort, allowed in settings.get("allowed_conditions_by_cohort", {}).items():
            cohort_rows = cells.cohort.eq(cohort).to_numpy()
            if not cohort_rows.any():
                continue
            if "condition" not in cells:
                raise ValueError(f"Condition is required for prespecified cohort roster: {cohort}")
            allowed = {str(value).strip().upper() for value in allowed}
            keep[cohort_rows] = cells.loc[cohort_rows, "condition"].str.strip().str.upper().isin(allowed)
    fields = ["cohort", "lineage", "donor_id"] + (["condition"] if "condition" in cells else [])
    excluded = cells.loc[~keep].groupby(fields, dropna=False).size().reset_index(name="n_excluded_cells")
    excluded["reason"] = "outside_prespecified_disease_population"
    return keep, excluded


def load_input(root, seed=42, cap=1000):
    cells, peaks = read_table(root / "cells.tsv"), read_table(root / "peaks.tsv")
    matrix_names = ["cell_counts.npz", "matrix.npz", "ATAC_counts.mtx.gz", "ATAC_counts.mtx"]
    if any((root / name).exists() for name in matrix_names):
        matrix = read_matrix(root, matrix_names)
    else:
        # Load at most cap selected rows per donor directly from sparse chunks.
        description = json.loads((root / "diagnostics.json").read_text())
        chosen = []
        for key, group in cells.groupby(["cohort", "lineage", "donor_id"], sort=True).indices.items():
            chosen.extend(np.random.default_rng(stable_seed(seed, key)).permutation(group)[:cap])
        chosen = np.array(sorted(chosen), dtype=int)
        pieces = []
        seen = []
        for chunk in sorted(description["chunks"], key=lambda x: int(x["row_start"])):
            start, n = int(chunk["row_start"]), int(chunk["rows"])
            rows = chosen[(chosen >= start) & (chosen < start + n)]
            if len(rows):
                pieces.append(sparse.load_npz(root / chunk["path"])[rows - start])
                seen.extend(rows.tolist())
        if seen != chosen.tolist():
            raise ValueError("Chunk metadata does not cover selected network cells")
        matrix = sparse.vstack(pieces, format="csr")
        cells = cells.iloc[chosen].reset_index(drop=True)
        cells["original_source_index"] = chosen
    cell_keys(cells)
    if matrix.shape != (len(cells), len(peaks)):
        raise ValueError("ATAC dimensions disagree with metadata")
    if "lineage" not in cells:
        raise ValueError("cells.tsv requires lineage")
    if "peak_coordinate" in peaks and "peak_id" not in peaks:
        peaks["peak_id"] = peaks.peak_coordinate
    if not {"chrom", "start0", "end"}.issubset(peaks):
        parsed = peaks.peak_id.str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
        if parsed.isna().any().any():
            raise ValueError("Malformed peak coordinates")
        peaks[["chrom", "start0", "end"]] = parsed.to_numpy()
    peaks[["start0", "end"]] = peaks[["start0", "end"]].astype(np.int64)
    canonical = peaks.chrom + ":" + peaks.start0.astype(str) + "-" + peaks.end.astype(str)
    if not peaks.peak_id.equals(canonical):
        declaration = root / "paired_export.json"
        system = json.loads(declaration.read_text()).get("peak_id_coordinate_system") if declaration.exists() else None
        native = peaks.peak_id.str.extract(r"^(chr[^-]+)-(\d+)-(\d+)$")
        if system != "GRanges_1based_closed" or native.isna().any().any():
            raise ValueError("peak_id disagrees with canonical coordinates without an explicit native-ID mapping")
        if not (native[0].equals(peaks.chrom) and
                (native[1].astype(np.int64) - 1).equals(peaks.start0) and
                native[2].astype(np.int64).equals(peaks.end)):
            raise ValueError("Declared native peak_id mapping disagrees with canonical coordinates")
        if canonical.duplicated().any():
            raise ValueError("Native peak identity mapping is not one-to-one")
    peaks["canonical_peak_id"] = canonical
    if peaks.peak_id.duplicated().any() or (peaks.start0 < 0).any() or (peaks.end <= peaks.start0).any():
        raise ValueError("Invalid or duplicate peak coordinates")
    if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise ValueError("ATAC counts must be finite and nonnegative")
    return matrix, cells, peaks


def build_metacells(matrix, cells, seed=42, size=50, cap=1000):
    """Donor-confined groups with source indices, equal downstream donor weights."""
    groups, meta, eligibility = [], [], []
    for key, values in cells.groupby(["cohort", "lineage", "donor_id"], sort=True).indices.items():
        idx = np.asarray(values, dtype=int)
        reason = "" if len(idx) >= 2 * size else "fewer_than_two_50_cell_aggregates"
        eligibility.append(dict(zip(["cohort", "lineage", "donor_id"], key), n_cells=len(idx), eligible=not reason, reason=reason))
        if reason:
            continue
        chosen = np.random.default_rng(stable_seed(seed, key)).permutation(idx)[: min(len(idx), cap)]
        for start in range(0, len(chosen) - size + 1, size):
            groups.append(chosen[start:start + size])
            row_meta = dict(zip(["cohort", "lineage", "donor_id"], key), metacell_id=f"{'_'.join(key)}_{start // size:03d}", n_cells=size)
            for field in ["lineage_state", "source_lineage", "source_label", "well_id"]:
                if field in cells:
                    row_meta[field] = ";".join(sorted(map(str, cells.iloc[groups[-1]][field].unique())))
            meta.append(row_meta)
    if not groups:
        return sparse.csr_matrix((0, matrix.shape[1])), pd.DataFrame(meta), groups, pd.DataFrame(eligibility)
    row = np.repeat(np.arange(len(groups)), size)
    col = np.concatenate(groups)
    assignment = sparse.csr_matrix((np.ones(len(col)), (row, col)), shape=(len(groups), len(cells)))
    aggregated = (assignment @ matrix).tocsr()
    if int(aggregated.sum()) != int(matrix[np.concatenate(groups)].sum()):
        raise ValueError("Metacell aggregation failed count conservation")
    return aggregated, pd.DataFrame(meta), groups, pd.DataFrame(eligibility)


def donor_folds(donors, n_splits, seed):
    values = np.array(sorted(set(map(str, donors))))
    if len(values) < n_splits:
        raise ValueError(f"Need at least {n_splits} biological donors")
    values = np.random.default_rng(seed).permutation(values)
    mapping = {donor: i % n_splits for i, donor in enumerate(values)}
    return np.array([mapping[str(d)] for d in donors], dtype=int)


def donor_weights(donors):
    names, counts = np.unique(np.asarray(donors, dtype=str), return_counts=True)
    lookup = dict(zip(names, counts))
    weights = np.array([1.0 / lookup[str(x)] for x in donors])
    return weights / weights.mean()


def normalize_sparse(matrix):
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    return (sparse.diags(1e4 / np.maximum(totals, 1)) @ matrix).tocsr()


def log_columns(matrix, indices):
    # CSR historically produced C-contiguous designs. Preserve that layout
    # after column caching so BLAS keeps identical floating-point arithmetic.
    return np.log1p(matrix[:, indices].toarray(order="C"))


def gene_peak_indices(peaks, gene, window=250000, promoter=2000):
    """TSS is strand-correct, 0-based; intervals are half-open."""
    tss = int(gene["tss0"])
    if gene["strand"] not in {"+", "-"}:
        raise ValueError("Gene strand must be + or -")
    same = np.asarray(peaks.chrom == gene["chrom"])
    prom = same & (np.asarray(peaks.start0) < tss + promoter + 1) & (np.asarray(peaks.end) > max(0, tss - promoter))
    cis = same & (np.asarray(peaks.start0) < tss + window + 1) & (np.asarray(peaks.end) > max(0, tss - window))
    return np.flatnonzero(prom), np.flatnonzero(cis & ~prom)


class PeakGeometry:
    """Exact interval queries bounded by chromosome and observed maximum width."""
    def __init__(self, peaks):
        self.chromosomes = {}
        for chrom, frame in peaks.groupby("chrom", sort=False):
            order = np.argsort(frame.start0.to_numpy(), kind="stable")
            starts, ends = frame.start0.to_numpy()[order], frame.end.to_numpy()[order]
            self.chromosomes[chrom] = (frame.index.to_numpy()[order], starts, ends, int(np.max(ends - starts)))

    def query(self, gene, window=250000, promoter=2000):
        if gene["strand"] not in {"+", "-"}:
            raise ValueError("Gene strand must be + or -")
        if gene["chrom"] not in self.chromosomes:
            return np.array([], dtype=int), np.array([], dtype=int)
        indices, starts, ends, width = self.chromosomes[gene["chrom"]]
        tss = int(gene["tss0"])
        low, high = max(0, tss - window), tss + window + 1
        left, right = np.searchsorted(starts, [low - width, high], side="left")
        selected = np.arange(left, right)
        selected = selected[ends[selected] > low]
        prom = (starts[selected] < tss + promoter + 1) & (ends[selected] > max(0, tss - promoter))
        return np.sort(indices[selected[prom]]), np.sort(indices[selected[~prom]])


def weighted_ridge(x, y, donors, alpha, predict):
    w = donor_weights(donors)
    mean = np.average(x, axis=0, weights=w) if x.shape[1] else np.empty(0)
    scale = np.sqrt(np.average((x - mean) ** 2, axis=0, weights=w)) if x.shape[1] else np.empty(0)
    scale = np.maximum(scale, 1e-8)
    z, test = (x - mean) / scale, (predict - mean) / scale
    ym = np.average(y, weights=w)
    if not x.shape[1]:
        return np.repeat(ym, len(predict))
    coef = np.linalg.solve(z.T @ (w[:, None] * z) + alpha * np.eye(z.shape[1]), z.T @ (w * (y - ym)))
    return test @ coef + ym


def holm(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p, kind="stable")
    out = np.empty(len(p))
    out[order] = np.minimum(1, np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    return out


def paired_inference(difference, seed=42, bootstrap_reps=2000, null_reps=10000):
    """Bootstrap donor CI and separate sign-flip null (symmetry assumption)."""
    d = np.asarray(difference, dtype=float)
    d = d[np.isfinite(d)]
    if len(d) < 5:
        return dict(n_donors=len(d), improvement=np.nan, ci_low=np.nan, ci_high=np.nan, pvalue=1.0, status="untestable")
    rng = np.random.default_rng(seed)
    boot = np.array([rng.choice(d, len(d), replace=True).mean() for _ in range(bootstrap_reps)])
    observed = d.mean()
    if len(d) <= 15:
        null = np.array([np.mean(d * (2 * ((code >> np.arange(len(d))) & 1) - 1)) for code in range(2 ** len(d))])
        pvalue = float(np.mean(null >= observed - 1e-15))
    else:
        exceed = sum(float(np.mean(d * rng.choice([-1, 1], size=len(d)))) >= observed - 1e-15 for _ in range(null_reps))
        pvalue = (exceed + 1) / (null_reps + 1)
    lo, hi = np.quantile(boot, [0.025, 0.975])
    return dict(n_donors=len(d), improvement=float(observed), ci_low=float(lo), ci_high=float(hi), pvalue=pvalue, status="evaluated")


def inventory():
    result = {"python": sys.version.split()[0]}
    for package in ["numpy", "scipy", "pandas", "xgboost", "threadpoolctl"]:
        try:
            result[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            result[package] = None
    return result


class CiceroGraph:
    """Sparse directed edge storage; query both orientations without zero-padding negatives."""
    def __init__(self, matrix):
        self.csr = matrix.tocsr()
        self.csc = matrix.tocsc()

    @classmethod
    def from_tsv(cls, path, peaks):
        mapping = dict(zip(peaks.peak_id, peaks.index))
        rows, columns, values = [], [], []
        finite_rows = 0
        for frame in pd.read_csv(path, sep="\t", chunksize=250000, dtype={"peak1": str, "peak2": str, "coaccess": float}):
            frame = frame[np.isfinite(frame.coaccess)]
            first, second = frame.peak1.map(mapping), frame.peak2.map(mapping)
            if first.isna().any() or second.isna().any():
                raise ValueError("Cicero output contains peaks outside its training universe")
            rows.append(first.to_numpy(dtype=np.int32))
            columns.append(second.to_numpy(dtype=np.int32))
            values.append(frame.coaccess.to_numpy(dtype=float))
            finite_rows += len(frame)
        size = int(max(peaks.index)) + 1 if len(peaks) else 0
        if finite_rows:
            matrix = sparse.coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))), shape=(size, size)).tocsr()
            if matrix.nnz != finite_rows:
                raise ValueError("Duplicate directed Cicero edges must be adjudicated before sparse scoring")
        else:
            matrix = sparse.csr_matrix((size, size))
        result = cls(matrix)
        metrics = {"finite_edges": finite_rows, "shape": [size, size],
                   "csr_csc_bytes": sum(x.data.nbytes + x.indices.nbytes + x.indptr.nbytes for x in [result.csr, result.csc])}
        Path(str(path) + ".sparse_graph.json").write_text(json.dumps(metrics, indent=2) + "\n")
        return result

    def promoter_scores(self, promoter_indices):
        scores = {}
        for index in promoter_indices:
            if index >= self.csr.shape[0]:
                continue
            for matrix in [self.csr, self.csc]:
                start, end = matrix.indptr[index:index + 2]
                for neighbor, value in zip(matrix.indices[start:end], matrix.data[start:end]):
                    scores[int(neighbor)] = max(scores.get(int(neighbor), -np.inf), float(value))
        return scores


class TableStream:
    """Append bounded groups of analytical rows with one stable TSV schema."""
    def __init__(self, path, columns):
        self.handle = Path(path).open("x", newline="")
        self.columns = columns
        self.writer = csv.DictWriter(self.handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        self.writer.writeheader()
        self.rows = 0

    def append(self, rows):
        self.writer.writerows(rows)
        self.rows += len(rows)

    def close(self):
        self.handle.close()


class ComparisonAccumulator:
    """Exact donor means over common held-out genes without retaining their rows."""
    def __init__(self, endpoints, methods):
        self.statistics = {}
        for endpoint in endpoints:
            baseline = "promoter_only" if endpoint == "rna" else "intercept"
            for method in methods:
                comparators = [baseline] + (["distance"] if method != "distance" else [])
                if method not in {"distance", "cicero"}:
                    comparators.append("cicero")
                for comparator in comparators:
                    self.statistics[(endpoint, method, comparator)] = {"donors": {}, "genes": set(), "rows": 0}

    def append(self, rows):
        # Caller supplies a single outer-fold gene, so each donor/method is unique.
        indexed = {(r["endpoint"], r["method"], r["donor_id"]): r for r in rows}
        if len(indexed) != len(rows):
            raise ValueError("Duplicate donor prediction within a held-out gene")
        donors = sorted({r["donor_id"] for r in rows})
        for (endpoint, method, baseline), summary in self.statistics.items():
            for donor in donors:
                candidate = indexed.get((endpoint, method, donor))
                reference = indexed.get((endpoint, baseline, donor))
                if candidate is None or reference is None:
                    continue
                if candidate["n_metacells"] != reference["n_metacells"]:
                    raise ValueError("Compared methods have different donor aggregate coverage")
                difference = reference["standardized_squared_error"] - candidate["standardized_squared_error"]
                total, count = summary["donors"].get(donor, (0.0, 0))
                summary["donors"][donor] = (total + difference, count + 1)
                summary["genes"].add(candidate["gene_id"])
                summary["rows"] += 1

    def finalize(self, seed, reps):
        rows, donor_rows = [], []
        for (endpoint, method, baseline), summary in self.statistics.items():
            differences = []
            for donor, (total, count) in sorted(summary["donors"].items()):
                differences.append(total / count)
                donor_rows.append(dict(endpoint=endpoint, method=method, baseline=baseline, donor_id=donor,
                                       mean_loss_difference=total / count, n_common_genes=count))
            statistics = paired_inference(differences, stable_seed(seed, (endpoint, method, baseline)), reps)
            rows.append(dict(endpoint=endpoint, method=method, baseline=baseline,
                             n_common_predictions=summary["rows"], n_common_genes=len(summary["genes"]), **statistics))
        return pd.DataFrame(rows), pd.DataFrame(donor_rows)


def worker_count(value, label):
    workers = int(value)
    if workers < 1:
        raise ValueError(f"{label} must be positive")
    allocation = os.environ.get("SLURM_CPUS_PER_TASK")
    if workers > 1 and (not allocation or workers > int(allocation)):
        raise ValueError(f"{label} exceeds SLURM_CPUS_PER_TASK")
    return workers


def cicero_scores(counts, peaks, output, rscript, seed, window, scoring_mask=None, distance_samples=100, donor_ids=None, workers=1):
    """Native Cicero on training-only, donor-confined aggregates."""
    workers = worker_count(workers, "cicero_workers")
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    with threadpool_limits(limits=1):
        mmwrite(output / "counts.mtx", counts.T)
    exported = peaks.copy()
    exported["score_peak"] = True if scoring_mask is None else scoring_mask
    exported.to_csv(output / "peaks.tsv", sep="\t", index=False)
    script = Path(__file__).with_name("network_cicero.R")
    command = [str(rscript), str(script), "--counts", str(output / "counts.mtx"), "--peaks", str(output / "peaks.tsv"), "--output", str(output / "links.tsv"), "--seed", str(seed), "--window", str(window), "--distance-samples", str(distance_samples), "--workers", str(workers)]
    if donor_ids is not None:
        if len(donor_ids) != counts.shape[0]:
            raise ValueError("Cicero donor identities do not match aggregate rows")
        pd.DataFrame(dict(metacell_id=[f"metacell_{i+1}" for i in range(len(donor_ids))],
                          donor_id=donor_ids)).to_csv(output / "donors.tsv", sep="\t", index=False)
        command.extend(["--donors", str(output / "donors.tsv")])
    with (output / "run.log").open("w") as handle:
        run = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if run.returncode:
        raise RuntimeError(f"Cicero failed ({run.returncode}); inspect {output / 'run.log'}")
    if donor_ids is not None and "covariance=donor_equal_reliability_weighted" not in Path(str(output / "links.tsv") + ".method.txt").read_text():
        raise RuntimeError("Cicero wrapper did not confirm donor-balanced calibration and covariance")
    graph = CiceroGraph.from_tsv(output / "links.tsv", peaks)
    (output / "timing.json").write_text(json.dumps({"elapsed_seconds": time.monotonic() - started,
        "training_metacells": counts.shape[0], "training_peaks": counts.shape[1],
        "scoring_peaks": len(peaks) if scoring_mask is None else int(np.sum(scoring_mask)),
        "directed_edges": graph.csr.nnz}, indent=2) + "\n")
    return graph


def selected_indices(method, candidates, promoter_indices, peaks, gene, train_x, train_y, train_donors, seed, max_links, cicero=None, external=None):
    if not len(candidates):
        return np.empty(0, dtype=int), np.empty(0)
    local = log_columns(train_x, candidates)
    weights = donor_weights(train_donors)
    mean = np.average(local, axis=0, weights=weights)
    valid = (np.average(local > 0, axis=0, weights=weights) >= 0.01) & (np.average((local - mean) ** 2, axis=0, weights=weights) > 1e-12)
    candidates, local = candidates[valid], local[:, valid]
    if not len(candidates):
        return candidates, np.empty(0)
    if method == "distance":
        score = -np.abs((peaks.iloc[candidates].start0.to_numpy() + peaks.iloc[candidates].end.to_numpy()) / 2 - int(gene["tss0"]))
    elif method == "directnet_adapted":
        from xgboost import XGBRegressor
        model = XGBRegressor(n_estimators=100, max_depth=3, learning_rate=0.1, objective="reg:squarederror", n_jobs=1, random_state=seed, importance_type="gain", subsample=1, colsample_bytree=1)
        model.fit(local, train_y, sample_weight=donor_weights(train_donors))
        score = model.feature_importances_
    elif method == "cicero":
        if cicero is None:
            return np.empty(0, dtype=int), np.empty(0)
        if isinstance(cicero, CiceroGraph):
            scores = cicero.promoter_scores(promoter_indices)
            score = np.array([scores.get(int(x), -np.inf) for x in candidates])
        else:
            prom_names = set(peaks.iloc[promoter_indices].peak_id)
            scores = {}
            for row in cicero.itertuples(index=False):
                if row.peak1 in prom_names:
                    scores[row.peak2] = max(scores.get(row.peak2, -np.inf), float(row.coaccess))
                if row.peak2 in prom_names:
                    scores[row.peak1] = max(scores.get(row.peak1, -np.inf), float(row.coaccess))
            score = np.array([scores.get(x, -np.inf) for x in peaks.iloc[candidates].peak_id])
    elif external is not None:
        gene_links = external[external.gene_id == gene["gene_id"]]
        lookup = dict(zip(gene_links.peak_id, gene_links.score))
        score = np.array([lookup.get(x, -np.inf) for x in peaks.iloc[candidates].peak_id])
    else:
        raise ValueError(method)
    keep = np.isfinite(score)
    if method != "distance":
        # An exact zero is no learned link, not a feature for padding the ridge
        # design. Signed nonzero coaccessibility remains an evaluable coupling.
        keep &= np.asarray(score) != 0
    candidates, score = candidates[keep], np.asarray(score)[keep]
    ranked = np.lexsort((np.asarray(peaks.iloc[candidates].peak_id), -score))[:max_links]
    return candidates[ranked], score[ranked]


def design_matrix(x, selected, promoter_values, endpoint):
    distal = log_columns(x, selected)
    return np.column_stack([promoter_values, distal]) if endpoint == "rna" else distal


def tune_alpha(design, y, donors, train, inner, grid):
    losses = {}
    for alpha in grid:
        values = []
        for k in sorted(set(inner)):
            it, iv = train[inner != k], train[inner == k]
            pred = weighted_ridge(design[it], y[it], donors[it], alpha, design[iv])
            values.extend(pd.DataFrame({"donor": donors[iv], "loss": (pred - y[iv]) ** 2}).groupby("donor").loss.mean().tolist())
        losses[alpha] = np.mean(values)
    return min(grid, key=lambda a: (losses[a], a))


def save_predictions(rows, fold, gene, endpoint, method, indices, metadata, y, prediction, variance, alpha):
    # Reduce to the biological replication unit immediately. Squaring the mean
    # residual would incorrectly discard within-donor prediction error.
    held_donors = metadata.iloc[indices].donor_id.to_numpy()
    observed = y[indices]
    for donor in sorted(set(held_donors)):
        selected = held_donors == donor
        mse = float(np.mean((prediction[selected] - observed[selected]) ** 2))
        rows.append(dict(outer_fold=fold, gene_id=gene, endpoint=endpoint, method=method, donor_id=donor, n_metacells=int(selected.sum()), observed=float(observed[selected].mean()), predicted=float(prediction[selected].mean()), squared_error=mse, standardized_squared_error=mse / variance, ridge_alpha=alpha))


def compare_predictions(predictions, methods, seed, reps, endpoints=None):
    rows = []
    for endpoint in (sorted(predictions.endpoint.unique()) if endpoints is None else endpoints):
        baseline = "promoter_only" if endpoint == "rna" else "intercept"
        for method in methods:
            comparators = [baseline] + (["distance"] if method != "distance" else [])
            if method == "directnet_adapted":
                comparators.append("cicero")
            for comparator in comparators:
                a = predictions[(predictions.endpoint == endpoint) & (predictions.method == method)]
                b = predictions[(predictions.endpoint == endpoint) & (predictions.method == comparator)]
                keys = ["outer_fold", "gene_id", "donor_id"]
                joined = a.merge(b, on=keys, suffixes=("_candidate", "_baseline"), validate="one_to_one")
                if len(joined):
                    if "n_metacells_candidate" in joined and not joined.n_metacells_candidate.equals(joined.n_metacells_baseline):
                        raise ValueError("Compared methods have different donor aggregate coverage")
                    joined["difference"] = joined.standardized_squared_error_baseline - joined.standardized_squared_error_candidate
                    differences = joined.groupby(["donor_id", "gene_id"]).difference.mean().groupby("donor_id").mean()
                    stats = paired_inference(differences, stable_seed(seed, (endpoint, method, comparator)), reps)
                else:
                    stats = paired_inference([], seed, reps)
                rows.append(dict(endpoint=endpoint, method=method, baseline=comparator, n_common_predictions=len(joined), n_common_genes=joined.gene_id.nunique(), **stats))
    return pd.DataFrame(rows)


def _fit_gene_readout(job, context):
    """Pure gene calculation; table writing and donor reduction stay in the parent."""
    gene_number, gene = job
    gene_started = time.monotonic()
    (x, feature_x, donors, metadata, peaks, geometry, train, test, inner, methods,
     external_settings, cicero_cache, external_cache, training_x, inner_training_x,
     alpha_grid, max_links, window, promoter_bp, seed, fold, rna_normal, rna_lookup) = (
        context[key] for key in ("x", "feature_x", "donors", "metadata", "peaks", "geometry",
            "train", "test", "inner", "methods", "external_settings", "cicero_cache",
            "external_cache", "training_x", "inner_training_x", "alpha_grid", "max_links",
            "window", "promoter_bp", "seed", "fold", "rna_normal", "rna_lookup"))
    predictions, links = [], []
    gene_id = gene["gene_id"]
    prom, candidates = geometry.query(gene, window, promoter_bp)
    result = dict(gene_number=gene_number, predictions=predictions, links=links, timing=None)
    coverage = dict(outer_fold=fold, gene_id=gene_id, n_candidates=len(candidates))
    if not len(prom) or not len(candidates):
        result["coverage"] = dict(coverage, status="untestable", reason="missing_promoter_or_distal_candidates")
        return result
    # Retain the original CSR sum for bit-identical promoter target arithmetic.
    promoter = np.log1p(np.asarray(x[:, prom].sum(axis=1)).ravel())
    if np.var(promoter[train]) <= 1e-12:
        result["coverage"] = dict(coverage, status="untestable", reason="constant_training_promoter")
        return result
    endpoints = {"promoter_accessibility": promoter}
    if rna_normal is not None and gene_id in rna_lookup:
        endpoints["rna"] = log_columns(rna_normal, [rna_lookup[gene_id]]).ravel()
    result["coverage"] = dict(coverage, status="testable", reason="")
    # Every method learns links against the same promoter target for both readouts.
    selected_cache = {}
    for endpoint, y in endpoints.items():
        variance = float(np.average((y[train] - np.average(y[train], weights=donor_weights(donors[train]))) ** 2, weights=donor_weights(donors[train])))
        if variance <= 1e-12:
            continue
        baseline_x = promoter[:, None] if endpoint == "rna" else np.empty((len(donors), 0))
        baseline_alpha = tune_alpha(baseline_x, y, donors, train, inner, alpha_grid)
        base_pred = weighted_ridge(baseline_x[train], y[train], donors[train], baseline_alpha, baseline_x[test])
        save_predictions(predictions, fold, gene_id, endpoint, "promoter_only" if endpoint == "rna" else "intercept", test, metadata, y, base_pred, variance, baseline_alpha)
        for method in methods:
            if method == "cicero" and any(v is None for v in cicero_cache.values()):
                continue
            if method in external_settings and any(external_cache[(method, tag)] is None for tag in ["outer"] + [f"inner{k}" for k in sorted(set(inner))]):
                continue
            if method not in selected_cache:
                chosen, scores = selected_indices(method, candidates, prom, peaks, gene, training_x, promoter[train], donors[train], stable_seed(seed, (fold, gene_id)), max_links, cicero_cache.get("outer"), external_cache.get((method, "outer")))
                inner_selected = {}
                if len(chosen):
                    for k in sorted(set(inner)):
                        it = train[inner != k]
                        inner_selected[k], _ = selected_indices(method, candidates, prom, peaks, gene, inner_training_x[k], promoter[it], donors[it], stable_seed(seed, (fold, k, gene_id)), max_links, cicero_cache.get(f"inner{k}"), external_cache.get((method, f"inner{k}")))
                selected_cache[method] = (chosen, scores, inner_selected)
            chosen, scores, inner_selected = selected_cache[method]
            if not len(chosen):
                continue
            validation_losses = {alpha: [] for alpha in alpha_grid}
            for k in sorted(set(inner)):
                it, iv = train[inner != k], train[inner == k]
                design = design_matrix(feature_x, inner_selected[k], promoter, endpoint)
                for alpha in alpha_grid:
                    pred = weighted_ridge(design[it], y[it], donors[it], alpha, design[iv])
                    loss = pd.DataFrame({"donor": donors[iv], "loss": (pred - y[iv]) ** 2}).groupby("donor").loss.mean()
                    validation_losses[alpha].extend(loss.tolist())
            best = min(alpha_grid, key=lambda a: (np.mean(validation_losses[a]), a))
            design = design_matrix(feature_x, chosen, promoter, endpoint)
            pred = weighted_ridge(design[train], y[train], donors[train], best, design[test])
            save_predictions(predictions, fold, gene_id, endpoint, method, test, metadata, y, pred, variance, best)
            if endpoint == "promoter_accessibility":
                for rank, (idx, score) in enumerate(zip(chosen, scores), start=1):
                    links.append(dict(outer_fold=fold, gene_id=gene_id, peak_id=peaks.iloc[idx].peak_id, method=method, rank=rank, score=float(score), score_unit="external_method_unscaled_score" if method in external_settings else "genomic_distance_negative_bp" if method == "distance" else "xgboost_gain_fraction" if method == "directnet_adapted" else "donor_weighted_Cicero_coaccessibility", evidence_state="candidate_accessibility_link"))
    result["timing"] = dict(outer_fold=fold, gene_id=gene_id, elapsed_seconds=time.monotonic()-gene_started,
        n_endpoints=len(endpoints), n_prediction_rows=len(predictions), n_link_rows=len(links))
    return result


_READOUT_CONTEXT = None
_READOUT_THREAD_LIMIT = None


def _readout_worker(job):
    try:
        return _fit_gene_readout(job, _READOUT_CONTEXT)
    except Exception as exc:
        raise RuntimeError(f"Gene readout failed at gene {job[1]['gene_id']}") from exc


def _limit_readout_threads():
    global _READOUT_THREAD_LIMIT
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[name] = "1"
    _READOUT_THREAD_LIMIT = threadpool_limits(limits=1)


def _ordered_gene_readouts(jobs, context, workers):
    """Bound work/results to four genes per worker and retain annotation order."""
    if workers == 1:
        for job in jobs:
            yield _fit_gene_readout(job, context)
        return
    global _READOUT_CONTEXT
    _READOUT_CONTEXT = MappingProxyType(context)
    # Analysis helpers also load this module via importlib without registering it.
    # A fork-inherited alias lets Pool pickle only this function, not the arrays.
    alias = f"_masld_readout_{id(_readout_worker)}"
    registry = ModuleType(alias)
    registry._readout_worker = _readout_worker
    _readout_worker.__module__ = alias
    sys.modules[alias] = registry
    try:
        with multiprocessing.get_context("fork").Pool(workers, initializer=_limit_readout_threads) as pool:
            for start in range(0, len(jobs), 4 * workers):
                yield from pool.imap(_readout_worker, jobs[start:start + 4 * workers], chunksize=4)
    finally:
        _READOUT_CONTEXT = None


def selected_folds(settings, available):
    available = sorted(set(map(int, available)))
    requested = settings.get("selected_outer_folds")
    if requested is None:
        return available
    if not isinstance(requested, list) or not requested or any(type(f) is not int for f in requested):
        raise ValueError("selected_outer_folds must be a nonempty list of integer fold IDs")
    if len(set(requested)) != len(requested) or set(requested) - set(available):
        raise ValueError("selected_outer_folds contains duplicate or unavailable fold IDs")
    return sorted(requested)


def scientific_settings(settings):
    return {key: value for key, value in settings.items()
            if key not in {"selected_outer_folds", "readout_workers", "cicero_workers"}}


def table_digest(frame):
    return hashlib.sha256(frame.to_csv(sep="\t", index=False).encode()).hexdigest()


def matrix_digest(matrix):
    digest = hashlib.sha256(str(matrix.shape).encode())
    for array in (matrix.data, matrix.indices, matrix.indptr):
        digest.update(str(array.dtype).encode())
        digest.update(memoryview(array).cast("B"))
    return digest.hexdigest()


def fit_group(counts, metadata, peaks, genes, output, settings, seed, rna=None, rna_genes=None):
    readout_workers = worker_count(settings.get("readout_workers", 1), "readout_workers")
    output.mkdir(parents=True, exist_ok=False)
    x = normalize_sparse(counts)
    feature_x = x.tocsc()
    donors = metadata.donor_id.to_numpy()
    if settings.get("outer_fold_map") is not None:
        mapping = settings["outer_fold_map"]
        missing = set(donors) - set(mapping)
        if missing:
            raise ValueError(f"Fixed outer fold map omits donors: {sorted(missing)}")
        outer = np.array([int(mapping[str(donor)]) for donor in donors], dtype=int)
        if len(set(outer)) < 2 or (outer < 0).any():
            raise ValueError("Fixed outer fold map needs at least two nonnegative folds")
    else:
        outer = donor_folds(donors, int(settings.get("outer_folds", 5)), seed)
    folds = metadata.copy()
    folds["outer_fold"] = outer
    requested_folds = selected_folds(settings, outer)
    folds.to_csv(output / "fold_assignments.tsv", sep="\t", index=False)
    geometry = PeakGeometry(peaks) if len(peaks) else None
    # A deterministic annotation order, never phenotype, specifies a bounded pilot.
    max_genes = int(settings.get("max_genes", 0))
    genes = genes.sort_values("gene_id")
    if max_genes > 0:
        # Geometry-only selection does not inspect held-out counts or disease labels.
        supported = []
        for gene in genes.to_dict("records"):
            promoter, distal = geometry.query(gene)
            if len(promoter) and len(distal):
                supported.append(gene)
                if len(supported) == max_genes:
                    break
        genes = pd.DataFrame(supported, columns=genes.columns)
    native_paired = rna is not None or ("cohort" in metadata and metadata.cohort.eq("GSE296875").all())
    estimated_rows = len(genes) * len(set(donors)) * 4 * (2 if native_paired else 1)
    maximum_rows = int(settings.get("max_prediction_rows", 0))
    projection = {"candidate_genes": len(genes), "donors": len(set(donors)),
                  "estimated_donor_prediction_rows_upper_bound": estimated_rows,
                  "max_prediction_rows": maximum_rows,
                  "estimated_prediction_working_bytes_upper_bound": len(set(donors)) * 8 * 2048,
                  "prediction_storage": "streamed_one_gene_fold_at_a_time",
                  "inference_unit": "donor", "loss_aggregation": "mean_of_squared_errors"}
    (output / "memory_projection.json").write_text(json.dumps(projection, indent=2) + "\n")
    if maximum_rows > 0 and estimated_rows > maximum_rows:
        raise MemoryError("Prediction output exceeds explicit max_prediction_rows disk-output policy")
    dependency = inventory()
    methods = ["distance"]
    statuses = [{"method": "distance", "status": "available", "reason": ""}]
    if dependency["xgboost"]:
        methods.append("directnet_adapted")
    statuses.append({"method": "directnet_adapted", "status": "available" if dependency["xgboost"] else "not_run", "reason": "" if dependency["xgboost"] else "xgboost_missing_in_active_environment"})
    rscript = settings.get("cicero_rscript", "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript")
    if rscript and Path(rscript).exists():
        methods.append("cicero")
    statuses.append({"method": "cicero", "status": "available" if "cicero" in methods else "not_run", "reason": "" if "cicero" in methods else "cicero_rscript_not_configured", "method_variant": "donor_balanced_Cicero_adaptation"})
    external_settings = settings.get("external_stage_links", {})
    if set(external_settings) & {"distance", "directnet_adapted", "cicero", "intercept", "promoter_only"}:
        raise ValueError("External methods require distinct labels")
    for method in external_settings:
        methods.append(method)
        statuses.append(dict(method=method, status="available", reason="external_training_stage_links_exploratory"))
    for method in ["directnet_native", "LINGER", "TRIPOD"]:
        statuses.append({"method": method, "status": "not_run", "reason": "conditional_dependency_and_paired_data_stage"})
    pd.DataFrame(statuses).to_csv(output / "method_status.tsv", sep="\t", index=False)
    (output / "dependencies.json").write_text(json.dumps(dependency, indent=2) + "\n")
    failures = []
    common = dict(cohort=metadata.cohort.iloc[0] if "cohort" in metadata else "fixture",
                  lineage=metadata.lineage.iloc[0] if "lineage" in metadata else "fixture",
                  lineage_state=";".join(sorted(metadata.lineage_state.unique())) if "lineage_state" in metadata else "source_label",
                  input_representation=settings.get("representation", "binary"))
    prediction_columns = ["outer_fold", "gene_id", "endpoint", "method", "donor_id", "n_metacells", "observed", "predicted", "squared_error", "standardized_squared_error", "ridge_alpha", "target_unit"]
    prediction_stream = TableStream(output / "predictions.tsv", prediction_columns + list(common))
    coverage_stream = TableStream(output / "coverage.tsv", ["outer_fold", "gene_id", "status", "reason", "n_candidates"] + list(common))
    link_stream = TableStream(output / "links.tsv", ["outer_fold", "gene_id", "peak_id", "method", "rank", "score", "score_unit", "evidence_state"] + list(common))
    timing_stream = TableStream(output / "gene_timing.tsv", ["outer_fold", "gene_id", "elapsed_seconds", "n_endpoints", "n_prediction_rows", "n_link_rows"])
    planned_endpoints = ["promoter_accessibility", "rna"] if native_paired else ["promoter_accessibility"]
    accumulator = ComparisonAccumulator(planned_endpoints, ["distance", "directnet_adapted", "cicero"])
    external_accumulator = ComparisonAccumulator(planned_endpoints, list(external_settings))
    for candidate, control in settings.get("external_method_controls", {}).items():
        if candidate not in external_settings or control not in external_settings:
            raise ValueError("External control contrasts require both configured method stages")
        for endpoint in planned_endpoints:
            external_accumulator.statistics[(endpoint, candidate, control)] = {"donors": {}, "genes": set(), "rows": 0}
    scored_counts, scored_folds = {}, {}
    def write_coverage(row):
        coverage_stream.append([{**row, **common}])
    rna_lookup = {} if rna_genes is None else {g: i for i, g in enumerate(rna_genes.gene_id)}
    rna_normal = None if rna is None else normalize_sparse(rna).tocsc()
    alpha_grid = [float(a) for a in settings.get("ridge_alphas", [0.1, 1, 10, 100])]
    max_links = int(settings.get("max_distal_peaks", settings.get("max_links", 10)))
    window = int(settings.get("cis_window", settings.get("cis_window_bp", 250000)))
    promoter_bp = int(settings.get("promoter_flank", settings.get("promoter_bp", 2000)))
    if promoter_bp > window:
        raise ValueError("Promoter flank must not exceed the cis window")
    cicero_universe = np.zeros(len(peaks), dtype=bool)
    for gene in genes.to_dict("records"):
        promoter_indices, distal = geometry.query(gene, window, promoter_bp)
        cicero_universe[promoter_indices] = True
        cicero_universe[distal] = True
    contract = dict(seed=seed, settings=scientific_settings(settings),
        available_folds=sorted(set(map(int, outer))), gene_ids=genes.gene_id.tolist(),
        counts_sha256=matrix_digest(counts), metadata_sha256=table_digest(metadata),
        peaks_sha256=table_digest(peaks), genes_sha256=table_digest(genes),
        rna_sha256=None if rna is None else matrix_digest(rna),
        rna_genes_sha256=None if rna_genes is None else table_digest(rna_genes),
        fold_assignments_sha256=table_digest(folds),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        planned_endpoints=planned_endpoints)
    (output / "fold_contract.json").write_text(json.dumps(contract, indent=2) + "\n")
    for fold in requested_folds:
        train, test = np.flatnonzero(outer != fold), np.flatnonzero(outer == fold)
        inner = donor_folds(donors[train], int(settings.get("inner_folds", 3)), stable_seed(seed, ("inner", fold)))
        # These exact row subsets are shared by every gene, method and endpoint.
        # Cache each once rather than copying the complete sparse peak matrix
        # before every local cis-feature query.
        training_x = x[train].tocsc()
        inner_training_x = {k: x[train[inner != k]].tocsc() for k in sorted(set(inner))}
        cicero_cache = {}
        external_cache = {}
        for method, paths in external_settings.items():
            for tag, fit_idx in [("outer", train)] + [(f"inner{k}", train[inner != k]) for k in sorted(set(inner))]:
                try:
                    path = Path(paths[f"{fold}:{tag}"])
                    run_root = Path(settings.get("external_stage_roots", {})[method])
                    for authority in {path.parent.resolve(), run_root.resolve()}:
                        if any(authority.glob("INVALIDATED*.json")):
                            raise ValueError("External stage or parent run is invalidated")
                        completion = json.loads((authority / "completion.json").read_text())
                        if completion.get("status") != "completed":
                            raise ValueError("External stage or parent run has not completed successfully")
                    if not path.resolve().is_relative_to(run_root.resolve()):
                        raise ValueError("External link stage is outside its declared parent run")
                    observed_donors = json.loads(path.with_suffix(".training_donors.json").read_text())["donor_ids"]
                    if set(map(str, observed_donors)) != set(donors[fit_idx]):
                        raise ValueError("External link donor set does not match this training stage")
                    frame = read_table(path)
                    if not {"gene_id", "peak_id", "score"}.issubset(frame):
                        raise ValueError("External links require gene_id, peak_id, score")
                    frame.score = frame.score.astype(float)
                    if frame.duplicated(["gene_id", "peak_id"]).any() or not np.isfinite(frame.score).all():
                        raise ValueError("External links must be unique and finite")
                    if not frame.peak_id.isin(peaks.peak_id).all():
                        raise ValueError("External peak IDs differ from the common peak universe")
                    external_cache[(method, tag)] = frame
                except Exception as exc:
                    external_cache[(method, tag)] = None
                    failures.append(dict(outer_fold=fold, method=method, stage=tag, reason=str(exc)))
        if "cicero" in methods:
            for tag, fit_idx in [("outer", train)] + [(f"inner{k}", train[inner != k]) for k in sorted(set(inner))]:
                try:
                    nonzero = np.flatnonzero(np.asarray(counts[fit_idx].sum(axis=0)).ravel() > 0)
                    print(f"Cicero fold={fold} stage={tag} training_donors={len(set(donors[fit_idx]))} peaks={len(nonzero)}", flush=True)
                    cicero_cache[tag] = cicero_scores(counts[fit_idx][:, nonzero], peaks.iloc[nonzero], output / f"cicero_fold{fold}_{tag}", rscript, stable_seed(seed, (fold, tag)), window * 2, cicero_universe[nonzero], int(settings.get("cicero_distance_samples", 100)), donors[fit_idx], worker_count(settings.get("cicero_workers", 1), "cicero_workers"))
                except Exception as exc:
                    cicero_cache[tag] = None
                    failures.append(dict(outer_fold=fold, method="cicero", stage=tag, reason=str(exc)))
        context = dict(x=x, feature_x=feature_x, donors=donors, metadata=metadata, peaks=peaks,
            geometry=geometry, train=train, test=test, inner=inner, methods=methods,
            external_settings=external_settings, cicero_cache=cicero_cache, external_cache=external_cache,
            training_x=training_x, inner_training_x=inner_training_x, alpha_grid=alpha_grid,
            max_links=max_links, window=window, promoter_bp=promoter_bp, seed=seed, fold=fold,
            rna_normal=rna_normal, rna_lookup=rna_lookup)
        jobs = list(enumerate(genes.to_dict("records"), start=1))
        for result in _ordered_gene_readouts(jobs, context, readout_workers):
            write_coverage(result["coverage"])
            if result["timing"] is None:
                continue
            gene_number = result["gene_number"]
            predictions, links = result["predictions"], result["links"]
            accumulator.append(predictions)
            external_accumulator.append(predictions)
            for row in predictions:
                row.update(common)
                raw_unit = "log1p_native_author_peak_count_events_per_10000" if common["cohort"] == "GSE296875" else "log1p_fragment_counts_per_10000"
                row["target_unit"] = "log1p_RNA_counts_per_10000" if row["endpoint"] == "rna" else "log1p_accessible_nucleus_counts_per_10000" if common["input_representation"] == "binary" else raw_unit
                method = row["method"]
                scored_counts[method] = scored_counts.get(method, 0) + 1
                scored_folds.setdefault(method, set()).add(fold)
            prediction_stream.append(predictions)
            link_stream.append([{**row, **common} for row in links])
            timing_stream.append([result["timing"]])
            if gene_number % 100 == 0 or gene_number == len(genes):
                print(f"Gene readout fold={fold} processed={gene_number}/{len(genes)} donor_prediction_rows={prediction_stream.rows}", flush=True)
        for stream in [prediction_stream, coverage_stream, link_stream, timing_stream]:
            stream.handle.flush()
        (output / f"fold{fold}_completed.json").write_text(json.dumps(dict(
            status="completed_fold", outer_fold=int(fold), n_genes=len(genes),
            failures=[row for row in failures if row["outer_fold"] == fold]), indent=2) + "\n")
    for stream in [prediction_stream, coverage_stream, link_stream, timing_stream]:
        stream.close()
    pd.DataFrame(failures, columns=["outer_fold", "method", "stage", "reason"]).to_csv(output / "failures.tsv", sep="\t", index=False)
    summary = []
    for row in statuses:
        method = row["method"]
        failures_for_method = [x for x in failures if x["method"] == method]
        summary.append(dict(method=method, eligibility=row["status"], n_prediction_rows=scored_counts.get(method, 0), n_scored_outer_folds=len(scored_folds.get(method, set())), n_failed_stages=len(failures_for_method), execution_status="partial" if failures_for_method else "completed" if scored_counts.get(method, 0) else "not_scored", reason=row["reason"]))
    pd.DataFrame(summary).to_csv(output / "method_execution_summary.tsv", sep="\t", index=False)
    comparisons, donor_errors = accumulator.finalize(seed, int(settings.get("bootstrap_reps", 2000)))
    if external_settings:
        exploratory, exploratory_donors = external_accumulator.finalize(seed, int(settings.get("bootstrap_reps", 2000)))
        exploratory["family_complete"] = False
        exploratory["predictive_improvement"] = False
        exploratory["evaluation_status"] = "exploratory_native_stage_comparison_not_in_prespecified_primary_family"
        exploratory["holm_pvalue"] = holm(exploratory.pvalue.to_numpy())
        exploratory["family"] = "exploratory_native_methods_all_endpoint_baseline_control_contrasts_within_context"
        exploratory.to_csv(output / "exploratory_comparisons.tsv", sep="\t", index=False)
        exploratory_donors.to_csv(output / "exploratory_donor_errors.tsv", sep="\t", index=False)
    donor_errors.to_csv(output / "donor_errors.tsv", sep="\t", index=False)
    comparisons.to_csv(output / "comparisons.tsv", sep="\t", index=False)
    if settings.get("selected_outer_folds") is not None:
        (output / "shard_completion.json").write_text(json.dumps(dict(
            status="completed" if not failures else "failed", selected_outer_folds=requested_folds,
            all_outer_folds=False, scientific_adoption=False), indent=2) + "\n")
    return comparisons


def native_directnet_diagnostic(counts, cells, peaks, genes, output, settings, seed, rscript, api_compatibility=False):
    """Bounded native call; its aggregates and optional splits are not donor held out."""
    keep, exclusions = condition_eligibility(cells, settings)
    exclusions.to_csv(output / "native_condition_exclusions.tsv", sep="\t", index=False)
    indices = np.flatnonzero(keep)
    if len(cells.iloc[indices][["cohort", "lineage"]].drop_duplicates()) != 1:
        raise ValueError("Native diagnostic requires one cohort and lineage per invocation")
    rng = np.random.default_rng(seed)
    donors = cells.iloc[indices].donor_id.to_numpy()
    # Balanced donor round-robin before native neighbor aggregation.
    buckets = [rng.permutation(indices[donors == donor]).tolist() for donor in sorted(set(donors))]
    chosen = []
    while buckets and len(chosen) < 1000:
        chosen.extend(bucket.pop() for bucket in buckets if bucket)
        buckets = [bucket for bucket in buckets if bucket]
    chosen = np.array(chosen[:1000], dtype=int)
    geometry = PeakGeometry(peaks)
    targets = []
    for gene in genes.sort_values("gene_id").to_dict("records"):
        promoter, distal = geometry.query(gene)
        if len(promoter) and len(distal):
            targets.append(gene)
        if len(targets) == 5:
            break
    if not targets or not {"start0", "end"}.issubset(genes):
        raise ValueError("Native diagnostic requires gene intervals and geometrically eligible targets")
    with threadpool_limits(limits=1):
        mmwrite(output / "native_counts.mtx", counts[chosen].T)
    cells.iloc[chosen].to_csv(output / "native_cells.tsv", sep="\t", index=False)
    peaks.to_csv(output / "native_peaks.tsv", sep="\t", index=False)
    pd.DataFrame(targets).to_csv(output / "native_genes.tsv", sep="\t", index=False)
    script = Path(__file__).with_name("network_directnet_native.R")
    command = [str(rscript), str(script), "--counts", str(output / "native_counts.mtx"),
               "--cells", str(output / "native_cells.tsv"), "--peaks", str(output / "native_peaks.tsv"),
               "--genes", str(output / "native_genes.tsv"), "--output", str(output / "native"), "--seed", str(seed)]
    if api_compatibility:
        command.append("--xgboost-api-compatibility")
    environment = os.environ.copy()
    prefix = Path(rscript).resolve().parent.parent
    environment["LD_LIBRARY_PATH"] = str(prefix / "lib") + (":" + environment["LD_LIBRARY_PATH"] if environment.get("LD_LIBRARY_PATH") else "")
    environment.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with (output / "native_run.log").open("w") as handle:
        run = subprocess.run(command, env=environment, stdout=handle, stderr=subprocess.STDOUT, check=False)
    native_status = output / "native" / "completion.json"
    status = json.loads(native_status.read_text()).get("status", "failed") if native_status.exists() else "failed"
    (output / "completion.json").write_text(json.dumps({"status": status if run.returncode == 0 else "failed",
        "interpretation": "native_technical_comparison_only", "donor_held_out": False,
        "xgboost_API_compatibility_shim": api_compatibility,
        "scientific_adoption": False, "returncode": run.returncode}, indent=2) + "\n")
    if run.returncode:
        raise RuntimeError("Native DIRECT-NET diagnostic failed; inspect native_run.log and native/completion.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--genes", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--rna-root", type=Path)
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--max-genes", type=int)
    parser.add_argument("--bootstrap-reps", type=int)
    parser.add_argument("--cicero-distance-samples", type=int)
    parser.add_argument("--native-directnet-only", action="store_true")
    parser.add_argument("--directnet-rscript", type=Path)
    parser.add_argument("--native-directnet-api-compat", action="store_true")
    args = parser.parse_args()
    require_compute()
    config = json.loads(args.config.read_text())
    settings = dict(config.get("network", {}))
    settings.setdefault("representation", "binary")
    settings.setdefault("cicero_distance_samples", 100)
    seed = int(config.get("seed", 42))
    source_pilot = False
    for name in ["diagnostics.json", "paired_export.json"]:
        evidence = args.input_root / name
        if evidence.exists():
            source = json.loads(evidence.read_text())
            source_pilot = source_pilot or bool(source.get("pilot", False)) or source.get("mode") == "pilot"
    technical_pilot = args.pilot or source_pilot
    if technical_pilot:
        maximum = int(settings.get("max_genes", 0))
        settings["max_genes"] = min(maximum, int(settings.get("pilot_max_genes", 25))) if maximum > 0 else int(settings.get("pilot_max_genes", 25))
        settings["cicero_distance_samples"] = int(settings.get("pilot_cicero_distance_samples", 10))
    if args.max_genes is not None:
        settings["max_genes"] = args.max_genes
    if args.bootstrap_reps is not None:
        settings["bootstrap_reps"] = args.bootstrap_reps
    if args.cicero_distance_samples is not None:
        if args.cicero_distance_samples < 1:
            raise ValueError("Cicero distance samples must be positive")
        settings["cicero_distance_samples"] = args.cicero_distance_samples
    if args.output.exists():
        raise FileExistsError(f"Refusing overwrite: {args.output}")
    counts, cells, peaks = load_input(args.input_root, seed, int(settings.get("max_cells_per_donor", 1000)))
    genes = read_table(args.genes)
    if not {"gene_id", "chrom", "tss0", "strand"}.issubset(genes) or genes.gene_id.duplicated().any():
        raise ValueError("genes.tsv requires unique gene_id/chrom/tss0/strand")
    genes.tss0 = genes.tss0.astype(np.int64)
    rna, rna_genes = None, None
    if args.rna_root:
        if set(cells.cohort) & {"GSE244832", "GSE281367"}:
            raise ValueError("These disease cohorts do not have qualified same-nucleus RNA/ATAC")
        rna_cells = read_table(args.rna_root / "cells.tsv")
        if not cell_keys(cells).equals(cell_keys(rna_cells)):
            raise ValueError("RNA/ATAC cell identities or order do not match exactly")
        rna = read_matrix(args.rna_root, ["rna_counts.npz", "matrix.npz", "RNA_counts.mtx.gz", "RNA_counts.mtx"])
        rna_genes = read_table(args.rna_root / "genes.tsv")
        if rna.shape != (len(cells), len(rna_genes)) or rna_genes.gene_id.duplicated().any():
            raise ValueError("RNA dimensions or gene identities are invalid")
    args.output.mkdir(parents=True, exist_ok=False)
    versions = {}
    snapshots = args.output / "source_snapshot"
    snapshots.mkdir(exist_ok=False)
    for source in [Path(__file__), Path(__file__).with_name("network_cicero.R"), Path(__file__).with_name("network_directnet_native.R"), args.config, args.genes]:
        content = source.read_bytes()
        versions[str(source.resolve())] = hashlib.sha256(content).hexdigest()
        if source.suffix in {".py", ".R", ".json"}:
            with (snapshots / source.name).open("xb") as handle:
                handle.write(content)
    (args.output / "source_sha256.json").write_text(json.dumps(versions, indent=2) + "\n")
    if args.native_directnet_only:
        if not args.directnet_rscript or not args.directnet_rscript.exists():
            raise ValueError("Native DIRECT-NET requires explicit --directnet-rscript from its approved environment")
        (args.output / "run_settings.json").write_text(json.dumps({"seed": seed, **settings,
            "pilot": True, "source_pilot": source_pilot, "analysis_scope": "native_technical_comparison_only",
            "donor_held_out": False, "native_aggregation_donor_mixing_allowed": True}, indent=2) + "\n")
        native_directnet_diagnostic(counts, cells, peaks, genes, args.output, settings, seed, args.directnet_rscript, args.native_directnet_api_compat)
        return
    (args.output / "run_settings.json").write_text(json.dumps({"seed": seed, **settings, "aggregation": "donor_confined_disjoint_random_50_cells", "pilot": technical_pilot, "source_pilot": source_pilot, "analysis_scope": "technical_pilot" if technical_pilot else "bounded_gene_development" if int(settings.get("max_genes", 0)) > 0 else "genome_wide", "inference_unit": "donor", "pvalue_assumption": "symmetric independent donor paired differences"}, indent=2) + "\n")
    condition_mask, excluded_conditions = condition_eligibility(cells, settings, technical_pilot)
    excluded_conditions.to_csv(args.output / "condition_exclusions.tsv", sep="\t", index=False)
    all_comparisons = []
    execution_qualifications = []
    for (cohort, lineage), indexes in cells.groupby(["cohort", "lineage"], sort=True).indices.items():
        if lineage not in settings.get("lineages", ["hepatocyte", "stellate", "mesenchymal"]):
            continue
        idx = np.asarray(indexes)
        idx = idx[condition_mask[idx]]
        if not len(idx):
            continue
        local_cells = cells.iloc[idx].reset_index(drop=True)
        local_counts = counts[idx]
        if technical_pilot:
            # Technical limit only; no inferential MASH claim follows from this subset.
            donors = sorted(local_cells.donor_id.unique())[:10]
            cap = min(1000, max(100, (1000 // max(len(donors), 1) // 50) * 50))
            selected = []
            for donor in donors:
                source = np.flatnonzero(local_cells.donor_id.to_numpy() == donor)
                selected.extend(np.random.default_rng(stable_seed(seed, donor)).permutation(source)[:cap])
            local_counts = local_counts[selected]
            idx = idx[selected]
            local_cells = local_cells.iloc[selected].reset_index(drop=True)
        representation = settings.get("representation", "binary")
        if representation == "binary":
            local_counts = local_counts.copy()
            local_counts.eliminate_zeros()
            local_counts.data[:] = 1
        elif representation != "fragment":
            raise ValueError("network.representation must be binary or fragment")
        aggregate, meta, groups, eligible = build_metacells(local_counts, local_cells, seed, int(settings.get("metacell_size", 50)), int(settings.get("max_cells_per_donor", 1000)))
        name = f"{cohort}_{lineage}"
        eligible.to_csv(args.output / f"{name}_eligibility.tsv", sep="\t", index=False)
        if meta.empty or meta.donor_id.nunique() < int(settings.get("outer_folds", 5)):
            (args.output / f"{name}_untestable.json").write_text(json.dumps({"status": "untestable", "reason": "insufficient_eligible_donors_for_grouped_folds"}) + "\n")
            continue
        membership = []
        for metacell_index, group in enumerate(groups):
            for row_index in group:
                original = int(local_cells.iloc[row_index].original_source_index) if "original_source_index" in local_cells else int(idx[row_index])
                membership.append(dict(cohort=cohort, lineage=lineage, donor_id=local_cells.iloc[row_index].donor_id, metacell_id=meta.iloc[metacell_index].metacell_id, cell_id=local_cells.iloc[row_index].cell_id, source_index=original))
        pd.DataFrame(membership).to_csv(args.output / f"{name}_metacell_membership.tsv", sep="\t", index=False)
        paired_aggregate = None
        if rna is not None:
            paired_aggregate = sparse.vstack([sparse.csr_matrix(rna[idx[group]].sum(axis=0)) for group in groups], format="csr")
        comparison = fit_group(aggregate, meta, peaks, genes, args.output / name, settings, seed, paired_aggregate, rna_genes)
        summary_path = args.output / name / "method_execution_summary.tsv"
        qualified = False
        if summary_path.exists():
            summary = pd.read_csv(summary_path, sep="\t")
            expected_methods = {"distance", "directnet_adapted", "cicero"} | set(settings.get("external_stage_links", {}))
            expected = summary[summary.method.isin(expected_methods)]
            qualified = bool(len(expected) == len(expected_methods) and expected.execution_status.eq("completed").all()
                             and expected.n_scored_outer_folds.eq(int(settings.get("outer_folds", 5))).all())
        execution_qualifications.append(dict(group=name, all_methods_all_outer_folds=qualified))
        if len(comparison):
            comparison["cohort"], comparison["lineage"] = cohort, lineage
            comparison["lineage_state"] = ";".join(sorted(local_cells.lineage_state.unique())) if "lineage_state" in local_cells else "source_label"
            comparison["source_lineage"] = ";".join(sorted(local_cells.source_lineage.unique())) if "source_lineage" in local_cells else lineage
            all_comparisons.append(comparison)
    if all_comparisons:
        comparisons = pd.concat(all_comparisons, ignore_index=True)
        comparisons["holm_pvalue"] = holm(comparisons.pvalue.to_numpy())
        comparisons["provisional_predictive_improvement"] = (comparisons.ci_low > 0) & (comparisons.holm_pvalue < 0.05)
        comparisons["predictive_improvement"] = False
        comparisons["family_complete"] = False
        comparisons["family"] = "provisional_within_invocation_all_planned_contrasts_global_pooling_required"
        comparisons.to_csv(args.output / "comparisons.tsv", sep="\t", index=False)
    qualified = bool(execution_qualifications and all(x["all_methods_all_outer_folds"] for x in execution_qualifications))
    if settings.get("selected_outer_folds") is not None:
        shard_markers = [json.loads(path.read_text()) for path in args.output.glob("*/shard_completion.json")]
        shard_ok = bool(shard_markers and all(marker.get("status") == "completed" for marker in shard_markers))
        (args.output / "shard_completion.json").write_text(json.dumps(dict(
            status="completed" if shard_ok else "failed", selected_outer_folds=sorted(settings["selected_outer_folds"]),
            all_outer_folds=False, scientific_adoption=False), indent=2) + "\n")
        qualified = False
    (args.output / "completion.json").write_text(json.dumps({"status": "completed" if qualified else "partial", "method_fold_qualification": execution_qualifications, "interpretation": "technical_pilot_no_research_calls" if technical_pilot else "bounded_gene_development" if int(settings.get("max_genes", 0)) > 0 else "candidate_accessibility_context_and_internal_prediction", "pilot": technical_pilot, "rna_available": rna is not None, "scientific_adoption": False}, indent=2) + "\n")


if __name__ == "__main__":
    main()
