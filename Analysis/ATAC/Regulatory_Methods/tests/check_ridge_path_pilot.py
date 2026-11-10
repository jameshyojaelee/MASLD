#!/usr/bin/env python3
"""Bounded arithmetic-equivalence and runtime comparison of ridge paths."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import socket
import time

import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def old_tune(core, designs, y, donors, train, inner, grid):
    losses = {alpha: [] for alpha in grid}
    for k in sorted(set(inner)):
        it, iv = train[inner != k], train[inner == k]
        design = designs[k]
        for alpha in grid:
            prediction = core.weighted_ridge(design[it], y[it], donors[it], alpha, design[iv])
            values = pd.DataFrame({"donor": donors[iv], "loss": (prediction - y[iv]) ** 2}).groupby("donor").loss.mean()
            losses[alpha].extend(values.tolist())
    return {alpha: float(np.mean(values)) for alpha, values in losses.items()}


def compare(core, path, inputs, grid):
    left = old_tune(core, *inputs, grid)
    right = path.tune_path(core, *inputs, grid)
    if not np.array_equal(list(left.values()), list(right.values())):
        raise ValueError(f"Ridge path losses differ: {left}, {right}")
    for subset in [grid[:4], grid]:
        a = min(subset, key=lambda alpha: (left[alpha], alpha))
        b = min(subset, key=lambda alpha: (right[alpha], alpha))
        if a != b:
            raise ValueError("Selected alpha differs")
    return max(abs(left[alpha] - right[alpha]) for alpha in grid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID") or "login" in socket.gethostname():
        raise RuntimeError("Runtime pilot requires a compute allocation")
    if args.output.exists():
        raise FileExistsError(args.output)
    root = Path(__file__).resolve().parents[1]
    core = load("path_original_core", args.prepared / "source_snapshot/network_compare.py")
    path = load("path_candidate", root / "scripts/ridge_path.py")
    grid = [.1, 1., 10., 100., 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
    threads = threadpool_limits(limits=1)
    rng = np.random.default_rng(42)
    donors = np.repeat(np.array(["a", "b", "c", "d", "e", "f"]), [2, 3, 1, 4, 2, 3])
    train = np.arange(len(donors))
    mapping = {donor: index % 3 for index, donor in enumerate(sorted(set(donors)))}
    inner = np.array([mapping[donor] for donor in donors])
    y = rng.normal(size=len(donors))
    fixtures = []
    for width in [0, 1, 4, 11]:
        designs = {k: rng.normal(size=(len(donors), max(0, width-k))) for k in range(3)}
        if width:
            for value in designs.values():
                if value.shape[1]:
                    value[:, 0] = 1  # Constant feature exercises the scale floor.
        delta = compare(core, path, (designs, y, donors, train, inner), grid)
        fixtures.append(dict(nominal_features=width, max_loss_delta=delta, uneven_donors=True, constant_features=bool(width)))
    metadata = core.read_table(args.prepared / "metadata.tsv")
    donors = metadata.donor_id.to_numpy()
    folds = core.read_table(args.prepared / "fold_assignments.tsv").outer_fold.to_numpy(int)
    train = np.flatnonzero(folds != 0)
    inner = core.donor_folds(donors[train], 3, core.stable_seed(42, ("inner", 0)))
    x = core.normalize_sparse(sparse.load_npz(args.prepared / "counts.npz"))
    feature_x = x.tocsc()
    rna = core.normalize_sparse(sparse.load_npz(args.prepared / "rna.npz")).tocsc()
    rna_genes = core.read_table(args.prepared / "rna_genes.tsv")
    lookup = {gene: index for index, gene in enumerate(rna_genes.gene_id)}
    genes = core.read_table(args.prepared / "genes.tsv")
    genes.tss0 = genes.tss0.astype(int)
    peaks = core.read_table(args.prepared / "peaks.tsv")
    peaks[["start0", "end"]] = peaks[["start0", "end"]].astype(int)
    geometry = core.PeakGeometry(peaks)
    training_x = {k: x[train[inner != k]].tocsc() for k in sorted(set(inner))}
    jobs, details = [], []
    for gene in genes.to_dict("records"):
        gene_id = gene["gene_id"]
        if gene_id not in lookup:
            continue
        promoter, candidates = geometry.query(gene, 250000, 2000)
        if not len(promoter) or not len(candidates):
            continue
        p = np.log1p(np.asarray(x[:, promoter].sum(axis=1)).ravel())
        target = core.log_columns(rna, [lookup[gene_id]]).ravel()
        if np.var(p[train]) <= 1e-12 or np.var(target[train]) <= 1e-12:
            continue
        cap = len(jobs) % 11
        designs = {}
        for k in sorted(set(inner)):
            it = train[inner != k]
            selected, _ = core.selected_indices("distance", candidates, promoter, peaks, gene,
                training_x[k], p[it], donors[it], core.stable_seed(42, (0, k, gene_id)), 10)
            selected = selected[:max(0, cap-k)]
            designs[k] = core.design_matrix(feature_x, selected, p, "rna")
        if max(design.shape[1] for design in designs.values()) < cap + 1:
            continue
        inputs = designs, target, donors, train, inner
        delta = compare(core, path, inputs, grid)
        jobs.append(inputs)
        details.append(dict(gene_id=gene_id, inner_feature_counts={str(k): value.shape[1] for k, value in designs.items()}, max_loss_delta=delta))
        if len(jobs) == 30:
            break
    if len(jobs) != 30:
        raise ValueError("Insufficient valid real-data pilot genes")
    timings = []
    for repeat in range(3):
        elapsed = {}
        order = ["original", "path"] if repeat % 2 == 0 else ["path", "original"]
        for mode in order:
            start = time.perf_counter()
            for inputs in jobs:
                function = old_tune if mode == "original" else path.tune_path
                function(core, *inputs, grid)
            elapsed[mode] = time.perf_counter() - start
        elapsed["speedup"] = elapsed["original"] / elapsed["path"]
        timings.append(elapsed)
    args.output.mkdir(parents=True)
    result = dict(status="completed", n_real_genes=len(jobs), n_synthetic_fixtures=len(fixtures),
                  bit_identical_losses=True, selected_original_and_wide_alpha_identical=True,
                  observed_feature_widths=sorted({n for row in details for n in row["inner_feature_counts"].values()}),
                  fixtures=fixtures, genes=details, timings_seconds=timings,
                  median_speedup=float(np.median([row["speedup"] for row in timings])),
                  threads=1, source_prepared=str(args.prepared.resolve()), outcome_interpretation=False)
    (args.output / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
