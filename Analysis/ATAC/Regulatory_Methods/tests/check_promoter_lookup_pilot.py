#!/usr/bin/env python3
"""Verify exact original CSR arrays and promoter sums on bounded real inputs."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
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


def same_subset(matrix, lookup, columns):
    original, indexed = matrix[:, columns], lookup.columns(columns)
    for name in ["indptr", "indices", "data"]:
        left, right = getattr(original, name), getattr(indexed, name)
        if left.dtype != right.dtype or not np.array_equal(left, right):
            raise ValueError(f"Original CSR subset {name} differs for columns {columns}")
    left = np.asarray(original.sum(axis=1)).ravel()
    right = np.asarray(indexed.sum(axis=1)).ravel()
    if not np.array_equal(left, right):
        raise ValueError("Original CSR row sums differ")
    if not np.array_equal(np.log1p(left), lookup.promoter(columns)):
        raise ValueError("Original promoter targets differ")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID") or "login" in socket.gethostname():
        raise RuntimeError("Actual-input pilot requires a compute allocation")
    if args.output.exists():
        raise FileExistsError(args.output)
    root = Path(__file__).resolve().parents[1]
    source = root / "scripts/promoter_lookup.py"
    module = load("promoter_lookup_candidate", source)
    core = load("promoter_original_core", args.prepared / "source_snapshot/network_compare.py")
    threads = threadpool_limits(limits=1)
    fixtures = [sparse.csr_matrix((np.array([1e16, 1., -1e16, 0., 2., 3.]),
                                  np.array([3, 1, 2, 0, 2, 1], dtype=np.int32),
                                  np.array([0, 4, 4, 6], dtype=np.int32)), shape=(3, 5)),
                sparse.csr_matrix((np.array([1., 2., 0., 4.]),
                                   np.array([2, 2, 1, 0], dtype=np.int32),
                                   np.array([0, 3, 4], dtype=np.int32)), shape=(2, 5)),
                sparse.csr_matrix((3, 5), dtype=float)]
    synthetic_checks = 0
    for matrix in fixtures:
        lookup = module.PromoterLookup(matrix)
        for columns in [np.array([], dtype=int), np.array([0]), np.array([1, 2, 3]), np.array([4, 1, 0])]:
            same_subset(matrix, lookup, columns)
            synthetic_checks += 1
    start = time.perf_counter()
    counts = sparse.load_npz(args.prepared / "counts.npz")
    x = core.normalize_sparse(counts)
    load_seconds = time.perf_counter() - start
    del counts
    peaks = core.read_table(args.prepared / "peaks.tsv")
    peaks[["start0", "end"]] = peaks[["start0", "end"]].astype(int)
    genes = core.read_table(args.prepared / "genes.tsv")
    genes.tss0 = genes.tss0.astype(int)
    rna_genes = set(core.read_table(args.prepared / "rna_genes.tsv").gene_id)
    geometry = core.PeakGeometry(peaks)
    eligible = []
    for gene in genes.to_dict("records"):
        if gene["gene_id"] in rna_genes:
            promoter, _ = geometry.query(gene, 250000, 2000)
            if len(promoter):
                eligible.append((gene, promoter))
    # Even annotation spacing plus promoter-size extrema covers both strands,
    # chromosomes, and small/large selections without inspecting RNA outcomes.
    chosen = sorted(set(np.linspace(0, len(eligible)-1, 96).astype(int).tolist() +
                        [int(np.argmin([len(value[1]) for value in eligible])),
                         int(np.argmax([len(value[1]) for value in eligible]))]))
    sample = [eligible[index] for index in chosen]
    start = time.perf_counter()
    lookup = module.PromoterLookup(x)
    build_seconds = time.perf_counter() - start
    details = []
    for gene, promoter in sample:
        same_subset(x, lookup, promoter)
        details.append(dict(gene_id=gene["gene_id"], chrom=gene["chrom"], strand=gene["strand"], n_promoter_peaks=len(promoter)))
    timings = []
    for repeat in range(3):
        elapsed, results = {}, {}
        for mode in (["original", "indexed"] if repeat % 2 == 0 else ["indexed", "original"]):
            start = time.perf_counter()
            if mode == "original":
                results[mode] = [np.log1p(np.asarray(x[:, promoter].sum(axis=1)).ravel()) for _, promoter in sample]
            else:
                results[mode] = [lookup.promoter(promoter) for _, promoter in sample]
            elapsed[mode] = time.perf_counter() - start
        if not all(np.array_equal(a, b) for a, b in zip(results["original"], results["indexed"])):
            raise ValueError("Timed promoter arrays differ")
        elapsed["speedup"] = elapsed["original"] / elapsed["indexed"]
        timings.append(elapsed)
    args.output.mkdir(parents=True)
    result = dict(status="completed", production_used=False, matrix_shape=list(x.shape), matrix_nnz=x.nnz,
                  original_sorted_indices=bool(x.has_sorted_indices), n_real_genes=len(sample), n_synthetic_checks=synthetic_checks,
                  bit_identical_csr_indptr_indices_data=True, bit_identical_row_sums=True, bit_identical_promoter_targets=True,
                  source_prepared=str(args.prepared.resolve()), helper_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  source_load_seconds=load_seconds, index_build_seconds=build_seconds, index_bytes=lookup.index_bytes,
                  maximum_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  timings_seconds=timings, median_speedup=float(np.median([row["speedup"] for row in timings])),
                  strands=sorted({row["strand"] for row in details}),
                  promoter_peak_range=[min(row["n_promoter_peaks"] for row in details), max(row["n_promoter_peaks"] for row in details)],
                  genes=details, threads=1, outcome_interpretation=False)
    (args.output / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
