#!/usr/bin/env python3
"""RNA-only ridge sensitivity using the original donor folds and training links.

The original grid is replayed on every evaluable gene as a numerical anchor.
Cached Cicero stages are read through new-directory links; original files are
never written. DIRECT-NET feature selection is replayed with its original seed.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import multiprocessing
import os
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits

ORIGINAL_GRID = [0.1, 1., 10., 100.]
WIDE_GRID = [0.1, 1., 10., 100., 1000., 10000., 100000., 1000000., 10000000., 100000000.]
METHODS = ["distance", "directnet_adapted", "cicero"]
KEYS = ["outer_fold", "gene_id", "endpoint", "method", "donor_id"]
CORE = None
CONTEXT = None
THREAD_LIMIT = None


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def same_frame(actual, expected, columns, label):
    if not actual[columns].astype(str).reset_index(drop=True).equals(expected[columns].astype(str).reset_index(drop=True)):
        raise ValueError(f"{label} differs from original")


def load_core(prepared):
    global CORE
    path = prepared / "source_snapshot/network_compare.py"
    receipt = json.loads((prepared / "preparation.json").read_text())
    if receipt.get("status") != "completed":
        raise ValueError("Completed input preparation required before source import")
    expected = [digest for name, digest in receipt["source_inventory"].items()
                if Path(name).name == "network_compare.py"]
    if expected != [sha256(path)]:
        raise ValueError("Prepared core snapshot differs from original source")
    spec = importlib.util.spec_from_file_location("original_ridge_core", path)
    CORE = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(CORE)
    CORE.require_compute()


def init_worker():
    global THREAD_LIMIT
    THREAD_LIMIT = threadpool_limits(limits=1)


def tune(designs, y, donors, train, inner, grid):
    """Preserve original donor loss averaging, alpha order and tie breaking."""
    losses = {alpha: [] for alpha in grid}
    for k in sorted(set(inner)):
        it, iv = train[inner != k], train[inner == k]
        design = designs[k]
        for alpha in grid:
            prediction = CORE.weighted_ridge(design[it], y[it], donors[it], alpha, design[iv])
            values = pd.DataFrame({"donor": donors[iv], "loss": (prediction - y[iv]) ** 2}).groupby("donor").loss.mean()
            losses[alpha].extend(values.tolist())
    return {alpha: float(np.mean(values)) for alpha, values in losses.items()}


def fit_gene(job):
    gene_number, gene = job
    c = CONTEXT
    gene_id = gene["gene_id"]
    rows, selections, tuning = [], [], []
    if gene_id not in c["rna_lookup"]:
        return gene_number, gene_id, rows, selections, tuning
    prom, candidates = c["geometry"].query(gene, c["window"], c["promoter_bp"])
    if not len(prom) or not len(candidates):
        return gene_number, gene_id, rows, selections, tuning
    lookup = c.get("promoter_lookup")
    promoter = (lookup.promoter(prom) if lookup is not None else
                np.log1p(np.asarray(c["x"][:, prom].sum(axis=1)).ravel()))
    train, test, inner, donors = (c[key] for key in ("train", "test", "inner", "donors"))
    if np.var(promoter[train]) <= 1e-12:
        return gene_number, gene_id, rows, selections, tuning
    y = CORE.log_columns(c["rna_normal"], [c["rna_lookup"][gene_id]]).ravel()
    weights = CORE.donor_weights(donors[train])
    variance = float(np.average((y[train] - np.average(y[train], weights=weights)) ** 2, weights=weights))
    if variance <= 1e-12:
        return gene_number, gene_id, rows, selections, tuning
    fold, seed = c["fold"], c["seed"]
    for method in ["promoter_only", *METHODS]:
        if method == "promoter_only":
            design = promoter[:, None]
            designs = {k: design for k in sorted(set(inner))}
        else:
            chosen, scores = CORE.selected_indices(method, candidates, prom, c["peaks"], gene,
                c["training_x"], promoter[train], donors[train], CORE.stable_seed(seed, (fold, gene_id)),
                c["max_links"], c["cicero"]["outer"])
            if not len(chosen):
                continue
            for rank, (idx, score) in enumerate(zip(chosen, scores), 1):
                selections.append(dict(outer_fold=fold, gene_id=gene_id, method=method,
                    rank=rank, peak_id=c["peaks"].iloc[idx].peak_id, score=float(score)))
            design = CORE.design_matrix(c["feature_x"], chosen, promoter, "rna")
            designs = {}
            for k in sorted(set(inner)):
                it = train[inner != k]
                selected, _ = CORE.selected_indices(method, candidates, prom, c["peaks"], gene,
                    c["inner_training_x"][k], promoter[it], donors[it],
                    CORE.stable_seed(seed, (fold, k, gene_id)), c["max_links"], c["cicero"][f"inner{k}"])
                designs[k] = CORE.design_matrix(c["feature_x"], selected, promoter, "rna")
        losses = tune(designs, y, donors, train, inner, WIDE_GRID)
        for grid_name, grid in [("original_replay", ORIGINAL_GRID), ("wide", WIDE_GRID)]:
            best = min(grid, key=lambda alpha: (losses[alpha], alpha))
            predicted = CORE.weighted_ridge(design[train], y[train], donors[train], best, design[test])
            local = []
            CORE.save_predictions(local, fold, gene_id, "rna", method, test, c["metadata"], y, predicted, variance, best)
            for row in local:
                row.update(grid=grid_name, target_unit="log1p_RNA_counts_per_10000", cohort="GSE296875",
                    lineage=c["lineage"], input_representation="binary")
            rows.extend(local)
            tuning.append(dict(outer_fold=fold, gene_id=gene_id, method=method, grid=grid_name,
                ridge_alpha=best, inner_mse=losses[best], lower_boundary_mse=losses[grid[0]],
                upper_boundary_mse=losses[grid[-1]], n_outer_features=design.shape[1]))
    return gene_number, gene_id, rows, selections, tuning


def verify_anchor(actual, expected):
    """Reject changed coverage, transformations, tuning or prediction errors."""
    actual = pd.DataFrame(actual)
    if actual.empty and expected.empty:
        return 0., 0
    if actual.empty or expected.empty:
        raise ValueError("Original replay changed evaluable gene coverage")
    actual = actual[actual.grid.eq("original_replay")].sort_values(KEYS).reset_index(drop=True)
    expected = expected.sort_values(KEYS).reset_index(drop=True)
    if actual.duplicated(KEYS).any() or expected.duplicated(KEYS).any():
        raise ValueError("Duplicate prediction keys")
    same_frame(actual, expected, KEYS + ["n_metacells", "input_representation", "target_unit"], "Original prediction keys")
    if not np.array_equal(actual.ridge_alpha.to_numpy(float), expected.ridge_alpha.to_numpy(float)):
        raise ValueError("Original selected alpha changed")
    columns = ["observed", "predicted", "squared_error", "standardized_squared_error"]
    a, b = actual[columns].to_numpy(float), expected[columns].to_numpy(float)
    if not np.isfinite(a).all() or not np.allclose(a, b, rtol=1e-10, atol=1e-11):
        raise ValueError(f"Original prediction numerical replay failed: {np.max(np.abs(a-b))}")
    return float(np.max(np.abs(a-b))), len(actual)


def cached_cicero(source, output, aggregate, peaks, metadata, fit_idx, calibration):
    if any(source.glob("INVALIDATED*")):
        raise ValueError("Invalidated Cicero source")
    method = (source / "links.tsv.method.txt").read_text()
    timing = json.loads((source / "timing.json").read_text())
    log = (source / "run.log").read_text()
    if re.findall(r"Models with errors:\s*(\d+)", log) != ["0"] or "window_models_complete" not in log or timing.get("directed_edges", 0) <= 0:
        raise ValueError("Original Cicero stage did not complete successfully")
    if "covariance=donor_equal_reliability_weighted\n" not in method or f"distance_samples={calibration}\n" not in method:
        raise ValueError("Original Cicero weighting/calibration changed")
    nonzero = np.flatnonzero(np.asarray(aggregate[fit_idx].sum(axis=0)).ravel() > 0)
    local_peaks = peaks.iloc[nonzero]
    saved_peaks = CORE.read_table(source / "peaks.tsv")
    same_frame(saved_peaks, local_peaks, list(peaks.columns), "Cicero training peaks")
    saved_counts = CORE.read_matrix(source, ["counts.mtx"])
    current = aggregate[fit_idx][:, nonzero].T.tocsr()
    if saved_counts.shape != current.shape or (saved_counts != current).nnz:
        raise ValueError("Cicero training counts changed")
    donors = metadata.iloc[fit_idx].donor_id.to_numpy()
    saved_donors = CORE.read_table(source / "donors.tsv")
    expected = pd.DataFrame(dict(metacell_id=[f"metacell_{i+1}" for i in range(len(donors))], donor_id=donors))
    same_frame(saved_donors, expected, ["metacell_id", "donor_id"], "Cicero donor order")
    weights = CORE.read_table(source / "links.tsv.donor_weights.tsv")
    same_frame(weights, expected, ["metacell_id", "donor_id"], "Cicero covariance donor order")
    n = pd.Series(donors).value_counts()
    if not np.allclose(weights.covariance_weight.to_numpy(float), [1. / n[d] for d in donors], rtol=0, atol=1e-14):
        raise ValueError("Cicero covariance does not give donors equal weight")
    output.mkdir()
    # from_tsv writes a sidecar, so its input path must be in the NEW directory.
    link = output / "links.tsv"
    link.symlink_to((source / "links.tsv").resolve())
    graph = CORE.CiceroGraph.from_tsv(link, local_peaks)
    receipt = dict(source=str(source.resolve()), n_training_donors=len(set(donors)),
        n_training_metacells=len(donors), n_training_peaks=len(local_peaks),
        exact_counts=True, exact_donor_order=True, equal_donor_weight=True,
        links_sha256=sha256(source / "links.tsv"))
    (output / "source_verification.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return graph, receipt


def main():
    global CONTEXT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--fold", required=True, type=int, choices=range(5))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--max-genes", type=int, default=0, help="Technical check only; zero evaluates all genes")
    parser.add_argument("--indexed-promoters", action="store_true",
                        help="Preserve original CSR order using an indexed column lookup")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    load_core(args.prepared)
    workers = CORE.worker_count(args.workers, "readout_workers")
    preparation = json.loads((args.prepared / "preparation.json").read_text())
    if preparation["status"] != "completed":
        raise ValueError("Completed original input reconstruction required")
    baseline = Path(preparation["baseline_root"])
    if args.output.resolve().is_relative_to(baseline.resolve()):
        raise ValueError("Sensitivity must be separate from original results")
    group = preparation["group"]
    lineage = group.removeprefix("GSE296875_")
    settings = json.loads((baseline / "run_settings.json").read_text())
    if settings["ridge_alphas"] != ORIGINAL_GRID:
        raise ValueError("Original ridge grid differs")
    if CORE.inventory() != json.loads((baseline / group / "dependencies.json").read_text()):
        raise ValueError("Original numerical environment changed")
    args.output.mkdir(parents=True)
    started = time.monotonic()
    config = dict(status="running", prepared_root=str(args.prepared.resolve()), baseline_root=str(baseline),
        group=group, lineage=lineage, outer_fold=args.fold, original_grid=ORIGINAL_GRID, wide_grid=WIDE_GRID,
        workers=workers, max_genes=args.max_genes, indexed_promoters=args.indexed_promoters,
        endpoint="rna", scientific_adoption=False,
        primary_family="six_RNA_distal_method_by_population_vs_promoter_only_Holm", seed=settings["seed"])
    (args.output / "run_settings.json").write_text(json.dumps(config, indent=2) + "\n")
    metadata = CORE.read_table(args.prepared / "metadata.tsv")
    metadata.n_cells = metadata.n_cells.astype(int)
    folds = CORE.read_table(args.prepared / "fold_assignments.tsv")
    same_frame(folds, metadata, list(metadata.columns), "Prepared folds/metadata")
    peaks = CORE.read_table(args.prepared / "peaks.tsv")
    peaks[["start0", "end"]] = peaks[["start0", "end"]].astype(np.int64)
    genes = CORE.read_table(args.prepared / "genes.tsv")
    genes.tss0 = genes.tss0.astype(np.int64)
    rna_genes = CORE.read_table(args.prepared / "rna_genes.tsv")
    aggregate = sparse.load_npz(args.prepared / "counts.npz")
    rna = sparse.load_npz(args.prepared / "rna.npz")
    contract = json.loads((baseline / group / "fold_contract.json").read_text())
    for key, digest in {"counts_sha256": CORE.matrix_digest(aggregate), "rna_sha256": CORE.matrix_digest(rna),
                        "metadata_sha256": CORE.table_digest(metadata), "peaks_sha256": CORE.table_digest(peaks),
                        "genes_sha256": CORE.table_digest(genes), "rna_genes_sha256": CORE.table_digest(rna_genes),
                        "fold_assignments_sha256": CORE.table_digest(folds.assign(outer_fold=folds.outer_fold.astype(int)))}.items():
        if digest != contract[key]:
            raise ValueError(f"Prepared input contract changed: {key}")
    seed = int(settings["seed"])
    outer, donors = folds.outer_fold.to_numpy(int), metadata.donor_id.to_numpy()
    if not np.array_equal(outer, CORE.donor_folds(donors, 5, seed)):
        raise ValueError("Original outer folds differ")
    train, test = np.flatnonzero(outer != args.fold), np.flatnonzero(outer == args.fold)
    inner = CORE.donor_folds(donors[train], 3, CORE.stable_seed(seed, ("inner", args.fold)))
    pd.DataFrame(dict(metacell_id=metadata.iloc[train].metacell_id, donor_id=donors[train], inner_fold=inner)).to_csv(args.output / "inner_fold_assignments.tsv", sep="\t", index=False)
    receipt = json.loads((baseline / "merge_receipt.json").read_text())
    source_info = receipt["sources"][str(args.fold)]
    source_group = Path(source_info["root"]) / group
    if source_info["legacy_snapshot"]:
        qualification = json.loads((baseline.parents[2] / "snapshot_qualifications" / f"{group}_fold{args.fold}.json").read_text())
        if (qualification["status"] != "qualified_snapshot" or qualification["outer_fold"] != args.fold
                or Path(qualification["source_root"]).resolve() != source_group.parent.resolve()
                or sha256(source_group.parent / "source_snapshot/network_compare.py") != qualification["legacy_source_sha256"]
                or not (source_group / f"cicero_fold{args.fold+1}_outer/counts.mtx").exists()):
            raise ValueError("Original early-fold qualification differs")
    else:
        marker = json.loads((source_group / f"fold{args.fold}_completed.json").read_text())
        if marker["status"] != "completed_fold" or marker["outer_fold"] != args.fold or marker["failures"]:
            raise ValueError("Original fold contains failures")
    cache, cache_receipts = {}, []
    for tag, indices in [("outer", train)] + [(f"inner{k}", train[inner != k]) for k in sorted(set(inner))]:
        print(f"Verify/load cached Cicero fold={args.fold} stage={tag}", flush=True)
        cache[tag], verified = cached_cicero(source_group / f"cicero_fold{args.fold}_{tag}", args.output / tag,
            aggregate, peaks, metadata, indices, int(settings["cicero_distance_samples"]))
        cache_receipts.append(verified)
    original_parts = []
    for chunk in pd.read_csv(source_group / "predictions.tsv", sep="\t", chunksize=200000, dtype={"donor_id": str}):
        original_parts.append(chunk[(chunk.outer_fold == args.fold) & chunk.endpoint.eq("rna")])
    original = pd.concat(original_parts, ignore_index=True)
    if original.empty or original.duplicated(KEYS).any():
        raise ValueError("Missing or duplicate original RNA predictions")
    original_groups = original.groupby("gene_id", sort=False).indices
    saved_links = pd.read_csv(source_group / "links.tsv", sep="\t")
    saved_links = saved_links[saved_links.outer_fold.eq(args.fold)]
    link_groups = saved_links.groupby("gene_id", sort=False).indices
    x = CORE.normalize_sparse(aggregate)
    promoter_lookup = None
    if args.indexed_promoters:
        from promoter_lookup import PromoterLookup
        lookup_started = time.monotonic()
        promoter_lookup = PromoterLookup(x)
        config.update(promoter_index_bytes=promoter_lookup.index_bytes,
                      promoter_index_seconds=time.monotonic()-lookup_started)
    CONTEXT = dict(x=x, feature_x=x.tocsc(), donors=donors, metadata=metadata, peaks=peaks,
        promoter_lookup=promoter_lookup,
        geometry=CORE.PeakGeometry(peaks), train=train, test=test, inner=inner, fold=args.fold, seed=seed,
        window=settings["cis_window"], promoter_bp=settings["promoter_flank"], max_links=settings["max_distal_peaks"],
        training_x=x[train].tocsc(), inner_training_x={k: x[train[inner != k]].tocsc() for k in sorted(set(inner))},
        cicero=cache, rna_normal=CORE.normalize_sparse(rna).tocsc(),
        rna_lookup={g: i for i, g in enumerate(rna_genes.gene_id)}, lineage=lineage)
    jobs = list(enumerate(genes.to_dict("records"), 1))
    if args.max_genes:
        jobs = [j for j in jobs if j[1]["gene_id"] in original_groups][:args.max_genes]
    columns = KEYS + ["n_metacells", "observed", "predicted", "squared_error", "standardized_squared_error", "ridge_alpha", "grid", "target_unit", "cohort", "lineage", "input_representation"]
    predictions = CORE.TableStream(args.output / "predictions.tsv", columns)
    tune_stream = CORE.TableStream(args.output / "tuning.tsv", ["outer_fold", "gene_id", "method", "grid", "ridge_alpha", "inner_mse", "lower_boundary_mse", "upper_boundary_mse", "n_outer_features"])
    maximum_difference, n_anchors, n_links, n_tested = 0., 0, 0, 0
    readout_started = time.monotonic()
    pool = multiprocessing.get_context("fork").Pool(workers, initializer=init_worker)
    try:
        for start in range(0, len(jobs), 4 * workers):
            for number, gene_id, rows, selected, tuning in pool.imap(fit_gene, jobs[start:start+4*workers], chunksize=4):
                expected = original.iloc[original_groups.get(gene_id, [])]
                difference, anchors = verify_anchor(rows, expected)
                maximum_difference = max(maximum_difference, difference)
                n_anchors += anchors
                if selected:
                    actual_links = pd.DataFrame(selected)
                    expected_links = saved_links.iloc[link_groups.get(gene_id, [])]
                    lkeys = ["outer_fold", "gene_id", "method", "rank", "peak_id"]
                    a = actual_links.sort_values(lkeys).reset_index(drop=True)
                    b = expected_links.sort_values(lkeys).reset_index(drop=True)
                    same_frame(a, b, lkeys, "Outer selected distal peaks")
                    if not np.allclose(a.score, b.score, rtol=1e-10, atol=1e-12):
                        raise ValueError("Replayed outer link scores differ")
                    n_links += len(a)
                n_tested += bool(rows)
                predictions.append(rows)
                tune_stream.append(tuning)
                if number % 100 == 0 or args.max_genes:
                    print(f"Readout fold={args.fold} gene={number}/{len(genes)} evaluated={n_tested} anchors={n_anchors} max_delta={maximum_difference:.3g}", flush=True)
        pool.close()
        pool.join()
    except BaseException:
        pool.terminate()
        pool.join()
        raise
    finally:
        predictions.close()
        tune_stream.close()
    if not args.max_genes and n_anchors != len(original):
        raise ValueError("Original RNA coverage was not fully replayed")
    config.update(status="completed", technical_check=bool(args.max_genes), n_evaluable_genes=n_tested,
        n_original_predictions_verified=n_anchors, n_outer_links_verified=n_links,
        maximum_original_numeric_difference=maximum_difference, n_prediction_rows=predictions.rows,
        readout_elapsed_seconds=time.monotonic()-readout_started,
        elapsed_seconds=time.monotonic()-started, cache_verification=cache_receipts,
        n_donors=len(set(donors)), n_held_donors=len(set(donors[test])), versions=CORE.inventory())
    (args.output / "completion.json").write_text(json.dumps(config, indent=2) + "\n")
    print(json.dumps({key: config[key] for key in ["status", "outer_fold", "n_prediction_rows", "elapsed_seconds", "maximum_original_numeric_difference"]}), flush=True)


if __name__ == "__main__":
    main()
