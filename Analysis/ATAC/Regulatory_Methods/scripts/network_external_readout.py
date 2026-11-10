#!/usr/bin/env python3
"""One-gene native-method readout against a verified completed balanced baseline.

Rebuild donor-confined aggregates, verify all cached outer-training ATAC inputs,
and reproduce distance/DIRECT predictions before reusing any Cicero predictions.
This analysis is exploratory and cannot establish multimodal improvement.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def completed(root):
    root = Path(root)
    if any(root.glob("INVALIDATED*.json")):
        raise ValueError(f"Invalidated baseline or native run: {root}")
    state = json.loads((root / "completion.json").read_text())
    if state.get("status") != "completed":
        raise ValueError(f"Completed baseline or native run required: {root}")
    return state


def frozen_core(baseline):
    inventory = json.loads((baseline / "source_sha256.json").read_text())
    for source, expected in inventory.items():
        snapshot = baseline / "source_snapshot" / Path(source).name
        if snapshot.exists() and sha256(snapshot) != expected:
            raise ValueError(f"Baseline source snapshot changed: {snapshot}")
    source = baseline / "source_snapshot/network_compare.py"
    if not any(Path(key).name == source.name and value == sha256(source) for key, value in inventory.items()):
        raise ValueError("Frozen core is missing its recorded source hash")
    spec = importlib.util.spec_from_file_location("verified_native_baseline_core", source)
    core = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(core)
    core.require_compute()
    return core, inventory


def same_frame(actual, expected, columns, label):
    if not actual[columns].astype(str).reset_index(drop=True).equals(expected[columns].astype(str).reset_index(drop=True)):
        raise ValueError(f"{label} differs from the completed baseline")


def verify_native_inputs(baseline, roots, paths):
    paths = {Path(path).resolve() for path in paths}
    baseline_time = (baseline / "run_settings.json").stat().st_mtime_ns
    for path in paths:
        if path.stat().st_mtime_ns > baseline_time:
            raise ValueError(f"Raw input changed after baseline initialization: {path}")
    observed = {}
    checks = []
    for root in sorted({Path(root).resolve() for root in roots}):
        manifest = json.loads((root / "input_sha256.json").read_text())
        expected = {Path(path).resolve(): value for path, value in manifest.items()}
        if not paths.issubset(expected):
            raise ValueError("Native source manifest does not cover all common-readout raw inputs")
        for path, digest in expected.items():
            if path not in observed:
                observed[path] = sha256(path)
            if observed[path] != digest:
                raise ValueError(f"Native source input hash changed: {path}")
        checks.append(dict(native_root=str(root), n_verified_files=len(expected),
                           manifest_sha256=sha256(root / "input_sha256.json")))
    return checks


def anchor_predictions(actual, expected, gene):
    methods = {"intercept", "promoter_only", "distance", "directnet_adapted"}
    keys = ["outer_fold", "gene_id", "endpoint", "method", "donor_id"]
    actual = actual[(actual.gene_id == gene) & actual.method.isin(methods)].sort_values(keys).reset_index(drop=True)
    expected = expected[(expected.gene_id == gene) & expected.method.isin(methods)].sort_values(keys).reset_index(drop=True)
    if actual.empty or set(actual.method) != methods or set(expected.method) != methods or actual.duplicated(keys).any() or expected.duplicated(keys).any():
        raise ValueError("Unique nonempty anchor predictions are required")
    same_frame(actual, expected, keys + ["n_metacells", "input_representation", "target_unit"], "Anchor prediction keys/units")
    columns = ["observed", "predicted", "squared_error", "standardized_squared_error", "ridge_alpha"]
    left, right = actual[columns].to_numpy(float), expected[columns].to_numpy(float)
    if not np.allclose(left, right, rtol=1e-10, atol=1e-12, equal_nan=False):
        raise ValueError("Distance/DIRECT/baseline numerical predictions do not reproduce the completed baseline")
    return dict(rows=len(actual), maximum_absolute_difference=float(np.max(np.abs(left - right))))


def verify_atac_caches(core, group_root, aggregate, metadata, peaks, folds, calibration_samples):
    checks = []
    outer = folds.outer_fold.to_numpy(int)
    for fold in sorted(set(outer)):
        train = np.flatnonzero(outer != fold)
        stage = group_root / f"cicero_fold{fold}_outer"
        if any(stage.glob("INVALIDATED*.json")):
            raise ValueError(f"Invalidated Cicero stage: {stage}")
        method = (stage / "links.tsv.method.txt").read_text()
        if "covariance=donor_equal_reliability_weighted" not in method or f"distance_samples={calibration_samples}\n" not in method:
            raise ValueError("Cached Cicero weighting/calibration settings differ")
        nonzero = np.flatnonzero(np.asarray(aggregate[train].sum(axis=0)).ravel() > 0)
        expected_peaks = peaks.iloc[nonzero]
        saved_peaks = core.read_table(stage / "peaks.tsv")
        same_frame(saved_peaks, expected_peaks, ["peak_id", "chrom", "start0", "end"], "Cached training peak order/coordinates")
        saved_counts = core.read_matrix(stage, ["counts.mtx"])
        current = aggregate[train][:, nonzero].T.tocsr()
        if saved_counts.shape != current.shape or (saved_counts != current).nnz:
            raise ValueError(f"Cached ATAC aggregate counts differ: {stage}")
        saved_donors = core.read_table(stage / "donors.tsv")
        expected_donors = pd.DataFrame(dict(metacell_id=[f"metacell_{i+1}" for i in range(len(train))],
                                            donor_id=metadata.iloc[train].donor_id.to_numpy()))
        same_frame(saved_donors, expected_donors, ["metacell_id", "donor_id"], "Cached donor order")
        weights = core.read_table(stage / "links.tsv.donor_weights.tsv")
        same_frame(weights, expected_donors, ["metacell_id", "donor_id"], "Cached covariance donor order")
        if not np.allclose(weights.covariance_weight.astype(float), .5, rtol=0, atol=1e-14):
            raise ValueError("Native secondary comparison requires two equally weighted aggregates per donor")
        checks.append(dict(outer_fold=int(fold), n_training_donors=metadata.iloc[train].donor_id.nunique(),
                           n_training_metacells=len(train), n_training_peaks=len(nonzero),
                           exact_count_match=True, counts_sha256=sha256(stage / "counts.mtx")))
    return checks


def reusable_predictions(baseline, gene, methods, nfolds):
    """Require each requested method on the exact distance donor/endpoint keys."""
    keys = ["outer_fold", "gene_id", "endpoint", "donor_id"]
    anchor = baseline[(baseline.gene_id == gene) & (baseline.method == "distance")]
    expected = anchor[keys].sort_values(keys).reset_index(drop=True)
    if expected.empty or expected.duplicated(keys).any() or anchor.outer_fold.nunique() != nfolds:
        raise ValueError("Baseline distance coverage is incomplete for requested reuse")
    rows = []
    for method in methods:
        selected = baseline[(baseline.gene_id == gene) & (baseline.method == method)]
        actual = selected[keys].sort_values(keys).reset_index(drop=True)
        if not actual.equals(expected) or selected[["squared_error", "standardized_squared_error"]].isna().any().any():
            raise ValueError(f"Missing or incomplete requested reusable method: {method}")
        rows.append(selected)
    return pd.concat(rows, ignore_index=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-root", required=True, type=Path)
    parser.add_argument("--external-config", required=True, type=Path)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--rna-root", required=True, type=Path)
    parser.add_argument("--genes", required=True, type=Path)
    parser.add_argument("--gene", default="ENSG00000000419")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    completed(args.baseline_root)
    core, source_inventory = frozen_core(args.baseline_root)
    settings = json.loads((args.baseline_root / "run_settings.json").read_text())
    seed = int(settings["seed"])
    group_paths = list(args.baseline_root.glob("*/fold_assignments.tsv"))
    if len(group_paths) != 1:
        raise ValueError("Exactly one completed baseline cohort/lineage is required")
    group_root = group_paths[0].parent
    expected_dependencies = json.loads((group_root / "dependencies.json").read_text())
    if core.inventory() != expected_dependencies:
        raise ValueError("Readout environment differs from the completed baseline")
    annotation_hashes = [value for key, value in source_inventory.items() if Path(key).name == args.genes.name]
    if annotation_hashes != [sha256(args.genes)]:
        raise ValueError("Gene annotation differs from the completed baseline")
    extension = json.loads(args.external_config.read_text())
    extension = extension.get("network", extension)
    for key in ["seed", "max_cells_per_donor", "metacell_size", "representation", "outer_folds", "inner_folds", "ridge_alphas", "cis_window", "promoter_flank", "max_distal_peaks"]:
        if key in extension and extension[key] != settings.get(key):
            raise ValueError(f"Native extension changes baseline setting: {key}")
    external = extension.get("external_stage_links", {})
    reused_methods = list(dict.fromkeys(["cicero", *extension.get("reuse_baseline_methods", [])]))
    if set(reused_methods) & set(external):
        raise ValueError("A method cannot be both freshly evaluated and reused")
    saved_external = settings.get("external_stage_links", {})
    saved_controls = settings.get("external_method_controls", {})
    for method in reused_methods:
        if method != "cicero":
            if method not in saved_external:
                raise ValueError(f"Requested reusable method was not configured in baseline: {method}")
            completed(Path(settings["external_stage_roots"][method]))
            for path in saved_external[method].values():
                completed(Path(path).parent)
    if not external:
        raise ValueError("Explicit native stage-link maps are required")
    for method in external:
        completed(Path(extension["external_stage_roots"][method]))
    atac_file = next((args.input_root / name for name in ["cell_counts.npz", "matrix.npz", "ATAC_counts.mtx.gz", "ATAC_counts.mtx"] if (args.input_root / name).exists()), None)
    rna_file = next((args.rna_root / name for name in ["rna_counts.npz", "matrix.npz", "RNA_counts.mtx.gz", "RNA_counts.mtx"] if (args.rna_root / name).exists()), None)
    if atac_file is None or rna_file is None:
        raise ValueError("Balanced native replay requires the explicit paired ATAC/RNA matrix files")
    native_input_checks = verify_native_inputs(args.baseline_root,
        [extension["external_stage_roots"][method] for method in external],
        [atac_file, rna_file, args.input_root / "cells.tsv", args.input_root / "peaks.tsv",
         args.rna_root / "cells.tsv", args.rna_root / "genes.tsv", args.genes])
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "external_stage_config.json").write_text(args.external_config.read_text())
    counts, cells, peaks = core.load_input(args.input_root, seed, int(settings["max_cells_per_donor"]))
    if len(cells[["cohort", "lineage"]].drop_duplicates()) != 1:
        raise ValueError("Exactly one paired cohort/lineage is required")
    keep, excluded = core.condition_eligibility(cells, settings)
    selected = np.flatnonzero(keep)
    local_cells = cells.iloc[selected].reset_index(drop=True)
    if settings["representation"] == "binary":
        counts.eliminate_zeros()
        counts.data[:] = 1
    aggregate, metadata, groups, eligibility = core.build_metacells(counts[selected], local_cells,
        seed, int(settings.get("metacell_size", 50)), int(settings["max_cells_per_donor"]))
    del counts
    if not metadata.groupby("donor_id").size().eq(2).all() or not metadata.n_cells.eq(50).all():
        raise ValueError("Balanced native readout requires exactly two 50-cell aggregates per donor")
    membership = []
    for i, group in enumerate(groups):
        for index in group:
            original = int(local_cells.iloc[index].original_source_index) if "original_source_index" in local_cells else int(selected[index])
            membership.append(dict(cohort=local_cells.iloc[index].cohort, lineage=local_cells.iloc[index].lineage,
                donor_id=local_cells.iloc[index].donor_id, metacell_id=metadata.iloc[i].metacell_id,
                cell_id=local_cells.iloc[index].cell_id, source_index=original))
    membership = pd.DataFrame(membership)
    saved_membership = core.read_table(args.baseline_root / f"{group_root.name}_metacell_membership.tsv")
    same_frame(membership, saved_membership, membership.columns.tolist(), "Raw cell aggregate membership")
    folds = core.read_table(group_paths[0])
    same_frame(metadata, folds, metadata.columns.tolist(), "Aggregate metadata/order")
    calibration_samples = int(settings["cicero_distance_samples"])
    cache_checks = verify_atac_caches(core, group_root, aggregate, metadata, peaks, folds, calibration_samples)
    rna_cells = core.read_table(args.rna_root / "cells.tsv")
    if not core.cell_keys(cells).equals(core.cell_keys(rna_cells)):
        raise ValueError("RNA/ATAC raw cell axes differ")
    rna = core.read_matrix(args.rna_root, ["rna_counts.npz", "matrix.npz", "RNA_counts.mtx.gz", "RNA_counts.mtx"])
    rna_genes = core.read_table(args.rna_root / "genes.tsv")
    if rna.shape != (len(rna_cells), len(rna_genes)) or rna_genes.gene_id.duplicated().any():
        raise ValueError("RNA matrix/gene axes are invalid")
    paired = sparse.vstack([sparse.csr_matrix(rna[selected[group]].sum(axis=0)) for group in groups], format="csr")
    del rna
    genes = core.read_table(args.genes)
    genes.tss0 = genes.tss0.astype(np.int64)
    genes = genes[genes.gene_id == args.gene]
    if len(genes) != 1:
        raise ValueError("Exactly one prespecified annotation gene is required")
    outer_map = folds.groupby("donor_id").outer_fold.agg(lambda values: set(values))
    if not outer_map.map(len).eq(1).all():
        raise ValueError("Baseline donor appears in multiple outer folds")
    # No Cicero refit. Its rows are admitted only after exact aggregate checks
    # above and independent numerical reconstruction of the anchor readouts.
    settings.update(max_genes=0, cicero_rscript="__verified_cached_Cicero_no_refit__",
        outer_fold_map={str(donor): int(next(iter(values))) for donor, values in outer_map.items()},
        external_stage_links=external, external_stage_roots=extension["external_stage_roots"],
        external_method_controls=extension.get("external_method_controls", {}))
    core.fit_group(aggregate, metadata, peaks, genes, args.output / "fresh_readout", settings, seed, paired, rna_genes)
    fresh = pd.read_csv(args.output / "fresh_readout/predictions.tsv", sep="\t", dtype={"donor_id": str})
    baseline = pd.read_csv(group_root / "predictions.tsv", sep="\t", dtype={"donor_id": str})
    anchors = anchor_predictions(fresh, baseline, args.gene)
    reused = reusable_predictions(baseline, args.gene, reused_methods, int(settings["outer_folds"]))
    merged = pd.concat([fresh, reused], ignore_index=True)
    keys = ["outer_fold", "gene_id", "endpoint", "method", "donor_id"]
    if merged.duplicated(keys).any():
        raise ValueError("Duplicate merged prediction rows")
    merged.to_csv(args.output / "predictions.tsv", sep="\t", index=False)
    endpoints = ["promoter_accessibility", "rna"]
    native_methods = list(external) + [method for method in reused_methods if method != "cicero"]
    accumulator = core.ComparisonAccumulator(endpoints, native_methods)
    controls = dict(settings["external_method_controls"])
    controls.update({method: control for method, control in saved_controls.items()
                     if method in reused_methods and control in reused_methods})
    comparison_pairs = list(controls.items()) + [tuple(pair) for pair in extension.get("external_comparison_pairs", [])]
    for candidate, control in comparison_pairs:
        if candidate not in native_methods or control not in native_methods:
            raise ValueError("Both native control methods must be configured")
        for endpoint in endpoints:
            accumulator.statistics[(endpoint, candidate, control)] = {"donors": {}, "genes": set(), "rows": 0}
    for _, frame in merged.groupby(["outer_fold", "gene_id"], sort=True):
        accumulator.append(frame.to_dict("records"))
    comparisons, donor_errors = accumulator.finalize(seed, int(settings["bootstrap_reps"]))
    comparisons["holm_pvalue"] = core.holm(comparisons.pvalue.to_numpy())
    comparisons["family"] = "exploratory_one_gene_native_endpoint_baseline_control_contrasts"
    comparisons["family_complete"] = False
    comparisons["predictive_improvement"] = False
    comparisons["evaluation_status"] = "bounded_one_gene_native_comparison_no_multimodal_improvement_claim"
    comparisons.to_csv(args.output / "exploratory_comparisons.tsv", sep="\t", index=False)
    contextual = comparisons.assign(cohort=metadata.cohort.iloc[0], lineage=metadata.lineage.iloc[0],
        lineage_state=metadata.lineage_state.iloc[0] if "lineage_state" in metadata else "unrecorded",
        gene_id=args.gene)
    contextual.to_csv(args.output / "exploratory_comparisons_context.tsv", sep="\t", index=False)
    donor_errors.to_csv(args.output / "donor_errors.tsv", sep="\t", index=False)
    receipt = dict(status="passed", baseline=str(args.baseline_root.resolve()), gene_id=args.gene,
        source_hashes=source_inventory, dependency_inventory=expected_dependencies,
        exact_cell_membership=True, exact_outer_folds=True, cached_training_ATAC=cache_checks,
        anchor_prediction_check=anchors, calibration_samples=calibration_samples,
        native_stage_config_sha256=sha256(args.external_config), RNA_units="log1p_counts_per_10000_common_readout",
        native_input_hash_checks=native_input_checks,
        reused_methods=reused_methods,
        reused_prediction_sha256=sha256(group_root / "predictions.tsv"),
        baseline_settings_sha256=sha256(args.baseline_root / "run_settings.json"),
        baseline_folds_sha256=sha256(group_root / "fold_assignments.tsv"),
        explicit_native_comparison_pairs=comparison_pairs,
        native_training_units="source_owned_native_stage_manifests", n_donors=metadata.donor_id.nunique(),
        n_metacells=len(metadata), scientific_adoption=False)
    (args.output / "compatibility_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    qualification = []
    for method in native_methods:
        rows = merged[merged.method == method]
        qualification.append(dict(method=method, n_genes=rows.gene_id.nunique(), n_donors=rows.donor_id.nunique(),
            n_outer_folds=rows.outer_fold.nunique(), n_rows=len(rows)))
    pd.DataFrame(qualification).to_csv(args.output / "native_coverage.tsv", sep="\t", index=False)
    qualified = all(item["n_outer_folds"] == int(settings["outer_folds"]) for item in qualification)
    (args.output / "completion.json").write_text(json.dumps(dict(status="completed" if qualified else "partial",
        analysis_scope="bounded_one_gene_native_comparison", scientific_adoption=False,
        matched_native_RNA_only_ablation="not_run", matched_native_prior_only_ablation="not_run",
        proxy_reference_controls={method: "completed" if method in native_methods else "not_run"
                                  for method in extension.get("proxy_reference_methods", [])},
        primary_improvement_claim=False), indent=2) + "\n")


if __name__ == "__main__":
    main()
