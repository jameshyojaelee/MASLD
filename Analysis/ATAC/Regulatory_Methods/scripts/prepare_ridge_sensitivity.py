#!/usr/bin/env python3
"""Reconstruct and verify paired aggregates for a separate ridge sensitivity."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from scipy import sparse


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_core(baseline):
    if any(baseline.glob("INVALIDATED*.json")):
        raise ValueError("Invalidated baseline")
    state = json.loads((baseline / "completion.json").read_text())
    if state.get("status") != "completed" or state.get("pilot", True):
        raise ValueError("Completed non-pilot baseline required")
    inventory = json.loads((baseline / "source_sha256.json").read_text())
    source = baseline / "source_snapshot/network_compare.py"
    expected = [value for key, value in inventory.items() if Path(key).name == source.name]
    if expected != [sha256(source)]:
        raise ValueError("Baseline core snapshot hash differs")
    for path, digest in inventory.items():
        snapshot = baseline / "source_snapshot" / Path(path).name
        if snapshot.exists() and sha256(snapshot) != digest:
            raise ValueError(f"Baseline source snapshot differs: {snapshot}")
    spec = importlib.util.spec_from_file_location("ridge_baseline_core", source)
    core = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(core)
    core.require_compute()
    return core, inventory


def same_frame(actual, expected, columns, label):
    if not actual[columns].astype(str).reset_index(drop=True).equals(
            expected[columns].astype(str).reset_index(drop=True)):
        raise ValueError(f"{label} differs from baseline")


def verify_contract(core, contract, counts, metadata, peaks, genes, rna, rna_genes, folds):
    objects = dict(counts_sha256=(core.matrix_digest, counts),
                   metadata_sha256=(core.table_digest, metadata),
                   peaks_sha256=(core.table_digest, peaks),
                   genes_sha256=(core.table_digest, genes),
                   rna_sha256=(core.matrix_digest, rna),
                   rna_genes_sha256=(core.table_digest, rna_genes),
                   fold_assignments_sha256=(core.table_digest, folds))
    observed = {key: function(value) for key, (function, value) in objects.items()}
    for key, value in observed.items():
        if contract.get(key) != value:
            raise ValueError(f"Original aggregate contract differs: {key}")
    return observed


def source_files(root, rna=False):
    names = (["rna_counts.npz", "matrix.npz", "RNA_counts.mtx.gz", "RNA_counts.mtx"]
             if rna else ["cell_counts.npz", "matrix.npz", "ATAC_counts.mtx.gz", "ATAC_counts.mtx"])
    matrix = next((root / name for name in names if (root / name).exists()), None)
    if matrix is None:
        raise ValueError("Explicit paired matrices required for ridge preparation")
    result = [matrix, root / "cells.tsv", root / ("genes.tsv" if rna else "peaks.tsv")]
    for name in ["paired_export.json", "input_sha256.json", "diagnostics.json"]:
        if (root / name).exists():
            result.append(root / name)
    return result


def prepare(baseline, input_root, rna_root, genes_path, output):
    baseline, input_root, rna_root, genes_path, output = map(
        lambda path: Path(path).resolve(), [baseline, input_root, rna_root, genes_path, output])
    if output.exists():
        raise FileExistsError(f"Refusing overwrite: {output}")
    if output.is_relative_to(baseline):
        raise ValueError("Sensitivity output must be separate from original baseline results")
    core, source_inventory = load_core(baseline)
    group_paths = list(baseline.glob("*/fold_assignments.tsv"))
    if len(group_paths) != 1:
        raise ValueError("Exactly one completed paired baseline group is required")
    group_root = group_paths[0].parent
    contract = json.loads((group_root / "fold_contract.json").read_text())
    settings = json.loads((baseline / "run_settings.json").read_text())
    if settings.get("source_pilot") or settings.get("pilot") or settings.get("max_genes", 0):
        raise ValueError("Genome-wide non-pilot baseline required")
    if core.inventory() != json.loads((group_root / "dependencies.json").read_text()):
        raise ValueError("Current dependency inventory differs from baseline")
    annotation = [value for path, value in source_inventory.items() if Path(path).name == genes_path.name]
    if annotation != [sha256(genes_path)]:
        raise ValueError("Gene annotation differs from baseline")
    seed = int(settings["seed"])
    counts, cells, peaks = core.load_input(input_root, seed, int(settings["max_cells_per_donor"]))
    populations = cells[["cohort", "lineage"]].drop_duplicates()
    if len(populations) != 1 or set(cells.cohort) != {"GSE296875"}:
        raise ValueError("Exactly one genuine paired GSE296875 population is required")
    if "_".join(populations.iloc[0]) != group_root.name:
        raise ValueError("Input population differs from baseline")
    rna_cells = core.read_table(rna_root / "cells.tsv")
    if not core.cell_keys(cells).equals(core.cell_keys(rna_cells)):
        raise ValueError("RNA/ATAC cell identities or order differ")
    rna = core.read_matrix(rna_root, ["rna_counts.npz", "matrix.npz", "RNA_counts.mtx.gz", "RNA_counts.mtx"])
    rna_genes = core.read_table(rna_root / "genes.tsv")
    if rna.shape != (len(cells), len(rna_genes)) or rna_genes.gene_id.duplicated().any():
        raise ValueError("Invalid RNA axes")
    if not np.isfinite(rna.data).all() or (rna.data < 0).any():
        raise ValueError("RNA counts must be finite and nonnegative")
    keep, excluded = core.condition_eligibility(cells, settings)
    indexes = np.flatnonzero(keep)
    local_cells = cells.iloc[indexes].reset_index(drop=True)
    local_counts = counts[indexes]
    if settings["representation"] == "binary":
        local_counts = local_counts.copy()
        local_counts.eliminate_zeros()
        local_counts.data[:] = 1
    elif settings["representation"] != "fragment":
        raise ValueError("Unknown baseline representation")
    aggregate, metadata, groups, eligibility = core.build_metacells(
        local_counts, local_cells, seed, int(settings["metacell_size"]), int(settings["max_cells_per_donor"]))
    paired = sparse.vstack([sparse.csr_matrix(rna[indexes[group]].sum(axis=0)) for group in groups], format="csr")
    membership = []
    for metacell_index, group in enumerate(groups):
        for row_index in group:
            row = local_cells.iloc[row_index]
            original = int(row.original_source_index) if "original_source_index" in local_cells else int(indexes[row_index])
            membership.append(dict(cohort=row.cohort, lineage=row.lineage, donor_id=row.donor_id,
                                   metacell_id=metadata.iloc[metacell_index].metacell_id,
                                   cell_id=row.cell_id, source_index=original))
    membership = pd.DataFrame(membership)
    saved_membership = core.read_table(baseline / f"{group_root.name}_metacell_membership.tsv")
    same_frame(membership, saved_membership, list(membership), "Exact nucleus aggregate membership")
    saved_folds = core.read_table(group_paths[0])
    same_frame(metadata, saved_folds, list(metadata), "Aggregate metadata and order")
    donors = metadata.donor_id.to_numpy()
    folds = metadata.copy()
    folds["outer_fold"] = core.donor_folds(donors, int(settings["outer_folds"]), seed)
    same_frame(folds, saved_folds, list(folds), "Original outer donor folds")
    genes = core.read_table(genes_path)
    genes.tss0 = genes.tss0.astype(np.int64)
    genes = genes.sort_values("gene_id")
    verified = verify_contract(core, contract, aggregate, metadata, peaks, genes, paired, rna_genes, folds)
    # Source files remain read-only; all reconstruction products go to a new root.
    paths = sorted(set(source_files(input_root) + source_files(rna_root, True) + [genes_path]))
    sources = {str(path): dict(size_bytes=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns,
                              sha256=sha256(path)) for path in paths}
    output.mkdir(parents=True, exist_ok=False)
    sparse.save_npz(output / "counts.npz", aggregate)
    sparse.save_npz(output / "rna.npz", paired)
    for name, frame in dict(metadata=metadata, peaks=peaks, genes=genes, rna_genes=rna_genes,
                            fold_assignments=folds, metacell_membership=membership,
                            eligibility=eligibility, condition_exclusions=excluded).items():
        frame.to_csv(output / f"{name}.tsv", sep="\t", index=False)
    snapshot = output / "source_snapshot"
    snapshot.mkdir()
    for path in (baseline / "source_snapshot").iterdir():
        if path.is_file():
            shutil.copyfile(path, snapshot / path.name)
    shutil.copyfile(Path(__file__), snapshot / Path(__file__).name)
    receipt = dict(status="completed", baseline_root=str(baseline), group=group_root.name,
                   input_root=str(input_root), rna_root=str(rna_root), genes_source=str(genes_path),
                   seed=seed, settings=settings, verified_contract=verified, source_files=sources,
                   source_inventory=source_inventory, dependencies=core.inventory(),
                   n_donors=metadata.donor_id.nunique(), n_metacells=len(metadata),
                   n_nuclei=len(membership), scientific_adoption=False)
    (output / "preparation.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({key: receipt[key] for key in ["status", "group", "n_donors", "n_metacells", "n_nuclei"]}), flush=True)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["baseline-root", "input-root", "rna-root", "genes", "output"]:
        parser.add_argument(f"--{name}", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.baseline_root, args.input_root, args.rna_root, args.genes, args.output)


if __name__ == "__main__":
    main()
