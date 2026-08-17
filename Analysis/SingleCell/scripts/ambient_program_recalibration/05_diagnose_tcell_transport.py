#!/usr/bin/env python3
"""Diagnose the frozen T-cell scoring substrate against complete-atlas counts."""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import hotspot
import numpy as np
import pandas as pd
import scipy.sparse as sp
from hotspot import modules as hotspot_modules

from common import HOTSPOT_ROOT, MEMBERSHIP, REGISTRY, refuse_existing, require


SCRIPT_ROOT = Path(__file__).resolve().parent
CANDIDATE = Path(os.environ["CAND_ROOT"])
OUTPUT = CANDIDATE / "diagnostics/tcell_transport_diagnostic.tsv"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    require(specification is not None and specification.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def correlation(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.corrcoef(np.asarray(left, dtype=float), np.asarray(right, dtype=float))[0, 1])


def main() -> None:
    refuse_existing(OUTPUT)
    rescore = load_module("ambient_rescore", SCRIPT_ROOT / "02_rescore_all117.py")
    producer = rescore.load_score_producer()
    lineage = "tcells"
    run_metadata = __import__("json").loads(
        (HOTSPOT_ROOT / lineage / "run_metadata.json").read_text()
    )
    atlas = producer.load_atlas(lineage, smoke=False)
    excluded = set(run_metadata.get("exclude_datasets") or [])
    if excluded:
        atlas = atlas[~atlas.obs["dataset"].astype(str).isin(excluded)].copy()
    atlas = producer.strip_confounders(atlas)
    atlas = producer.filter_detected(atlas, min_frac=0.01)
    producer.ensure_raw_layer(atlas)
    latent = producer.resolve_latent(atlas)
    hotspot_object = hotspot.Hotspot(
        atlas,
        layer_key="counts",
        model="danb",
        latent_obsm_key=latent,
        umi_counts_obs_key="n_counts",
    )
    hotspot_object.create_knn_graph(
        weighted_graph=False, n_neighbors=int(run_metadata["params"]["n_neighbors"])
    )

    score_genes = (CANDIDATE / "work/score_genes.txt").read_text().splitlines()
    gene_position = {gene: index for index, gene in enumerate(score_genes)}
    cell_index = pd.Index(atlas.obs_names.astype(str))
    raw_counts, corrected_counts, _ = rescore.assemble_lineage_counts(
        lineage, cell_index, score_genes
    )
    registry = pd.read_csv(REGISTRY, sep="\t")
    registry = registry[registry["cell_type"] == lineage].set_index("module")
    membership = pd.read_csv(MEMBERSHIP, sep="\t")
    membership = membership[membership["cell_type"] == lineage]
    stored = pd.read_parquet(HOTSPOT_ROOT / lineage / "cell_scores.parquet")
    stored["module"] = stored["module"].astype(int)
    arguments = (
        hotspot_object.model,
        hotspot_object.umi_counts.values,
        hotspot_object.neighbors.values,
        hotspot_object.weights.values,
    )

    rows = []
    for module in sorted(int(value) for value in registry.index):
        genes = membership.loc[
            membership["module"].astype(int) == module, "source_gene"
        ].astype(str).tolist()
        gene_indices = np.array([gene_position[gene] for gene in genes])
        native = atlas[:, genes].layers["counts"]
        if sp.issparse(native):
            native = native.toarray()
        native = np.asarray(native, dtype=np.float64).T
        complete_raw = np.asarray(raw_counts[gene_indices, :].todense(), dtype=np.float64)
        complete_corrected = np.asarray(
            corrected_counts[gene_indices, :].todense(), dtype=np.float64
        )
        native_score = hotspot_modules.compute_scores(native, *arguments)
        complete_raw_score = hotspot_modules.compute_scores(complete_raw, *arguments)
        complete_corrected_score = hotspot_modules.compute_scores(
            complete_corrected, *arguments
        )
        stored_score = (
            stored[stored["module"] == module]
            .set_index("cell_id")["score"]
            .reindex(cell_index)
            .to_numpy(float)
        )
        require(np.isfinite(stored_score).all(), f"stored score join failed: {module}")
        rows.append(
            {
                "program_uid": registry.loc[module, "program_uid"],
                "module": module,
                "module_name": registry.loc[module, "module_name"],
                "n_cells": len(cell_index),
                "n_genes": len(genes),
                "native_values_all_integer": bool(np.all(native == np.round(native))),
                "native_reconstruction_r": correlation(native_score, stored_score),
                "complete_integer_raw_reconstruction_r": correlation(
                    complete_raw_score, stored_score
                ),
                "native_vs_complete_raw_score_r": correlation(
                    native_score, complete_raw_score
                ),
                "complete_raw_vs_corrected_score_r": correlation(
                    complete_raw_score, complete_corrected_score
                ),
            }
        )
    pd.DataFrame(rows).to_csv(OUTPUT, sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
