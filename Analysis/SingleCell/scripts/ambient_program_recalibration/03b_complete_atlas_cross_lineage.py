#!/usr/bin/env python3
"""Run the fixed two-program localization control on the complete six-lineage atlas."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

from common import OLD_RAW_EXPORT, build_sample_to_donor, refuse_existing, require, sha256


SCRIPT_ROOT = Path(__file__).resolve().parent
AMBIENT_INPUT = Path(os.environ["AMBIENT_CAND_ROOT"])
CROSS_CANDIDATE = Path(os.environ["CROSS_CAND_ROOT"])
COMPLETE_ATLAS = (
    Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
    / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
)


def load_prototype():
    path = SCRIPT_ROOT / "03_cross_lineage.py"
    specification = importlib.util.spec_from_file_location("retired_cross_prototype", path)
    require(
        specification is not None and specification.loader is not None,
        "cannot load cross-lineage prototype",
    )
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


prototype = load_prototype()


def read_raw_chunk(
    directory: Path, chunk: dict[str, int], n_genes: int
) -> sp.csc_matrix:
    chunk_id = int(chunk["chunk"])
    n_cells = int(chunk["ncells"])
    n_nonzero = int(chunk["nnz"])
    indices = np.fromfile(directory / f"chunk{chunk_id}_indices.bin", dtype="<i4")
    indptr = np.fromfile(directory / f"chunk{chunk_id}_indptr.bin", dtype="<i4")
    values = np.fromfile(directory / f"chunk{chunk_id}_data.bin", dtype="<i4")
    require(
        len(indices) == n_nonzero
        and len(values) == n_nonzero
        and len(indptr) == n_cells + 1,
        f"raw chunk shape drift: {directory.name} chunk {chunk_id}",
    )
    return sp.csc_matrix(
        (values.astype(np.float64), indices, indptr),
        shape=(n_genes, n_cells),
    )


def load_complete_annotation() -> tuple[pd.Series, pd.Series]:
    with h5py.File(COMPLETE_ATLAS, "r") as handle:
        cell_ids = pd.Index(read_elem(handle["obs/_index"]).astype(str))
        confidence = pd.Series(
            np.asarray(read_elem(handle["obs/cell_type_conf"]), dtype=float),
            index=cell_ids,
        )
        cell_type = pd.Series(
            np.asarray(read_elem(handle["obs/cell_type"]), dtype=str),
            index=cell_ids,
        )
    require(confidence.index.is_unique, "complete-atlas cell identifiers are not unique")
    return confidence, cell_type


def aggregate_complete_atlas(
    programs: dict[str, dict[str, float]],
) -> tuple[dict[str, tuple[pd.DataFrame, pd.DataFrame]], list[str]]:
    sample_to_donor = build_sample_to_donor()
    manifest = json.loads((AMBIENT_INPUT / "work/manifest.json").read_text())
    raw_work = OLD_RAW_EXPORT / "work"
    raw_genes = (raw_work / "genes.txt").read_text().splitlines()
    raw_position = {gene: index for index, gene in enumerate(raw_genes)}
    genes = sorted(set().union(*[set(weights) for weights in programs.values()]))
    require(set(genes) <= set(raw_position), "hero genes are absent from raw transport")
    gene_indices = np.array([raw_position[gene] for gene in genes], dtype=np.int64)
    confidence, atlas_cell_type = load_complete_annotation()

    expression_parts: dict[str, list[pd.DataFrame]] = {
        annotation_filter: [] for annotation_filter in prototype.ANNOTATION_FILTERS
    }
    metadata_parts: dict[str, list[pd.DataFrame]] = {
        annotation_filter: [] for annotation_filter in prototype.ANNOTATION_FILTERS
    }

    for dataset in sorted(manifest["datasets"]):
        dataset_work = raw_work / dataset
        dimensions = json.loads((dataset_work / "dims.json").read_text())
        meta = pd.read_csv(AMBIENT_INPUT / "work" / dataset / "meta.csv.gz", dtype=str)
        require(len(meta) == int(dimensions["ncells"]), f"metadata census drift: {dataset}")
        annotation_index = pd.Index(meta["cell_id"].astype(str))
        cell_confidence = confidence.reindex(annotation_index)
        joined_cell_type = atlas_cell_type.reindex(annotation_index)
        require(not cell_confidence.isna().any(), f"confidence join failed: {dataset}")
        require(not joined_cell_type.isna().any(), f"cell-type join failed: {dataset}")
        require(
            (joined_cell_type.to_numpy(str) == meta["cell_type"].to_numpy(str)).all(),
            f"cell-type annotation drift: {dataset}",
        )
        meta["cell_type_conf"] = cell_confidence.to_numpy(float)
        meta["donor"] = meta["sample"].map(sample_to_donor).fillna(meta["sample"])
        relevant = meta["cell_type"].isin(prototype.LINEAGES).to_numpy()
        require(relevant.any(), f"no requested lineages: {dataset}")
        keys = meta.loc[relevant, "donor"] + "||" + meta.loc[relevant, "cell_type"]
        levels = pd.Index(sorted(keys.unique()))
        key_to_code = pd.Series(np.arange(len(levels)), index=levels)
        full_keys = meta["donor"] + "||" + meta["cell_type"]
        codes = full_keys.map(key_to_code).fillna(-1).astype(int).to_numpy()
        sums = {
            annotation_filter: np.zeros((len(levels), len(genes)), dtype=np.float64)
            for annotation_filter in prototype.ANNOTATION_FILTERS
        }
        counts = {
            annotation_filter: np.zeros(len(levels), dtype=np.int64)
            for annotation_filter in prototype.ANNOTATION_FILTERS
        }

        offset = 0
        for chunk in dimensions["chunks"]:
            raw = read_raw_chunk(dataset_work, chunk, len(raw_genes))
            end = offset + raw.shape[1]
            totals = np.asarray(raw.sum(axis=0)).ravel()
            selected = raw[gene_indices, :].astype(np.float64).tocsc()
            scale = np.divide(
                1e4,
                totals,
                out=np.zeros_like(totals, dtype=np.float64),
                where=totals > 0,
            )
            selected = selected @ sp.diags(scale)
            selected.data = np.log1p(selected.data)
            block_codes = codes[offset:end]
            block_confidence = meta["cell_type_conf"].to_numpy(float)[offset:end]
            for annotation_filter, confidence_keep in (
                ("all_annotated_cells", np.ones(raw.shape[1], dtype=bool)),
                ("cell_type_conf_ge_0.90", block_confidence >= 0.90),
            ):
                keep = (block_codes >= 0) & confidence_keep
                if not keep.any():
                    continue
                kept_codes = block_codes[keep]
                indicator = sp.csr_matrix(
                    (
                        np.ones(len(kept_codes)),
                        (kept_codes, np.arange(len(kept_codes))),
                    ),
                    shape=(len(levels), len(kept_codes)),
                )
                sums[annotation_filter] += (
                    indicator @ selected[:, keep].T
                ).toarray()
                counts[annotation_filter] += np.bincount(
                    kept_codes, minlength=len(levels)
                )
            offset = end
            del raw, selected
        require(offset == len(meta), f"raw chunk cell census drift: {dataset}")

        for annotation_filter in prototype.ANNOTATION_FILTERS:
            eligible = counts[annotation_filter] > 0
            result_meta = pd.DataFrame(
                {
                    "group_id": levels[eligible],
                    "n_cells": counts[annotation_filter][eligible],
                }
            )
            result_meta[["donor", "lineage"]] = result_meta["group_id"].str.split(
                "||", n=1, expand=True, regex=False
            )
            result_meta["dataset"] = dataset
            result_expression = pd.DataFrame(
                sums[annotation_filter][eligible, :]
                / counts[annotation_filter][eligible, None],
                columns=genes,
            )
            expression_parts[annotation_filter].append(result_expression)
            metadata_parts[annotation_filter].append(result_meta)

    outputs = {}
    for annotation_filter in prototype.ANNOTATION_FILTERS:
        expression = pd.concat(expression_parts[annotation_filter], ignore_index=True)
        metadata = pd.concat(metadata_parts[annotation_filter], ignore_index=True)
        require(
            not metadata[["dataset", "donor", "lineage"]].duplicated().any(),
            f"duplicate donor-lineage profile: {annotation_filter}",
        )
        outputs[annotation_filter] = (expression, metadata)
    return outputs, genes


def main() -> None:
    prototype.CANDIDATE = CROSS_CANDIDATE
    prototype.RESULTS = CROSS_CANDIDATE / "results"
    prototype.ATLAS = COMPLETE_ATLAS
    prototype.aggregate_donor_lineage = aggregate_complete_atlas
    prototype.main()
    extra_provenance_path = CROSS_CANDIDATE / "results/complete_atlas_input_provenance.tsv"
    refuse_existing(extra_provenance_path)
    prototype_provenance = json.loads(
        (CROSS_CANDIDATE / "results/hero_lineage_provenance.json").read_text()
    )
    pd.DataFrame(
        [
            {
                "role": "complete_atlas_annotations",
                "path": str(COMPLETE_ATLAS.resolve()),
                "sha256": prototype_provenance["atlas_sha256"],
            },
            {
                "role": "checksum_verified_raw_transport_manifest",
                "path": str((AMBIENT_INPUT / "work/manifest.json").resolve()),
                "sha256": sha256(AMBIENT_INPUT / "work/manifest.json"),
            },
            {
                "role": "raw_gene_axis",
                "path": str((OLD_RAW_EXPORT / "work/genes.txt").resolve()),
                "sha256": sha256(OLD_RAW_EXPORT / "work/genes.txt"),
            },
        ]
    ).to_csv(extra_provenance_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
