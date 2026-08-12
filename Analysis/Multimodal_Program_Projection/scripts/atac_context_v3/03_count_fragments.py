#!/usr/bin/env python3
"""Count exact fragments in lineage consensus peaks and donor-pseudobulk them."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
import snapatac2 as snap

from atac_context_v3_lib import (
    LINEAGES,
    MIN_CELLS,
    RELEASE_ID,
    ContractError,
    candidate_root,
    default_candidate_root,
    require_new_path,
)
from snap_helpers import build_dataset, close_dataset, metadata_conditions


LEGACY_PEAKS = {
    "hepatocyte": "Hepatocytes_peaks.bed",
    "stellate": "Fibroblasts_peaks.bed",
    "macrophage": "Macrophages_peaks.bed",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", required=True, choices=["GSE244832", "GSE281367"])
    parser.add_argument("--candidate-root", type=Path, default=default_candidate_root())
    parser.add_argument("--tempdir", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=500)
    parser.add_argument(
        "--peak-space",
        choices=["consensus", "legacy_gse244832"],
        default="consensus",
    )
    return parser.parse_args()


def write_sparse(path: Path, matrix: sp.spmatrix) -> None:
    require_new_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as handle:
        mmwrite(handle, matrix.tocsr(), field="integer")


def main() -> None:
    args = parse_args()
    out = candidate_root(args.candidate_root)
    if not (out / "consensus_peak_manifest.tsv").is_file():
        raise ContractError("consensus peaks must complete before fragment counting")
    cohort_target = out / (
        "counts" if args.peak_space == "consensus" else "counts_legacy_gse244832"
    ) / args.cohort
    require_new_path(cohort_target)
    cohort_target.mkdir(parents=True)
    args.tempdir.mkdir(parents=True, exist_ok=True)

    dataset_file = args.tempdir / f"{args.cohort}.counts.h5ads"
    dataset, adatas = build_dataset(args.cohort, dataset_file)
    conditions = metadata_conditions(args.cohort)
    try:
        lineages = LINEAGES if args.peak_space == "consensus" else tuple(LEGACY_PEAKS)
        for lineage in lineages:
            if args.peak_space == "consensus":
                peak_bed = out / "peaks/consensus" / f"{lineage}.bed"
            else:
                peak_bed = (
                    Path(__file__).resolve().parents[4]
                    / "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2"
                    / LEGACY_PEAKS[lineage]
                )
            if not peak_bed.is_file():
                raise ContractError(f"missing {args.peak_space} peak BED: {peak_bed}")
            cell_matrix_file = args.tempdir / (
                f"{args.cohort}.{args.peak_space}.{lineage}.peak_matrix.h5ad"
            )
            peak_matrix = snap.pp.make_peak_matrix(
                dataset,
                peak_file=peak_bed,
                counting_strategy="fragment",
                chunk_size=args.chunk_size,
                file=cell_matrix_file,
                inplace=False,
            )
            try:
                obs = peak_matrix.obs
                try:
                    lineage_values = np.asarray(obs["lineage_v3"]).astype(str)
                    donor_values = np.asarray(obs["donor_id"]).astype(str)
                except (KeyError, IndexError) as error:
                    raise ContractError("peak matrix lost lineage/donor annotations") from error
                if len(lineage_values) != peak_matrix.n_obs or len(donor_values) != peak_matrix.n_obs:
                    raise ContractError("peak matrix lost lineage/donor annotations")
                mask = lineage_values == lineage
                donors = sorted(np.unique(donor_values[mask]).tolist())
                rows = []
                donor_meta = []
                matrix = peak_matrix.X
                for donor in donors:
                    donor_mask = mask & (donor_values == donor)
                    n_cells = int(donor_mask.sum())
                    summed = np.asarray(matrix[donor_mask, :].sum(axis=0)).ravel()
                    summed = np.rint(summed).astype(np.int64)
                    rows.append(sp.csr_matrix(summed.reshape(1, -1)))
                    donor_meta.append(
                        {
                            "release_id": RELEASE_ID,
                            "cohort": args.cohort,
                            "donor_id": donor,
                            "condition": conditions[donor],
                            "lineage": lineage,
                            "n_cells": n_cells,
                            "total_counted_fragments": int(summed.sum()),
                            "contrast_eligible": str(n_cells >= MIN_CELLS).upper(),
                            "exclusion_reason": "" if n_cells >= MIN_CELLS else "fewer_than_20_cells",
                        }
                    )
                pseudobulk = sp.vstack(rows, format="csr") if rows else sp.csr_matrix((0, peak_matrix.n_vars))
                write_sparse(cohort_target / f"{lineage}.mtx.gz", pseudobulk)
                pd.DataFrame(donor_meta).to_csv(
                    cohort_target / f"{lineage}.donors.tsv", sep="\t", index=False
                )
                pd.DataFrame(
                    {
                        "column_index_1based": np.arange(1, peak_matrix.n_vars + 1),
                        "peak_coordinate": [str(value) for value in peak_matrix.var_names],
                    }
                ).to_csv(cohort_target / f"{lineage}.peaks.tsv", sep="\t", index=False)
            finally:
                peak_matrix.close()
                if cell_matrix_file.exists():
                    cell_matrix_file.unlink()
    finally:
        close_dataset(dataset, adatas)
        if dataset_file.exists():
            dataset_file.unlink()
    print(f"Wrote exact donor fragment counts ({args.peak_space}): {cohort_target}")


if __name__ == "__main__":
    main()
