#!/usr/bin/env python3
"""Prepare five donor-held MultiVI folds with query ATAC withheld."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


LINEAGES = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)
N_SELECTED_PEAKS = 10_000


class MultiVIPrepareError(ValueError):
    """Raised when same-nucleus MultiVI preparation differs."""


def _decode(values: Any) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def _read_csr(group: Any) -> Any:
    import numpy as np
    from scipy import sparse

    matrix = sparse.csr_matrix(
        (
            np.asarray(group["data"][:]),
            np.asarray(group["indices"][:], dtype=np.int64),
            np.asarray(group["indptr"][:], dtype=np.int64),
        ),
        shape=tuple(int(value) for value in group["shape"][:]),
    )
    if matrix.data.size and (
        (matrix.data < 0).any() or (matrix.data != matrix.data.astype(int)).any()
    ):
        raise MultiVIPrepareError("raw count matrix differs")
    return matrix


def fold_index(donor: str, seed: int = 20260821) -> int:
    return int.from_bytes(sha256(f"{seed}\0{donor}".encode()).digest()[:8], "big") % 5


def deterministic_peak_indices(peak_ids: Sequence[str], n_peaks: int) -> tuple[int, ...]:
    ranked = sorted(
        range(len(peak_ids)),
        key=lambda index: (
            sha256(f"masld-rna-atac-smoke-v1\0{peak_ids[index]}".encode()).digest(),
            peak_ids[index],
            index,
        ),
    )
    return tuple(sorted(ranked[:n_peaks]))


def training_hvgs(matrix: Any, n_genes: int) -> tuple[int, ...]:
    import numpy as np

    depth = np.asarray(matrix.sum(axis=1)).ravel()
    if (depth <= 0).any():
        raise MultiVIPrepareError("training RNA contains an empty row")
    normalized = matrix.multiply((10_000.0 / depth)[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)
    mean = np.asarray(normalized.mean(axis=0)).ravel()
    variance = np.maximum(
        np.asarray(normalized.power(2).mean(axis=0)).ravel() - mean * mean, 0.0
    )
    score = variance / np.maximum(mean, 1.0e-12)
    ranked = sorted(range(matrix.shape[1]), key=lambda index: (-score[index], index))
    return tuple(sorted(ranked[:n_genes]))


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def _save_sparse(path: Path, matrix: Any) -> None:
    from scipy import sparse

    with path.open("xb") as handle:
        sparse.save_npz(handle, sparse.csr_matrix(matrix), compressed=True)


def prepare(source: Path, output: Path) -> dict[str, Any]:
    import h5py
    import numpy as np

    if output.exists():
        raise MultiVIPrepareError("output exists")
    with h5py.File(source, "r") as handle:
        if (
            handle.attrs.get("schema_version") != "masld-bench-multimodal-h5-v1"
            or handle.attrs.get("pairing_level") != "same_nucleus"
        ):
            raise MultiVIPrepareError("same-nucleus HDF5 contract differs")
        cell_ids = _decode(handle["obs/cell_id"][:])
        donors = _decode(handle["obs/donor_id"][:])
        lineages = _decode(handle["obs/broad_label"][:])
        rna_states = _decode(handle["obs/rna_status"][:])
        atac_states = _decode(handle["obs/atac_status"][:])
        rna = _read_csr(handle["rna/counts_csr"])
        atac = _read_csr(handle["atac/counts_csr"])
        gene_ids = _decode(handle["rna/ensembl_id"][:])
        gene_names = _decode(handle["rna/gene_name"][:])
        peak_ids = _decode(handle["atac/peak_id"][:])
        chromosomes = _decode(handle["atac/chromosome"][:])
        starts = np.asarray(handle["atac/bed_start_0based"][:], dtype=np.int64)
        ends = np.asarray(handle["atac/bed_end_half_open"][:], dtype=np.int64)
    if (
        len(cell_ids) != 1000
        or len(set(cell_ids)) != 1000
        or len(set(donors)) != 39
        or set(lineages) != set(LINEAGES)
        or set(rna_states) != {"observed"}
        or set(atac_states) != {"observed"}
        or rna.shape != (1000, 36_601)
        or atac.shape != (1000, 306_706)
    ):
        raise MultiVIPrepareError("source row, modality, or missingness contract differs")
    selected_peaks = deterministic_peak_indices(peak_ids, N_SELECTED_PEAKS)
    output.mkdir(parents=True, mode=0o750)
    folds: dict[str, Any] = {}
    assignments = np.asarray([fold_index(donor) for donor in donors], dtype=np.int8)
    for donor in set(donors):
        if len({int(assignments[i]) for i, value in enumerate(donors) if value == donor}) != 1:
            raise MultiVIPrepareError("donor crosses outer folds")
    for fold in range(5):
        root = output / f"fold_{fold}"
        root.mkdir(mode=0o750)
        train_mask = assignments != fold
        query_mask = ~train_mask
        selected_genes = training_hvgs(rna[train_mask], 2000)
        train_rna = rna[train_mask][:, selected_genes]
        train_atac = atac[train_mask][:, selected_peaks]
        query_rna = rna[query_mask][:, selected_genes]
        if (
            train_rna.shape[0] == 0
            or query_rna.shape[0] == 0
            or np.any(np.asarray(train_rna.sum(axis=1)).ravel() <= 0)
            or np.any(np.asarray(train_atac.sum(axis=1)).ravel() <= 0)
            or np.any(np.asarray(query_rna.sum(axis=1)).ravel() <= 0)
        ):
            raise MultiVIPrepareError("fold contains an empty observed assay")
        _save_sparse(root / "training_rna.npz", train_rna)
        _save_sparse(root / "training_atac.npz", train_atac)
        _save_sparse(root / "query_rna.npz", query_rna)
        row_fields = (
            "cell_id",
            "donor_id",
            "lineage",
            "rna_state",
            "atac_state",
        )
        write_tsv(
            root / "training_rows.tsv",
            row_fields,
            (
                {
                    "cell_id": cell_ids[index],
                    "donor_id": donors[index],
                    "lineage": lineages[index],
                    "rna_state": "observed",
                    "atac_state": "observed",
                }
                for index in np.flatnonzero(train_mask)
            ),
        )
        write_tsv(
            root / "query_rows.tsv",
            row_fields,
            (
                {
                    "cell_id": cell_ids[index],
                    "donor_id": donors[index],
                    "lineage": lineages[index],
                    "rna_state": "observed",
                    "atac_state": "structurally_missing",
                }
                for index in np.flatnonzero(query_mask)
            ),
        )
        write_tsv(
            root / "selected_genes.tsv",
            ("selected_index", "source_index", "ensembl_id", "gene_name"),
            (
                {
                    "selected_index": selected_index,
                    "source_index": source_index,
                    "ensembl_id": gene_ids[source_index],
                    "gene_name": gene_names[source_index],
                }
                for selected_index, source_index in enumerate(selected_genes)
            ),
        )
        write_tsv(
            root / "selected_peaks.tsv",
            (
                "selected_index",
                "source_index",
                "peak_id",
                "chromosome",
                "bed_start",
                "bed_end",
            ),
            (
                {
                    "selected_index": selected_index,
                    "source_index": source_index,
                    "peak_id": peak_ids[source_index],
                    "chromosome": chromosomes[source_index],
                    "bed_start": int(starts[source_index]),
                    "bed_end": int(ends[source_index]),
                }
                for selected_index, source_index in enumerate(selected_peaks)
            ),
        )
        folds[str(fold)] = {
            "training_nuclei": int(train_mask.sum()),
            "query_nuclei": int(query_mask.sum()),
            "training_donors": len({donors[i] for i in np.flatnonzero(train_mask)}),
            "query_donors": len({donors[i] for i in np.flatnonzero(query_mask)}),
            "training_rna_state": "observed",
            "training_atac_state": "observed",
            "query_rna_state": "observed",
            "query_atac_state": "structurally_missing",
            "query_atac_exported": False,
            "selected_genes": 2000,
            "selected_peaks": N_SELECTED_PEAKS,
        }
    receipt = {
        "schema_version": "masld-bench-multivi-smoke-prepare-v1",
        "status": "pass",
        "dataset_id": "gse296875",
        "view_id": "gse296875_rna_atac_smoke_1000_v1",
        "pairing_topology": "same_nucleus",
        "outer_unit": "donor",
        "split_seed": 20260821,
        "folds": folds,
        "peak_selection": "identifier_only_sha256_common_10000",
        "gene_selection": "training_fold_only_log_cpm_dispersion",
        "held_atac_read_in_prepare": True,
        "held_atac_exported_from_prepare": False,
        "fit_predict_source_hdf5_required": False,
        "gse244832_false_pairing": False,
        "gse281367_false_pairing": False,
        "outcomes_read": False,
        "smoke_only": True,
        "champion_claim_allowed": False,
    }
    (output / "prepare_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    print(json.dumps(prepare(arguments.source, arguments.output), sort_keys=True))


if __name__ == "__main__":
    main()
