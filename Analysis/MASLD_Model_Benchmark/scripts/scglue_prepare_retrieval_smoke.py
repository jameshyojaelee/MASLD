#!/usr/bin/env python3
"""Prepare donor-held, identity-blinded same-nucleus scGLUE retrieval folds."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


class SCGLUERetrievalPrepareError(ValueError):
    """Raised when the identity-blinded scGLUE preparation differs."""


def decode(values: Any) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def read_csr(group: Any) -> Any:
    import numpy as np
    from scipy import sparse

    return sparse.csr_matrix(
        (
            np.asarray(group["data"][:]),
            np.asarray(group["indices"][:], dtype=np.int64),
            np.asarray(group["indptr"][:], dtype=np.int64),
        ),
        shape=tuple(int(value) for value in group["shape"][:]),
    )


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SCGLUERetrievalPrepareError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


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


def save_sparse(path: Path, matrix: Any) -> None:
    from scipy import sparse

    with path.open("xb") as handle:
        sparse.save_npz(handle, sparse.csr_matrix(matrix), compressed=True)


def pseudonym(kind: str, identifier: str) -> str:
    return sha256(f"masld-scglue-retrieval-v1\0{kind}\0{identifier}".encode()).hexdigest()


def parse_attributes(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in value.rstrip(";").split("; "):
        if " " not in item:
            continue
        key, raw = item.split(" ", 1)
        result[key] = raw.strip('"')
    return result


def gene_tss(gtf: Path, selected: set[str]) -> dict[str, tuple[str, int]]:
    result: dict[str, tuple[str, int]] = {}
    with gzip.open(gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = parse_attributes(fields[8])
            gene_id = attributes.get("gene_id", "").split(".", 1)[0]
            if gene_id not in selected:
                continue
            start = int(fields[3]) - 1
            end = int(fields[4])
            tss = start if fields[6] == "+" else end - 1
            if gene_id in result:
                raise SCGLUERetrievalPrepareError(f"gene TSS duplicated: {gene_id}")
            result[gene_id] = (fields[0], tss)
    if len(result) < int(0.90 * len(selected)):
        raise SCGLUERetrievalPrepareError("fewer than 90% of selected genes have GTF TSS")
    return result


def guidance_edges(
    genes: Sequence[Mapping[str, str]],
    peaks: Sequence[Mapping[str, str]],
    tss: Mapping[str, tuple[str, int]],
    max_distance: int = 150_000,
) -> list[dict[str, Any]]:
    by_chromosome: dict[str, list[tuple[int, str]]] = {}
    for row in genes:
        gene_id = row["ensembl_id"].split(".", 1)[0]
        if gene_id not in tss:
            continue
        chromosome, position = tss[gene_id]
        by_chromosome.setdefault(chromosome, []).append((position, gene_id))
    edges: list[dict[str, Any]] = []
    vertices = [row["ensembl_id"].split(".", 1)[0] for row in genes] + [
        row["peak_id"] for row in peaks
    ]
    if len(vertices) != len(set(vertices)):
        raise SCGLUERetrievalPrepareError("guidance vertices are not unique")
    for vertex in vertices:
        edges.append({"source": vertex, "target": vertex, "weight": "1", "sign": "1"})
    regulatory = 0
    for peak in peaks:
        chromosome = peak["chromosome"]
        midpoint = (int(peak["bed_start"]) + int(peak["bed_end"])) // 2
        for position, gene_id in by_chromosome.get(chromosome, []):
            distance = abs(midpoint - position)
            if distance > max_distance:
                continue
            weight = ((distance + 1000.0) / 1000.0) ** -0.75
            formatted = format(weight, ".17g")
            edges.append(
                {"source": gene_id, "target": peak["peak_id"], "weight": formatted, "sign": "1"}
            )
            edges.append(
                {"source": peak["peak_id"], "target": gene_id, "weight": formatted, "sign": "1"}
            )
            regulatory += 2
    if regulatory == 0 or any(
        not math.isfinite(float(row["weight"])) or float(row["weight"]) <= 0
        for row in edges
    ):
        raise SCGLUERetrievalPrepareError("guidance graph has no valid regulatory edge")
    return edges


def prepare(source: Path, prepared: Path, gtf: Path, output: Path) -> dict[str, Any]:
    import h5py
    import numpy as np

    if output.exists():
        raise SCGLUERetrievalPrepareError("output exists")
    with h5py.File(source, "r") as handle:
        if handle.attrs.get("pairing_level") != "same_nucleus":
            raise SCGLUERetrievalPrepareError("source is not same-nucleus")
        cell_ids = decode(handle["obs/cell_id"][:])
        donors = decode(handle["obs/donor_id"][:])
        rna = read_csr(handle["rna/counts_csr"])
        atac = read_csr(handle["atac/counts_csr"])
        source_gene_ids = decode(handle["rna/ensembl_id"][:])
        source_peak_ids = decode(handle["atac/peak_id"][:])
    if len(cell_ids) != 1000 or len(set(cell_ids)) != 1000 or len(set(donors)) != 39:
        raise SCGLUERetrievalPrepareError("source cell or donor census differs")
    cell_index = {cell_id: index for index, cell_id in enumerate(cell_ids)}
    selected_gene_union: set[str] = set()
    for fold in range(5):
        gene_fields, fold_genes = read_tsv(
            prepared / f"fold_{fold}/selected_genes.tsv"
        )
        if gene_fields != (
            "selected_index",
            "source_index",
            "ensembl_id",
            "gene_name",
        ):
            raise SCGLUERetrievalPrepareError("selected gene schema differs")
        selected_gene_union.update(
            row["ensembl_id"].split(".", 1)[0] for row in fold_genes
        )
    all_tss = gene_tss(gtf, selected_gene_union)
    output.mkdir(parents=True, mode=0o750)
    model_root = output / "model_inputs"
    evaluator_root = output / "evaluator_only"
    model_root.mkdir(mode=0o750)
    evaluator_root.mkdir(mode=0o750)
    folds: dict[str, Any] = {}
    for fold in range(5):
        source_root = prepared / f"fold_{fold}"
        _train_fields, train_rows = read_tsv(source_root / "training_rows.tsv")
        _query_fields, query_rows = read_tsv(source_root / "query_rows.tsv")
        gene_fields, genes = read_tsv(source_root / "selected_genes.tsv")
        peak_fields, peaks = read_tsv(source_root / "selected_peaks.tsv")
        if (
            gene_fields != ("selected_index", "source_index", "ensembl_id", "gene_name")
            or peak_fields
            != (
                "selected_index",
                "source_index",
                "peak_id",
                "chromosome",
                "bed_start",
                "bed_end",
            )
            or len(genes) != 2000
            or len(peaks) != 10_000
        ):
            raise SCGLUERetrievalPrepareError("selected feature roster differs")
        gene_indices = np.asarray([int(row["source_index"]) for row in genes], dtype=np.int64)
        peak_indices = np.asarray([int(row["source_index"]) for row in peaks], dtype=np.int64)
        gene_ids = [row["ensembl_id"].split(".", 1)[0] for row in genes]
        if [source_gene_ids[index].split(".", 1)[0] for index in gene_indices] != gene_ids:
            raise SCGLUERetrievalPrepareError("gene source join differs")
        if [source_peak_ids[index] for index in peak_indices] != [row["peak_id"] for row in peaks]:
            raise SCGLUERetrievalPrepareError("peak source join differs")
        train_indices = np.asarray([cell_index[row["cell_id"]] for row in train_rows], dtype=np.int64)
        query_indices = np.asarray([cell_index[row["cell_id"]] for row in query_rows], dtype=np.int64)
        if {donors[index] for index in train_indices} & {donors[index] for index in query_indices}:
            raise SCGLUERetrievalPrepareError("donor crosses retrieval fold")
        train_rna = rna[train_indices][:, gene_indices]
        train_atac = atac[train_indices][:, peak_indices]
        query_rna = rna[query_indices][:, gene_indices]
        order = np.asarray(
            sorted(
                range(len(query_indices)),
                key=lambda offset: pseudonym("atac_order", cell_ids[query_indices[offset]]),
            ),
            dtype=np.int64,
        )
        query_atac = atac[query_indices[order]][:, peak_indices]
        for matrix in (train_rna, train_atac, query_rna, query_atac):
            if np.any(np.asarray(matrix.sum(axis=1)).ravel() <= 0):
                raise SCGLUERetrievalPrepareError("observed selected assay row is empty")
        fold_model = model_root / f"fold_{fold}"
        fold_eval = evaluator_root / f"fold_{fold}"
        fold_model.mkdir(mode=0o750)
        fold_eval.mkdir(mode=0o750)
        save_sparse(fold_model / "training_rna.npz", train_rna)
        save_sparse(fold_model / "training_atac.npz", train_atac)
        save_sparse(fold_model / "query_rna.npz", query_rna)
        save_sparse(fold_model / "query_atac.npz", query_atac)
        training_ids = [pseudonym("train", cell_ids[index]) for index in train_indices]
        rna_query_ids = [pseudonym("rna_query", cell_ids[index]) for index in query_indices]
        atac_query_ids = [
            pseudonym("atac_query", cell_ids[query_indices[offset]]) for offset in order
        ]
        write_tsv(
            fold_model / "training_rows.tsv",
            ("training_id", "rna_state", "atac_state"),
            (
                {"training_id": identifier, "rna_state": "observed", "atac_state": "observed"}
                for identifier in training_ids
            ),
        )
        write_tsv(
            fold_model / "query_rna_rows.tsv",
            ("rna_query_id", "rna_state"),
            ({"rna_query_id": identifier, "rna_state": "observed"} for identifier in rna_query_ids),
        )
        write_tsv(
            fold_model / "query_atac_rows.tsv",
            ("atac_query_id", "atac_state"),
            ({"atac_query_id": identifier, "atac_state": "observed"} for identifier in atac_query_ids),
        )
        write_tsv(
            fold_model / "genes.tsv",
            ("feature_index", "gene_id"),
            ({"feature_index": index, "gene_id": gene_id} for index, gene_id in enumerate(gene_ids)),
        )
        write_tsv(
            fold_model / "peaks.tsv",
            ("feature_index", "peak_id"),
            ({"feature_index": index, "peak_id": row["peak_id"]} for index, row in enumerate(peaks)),
        )
        tss = {gene_id: all_tss[gene_id] for gene_id in gene_ids if gene_id in all_tss}
        edges = guidance_edges(genes, peaks, tss)
        write_tsv(fold_model / "guidance_edges.tsv", ("source", "target", "weight", "sign"), edges)
        pair_rows = []
        atac_id_by_cell = {
            cell_ids[query_indices[offset]]: atac_query_ids[position]
            for position, offset in enumerate(order)
        }
        for index, rna_id in zip(query_indices, rna_query_ids, strict=True):
            cell_id = cell_ids[index]
            pair_rows.append(
                {
                    "pair_hash": pseudonym("pair", cell_id),
                    "donor_hash": pseudonym("donor", donors[index]),
                    "rna_query_id": rna_id,
                    "atac_query_id": atac_id_by_cell[cell_id],
                }
            )
        write_tsv(
            fold_eval / "hidden_pairs.tsv",
            ("pair_hash", "donor_hash", "rna_query_id", "atac_query_id"),
            pair_rows,
        )
        folds[str(fold)] = {
            "training_nuclei": len(train_indices),
            "query_nuclei": len(query_indices),
            "training_donors": len({donors[index] for index in train_indices}),
            "query_donors": len({donors[index] for index in query_indices}),
            "genes": len(genes),
            "peaks": len(peaks),
            "guidance_edges": len(edges),
            "genes_with_gtf_tss": len(tss),
            "genes_without_gtf_tss_self_loop_only": len(genes) - len(tss),
            "query_atac_permuted": True,
        }
    receipt = {
        "schema_version": "masld-bench-scglue-retrieval-prepare-v1",
        "status": "pass",
        "dataset_id": "gse296875",
        "pairing_topology": "same_nucleus",
        "outer_unit": "donor",
        "folds": folds,
        "model_input_subtree": "model_inputs",
        "evaluator_only_subtree": "evaluator_only",
        "hidden_pair_map_available_to_model": False,
        "query_modality_ids_disjoint": True,
        "query_atac_row_order_permuted": True,
        "query_rna_state": "observed",
        "query_atac_state": "observed",
        "outcomes_read": False,
        "retrieval_metrics_calculated": False,
        "smoke_only": True,
        "champion_claim_allowed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    prepare(arguments.source, arguments.prepared, arguments.gtf, arguments.output)


if __name__ == "__main__":
    main()
