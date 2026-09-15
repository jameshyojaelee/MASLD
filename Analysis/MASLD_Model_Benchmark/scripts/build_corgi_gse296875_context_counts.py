#!/usr/bin/env python3
"""Build masked donor-by-lineage raw RNA contexts for regular Corgi."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re

import h5py
import numpy as np
from scipy import sparse


MEMBERSHIP_SHA256 = "582387b78ec25215d45820c24aace12f130b826e91755bf19ac4581730611e64"
SOURCE_LOCK_SHA256 = "6be686b4187309f038097cf73d3f2a8dd34e72d666bc92733a53e99a1ccbb620"
CORGI_SHA256 = "73aad9ab9ad3e080a2ccaa0cee1a1d727d126a454f0e69f74f343a900b168c0f"
GTF_SHA256 = "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
ATTRIBUTE = re.compile(r'(?P<key>[A-Za-z_]+) "(?P<value>[^"]*)";')


class CorgiContextError(RuntimeError):
    """Raised when a source, mapping, or raw-count requirement differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def decode(values: np.ndarray) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def read_tsv(path: Path, *, compressed: bool = False) -> list[dict[str, str]]:
    opener = gzip.open if compressed else Path.open
    kwargs = {"mode": "rt", "encoding": "utf-8", "newline": ""} if compressed else {
        "mode": "r", "encoding": "utf-8", "newline": ""
    }
    with opener(path, **kwargs) as handle:  # type: ignore[arg-type]
        return list(csv.DictReader(handle, delimiter="\t"))


def resolve_roster(
    roster_ids: list[str],
    roster_symbols: list[str],
    source_ids: list[str],
    source_symbols: list[str],
) -> list[dict[str, object]]:
    if len(roster_ids) != 2_891 or len(roster_symbols) != 2_891:
        raise CorgiContextError("Corgi roster length differs")
    if len(set(roster_ids)) != len(roster_ids):
        raise CorgiContextError("Corgi Ensembl roster is duplicated")
    by_id = {value.split(".", 1)[0]: index for index, value in enumerate(source_ids)}
    if len(by_id) != len(source_ids):
        raise CorgiContextError("source stable Ensembl IDs are duplicated")
    by_symbol: dict[str, list[int]] = {}
    for index, symbol in enumerate(source_symbols):
        by_symbol.setdefault(symbol, []).append(index)
    output: list[dict[str, object]] = []
    used_source_indices: set[int] = set()
    for position, (gene_id, symbol) in enumerate(zip(roster_ids, roster_symbols, strict=True)):
        source_index = by_id.get(gene_id)
        state = "exact_stable_id"
        if source_index is None:
            candidates = by_symbol.get(symbol, [])
            if len(candidates) == 1:
                source_index = candidates[0]
                state = "unique_symbol_rescue"
            elif len(candidates) > 1:
                state = "ambiguous_symbol_structurally_missing"
            else:
                state = "absent_source_annotation_structurally_missing"
        if source_index is not None:
            if source_index in used_source_indices:
                raise CorgiContextError("one source feature maps to multiple Corgi genes")
            used_source_indices.add(source_index)
        output.append(
            {
                "corgi_position": position,
                "corgi_ensembl_stable_id": gene_id,
                "corgi_hgnc_symbol": symbol,
                "source_feature_index": source_index,
                "source_ensembl_stable_id": (
                    source_ids[source_index].split(".", 1)[0] if source_index is not None else ""
                ),
                "source_symbol": source_symbols[source_index] if source_index is not None else "",
                "mapping_state": state,
                "observed_feature": source_index is not None,
            }
        )
    return output


def gencode_v49_lengths(gtf: Path, roster: list[dict[str, object]]) -> dict[str, int]:
    wanted = {str(row["corgi_ensembl_stable_id"]) for row in roster}
    wanted.update(
        str(row["source_ensembl_stable_id"])
        for row in roster
        if row["source_ensembl_stable_id"]
    )
    intervals: dict[str, list[tuple[str, int, int]]] = {gene: [] for gene in wanted}
    with gzip.open(gtf, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = {match["key"]: match["value"] for match in ATTRIBUTE.finditer(fields[8])}
            stable = attrs.get("gene_id", "").split(".", 1)[0]
            if stable in intervals:
                intervals[stable].append((fields[0], int(fields[3]) - 1, int(fields[4])))
    lengths: dict[str, int] = {}
    for gene, pieces in intervals.items():
        by_contig: dict[str, list[tuple[int, int]]] = {}
        for contig, start, end in pieces:
            by_contig.setdefault(contig, []).append((start, end))
        total = 0
        for contig_pieces in by_contig.values():
            current_start = current_end = -1
            for start, end in sorted(contig_pieces):
                if start > current_end:
                    total += max(0, current_end - current_start)
                    current_start, current_end = start, end
                else:
                    current_end = max(current_end, end)
            total += max(0, current_end - current_start)
        if total:
            lengths[gene] = total
    return lengths


def aggregate_well(
    raw_path: Path,
    rows: list[dict[str, str]],
    unit_index: dict[tuple[str, str], int],
    roster: list[dict[str, object]] | None,
) -> tuple[np.ndarray, list[dict[str, object]], list[str], list[str]]:
    with h5py.File(raw_path, "r") as handle:
        matrix = handle["matrix"]
        feature_types = matrix["features/feature_type"][:]
        rna_indices = np.flatnonzero(feature_types == b"Gene Expression")
        source_ids = decode(matrix["features/id"][:][rna_indices])
        source_symbols = decode(matrix["features/name"][:][rna_indices])
        if len(source_ids) != 36_601:
            raise CorgiContextError("source RNA feature count differs")
        if roster is None:
            raise CorgiContextError("resolved roster is required")
        source_positions = [
            int(row["source_feature_index"])
            for row in roster
            if row["observed_feature"]
        ]
        roster_positions = [
            int(row["corgi_position"])
            for row in roster
            if row["observed_feature"]
        ]
        barcodes = decode(matrix["barcodes"][:])
        barcode_to_column = {barcode: index for index, barcode in enumerate(barcodes)}
        if len(barcode_to_column) != len(barcodes):
            raise CorgiContextError("source RNA barcodes are duplicated")
        group_rows: list[int] = []
        group_columns: list[int] = []
        for row in rows:
            column = barcode_to_column.get(row["raw_barcode"])
            if column is None:
                raise CorgiContextError("membership barcode is absent from raw RNA")
            group_rows.append(column)
            group_columns.append(unit_index[(row["donor_id"], row["lineage_id"])])
        group = sparse.csc_matrix(
            (np.ones(len(group_rows), dtype=np.uint8), (group_rows, group_columns)),
            shape=(len(barcodes), len(unit_index)),
        )
        raw_data = matrix["data"][:]
        if np.any(raw_data < 0) or np.any(raw_data != np.floor(raw_data)):
            raise CorgiContextError("source RNA values are not nonnegative integer counts")
        source_matrix = sparse.csc_matrix(
            (raw_data, matrix["indices"][:], matrix["indptr"][:]),
            shape=tuple(int(value) for value in matrix["shape"][:]),
        )
        if source_matrix.shape != (len(feature_types), len(barcodes)):
            raise CorgiContextError("source RNA sparse matrix shape differs")
        selected = source_matrix[rna_indices[np.asarray(source_positions)], :]
        aggregated = selected @ group
        dense = np.asarray(aggregated.todense(), dtype=np.uint64)
        output = np.zeros((len(unit_index), 2_891), dtype=np.uint64)
        output[:, np.asarray(roster_positions)] = dense.T
        return output, roster, source_ids, source_symbols


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--membership-root", type=Path, required=True)
    parser.add_argument("--source-lock", type=Path, required=True)
    parser.add_argument("--corgi-root", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    checks = (
        (args.membership_root / "ARTIFACTS.json", MEMBERSHIP_SHA256),
        (args.source_lock, SOURCE_LOCK_SHA256),
        (args.corgi_root / "ARTIFACTS.json", CORGI_SHA256),
        (args.gtf, GTF_SHA256),
    )
    for path, expected in checks:
        if digest(path) != expected:
            raise CorgiContextError(f"input checksum differs: {path}")
    acquisition_source = next((args.corgi_root / "source").glob("corgi-*/data"))
    roster_ids = (acquisition_source / "trans_regulators_final_ensg.txt").read_text().splitlines()
    roster_symbols = (acquisition_source / "trans_regulators_final_hgnc.txt").read_text().splitlines()
    reference = np.load(acquisition_source / "tf_reference.npy", allow_pickle=False)
    if reference.shape != (2_891,) or not np.isfinite(reference).all():
        raise CorgiContextError("Corgi reference vector differs")
    membership = read_tsv(args.membership_root / "cell_membership.tsv.gz", compressed=True)
    membership = [row for row in membership if row["analysis_role"] == "primary"]
    if not membership or any(row["lineage_id"] not in LINEAGES for row in membership):
        raise CorgiContextError("primary membership lineage differs")
    donor_folds = read_tsv(args.membership_root / "donor_folds.tsv")
    donors = [row["donor_id"] for row in donor_folds]
    if len(donors) != 39 or len(set(donors)) != 39:
        raise CorgiContextError("donor fold roster differs")
    units = [(donor, lineage) for donor in donors for lineage in LINEAGES]
    unit_index = {unit: index for index, unit in enumerate(units)}
    fold_by_donor = {row["donor_id"]: int(row["outer_fold"]) for row in donor_folds}
    source_lock = json.loads(args.source_lock.read_text())
    raw_sources = source_lock["raw_rna_sources"]
    by_well: dict[str, list[dict[str, str]]] = {well: [] for well in raw_sources}
    for row in membership:
        by_well[row["well_id"]].append(row)
    first_well = sorted(raw_sources)[0]
    first_raw = (args.source_lock.parent / raw_sources[first_well]["path"]).resolve(strict=True)
    if digest(first_raw) != raw_sources[first_well]["sha256"]:
        raise CorgiContextError("first raw RNA checksum differs")
    with h5py.File(first_raw, "r") as handle:
        mask = handle["matrix/features/feature_type"][:] == b"Gene Expression"
        source_ids = decode(handle["matrix/features/id"][:][mask])
        source_symbols = decode(handle["matrix/features/name"][:][mask])
    roster = resolve_roster(roster_ids, roster_symbols, source_ids, source_symbols)
    counts = np.zeros((len(units), 2_891), dtype=np.uint64)
    source_feature_hash: str | None = None
    for well in sorted(raw_sources):
        raw = (args.source_lock.parent / raw_sources[well]["path"]).resolve(strict=True)
        if raw.stat().st_size != int(raw_sources[well]["bytes"]) or digest(raw) != raw_sources[well]["sha256"]:
            raise CorgiContextError(f"raw RNA identity differs: {well}")
        well_counts, _, ids, symbols = aggregate_well(raw, by_well[well], unit_index, roster)
        feature_hash = sha256("\n".join(f"{a}\t{b}" for a, b in zip(ids, symbols, strict=True)).encode()).hexdigest()
        if source_feature_hash is None:
            source_feature_hash = feature_hash
        elif source_feature_hash != feature_hash:
            raise CorgiContextError("RNA feature axis differs across wells")
        counts += well_counts
    nuclei = Counter((row["donor_id"], row["lineage_id"]) for row in membership)
    unit_rows: list[dict[str, object]] = []
    for index, (donor, lineage) in enumerate(units):
        n = nuclei.get((donor, lineage), 0)
        unit_rows.append(
            {
                "unit_index": index,
                "donor_id": donor,
                "outer_fold": fold_by_donor[donor],
                "lineage_id": lineage,
                "nuclei": n,
                "rna_state": "observed" if n else "structurally_missing",
                "library_umis_observed_features": int(counts[index].sum()),
                "biological_unit": "donor_by_lineage_pseudobulk",
            }
        )
    lengths = gencode_v49_lengths(args.gtf, roster)
    for row in roster:
        gene = str(row["source_ensembl_stable_id"] or row["corgi_ensembl_stable_id"])
        row["gencode_v49_exon_union_length"] = lengths.get(gene, "")
        row["length_state"] = "observed" if gene in lengths else "unavailable"
    observed_mask = np.asarray([bool(row["observed_feature"]) for row in roster], dtype=np.bool_)
    if int(observed_mask.sum()) != 2_857:
        raise CorgiContextError("observed Corgi gene census differs")
    args.output.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(
        args.output / "raw_context_counts.npz",
        counts=counts,
        observed_gene_mask=observed_mask,
        tf_reference=np.asarray(reference, dtype=np.float32),
    )
    write_tsv(args.output / "units.tsv", unit_rows)
    write_tsv(args.output / "corgi_gene_crosswalk.tsv", roster)
    audit = {
        "schema_version": "masld-bench-corgi-gse296875-context-counts-v1",
        "status": "pass_raw_context_counts_mapper_fit_blocked",
        "dataset_id": "gse296875",
        "biological_unit": "donor_by_lineage_pseudobulk",
        "donors": 39,
        "lineages": list(LINEAGES),
        "units": len(units),
        "observed_units": sum(row["rna_state"] == "observed" for row in unit_rows),
        "nuclei": len(membership),
        "corgi_genes": 2_891,
        "observed_source_genes": int(observed_mask.sum()),
        "structurally_missing_source_genes": int((~observed_mask).sum()),
        "mapping_counts": dict(sorted(Counter(str(row["mapping_state"]) for row in roster).items())),
        "missing_encoded_as_biological_zero": False,
        "counts_storage_zero_requires_observed_gene_mask": True,
        "model_training_activated": False,
        "outcomes_read": False,
        "sealed_features_or_labels_read": False,
        "remaining_blockers": [
            "fit_and_freeze_mapper_inside_each_outer_training_partition",
            "define_nonzero_imputation_for_34_structurally_missing_genes",
            "freeze_assay_native_ATAC_endpoint_and_sequence_windows",
            "resolve_Corgi_license_native_reference_and_numeric_parity_gates",
        ],
    }
    (args.output / "audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(json.dumps(audit, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
