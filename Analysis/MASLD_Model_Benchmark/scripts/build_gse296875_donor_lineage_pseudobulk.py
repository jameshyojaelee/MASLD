#!/usr/bin/env python3
"""Build the donor-by-lineage RNA pseudobulk fixture for GSE296875.

Units are donor by lineage over the frozen fragment-membership partition: five
primary lineages that are confirmatory-family members, plus two secondary
lineages that are reported and never conditional on.  The inferential unit remains
the donor; the unit count is a storage fact, not a sample size.

A unit below the frozen minimum-cell threshold is marked ``insufficient_cells``
and keeps its observed nuclei count.  It is never written as a zero row and
never dropped, because a lineage that was not sampled deeply enough is missing,
not silent.

This module deliberately does not import the benchmark package.  It runs on the
scanpy environment, which has NumPy, SciPy, and h5py but is Python 3.10 and
therefore cannot import ``masld_bench``.  The scoring module is unimportable
from here, so no outcome can be read on this path by construction.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse


class PseudobulkError(RuntimeError):
    """Raised when a frozen input or a structural expectation does not hold."""


PRIMARY_LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
SECONDARY_LINEAGES = ("endothelial_cell", "b_cell")
LINEAGES = PRIMARY_LINEAGES + SECONDARY_LINEAGES

MINIMUM_CELLS_PER_UNIT = 20
EXPECTED_DONORS = 39
EXPECTED_NUCLEI = 68_398
EXPECTED_RNA_FEATURES = 36_601

BOUND = {
    "fragment_membership": "582387b78ec25215d45820c24aace12f130b826e91755bf19ac4581730611e64",
    "barcode_donor_join": "ead9a54099fe2dfb47517497c072feec53fa5ebf9c0392e38ee34a164f8d64e8",
    "phenotype_endpoints": "6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140",
    "donor_folds": "e6bc9161acaab508c298c49036af1ad10a720949178a14cba42282933195d416",
    "corgi_context_counts": "6946c329e74cd1a6cd72fe37faaf805183ff24c729ef9e549da11c49ca02fef8",
}
SOURCE_LOCK_SHA256 = "6be686b4187309f038097cf73d3f2a8dd34e72d666bc92733a53e99a1ccbb620"


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def decode(values: np.ndarray) -> list[str]:
    return [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values]


def read_tsv(path: Path, *, compressed: bool = False) -> list[dict[str, str]]:
    if compressed:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def verify_bound(root: Path, corgi_root: Path, source_lock: Path) -> None:
    paths = {
        "fragment_membership": root / "executions/gse296875-fragment-membership-21063829",
        "barcode_donor_join": root / "executions/gse296875-barcode-donor-join-20260825",
        "phenotype_endpoints": root / "executions/gse296875-phenotype-endpoints-20260825",
        "donor_folds": root / "executions/gse296875-donor-folds-20260825",
        "corgi_context_counts": corgi_root,
    }
    for name, expected in BOUND.items():
        observed = digest(paths[name] / "ARTIFACTS.json")
        if observed != expected:
            raise PseudobulkError(f"bound fixture changed: {name}")
    if digest(source_lock) != SOURCE_LOCK_SHA256:
        raise PseudobulkError("raw RNA source lock changed")


def aggregate_well(
    raw_path: Path,
    rows: list[dict[str, str]],
    unit_index: dict[tuple[str, str], int],
    n_units: int,
) -> tuple[np.ndarray, str]:
    """Sum raw UMIs over the nuclei of each donor-by-lineage unit in one well."""

    with h5py.File(raw_path, "r") as handle:
        matrix = handle["matrix"]
        feature_types = matrix["features/feature_type"][:]
        rna_indices = np.flatnonzero(feature_types == b"Gene Expression")
        if rna_indices.size != EXPECTED_RNA_FEATURES:
            raise PseudobulkError("source RNA feature count differs")
        ids = decode(matrix["features/id"][:][rna_indices])
        symbols = decode(matrix["features/name"][:][rna_indices])
        feature_hash = sha256(
            "\n".join(f"{a}\t{b}" for a, b in zip(ids, symbols)).encode("utf-8")
        ).hexdigest()

        barcodes = decode(matrix["barcodes"][:])
        column_of = {barcode: index for index, barcode in enumerate(barcodes)}
        if len(column_of) != len(barcodes):
            raise PseudobulkError("source RNA barcodes are duplicated")

        rows_index: list[int] = []
        columns_index: list[int] = []
        for row in rows:
            column = column_of.get(row["raw_barcode"])
            if column is None:
                raise PseudobulkError(
                    f"membership barcode absent from raw RNA: {row['cell_id']}"
                )
            rows_index.append(column)
            columns_index.append(unit_index[(row["donor_id"], row["lineage_id"])])
        group = sparse.csc_matrix(
            (
                np.ones(len(rows_index), dtype=np.uint8),
                (rows_index, columns_index),
            ),
            shape=(len(barcodes), n_units),
        )

        data = matrix["data"][:]
        if np.any(data < 0) or np.any(data != np.floor(data)):
            raise PseudobulkError("source RNA values are not nonnegative integers")
        source = sparse.csc_matrix(
            (data, matrix["indices"][:], matrix["indptr"][:]),
            shape=tuple(int(v) for v in matrix["shape"][:]),
        )
        if source.shape != (len(feature_types), len(barcodes)):
            raise PseudobulkError("source RNA sparse matrix shape differs")
        aggregated = source[rna_indices, :] @ group
        return np.asarray(aggregated.todense(), dtype=np.uint64).T, feature_hash


def corgi_cross_check(
    counts: np.ndarray,
    units: list[tuple[str, str]],
    gene_ids: list[str],
    corgi_root: Path,
) -> dict[str, object]:
    """Reconcile against the frozen Corgi context counts on their shared cells.

    The Corgi output file was built for a different task over a 2,891-gene roster
    and the five primary lineages.  Where the two overlap they must agree
    exactly, otherwise one of the two aggregations is wrong.
    """

    corgi_units = read_tsv(corgi_root / "units.tsv")
    crosswalk = read_tsv(corgi_root / "corgi_gene_crosswalk.tsv")
    with np.load(corgi_root / "raw_context_counts.npz") as payload:
        corgi_counts = payload["counts"]

    unit_position = {unit: index for index, unit in enumerate(units)}
    stable_of = {value.split(".", 1)[0]: index for index, value in enumerate(gene_ids)}

    gene_pairs = [
        (int(row["corgi_position"]), stable_of[row["source_ensembl_stable_id"]])
        for row in crosswalk
        if row["observed_feature"] == "True"
        and row["source_ensembl_stable_id"] in stable_of
    ]
    if len(gene_pairs) != 2_857:
        raise PseudobulkError(
            f"Corgi observed-gene overlap differs: {len(gene_pairs)}"
        )

    compared = 0
    for row in corgi_units:
        position = unit_position[(row["donor_id"], row["lineage_id"])]
        mine = counts[position]
        theirs = corgi_counts[int(row["unit_index"])]
        for corgi_position, my_index in gene_pairs:
            if int(theirs[corgi_position]) != int(mine[my_index]):
                raise PseudobulkError(
                    "Corgi cross-check disagrees at "
                    f"{row['donor_id']}/{row['lineage_id']} gene {corgi_position}"
                )
        compared += 1
    return {
        "corgi_units_compared": compared,
        "corgi_genes_compared": len(gene_pairs),
        "corgi_cell_values_compared": compared * len(gene_pairs),
        "corgi_agreement": "exact",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--source-lock", type=Path, required=True)
    parser.add_argument("--corgi-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    root = arguments.root
    verify_bound(root, arguments.corgi_root, arguments.source_lock)

    membership_root = root / "executions/gse296875-fragment-membership-21063829"
    membership = read_tsv(membership_root / "cell_membership.tsv.gz", compressed=True)
    if len(membership) != EXPECTED_NUCLEI:
        raise PseudobulkError("frozen nucleus roster size differs")
    if {row["lineage_id"] for row in membership} != set(LINEAGES):
        raise PseudobulkError("frozen lineage set differs")

    folds = read_tsv(
        root / "executions/gse296875-donor-folds-20260825/split_lock/donor_folds.tsv"
    )
    donors = [row["donor_id"] for row in folds]
    if len(donors) != EXPECTED_DONORS or len(set(donors)) != EXPECTED_DONORS:
        raise PseudobulkError("donor roster differs")
    fold_of = {row["donor_id"]: int(row["outer_fold"]) for row in folds}
    well_of = {row["donor_id"]: row["well_id"] for row in folds}
    adult_of = {row["donor_id"]: row["is_adult"] == "true" for row in folds}

    units = [(donor, lineage) for donor in donors for lineage in LINEAGES]
    unit_index = {unit: index for index, unit in enumerate(units)}

    source_lock = json.loads(arguments.source_lock.read_text())
    raw_sources = source_lock["raw_rna_sources"]
    by_well: dict[str, list[dict[str, str]]] = {well: [] for well in raw_sources}
    for row in membership:
        by_well[row["well_id"]].append(row)

    counts = np.zeros((len(units), EXPECTED_RNA_FEATURES), dtype=np.uint64)
    feature_hash: str | None = None
    gene_ids: list[str] = []
    gene_symbols: list[str] = []
    for well in sorted(raw_sources):
        record = raw_sources[well]
        path = (arguments.source_lock.parent / record["path"]).resolve(strict=True)
        if path.stat().st_size != int(record["bytes"]) or digest(path) != record["sha256"]:
            raise PseudobulkError(f"raw RNA identity differs: {well}")
        well_counts, observed_hash = aggregate_well(
            path, by_well[well], unit_index, len(units)
        )
        if feature_hash is None:
            feature_hash = observed_hash
            with h5py.File(path, "r") as handle:
                matrix = handle["matrix"]
                mask = matrix["features/feature_type"][:] == b"Gene Expression"
                gene_ids = decode(matrix["features/id"][:][mask])
                gene_symbols = decode(matrix["features/name"][:][mask])
        elif feature_hash != observed_hash:
            raise PseudobulkError("RNA feature axis differs across wells")
        counts += well_counts

    nuclei = Counter((row["donor_id"], row["lineage_id"]) for row in membership)
    unit_rows: list[dict[str, object]] = []
    for index, (donor, lineage) in enumerate(units):
        observed_nuclei = nuclei.get((donor, lineage), 0)
        sufficient = observed_nuclei >= MINIMUM_CELLS_PER_UNIT
        unit_rows.append(
            {
                "unit_index": index,
                "donor_id": donor,
                "well_id": well_of[donor],
                "outer_fold": fold_of[donor],
                "is_adult": str(adult_of[donor]).lower(),
                "lineage_id": lineage,
                "analysis_role": "primary" if lineage in PRIMARY_LINEAGES else "secondary",
                "in_confirmatory_family": str(lineage in PRIMARY_LINEAGES).lower(),
                "nuclei": observed_nuclei,
                "minimum_nuclei_per_unit": MINIMUM_CELLS_PER_UNIT,
                "rna_state": "observed" if sufficient else "insufficient_cells",
                "library_umis": int(counts[index].sum()),
                "detected_genes": int((counts[index] > 0).sum()),
                "biological_unit": "donor_by_lineage_pseudobulk",
            }
        )

    gene_rows = [
        {
            "gene_index": index,
            "ensembl_id": gene_id,
            "ensembl_stable_id": gene_id.split(".", 1)[0],
            "symbol": symbol,
        }
        for index, (gene_id, symbol) in enumerate(zip(gene_ids, gene_symbols))
    ]

    reconciliation = corgi_cross_check(counts, units, gene_ids, arguments.corgi_root)

    primary_rows = [r for r in unit_rows if r["analysis_role"] == "primary"]
    secondary_rows = [r for r in unit_rows if r["analysis_role"] == "secondary"]
    masked_by_lineage = {
        lineage: sum(
            1
            for r in unit_rows
            if r["lineage_id"] == lineage and r["rna_state"] == "insufficient_cells"
        )
        for lineage in LINEAGES
    }

    arguments.output.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(
        arguments.output / "donor_lineage_counts.npz",
        counts=counts,
        unit_observed=np.asarray(
            [r["rna_state"] == "observed" for r in unit_rows], dtype=np.bool_
        ),
        nuclei=np.asarray([r["nuclei"] for r in unit_rows], dtype=np.int64),
    )
    write_tsv(arguments.output / "units.tsv", unit_rows)
    write_tsv(arguments.output / "genes.tsv", gene_rows)

    audit = {
        "schema_version": "masld-bench-gse296875-donor-lineage-pseudobulk-v1",
        "dataset_id": "gse296875",
        "biological_unit": "donor_by_lineage_pseudobulk",
        "unit_of_inference": "donor",
        "inferential_n": EXPECTED_DONORS,
        "unit_count_is_not_a_sample_size": True,
        "donors": EXPECTED_DONORS,
        "nuclei": len(membership),
        "wells": len(raw_sources),
        "genes": EXPECTED_RNA_FEATURES,
        "counts": "raw_cell_ranger_arc_gene_expression_umis",
        "decontaminated": False,
        "lineage_partition": "frozen_fragment_membership",
        "primary_lineages": list(PRIMARY_LINEAGES),
        "secondary_lineages": list(SECONDARY_LINEAGES),
        "secondary_lineages_are_not_confirmatory_family_members": True,
        "primary_nuclei": sum(r["nuclei"] for r in primary_rows),
        "secondary_nuclei": sum(r["nuclei"] for r in secondary_rows),
        "minimum_nuclei_per_unit": MINIMUM_CELLS_PER_UNIT,
        "primary_units_total": len(primary_rows),
        "primary_units_observed": sum(
            1 for r in primary_rows if r["rna_state"] == "observed"
        ),
        "primary_units_masked": sum(
            1 for r in primary_rows if r["rna_state"] == "insufficient_cells"
        ),
        "secondary_units_total": len(secondary_rows),
        "secondary_units_observed": sum(
            1 for r in secondary_rows if r["rna_state"] == "observed"
        ),
        "masked_by_lineage": masked_by_lineage,
        "masked_state": "insufficient_cells",
        "missing_encoded_as_biological_zero": False,
        "masked_units_retain_observed_cell_count": True,
        "rna_feature_axis_sha256": feature_hash,
        "reconciliation": reconciliation,
        "outcomes_read": False,
        "model_training_activated": False,
        "labels_are_evaluator_only": True,
        "status": "pass_donor_lineage_pseudobulk",
    }
    (arguments.output / "audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
