#!/usr/bin/env python3
"""Independently verify the GSE296875 RNA cell-state transfer view."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


VIEW_ID = "gse296875_rna_cell_state_external_development_7500_v1"
SEED = 20260825
ROWS_PER_CLASS = 1_500
EXPECTED_ROWS = 7_500
EXPECTED_DONORS = 39
ROW_NAMESPACE = "gse296875:rna_cell_state_external_development:row:v1"
DONOR_NAMESPACE = "gse296875:rna_cell_state_external_development:donor:v1"
SOURCE_TO_BROAD = {
    "B cells": "immune",
    "Cholangiocytes": "cholangiocyte",
    "Hepatocytes": "hepatocyte",
    "Kupffer": "immune",
    "LSEC": "endothelial",
    "Mesenchymal": "mesenchymal_stromal",
    "NK-T": "immune",
}
CLASSES = (
    "cholangiocyte",
    "endothelial",
    "hepatocyte",
    "immune",
    "mesenchymal_stromal",
)


class VerificationError(ValueError):
    """Raised when independent verification differs from the frozen view."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def identifier(namespace: str, raw_value: str) -> str:
    return hashlib.sha256(f"{namespace}\0{raw_value}".encode()).hexdigest()


def _priority(seed: int, *parts: str) -> int:
    text = "\x1f".join((str(seed), *parts))
    return int.from_bytes(hashlib.sha256(text.encode()).digest(), "big")


def independent_select(
    mapping_rows: Sequence[Mapping[str, str]],
    label_rows: Sequence[Mapping[str, str]],
    *,
    per_class: int = ROWS_PER_CLASS,
    seed: int = SEED,
) -> list[dict[str, str]]:
    mapping = {row["cell_id"]: dict(row) for row in mapping_rows}
    if len(mapping) != len(mapping_rows) or set(mapping) != {
        row["cell_id"] for row in label_rows
    }:
        raise VerificationError("mapping and label row identities differ")
    donors = {row["donor_id"] for row in mapping_rows}
    pools: dict[str, dict[str, list[dict[str, str]]]] = {
        broad: {} for broad in CLASSES
    }
    for label_row in label_rows:
        source_label = label_row["author_label"]
        if source_label not in SOURCE_TO_BROAD:
            raise VerificationError("source label roster differs")
        value = mapping[label_row["cell_id"]]
        well = value["well_id"]
        encoded = value["raw_barcode"]
        prefix = f"{well}_"
        if not encoded.startswith(prefix):
            raise VerificationError("well-barcode namespace differs")
        record = {
            "cell_id": value["cell_id"],
            "donor_id": value["donor_id"],
            "well_id": well,
            "raw_barcode": encoded[len(prefix):],
            "source_label": source_label,
            "broad_label": SOURCE_TO_BROAD[source_label],
        }
        pools[record["broad_label"]].setdefault(record["donor_id"], []).append(record)

    selected: list[dict[str, str]] = []
    for broad in CLASSES:
        grouped = pools[broad]
        if set(grouped) != donors:
            raise VerificationError(f"{broad} donor coverage differs")
        for donor, rows in grouped.items():
            rows.sort(
                key=lambda row: (
                    _priority(seed, "cell", broad, donor, row["cell_id"]),
                    row["cell_id"],
                )
            )
        donor_order = sorted(
            grouped,
            key=lambda donor: (_priority(seed, "donor", broad, donor), donor),
        )
        chosen: list[dict[str, str]] = []
        depth = 0
        while len(chosen) < per_class:
            before = len(chosen)
            for donor in donor_order:
                if depth < len(grouped[donor]):
                    chosen.append(grouped[donor][depth])
                    if len(chosen) == per_class:
                        break
            if len(chosen) == before:
                raise VerificationError(f"{broad} selection capacity differs")
            depth += 1
        selected.extend(chosen)
    return selected


def _read_tsv(path: Path, expected: Sequence[str]) -> list[dict[str, str]]:
    if path.suffix == ".gz":
        handle_context = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        handle_context = path.open("r", encoding="utf-8", newline="")
    with handle_context as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(expected):
            raise VerificationError(f"{path.name} columns differ")
        return [dict(row) for row in reader]


def _decode(values: Any) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def verify_raw_counts(
    *,
    feature_h5: Path,
    selected: Sequence[Mapping[str, str]],
    source_lock_path: Path,
) -> tuple[int, list[str]]:
    import h5py
    import numpy as np

    lock = json.loads(source_lock_path.read_text())
    raw_sources = lock.get("raw_rna_sources")
    if not isinstance(raw_sources, dict) or set(raw_sources) != {
        f"well{index}" for index in range(1, 9)
    }:
        raise VerificationError("raw RNA source roster differs")

    with h5py.File(feature_h5, "r") as observed:
        if set(observed["obs"].keys()) != {
            "row_id", "rna_status", "rna_observed_mask", "atac_observed_mask"
        }:
            raise VerificationError("query observation fields differ")
        if observed.attrs.get("label_visibility") != "evaluator_only":
            raise VerificationError("query label firewall differs")
        output_data = observed["rna/counts_csr/data"][:]
        output_indices = observed["rna/counts_csr/indices"][:]
        output_indptr = observed["rna/counts_csr/indptr"][:]
        output_shape = tuple(int(v) for v in observed["rna/counts_csr/shape"][:])
        output_gene_ids = _decode(observed["rna/ensembl_id"][:])
        output_gene_names = _decode(observed["rna/gene_name"][:])
        observed_rows = _decode(observed["obs/row_id"][:])
        observed_rna_mask = observed["obs/rna_observed_mask"][:]
        observed_atac_mask = observed["obs/atac_observed_mask"][:]

    expected_rows = [identifier(ROW_NAMESPACE, row["cell_id"]) for row in selected]
    if observed_rows != expected_rows:
        raise VerificationError("query row order differs")
    if output_shape != (EXPECTED_ROWS, 36_601):
        raise VerificationError("query RNA shape differs")
    if not np.all(observed_rna_mask == 1) or not np.all(observed_atac_mask == 0):
        raise VerificationError("query missingness masks differ")
    if len(output_indptr) != EXPECTED_ROWS + 1 or output_indptr[-1] != len(output_data):
        raise VerificationError("query CSR structure differs")

    selected_by_well: dict[str, list[tuple[int, Mapping[str, str]]]] = {
        well: [] for well in raw_sources
    }
    for row_index, row in enumerate(selected):
        selected_by_well[row["well_id"]].append((row_index, row))

    for well in sorted(raw_sources):
        binding = raw_sources[well]
        source = (source_lock_path.parent / binding["path"]).resolve(strict=True)
        if source.stat().st_size != int(binding["bytes"]):
            raise VerificationError(f"{well} raw source size differs")
        if sha256_file(source) != binding["sha256"]:
            raise VerificationError(f"{well} raw source hash differs")
        with h5py.File(source, "r") as handle:
            matrix = handle["matrix"]
            features = matrix["features"]
            feature_types = features["feature_type"][:]
            gene_source = np.flatnonzero(feature_types == b"Gene Expression")
            gene_ids = _decode(features["id"][:][gene_source])
            gene_names = _decode(features["name"][:][gene_source])
            if gene_ids != output_gene_ids or gene_names != output_gene_names:
                raise VerificationError(f"{well} gene axis differs")
            source_to_gene = np.full(len(feature_types), -1, dtype=np.int64)
            source_to_gene[gene_source] = np.arange(len(gene_source), dtype=np.int64)
            barcodes = _decode(matrix["barcodes"][:])
            barcode_to_index = {barcode: index for index, barcode in enumerate(barcodes)}
            if len(barcode_to_index) != len(barcodes):
                raise VerificationError(f"{well} barcode axis is not unique")
            for output_row, row in selected_by_well[well]:
                column = barcode_to_index.get(row["raw_barcode"])
                if column is None:
                    raise VerificationError(f"{well} selected barcode is absent")
                source_start = int(matrix["indptr"][column])
                source_end = int(matrix["indptr"][column + 1])
                raw_indices = matrix["indices"][source_start:source_end].astype(np.int64)
                raw_values = matrix["data"][source_start:source_end]
                mapped = source_to_gene[raw_indices]
                keep = mapped >= 0
                expected_indices = mapped[keep]
                expected_values = raw_values[keep]
                observed_start = int(output_indptr[output_row])
                observed_end = int(output_indptr[output_row + 1])
                if not np.array_equal(output_indices[observed_start:observed_end], expected_indices):
                    raise VerificationError(f"query indices differ at row {output_row}")
                if not np.array_equal(output_data[observed_start:observed_end], expected_values):
                    raise VerificationError(f"query counts differ at row {output_row}")
    return len(output_data), output_gene_ids


def verify(args: argparse.Namespace) -> dict[str, Any]:
    producer = args.producer.resolve(strict=True)
    feature_h5 = producer / "view/features/rna_query.h5"
    evaluator_path = producer / "view/evaluator/labels.tsv.gz"
    mapping_rows = _read_tsv(
        args.mapping_metadata.resolve(strict=True),
        ("cell_id", "donor_id", "well_id", "raw_barcode"),
    )
    label_rows = _read_tsv(
        args.author_labels.resolve(strict=True),
        ("cell_id", "author_label"),
    )
    if len(mapping_rows) != 68_398 or len(label_rows) != 68_398:
        raise VerificationError("source row census differs")
    selected = independent_select(mapping_rows, label_rows)
    if len(selected) != EXPECTED_ROWS:
        raise VerificationError("selected row census differs")

    evaluator = _read_tsv(evaluator_path, ("row_id", "donor_id", "broad_label"))
    expected_evaluator = [
        {
            "row_id": identifier(ROW_NAMESPACE, row["cell_id"]),
            "donor_id": identifier(DONOR_NAMESPACE, row["donor_id"]),
            "broad_label": row["broad_label"],
        }
        for row in selected
    ]
    if evaluator != expected_evaluator:
        raise VerificationError("evaluator identities or labels differ")
    donors = {row["donor_id"] for row in evaluator}
    if len(donors) != EXPECTED_DONORS:
        raise VerificationError("evaluator donor census differs")
    for donor in donors:
        if {row["broad_label"] for row in evaluator if row["donor_id"] == donor} != set(CLASSES):
            raise VerificationError("evaluator donor-by-class coverage differs")

    nnz, gene_ids = verify_raw_counts(
        feature_h5=feature_h5,
        selected=selected,
        source_lock_path=args.source_lock.resolve(strict=True),
    )
    campaign = json.loads((producer / "view/campaign_receipt.json").read_text())
    if (
        campaign.get("sealed_outcomes_read") is not False
        or campaign.get("model_fit") is not False
        or campaign.get("metric_calculated") is not False
        or campaign.get("external_or_champion_claim_allowed") is not False
    ):
        raise VerificationError("campaign claim firewall differs")

    receipt = {
        "schema_version": "masld-bench-gse296875-cell-state-external-development-verification-v1",
        "view_id": VIEW_ID,
        "producer": producer.as_posix(),
        "producer_artifacts_sha256": args.expected_producer_sha256,
        "rows": EXPECTED_ROWS,
        "donors": EXPECTED_DONORS,
        "classes": list(CLASSES),
        "class_counts": {
            label: sum(row["broad_label"] == label for row in evaluator)
            for label in CLASSES
        },
        "rna_genes": len(gene_ids),
        "rna_nnz": nnz,
        "raw_counts_rederived_exactly": True,
        "selection_rederived_independently": True,
        "evaluator_join_rederived_independently": True,
        "feature_labels_present": False,
        "feature_donor_ids_present": False,
        "atac_values_read_from_query_artifact": False,
        "donor_is_biological_unit": True,
        "external_or_champion_claim_allowed": False,
        "sealed_outcomes_read": False,
        "model_fit": False,
        "metric_calculated": False,
        "status": "pass_external_development_view_only",
        "next_gate": "atlas_fit_model_inference_then_prediction_hash_commit_then_evaluator_join",
    }
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / "verification_receipt.json").open("x", encoding="utf-8") as handle:
        json.dump(receipt, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--producer", type=Path, required=True)
    parser.add_argument("--expected-producer-sha256", required=True)
    parser.add_argument("--source-lock", type=Path, required=True)
    parser.add_argument("--mapping-metadata", type=Path, required=True)
    parser.add_argument("--author-labels", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    print(json.dumps(verify(parse_args()), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
