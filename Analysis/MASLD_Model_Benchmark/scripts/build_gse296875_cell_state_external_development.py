#!/usr/bin/env python3
"""Build a label-separated GSE296875 RNA cell-state transfer view."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from scripts import build_gse296875_rna_atac_smoke as source_builder


VIEW_ID = "gse296875_rna_cell_state_external_development_7500_v1"
DATASET_ID = "gse296875"
SELECTION_SEED = 20260825
CELLS_PER_CLASS = 1_500
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


class ExternalCellStateError(ValueError):
    """Raised when the transfer-view requirements are not met."""


def canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def identifier(namespace: str, raw_value: str) -> str:
    return hashlib.sha256(f"{namespace}\0{raw_value}".encode("utf-8")).hexdigest()


def priority(seed: int, *parts: object) -> int:
    value = "\x1f".join(map(str, (seed, *parts)))
    return int.from_bytes(hashlib.sha256(value.encode("utf-8")).digest(), "big")


def select_records(
    joined: Sequence[Mapping[str, str]],
    *,
    cells_per_class: int = CELLS_PER_CLASS,
    selection_seed: int = SELECTION_SEED,
) -> list[dict[str, str]]:
    """Select an equal class census by donor-round-robin without count values."""

    if cells_per_class < 1:
        raise ExternalCellStateError("cells_per_class must be positive")
    pools: dict[str, dict[str, list[dict[str, str]]]] = {
        label: {} for label in CLASSES
    }
    for source in joined:
        source_label = str(source["source_label"])
        broad = SOURCE_TO_BROAD.get(source_label)
        if broad is None:
            raise ExternalCellStateError(f"unexpected source label: {source_label}")
        record = {str(key): str(value) for key, value in source.items()}
        record["broad_label"] = broad
        pools[broad].setdefault(record["donor_id"], []).append(record)

    donor_universe = {str(row["donor_id"]) for row in joined}
    if not donor_universe:
        raise ExternalCellStateError("joined metadata is empty")
    selected: list[dict[str, str]] = []
    for broad in CLASSES:
        donor_rows = pools[broad]
        if set(donor_rows) != donor_universe:
            raise ExternalCellStateError(
                f"{broad} is not observed in every donor"
            )
        for donor, rows in donor_rows.items():
            rows.sort(
                key=lambda row: (
                    priority(selection_seed, "cell", broad, donor, row["cell_id"]),
                    row["cell_id"],
                )
            )
        donor_order = sorted(
            donor_rows,
            key=lambda donor: (
                priority(selection_seed, "donor", broad, donor),
                donor,
            ),
        )
        class_rows: list[dict[str, str]] = []
        offset = 0
        while len(class_rows) < cells_per_class:
            progressed = False
            for donor in donor_order:
                rows = donor_rows[donor]
                if offset < len(rows):
                    class_rows.append(rows[offset])
                    progressed = True
                    if len(class_rows) == cells_per_class:
                        break
            if not progressed:
                raise ExternalCellStateError(
                    f"{broad} cannot supply {cells_per_class} cells"
                )
            offset += 1
        selected.extend(class_rows)

    if len(selected) != cells_per_class * len(CLASSES):
        raise ExternalCellStateError("selected row census differs")
    cell_ids = [row["cell_id"] for row in selected]
    if len(set(cell_ids)) != len(cell_ids):
        raise ExternalCellStateError("selected cells are not unique")
    return selected


def _write_strings(group: Any, name: str, values: Iterable[str]) -> None:
    import h5py

    group.create_dataset(
        name,
        data=list(values),
        dtype=h5py.string_dtype(encoding="utf-8"),
        compression="gzip",
    )


def write_feature_h5(
    path: Path,
    *,
    row_ids: Sequence[str],
    rna: Any,
    gene_ids: Sequence[str],
    gene_names: Sequence[str],
) -> dict[str, Any]:
    """Write only query-time RNA features; labels and donor IDs are forbidden."""

    import h5py
    import numpy as np

    if len(row_ids) != rna.shape[0] or len(set(row_ids)) != len(row_ids):
        raise ExternalCellStateError("feature row identities differ")
    if rna.shape[1] != len(gene_ids) or len(gene_ids) != len(gene_names):
        raise ExternalCellStateError("feature gene axis differs")
    with h5py.File(path, "x") as handle:
        handle.attrs["schema_version"] = "masld-bench-rna-cell-state-query-v1"
        handle.attrs["view_id"] = VIEW_ID
        handle.attrs["dataset_id"] = DATASET_ID
        handle.attrs["label_visibility"] = "evaluator_only"
        handle.attrs["native_reference"] = "10x_GRCh38_2020-A"
        handle.attrs["native_gene_annotation"] = "GENCODE_v32_filtered"
        handle.attrs["atac_status"] = "structurally_missing_for_this_task"
        obs = handle.create_group("obs")
        _write_strings(obs, "row_id", row_ids)
        _write_strings(obs, "rna_status", ("observed" for _ in row_ids))
        obs.create_dataset(
            "rna_observed_mask",
            data=np.ones(len(row_ids), dtype=np.uint8),
            compression="gzip",
        )
        obs.create_dataset(
            "atac_observed_mask",
            data=np.zeros(len(row_ids), dtype=np.uint8),
            compression="gzip",
        )
        rna_group = handle.create_group("rna")
        identity = source_builder._write_csr(rna_group.create_group("counts_csr"), rna)
        _write_strings(rna_group, "ensembl_id", gene_ids)
        _write_strings(rna_group, "gene_name", gene_names)
        rna_group.attrs["source_scale"] = "raw_integer_umi_counts"

    with h5py.File(path, "r") as handle:
        forbidden = {"donor_id", "source_label", "broad_label", "cell_id"}
        if forbidden.intersection(handle["obs"].keys()):
            raise ExternalCellStateError("feature artifact contains evaluator fields")
        if tuple(handle["rna/counts_csr/shape"][:]) != tuple(rna.shape):
            raise ExternalCellStateError("feature HDF5 failed readback")
    return identity


def _write_tsv_gzip(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    fields = ("row_id", "donor_id", "broad_label")
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)


def _write_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")


def build(args: argparse.Namespace) -> dict[str, Any]:
    source_lock_path = args.source_lock.resolve(strict=True)
    source_lock = source_builder._source_lock(source_lock_path)
    source_lock["_path"] = source_lock_path.as_posix()
    joined = source_builder._joined_metadata(args.mapping_metadata, args.author_labels)
    if set(row["source_label"] for row in joined) != set(SOURCE_TO_BROAD):
        raise ExternalCellStateError("author source-label roster differs")
    selected = select_records(joined)
    if len(selected) != EXPECTED_ROWS:
        raise ExternalCellStateError("production selection row census differs")

    atlas_ontology = json.loads(args.atlas_ontology.read_text(encoding="utf-8"))
    if (
        atlas_ontology.get("ontology_id") != "broad_liver_five_v1"
        or tuple(atlas_ontology.get("classes", {})) != CLASSES
    ):
        raise ExternalCellStateError("Atlas five-class ontology differs")

    args.output.mkdir(parents=True, exist_ok=False)
    feature_root = args.output / "features"
    evaluator_root = args.output / "evaluator"
    feature_root.mkdir()
    evaluator_root.mkdir()

    rna, gene_ids, gene_names, raw_source_records = source_builder._rna_matrix(
        selected=selected,
        source_lock=source_lock,
    )
    row_ids = [identifier(ROW_NAMESPACE, row["cell_id"]) for row in selected]
    feature_h5 = feature_root / "rna_query.h5"
    matrix_identity = write_feature_h5(
        feature_h5,
        row_ids=row_ids,
        rna=rna,
        gene_ids=gene_ids,
        gene_names=gene_names,
    )

    evaluator_rows = [
        {
            "row_id": row_id,
            "donor_id": identifier(DONOR_NAMESPACE, row["donor_id"]),
            "broad_label": row["broad_label"],
        }
        for row_id, row in zip(row_ids, selected, strict=True)
    ]
    evaluator_labels = evaluator_root / "labels.tsv.gz"
    _write_tsv_gzip(evaluator_labels, evaluator_rows)

    class_counts = {
        label: sum(row["broad_label"] == label for row in selected)
        for label in CLASSES
    }
    source_label_counts = {
        label: sum(row["source_label"] == label for row in selected)
        for label in sorted(SOURCE_TO_BROAD)
    }
    donor_ids = {row["donor_id"] for row in evaluator_rows}
    donor_class_counts = {
        donor: {
            label: sum(
                row["donor_id"] == donor and row["broad_label"] == label
                for row in evaluator_rows
            )
            for label in CLASSES
        }
        for donor in sorted(donor_ids)
    }
    if len(donor_ids) != EXPECTED_DONORS or any(
        any(value == 0 for value in counts.values())
        for counts in donor_class_counts.values()
    ):
        raise ExternalCellStateError("donor-by-class evaluator coverage differs")

    crosswalk = {
        label: sorted(
            source for source, broad in SOURCE_TO_BROAD.items() if broad == label
        )
        for label in CLASSES
    }
    ontology = {
        "schema_version": "masld-bench-cell-ontology-crosswalk-v1",
        "ontology_id": "gse296875_to_broad_liver_five_v1",
        "target_ontology_id": "broad_liver_five_v1",
        "target_ontology_sha256": source_builder.sha256_file(args.atlas_ontology),
        "source_field": "author_label",
        "classes": crosswalk,
        "labels_used_only_for_prespecified_balanced_sampling_and_evaluation": True,
        "sealed_outcomes_used": False,
    }
    ontology["ontology_sha256"] = canonical_sha256(ontology)
    _write_json(evaluator_root / "ontology.json", ontology)

    feature_receipt = {
        "schema_version": "masld-bench-rna-cell-state-query-receipt-v1",
        "view_id": VIEW_ID,
        "dataset_id": DATASET_ID,
        "role": "project_exposed_external_development",
        "rows": EXPECTED_ROWS,
        "native_genes": len(gene_ids),
        "matrix": matrix_identity,
        "row_ids_sha256": canonical_sha256(row_ids),
        "query_time_inputs": ["raw_integer_umi_counts", "rna_observed_mask"],
        "forbidden_query_time_inputs": [
            "author_cell_label",
            "broad_cell_label",
            "donor_id",
            "ATAC",
            "histology",
            "MASLD_phenotype",
        ],
        "labels_present": False,
        "donor_ids_present": False,
        "atac_values_read": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "model_fit": False,
        "metric_calculated": False,
    }
    _write_json(feature_root / "feature_receipt.json", feature_receipt)

    evaluator_receipt = {
        "schema_version": "masld-bench-rna-cell-state-evaluator-receipt-v1",
        "view_id": VIEW_ID,
        "dataset_id": DATASET_ID,
        "role": "project_exposed_external_development",
        "rows": EXPECTED_ROWS,
        "donors": EXPECTED_DONORS,
        "classes": list(CLASSES),
        "class_counts": class_counts,
        "source_label_counts": source_label_counts,
        "every_donor_has_every_class": True,
        "biological_replication_unit": "donor",
        "primary_metric": "donor_balanced_macro_f1",
        "uncertainty": "paired_donor_cluster_bootstrap_10000",
        "calibration_metrics": ["multiclass_brier", "ece"],
        "labels_exposed_to_model_adapter": False,
        "donor_ids_exposed_to_model_adapter": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "champion_claim_allowed": False,
        "external_claim_allowed": False,
        "reason": "project-exposed development cohort; donor fingerprint overlap audit not complete",
    }
    _write_json(evaluator_root / "evaluator_receipt.json", evaluator_receipt)

    feature_sha = source_builder._freeze_smoke_tree(
        feature_root,
        {
            "artifact_class": "gse296875_rna_cell_state_external_development_features",
            "view_id": VIEW_ID,
            "rows": EXPECTED_ROWS,
            "labels_present": False,
            "donor_ids_present": False,
            "sealed_outcomes_accessed": False,
        },
    )
    evaluator_sha = source_builder._freeze_smoke_tree(
        evaluator_root,
        {
            "artifact_class": "gse296875_rna_cell_state_external_development_evaluator",
            "view_id": VIEW_ID,
            "rows": EXPECTED_ROWS,
            "donors": EXPECTED_DONORS,
            "champion_claim_allowed": False,
            "sealed_outcomes_accessed": False,
        },
    )

    campaign_receipt = {
        "schema_version": "masld-bench-gse296875-cell-state-external-development-v1",
        "view_id": VIEW_ID,
        "dataset_id": DATASET_ID,
        "rows": EXPECTED_ROWS,
        "donors": EXPECTED_DONORS,
        "classes": list(CLASSES),
        "class_counts": class_counts,
        "selection_seed": SELECTION_SEED,
        "selection_policy": "class_then_donor_round_robin_then_sha256_cell_priority",
        "endpoint_labels_used_for_prespecified_balanced_sampling": True,
        "feature_artifacts_sha256": feature_sha,
        "evaluator_artifacts_sha256": evaluator_sha,
        "source_lock_sha256": source_builder.sha256_file(source_lock_path),
        "mapping_metadata_sha256": source_builder.sha256_file(args.mapping_metadata),
        "author_labels_sha256": source_builder.sha256_file(args.author_labels),
        "atlas_ontology_sha256": source_builder.sha256_file(args.atlas_ontology),
        "raw_matrix_roster_sha256": canonical_sha256(raw_source_records),
        "cohort_family_relation_to_atlas": "distinct_accession_no_alias_match",
        "donor_fingerprint_overlap_status": "not_run",
        "role": "project_exposed_external_development",
        "not_a_masld_phenotype_cohort": True,
        "external_or_champion_claim_allowed": False,
        "sealed_outcomes_read": False,
        "model_fit": False,
        "metric_calculated": False,
        "next_gate": "atlas_fit_model_inference_then_prediction_hash_commit_then_evaluator_join",
    }
    _write_json(args.output / "campaign_receipt.json", campaign_receipt)
    root_sha = source_builder._freeze_smoke_tree(
        args.output,
        {
            "artifact_class": "gse296875_rna_cell_state_external_development_campaign",
            "view_id": VIEW_ID,
            "rows": EXPECTED_ROWS,
            "donors": EXPECTED_DONORS,
            "feature_artifacts_sha256": feature_sha,
            "evaluator_artifacts_sha256": evaluator_sha,
            "external_or_champion_claim_allowed": False,
            "sealed_outcomes_accessed": False,
        },
    )
    return {**campaign_receipt, "artifacts_sha256": root_sha}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-lock", type=Path, required=True)
    parser.add_argument("--mapping-metadata", type=Path, required=True)
    parser.add_argument("--author-labels", type=Path, required=True)
    parser.add_argument("--atlas-ontology", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    receipt = build(parse_args())
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
