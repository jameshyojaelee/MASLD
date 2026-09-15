#!/usr/bin/env python3
"""Consolidate label-free GSE105127 RSEM outputs without normalization."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


class GSE105127ConsolidationError(RuntimeError):
    """Raised when participant-zone RNA quantifications cannot be consolidated."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_plan(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        forbidden = {"phenotype", "label", "outcome", "disease", "fibrosis", "nas", "sex", "age", "bmi", "outer_fold"}
        if forbidden & set(reader.fieldnames or ()):
            raise GSE105127ConsolidationError("RNA plan violates the label firewall")
        rows = list(reader)
    if len(rows) != 57 or len({row["row_id"] for row in rows}) != 57:
        raise GSE105127ConsolidationError("RNA plan row census differs")
    return sorted(rows, key=lambda row: (row["participant_group_id"], row["zone"]))


def read_quant(path: Path) -> tuple[list[tuple[str, str]], np.ndarray]:
    features: list[tuple[str, str]] = []
    values: list[tuple[float, float, float, float]] = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"gene_id", "transcript_id(s)", "length", "effective_length", "expected_count", "TPM"} <= set(reader.fieldnames or ()):
            raise GSE105127ConsolidationError("RSEM gene-results schema differs")
        for row in reader:
            features.append((row["gene_id"], row["transcript_id(s)"]))
            try:
                item = tuple(float(row[key]) for key in ("expected_count", "TPM", "effective_length", "length"))
            except ValueError as error:
                raise GSE105127ConsolidationError("RSEM value is nonnumeric") from error
            if not all(math.isfinite(value) and value >= 0 for value in item):
                raise GSE105127ConsolidationError("RSEM value is invalid")
            values.append(item)
    if not features or len(set(features)) != len(features):
        raise GSE105127ConsolidationError("RSEM feature axis is empty or duplicated")
    return features, np.asarray(values, dtype=np.float64)


def consolidate(*, plan_root: Path, quant_root: Path, output: Path) -> dict[str, object]:
    if output.exists():
        raise GSE105127ConsolidationError(f"refusing to overwrite consolidation: {output}")
    verify_frozen_tree(plan_root)
    rows = read_plan(plan_root / "rna_rows.tsv")
    matrices = []
    features: list[tuple[str, str]] | None = None
    receipts = []
    for row in rows:
        source = quant_root / "participants" / row["row_id"]
        verify_frozen_tree(source)
        receipt = json.loads((source / "receipt.json").read_text(encoding="utf-8"))
        if (
            receipt.get("status") != "passed_raw_scale"
            or receipt.get("row_id") != row["row_id"]
            or receipt.get("participant_group_id") != row["participant_group_id"]
            or receipt.get("zone") != row["zone"]
            or receipt.get("labels_accessed") is not False
            or receipt.get("fit_or_score_performed") is not False
        ):
            raise GSE105127ConsolidationError("RNA participant receipt differs")
        current_features, current = read_quant(source / "rsem.genes.results.gz")
        if features is None:
            features = current_features
        elif current_features != features:
            raise GSE105127ConsolidationError("RSEM feature order differs across rows")
        matrices.append(current)
        receipts.append({
            "row_id": row["row_id"],
            "source_artifacts_sha256": sha256_file(source / "ARTIFACTS.json"),
            "source_gene_results_sha256": receipt["gene_results_sha256"],
        })
    assert features is not None
    matrix = np.stack(matrices, axis=0)
    if matrix.shape != (57, len(features), 4) or not np.isfinite(matrix).all() or np.any(matrix < 0):
        raise GSE105127ConsolidationError("consolidated RNA tensor differs")
    if not np.any(matrix[:, :, 0] != np.floor(matrix[:, :, 0])):
        raise GSE105127ConsolidationError("continuous expected counts were lost")
    output.mkdir(parents=True)
    names = ("expected_counts", "tpm", "effective_lengths", "lengths")
    for index, name in enumerate(names):
        np.save(output / f"{name}.npy", matrix[:, :, index], allow_pickle=False)
    with (output / "row_axis.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("row_index", "row_id", "participant_group_id", "zone", "pairing_topology"),
            delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            {
                "row_index": index,
                "row_id": row["row_id"],
                "participant_group_id": row["participant_group_id"],
                "zone": row["zone"],
                "pairing_topology": "adjacent_section",
            }
            for index, row in enumerate(rows)
        )
    with (output / "gene_axis.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("gene_index", "gencode_v49_gene_id", "transcript_ids"))
        for index, (gene, transcripts) in enumerate(features):
            writer.writerow((index, gene, transcripts))
    (output / "source_receipts.json").write_text(json.dumps(receipts, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    receipt = {
        "schema_version": "masld-bench-gse105127-rna-consolidation-v1",
        "status": "passed_label_free_raw_scale",
        "participant_zone_rows": 57,
        "participants": 19,
        "genes": len(features),
        "matrix_shape": [57, len(features)],
        "measurements": list(names),
        "dtype": "float64",
        "expected_counts_continuous": True,
        "rounding_applied": False,
        "rsem_native_tpm_retained": True,
        "post_quantification_normalization_or_feature_filter_applied": False,
        "pairing_topology": "adjacent_section",
        "labels_accessed": False,
        "fit_or_score_performed": False,
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    freeze_tree(output, {"artifact_class": "gse105127_rna_continuous_matrices", "status": "passed"})
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--quant-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(consolidate(plan_root=args.plan_root, quant_root=args.quant_root, output=args.output), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
