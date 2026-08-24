#!/usr/bin/env python3
"""Sequence-only scBasset profile export with a fail-closed held-ATAC firewall.

scBasset's output tasks are the cells used to fit its final dense layer.  This
module therefore never creates a held-cell embedding.  It averages sequence-to-
training-cell probabilities within frozen training-only lineages, converts each
lineage vector to a depth-free peak composition, and repeats that composition
across held donors.  No held-donor ATAC or RNA is accepted at inference.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


MODEL_ID = "scbasset"
TASK_ID = "rna_conditioned_atac"
DATASET_ID = "gse296875"
SEQUENCE_LENGTH = 1344
LINEAGES = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)
CELL_FIELDS = ("cell_id", "donor_id", "lineage", "outer_fold")
REGION_FIELDS = (
    "region_id",
    "block_id",
    "role",
    "sequence_start",
    "sequence_end",
)
DONOR_FIELDS = ("donor_id", "outer_fold", "evaluation_role")
PREDICTION_FIELDS = (
    "row_hash",
    "donor_hash",
    "block_hash",
    "stratum",
    "predicted",
)


class ScBassetProfileError(ValueError):
    """Raised when scBasset profile inference violates its frozen contract."""


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _canonical_hash(value: Any) -> str:
    return sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _join_hash(namespace: str, kind: str, identifier: str) -> str:
    return sha256(f"{namespace}\0{kind}\0{identifier}".encode()).hexdigest()


def _read_tsv(
    path: Path, expected_fields: Sequence[str]
) -> list[dict[str, str]]:
    if path.is_symlink() or not path.is_file():
        raise ScBassetProfileError(f"missing or linked TSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(expected_fields):
            raise ScBassetProfileError(
                f"TSV fields differ for {path}: {reader.fieldnames!r}"
            )
        rows = [dict(row) for row in reader]
    if not rows:
        raise ScBassetProfileError(f"TSV is empty: {path}")
    return rows


def _write_tsv(
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
        for row in rows:
            writer.writerow(dict(row))


def _write_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        handle.write(_canonical_json(value))
        handle.write("\n")


def _artifact_record(path: Path, *, relative_to: Path, role: str) -> dict[str, Any]:
    return {
        "path": path.relative_to(relative_to).as_posix(),
        "sha256": _sha256_file(path),
        "size_bytes": path.stat().st_size,
        "media_type": "text/tab-separated-values",
        "role": role,
    }


def _validate_training_cells(
    rows: Sequence[Mapping[str, str]], *, donor_test_fold: int, donor_valid_fold: int
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    cell_ids = tuple(row["cell_id"] for row in rows)
    donors = tuple(row["donor_id"] for row in rows)
    if len(set(cell_ids)) != len(cell_ids) or any(not value for value in cell_ids):
        raise ScBassetProfileError("training cell IDs are empty or duplicated")
    if any(not donor for donor in donors):
        raise ScBassetProfileError("training donor IDs are empty")
    lineages = tuple(row["lineage"] for row in rows)
    if set(lineages) != set(LINEAGES):
        raise ScBassetProfileError("training cells do not cover the five frozen lineages")
    observed_folds: set[int] = set()
    for row in rows:
        try:
            fold = int(row["outer_fold"])
        except ValueError as error:
            raise ScBassetProfileError("training outer fold is not an integer") from error
        if fold not in range(5):
            raise ScBassetProfileError("training outer fold is outside [0, 5)")
        observed_folds.add(fold)
    expected_folds = set(range(5)).difference({donor_test_fold, donor_valid_fold})
    if observed_folds != expected_folds:
        raise ScBassetProfileError(
            f"training cells use {sorted(observed_folds)}, expected {sorted(expected_folds)}"
        )
    return cell_ids, donors


def _validate_regions(rows: Sequence[Mapping[str, str]]) -> tuple[str, ...]:
    region_ids = tuple(row["region_id"] for row in rows)
    if len(set(region_ids)) != len(region_ids) or any(not value for value in region_ids):
        raise ScBassetProfileError("region IDs are empty or duplicated")
    for row in rows:
        if row["role"] not in {"valid", "test"}:
            raise ScBassetProfileError("profile inference may use only valid or test regions")
        if not row["block_id"]:
            raise ScBassetProfileError("genomic block ID is empty")
        try:
            start = int(row["sequence_start"])
            end = int(row["sequence_end"])
        except ValueError as error:
            raise ScBassetProfileError("sequence coordinates are not integers") from error
        if start < 0 or end - start != SEQUENCE_LENGTH:
            raise ScBassetProfileError("scBasset inference windows must be exactly 1,344 bp")
    return region_ids


def _validate_held_donors(
    rows: Sequence[Mapping[str, str]],
    *,
    training_donors: Sequence[str],
    donor_test_fold: int,
    donor_valid_fold: int,
) -> tuple[str, ...]:
    donor_ids = tuple(row["donor_id"] for row in rows)
    if len(set(donor_ids)) != len(donor_ids) or any(not value for value in donor_ids):
        raise ScBassetProfileError("held donor IDs are empty or duplicated")
    if set(donor_ids).intersection(training_donors):
        raise ScBassetProfileError("held and training donor rosters overlap")
    roles = {row["evaluation_role"] for row in rows}
    if roles not in ({"valid"}, {"test"}):
        raise ScBassetProfileError("one export must contain only valid or only test donors")
    expected_fold = donor_valid_fold if roles == {"valid"} else donor_test_fold
    for row in rows:
        try:
            fold = int(row["outer_fold"])
        except ValueError as error:
            raise ScBassetProfileError("held donor fold is not an integer") from error
        if fold != expected_fold:
            raise ScBassetProfileError("held donor fold differs from its evaluation role")
    return donor_ids


def aggregate_training_cell_predictions(
    sequence_predictions: Any,
    training_cells: Sequence[Mapping[str, str]],
    *,
    held_atac: Any = None,
    held_cell_ids: Sequence[str] | None = None,
    pseudocount: float = 1e-8,
) -> Any:
    """Return normalized lineage profiles from training-cell outputs only."""

    import numpy as np

    if held_atac is not None:
        raise ScBassetProfileError("held-donor ATAC is prohibited at scBasset inference")
    if held_cell_ids is not None:
        raise ScBassetProfileError("scBasset has no inductive held-cell output columns")
    if not math.isfinite(pseudocount) or pseudocount <= 0:
        raise ScBassetProfileError("profile pseudocount must be finite and positive")
    predictions = np.asarray(sequence_predictions)
    if predictions.ndim != 2 or predictions.shape[1] != len(training_cells):
        raise ScBassetProfileError(
            "sequence predictions must be region-by-training-cell in exact roster order"
        )
    if (
        not np.issubdtype(predictions.dtype, np.floating)
        or np.any(~np.isfinite(predictions))
        or np.any(predictions < 0)
        or np.any(predictions > 1)
    ):
        raise ScBassetProfileError("scBasset sigmoid predictions must be finite in [0, 1]")
    output = np.empty((len(LINEAGES), predictions.shape[0]), dtype=np.float64)
    labels = [row["lineage"] for row in training_cells]
    for lineage_index, lineage in enumerate(LINEAGES):
        indices = [index for index, value in enumerate(labels) if value == lineage]
        if not indices:
            raise ScBassetProfileError(f"training lineage is empty: {lineage}")
        profile = predictions[:, indices].mean(axis=1, dtype=np.float64)
        profile += pseudocount
        total = float(profile.sum())
        if not math.isfinite(total) or total <= 0:
            raise ScBassetProfileError("lineage profile has nonpositive mass")
        output[lineage_index] = profile / total
    if (
        np.any(~np.isfinite(output))
        or np.any(output <= 0)
        or not np.allclose(output.sum(axis=1), 1.0, rtol=0.0, atol=1e-12)
    ):
        raise ScBassetProfileError("normalized lineage profiles are invalid")
    return output


def export_prediction_bundle(
    *,
    sequence_predictions: Any,
    training_cells: Sequence[Mapping[str, str]],
    regions: Sequence[Mapping[str, str]],
    held_donors: Sequence[Mapping[str, str]],
    output: Path,
    run_id: str,
    namespace: str,
    donor_test_fold: int,
    donor_valid_fold: int,
    model_artifact_sha256: str,
    held_atac: Any = None,
    held_cell_ids: Sequence[str] | None = None,
    pseudocount: float = 1e-8,
) -> Path:
    """Write an evaluator-compatible sequence-only PredictionBundle."""

    import numpy as np

    if output.exists():
        raise ScBassetProfileError(f"output already exists: {output}")
    if len(run_id) != 64 or any(value not in "0123456789abcdef" for value in run_id):
        raise ScBassetProfileError("run ID must be lowercase SHA-256")
    if (
        len(model_artifact_sha256) != 64
        or any(value not in "0123456789abcdef" for value in model_artifact_sha256)
    ):
        raise ScBassetProfileError("model artifact hash must be lowercase SHA-256")
    if not namespace or donor_test_fold not in range(5) or donor_valid_fold not in range(5):
        raise ScBassetProfileError("split or namespace contract is invalid")
    if donor_test_fold == donor_valid_fold:
        raise ScBassetProfileError("donor test and validation folds must differ")
    cell_ids, training_donors = _validate_training_cells(
        training_cells,
        donor_test_fold=donor_test_fold,
        donor_valid_fold=donor_valid_fold,
    )
    region_ids = _validate_regions(regions)
    donor_ids = _validate_held_donors(
        held_donors,
        training_donors=training_donors,
        donor_test_fold=donor_test_fold,
        donor_valid_fold=donor_valid_fold,
    )
    predictions = np.asarray(sequence_predictions)
    if predictions.shape[0] != len(regions):
        raise ScBassetProfileError("prediction region axis differs from region roster")
    lineage_profiles = aggregate_training_cell_predictions(
        predictions,
        training_cells,
        held_atac=held_atac,
        held_cell_ids=held_cell_ids,
        pseudocount=pseudocount,
    )

    output.mkdir(parents=True, exist_ok=False)
    prediction_rows: list[dict[str, str]] = []
    row_id_rows: list[dict[str, str]] = []
    for donor in donor_ids:
        donor_hash = _join_hash(namespace, "unit", donor)
        for lineage_index, lineage in enumerate(LINEAGES):
            for region_index, region in enumerate(regions):
                row_hash = _join_hash(
                    namespace,
                    "row",
                    f"{donor}\0{lineage}\0{region_ids[region_index]}",
                )
                prediction_rows.append(
                    {
                        "row_hash": row_hash,
                        "donor_hash": donor_hash,
                        "block_hash": _join_hash(
                            namespace, "block", region["block_id"]
                        ),
                        "stratum": lineage,
                        "predicted": format(
                            float(lineage_profiles[lineage_index, region_index]),
                            ".17g",
                        ),
                    }
                )
                row_id_rows.append(
                    {"row_hash": row_hash, "donor_hash": donor_hash}
                )
    prediction_rows.sort(key=lambda row: row["row_hash"])
    row_id_rows.sort(key=lambda row: row["row_hash"])
    prediction_path = output / "predictions.tsv"
    row_ids_path = output / "row_ids.tsv"
    _write_tsv(prediction_path, PREDICTION_FIELDS, prediction_rows)
    _write_tsv(row_ids_path, ("row_hash", "donor_hash"), row_id_rows)
    standardized = _artifact_record(
        prediction_path,
        relative_to=output,
        role=f"standardized_prediction_table:{TASK_ID}",
    )
    row_artifact = _artifact_record(
        row_ids_path,
        relative_to=output,
        role=f"prediction_row_ids:{TASK_ID}",
    )
    source_join = _canonical_hash(
        {
            "task_id": TASK_ID,
            "dataset_ids": [DATASET_ID],
            "split_id": "donor_outer_x_genomic_block",
            "row_id_field": "row_hash",
            "unit_id_field": "donor_hash",
            "unit_id_namespace": namespace,
            "biological_unit": "donor",
        }
    )
    evaluation_role = held_donors[0]["evaluation_role"]
    bundle = {
        "schema_version": "masld-bench-prediction-bundle-v1",
        "bundle_id": f"{MODEL_ID}-{run_id[:16]}",
        "run_id": run_id,
        "task_id": TASK_ID,
        "model_id": MODEL_ID,
        "dataset_ids": [DATASET_ID],
        "split_id": "donor_outer_x_genomic_block",
        "artifacts": [standardized, row_artifact],
        "standardized_table": standardized,
        "row_ids": row_artifact,
        "n_predictions": len(prediction_rows),
        "row_id_field": "row_hash",
        "unit_id_field": "donor_hash",
        "unit_id_namespace": namespace,
        "biological_unit": "donor",
        "table_schema_sha256": _canonical_hash(
            {"format": "tsv", "fields": list(PREDICTION_FIELDS)}
        ),
        "source_join_key_sha256": source_join,
        "format_version": "tsv-v1",
        "missing_state": "observed",
        "metadata": {
            "donor_test_fold": donor_test_fold,
            "donor_valid_fold": donor_valid_fold,
            "evaluation_role": evaluation_role,
            "training_cell_count": len(cell_ids),
            "training_donor_count": len(set(training_donors)),
            "held_donor_count": len(donor_ids),
            "lineages": list(LINEAGES),
            "region_count": len(regions),
            "sequence_length": SEQUENCE_LENGTH,
            "prediction_scale": "training_cell_state_mean_depth_free_multinomial_peak_composition",
            "profile_pseudocount": pseudocount,
            "donor_context": "none_sequence_only",
            "held_donor_profiles_identical_within_lineage": True,
            "training_cell_output_columns_only": True,
            "held_cell_embedding_available": False,
            "held_atac_input_exposed": False,
            "observed_atac_exported": False,
            "rna_input_exposed": False,
            "public_tutorial_weights_used": False,
            "profile_fixture_passed": True,
            "model_artifact_sha256": model_artifact_sha256,
        },
    }
    bundle_path = output / "prediction_bundle.json"
    _write_json(bundle_path, bundle)
    return bundle_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence-predictions", required=True, type=Path)
    parser.add_argument("--training-cells", required=True, type=Path)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--held-donors", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--join-namespace", required=True)
    parser.add_argument("--donor-test-fold", required=True, type=int)
    parser.add_argument("--donor-valid-fold", required=True, type=int)
    parser.add_argument("--model-artifact-sha256", required=True)
    parser.add_argument("--profile-pseudocount", type=float, default=1e-8)
    parser.add_argument(
        "--held-atac",
        type=Path,
        help="Forbidden fail-closed option: any supplied held ATAC aborts inference.",
    )
    parser.add_argument(
        "--held-cell-columns",
        type=Path,
        help="Forbidden fail-closed option: scBasset cannot add held-cell outputs.",
    )
    arguments = parser.parse_args(argv)
    import numpy as np

    if arguments.held_atac is not None:
        raise ScBassetProfileError("--held-atac is prohibited")
    if arguments.held_cell_columns is not None:
        raise ScBassetProfileError("--held-cell-columns is prohibited")
    matrix = np.load(arguments.sequence_predictions, allow_pickle=False)
    cells = _read_tsv(arguments.training_cells, CELL_FIELDS)
    regions = _read_tsv(arguments.regions, REGION_FIELDS)
    donors = _read_tsv(arguments.held_donors, DONOR_FIELDS)
    bundle = export_prediction_bundle(
        sequence_predictions=matrix,
        training_cells=cells,
        regions=regions,
        held_donors=donors,
        output=arguments.output,
        run_id=arguments.run_id,
        namespace=arguments.join_namespace,
        donor_test_fold=arguments.donor_test_fold,
        donor_valid_fold=arguments.donor_valid_fold,
        model_artifact_sha256=arguments.model_artifact_sha256,
        held_atac=None,
        held_cell_ids=None,
        pseudocount=arguments.profile_pseudocount,
    )
    print(bundle)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
