#!/usr/bin/env python3
"""Evaluate frozen scBasset validation predictions without reading test ATAC."""

from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


LINEAGES = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)
PREDICTION_FIELDS = (
    "row_hash",
    "donor_hash",
    "block_hash",
    "stratum",
    "predicted",
)
REGION_FIELDS = (
    "region_id",
    "block_id",
    "role",
    "sequence_start",
    "sequence_end",
)
CELL_FIELDS = ("cell_id", "donor_id", "lineage", "outer_fold")
DONOR_FIELDS = ("donor_id", "outer_fold", "evaluation_role")
CCRE_FIELDS = (
    "contig",
    "output_start",
    "output_end",
    "window_id",
    "genomic_fold",
    "window_class",
    "ccre_class",
    "ccre_id",
    "ccre_start",
    "ccre_end",
    "input_start",
    "input_end",
    "selection_hash",
)
BIGWIG_FIELDS = (
    "donor_id",
    "outer_fold",
    "lineage_id",
    "analysis_role",
    "nuclei",
    "unique_fragments",
    "read_support",
    "tn5_insertions",
    "nonzero_positions",
    "max_pending_positions",
    "path",
    "size_bytes",
    "sha256",
)
EPSILON = 1.0e-6


class ScBassetEvaluationError(ValueError):
    """Raised when the validation-only evaluator requirement differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def join_hash(namespace: str, kind: str, identifier: str) -> str:
    return sha256(f"{namespace}\0{kind}\0{identifier}".encode()).hexdigest()


def read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    if path.is_symlink() or not path.is_file():
        raise ScBassetEvaluationError(f"missing or linked TSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise ScBassetEvaluationError(f"TSV fields differ: {path}")
        return [dict(row) for row in reader]


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


def deviance_per_insertion(observed: Any, predicted: Any) -> float:
    import numpy as np

    truth = np.asarray(observed, dtype=np.float64)
    score = np.asarray(predicted, dtype=np.float64)
    if truth.shape != score.shape or truth.size == 0 or np.any(truth < 0):
        raise ScBassetEvaluationError("multinomial axes or truth differ")
    total = float(truth.sum())
    if total <= 0 or np.any(~np.isfinite(score)) or np.any(score < 0):
        raise ScBassetEvaluationError("multinomial mass differs")
    score_total = float(score.sum())
    probability = (
        np.full(score.shape, 1.0 / score.size)
        if score_total <= 0
        else score / score_total
    )
    probability = (1.0 - EPSILON) * probability + EPSILON / score.size
    positive = truth > 0
    expected = total * probability[positive]
    return float(
        2.0 * np.sum(truth[positive] * np.log(truth[positive] / expected)) / total
    )


def average_precision(observed: Any, scores: Any) -> float:
    from sklearn.metrics import average_precision_score

    labels = observed > 0
    if not labels.any():
        raise ScBassetEvaluationError("AUPRC has no positive region")
    return float(average_precision_score(labels, scores))


def spearman_or_zero(observed: Any, predicted: Any) -> float:
    import numpy as np
    from scipy.stats import spearmanr

    if np.all(observed == observed[0]) or np.all(predicted == predicted[0]):
        return 0.0
    result = spearmanr(observed, predicted)
    value = float(getattr(result, "statistic", result[0]))
    return value if math.isfinite(value) else 0.0


def _extract_bigwig(arguments: Mapping[str, Any]) -> dict[str, Any]:
    import numpy as np
    import pyBigWig

    path = Path(arguments["path"])
    if (
        path.is_symlink()
        or path.stat().st_size != int(arguments["size_bytes"])
        or sha256_file(path) != arguments["sha256"]
    ):
        raise ScBassetEvaluationError("donor-lineage bigWig differs")
    track = pyBigWig.open(str(path))
    try:
        values = [
            track.stats(contig, start, end, type="sum", exact=True)[0]
            for contig, start, end in arguments["windows"]
        ]
    finally:
        track.close()
    counts = np.asarray([0.0 if value is None else value for value in values])
    rounded = np.rint(counts)
    if (
        np.any(~np.isfinite(counts))
        or np.any(counts < 0)
        or not np.allclose(counts, rounded, rtol=0.0, atol=1.0e-6)
    ):
        raise ScBassetEvaluationError("donor-lineage counts are not integers")
    return {
        "donor_id": arguments["donor_id"],
        "lineage": arguments["lineage_id"],
        "counts": rounded.astype(np.uint32),
    }


def _verify_root(root: Path, expected: str) -> Path:
    value = root.resolve(strict=True)
    if sha256_file(value / "ARTIFACTS.json") != expected:
        raise ScBassetEvaluationError(f"ARTIFACTS differs: {value}")
    return value


def _load_bundle(root: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    bundle_path = root / "predictions/bundle/prediction_bundle.json"
    bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
    metadata = bundle.get("metadata", {})
    if (
        bundle.get("task_id") != "rna_conditioned_atac"
        or bundle.get("model_id") != "scbasset"
        or bundle.get("n_predictions") != 720_000
        or metadata.get("evaluation_role") != "valid"
        or metadata.get("held_donor_count") != 9
        or metadata.get("region_count") != 16_000
        or metadata.get("held_atac_input_exposed") is not False
        or metadata.get("observed_atac_exported") is not False
        or metadata.get("rna_input_exposed") is not False
    ):
        raise ScBassetEvaluationError("prediction bundle contract differs")
    record = bundle["standardized_table"]
    path = bundle_path.parent / record["path"]
    if (
        path.stat().st_size != record["size_bytes"]
        or sha256_file(path) != record["sha256"]
    ):
        raise ScBassetEvaluationError("prediction table differs")
    rows = read_tsv(path, PREDICTION_FIELDS)
    if len(rows) != 720_000:
        raise ScBassetEvaluationError("prediction row count differs")
    return bundle, rows


def evaluate(
    *,
    prediction_root: Path,
    prediction_artifacts_sha256: str,
    input_root: Path,
    input_artifacts_sha256: str,
    donor_bigwigs: Path,
    donor_bigwigs_artifacts_sha256: str,
    split_root: Path,
    split_artifacts_sha256: str,
    output: Path,
    workers: int,
) -> dict[str, Any]:
    import numpy as np
    from scipy import sparse

    if output.exists() or not 1 <= workers <= 16:
        raise ScBassetEvaluationError("output or worker contract differs")
    prediction_root = _verify_root(prediction_root, prediction_artifacts_sha256)
    input_root = _verify_root(input_root, input_artifacts_sha256)
    donor_bigwigs = _verify_root(donor_bigwigs, donor_bigwigs_artifacts_sha256)
    split_root = _verify_root(split_root, split_artifacts_sha256)
    bundle, prediction_rows = _load_bundle(prediction_root)
    namespace = str(bundle["unit_id_namespace"])
    inputs = input_root / "inputs"
    regions = read_tsv(inputs / "valid_regions.tsv", REGION_FIELDS)
    cells = read_tsv(inputs / "training_cells.tsv", CELL_FIELDS)
    donors = read_tsv(inputs / "held_valid_donors.tsv", DONOR_FIELDS)
    if (
        len(regions) != 16_000
        or len(cells) != 27_377
        or len(donors) != 9
        or any(row["role"] != "valid" for row in regions)
        or any(row["evaluation_role"] != "valid" for row in donors)
        or any(int(row["outer_fold"]) != 1 for row in donors)
    ):
        raise ScBassetEvaluationError("validation roster differs")
    donor_ids = [row["donor_id"] for row in donors]
    donor_index = {value: index for index, value in enumerate(donor_ids)}
    lineage_index = {value: index for index, value in enumerate(LINEAGES)}
    region_index = {row["region_id"]: index for index, row in enumerate(regions)}

    split_windows = {
        row["window_id"]: row
        for row in read_tsv(split_root / "ccre_evaluation_windows.tsv", CCRE_FIELDS)
        if int(row["genomic_fold"]) == 1
    }
    if set(split_windows) != set(region_index):
        raise ScBassetEvaluationError("scBasset and split validation cCREs differ")
    for row in regions:
        window = split_windows[row["region_id"]]
        midpoint = (int(window["output_start"]) + int(window["output_end"])) // 2
        if (
            row["block_id"] != window["contig"]
            or int(row["sequence_start"]) != midpoint - 672
            or int(row["sequence_end"]) != midpoint + 672
        ):
            raise ScBassetEvaluationError("scBasset sequence/output window join differs")

    expected: dict[str, tuple[int, int, int]] = {}
    for donor in donor_ids:
        for lineage in LINEAGES:
            for region in regions:
                row_hash = join_hash(
                    namespace,
                    "row",
                    f"{donor}\0{lineage}\0{region['region_id']}",
                )
                expected[row_hash] = (
                    donor_index[donor],
                    lineage_index[lineage],
                    region_index[region["region_id"]],
                )
    candidate = np.full((9, 5, 16_000), np.nan, dtype=np.float64)
    for row in prediction_rows:
        position = expected.pop(row["row_hash"], None)
        if position is None:
            raise ScBassetEvaluationError("prediction row identity is unexpected")
        donor = donor_ids[position[0]]
        lineage = LINEAGES[position[1]]
        region = regions[position[2]]
        if (
            row["donor_hash"] != join_hash(namespace, "unit", donor)
            or row["block_hash"]
            != join_hash(namespace, "block", region["block_id"])
            or row["stratum"] != lineage
        ):
            raise ScBassetEvaluationError("prediction join fields differ")
        value = float(row["predicted"])
        if not math.isfinite(value) or value <= 0:
            raise ScBassetEvaluationError("prediction is not finite and positive")
        candidate[position] = value
    if expected or np.any(~np.isfinite(candidate)):
        raise ScBassetEvaluationError("validation prediction universe is incomplete")

    training = sparse.load_npz(inputs / "m_valid.npz").tocsr()
    if training.shape != (16_000, 27_377):
        raise ScBassetEvaluationError("training-only validation matrix axes differ")
    training_labels = np.asarray([row["lineage"] for row in cells])
    lineage_mean = np.empty((5, 16_000), dtype=np.float64)
    for lineage, index in lineage_index.items():
        positions = np.flatnonzero(training_labels == lineage)
        lineage_mean[index] = np.asarray(training[:, positions].mean(axis=1)).ravel()
    global_mean = np.asarray(training.mean(axis=1)).ravel()
    lineage_mean += 1.0e-8
    global_mean += 1.0e-8
    lineage_mean /= lineage_mean.sum(axis=1, keepdims=True)
    global_mean /= global_mean.sum()

    manifest = read_tsv(donor_bigwigs / "bigwig_manifest.tsv", BIGWIG_FIELDS)
    manifest_index = {
        (row["donor_id"], row["lineage_id"]): row
        for row in manifest
        if row["donor_id"] in donor_index and row["lineage_id"] in lineage_index
    }
    if len(manifest_index) != 45:
        raise ScBassetEvaluationError("held-valid donor-lineage bigWig census differs")
    windows = [
        (
            split_windows[row["region_id"]]["contig"],
            int(split_windows[row["region_id"]]["output_start"]),
            int(split_windows[row["region_id"]]["output_end"]),
        )
        for row in regions
    ]
    extraction = [
        {
            **row,
            "lineage_id": lineage,
            "path": str(donor_bigwigs / row["path"]),
            "windows": windows,
        }
        for (donor, lineage), row in sorted(manifest_index.items())
    ]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        extracted = list(executor.map(_extract_bigwig, extraction))
    observed = np.zeros((9, 5, 16_000), dtype=np.uint32)
    for record in extracted:
        observed[
            donor_index[record["donor_id"]], lineage_index[record["lineage"]]
        ] = record["counts"]
    if np.any(observed.sum(axis=2) <= 0):
        raise ScBassetEvaluationError("held-valid donor-lineage outcome is empty")

    block_positions = {
        block: np.asarray(
            [index for index, region in enumerate(regions) if region["block_id"] == block],
            dtype=np.int64,
        )
        for block in sorted({row["block_id"] for row in regions})
    }
    models = ("scbasset", "training_lineage_mean", "training_global_mean")
    unit_rows: list[dict[str, Any]] = []
    full_rows: list[dict[str, Any]] = []
    for donor_offset, donor in enumerate(donor_ids):
        for lineage_offset, lineage in enumerate(LINEAGES):
            truth = observed[donor_offset, lineage_offset].astype(np.float64)
            vectors = {
                "scbasset": candidate[donor_offset, lineage_offset],
                "training_lineage_mean": lineage_mean[lineage_offset],
                "training_global_mean": global_mean,
            }
            for model_id, vector in vectors.items():
                full_rows.append(
                    {
                        "model_id": model_id,
                        "donor_hash": join_hash(namespace, "unit", donor),
                        "stratum": lineage,
                        "peak_auprc": format(average_precision(truth, vector), ".17g"),
                        "profile_spearman": format(spearman_or_zero(truth, vector), ".17g"),
                    }
                )
                for block, positions in block_positions.items():
                    block_truth = truth[positions]
                    if block_truth.sum() <= 0:
                        continue
                    unit_rows.append(
                        {
                            "model_id": model_id,
                            "donor_hash": join_hash(namespace, "unit", donor),
                            "block_hash": join_hash(namespace, "block", block),
                            "stratum": lineage,
                            "observed_insertions": int(block_truth.sum()),
                            "regions": len(positions),
                            "deviance_per_insertion": format(
                                deviance_per_insertion(block_truth, vector[positions]),
                                ".17g",
                            ),
                        }
                    )

    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in unit_rows:
        grouped[(row["model_id"], row["stratum"])].append(
            float(row["deviance_per_insertion"])
        )
    model_mean = {
        model: sum(
            float(row["deviance_per_insertion"])
            for row in unit_rows
            if row["model_id"] == model
        )
        / sum(row["model_id"] == model for row in unit_rows)
        for model in models
    }
    strongest = min(
        ("training_lineage_mean", "training_global_mean"),
        key=lambda model: (model_mean[model], model),
    )
    relative = (model_mean[strongest] - model_mean["scbasset"]) / model_mean[strongest]
    lineage_gain = {}
    for lineage in LINEAGES:
        baseline = sum(grouped[(strongest, lineage)]) / len(grouped[(strongest, lineage)])
        model = sum(grouped[("scbasset", lineage)]) / len(grouped[("scbasset", lineage)])
        lineage_gain[lineage] = (baseline - model) / baseline
    output.mkdir(parents=True, mode=0o750)
    write_tsv(
        output / "donor_lineage_block_metrics.tsv",
        (
            "model_id",
            "donor_hash",
            "block_hash",
            "stratum",
            "observed_insertions",
            "regions",
            "deviance_per_insertion",
        ),
        sorted(
            unit_rows,
            key=lambda row: (
                row["model_id"], row["donor_hash"], row["stratum"], row["block_hash"]
            ),
        ),
    )
    write_tsv(
        output / "donor_lineage_secondary_metrics.tsv",
        ("model_id", "donor_hash", "stratum", "peak_auprc", "profile_spearman"),
        sorted(full_rows, key=lambda row: (row["model_id"], row["donor_hash"], row["stratum"])),
    )
    result = {
        "schema_version": "masld-bench-scbasset-valid-evaluation-v1",
        "status": "pass",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "split_id": "donor0_genomic0",
        "evaluation_role": "valid",
        "seed": 11,
        "n_donors": 9,
        "n_lineages": 5,
        "n_regions": 16_000,
        "n_genomic_blocks": len(block_positions),
        "mean_deviance_per_insertion": model_mean,
        "strongest_training_only_baseline": strongest,
        "scbasset_relative_deviance_reduction": relative,
        "lineage_relative_deviance_reduction": lineage_gain,
        "development_gate": {
            "minimum_relative_deviance_reduction": 0.05,
            "minimum_improved_lineages": 4,
            "maximum_allowed_lineage_worsening": -0.02,
            "overall_threshold_passed": relative >= 0.05,
            "improved_lineages": sum(value > 0 for value in lineage_gain.values()),
            "no_lineage_worse_than_threshold": min(lineage_gain.values()) >= -0.02,
        },
        "one_seed_one_crossed_fold_screen": True,
        "promotion_gate_finalized": False,
        "champion_claim_allowed": False,
        "biological_unit": "donor",
        "cells_used_as_independent_replicates": False,
        "prediction_frozen_before_outcomes": True,
        "held_valid_atac_read_by_evaluator_only": True,
        "held_valid_atac_exposed_to_model": False,
        "test_regions_read": False,
        "test_atac_read": False,
        "test_predictions_read": False,
        "external_evaluation": False,
        "prediction_artifacts_sha256": prediction_artifacts_sha256,
        "input_artifacts_sha256": input_artifacts_sha256,
        "donor_bigwigs_artifacts_sha256": donor_bigwigs_artifacts_sha256,
        "split_artifacts_sha256": split_artifacts_sha256,
        "smoothing_epsilon": EPSILON,
    }
    (output / "evaluation.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction-root", type=Path, required=True)
    parser.add_argument("--prediction-artifacts-sha256", required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--input-artifacts-sha256", required=True)
    parser.add_argument("--donor-bigwigs", type=Path, required=True)
    parser.add_argument("--donor-bigwigs-artifacts-sha256", required=True)
    parser.add_argument("--split-root", type=Path, required=True)
    parser.add_argument("--split-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    arguments = parser.parse_args()
    result = evaluate(**vars(arguments))
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
