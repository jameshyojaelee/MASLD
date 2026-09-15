#!/usr/bin/env python3
"""Meta-aggregate frozen scBasset validation metrics without reopening outcomes."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping


MODELS = ("scbasset", "training_lineage_mean", "training_global_mean")
BASELINES = ("training_lineage_mean", "training_global_mean")
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
EXPECTED_EVALUATOR_SHA = "61018950b7772b306935f38ef0fa074ef5552aafa8fa71ce03bfa4c162676b3c"
EXPECTED_BIGWIG_SHA = "01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f"
EXPECTED_SPLIT_SHA = "00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f"


class ScBassetAllFoldAggregateError(ValueError):
    """Raised when frozen fold results are not meta-aggregate compatible."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ScBassetAllFoldAggregateError(f"missing or linked {label}: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScBassetAllFoldAggregateError(f"invalid {label}: {path}") from error
    if not isinstance(value, dict):
        raise ScBassetAllFoldAggregateError(f"{label} must be an object")
    return value


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    if path.is_symlink() or not path.is_file():
        raise ScBassetAllFoldAggregateError(f"missing or linked TSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        return fields, [dict(row) for row in reader]


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def source_evaluator_sha(path: Path) -> str:
    records = []
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split(maxsplit=1)
        if len(parts) == 2 and parts[1].endswith("/scripts/evaluate_scbasset_valid.py"):
            records.append(parts[0])
    if records != [EXPECTED_EVALUATOR_SHA]:
        raise ScBassetAllFoldAggregateError("base evaluator source hash differs")
    return records[0]


def summarize(rows: list[Mapping[str, Any]]) -> dict[str, Any]:
    model_rows: dict[str, list[float]] = defaultdict(list)
    lineage_rows: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        model = str(row["model_id"])
        lineage = str(row["stratum"])
        value = float(row["deviance_per_insertion"])
        if model not in MODELS or lineage not in LINEAGES or not math.isfinite(value) or value < 0:
            raise ScBassetAllFoldAggregateError("invalid unit metric value")
        model_rows[model].append(value)
        lineage_rows[(model, lineage)].append(value)
    if set(model_rows) != set(MODELS):
        raise ScBassetAllFoldAggregateError("unit metric model roster differs")
    model_mean = {model: sum(model_rows[model]) / len(model_rows[model]) for model in MODELS}
    strongest = min(BASELINES, key=lambda model: (model_mean[model], model))
    relative = (model_mean[strongest] - model_mean["scbasset"]) / model_mean[strongest]
    lineage_gain = {}
    for lineage in LINEAGES:
        baseline = lineage_rows[(strongest, lineage)]
        candidate = lineage_rows[("scbasset", lineage)]
        if not baseline or len(baseline) != len(candidate):
            raise ScBassetAllFoldAggregateError("lineage unit roster differs")
        baseline_mean = sum(baseline) / len(baseline)
        candidate_mean = sum(candidate) / len(candidate)
        lineage_gain[lineage] = (baseline_mean - candidate_mean) / baseline_mean
    return {
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
    }


def aggregate(
    *,
    fold0_root: Path,
    fold0_artifacts_sha256: str,
    fold1to4_root: Path,
    fold1to4_artifacts_sha256: str,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise ScBassetAllFoldAggregateError("output already exists")
    roots = [fold0_root.resolve(strict=True), fold1to4_root.resolve(strict=True)]
    expected_hashes = [fold0_artifacts_sha256, fold1to4_artifacts_sha256]
    for root, expected in zip(roots, expected_hashes):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise ScBassetAllFoldAggregateError("source ARTIFACTS hash differs")
    fold0 = read_json(roots[0] / "evaluation/evaluation.json", "fold-0 evaluation")
    fold1to4 = read_json(roots[1] / "evaluation/evaluation.json", "fold-1-to-4 evaluation")
    source_hashes = [source_evaluator_sha(root / "source.sha256") for root in roots]
    if (
        fold0.get("schema_version") != "masld-bench-scbasset-valid-evaluation-v1"
        or fold1to4.get("schema_version")
        != "masld-bench-scbasset-crossfold-valid-evaluation-v1"
        or fold0.get("status") != "pass"
        or fold1to4.get("status") != "pass"
        or fold0.get("split_id") != "donor0_genomic0"
        or fold1to4.get("outer_folds") != [1, 2, 3, 4]
        or fold0.get("seed") != 11
        or fold1to4.get("seed") != 20260824
        or fold0.get("n_donors") != 9
        or fold1to4.get("n_donors") != 30
        or fold0.get("donor_bigwigs_artifacts_sha256") != EXPECTED_BIGWIG_SHA
        or fold1to4.get("donor_bigwigs_artifacts_sha256") != EXPECTED_BIGWIG_SHA
        or fold0.get("split_artifacts_sha256") != EXPECTED_SPLIT_SHA
        or fold1to4.get("split_artifacts_sha256") != EXPECTED_SPLIT_SHA
        or fold0.get("test_atac_read")
        or fold1to4.get("test_atac_read")
        or fold0.get("champion_claim_allowed")
        or fold1to4.get("champion_claim_allowed")
    ):
        raise ScBassetAllFoldAggregateError("evaluation compatibility contract differs")

    unit_fields = (
        "outer_fold",
        "model_id",
        "donor_hash",
        "block_hash",
        "stratum",
        "observed_insertions",
        "regions",
        "deviance_per_insertion",
    )
    combined_units: list[dict[str, Any]] = []
    donor_sets: dict[int, set[str]] = {}
    block_sets: dict[int, set[str]] = {}
    for root, expected_folds in ((roots[0], (0,)), (roots[1], (1, 2, 3, 4))):
        fields, rows = read_tsv(root / "evaluation/donor_lineage_block_metrics.tsv")
        expected_fields = unit_fields[1:] if expected_folds == (0,) else unit_fields
        if fields != expected_fields:
            raise ScBassetAllFoldAggregateError("unit metric schema differs")
        for row in rows:
            fold = 0 if expected_folds == (0,) else int(row["outer_fold"])
            if fold not in expected_folds:
                raise ScBassetAllFoldAggregateError("unit metric fold differs")
            combined_units.append({**row, "outer_fold": fold})
            donor_sets.setdefault(fold, set()).add(row["donor_hash"])
            block_sets.setdefault(fold, set()).add(row["block_hash"])

    if [len(donor_sets.get(fold, set())) for fold in range(5)] != [9, 7, 4, 8, 11]:
        raise ScBassetAllFoldAggregateError("donor fold census differs")
    all_donors = set().union(*donor_sets.values())
    if len(all_donors) != 39 or sum(map(len, donor_sets.values())) != 39:
        raise ScBassetAllFoldAggregateError("donors overlap across validation folds")
    all_blocks = set().union(*block_sets.values())
    if len(all_blocks) != sum(map(len, block_sets.values())):
        raise ScBassetAllFoldAggregateError("genomic blocks overlap across validation folds")
    grouped: dict[tuple[int, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in combined_units:
        grouped[(int(row["outer_fold"]), row["donor_hash"], row["block_hash"], row["stratum"])].append(row)
    for rows in grouped.values():
        if {row["model_id"] for row in rows} != set(MODELS) or len(rows) != len(MODELS):
            raise ScBassetAllFoldAggregateError("donor-block model roster differs")
        if len({(row["observed_insertions"], row["regions"]) for row in rows}) != 1:
            raise ScBassetAllFoldAggregateError("donor-block observed axes differ")

    secondary_fields = (
        "outer_fold",
        "model_id",
        "donor_hash",
        "stratum",
        "peak_auprc",
        "profile_spearman",
    )
    combined_secondary: list[dict[str, Any]] = []
    for root, expected_folds in ((roots[0], (0,)), (roots[1], (1, 2, 3, 4))):
        fields, rows = read_tsv(root / "evaluation/donor_lineage_secondary_metrics.tsv")
        expected_fields = secondary_fields[1:] if expected_folds == (0,) else secondary_fields
        if fields != expected_fields:
            raise ScBassetAllFoldAggregateError("secondary metric schema differs")
        for row in rows:
            fold = 0 if expected_folds == (0,) else int(row["outer_fold"])
            if fold not in expected_folds:
                raise ScBassetAllFoldAggregateError("secondary metric fold differs")
            combined_secondary.append({**row, "outer_fold": fold})
    secondary_means = {}
    for metric in ("peak_auprc", "profile_spearman"):
        secondary_means[metric] = {}
        for model in MODELS:
            values = [float(row[metric]) for row in combined_secondary if row["model_id"] == model]
            if len(values) != 39 * len(LINEAGES) or any(not math.isfinite(value) for value in values):
                raise ScBassetAllFoldAggregateError("secondary metric census differs")
            secondary_means[metric][model] = sum(values) / len(values)

    output.mkdir(mode=0o750)
    write_tsv(output / "donor_lineage_block_metrics.tsv", unit_fields, combined_units)
    write_tsv(output / "donor_lineage_secondary_metrics.tsv", secondary_fields, combined_secondary)
    result = {
        "schema_version": "masld-bench-scbasset-all5-valid-meta-aggregate-v1",
        "status": "pass",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "evaluation_role": "valid",
        "outer_folds": [0, 1, 2, 3, 4],
        "n_donors": 39,
        "n_lineages": 5,
        "n_genomic_blocks": len(all_blocks),
        "donors_by_outer_fold": {str(fold): len(donor_sets[fold]) for fold in range(5)},
        "blocks_by_outer_fold": {str(fold): len(block_sets[fold]) for fold in range(5)},
        "seed_by_outer_fold": {"0": 11, "1": 20260824, "2": 20260824, "3": 20260824, "4": 20260824},
        "mixed_seed_meta_aggregate": True,
        "fixed_seed_five_fold_cv": False,
        "unit_weighting": "equal donor x lineage x genomic-block rows",
        "cells_used_as_independent_replicates": False,
        "folds_used_as_independent_biological_replicates": False,
        "base_evaluator_sha256": source_hashes[0],
        "base_evaluator_identical_across_sources": len(set(source_hashes)) == 1,
        "source_evaluations": [
            {"outer_folds": [0], "path": roots[0].as_posix(), "artifacts_sha256": fold0_artifacts_sha256},
            {"outer_folds": [1, 2, 3, 4], "path": roots[1].as_posix(), "artifacts_sha256": fold1to4_artifacts_sha256},
        ],
        "donor_bigwigs_artifacts_sha256": EXPECTED_BIGWIG_SHA,
        "split_artifacts_sha256": EXPECTED_SPLIT_SHA,
        "frozen_metric_rows_only": True,
        "outcome_files_opened_or_recomputed": False,
        "test_regions_read": False,
        "test_atac_read": False,
        "test_predictions_read": False,
        "champion_claim_allowed": False,
        "secondary_donor_lineage_means": secondary_means,
        "fold_summaries": [{"outer_fold": 0, **{key: fold0[key] for key in (
            "n_donors", "strongest_training_only_baseline", "mean_deviance_per_insertion",
            "scbasset_relative_deviance_reduction", "lineage_relative_deviance_reduction", "development_gate",
        )}}] + fold1to4["fold_summaries"],
        **summarize(combined_units),
    }
    (output / "evaluation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fold0-root", type=Path, required=True)
    parser.add_argument("--fold0-artifacts-sha256", required=True)
    parser.add_argument("--fold1to4-root", type=Path, required=True)
    parser.add_argument("--fold1to4-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    print(json.dumps(aggregate(**vars(parser.parse_args())), sort_keys=True))


if __name__ == "__main__":
    main()
