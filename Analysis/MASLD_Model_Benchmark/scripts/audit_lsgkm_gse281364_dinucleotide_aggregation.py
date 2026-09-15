#!/usr/bin/env python3
"""Independently rederive the LS-GKM raw aggregation without control or outcomes."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping


ROW_FIELDS = (
    "seed", "row_hash", "unit_hash", "block_hash", "stratum", "outer_fold",
    "study_id", "assay_context_id", "element_id", "source_locus_group_id",
    "long_range_block_id",
)
OUTPUT_FIELDS = (
    *ROW_FIELDS, "model_id", "head_id", "prediction", "prediction_unit",
    "calibration_status", "context_specific_prediction", "experimental_replicates",
    "biological_donors", "outcome_role",
)
SCORE_FIELDS = (
    "split_id", "model_seed", "model_state_id", "element_id",
    "source_locus_group_id", "long_range_block_id", "outer_fold", "contig",
    "variant_pos0", "ref", "alt", "direct_gkmsvm_alt_minus_ref",
    "deltasvm_alt_minus_ref",
)
READOUT_FIELD = {
    "direct_gkmsvm": "direct_gkmsvm_alt_minus_ref",
    "deltasvm": "deltasvm_alt_minus_ref",
}


class AggregationAuditError(RuntimeError):
    """Raised when independent aggregation verification fails."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise AggregationAuditError(f"JSON object differs: {path}")
    return value


def project_file(root: Path, relative_text: str, expected: str) -> Path:
    relative = Path(relative_text)
    if relative.is_absolute() or ".." in relative.parts:
        raise AggregationAuditError("unsafe project-relative path")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    if path.is_symlink() or not path.is_file() or sha256(path) != expected:
        raise AggregationAuditError(f"file authority differs: {relative_text}")
    return path


def manifest(root: Path, authority: Mapping[str, Any]) -> tuple[Path, dict[str, dict[str, Any]]]:
    path = project_file(root, authority["artifacts_path"], authority["artifacts_sha256"])
    value = load_json(path)
    roster = {
        row["path"]: row
        for row in value.get("artifacts", [])
        if isinstance(row, dict) and isinstance(row.get("path"), str)
    }
    if value.get("schema_version") != "masld-bench-artifacts-v1" or len(roster) != len(value.get("artifacts", [])):
        raise AggregationAuditError("artifact manifest differs")
    return path, roster


def member(authority: Path, roster: Mapping[str, Mapping[str, Any]], relative: str) -> Path:
    record = roster.get(relative)
    if record is None or not isinstance(record.get("sha256"), str):
        raise AggregationAuditError(f"member absent: {relative}")
    path = (authority.parent / relative).resolve(strict=True)
    path.relative_to(authority.parent)
    if path.is_symlink() or not path.is_file() or sha256(path) != record["sha256"]:
        raise AggregationAuditError(f"member bytes differ: {relative}")
    return path


def tsv(path: Path, fields: tuple[str, ...], compressed: bool) -> list[dict[str, str]]:
    opener = gzip.open if compressed else Path.open
    if compressed:
        handle_context = opener(path, "rt", encoding="utf-8", newline="")
    else:
        handle_context = opener(path, "r", encoding="utf-8", newline="")
    with handle_context as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != fields:
            raise AggregationAuditError(f"TSV header differs: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise AggregationAuditError(f"TSV empty: {path}")
    return rows


def audit(root: Path, config_path: Path, aggregation_artifacts: Path, aggregation_sha: str, output: Path) -> dict[str, Any]:
    config = load_json(config_path)
    if config.get("status") != "aggregate_25_frozen_fits_without_control_or_outcome_values":
        raise AggregationAuditError("aggregation config differs")
    if output.exists():
        raise AggregationAuditError("audit output exists")
    output.mkdir(parents=True, mode=0o750)
    aggregation_relative = aggregation_artifacts.resolve(strict=True).relative_to(root).as_posix()
    aggregation_authority = project_file(root, aggregation_relative, aggregation_sha)
    aggregation_manifest = load_json(aggregation_authority)
    aggregation_roster = {
        row["path"]: row
        for row in aggregation_manifest.get("artifacts", [])
        if isinstance(row, dict) and isinstance(row.get("path"), str)
    }
    raw_path = member(aggregation_authority, aggregation_roster, "aggregation/raw_predictions.tsv.gz")
    receipt_path = member(aggregation_authority, aggregation_roster, "aggregation/aggregation_receipt.json")
    aggregate_receipt = load_json(receipt_path)
    if (
        aggregate_receipt.get("status") != "pass_outcome_blind_raw_score_aggregation"
        or aggregate_receipt.get("raw_predictions_sha256") != sha256(raw_path)
        or aggregate_receipt.get("shared_control_prediction_member_opened") is not False
        or aggregate_receipt.get("outcomes_read") is not False
    ):
        raise AggregationAuditError("aggregation receipt differs")

    row_authority, row_roster = manifest(root, config["row_universe"])
    row_path = member(row_authority, row_roster, config["row_universe"]["member"])
    if sha256(row_path) != config["row_universe"]["member_sha256"]:
        raise AggregationAuditError("row member differs")
    rows = tsv(row_path, ROW_FIELDS, False)
    if len(rows) != 10330:
        raise AggregationAuditError("row denominator differs")

    scores: dict[tuple[int, str], dict[str, str]] = {}
    fit_authorities = [
        (config["scale_fit"], "donor0_genomic0", (1103,), True),
        *[
            (bundle, bundle["split_id"], tuple(bundle["model_seeds"]), False)
            for bundle in config["bundle_fits"]
        ],
    ]
    fit_count = 0
    for authority_config, split_id, seeds, scale in fit_authorities:
        authority, roster = manifest(root, authority_config)
        for seed in seeds:
            prefix = "probe" if scale else f"fits/{split_id}.seed{seed}"
            score_path = member(authority, roster, f"{prefix}/predictions/held_element_scores.tsv.gz")
            for row in tsv(score_path, SCORE_FIELDS, True):
                key = (seed, row["element_id"])
                if key in scores or row["split_id"] != split_id or int(row["model_seed"]) != seed:
                    raise AggregationAuditError("source score identity differs")
                scores[key] = row
            fit_count += 1
    if fit_count != 25 or len(scores) != 5165:
        raise AggregationAuditError("source score rectangle differs")

    prediction_rows = tsv(raw_path, OUTPUT_FIELDS, True)
    identities: set[tuple[str, int, str]] = set()
    context_values: dict[tuple[str, int, str], dict[str, float]] = {}
    row_lookup = {(int(row["seed"]), row["row_hash"]): row for row in rows}
    for row in prediction_rows:
        seed = int(row["seed"])
        identity = (row["head_id"], seed, row["row_hash"])
        authority = row_lookup.get((seed, row["row_hash"]))
        source = scores.get((seed, row["element_id"]))
        field = READOUT_FIELD.get(row["head_id"])
        try:
            predicted = float(row["prediction"])
            expected = float(source[field]) if source is not None and field is not None else math.nan
        except ValueError as error:
            raise AggregationAuditError("prediction is not numeric") from error
        if (
            identity in identities
            or authority is None
            or source is None
            or field is None
            or not math.isfinite(predicted)
            or predicted != expected
            or any(row[name] != authority[name] for name in ROW_FIELDS if name != "seed")
            or row["model_id"] != "lsgkm_exact_dinucleotide_null"
            or row["calibration_status"] != "uncalibrated_static_rank_score"
            or row["context_specific_prediction"] != "false"
        ):
            raise AggregationAuditError("raw aggregation row differs")
        identities.add(identity)
        context_values.setdefault((row["head_id"], seed, row["element_id"]), {})[
            row["assay_context_id"]
        ] = predicted
    if len(prediction_rows) != 20660 or len(identities) != 20660:
        raise AggregationAuditError("raw prediction coverage differs")
    if len(context_values) != 10330 or any(
        set(values) != {"HepG2_control", "HepG2_PAOA"}
        or values["HepG2_control"] != values["HepG2_PAOA"]
        for values in context_values.values()
    ):
        raise AggregationAuditError("static context-copy invariant differs")

    receipt = {
        "schema_version": "masld-bench-lsgkm-gse281364-dinucleotide-aggregation-audit-v1",
        "status": "pass_independent_raw_aggregation_rederivation",
        "aggregation_artifacts_sha256": aggregation_sha,
        "raw_predictions_sha256": sha256(raw_path),
        "fit_count": fit_count,
        "seed_by_element_scores": len(scores),
        "prediction_rows": len(prediction_rows),
        "context_copy_groups": len(context_values),
        "control_prediction_member_opened": False,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "metrics_calculated": False,
        "sealed_assets_read": False,
    }
    (output / "aggregation_audit.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--aggregation-artifacts", type=Path, required=True)
    parser.add_argument("--aggregation-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(
        args.project_root.resolve(strict=True),
        args.config.resolve(strict=True),
        args.aggregation_artifacts,
        args.aggregation_artifacts_sha256,
        args.output,
    )


if __name__ == "__main__":
    main()
