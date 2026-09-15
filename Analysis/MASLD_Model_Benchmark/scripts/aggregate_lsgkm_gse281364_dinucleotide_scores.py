#!/usr/bin/env python3
"""Aggregate LS-GKM held raw scores onto the outcome-blind canonical row universe."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCHEMA = "masld-bench-lsgkm-gse281364-dinucleotide-aggregation-input-v1"
STATUS = "aggregate_25_frozen_fits_without_control_or_outcome_values"
DESIGN_ID = "lsgkm_hepatocyte_atac_peak_exact_dinucleotide_null_v1"
SEEDS = (1103, 2909, 4721, 6673, 8111)
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
SPLIT_SEEDS = {
    "donor0_genomic0": (2909, 4721, 6673, 8111),
    "donor1_genomic1": SEEDS,
    "donor2_genomic2": SEEDS,
    "donor3_genomic3": SEEDS,
    "donor4_genomic4": SEEDS,
}
REUSED_FIT = ("donor0_genomic0", 1103)
ROW_FIELDS = (
    "seed",
    "row_hash",
    "unit_hash",
    "block_hash",
    "stratum",
    "outer_fold",
    "study_id",
    "assay_context_id",
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
)
SCORE_FIELDS = (
    "split_id",
    "model_seed",
    "model_state_id",
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
    "outer_fold",
    "contig",
    "variant_pos0",
    "ref",
    "alt",
    "direct_gkmsvm_alt_minus_ref",
    "deltasvm_alt_minus_ref",
)
READOUTS = {
    "direct_gkmsvm": (
        "direct_gkmsvm_alt_minus_ref",
        "native_gkmpredict_margin_ALT_minus_REF",
    ),
    "deltasvm": (
        "deltasvm_alt_minus_ref",
        "canonical_local_11mer_score_sum_ALT_minus_REF",
    ),
}
OUTPUT_FIELDS = (
    *ROW_FIELDS,
    "model_id",
    "head_id",
    "prediction",
    "prediction_unit",
    "calibration_status",
    "context_specific_prediction",
    "experimental_replicates",
    "biological_donors",
    "outcome_role",
)


class LSGKMAggregationError(RuntimeError):
    """Raised when an authority, fit rectangle, or row identity differs."""


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise LSGKMAggregationError(f"JSON object differs: {path}")
    return value


def safe_project_file(root: Path, relative_text: str, expected_sha256: str) -> Path:
    relative = Path(relative_text)
    if relative.is_absolute() or ".." in relative.parts:
        raise LSGKMAggregationError("unsafe project-relative path")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    if path.is_symlink() or not path.is_file() or file_sha256(path) != expected_sha256:
        raise LSGKMAggregationError(f"file authority differs: {relative_text}")
    return path


def bind_manifest(
    root: Path, authority: Mapping[str, Any]
) -> tuple[Path, dict[str, Any], dict[str, dict[str, Any]]]:
    path = safe_project_file(
        root, str(authority["artifacts_path"]), str(authority["artifacts_sha256"])
    )
    manifest = load_json(path)
    if manifest.get("schema_version") != "masld-bench-artifacts-v1":
        raise LSGKMAggregationError("artifact schema differs")
    roster: dict[str, dict[str, Any]] = {}
    for member in manifest.get("artifacts", []):
        if not isinstance(member, dict) or not isinstance(member.get("path"), str):
            raise LSGKMAggregationError("artifact member differs")
        if member["path"] in roster:
            raise LSGKMAggregationError("duplicate artifact member")
        roster[member["path"]] = member
    return path, manifest, roster


def bind_openable_member(
    authority_path: Path,
    roster: Mapping[str, Mapping[str, Any]],
    relative: str,
    expected_sha256: str | None = None,
) -> Path:
    member = roster.get(relative)
    if member is None or not isinstance(member.get("sha256"), str):
        raise LSGKMAggregationError(f"artifact member is absent: {relative}")
    if expected_sha256 is not None and member["sha256"] != expected_sha256:
        raise LSGKMAggregationError(f"artifact member hash differs: {relative}")
    path = (authority_path.parent / relative).resolve(strict=True)
    path.relative_to(authority_path.parent)
    if path.is_symlink() or not path.is_file() or file_sha256(path) != member["sha256"]:
        raise LSGKMAggregationError(f"artifact member bytes differ: {relative}")
    return path


def bind_manifest_record_only(
    roster: Mapping[str, Mapping[str, Any]], relative: str, expected_sha256: str
) -> dict[str, Any]:
    """Bind a member record without resolving, hashing, or opening member bytes."""
    member = roster.get(relative)
    if member is None or member.get("sha256") != expected_sha256:
        raise LSGKMAggregationError(f"manifest-only member differs: {relative}")
    return dict(member)


def read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise LSGKMAggregationError(f"TSV header differs: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise LSGKMAggregationError(f"TSV is empty: {path}")
    return rows


def read_gzip_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise LSGKMAggregationError(f"gzip TSV header differs: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise LSGKMAggregationError(f"gzip TSV is empty: {path}")
    return rows


def write_gzip_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> int:
    if path.exists():
        raise LSGKMAggregationError(f"refusing to overwrite: {path}")
    count = 0
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                for row in rows:
                    writer.writerow({field: row[field] for field in fields})
                    count += 1
    return count


def expected_remaining_fits() -> set[tuple[str, int]]:
    return {
        (split_id, seed)
        for split_id, seeds in SPLIT_SEEDS.items()
        for seed in seeds
    }


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != STATUS
        or config.get("design_id") != DESIGN_ID
    ):
        raise LSGKMAggregationError("aggregation identity differs")
    fit_grid = config.get("fit_grid", {})
    if (
        tuple(fit_grid.get("fixed_seeds", ())) != SEEDS
        or fit_grid.get("total_fit_count") != 25
        or fit_grid.get("reused_fit")
        != {"split_id": REUSED_FIT[0], "model_seed": REUSED_FIT[1]}
        or fit_grid.get("remaining_fit_count") != 24
        or fit_grid.get("expected_elements_per_seed") != 1033
        or fit_grid.get("expected_rows_per_readout") != 10330
    ):
        raise LSGKMAggregationError("fit grid differs")
    bundles = config.get("bundle_fits", [])
    bundle_map = {
        row.get("split_id"): tuple(row.get("model_seeds", ()))
        for row in bundles
        if isinstance(row, dict)
    }
    if len(bundles) != 5 or bundle_map != SPLIT_SEEDS:
        raise LSGKMAggregationError("bundle roster differs")
    firewall = config.get("action_firewall", {})
    if (
        firewall.get("lsgkm_prediction_value_read_authorized") is not True
        or firewall.get("row_universe_metadata_read_authorized") is not True
        or firewall.get("raw_score_aggregation_authorized") is not True
        or firewall.get("control_manifest_identity_bind_authorized") is not True
        or any(
            firewall.get(field) is not False
            for field in (
                "control_prediction_member_open_authorized",
                "outcome_access_authorized",
                "reporter_count_access_authorized",
                "benchmark_metric_calculation_authorized",
                "calibration_or_head_fit_authorized",
                "candidate_selection_authorized",
                "sealed_asset_access_authorized",
            )
        )
    ):
        raise LSGKMAggregationError("aggregation action firewall differs")
    output = config.get("output_contract", {})
    if (
        tuple(output.get("readout_ids", ())) != tuple(READOUTS)
        or tuple(output.get("contexts", ())) != CONTEXTS
        or output.get("rows_per_readout") != 10330
        or output.get("total_rows") != 20660
        or output.get("static_score_is_context_specific") is not False
        or output.get("static_native_score_directly_stackable") is not False
        or output.get("calibration_status") != "uncalibrated_static_rank_score"
    ):
        raise LSGKMAggregationError("output contract differs")


def finite_number(text: str, label: str) -> float:
    try:
        value = float(text)
    except ValueError as error:
        raise LSGKMAggregationError(f"non-numeric {label}") from error
    if not math.isfinite(value):
        raise LSGKMAggregationError(f"non-finite {label}")
    return value


def collect_fit_scores(
    authority_path: Path,
    manifest: Mapping[str, Any],
    roster: Mapping[str, Mapping[str, Any]],
    split_id: str,
    model_seeds: Sequence[int],
    *,
    scale_probe: bool,
) -> tuple[dict[tuple[int, str], dict[str, str]], list[dict[str, Any]]]:
    metadata = manifest.get("metadata", {})
    expected_status = (
        "pass_one_fit_production_command_scale_probe"
        if scale_probe
        else "pass_authorized_split_bundle"
    )
    expected_receipt_status = (
        "pass_one_fit_production_command_scale_probe"
        if scale_probe
        else "pass_one_authorized_full_campaign_fit"
    )
    if metadata.get("status") != expected_status:
        raise LSGKMAggregationError("fit authority status differs")
    if scale_probe:
        if metadata.get("split_id") != split_id or metadata.get("model_seed") != model_seeds[0]:
            raise LSGKMAggregationError("scale fit identity differs")
    else:
        if (
            metadata.get("split_id") != split_id
            or tuple(metadata.get("model_seeds", ())) != tuple(model_seeds)
            or metadata.get("fit_count") != len(model_seeds)
        ):
            raise LSGKMAggregationError("bundle fit identity differs")
    result: dict[tuple[int, str], dict[str, str]] = {}
    receipts: list[dict[str, Any]] = []
    for seed in model_seeds:
        prefix = "probe" if scale_probe else f"fits/{split_id}.seed{seed}"
        receipt_member = f"{prefix}/{'scale_probe_receipt.json' if scale_probe else 'fit_receipt.json'}"
        score_member = f"{prefix}/predictions/held_element_scores.tsv.gz"
        receipt_path = bind_openable_member(authority_path, roster, receipt_member)
        score_path = bind_openable_member(authority_path, roster, score_member)
        receipt = load_json(receipt_path)
        receipt_split = (
            receipt.get("active_fit", {}).get("split_id")
            if scale_probe
            else receipt.get("split_id")
        )
        receipt_seed = (
            receipt.get("active_fit", {}).get("model_seed")
            if scale_probe
            else receipt.get("model_seed")
        )
        if (
            receipt.get("status") != expected_receipt_status
            or receipt_split != split_id
            or receipt_seed != seed
            or receipt.get("benchmark_metrics_calculated") is not False
            or any(
                receipt.get(field) is not False
                for field in (
                    "outcomes_read",
                    "reporter_counts_read",
                    "control_prediction_values_read",
                    "sealed_assets_read",
                )
            )
        ):
            raise LSGKMAggregationError("fit receipt or firewall differs")
        rows = read_gzip_tsv(score_path, SCORE_FIELDS)
        if len(rows) != receipt.get("held_elements_scored"):
            raise LSGKMAggregationError("score denominator differs")
        for row in rows:
            key = (seed, row["element_id"])
            if (
                key in result
                or row["split_id"] != split_id
                or int(row["model_seed"]) != seed
                or row["model_state_id"] != "hepatocyte"
            ):
                raise LSGKMAggregationError("score identity differs")
            for score_field, _unit in READOUTS.values():
                finite_number(row[score_field], score_field)
            result[key] = row
        receipts.append(
            {
                "split_id": split_id,
                "model_seed": seed,
                "model_sha256": receipt["model_receipt"]["sha256"],
                "scores_sha256": receipt["held_element_scores_sha256"],
                "held_elements_scored": receipt["held_elements_scored"],
            }
        )
    return result, receipts


def aggregate(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    config = load_json(config_path)
    validate_config(config)
    if output.exists():
        raise LSGKMAggregationError("aggregation output exists")
    output.mkdir(parents=True, mode=0o750)

    incident = config["protocol_incident"]
    incident_path = safe_project_file(root, incident["path"], incident["sha256"])
    incident_record = load_json(incident_path)
    if (
        incident_record.get("status")
        != "contained_narrow_control_value_exposure_no_label_or_metric_leakage"
        or incident_record.get("containment", {}).get("incident_artifact_contains_no_prediction_value")
        is not True
    ):
        raise LSGKMAggregationError("protocol incident containment differs")
    validation_path, validation_manifest, validation_roster = bind_manifest(
        root, incident["validation_authority"]
    )
    validation_receipt = load_json(
        bind_openable_member(
            validation_path, validation_roster, "validation_receipt.json"
        )
    )
    if (
        validation_receipt.get("status") != "pass_structure_only_incident_validation"
        or validation_receipt.get("incident_sha256") != incident["sha256"]
        or validation_receipt.get("referenced_prediction_member_opened") is not False
    ):
        raise LSGKMAggregationError("incident validation differs")

    campaign_path = safe_project_file(
        root, config["campaign_config"]["path"], config["campaign_config"]["sha256"]
    )
    campaign = load_json(campaign_path)
    if campaign.get("design_id") != DESIGN_ID:
        raise LSGKMAggregationError("campaign config differs")
    authorization_path, authorization_manifest, authorization_roster = bind_manifest(
        root, config["authorization"]
    )
    authorization = load_json(
        bind_openable_member(
            authorization_path,
            authorization_roster,
            "authorization/full_campaign_authorization.json",
        )
    )
    if (
        authorization.get("status") != "authorized_remaining_24_outcome_blind_fits"
        or len(authorization.get("authorized_run_fits", ())) != 24
    ):
        raise LSGKMAggregationError("campaign authorization differs")

    row_path, row_manifest, row_roster = bind_manifest(root, config["row_universe"])
    row_member = config["row_universe"]["member"]
    row_file = bind_openable_member(
        row_path, row_roster, row_member, config["row_universe"]["member_sha256"]
    )
    row_rows = read_tsv(row_file, ROW_FIELDS)
    if len(row_rows) != 10330:
        raise LSGKMAggregationError("canonical row denominator differs")

    control_path, control_manifest, control_roster = bind_manifest(
        root, config["shared_control"]
    )
    control_member = bind_manifest_record_only(
        control_roster,
        config["shared_control"]["prediction_member"],
        config["shared_control"]["prediction_member_sha256"],
    )
    if control_path.parent / config["shared_control"]["prediction_member"] == row_file:
        raise LSGKMAggregationError("control member alias differs")

    scale_path, scale_manifest, scale_roster = bind_manifest(root, config["scale_fit"])
    score_lookup, fit_receipts = collect_fit_scores(
        scale_path,
        scale_manifest,
        scale_roster,
        REUSED_FIT[0],
        (REUSED_FIT[1],),
        scale_probe=True,
    )
    observed_fits = {REUSED_FIT}
    for bundle in config["bundle_fits"]:
        bundle_path, bundle_manifest, bundle_roster = bind_manifest(root, bundle)
        bundle_scores, bundle_receipts = collect_fit_scores(
            bundle_path,
            bundle_manifest,
            bundle_roster,
            bundle["split_id"],
            tuple(bundle["model_seeds"]),
            scale_probe=False,
        )
        if set(score_lookup) & set(bundle_scores):
            raise LSGKMAggregationError("duplicate score key across fit authorities")
        score_lookup.update(bundle_scores)
        fit_receipts.extend(bundle_receipts)
        observed_fits.update(
            (bundle["split_id"], seed) for seed in bundle["model_seeds"]
        )
    if observed_fits != expected_remaining_fits() | {REUSED_FIT} or len(fit_receipts) != 25:
        raise LSGKMAggregationError("25-fit rectangle differs")
    if len(score_lookup) != 5165:
        raise LSGKMAggregationError("seed by element score rectangle differs")

    row_identities: set[tuple[int, str]] = set()
    contexts_by_identity: dict[tuple[int, str], set[str]] = {}
    for row in row_rows:
        seed = int(row["seed"])
        key = (seed, row["element_id"])
        score = score_lookup.get(key)
        if (
            score is None
            or row["source_locus_group_id"] != score["source_locus_group_id"]
            or row["long_range_block_id"] != score["long_range_block_id"]
            or row["outer_fold"] != score["outer_fold"]
            or row["assay_context_id"] not in CONTEXTS
        ):
            raise LSGKMAggregationError("row-to-score identity differs")
        row_identities.add((seed, row["row_hash"]))
        contexts_by_identity.setdefault(key, set()).add(row["assay_context_id"])
    if (
        len(row_identities) != 10330
        or len(contexts_by_identity) != 5165
        or any(contexts != set(CONTEXTS) for contexts in contexts_by_identity.values())
    ):
        raise LSGKMAggregationError("canonical context duplication differs")

    output_rows: list[dict[str, object]] = []
    for head_id, (score_field, unit) in READOUTS.items():
        for row in row_rows:
            score = score_lookup[(int(row["seed"]), row["element_id"])]
            output_rows.append(
                {
                    **row,
                    "model_id": "lsgkm_exact_dinucleotide_null",
                    "head_id": head_id,
                    "prediction": score[score_field],
                    "prediction_unit": unit,
                    "calibration_status": "uncalibrated_static_rank_score",
                    "context_specific_prediction": "false",
                    "experimental_replicates": "4",
                    "biological_donors": "0",
                    "outcome_role": "exposed_development_MPRA_only",
                }
            )
    prediction_path = output / "raw_predictions.tsv.gz"
    prediction_count = write_gzip_tsv(prediction_path, OUTPUT_FIELDS, output_rows)
    if prediction_count != 20660:
        raise LSGKMAggregationError("raw prediction denominator differs")

    receipt = {
        "schema_version": "masld-bench-lsgkm-gse281364-dinucleotide-aggregation-v1",
        "status": "pass_outcome_blind_raw_score_aggregation",
        "design_id": DESIGN_ID,
        "config_sha256": file_sha256(config_path),
        "fit_count": 25,
        "fit_receipts": fit_receipts,
        "seed_by_element_scores": len(score_lookup),
        "row_universe_rows": len(row_rows),
        "prediction_rows_per_readout": 10330,
        "total_prediction_rows": prediction_count,
        "raw_predictions_sha256": file_sha256(prediction_path),
        "readouts": list(READOUTS),
        "readouts_are_independent_model_families": False,
        "static_score_is_context_specific": False,
        "static_native_score_directly_stackable": False,
        "calibration_status": "uncalibrated_static_rank_score",
        "shared_control_artifacts_sha256": config["shared_control"]["artifacts_sha256"],
        "shared_control_prediction_member": config["shared_control"]["prediction_member"],
        "shared_control_prediction_member_sha256": control_member["sha256"],
        "shared_control_prediction_member_opened": False,
        "protocol_incident_id": incident_record["incident_id"],
        "benchmark_metrics_calculated": False,
        "calibration_or_head_fit_executed": False,
        "candidate_selection_executed": False,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "control_prediction_values_read": False,
        "sealed_assets_read": False,
    }
    (output / "aggregation_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    aggregate(
        arguments.project_root.resolve(strict=True),
        arguments.config.resolve(strict=True),
        arguments.output,
    )


if __name__ == "__main__":
    main()
