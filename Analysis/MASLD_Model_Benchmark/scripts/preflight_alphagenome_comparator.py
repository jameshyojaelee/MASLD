#!/usr/bin/env python3
"""Freeze AlphaGenome's current restricted, non-executable comparator disposition."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import tomllib
from typing import Any, Iterable, Mapping

from masld_bench.artifacts import (
    ArtifactError,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)


CONFIG_SCHEMA = "masld-bench-alphagenome-comparator-preflight-v1"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
RELEVANT_CURIES = (
    "CL:0000115",
    "CL:0000182",
    "CL:0000235",
    "CL:0000632",
    "UBERON:0001114",
    "UBERON:0001115",
    "UBERON:0002107",
)
ROSTER_FIELDS = (
    "output_type",
    "ontology_curie",
    "track_name",
    "strand",
    "biosample_name",
    "biosample_type",
    "biosample_stage",
    "assay",
    "data_source",
    "gtex_tissue",
    "histone_mark",
    "transcription_factor",
)


class AlphaGenomePreflightError(RuntimeError):
    """Raised when a frozen AlphaGenome authority differs."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise AlphaGenomePreflightError(f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AlphaGenomePreflightError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise AlphaGenomePreflightError(f"{label} must be a JSON object")
    return value


def project_path(root: Path, relative: str, *, must_exist: bool) -> Path:
    value = Path(relative)
    if value.is_absolute() or not value.parts or ".." in value.parts:
        raise AlphaGenomePreflightError(f"unsafe project path: {relative}")
    path = reject_symlink_components(root / value, label="AlphaGenome authority")
    if must_exist:
        try:
            path.resolve(strict=True).relative_to(root)
        except (OSError, ValueError) as error:
            raise AlphaGenomePreflightError(
                f"authority is missing or escapes project root: {relative}"
            ) from error
    return path


def validate_authorities(
    root: Path, config: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    expected_names = {
        "checkpoint_contract",
        "exposure_audit",
        "development_crosswalk",
        "registry",
        "admission_contract",
        "admission_sources",
        "sequence_fixture",
        "project_reference",
    }
    raw = config.get("frozen_authorities")
    if not isinstance(raw, dict) or set(raw) != expected_names:
        raise AlphaGenomePreflightError("frozen authority census differs")
    records: dict[str, dict[str, Any]] = {}
    for name, binding in raw.items():
        if not isinstance(binding, dict):
            raise AlphaGenomePreflightError(f"authority binding differs: {name}")
        expected_hash = str(binding.get("sha256", ""))
        if not SHA256.fullmatch(expected_hash):
            raise AlphaGenomePreflightError(f"authority SHA-256 differs: {name}")
        path = project_path(root, str(binding.get("path", "")), must_exist=True)
        kind = binding.get("kind")
        if kind == "file":
            if path.is_symlink() or not path.is_file():
                raise AlphaGenomePreflightError(f"authority is not a file: {name}")
            observed = sha256_file(path)
        elif kind == "tree":
            try:
                verify_frozen_tree(path)
            except ArtifactError as error:
                raise AlphaGenomePreflightError(
                    f"frozen authority tree differs: {name}: {error}"
                ) from error
            observed = sha256_file(path / "ARTIFACTS.json")
        else:
            raise AlphaGenomePreflightError(f"authority kind differs: {name}")
        if observed != expected_hash:
            raise AlphaGenomePreflightError(f"frozen authority changed: {name}")
        records[name] = {
            "kind": kind,
            "path": str(path.relative_to(root)),
            "sha256": observed,
        }
    return records


def quoted_field(block: str, field: str, *, indent: int = 6) -> str:
    prefix = re.escape(" " * indent)
    match = re.search(
        rf"^{prefix}{re.escape(field)}: \"([^\"]*)\"$", block, re.M
    )
    return match.group(1) if match else ""


def scalar_field(block: str, field: str, *, indent: int = 6) -> str:
    prefix = re.escape(" " * indent)
    match = re.search(rf"^{prefix}{re.escape(field)}: ([^\s]+)$", block, re.M)
    return match.group(1) if match else ""


def parse_metadata_block(output_type: str, block: str) -> dict[str, str]:
    ontology_type_match = re.search(
        r"^        ontology_type: ONTOLOGY_TYPE_([A-Z_]+)$", block, re.M
    )
    ontology_id_match = re.search(r"^        id: (\d+)$", block, re.M)
    ontology_curie = ""
    if ontology_type_match and ontology_id_match:
        ontology_curie = (
            f"{ontology_type_match.group(1)}:{int(ontology_id_match.group(1)):07d}"
        )
    biosample_match = re.search(r"^      biosample \{\n(.*?)^      \}$", block, re.M | re.S)
    biosample = biosample_match.group(1) if biosample_match else ""
    return {
        "output_type": output_type.removeprefix("OUTPUT_TYPE_"),
        "ontology_curie": ontology_curie,
        "track_name": quoted_field(block, "name"),
        "strand": scalar_field(block, "strand"),
        "biosample_name": quoted_field(biosample, "name", indent=8),
        "biosample_type": scalar_field(biosample, "type", indent=8).removeprefix(
            "BIOSAMPLE_TYPE_"
        ),
        "biosample_stage": quoted_field(biosample, "stage", indent=8),
        "assay": quoted_field(block, "assay"),
        "data_source": quoted_field(block, "data_source"),
        "gtex_tissue": quoted_field(block, "gtex_tissue"),
        "histone_mark": quoted_field(block, "histone_mark_code"),
        "transcription_factor": quoted_field(block, "transcription_factor_code"),
    }


def iter_metadata_blocks(lines: Iterable[str]) -> Iterable[tuple[str, str]]:
    output_type = ""
    collecting = False
    depth = 0
    block: list[str] = []
    for raw_line in lines:
        line = raw_line.rstrip("\n")
        if line.startswith("  output_type: OUTPUT_TYPE_"):
            output_type = line.split(":", 1)[1].strip()
        if not collecting and line == "    metadata {":
            if not output_type:
                raise AlphaGenomePreflightError("metadata precedes output type")
            collecting = True
            depth = 1
            block = [line]
            continue
        if collecting:
            block.append(line)
            depth += line.count("{") - line.count("}")
            if depth == 0:
                yield output_type, "\n".join(block)
                collecting = False
                block = []
    if collecting:
        raise AlphaGenomePreflightError("unterminated output metadata block")


def audit_metadata(path: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    contract = config["metadata_contract"]
    if sha256_file(path) != contract["source_sha256"]:
        raise AlphaGenomePreflightError("human output metadata changed")
    rows = [
        parse_metadata_block(output_type, block)
        for output_type, block in iter_metadata_blocks(
            path.open(encoding="utf-8")
        )
    ]
    padding = sum(row["track_name"] == "Padding" for row in rows)
    nonpadding_rows = [row for row in rows if row["track_name"] != "Padding"]
    splice_junctions = sum(
        row["output_type"] == "SPLICE_JUNCTIONS" for row in nonpadding_rows
    )
    relevant_rows = [row for row in nonpadding_rows if row["ontology_curie"] in RELEVANT_CURIES]
    counts = {
        curie: sum(row["ontology_curie"] == curie for row in relevant_rows)
        for curie in RELEVANT_CURIES
    }
    cholangiocyte = sum(
        "cholangiocyte" in (row["track_name"] + " " + row["biosample_name"]).casefold()
        for row in nonpadding_rows
    )
    kupffer = sum(
        "kupffer" in (row["track_name"] + " " + row["biosample_name"]).casefold()
        for row in nonpadding_rows
    )
    if (
        len(rows) != contract["encoded_output_slots"]
        or padding != contract["padding_slots"]
        or len(nonpadding_rows) != contract["nonpadding_metadata_rows"]
        or len(nonpadding_rows) + splice_junctions
        != contract["active_tracks_with_second_splice_junction_strands"]
        or counts != contract["relevant_ontology_expected_counts"]
        or cholangiocyte != contract["cholangiocyte_exact_name_tracks"]
        or kupffer != contract["kupffer_exact_name_tracks"]
    ):
        raise AlphaGenomePreflightError("AlphaGenome output ontology census differs")
    by_output: dict[str, dict[str, int]] = {}
    for curie in RELEVANT_CURIES:
        by_output[curie] = {}
        for row in relevant_rows:
            if row["ontology_curie"] == curie:
                output = row["output_type"]
                by_output[curie][output] = by_output[curie].get(output, 0) + 1
    expected_cell_outputs = contract["cell_ontology_expected_counts_by_output"]
    observed_cell_outputs = {
        curie: by_output[curie] for curie in expected_cell_outputs
    }
    if observed_cell_outputs != expected_cell_outputs:
        raise AlphaGenomePreflightError("AlphaGenome cell-track assay census differs")
    return {
        "encoded_output_slots": len(rows),
        "padding_slots": padding,
        "nonpadding_metadata_rows": len(nonpadding_rows),
        "splice_junction_metadata_rows": splice_junctions,
        "active_tracks_with_second_splice_junction_strands": len(nonpadding_rows)
        + splice_junctions,
        "relevant_ontology_counts": counts,
        "relevant_ontology_counts_by_output": by_output,
        "cholangiocyte_exact_name_tracks": cholangiocyte,
        "kupffer_exact_name_tracks": kupffer,
        "relevant_rows": relevant_rows,
    }


def registry_row(path: Path) -> Mapping[str, Any]:
    raw = tomllib.loads(path.read_text(encoding="utf-8"))
    matches = [row for row in raw.get("models", []) if row.get("model_id") == "alphagenome"]
    if len(matches) != 1:
        raise AlphaGenomePreflightError("AlphaGenome registry row differs")
    return matches[0]


def validate_scientific_contracts(root: Path, config: Mapping[str, Any]) -> None:
    authorities = config["frozen_authorities"]
    checkpoint = load_json(
        root / authorities["checkpoint_contract"]["path"], label="checkpoint contract"
    )
    exposure = load_json(
        root / authorities["exposure_audit"]["path"], label="exposure audit"
    )
    crosswalk = load_json(
        root / authorities["development_crosswalk"]["path"],
        label="development crosswalk",
    )
    source_receipt = load_json(
        root
        / authorities["admission_sources"]["path"]
        / "review/source_admission_receipt.json",
        label="source admission receipt",
    )
    fixture_summary = load_json(
        root / authorities["sequence_fixture"]["path"] / "fixture/summary.json",
        label="sequence fixture summary",
    )
    firewall = load_json(
        root / authorities["sequence_fixture"]["path"] / "fixture/outcome_firewall.json",
        label="sequence fixture firewall",
    )
    registry = registry_row(root / authorities["registry"]["path"])
    if (
        checkpoint.get("weight_license") != "NONCOMMERCIAL_GATED"
        or checkpoint.get("terms_audit", {}).get("open_champion_disposition")
        != "ineligible_under_current_model_terms"
        or checkpoint.get("native_checkpoint", {}).get("status")
        != "exact_repository_revision_registered_but_gated_bytes_not_acquired"
        or checkpoint.get("native_checkpoint", {}).get("huggingface_revision")
        != config["local_checkpoint_route"]["checkpoint_revision"]
    ):
        raise AlphaGenomePreflightError("checkpoint or terms contract differs")
    if (
        exposure.get("checkpoint_findings", {})
        .get("alphagenome_all_folds", {})
        .get("exposure_state")
        != "target_label_unexposed"
        or exposure.get("sealed_output_gate", {}).get("primary_roster_complete")
        is not False
        or crosswalk.get("findings", {}).get("gse281364", {}).get("exposure_state")
        != "clean_declared"
        or crosswalk.get("findings", {}).get("gse289173", {}).get("exposure_state")
        != "target_label_unexposed"
    ):
        raise AlphaGenomePreflightError("exposure or lineage-roster contract differs")
    alpha_disposition = source_receipt.get("model_dispositions", {}).get(
        "alphagenome_all_folds", {}
    )
    if (
        source_receipt.get("status") != "pass"
        or source_receipt.get("checkpoint_bytes_downloaded") is not False
        or source_receipt.get("model_forward_executed") is not False
        or source_receipt.get("observed_outcomes_loaded") is not False
        or alpha_disposition.get("next_model_job_allowed") is not False
        or alpha_disposition.get("open_champion_eligible") is not False
    ):
        raise AlphaGenomePreflightError("source admission disposition differs")
    if (
        fixture_summary.get("status") != "pass"
        or fixture_summary.get("anchors") != 3
        or fixture_summary.get("reference_build") != "GRCh38.p14"
        or fixture_summary.get("checkpoint_bytes_loaded") is not False
        or fixture_summary.get("runtime_forward_executed") is not False
        or fixture_summary.get("observed_outcomes_loaded") is not False
        or firewall.get("variant_or_regulatory_outcomes_loaded") is not False
        or firewall.get("sealed_labels_loaded") is not False
    ):
        raise AlphaGenomePreflightError("outcome-blind fixture differs")
    if (
        registry.get("status") != "restricted_comparator"
        or registry.get("admission_blocking") is not True
        or "noncommercial" not in str(registry.get("license_status", "")).casefold()
        or registry.get("checkpoint_sha256") != "UNRESOLVED"
    ):
        raise AlphaGenomePreflightError("restricted registry disposition differs")


def audit_registered_artifacts(root: Path, config: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    results: dict[str, Any] = {}
    blockers: list[str] = []
    for route_name in ("local_checkpoint_route", "api_route"):
        route = config[route_name]
        route_results: dict[str, Any] = {}
        for artifact_name, binding in route["registered_artifacts"].items():
            path = project_path(root, binding["path"], must_exist=False)
            expected = binding["artifacts_sha256"]
            record = {
                "path": str(path.relative_to(root)),
                "expected_artifacts_sha256": expected,
                "exists": path.is_dir(),
                "passed": False,
            }
            if expected == "UNRESOLVED":
                blockers.append(f"artifact_hash_unresolved:{route_name}:{artifact_name}")
                if not path.is_dir():
                    blockers.append(f"artifact_missing:{route_name}:{artifact_name}")
                record["status"] = "blocked_unbound"
            else:
                if not SHA256.fullmatch(expected):
                    raise AlphaGenomePreflightError(
                        f"registered artifact SHA-256 differs: {artifact_name}"
                    )
                blockers.append(
                    f"artifact_activation_not_implemented:{route_name}:{artifact_name}"
                )
                record["status"] = "blocked_requires_new_preflight_revision"
            route_results[artifact_name] = record
        results[route_name] = route_results
    blockers.extend(
        [
            "restricted_noncommercial_never_open_champion",
            "four_lineage_signed_RNA_roster_incomplete",
            "api_protocol_does_not_return_immutable_server_build_identity",
            "api_metadata_request_is_not_model_version_specific",
            "current_execution_policy_closed",
        ]
    )
    return results, blockers


def preflight(*, root: Path, config_path: Path) -> dict[str, Any]:
    root = reject_symlink_components(root, label="benchmark root").resolve(strict=True)
    config_path = reject_symlink_components(
        config_path, label="AlphaGenome comparator config"
    ).resolve(strict=True)
    try:
        config_path.relative_to(root)
    except ValueError as error:
        raise AlphaGenomePreflightError("config escapes benchmark root") from error
    config = load_json(config_path, label="AlphaGenome comparator config")
    if (
        config.get("schema_version") != CONFIG_SCHEMA
        or config.get("model_id") != "alphagenome"
        or config.get("checkpoint_id") != "alphagenome_all_folds"
    ):
        raise AlphaGenomePreflightError("AlphaGenome comparator config identity differs")
    if (
        config.get("terms_contract", {}).get("open_champion_eligible") is not False
        or config.get("terms_contract", {}).get("conditional_open_model_training_allowed")
        is not False
        or config.get("local_checkpoint_route", {}).get("current_execution_allowed")
        is not False
        or config.get("api_route", {}).get("current_execution_allowed") is not False
        or config.get("api_route", {}).get("api_key_inspection_allowed") is not False
    ):
        raise AlphaGenomePreflightError("restricted execution firewall differs")
    authorities = validate_authorities(root, config)
    validate_scientific_contracts(root, config)
    metadata_path = (
        root
        / config["frozen_authorities"]["admission_sources"]["path"]
        / "sources/alphagenome_human_metadata"
    )
    metadata = audit_metadata(metadata_path, config)
    artifact_gates, blockers = audit_registered_artifacts(root, config)
    sdk_records = config["api_route"]["pinned_sdk_source_records"]
    if (
        len(sdk_records) != 9
        or any(not SHA256.fullmatch(str(record.get("sha256", ""))) for record in sdk_records)
        or config["api_route"]["logical_model_version"] != "ALL_FOLDS"
        or config["api_route"]["default_model_version_allowed"] is not False
    ):
        raise AlphaGenomePreflightError("API SDK identity or version policy differs")
    rows = metadata.pop("relevant_rows")
    return {
        "schema_version": "masld-bench-alphagenome-comparator-preflight-receipt-v1",
        "status": config["status"],
        "model_id": "alphagenome",
        "checkpoint_id": "alphagenome_all_folds",
        "config_sha256": sha256_file(config_path),
        "executable": False,
        "terminal_current_release": True,
        "fixed_authorities": authorities,
        "metadata_audit": metadata,
        "relevant_track_rows": rows,
        "artifact_gates": artifact_gates,
        "blockers": blockers,
        "local_checkpoint_route": {
            "checkpoint_downloaded": False,
            "checkpoint_deserialized": False,
            "model_forward_executed": False,
        },
        "api_route": {
            "sdk_source_locally_frozen": False,
            "api_key_inspected": False,
            "api_connection_attempted": False,
            "metadata_endpoint_called": False,
            "prediction_endpoint_called": False,
            "logical_model_version": "ALL_FOLDS",
            "immutable_server_build_identity_available": False,
        },
        "observed_RNA_read": False,
        "observed_ATAC_read": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "metrics_calculated": False,
        "open_champion_eligible": False,
        "conditional_open_model_training_allowed": False,
        "disposition": config["terminal_disposition"],
    }


def write_outputs(output: Path, receipt: Mapping[str, Any]) -> None:
    if output.exists() or output.is_symlink():
        raise AlphaGenomePreflightError(f"refusing to overwrite {output}")
    output.mkdir(parents=True, mode=0o750)
    rows = receipt["relevant_track_rows"]
    with (output / "relevant_track_roster.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=ROSTER_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    frozen = dict(receipt)
    del frozen["relevant_track_rows"]
    frozen["relevant_track_roster_rows"] = len(rows)
    frozen["relevant_track_roster_sha256"] = sha256_file(
        output / "relevant_track_roster.tsv"
    )
    write_json_exclusive(output / "receipt.json", frozen, mode=0o640)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = preflight(root=arguments.root, config_path=arguments.config)
    write_outputs(arguments.output, receipt)
    printable = dict(receipt)
    printable["relevant_track_rows"] = len(receipt["relevant_track_rows"])
    print(json.dumps(printable, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
