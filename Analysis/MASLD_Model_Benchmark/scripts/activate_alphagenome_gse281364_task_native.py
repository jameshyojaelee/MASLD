#!/usr/bin/env python3
"""Freeze an outcome-blind AlphaGenome GSE281364 native-comparator plan."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping

from masld_bench.artifacts import (
    ArtifactError,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)


SCHEMA = "masld-bench-alphagenome-gse281364-task-native-activation-v1"
RECEIPT_SCHEMA = "masld-bench-alphagenome-gse281364-task-native-activation-receipt-v1"
AUTOSOMES = tuple(f"chr{index}" for index in range(1, 23))


class AlphaGenomeActivationError(RuntimeError):
    """Raised when a frozen authority or outcome separation differs."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AlphaGenomeActivationError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    require(path.is_file() and not path.is_symlink(), f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AlphaGenomeActivationError(f"invalid {label}: {error}") from error
    require(isinstance(value, dict), f"{label} must be a JSON object")
    return value


def validate_file(root: Path, record: Mapping[str, Any], *, label: str) -> Path:
    path = reject_symlink_components(root / record["path"], label=label).resolve(strict=True)
    require(root in path.parents, f"{label} escapes benchmark root")
    require(path.is_file() and not path.is_symlink(), f"{label} is not a regular file")
    require(sha256_file(path) == record["sha256"], f"{label} hash differs")
    return path


def validate_tree(root: Path, record: Mapping[str, Any], *, label: str) -> Path:
    path = reject_symlink_components(root / record["path"], label=label).resolve(strict=True)
    require(root in path.parents, f"{label} escapes benchmark root")
    try:
        verify_frozen_tree(path)
    except ArtifactError as error:
        raise AlphaGenomeActivationError(f"{label} differs: {error}") from error
    require(
        sha256_file(path / "ARTIFACTS.json") == record["artifacts_sha256"],
        f"{label} artifact identity differs",
    )
    return path


def parse_assembly_report(path: Path) -> dict[str, dict[str, str]]:
    """Read assembled chromosome identities from an NCBI assembly report."""
    records: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            require(len(fields) == 10, f"invalid assembly-report row in {path.name}")
            (
                sequence_name,
                sequence_role,
                assigned_molecule,
                assigned_type,
                genbank_accession,
                relationship,
                refseq_accession,
                assembly_unit,
                sequence_length,
                ucsc_name,
            ) = fields
            if sequence_role != "assembled-molecule" or assigned_type != "Chromosome":
                continue
            records[ucsc_name] = {
                "sequence_name": sequence_name,
                "assigned_molecule": assigned_molecule,
                "genbank_accession": genbank_accession,
                "relationship": relationship,
                "refseq_accession": refseq_accession,
                "assembly_unit": assembly_unit,
                "sequence_length": sequence_length,
            }
    require(all(contig in records for contig in AUTOSOMES), f"autosome roster incomplete in {path.name}")
    return records


def compare_primary_assemblies(
    p13_path: Path, p14_path: Path, output_path: Path
) -> dict[str, Any]:
    p13 = parse_assembly_report(p13_path)
    p14 = parse_assembly_report(p14_path)
    rows: list[dict[str, str]] = []
    for contig in AUTOSOMES:
        left = p13[contig]
        right = p14[contig]
        require(
            left["genbank_accession"] == right["genbank_accession"]
            and left["refseq_accession"] == right["refseq_accession"]
            and left["sequence_length"] == right["sequence_length"],
            f"GRCh38.p13/p14 primary sequence identity differs for {contig}",
        )
        rows.append({
            "contig": contig,
            "genbank_accession": left["genbank_accession"],
            "refseq_accession": left["refseq_accession"],
            "sequence_length": left["sequence_length"],
            "p13_p14_primary_sequence_accession_equal": "true",
        })
    write_tsv(output_path, rows, tuple(rows[0]))
    return {
        "compared_autosomes": len(rows),
        "all_used_primary_sequence_accessions_equal": True,
        "p13_report_sha256": sha256_file(p13_path),
        "p14_report_sha256": sha256_file(p14_path),
    }


def write_tsv(path: Path, rows: Iterable[Mapping[str, Any]], fields: tuple[str, ...]) -> None:
    require(not path.exists() and not path.is_symlink(), f"refusing to overwrite {path}")
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o750)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    path.chmod(0o640)


def parse_output_metadata(
    metadata_path: Path, sdk_root: Path, expected: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Parse the frozen official all-fold output metadata with the pinned SDK."""
    sys.path.insert(0, str(sdk_root))
    for name in tuple(sys.modules):
        if name == "alphagenome" or name.startswith("alphagenome."):
            del sys.modules[name]
    text_format = importlib.import_module("google.protobuf.text_format")
    model_pb2 = importlib.import_module("alphagenome.protos.dna_model_pb2")
    service_pb2 = importlib.import_module("alphagenome.protos.dna_model_service_pb2")
    message = service_pb2.MetadataResponse()
    text_format.Parse(metadata_path.read_text(encoding="utf-8"), message)

    target_curie = expected["ontology_curie"]
    target_biosample = expected["biosample_name"]
    rows: list[dict[str, Any]] = []
    counts: dict[str, int] = {}
    for output in message.output_metadata:
        output_type = model_pb2.OutputType.Name(output.output_type).removeprefix("OUTPUT_TYPE_")
        payload = output.WhichOneof("payload")
        track_rows = output.tracks.metadata if payload == "tracks" else output.junctions.metadata
        for track_index, track in enumerate(track_rows):
            ontology_term = track.ontology_term
            ontology_type = model_pb2.OntologyType.Name(ontology_term.ontology_type).removeprefix(
                "ONTOLOGY_TYPE_"
            )
            curie = f"{ontology_type}:{ontology_term.id:07d}"
            biosample = track.biosample.name
            if curie != target_curie or biosample != target_biosample:
                continue
            counts[output_type] = counts.get(output_type, 0) + 1
            strand = "not_applicable"
            histone_mark = "not_applicable"
            tf = "not_applicable"
            if payload == "tracks":
                strand = model_pb2.Strand.Name(track.strand).removeprefix("STRAND_")
                histone_mark = track.histone_mark_code or "not_applicable"
                tf = track.transcription_factor_code or "not_applicable"
            rows.append({
                "output_type": output_type,
                "output_track_index0": track_index,
                "track_name": track.name,
                "strand": strand,
                "ontology_curie": curie,
                "biosample_name": biosample,
                "assay": track.assay or "not_available",
                "data_source": track.data_source or "not_available",
                "histone_mark": histone_mark,
                "transcription_factor": tf,
                "primary_zero_shot": str(output_type == expected["primary_zero_shot_output"]).lower(),
            })
    require(counts == expected["counts_by_output_type"], "HepG2 fixed-track census differs")
    return rows, counts


def read_fai(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            require(len(fields) >= 2, "invalid FASTA index row")
            lengths[fields[0]] = int(fields[1])
    return lengths


def build_request_manifest(
    fixture_manifest: Path,
    fai_path: Path,
    task: Mapping[str, Any],
    output_path: Path,
) -> dict[str, Any]:
    lengths = read_fai(fai_path)
    input_length = int(task["input_length_bp"])
    variant_index0 = int(task["variant_index0"])
    require(input_length in {16384, 131072, 524288, 1048576}, "unsupported AlphaGenome length")
    require(variant_index0 * 2 == input_length, "variant is not centered")
    rows: list[dict[str, Any]] = []
    seen_elements: set[str] = set()
    seen_groups: set[str] = set()
    with fixture_manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "fixture_id", "element_id", "outer_locus_sequence_group_id", "outer_fold",
            "contig", "variant_pos0", "variant_pos1", "variant_index0", "ref", "alt",
        }
        require(reader.fieldnames is not None and required.issubset(reader.fieldnames), "fixture columns differ")
        for source in reader:
            contig = source["contig"]
            pos0 = int(source["variant_pos0"])
            pos1 = int(source["variant_pos1"])
            require(contig in AUTOSOMES and contig in lengths, "non-autosomal fixture coordinate")
            require(pos1 == pos0 + 1, "fixture coordinate convention differs")
            require(source["ref"] in "ACGT" and source["alt"] in "ACGT", "fixture is not a SNV")
            require(source["ref"] != source["alt"], "REF and ALT are equal")
            start0 = pos0 - variant_index0
            end0 = start0 + input_length
            require(start0 >= 0 and end0 <= lengths[contig], "centered 524-kb request crosses a contig boundary")
            element = source["element_id"]
            group = source["outer_locus_sequence_group_id"]
            require(element not in seen_elements and group not in seen_groups, "fixture is not one per group")
            seen_elements.add(element)
            seen_groups.add(group)
            rows.append({
                "fixture_id": source["fixture_id"],
                "element_id": element,
                "outer_locus_sequence_group_id": group,
                "outer_fold": source["outer_fold"],
                "contig": contig,
                "interval_start0": start0,
                "interval_end0": end0,
                "variant_position1": pos1,
                "variant_index0": pos0 - start0,
                "ref": source["ref"],
                "alt": source["alt"],
                "logical_model_version": "ALL_FOLDS",
                "allele_effect_sign": task["allele_effect_sign"],
                "prediction_context_policy": "one_static_prediction_reused_for_HepG2_control_and_HepG2_PAOA",
            })
    require(len(rows) == 1033 and len(seen_groups) == 1033, "fixture row universe differs")
    require(set(row["outer_fold"] for row in rows) == {"0", "1", "2", "3", "4"}, "fold roster differs")
    write_tsv(output_path, rows, tuple(rows[0]))
    return {
        "requests": len(rows),
        "outer_locus_sequence_groups": len(seen_groups),
        "input_length_bp": input_length,
        "variant_index0": variant_index0,
        "contigs": sorted({row["contig"] for row in rows}),
        "outcomes_read": False,
    }


def validate_snapshot_sources(source_dir: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    sdk = config["frozen_sdk"]
    checkpoint = config["frozen_checkpoint"]
    tag = load_json(source_dir / "github_tag_v0.8.0.json", label="GitHub tag snapshot")
    pypi = load_json(source_dir / "pypi_0.8.0.json", label="PyPI snapshot")
    hf = load_json(source_dir / "hf_checkpoint_revision.json", label="Hugging Face snapshot")
    require((tag.get("object") or {}).get("sha") == sdk["revision"], "GitHub tag revision differs")
    require((pypi.get("info") or {}).get("version") == sdk["version"], "PyPI release differs")
    require(hf.get("sha") == checkpoint["revision"], "Hugging Face checkpoint revision differs")
    require(hf.get("gated") == "auto" and not hf.get("private"), "Hugging Face gate differs")

    records: dict[str, Any] = {}
    for name in (
        "github_tag_v0.8.0.json",
        "pypi_0.8.0.json",
        "hf_checkpoint_revision.json",
        "api_terms.html",
        "api_terms.headers.txt",
        "model_terms.html",
        "model_terms.headers.txt",
        "GRCh38.p13_assembly_report.txt",
        "GRCh38.p14_assembly_report.txt",
    ):
        path = source_dir / name
        require(path.is_file() and not path.is_symlink() and path.stat().st_size > 0, f"source absent: {name}")
        records[name] = {"sha256": sha256_file(path), "size_bytes": path.stat().st_size}
    api_terms_text = (source_dir / "api_terms.html").read_text(encoding="utf-8", errors="replace").lower()
    model_terms_text = (source_dir / "model_terms.html").read_text(encoding="utf-8", errors="replace").lower()
    return {
        "records": records,
        "api_terms_noncommercial_marker_visible_in_raw_html": "non-commercial" in api_terms_text,
        "model_terms_noncommercial_marker_visible_in_raw_html": "non-commercial" in model_terms_text,
        "terms_are_mutable_timestamped_web_snapshots": True,
        "terms_semantics_require_authoritative_human_review": True,
    }


def credential_presence(config: Mapping[str, Any]) -> dict[str, bool]:
    policy = config["credential_policy"]
    require(policy["presence_boolean_only"], "credential policy differs")
    require(policy["values_hashes_lengths_and_prefixes_forbidden"], "secret firewall differs")
    api_present = bool(os.environ.get(policy["api_environment_variable"]))
    hf_present = any(bool(os.environ.get(name)) for name in policy["huggingface_environment_variables"])
    return {"api_credential_present": api_present, "huggingface_credential_present": hf_present}


def activate(root: Path, config_path: Path, source_dir: Path, output: Path) -> dict[str, Any]:
    root = reject_symlink_components(root, label="benchmark root").resolve(strict=True)
    config_path = reject_symlink_components(config_path, label="activation config").resolve(strict=True)
    source_dir = reject_symlink_components(source_dir, label="source snapshot").resolve(strict=True)
    require(root in config_path.parents and root in source_dir.parents, "activation input escapes benchmark root")
    require(not output.exists() and not output.is_symlink(), "refusing to overwrite activation output")
    config = load_json(config_path, label="activation config")
    require(config.get("schema_version") == SCHEMA, "activation schema differs")
    disposition = config["activation_disposition"]
    require(
        not disposition["api_metadata_call_allowed"]
        and not disposition["api_prediction_call_allowed"]
        and not disposition["checkpoint_download_allowed"]
        and not disposition["outcome_access_allowed"]
        and not disposition["sealed_access_allowed"]
        and not disposition["open_champion_eligible"],
        "closed activation disposition differs",
    )

    authorities = config["authorities"]
    for key in ("checkpoint_contract", "development_crosswalk", "exposure_audit"):
        validate_file(root, authorities[key], label=key)
    sdk_artifact = validate_tree(
        root,
        {"path": config["frozen_sdk"]["prior_artifact_path"], "artifacts_sha256": config["frozen_sdk"]["prior_artifacts_sha256"]},
        label="pinned SDK artifact",
    )
    official_source = validate_tree(root, authorities["official_source_artifact"], label="official source artifact")
    fixture = validate_tree(root, authorities["outcome_blind_fixture"], label="outcome-blind fixture")
    split = validate_tree(root, authorities["outcome_blind_split"], label="outcome-blind split")
    reference = validate_tree(root, authorities["project_reference"], label="project reference")
    metadata_path = validate_file(root, authorities["official_output_metadata"], label="official output metadata")
    fixture_receipt = load_json(fixture / "fixture/receipt.json", label="fixture receipt")
    split_receipt = load_json(split / "split/receipt.json", label="split receipt")
    require(
        fixture_receipt["outcomes_read"] is False
        and fixture_receipt["reporter_counts_read"] is False
        and fixture_receipt["sealed_outcomes_read"] is False
        and split_receipt["outcomes_read"] is False
        and split_receipt["reporter_counts_read"] is False
        and split_receipt["sealed_outcomes_read"] is False,
        "outcome-blind authority opened outcomes",
    )

    output.mkdir(parents=True, mode=0o750)
    snapshot = validate_snapshot_sources(source_dir, config)
    assembly = compare_primary_assemblies(
        source_dir / "GRCh38.p13_assembly_report.txt",
        source_dir / "GRCh38.p14_assembly_report.txt",
        output / "reference/primary_assembly_crosswalk.tsv",
    )
    track_rows, track_counts = parse_output_metadata(
        metadata_path, sdk_artifact / "audit/sdk", config["expected_hepg2_tracks"]
    )
    write_tsv(
        output / "task/hepg2_fixed_tracks.tsv",
        track_rows,
        (
            "output_type", "output_track_index0", "track_name", "strand", "ontology_curie",
            "biosample_name", "assay", "data_source", "histone_mark", "transcription_factor",
            "primary_zero_shot",
        ),
    )
    scorer_rows = config["task_contract"]["native_scorers"]
    write_tsv(
        output / "task/native_scorers.tsv",
        scorer_rows,
        ("scorer_id", "output_type", "width_bp", "aggregation", "role"),
    )
    request_manifest = build_request_manifest(
        fixture / "fixture/manifest.tsv",
        reference / "GRCh38.p14.sequence_model.fa.fai",
        config["task_contract"],
        output / "task/request_manifest.tsv",
    )
    credentials = credential_presence(config)
    if credentials["api_credential_present"]:
        credential_status = "api_parity_campaign_ready_but_not_run"
    elif credentials["huggingface_credential_present"]:
        credential_status = "local_access_audit_ready_but_terms_acceptance_unresolved"
    else:
        credential_status = "execution_blocked_no_available_credential"
    status = "pass_outcome_blind_task_native_activation_" + credential_status
    receipt = {
        "schema_version": RECEIPT_SCHEMA,
        "status": status,
        "config_sha256": sha256_file(config_path),
        "sdk_version": config["frozen_sdk"]["version"],
        "sdk_revision": config["frozen_sdk"]["revision"],
        "logical_model_version": config["frozen_sdk"]["logical_model_version"],
        "checkpoint_revision": config["frozen_checkpoint"]["revision"],
        "source_snapshot": snapshot,
        "reference_crosswalk": assembly,
        "request_manifest": request_manifest,
        "hepg2_fixed_track_counts": track_counts,
        "native_scorers": scorer_rows,
        "prespecified_comparators": config["task_contract"]["prespecified_comparators"],
        "credential_presence": credentials,
        "credential_status": credential_status,
        "task_role": config["task_contract"]["role"],
        "contexts": config["task_contract"]["contexts"],
        "sequence_prediction_reused_identically_across_contexts": True,
        "api_metadata_called": False,
        "api_prediction_called": False,
        "checkpoint_bytes_downloaded": False,
        "checkpoint_deserialized": False,
        "model_forward_executed": False,
        "observed_outcomes_loaded": False,
        "sealed_assets_loaded": False,
        "metrics_calculated": False,
        "open_champion_eligible": False,
        "remaining_blockers": [
            "noncommercial_model_and_API_terms_restrict_release_role",
            "no_available_AlphaGenome_API_or_Hugging_Face_credential_in_job_environment" if not any(credentials.values()) else "credential_specific_authority_and_parity_campaign_not_yet_run",
            "API_response_does_not_echo_immutable_server_checkpoint_identity",
            "API_metadata_request_is_not_logical_model_version_specific",
            "no_frozen_API_metadata_and_repeatability_fixture",
        ],
        "next_allowed_action": (
            "run_immutable_API_metadata_and_repeated_synthetic_score_parity_campaign"
            if credentials["api_credential_present"]
            else "provide_AlphaGenome_API_key_via_ALPHA_GENOME_API_KEY_without_logging_its_value"
        ),
    }
    write_json_exclusive(output / "receipt.json", receipt, mode=0o640)
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    receipt = activate(args.root, args.config, args.source_dir, args.output)
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
