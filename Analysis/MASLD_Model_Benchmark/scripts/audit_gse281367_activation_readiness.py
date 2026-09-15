#!/usr/bin/env python3
"""Freeze public metadata and task boundaries for ATAC-only GSE281367."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import time
import tomllib
import urllib.error
import urllib.parse
import urllib.request


MAX_SOURCE_BYTES = 32 * 1024 * 1024
TITLE_PATTERN = re.compile(r"^(MASH|Normal), replicate([1-6]), snATAC-seq$")
REQUIRED_GEO_COUNTS = {"MASH": 6, "Normal": 6}
REQUIRED_PHENOTYPES = {
    "source_condition",
    "mash_adjudication",
    "masld_diagnostic_criteria",
    "fibrosis_stage",
    "nas_total",
    "histology_components",
    "sex",
    "age",
    "bmi",
    "alcohol_and_other_etiology_exclusions",
}
REQUIRED_ACTIVE_TASKS = {"atac_fixed_window_donor_lineage_transport"}
REQUIRED_CONDITIONAL_TASKS = {
    "observed_atac_profile_and_count_tournament",
    "atac_embedding_and_clustering",
    "source_label_group_contrast",
}
REQUIRED_PROHIBITED_TASKS = {
    "rna_conditioned_atac_or_trans_expression_context",
    "paired_rna_atac_integration_or_translation",
    "histology_or_clinical_metadata_conditioning",
    "independent_cell_type_accuracy_against_source_labels",
    "external_sealed_universal_or_masld_champion",
}
REQUIRED_OVERLAPS = {
    "same_publication_gse281364_gse281160",
    "gse189600_annotation_reference",
    "other_project_liver_cohorts",
}


class ActivationReadinessError(RuntimeError):
    """Raised when public evidence or a fail-closed boundary changes."""


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def normalized_text(text: str) -> str:
    return " ".join(text.replace("\u00a0", " ").split())


def validate_contract(
    contract: dict[str, object],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    for key in (
        "automatic_activation",
        "sealed_outcomes_read",
        "biological_matrices_downloaded",
        "controlled_data_accessed",
        "rna_observed",
        "rna_invention_allowed",
        "same_cell_pairing_allowed",
        "same_nucleus_pairing_allowed",
        "cells_are_biological_replicates",
        "geo_mash_label_is_histology_adjudicated",
        "normal_label_is_adjudicated_masld_negative",
        "external_or_sealed_evaluation",
        "external_masld_champion_eligible",
        "cross_modality_donor_grouped_split_available",
    ):
        if contract.get(key) is not False:
            raise ActivationReadinessError(f"{key} must remain false")
    if contract.get("assay_native_donor_grouped_split_available") is not True:
        raise ActivationReadinessError("assay-native donor grouping must remain available")
    expected_counts = {
        "declared_participant_n": 12,
        "geo_assay_record_n": 12,
        "atac_assay_record_n": 12,
        "sra_run_n": 12,
        "sra_experiment_n": 12,
        "source_reported_high_quality_nucleus_n": 69595,
        "project_audited_atac_cell_n": 226224,
    }
    for key, expected in expected_counts.items():
        if contract.get(key) != expected:
            raise ActivationReadinessError(f"{key} must equal {expected}")
    if contract.get("role") != "project_exposed_development":
        raise ActivationReadinessError("role must remain project-exposed development")
    if contract.get("pairing_topology") != "same_study_unpaired_atac_only":
        raise ActivationReadinessError("ATAC-only topology changed")
    if contract.get("source_label_semantics") != "CONFLICT_GEO_MASH_ARTICLE_MASLD":
        raise ActivationReadinessError("source label conflict must remain explicit")
    sources = contract.get("source_condition")
    if not isinstance(sources, list) or len(sources) != 2:
        raise ActivationReadinessError("source condition comparison missing")
    by_source = {str(row["source"]): row for row in sources}
    if by_source.get("GEO", {}).get("case_label") != "MASH":
        raise ActivationReadinessError("GEO MASH label changed")
    if by_source.get("preprint_and_final_article", {}).get("case_label") != "MASLD":
        raise ActivationReadinessError("article MASLD label changed")
    if any(row.get("diagnostic_criteria_available") is not False for row in sources):
        raise ActivationReadinessError("source labels are not diagnostically adjudicated")
    phenotypes = contract.get("phenotype_field")
    if not isinstance(phenotypes, list) or {row.get("field_id") for row in phenotypes} != REQUIRED_PHENOTYPES:
        raise ActivationReadinessError("phenotype inventory changed")
    by_field = {str(row["field_id"]): row for row in phenotypes}
    if by_field["source_condition"].get("missingness_state") != "observed_conflicting_semantics":
        raise ActivationReadinessError("source label conflict was erased")
    for field_id in REQUIRED_PHENOTYPES - {"source_condition"}:
        if by_field[field_id].get("missingness_state") != "structurally_missing":
            raise ActivationReadinessError(f"{field_id} must remain structurally missing")
    tasks = contract.get("task_disposition")
    if not isinstance(tasks, list):
        raise ActivationReadinessError("task disposition missing")
    by_status = {
        status: {str(row["task_id"]) for row in tasks if row.get("eligibility") == status}
        for status in ("active_development", "conditional", "prohibited")
    }
    if by_status["active_development"] != REQUIRED_ACTIVE_TASKS:
        raise ActivationReadinessError("active task census changed")
    if by_status["conditional"] != REQUIRED_CONDITIONAL_TASKS:
        raise ActivationReadinessError("conditional task census changed")
    if by_status["prohibited"] != REQUIRED_PROHIBITED_TASKS:
        raise ActivationReadinessError("prohibited task census changed")
    overlaps = contract.get("cohort_overlap")
    if not isinstance(overlaps, list) or {row.get("overlap_id") for row in overlaps} != REQUIRED_OVERLAPS:
        raise ActivationReadinessError("cohort overlap inventory changed")
    if any(row.get("independent_confirmation_allowed") is True for row in overlaps):
        raise ActivationReadinessError("same-publication evidence cannot become independent confirmation")
    reference = contract.get("reference")
    if not isinstance(reference, dict):
        raise ActivationReadinessError("reference contract missing")
    if reference.get("source_native_peak_model_input_ready") is not False:
        raise ActivationReadinessError("source-native coordinates must fail closed")
    if reference.get("fixed_project_window_evaluator_ready") is not True:
        raise ActivationReadinessError("frozen evaluator readiness changed")
    decision = contract.get("decision")
    if not isinstance(decision, dict):
        raise ActivationReadinessError("decision block missing")
    if decision.get("fixed_window_transport_active") is not True:
        raise ActivationReadinessError("frozen transport must remain active")
    for key in (
        "source_native_model_input_active",
        "source_label_disease_task_active",
        "histology_conditioning_active",
        "rna_conditioning_active",
        "external_or_sealed_champion_available",
    ):
        if decision.get(key) is not False:
            raise ActivationReadinessError(f"decision {key} must remain false")
    if "gse289173" in json.dumps(contract).lower():
        raise ActivationReadinessError("sealed accession must not enter this audit")
    return phenotypes, tasks, overlaps


def fetch_source(source_id: str, url: str) -> tuple[bytes, dict[str, object]]:
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme != "https":
        raise ActivationReadinessError(f"non-HTTPS source rejected: {url}")
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json,text/csv,application/xml,text/plain,application/gzip,*/*",
            "User-Agent": "masld-bench-gse281367-public-metadata-audit/1.0",
        },
    )
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read(MAX_SOURCE_BYTES + 1)
                if len(payload) > MAX_SOURCE_BYTES:
                    raise ActivationReadinessError(f"source exceeds metadata cap: {source_id}")
                return payload, {
                    "source_id": source_id,
                    "requested_url": url,
                    "resolved_url": response.geturl(),
                    "http_status": response.status,
                    "retrieved_at_utc": utc_now(),
                    "http_date": response.headers.get("Date", ""),
                    "content_type": response.headers.get("Content-Type", ""),
                    "last_modified": response.headers.get("Last-Modified", ""),
                    "etag": response.headers.get("ETag", ""),
                    "bytes": len(payload),
                    "sha256": sha256_bytes(payload),
                    "attempt": attempt,
                }
        except (urllib.error.URLError, TimeoutError, ConnectionError) as error:
            last_error = error
            if attempt < 3:
                time.sleep(2 ** (attempt - 1))
    raise ActivationReadinessError(f"failed to fetch {source_id}: {last_error}") from last_error


def decode_source(payload: bytes, encoding: str) -> str:
    if encoding == "gzip_text":
        try:
            return gzip.decompress(payload).decode("utf-8")
        except (OSError, EOFError, UnicodeDecodeError) as error:
            raise ActivationReadinessError("invalid public gzip text") from error
    return payload.decode("utf-8", errors="replace")


def safe_filename(source_id: str, encoding: str) -> str:
    suffix = {"gzip_text": ".soft.gz", "text": ".txt", "json_text": ".json"}[encoding]
    return re.sub(r"[^a-z0-9_.-]+", "_", source_id.lower()) + suffix


def parse_geo_soft(payload: bytes) -> tuple[list[dict[str, object]], list[str]]:
    try:
        text = gzip.decompress(payload).decode("utf-8")
    except (OSError, EOFError, UnicodeDecodeError) as error:
        raise ActivationReadinessError("GEO SOFT is not valid gzip text") from error
    records: list[dict[str, object]] = []
    series_files: list[str] = []
    current: dict[str, object] | None = None
    for line in text.splitlines():
        if line.startswith("!Series_supplementary_file = "):
            series_files.append(line.split(" = ", 1)[1])
        if line.startswith("^SAMPLE = "):
            if current is not None:
                records.append(current)
            current = {
                "gsm": line.split(" = ", 1)[1],
                "relations": [],
                "supplementary_files": [],
            }
        elif current is not None and line.startswith("!Sample_title = "):
            current["title"] = line.split(" = ", 1)[1]
        elif current is not None and line.startswith("!Sample_library_strategy = "):
            current["library_strategy"] = line.split(" = ", 1)[1]
        elif current is not None and line.startswith("!Sample_relation = "):
            current["relations"].append(line.split(" = ", 1)[1])
        elif current is not None and re.match(r"!Sample_supplementary_file(?:_\d+)? = ", line):
            current["supplementary_files"].append(line.split(" = ", 1)[1])
    if current is not None:
        records.append(current)
    if len(records) != 12:
        raise ActivationReadinessError(f"expected 12 GEO records, found {len(records)}")
    output: list[dict[str, object]] = []
    for record in records:
        title = str(record.get("title", ""))
        match = TITLE_PATTERN.fullmatch(title)
        if match is None:
            raise ActivationReadinessError(f"unexpected GEO title: {title}")
        condition, replicate = match.groups()
        if record.get("library_strategy") != "ATAC-seq":
            raise ActivationReadinessError(f"unexpected assay strategy: {record['gsm']}")
        relations = list(record["relations"])
        biosamples = [value.rsplit("/", 1)[-1] for value in relations if "biosample/" in value.lower()]
        experiments = [value.split("term=", 1)[-1] for value in relations if "term=SRX" in value]
        if len(biosamples) != 1 or len(experiments) != 1:
            raise ActivationReadinessError(f"BioSample/SRX unresolved: {record['gsm']}")
        files = list(record["supplementary_files"])
        suffixes = {
            "barcodes": "_barcodes.tsv.gz",
            "fragments": "_fragments.tsv.gz",
            "matrix": "_matrix.mtx.gz",
            "peaks": "_peaks.bed.gz",
        }
        by_type = {
            file_type: [url for url in files if url.endswith(suffix)]
            for file_type, suffix in suffixes.items()
        }
        if any(len(matches) != 1 for matches in by_type.values()) or len(files) != 4:
            raise ActivationReadinessError(f"processed file quartet unresolved: {record['gsm']}")
        output.append(
            {
                "gsm": record["gsm"],
                "source_title": title,
                "geo_condition": condition,
                "replicate_number": int(replicate),
                "participant_unit": "one_distinct_article_declared_individual",
                "biosample": biosamples[0],
                "experiment": experiments[0],
                "processed_files": {key: values[0] for key, values in by_type.items()},
                "rna_state": "structurally_missing",
            }
        )
    counts = {condition: sum(row["geo_condition"] == condition for row in output) for condition in REQUIRED_GEO_COUNTS}
    if counts != REQUIRED_GEO_COUNTS:
        raise ActivationReadinessError(f"GEO condition counts changed: {counts}")
    for condition in REQUIRED_GEO_COUNTS:
        if {row["replicate_number"] for row in output if row["geo_condition"] == condition} != set(range(1, 7)):
            raise ActivationReadinessError(f"replicate numbering changed: {condition}")
    if len(set(series_files)) != 2:
        raise ActivationReadinessError("expected two GEO series archives")
    if not any(url.endswith("GSE281367_RAW.tar") for url in series_files):
        raise ActivationReadinessError("GEO RAW tar missing")
    if not any(url.endswith("GSE281367_seurat_clustered.rds.gz") for url in series_files):
        raise ActivationReadinessError("GEO clustered RDS missing")
    return output, sorted(series_files)


def parse_sra_runinfo(text: str, geo_records: list[dict[str, object]]) -> tuple[list[dict[str, str]], dict[str, object]]:
    rows = list(csv.DictReader(io.StringIO(text)))
    if len(rows) != 12:
        raise ActivationReadinessError(f"expected 12 SRA runs, found {len(rows)}")
    if len({row.get("Run") for row in rows}) != 12:
        raise ActivationReadinessError("SRA run IDs duplicated")
    if {row.get("Experiment") for row in rows} != {row["experiment"] for row in geo_records}:
        raise ActivationReadinessError("SRA experiments do not match GEO")
    if {row.get("BioSample") for row in rows} != {row["biosample"] for row in geo_records}:
        raise ActivationReadinessError("SRA BioSamples do not match GEO")
    if any(row.get("LibraryStrategy") != "ATAC-seq" for row in rows):
        raise ActivationReadinessError("non-ATAC SRA run detected")
    if any(row.get("LibraryLayout") != "PAIRED" for row in rows):
        raise ActivationReadinessError("non-paired SRA run detected")
    if any(row.get("Consent") != "public" for row in rows):
        raise ActivationReadinessError("non-public SRA consent state detected")
    for field in ("Sex", "Disease", "Histological_Type"):
        if any(row.get(field, "") for row in rows):
            raise ActivationReadinessError(f"SRA {field} is no longer blank; revise metadata audit")
    return rows, {
        "run_n": 12,
        "experiment_n": 12,
        "biosample_n": 12,
        "library_strategy": "ATAC-seq",
        "library_layout": "PAIRED",
        "all_consent_public": True,
        "reported_size_mb": sum(int(row["size_MB"]) for row in rows),
        "reported_spots": sum(int(row["spots"]) for row in rows),
        "reported_bases": sum(int(row["bases"]) for row in rows),
        "sex_disease_histology_fields_blank": True,
        "raw_reads_downloaded": False,
    }


def parse_bioc(payload: bytes) -> dict[str, object]:
    try:
        collection = json.loads(payload)
        document = collection[0]["documents"][0]
    except (json.JSONDecodeError, IndexError, KeyError, TypeError) as error:
        raise ActivationReadinessError("invalid preprint BioC record") from error
    if document.get("id") != "PMC12633503" or document.get("infons", {}).get("license") != "CC BY":
        raise ActivationReadinessError("preprint identity or article license changed")
    passages = document.get("passages", [])
    texts = [str(row.get("text", "")) for row in passages]
    joined = normalized_text(" ".join(texts))
    required = [
        "liver nuclei from 12 individuals (6 Normal and 6 MASLD)",
        "69,595 high-quality nuclei",
        "Samples from different patients were processed separately",
        "downloaded from GSE189600",
        "GSE281367, GSE281364, and GSE281160",
    ]
    missing = [token for token in required if normalized_text(token).lower() not in joined.lower()]
    if missing:
        raise ActivationReadinessError(f"BioC mechanistic evidence missing: {missing}")
    return {
        "document_id": document["id"],
        "preprint_article_license": "CC_BY_4.0",
        "participant_n": 12,
        "article_condition_counts": {"Normal": 6, "MASLD": 6},
        "source_reported_high_quality_nucleus_n": 69595,
        "samples_from_different_patients_processed_separately": True,
        "human_liver_sources": ["Liver Tissue Cell Distribution System", "SEKISUI XenoTech"],
        "participant_diagnostic_criteria_reported": False,
        "participant_histology_reported": False,
        "rna_annotation_reference": "GSE189600",
        "same_publication_accessions": ["GSE281367", "GSE281364", "GSE281160"],
    }


def parse_crossref(payload: bytes) -> dict[str, object]:
    try:
        message = json.loads(payload)["message"]
    except (json.JSONDecodeError, KeyError, TypeError) as error:
        raise ActivationReadinessError("invalid Crossref record") from error
    relation = message.get("relation", {}).get("has-preprint", [])
    ids = {row.get("id") for row in relation}
    if "10.21203/rs.3.rs-6984670/v1" not in ids:
        raise ActivationReadinessError("final-to-preprint relation missing")
    return {
        "final_doi": "10.1038/s41588-026-02617-8",
        "preprint_doi": "10.21203/rs.3.rs-6984670/v1",
        "preprint_relation_frozen": True,
        "crossref_license_entries": message.get("license", []),
    }


def validate_project_artifacts(repo_root: Path, contract: dict[str, object]) -> dict[str, object]:
    spec = contract["project_artifact"]
    source_root = repo_root / str(spec["source_audit_path"])
    membership_root = repo_root / str(spec["membership_path"])
    outcome_root = repo_root / str(spec["outcome_path"])
    checks = [
        (source_root / "ARTIFACTS.json", str(spec["source_audit_artifacts_sha256"])),
        (source_root / "audit.json", str(spec["source_audit_json_sha256"])),
        (membership_root / "ARTIFACTS.json", str(spec["membership_artifacts_sha256"])),
        (membership_root / "membership.json", str(spec["membership_json_sha256"])),
        (outcome_root / "ARTIFACTS.json", str(spec["outcome_artifacts_sha256"])),
        (outcome_root / "outcome_summary.json", str(spec["outcome_summary_sha256"])),
    ]
    for path, expected in checks:
        if not path.is_file() or sha256_path(path) != expected:
            raise ActivationReadinessError(f"frozen project artifact mismatch: {path}")
    audit = json.loads((source_root / "audit.json").read_text(encoding="utf-8"))
    membership = json.loads((membership_root / "membership.json").read_text(encoding="utf-8"))
    outcome = json.loads((outcome_root / "outcome_summary.json").read_text(encoding="utf-8"))
    atac = audit.get("h5ad", {}).get("gse281367", {})
    if atac.get("donors") != 12 or atac.get("cells") != 226224:
        raise ActivationReadinessError("project ATAC census changed")
    if membership.get("pairing_by_dataset", {}).get("gse281367") != "same_study_unpaired":
        raise ActivationReadinessError("project pairing topology changed")
    if outcome.get("donors", {}).get("gse281367") != 12 or outcome.get("windows") != 32000:
        raise ActivationReadinessError("fixed-window outcome authority changed")
    coordinate = outcome.get("coordinate_contract_by_dataset", {}).get("gse281367", {})
    if coordinate.get("cut_sites") != ["fragment_start", "fragment_end_minus_1"]:
        raise ActivationReadinessError("Cell Ranger fragment semantics changed")
    return {
        "source_audit_artifacts_sha256": spec["source_audit_artifacts_sha256"],
        "membership_artifacts_sha256": spec["membership_artifacts_sha256"],
        "outcome_artifacts_sha256": spec["outcome_artifacts_sha256"],
        "project_audited_atac_cell_n": atac["cells"],
        "project_participant_n": atac["donors"],
        "project_condition_labels": atac["condition_counts"],
        "pairing_topology": membership["pairing_by_dataset"]["gse281367"],
        "fixed_window_n": outcome["windows"],
        "fixed_window_transport_active": outcome["dataset_activated_for_development_transport"],
        "cell_ranger_fragments_already_tn5_adjusted": True,
        "biological_matrix_read_by_this_audit": False,
    }


def write_run_inventory(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "Run", "Experiment", "BioSample", "SampleName", "LibraryStrategy",
        "LibraryLayout", "Platform", "Model", "spots", "bases", "size_MB", "Consent",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    contract_path = args.contract.resolve(strict=True)
    output = args.output.resolve()
    if output.exists():
        raise ActivationReadinessError(f"refusing to overwrite output: {output}")
    contract = tomllib.loads(contract_path.read_text(encoding="utf-8"))
    phenotypes, tasks, overlaps = validate_contract(contract)
    repo_root = Path(__file__).resolve().parents[1]
    project_authority = validate_project_artifacts(repo_root, contract)
    metadata_dir = output / "public_metadata"
    metadata_dir.mkdir(parents=True)

    payloads: dict[str, bytes] = {}
    texts: dict[str, str] = {}
    receipts: list[dict[str, object]] = []
    for spec in contract["public_source"]:
        source_id = str(spec["source_id"])
        encoding = str(spec["encoding"])
        payload, receipt = fetch_source(source_id, str(spec["url"]))
        text = decode_source(payload, encoding)
        evidence = normalized_text(text).lower()
        missing = [
            str(token) for token in spec["required_tokens"]
            if normalized_text(str(token)).lower() not in evidence
        ]
        if missing:
            raise ActivationReadinessError(f"{source_id} missing required evidence: {missing}")
        target = metadata_dir / safe_filename(source_id, encoding)
        target.write_bytes(payload)
        receipt["payload_stored"] = True
        receipt["local_path"] = str(target.relative_to(output))
        receipt["required_evidence_pass"] = True
        receipt["required_tokens"] = list(spec["required_tokens"])
        payloads[source_id] = payload
        texts[source_id] = text
        receipts.append(receipt)

    geo_records, series_files = parse_geo_soft(payloads["geo_soft"])
    sra_rows, sra_summary = parse_sra_runinfo(texts["sra_runinfo"], geo_records)
    article_evidence = parse_bioc(payloads["preprint_bioc"])
    publication_relation = parse_crossref(payloads["final_crossref"])
    write_run_inventory(output / "raw_run_inventory.tsv", sra_rows)

    topology = {
        "schema_version": "masld-bench-gse281367-participant-assay-topology-v1",
        "participant_n": 12,
        "geo_assay_record_n": 12,
        "sra_run_n": 12,
        "modality": "single_nucleus_atac",
        "pairing_topology": "same_study_unpaired_atac_only",
        "samples_from_different_patients_processed_separately": True,
        "assay_native_donor_grouped_split_available": True,
        "rna_state": "structurally_missing",
        "same_cell_or_nucleus_cross_modality_pairing": False,
        "cells_are_biological_replicates": False,
        "geo_records": geo_records,
        "article_evidence": article_evidence,
        "project_authority": project_authority,
        "cell_census_discrepancy": {
            "article_high_quality_nuclei": 69595,
            "project_audited_cells": 226224,
            "state": "distinct_qc_or_membership_authorities_not_interchangeable",
            "source_final_high_quality_barcode_join_frozen": False,
        },
    }
    phenotype_inventory = {
        "schema_version": "masld-bench-gse281367-phenotype-histology-v1",
        "participant_n": 12,
        "source_label_conflict": {
            "GEO": {"Normal": 6, "MASH": 6},
            "preprint_and_final_article": {"Normal": 6, "MASLD": 6},
            "resolved": False,
            "mash_histology_adjudicated": False,
            "normal_adjudicated_masld_negative": False,
        },
        "fields": phenotypes,
        "sra_sex_disease_histology_fields_blank": True,
        "disease_group_may_proxy_fibrosis_nas_or_histology": False,
        "missing_metadata_encoded_as_zero": False,
        "biological_matrices_read": False,
        "sealed_outcomes_read": False,
    }
    file_rights_reference = {
        "schema_version": "masld-bench-gse281367-files-rights-reference-v1",
        "raw_archive": sra_summary,
        "processed_files": {
            "per_sample_file_n": sum(len(row["processed_files"]) for row in geo_records),
            "per_sample_file_type_n": {key: 12 for key in ("barcodes", "fragments", "matrix", "peaks")},
            "series_archives": [{"url": url, "requested": False} for url in series_files],
            "per_sample_files": [
                {"gsm": row["gsm"], "files": row["processed_files"], "requested": False}
                for row in geo_records
            ],
            "biological_files_requested_by_this_audit": False,
        },
        "rights": contract["rights"],
        "reference": contract["reference"],
        "source_native_coordinate_model_input_ready": False,
        "fixed_project_window_evaluator_ready": True,
    }
    overlap = {
        "schema_version": "masld-bench-gse281367-cohort-family-overlap-v1",
        "cohort_family_id": contract["cohort_family_id"],
        "publication_family_id": contract["publication_family_id"],
        "publication_relation": publication_relation,
        "overlaps": overlaps,
        "participant_overlap_with_other_project_liver_cohorts": "UNRESOLVED_PUBLIC_METADATA",
        "cross_cohort_genotype_or_expression_fingerprint_run": False,
        "may_claim_independent_donor_family": False,
        "gse189600_source_label_dependency": True,
        "same_publication_evidence_counts_as_independent_confirmation": False,
    }
    task_eligibility = {
        "schema_version": "masld-bench-gse281367-task-eligibility-v1",
        "current_role": "project_exposed_atac_only_development",
        "task_dispositions": tasks,
        "fixed_window_transport_active": True,
        "source_native_model_input_active": False,
        "source_label_disease_task_active": False,
        "histology_conditioning_active": False,
        "rna_conditioning_active": False,
        "external_or_sealed_champion_available": False,
        "next_gate": contract["decision"]["next_gate"],
    }

    (output / "source_receipts.json").write_bytes(canonical_json({"sources": receipts}))
    (output / "participant_assay_topology.json").write_bytes(canonical_json(topology))
    (output / "phenotype_histology_inventory.json").write_bytes(canonical_json(phenotype_inventory))
    (output / "file_rights_reference_readiness.json").write_bytes(canonical_json(file_rights_reference))
    (output / "cohort_family_overlap.json").write_bytes(canonical_json(overlap))
    (output / "task_eligibility.json").write_bytes(canonical_json(task_eligibility))
    (output / "DECISION.md").write_text(
        "# GSE281367 activation-readiness decision\n\n"
        "**Decision: retain the frozen ATAC transport endpoint; block disease-label "
        "and source-native model-input activation.** Public evidence confirms 12 "
        "distinct individuals and 12 ATAC libraries. There is no RNA assay, paired "
        "modality, histology table, or participant clinical metadata.\n\n"
        "The phenotype is not resolved: GEO labels six samples `MASH` and six "
        "`Normal`, whereas the linked preprint and final article describe six "
        "`MASLD` and six `Normal`. Neither source gives participant diagnostic, "
        "fibrosis, NAS, histology, alcohol-exclusion, age, sex, or BMI criteria. "
        "These labels can support only explicitly source-labeled descriptive work; "
        "they cannot be treated as adjudicated MASH or MASLD-negative controls.\n\n"
        "The paper reports 69,595 high-quality nuclei, while the frozen project "
        "authority contains 226,224 cells. Those are distinct QC/membership "
        "authorities and cannot be interchanged until an exact source-final barcode "
        "join is frozen. The existing 12-participant, 32,000-window evaluator-only "
        "transport remains active. Source-native training waits for the exact "
        "Cell Ranger reference, coordinate crosswalk, QC membership, and label audit.\n\n"
        "GSE281367 shares a publication family with GSE281364 and GSE281160, and its "
        "cell labels were transferred from GSE189600 RNA. These are within-study and "
        "annotation dependencies, not independent confirmation. Participant overlap "
        "with other project liver cohorts remains unresolved because no stable public "
        "donor aliases or genotype fingerprints are available. No biological matrix, "
        "controlled data, or sealed outcome was accessed.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
