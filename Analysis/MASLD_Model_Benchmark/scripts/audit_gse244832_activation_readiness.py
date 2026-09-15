#!/usr/bin/env python3
"""Freeze public metadata and task boundaries for GSE244832 without matrices."""

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
import shutil
import subprocess
import time
import tomllib
import urllib.error
import urllib.parse
import urllib.request


MAX_SOURCE_BYTES = 32 * 1024 * 1024
ATAC_TITLE = re.compile(r"^(MM_\d+), snATACseq$")
RNA_TITLE = re.compile(r"^(JB_\d+), snRNAseq$")
REQUIRED_CONDITIONS = {"NORMAL": 5, "MASL": 4, "MASH": 9}
REQUIRED_PHENOTYPES = {
    "disease_group",
    "fibrosis_stage",
    "crn_score_or_diagnostic_criteria",
    "nas_total",
    "histology_components",
    "sex",
    "age",
    "bmi_or_body_composition",
}
REQUIRED_ACTIVE_TASKS = {"atac_fixed_window_donor_lineage_transport"}
REQUIRED_CONDITIONAL_TASKS = {
    "assay_native_donor_pseudobulk",
    "rna_cell_state_or_group_mapping",
    "cross_assay_donor_pseudobulk_late_fusion",
}
REQUIRED_PROHIBITED_TASKS = {
    "same_cell_or_same_nucleus_rna_atac_pairing",
    "histology_conditioned_multimodal_training",
    "external_or_sealed_masld_champion",
}


class ActivationReadinessError(RuntimeError):
    """Raised when public evidence or the fail-closed requirements change."""


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


def validate_contract(contract: dict[str, object]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    for key in (
        "automatic_activation",
        "sealed_outcomes_read",
        "biological_matrices_downloaded",
        "controlled_data_accessed",
        "same_cell_pairing_allowed",
        "same_nucleus_pairing_allowed",
        "infer_cross_assay_join_by_record_order",
        "cells_are_biological_replicates",
        "external_or_sealed_evaluation",
        "external_masld_champion_eligible",
        "cross_assay_donor_grouped_split_available",
    ):
        if contract.get(key) is not False:
            raise ActivationReadinessError(f"{key} must remain false")
    for key in ("assay_native_donor_grouped_split_available",):
        if contract.get(key) is not True:
            raise ActivationReadinessError(f"{key} must remain true")
    expected_counts = {
        "declared_participant_n": 18,
        "geo_assay_record_n": 36,
        "rna_assay_record_n": 18,
        "atac_assay_record_n": 18,
        "sra_run_n": 135,
        "sra_experiment_n": 36,
        "rna_sra_run_n": 117,
        "atac_sra_run_n": 18,
        "project_atac_participant_n": 18,
        "project_atac_cell_n": 88814,
    }
    for key, expected in expected_counts.items():
        if contract.get(key) != expected:
            raise ActivationReadinessError(f"{key} must equal {expected}")
    if contract.get("role") != "project_exposed_development":
        raise ActivationReadinessError("GSE244832 must remain project-exposed development")
    if contract.get("pairing_topology") != "same_sample_different_aliquot":
        raise ActivationReadinessError("pairing topology changed")
    if contract.get("study_level_same_liver_relationship") != "DECLARED":
        raise ActivationReadinessError("study-level same-liver declaration is required")
    if contract.get("machine_readable_cross_assay_participant_join") != "UNRESOLVED":
        raise ActivationReadinessError("cross-assay participant join must fail closed")
    conditions = contract.get("condition")
    if not isinstance(conditions, list):
        raise ActivationReadinessError("condition census missing")
    if {row.get("condition_id"): row.get("participant_n") for row in conditions} != REQUIRED_CONDITIONS:
        raise ActivationReadinessError("condition census must be NORMAL=5, MASL=4, MASH=9")
    phenotypes = contract.get("phenotype_field")
    if not isinstance(phenotypes, list):
        raise ActivationReadinessError("phenotype inventory missing")
    if {row.get("field_id") for row in phenotypes} != REQUIRED_PHENOTYPES:
        raise ActivationReadinessError("phenotype inventory changed")
    by_field = {str(row["field_id"]): row for row in phenotypes}
    if by_field["disease_group"].get("missingness_state") != "observed":
        raise ActivationReadinessError("disease group must be observed")
    for field_id in ("nas_total", "histology_components", "sex", "age", "bmi_or_body_composition"):
        if by_field[field_id].get("missingness_state") != "structurally_missing":
            raise ActivationReadinessError(f"{field_id} must not be imputed or encoded as zero")
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
    cross_assay = next(row for row in tasks if row["task_id"] == "cross_assay_donor_pseudobulk_late_fusion")
    if cross_assay.get("model_selection_allowed") is not False:
        raise ActivationReadinessError("unjoined cross-assay fusion cannot select models")
    partial = contract.get("partial_fibrosis_evidence")
    if not isinstance(partial, dict):
        raise ActivationReadinessError("partial fibrosis evidence missing")
    stage_3 = set(partial.get("exact_stage_3_donor_ids", []))
    stage_4 = set(partial.get("exact_stage_4_donor_ids", []))
    ambiguous = set(partial.get("ambiguous_stage_2_or_3_donor_ids", []))
    if stage_3 != {"D14", "D16"} or stage_4 != {"D11", "D12", "D15", "D18"}:
        raise ActivationReadinessError("exact fibrosis evidence changed")
    if ambiguous != {"D10", "D13", "D17"}:
        raise ActivationReadinessError("ambiguous fibrosis evidence changed")
    if stage_3 & stage_4 or stage_3 & ambiguous or stage_4 & ambiguous:
        raise ActivationReadinessError("fibrosis donor groups overlap")
    if partial.get("may_be_joined_to_geo_assay_rows") is not False:
        raise ActivationReadinessError("article donor IDs cannot be joined to GEO rows by assumption")
    decision = contract.get("decision")
    if not isinstance(decision, dict):
        raise ActivationReadinessError("decision block missing")
    for key in (
        "cross_assay_donor_fusion_available",
        "same_cell_or_nucleus_pairing_available",
        "histology_conditioning_available",
        "external_or_sealed_champion_available",
    ):
        if decision.get(key) is not False:
            raise ActivationReadinessError(f"decision {key} must remain false")
    if "gse289173" in json.dumps(contract).lower():
        raise ActivationReadinessError("sealed accession must not enter this audit")
    return phenotypes, tasks


def fetch_source(source_id: str, url: str) -> tuple[bytes, dict[str, object]]:
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme != "https":
        raise ActivationReadinessError(f"non-HTTPS public source rejected: {url}")
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/pdf,text/csv,application/xml,text/plain,application/gzip,*/*",
            "User-Agent": "masld-bench-gse244832-public-metadata-audit/1.0",
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


def decode_source(payload: bytes, encoding: str) -> tuple[str, dict[str, object]]:
    if encoding == "gzip_text":
        try:
            return gzip.decompress(payload).decode("utf-8"), {}
        except (OSError, EOFError, UnicodeDecodeError) as error:
            raise ActivationReadinessError("invalid public gzip text") from error
    if encoding == "pdf_text":
        binary = shutil.which("pdftotext")
        if binary is None:
            raise ActivationReadinessError("pdftotext is required for article evidence")
        completed = subprocess.run(
            [binary, "-layout", "-", "-"],
            input=payload,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if completed.returncode != 0:
            raise ActivationReadinessError(
                f"pdftotext failed: {completed.stderr.decode(errors='replace').strip()}"
            )
        text = completed.stdout.decode("utf-8", errors="replace")
        version = subprocess.run(
            [binary, "-v"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False
        )
        version_text = (version.stderr or version.stdout).decode(errors="replace").splitlines()
        return text, {
            "text_extractor": binary,
            "text_extractor_version": version_text[0] if version_text else "",
            "derived_text_sha256": sha256_bytes(completed.stdout),
            "derived_text_stored": False,
        }
    return payload.decode("utf-8", errors="replace"), {}


def safe_filename(source_id: str, encoding: str) -> str:
    suffix = {"gzip_text": ".soft.gz", "text": ".txt", "pdf_text": ".pdf"}[encoding]
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
        elif current is not None and line.startswith("!Sample_description = "):
            current["condition"] = line.split(" = ", 1)[1]
        elif current is not None and line.startswith("!Sample_library_strategy = "):
            current["library_strategy"] = line.split(" = ", 1)[1]
        elif current is not None and line.startswith("!Sample_relation = "):
            current["relations"].append(line.split(" = ", 1)[1])
        elif current is not None and re.match(r"!Sample_supplementary_file(?:_\d+)? = ", line):
            current["supplementary_files"].append(line.split(" = ", 1)[1])
    if current is not None:
        records.append(current)
    if len(records) != 36:
        raise ActivationReadinessError(f"expected 36 GEO records, found {len(records)}")
    output: list[dict[str, object]] = []
    for record in records:
        strategy = str(record.get("library_strategy", ""))
        title = str(record.get("title", ""))
        if strategy == "ATAC-seq":
            match = ATAC_TITLE.fullmatch(title)
            modality = "ATAC"
        elif strategy == "RNA-Seq":
            match = RNA_TITLE.fullmatch(title)
            modality = "RNA"
        else:
            raise ActivationReadinessError(f"unexpected library strategy: {strategy}")
        if match is None:
            raise ActivationReadinessError(f"unexpected assay title: {title}")
        condition = str(record.get("condition", ""))
        if condition not in REQUIRED_CONDITIONS:
            raise ActivationReadinessError(f"unexpected condition: {condition}")
        relations = list(record["relations"])
        biosamples = [value.rsplit("/", 1)[-1] for value in relations if "biosample/" in value.lower()]
        experiments = [value.split("term=", 1)[-1] for value in relations if "term=SRX" in value]
        if len(biosamples) != 1 or len(experiments) != 1:
            raise ActivationReadinessError(f"BioSample/SRX relation unresolved: {record['gsm']}")
        supplementary = list(record["supplementary_files"])
        if modality == "ATAC":
            if len(supplementary) != 1 or not supplementary[0].endswith(".sorted.bed.gz"):
                raise ActivationReadinessError(f"ATAC processed file unresolved: {record['gsm']}")
            processed_file = supplementary[0]
        else:
            if supplementary not in ([], ["NONE"]):
                raise ActivationReadinessError(f"unexpected RNA sample file: {record['gsm']}")
            processed_file = None
        output.append(
            {
                "gsm": record["gsm"],
                "source_title": title,
                "source_alias": match.group(1),
                "condition": condition,
                "modality": modality,
                "library_strategy": strategy,
                "biosample": biosamples[0],
                "experiment": experiments[0],
                "processed_file": processed_file,
                "cross_assay_participant_id": None,
                "cross_assay_join_state": "unresolved_not_inferred_from_record_order",
            }
        )
    for modality in ("ATAC", "RNA"):
        subset = [row for row in output if row["modality"] == modality]
        if len(subset) != 18:
            raise ActivationReadinessError(f"expected 18 {modality} records")
        counts = {condition: sum(row["condition"] == condition for row in subset) for condition in REQUIRED_CONDITIONS}
        if counts != REQUIRED_CONDITIONS:
            raise ActivationReadinessError(f"{modality} condition counts changed: {counts}")
    if len(set(series_files)) != 2:
        raise ActivationReadinessError("expected two GEO series supplementary archives")
    if not any(path.endswith("GSE244832_RAW.tar") for path in series_files):
        raise ActivationReadinessError("GEO RAW supplementary archive missing")
    if not any(path.endswith("GSE244832_hLIVER_processed_files.tar.gz") for path in series_files):
        raise ActivationReadinessError("GEO processed RNA archive missing")
    if {row["source_alias"] for row in output if row["modality"] == "ATAC"} & {
        row["source_alias"] for row in output if row["modality"] == "RNA"
    }:
        raise ActivationReadinessError("unexpected exact assay alias overlap; revise join audit")
    return output, sorted(series_files)


def parse_sra_runinfo(text: str, geo_records: list[dict[str, object]]) -> tuple[list[dict[str, str]], dict[str, object]]:
    rows = list(csv.DictReader(io.StringIO(text)))
    if len(rows) != 135:
        raise ActivationReadinessError(f"expected 135 SRA runs, found {len(rows)}")
    runs = [row.get("Run", "") for row in rows]
    if len(set(runs)) != 135 or any(not run.startswith("SRR") for run in runs):
        raise ActivationReadinessError("SRA run IDs are missing or duplicated")
    experiments = {row.get("Experiment", "") for row in rows}
    biosamples = {row.get("BioSample", "") for row in rows}
    if len(experiments) != 36 or len(biosamples) != 36:
        raise ActivationReadinessError("SRA experiment/BioSample census changed")
    if experiments != {str(row["experiment"]) for row in geo_records}:
        raise ActivationReadinessError("SRA experiments do not match GEO")
    if biosamples != {str(row["biosample"]) for row in geo_records}:
        raise ActivationReadinessError("SRA BioSamples do not match GEO")
    strategy_counts = {
        strategy: sum(row.get("LibraryStrategy") == strategy for row in rows)
        for strategy in ("ATAC-seq", "RNA-Seq")
    }
    if strategy_counts != {"ATAC-seq": 18, "RNA-Seq": 117}:
        raise ActivationReadinessError(f"SRA strategy counts changed: {strategy_counts}")
    if any(row.get("Consent") != "public" for row in rows):
        raise ActivationReadinessError("not every SRA run is marked public")
    if any(row.get("LibraryLayout") != "PAIRED" for row in rows):
        raise ActivationReadinessError("not every SRA run is paired-end")
    summary = {
        "run_n": len(rows),
        "experiment_n": len(experiments),
        "biosample_n": len(biosamples),
        "run_n_by_library_strategy": strategy_counts,
        "reported_size_mb_by_library_strategy": {
            strategy: sum(int(row["size_MB"]) for row in rows if row["LibraryStrategy"] == strategy)
            for strategy in ("ATAC-seq", "RNA-Seq")
        },
        "all_consent_public": True,
        "all_library_layout_paired": True,
        "raw_reads_downloaded": False,
    }
    return rows, summary


def validate_project_artifacts(repo_root: Path, contract: dict[str, object]) -> dict[str, object]:
    spec = contract["project_artifact"]
    source_root = repo_root / str(spec["source_audit_path"])
    membership_root = repo_root / str(spec["membership_path"])
    checks = [
        (source_root / "ARTIFACTS.json", str(spec["source_audit_artifacts_sha256"])),
        (source_root / "audit.json", str(spec["source_audit_json_sha256"])),
        (membership_root / "ARTIFACTS.json", str(spec["membership_artifacts_sha256"])),
        (membership_root / "membership.json", str(spec["membership_json_sha256"])),
    ]
    for path, expected in checks:
        if not path.is_file() or sha256_path(path) != expected:
            raise ActivationReadinessError(f"frozen project artifact mismatch: {path}")
    audit = json.loads((source_root / "audit.json").read_text(encoding="utf-8"))
    membership = json.loads((membership_root / "membership.json").read_text(encoding="utf-8"))
    atac = audit.get("h5ad", {}).get("gse244832", {})
    if atac.get("donors") != 18 or atac.get("cells") != 88814:
        raise ActivationReadinessError("frozen ATAC audit census changed")
    if membership.get("pairing_by_dataset", {}).get("gse244832") != "same_sample_different_aliquot":
        raise ActivationReadinessError("frozen membership topology changed")
    return {
        "source_audit_path": str(spec["source_audit_path"]),
        "source_audit_artifacts_sha256": spec["source_audit_artifacts_sha256"],
        "membership_path": str(spec["membership_path"]),
        "membership_artifacts_sha256": spec["membership_artifacts_sha256"],
        "atac_participant_n": atac["donors"],
        "atac_cell_n": atac["cells"],
        "pairing_topology": membership["pairing_by_dataset"]["gse244832"],
        "biological_matrix_read_by_this_audit": False,
    }


def write_run_inventory(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "Run",
        "Experiment",
        "BioSample",
        "SampleName",
        "LibraryStrategy",
        "LibraryLayout",
        "Platform",
        "Model",
        "spots",
        "bases",
        "size_MB",
        "Consent",
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
    phenotypes, tasks = validate_contract(contract)
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
        source_text, extractor = decode_source(payload, encoding)
        evidence_text = normalized_text(source_text).lower()
        missing = [
            str(token)
            for token in spec["required_tokens"]
            if normalized_text(str(token)).lower() not in evidence_text
        ]
        if missing:
            raise ActivationReadinessError(f"{source_id} missing required evidence: {missing}")
        receipt.update(extractor)
        receipt["required_evidence_pass"] = True
        receipt["required_tokens"] = list(spec["required_tokens"])
        if spec["store_payload"] is True:
            target = metadata_dir / safe_filename(source_id, encoding)
            target.write_bytes(payload)
            receipt["payload_stored"] = True
            receipt["local_path"] = str(target.relative_to(output))
        else:
            receipt["payload_stored"] = False
            receipt["local_path"] = None
        payloads[source_id] = payload
        texts[source_id] = source_text
        receipts.append(receipt)

    geo_records, series_files = parse_geo_soft(payloads["geo_soft"])
    sra_rows, sra_summary = parse_sra_runinfo(texts["sra_runinfo"], geo_records)
    write_run_inventory(output / "raw_run_inventory.tsv", sra_rows)

    topology = {
        "schema_version": "masld-bench-gse244832-participant-assay-topology-v1",
        "declared_participant_n": 18,
        "condition_participant_n": REQUIRED_CONDITIONS,
        "geo_assay_record_n": 36,
        "rna_assay_record_n": 18,
        "atac_assay_record_n": 18,
        "pairing_topology": "same_sample_different_aliquot",
        "study_level_same_liver_relationship": "declared_by_series_and_primary_article",
        "machine_readable_cross_assay_participant_join": "unresolved",
        "assay_title_namespaces": {"ATAC": "MM_*", "RNA": "JB_*"},
        "record_order_may_define_join": False,
        "same_cell_pairing": False,
        "same_nucleus_pairing": False,
        "assay_native_donor_grouped_split_available": True,
        "cross_assay_donor_grouped_split_available": False,
        "cells_are_biological_replicates": False,
        "geo_records": geo_records,
        "frozen_project_atac_authority": project_authority,
    }
    phenotype_inventory = {
        "schema_version": "masld-bench-gse244832-phenotype-histology-inventory-v1",
        "participant_n": 18,
        "condition_participant_n": REQUIRED_CONDITIONS,
        "fields": phenotypes,
        "partial_fibrosis_evidence": contract["partial_fibrosis_evidence"],
        "source_interpretation": {
            "normal": "source-defined CRN criteria below 3, not a general population control claim",
            "masl": "source-defined steatosis group",
            "mash": "source-defined MASH; nine participants span fibrosis stages 2-4",
            "nas_rule": "do not reinterpret the source CRN wording as observed participant NAS",
            "missing_rule": "unavailable metadata remains structurally missing or join-unresolved, never zero",
        },
        "biological_matrices_read": False,
        "supplementary_article_files_downloaded": False,
    }
    files_rights_reference = {
        "schema_version": "masld-bench-gse244832-files-rights-reference-v1",
        "raw_archive": sra_summary,
        "processed_files": {
            "atac_per_sample_bed_n": sum(row["processed_file"] is not None for row in geo_records),
            "rna_combined_processed_archive_n": 1,
            "geo_series_supplementary_archives": [
                {"url": url, "requested": False} for url in series_files
            ],
            "atac_per_sample_urls": [
                {"gsm": row["gsm"], "url": row["processed_file"], "requested": False}
                for row in geo_records
                if row["modality"] == "ATAC"
            ],
            "biological_files_requested_by_this_audit": False,
        },
        "rights": contract["rights"],
        "reference": contract["reference"],
        "reference_decision": "high_level_hg38_only_exact_reference_and_annotation_crosswalks_unresolved",
        "sequence_or_coordinate_task_ready": False,
    }
    task_eligibility = {
        "schema_version": "masld-bench-gse244832-task-eligibility-v1",
        "current_role": "project_exposed_development",
        "task_dispositions": tasks,
        "assay_native_donor_tasks_available": True,
        "cross_assay_donor_fusion_available": False,
        "same_cell_or_nucleus_pairing_available": False,
        "histology_conditioning_available": False,
        "external_or_sealed_champion_available": False,
        "next_gate": contract["decision"]["next_gate"],
        "sealed_outcomes_read": False,
        "biological_matrices_downloaded": False,
    }

    (output / "source_receipts.json").write_bytes(canonical_json({"sources": receipts}))
    (output / "participant_assay_topology.json").write_bytes(canonical_json(topology))
    (output / "phenotype_histology_inventory.json").write_bytes(canonical_json(phenotype_inventory))
    (output / "file_rights_reference_readiness.json").write_bytes(canonical_json(files_rights_reference))
    (output / "task_eligibility.json").write_bytes(canonical_json(task_eligibility))
    (output / "DECISION.md").write_text(
        "# GSE244832 activation-readiness decision\n\n"
        "**Decision: admit assay-native participant-level development tasks; keep "
        "cross-assay fusion blocked.** Public metadata confirms 18 livers: 5 NORMAL, "
        "4 MASL, and 9 MASH. Each liver contributed separate RNA and ATAC nuclei, so "
        "the topology is same-sample different-aliquot, never same-cell or same-nucleus.\n\n"
        "GEO exposes 18 `MM_*` ATAC records and 18 `JB_*` RNA records but no "
        "authoritative row-level crosswalk between those namespaces. Record order is "
        "not a donor join. The existing 18-participant, 88,814-cell ATAC transport "
        "authority remains valid for project-exposed development; RNA-only donor "
        "pseudobulk can be admitted separately after its matrix and reference gates. "
        "Cross-assay participant fusion waits for a frozen authoritative join.\n\n"
        "Disease group is observed. The primary article gives partial MASH fibrosis "
        "stage evidence, but complete participant fibrosis, NAS, histology components, "
        "sex, age, and BMI are not established in the audited sources and cannot be "
        "encoded as zero or inferred from group. The source is project-exposed and "
        "cannot support sealed, external, universal, diagnostic, or clinical claims. "
        "No biological matrix, controlled data, sealed outcome, or supplementary "
        "article file was accessed.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
