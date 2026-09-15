#!/usr/bin/env python3
"""Freeze the public pooling topology and task boundary for GSE281574."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import time
import tomllib
import urllib.error
import urllib.parse
import urllib.request


MAX_SOURCE_BYTES = 32 * 1024 * 1024
TITLE_PATTERN = re.compile(r"^(Normal|SAH|AC)_pool([12])_(ATAC|RNA)$")
REQUIRED_CONDITIONS = {"Normal", "SAH", "AC"}
REQUIRED_ALLOWED_TASKS = {
    "same_nucleus_rna_atac_pairing",
    "leave_one_library_out_fixed_reconstruction",
    "alcohol_etiology_ood_mapping",
    "library_level_condition_contrast",
}
REQUIRED_PROHIBITED_TASKS = {
    "donor_level_training_or_evaluation",
    "masld_case_control_or_negative_control",
    "histology_or_metadata_conditioning",
    "external_masld_champion_or_sealed_confirmation",
}


class PoolingTopologyError(RuntimeError):
    """Raised when the public pooling requirements or evidence is inconsistent."""


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def validate_contract(
    contract: dict[str, object],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    for key in (
        "automatic_activation",
        "sealed_outcomes_read",
        "biological_matrices_downloaded",
        "controlled_data_accessed",
        "pooled_libraries_are_donors",
        "cells_are_biological_replicates",
        "normal_stratum_is_masld_negative_control",
        "donor_safe_split_available",
        "external_masld_champion_eligible",
    ):
        if contract.get(key) is not False:
            raise PoolingTopologyError(f"{key} must remain false")
    expected_counts = {
        "reported_person_n": 15,
        "condition_n": 3,
        "persons_per_condition": 5,
        "physical_multiome_library_n": 6,
        "libraries_per_condition": 2,
        "geo_assay_record_n": 12,
        "assay_records_per_physical_library": 2,
    }
    for key, expected in expected_counts.items():
        if contract.get(key) != expected:
            raise PoolingTopologyError(f"{key} must equal {expected}")
    if sorted(contract.get("pool_size_multiset_per_condition", [])) != [2, 3]:
        raise PoolingTopologyError("each condition must retain one 2-person and one 3-person pool")
    if contract.get("numbered_pool_size_assignment") != "UNRESOLVED_PUBLIC_METADATA":
        raise PoolingTopologyError("numbered pool sizes cannot be inferred")
    if contract.get("individual_to_pool_join") != "UNAVAILABLE_PUBLIC_METADATA":
        raise PoolingTopologyError("individual-to-pool join cannot be asserted")
    if contract.get("nucleus_to_individual_join") != "UNAVAILABLE_PUBLIC_METADATA":
        raise PoolingTopologyError("nucleus-to-person join cannot be asserted")
    conditions = contract.get("condition")
    if not isinstance(conditions, list) or {row.get("condition_id") for row in conditions} != REQUIRED_CONDITIONS:
        raise PoolingTopologyError("condition census must be Normal, SAH, and AC")
    if any(row.get("person_n") != 5 or row.get("physical_library_n") != 2 for row in conditions):
        raise PoolingTopologyError("every condition must have five people in two pools")
    if any(row.get("masld_negative_control") is not False for row in conditions):
        raise PoolingTopologyError("no alcohol-study stratum is a MASLD-negative control")
    pairs = contract.get("assay_record_pair")
    if not isinstance(pairs, list) or len(pairs) != 6:
        raise PoolingTopologyError("six physical-library assay pairs are required")
    library_ids = [row.get("library_id") for row in pairs]
    if len(library_ids) != len(set(library_ids)):
        raise PoolingTopologyError("physical library IDs must be unique")
    gsm_ids = [str(row[key]) for row in pairs for key in ("atac_gsm", "rna_gsm")]
    if len(gsm_ids) != 12 or len(gsm_ids) != len(set(gsm_ids)):
        raise PoolingTopologyError("12 unique GEO assay records are required")
    if any(row.get("person_n") != -1 for row in pairs):
        raise PoolingTopologyError("numbered pool person counts must remain unresolved")
    tasks = contract.get("task_disposition")
    if not isinstance(tasks, list):
        raise PoolingTopologyError("task disposition census is missing")
    allowed = {row["task_id"] for row in tasks if row.get("allowed") is True}
    prohibited = {row["task_id"] for row in tasks if row.get("allowed") is False}
    if allowed != REQUIRED_ALLOWED_TASKS or prohibited != REQUIRED_PROHIBITED_TASKS:
        raise PoolingTopologyError("task disposition census changed")
    if any(row.get("model_selection_allowed") is not False for row in tasks):
        raise PoolingTopologyError("no GSE281574 task may select models")
    decision = contract.get("decision")
    if not isinstance(decision, dict):
        raise PoolingTopologyError("decision block is missing")
    if decision.get("current_role") != "alcohol_etiology_ood_and_descriptive_only":
        raise PoolingTopologyError("GSE281574 role must remain ALD OOD/descriptive")
    if decision.get("normal_stratum_may_be_used_as_masld_negative_controls") is not False:
        raise PoolingTopologyError("Normal pools cannot become MASLD-negative controls")
    if "gse289173" in json.dumps(contract).lower():
        raise PoolingTopologyError("sealed accession must not enter this audit")
    sources = contract.get("public_source")
    if not isinstance(sources, list) or len(sources) != 3:
        raise PoolingTopologyError("three public evidence sources are required")
    return pairs, tasks


def parse_geo_soft(payload: bytes) -> list[dict[str, object]]:
    try:
        text = gzip.decompress(payload).decode("utf-8")
    except (OSError, EOFError, UnicodeDecodeError) as error:
        raise PoolingTopologyError("GEO SOFT is not valid gzip text") from error
    records: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE = "):
            if current is not None:
                records.append(current)
            current = {
                "gsm": line.split(" = ", 1)[1],
                "relations": [],
                "characteristics": [],
            }
        elif current is not None and line.startswith("!Sample_title = "):
            current["title"] = line.split(" = ", 1)[1]
        elif current is not None and line.startswith("!Sample_relation = "):
            current["relations"].append(line.split(" = ", 1)[1])
        elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
            current["characteristics"].append(line.split(" = ", 1)[1])
    if current is not None:
        records.append(current)
    if len(records) != 12:
        raise PoolingTopologyError(f"expected 12 GEO assay records, found {len(records)}")
    for record in records:
        match = TITLE_PATTERN.fullmatch(str(record.get("title", "")))
        if match is None:
            raise PoolingTopologyError(f"unexpected GEO sample title: {record.get('title')}")
        condition, pool_number, modality = match.groups()
        record["condition_id"] = condition
        record["pool_number"] = int(pool_number)
        record["modality"] = modality
        record["library_id"] = f"{condition}_pool{pool_number}"
        relations = list(record["relations"])
        biosamples = [value.rsplit("/", 1)[-1] for value in relations if "biosample/" in value.lower()]
        srx = [value.split("term=", 1)[-1] for value in relations if "term=SRX" in value]
        if len(biosamples) != 1 or len(srx) != 1:
            raise PoolingTopologyError(f"BioSample/SRX relation unresolved: {record['gsm']}")
        record["biosample"] = biosamples[0]
        record["srx"] = srx[0]
    return records


def build_library_pairs(
    records: list[dict[str, object]], expected_pairs: list[dict[str, object]]
) -> list[dict[str, object]]:
    by_library: dict[str, dict[str, dict[str, object]]] = {}
    for record in records:
        library_id = str(record["library_id"])
        modality = str(record["modality"])
        if modality in by_library.setdefault(library_id, {}):
            raise PoolingTopologyError(f"duplicate modality record: {library_id} {modality}")
        by_library[library_id][modality] = record
    expected_by_library = {str(row["library_id"]): row for row in expected_pairs}
    if set(by_library) != set(expected_by_library):
        raise PoolingTopologyError("GEO physical-library roster differs from contract")
    output: list[dict[str, object]] = []
    for library_id in sorted(by_library):
        assays = by_library[library_id]
        if set(assays) != {"ATAC", "RNA"}:
            raise PoolingTopologyError(f"RNA-ATAC pair incomplete: {library_id}")
        expected = expected_by_library[library_id]
        if assays["ATAC"]["gsm"] != expected["atac_gsm"] or assays["RNA"]["gsm"] != expected["rna_gsm"]:
            raise PoolingTopologyError(f"GEO assay-pair identity changed: {library_id}")
        output.append(
            {
                "library_id": library_id,
                "condition_id": assays["ATAC"]["condition_id"],
                "pool_number": assays["ATAC"]["pool_number"],
                "reported_pool_person_n": None,
                "reported_pool_person_n_state": "one_of_2_or_3_unresolved_numbered_assignment",
                "same_nucleus_rna_atac": True,
                "person_identity_available": False,
                "nucleus_to_person_join_available": False,
                "atac": {
                    "gsm": assays["ATAC"]["gsm"],
                    "biosample": assays["ATAC"]["biosample"],
                    "srx": assays["ATAC"]["srx"],
                },
                "rna": {
                    "gsm": assays["RNA"]["gsm"],
                    "biosample": assays["RNA"]["biosample"],
                    "srx": assays["RNA"]["srx"],
                },
                "deposition_interpretation": "two_assay_records_one_physical_multiome_library_not_two_donors",
            }
        )
    return output


def task_decision(tasks: list[dict[str, object]]) -> dict[str, object]:
    allowed = [row for row in tasks if row["allowed"] is True]
    prohibited = [row for row in tasks if row["allowed"] is False]
    return {
        "schema_version": "masld-bench-gse281574-task-decision-v1",
        "current_role": "alcohol_etiology_ood_and_descriptive_only",
        "reported_person_n": 15,
        "physical_multiome_library_n": 6,
        "geo_assay_record_n": 12,
        "minimum_independent_molecular_unit": "physical_multiome_library_pool",
        "donor_safe_split_available": False,
        "library_grouped_technical_split_available": True,
        "condition_level_independent_unit_n": 2,
        "nested_model_selection_supported": False,
        "model_selection_allowed": False,
        "external_masld_champion_eligible": False,
        "normal_stratum_is_masld_negative_control": False,
        "allowed_descriptive_tasks": allowed,
        "prohibited_tasks": prohibited,
        "cells_are_biological_replicates": False,
        "biological_matrices_downloaded": False,
        "sealed_outcomes_read": False,
    }


def fetch_source(source_id: str, url: str) -> tuple[bytes, dict[str, object]]:
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme != "https":
        raise PoolingTopologyError(f"non-HTTPS public source rejected: {url}")
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json,application/xml,text/plain,text/html,application/gzip,*/*",
            "User-Agent": "masld-bench-gse281574-public-metadata-audit/1.0",
        },
    )
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read(MAX_SOURCE_BYTES + 1)
                if len(payload) > MAX_SOURCE_BYTES:
                    raise PoolingTopologyError(f"source exceeds public metadata cap: {source_id}")
                return payload, {
                    "source_id": source_id,
                    "url": response.geturl(),
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
    raise PoolingTopologyError(f"failed to fetch {source_id}: {last_error}") from last_error


def decode_source(payload: bytes, encoding: str) -> str:
    if encoding == "gzip_text":
        try:
            return gzip.decompress(payload).decode("utf-8")
        except (OSError, EOFError, UnicodeDecodeError) as error:
            raise PoolingTopologyError("invalid public gzip text") from error
    return payload.decode("utf-8", errors="replace")


def safe_filename(source_id: str, encoding: str) -> str:
    suffix = ".soft.gz" if encoding == "gzip_text" else ".txt"
    return re.sub(r"[^a-z0-9_.-]+", "_", source_id.lower()) + suffix


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
        raise PoolingTopologyError(f"refusing to overwrite output: {output}")
    contract = tomllib.loads(contract_path.read_text(encoding="utf-8"))
    expected_pairs, tasks = validate_contract(contract)
    metadata_dir = output / "public_metadata"
    metadata_dir.mkdir(parents=True)

    receipts: list[dict[str, object]] = []
    payloads: dict[str, bytes] = {}
    for spec in contract["public_source"]:
        source_id = str(spec["source_id"])
        encoding = str(spec["encoding"])
        payload, receipt = fetch_source(source_id, str(spec["url"]))
        text = decode_source(payload, encoding).lower()
        missing = [token for token in spec["required_tokens"] if str(token).lower() not in text]
        if missing:
            raise PoolingTopologyError(f"{source_id} missing required public evidence: {missing}")
        target = metadata_dir / safe_filename(source_id, encoding)
        target.write_bytes(payload)
        receipt["local_path"] = str(target.relative_to(output))
        receipt["required_evidence_pass"] = True
        receipts.append(receipt)
        payloads[source_id] = payload

    records = parse_geo_soft(payloads["geo_soft"])
    library_pairs = build_library_pairs(records, expected_pairs)
    topology = {
        "schema_version": "masld-bench-gse281574-pooling-topology-audit-v1",
        "reported_person_n": 15,
        "condition_person_counts": {"Normal": 5, "SAH": 5, "AC": 5},
        "physical_multiome_library_n": 6,
        "geo_assay_record_n": 12,
        "pool_size_multiset_per_condition": [2, 3],
        "numbered_pool_size_assignment": "UNRESOLVED_PUBLIC_METADATA",
        "individual_to_pool_join": "UNAVAILABLE_PUBLIC_METADATA",
        "nucleus_to_individual_join": "UNAVAILABLE_PUBLIC_METADATA",
        "same_nucleus_pairing_scope": "RNA_and_ATAC_within_each_physical_multiome_library",
        "cross_library_or_donor_pairing": "prohibited",
        "library_pairs": library_pairs,
    }
    decision = task_decision(tasks)
    rights_reference = {
        "schema_version": "masld-bench-gse281574-rights-reference-v1",
        "source_high_level_build": contract["reference"]["source_high_level_build"],
        "exact_reference_ready": contract["reference"]["sequence_or_coordinate_task_ready"],
        "exact_cell_ranger_arc_version": contract["reference"]["exact_cell_ranger_arc_version"],
        "geo_download_is_public": contract["rights"]["geo_download_is_public"],
        "processed_object_redistribution_terms": contract["rights"]["processed_object_redistribution_terms"],
        "derivative_weight_redistribution_terms": contract["rights"]["derivative_weight_redistribution_terms"],
        "biological_matrices_downloaded": False,
    }
    (output / "source_receipts.json").write_bytes(canonical_json({"sources": receipts}))
    (output / "pooling_topology.json").write_bytes(canonical_json(topology))
    (output / "task_eligibility.json").write_bytes(canonical_json(decision))
    (output / "rights_reference_readiness.json").write_bytes(canonical_json(rights_reference))
    (output / "DECISION.md").write_text(
        "# GSE281574 pooling and task decision\n\n"
        "**Decision: ALD OOD/descriptive only.** Fifteen people contributed to six "
        "physical same-nucleus RNA-ATAC libraries, with two pools per condition and "
        "pool sizes `{3,2}`. GEO exposes 12 assay records because RNA and ATAC are "
        "deposited separately; those records are not donors.\n\n"
        "Public metadata does not map numbered pools or nuclei to individuals. A "
        "six-library grouped technical split is possible, but a donor-safe split is "
        "not. With only two libraries per condition, this source cannot select models, "
        "support confirmatory inference, serve as MASLD-negative controls, or validate "
        "a MASLD champion. No biological matrix or sealed outcome was accessed.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
