#!/usr/bin/env python3
"""Freeze an outcome-blind public-metadata audit of liver RNA-ATAC cohorts."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import tomllib
import urllib.error
import urllib.request


MAX_SOURCE_BYTES = 32 * 1024 * 1024
VALID_TOPOLOGIES = {
    "same_nucleus",
    "same_nucleus_donor_join_unresolved",
    "same_sample_different_aliquot",
}
REQUIRED_CATEGORIES = {
    "masld_same_nucleus",
    "alcohol_ood_same_nucleus",
    "masld_same_sample_unpaired",
    "healthy_reference_same_nucleus",
    "masld_same_donor_unpaired",
    "controlled_pediatric_liver_ood_same_nucleus",
}

REMOTE_SOURCES = {
    "gse281574_geo_soft": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281574/soft/GSE281574_family.soft.gz",
        "encoding": "gzip_text",
        "required": ["^SERIES = GSE281574", "Normal_pool1_ATAC", "AC_pool2_RNA", "Assembly: hg38"],
    },
    "gse286690_geo_soft": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286690/soft/GSE286690_family.soft.gz",
        "encoding": "gzip_text",
        "required": ["^SERIES = GSE286690", "ENCSR367YUX", "ENCSR630YEA"],
    },
    "gse244832_geo_soft": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244832/soft/GSE244832_family.soft.gz",
        "encoding": "gzip_text",
        "required": ["^SERIES = GSE244832", "MASH, MASL, NORMAL", "18 single nucleus ATACseq samples"],
    },
    "encode_rna_experiment": {
        "url": "https://www.encodeproject.org/experiments/ENCSR367YUX/?format=json&frame=object",
        "encoding": "text",
        "required": ["ENCSR367YUX", "ENCSR630YEA", "ENCSR788DZO", "released"],
    },
    "encode_atac_experiment": {
        "url": "https://www.encodeproject.org/experiments/ENCSR630YEA/?format=json&frame=object",
        "encoding": "text",
        "required": ["ENCSR630YEA", "ENCSR367YUX", "ENCSR788DZO", "released"],
    },
    "encode_data_use_policy": {
        "url": "https://www.encodeproject.org/about/data-use-policy/",
        "encoding": "text",
        "required": ["unrestricted use", "immediately upon release"],
    },
    "fnih_pubmed": {
        "url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=40385416&retmode=xml",
        "encoding": "text",
        "required": ["Single cell multiomics reveals drivers", "Preprint"],
    },
    "fnih_zenodo": {
        "url": "https://zenodo.org/api/records/15298484",
        "encoding": "text",
        "required": ["15298484"],
    },
    "fnih_processed_portal": {
        "url": "https://epigenome.wustl.edu/MASLD/",
        "encoding": "text",
        "required": ["MASLD"],
        "allow_unavailable": True,
    },
    "emtab13130_biostudies": {
        "url": "https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-13130",
        "encoding": "text",
        "required": ["E-MTAB-13130"],
    },
    "emtab13131_biostudies": {
        "url": "https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-13131",
        "encoding": "text",
        "required": ["E-MTAB-13131"],
    },
    "hepatoblastoma_ega_multiome": {
        "url": "https://metadata.ega-archive.org/datasets/EGAD00001010049",
        "encoding": "text",
        "required": ["EGAD00001010049", "controlled"],
    },
    "hepatoblastoma_ega_study": {
        "url": "https://metadata.ega-archive.org/studies/EGAS00001006932",
        "encoding": "text",
        "required": ["EGAS00001006932"],
    },
    "hepatoblastoma_ega_multiome_study": {
        "url": "https://metadata.ega-archive.org/datasets/EGAD00001010049/studies",
        "encoding": "text",
        "required": ["EGAS00001006932"],
    },
}


class CohortAuditError(RuntimeError):
    """Raised when metadata or the frozen cohort requirements are invalid."""


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def validate_contract(contract: dict[str, object]) -> list[dict[str, object]]:
    if contract.get("sealed_outcomes_read") is not False:
        raise CohortAuditError("sealed_outcomes_read must remain false")
    if contract.get("biological_matrices_downloaded") is not False:
        raise CohortAuditError("biological_matrices_downloaded must remain false")
    if contract.get("candidate_activation_is_automatic") is not False:
        raise CohortAuditError("candidate activation cannot be automatic")
    if contract.get("anchor_excluded_from_ranking") != "gse296875":
        raise CohortAuditError("GSE296875 must remain the excluded topology anchor")
    candidates = contract.get("candidate")
    if not isinstance(candidates, list) or not candidates:
        raise CohortAuditError("candidate roster is empty")
    ranks = [row.get("rank") for row in candidates]
    if ranks != list(range(1, len(candidates) + 1)):
        raise CohortAuditError("candidate ranks must be consecutive and ordered")
    ids = [row.get("candidate_id") for row in candidates]
    if len(ids) != len(set(ids)):
        raise CohortAuditError("candidate IDs must be unique")
    categories = {str(row.get("category")) for row in candidates}
    if categories != REQUIRED_CATEGORIES:
        raise CohortAuditError(f"category census mismatch: {sorted(categories)}")
    for row in candidates:
        if not isinstance(row.get("exact_donor_n"), int) or int(row["exact_donor_n"]) < 1:
            raise CohortAuditError(f"invalid donor n: {row.get('candidate_id')}")
        if row.get("topology") not in VALID_TOPOLOGIES:
            raise CohortAuditError(f"invalid topology: {row.get('candidate_id')}")
        if row.get("eligible_for_masld_external_champion") is not False:
            raise CohortAuditError("no candidate is externally champion-eligible")
        for required_key in (
            "accessions",
            "phenotype_fields",
            "native_reference_build",
            "permissions",
            "contamination_role",
            "next_acquisition_job",
            "source_urls",
        ):
            if not row.get(required_key):
                raise CohortAuditError(
                    f"{row.get('candidate_id')} missing {required_key}"
                )
    decision = contract.get("decision")
    if not isinstance(decision, dict):
        raise CohortAuditError("decision block missing")
    if decision.get("paired_masld_external_seal_available_now") is not False:
        raise CohortAuditError("paired MASLD external seal must remain unavailable")
    if decision.get("observed_multiome_external_champion_available_now") is not False:
        raise CohortAuditError("observed multiome champion must remain unavailable")
    if "gse289173" in json.dumps(contract).lower():
        raise CohortAuditError("sealed accession must not enter this audit")
    return candidates


def fetch_source(source_id: str, spec: dict[str, object]) -> tuple[bytes | None, dict[str, object]]:
    url = str(spec["url"])
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json,text/plain,text/html,application/xml,*/*",
            "User-Agent": "masld-bench-outcome-blind-metadata-audit/1.0",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = response.read(MAX_SOURCE_BYTES + 1)
            if len(payload) > MAX_SOURCE_BYTES:
                raise CohortAuditError(f"source exceeds metadata cap: {source_id}")
            receipt: dict[str, object] = {
                "source_id": source_id,
                "url": url,
                "http_status": response.status,
                "content_type": response.headers.get("Content-Type", ""),
                "last_modified": response.headers.get("Last-Modified", ""),
                "etag": response.headers.get("ETag", ""),
                "bytes": len(payload),
                "sha256": sha256_bytes(payload),
                "available": True,
            }
    except (urllib.error.URLError, TimeoutError) as error:
        if spec.get("allow_unavailable") is True:
            return None, {
                "source_id": source_id,
                "url": url,
                "available": False,
                "error_class": type(error).__name__,
                "error": str(error),
            }
        raise CohortAuditError(f"failed to fetch {source_id}: {error}") from error
    if spec["encoding"] == "gzip_text":
        try:
            text = gzip.decompress(payload).decode("utf-8")
        except (OSError, EOFError, UnicodeDecodeError) as error:
            raise CohortAuditError(f"invalid gzip text: {source_id}") from error
    else:
        text = payload.decode("utf-8", errors="replace")
    missing = [pattern for pattern in spec["required"] if pattern.lower() not in text.lower()]
    if missing:
        raise CohortAuditError(f"{source_id} missing required evidence: {missing}")
    receipt["required_evidence_pass"] = True
    receipt["required_patterns"] = list(spec["required"])
    if source_id == "fnih_zenodo":
        record = json.loads(text)
        metadata = record.get("metadata", {})
        files = record.get("files", [])
        receipt["access_right"] = metadata.get("access_right")
        receipt["license_id"] = (metadata.get("license") or {}).get("id")
        receipt["file_count"] = len(files)
        receipt["publication_date"] = metadata.get("publication_date")
        if receipt["access_right"] != "restricted" or receipt["file_count"] != 0:
            raise CohortAuditError(
                "FNIH Zenodo access state changed; revise the activation contract"
            )
    if source_id == "hepatoblastoma_ega_multiome":
        record = json.loads(text)
        receipt["access_type"] = record.get("access_type")
        if receipt["access_type"] != "controlled":
            raise CohortAuditError(
                "hepatoblastoma multiome access state changed; revise the contract"
            )
    return payload, receipt


def safe_source_name(source_id: str, spec: dict[str, object]) -> str:
    suffix = ".soft.gz" if spec["encoding"] == "gzip_text" else ".txt"
    return re.sub(r"[^a-z0-9_.-]+", "_", source_id.lower()) + suffix


def write_tsv(path: Path, candidates: list[dict[str, object]]) -> None:
    fields = [
        "rank",
        "candidate_id",
        "category",
        "accessions",
        "exact_donor_n",
        "topology",
        "phenotype_fields",
        "native_reference_build",
        "permissions",
        "availability",
        "contamination_role",
        "activation_status",
        "eligible_for_masld_external_champion",
        "next_acquisition_job",
    ]
    with path.open("x", encoding="utf-8") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in candidates:
            values = []
            for key in fields:
                value = row[key]
                if isinstance(value, list):
                    value = "|".join(str(item) for item in value)
                values.append(str(value).replace("\t", " ").replace("\n", " "))
            handle.write("\t".join(values) + "\n")


def build_summary(candidates: list[dict[str, object]], receipts: list[dict[str, object]]) -> dict[str, object]:
    same_nucleus = [row for row in candidates if str(row["topology"]).startswith("same_nucleus")]
    receipt_by_id = {str(row["source_id"]): row for row in receipts}
    fnih_zenodo = receipt_by_id.get("fnih_zenodo", {})
    return {
        "schema_version": "masld-bench-external-observed-multiome-audit-summary-v1",
        "outcome_blind": True,
        "sealed_outcomes_read": False,
        "biological_matrices_downloaded": False,
        "gse296875_used_only_as_excluded_anchor": True,
        "candidate_count": len(candidates),
        "same_nucleus_candidate_count": len(same_nucleus),
        "masld_same_nucleus_candidate_count": sum(
            row["category"] == "masld_same_nucleus" for row in candidates
        ),
        "healthy_reference_candidate_count": sum(
            row["category"] == "healthy_reference_same_nucleus" for row in candidates
        ),
        "same_donor_or_sample_unpaired_candidate_count": sum(
            row["topology"] == "same_sample_different_aliquot" for row in candidates
        ),
        "public_metadata_sources_requested": len(receipts),
        "public_metadata_sources_available": sum(bool(row["available"]) for row in receipts),
        "fnih_zenodo_access_right": fnih_zenodo.get("access_right", "not_checked"),
        "fnih_zenodo_file_count": fnih_zenodo.get("file_count", "not_checked"),
        "paired_masld_external_seal_available_now": False,
        "observed_multiome_external_champion_available_now": False,
        "recommended_next_job": "fnih_86_portal_inventory_rights_and_donor_join_v1",
        "ranked_candidate_ids": [row["candidate_id"] for row in candidates],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    if output.exists():
        raise CohortAuditError(f"refusing to overwrite output: {output}")
    contract_path = args.contract.resolve(strict=True)
    contract = tomllib.loads(contract_path.read_text(encoding="utf-8"))
    candidates = validate_contract(contract)
    (output / "public_metadata").mkdir(parents=True)
    receipts: list[dict[str, object]] = []
    for source_id, spec in REMOTE_SOURCES.items():
        payload, receipt = fetch_source(source_id, spec)
        receipts.append(receipt)
        if payload is not None:
            target = output / "public_metadata" / safe_source_name(source_id, spec)
            with target.open("xb") as handle:
                handle.write(payload)
            receipt["local_path"] = str(target.relative_to(output))
    with (output / "ranked_candidates.json").open("xb") as handle:
        handle.write(canonical_json({"candidates": candidates}))
    write_tsv(output / "ranked_candidates.tsv", candidates)
    with (output / "source_receipts.json").open("xb") as handle:
        handle.write(canonical_json({"sources": receipts}))
    summary = build_summary(candidates, receipts)
    with (output / "audit_summary.json").open("xb") as handle:
        handle.write(canonical_json(summary))
    decision = contract["decision"]
    (output / "DECISION.md").write_text(
        "# External observed-multiome cohort expansion audit\n\n"
        "No topology-matched external MASLD champion source is active. "
        "FNIH 86-liver is the first activation target, but remains blocked on an "
        "open donor-resolved release, exact reference identity, joins, and reuse terms.\n\n"
        f"Next metadata acquisition job: `{decision['immediate_metadata_job']}`.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
