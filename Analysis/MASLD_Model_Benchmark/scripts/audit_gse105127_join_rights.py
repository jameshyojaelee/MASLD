#!/usr/bin/env python3
"""Freeze the authoritative GSE105127 donor/zone/assay join and reuse requirements."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re
import urllib.request
import xml.etree.ElementTree as ET


ARTICLE_URL = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6175862/fullTextXML"
)
ARTICLE_DOI = "10.1038/s41467-018-06611-5"
EXPECTED_METADATA_ARTIFACTS_SHA256 = (
    "bac433336705b74dcb2a6ef6b6594af10fef362769d564adfea9788fe546e88e"
)
MAX_ARTICLE_BYTES = 16 * 1024 * 1024
TITLE_PATTERN = re.compile(
    r"^(?P<participant>[0-9]+)_(?P<zone>CV|IZ|PP)_(?P<assay>RNA|RRBS)$"
)
GROUP_PATTERN = re.compile(
    r"(?P<label>NC|HO|STEA|EARLY)\s*=\s*(?P<ids>[0-9, ]+)"
)
EXPECTED_GROUP_COUNTS = {"NC": 4, "HO": 5, "STEA": 5, "EARLY": 5}
GROUP_NAMES = {
    "NC": "normal_control",
    "HO": "healthy_obese",
    "STEA": "bland_steatosis",
    "EARLY": "early_nash",
}
ZONE_CHARACTERISTICS = {
    "CV": "hepatic zone: pericentral (CV)",
    "IZ": "hepatic zone: intermediate (IZ)",
    "PP": "hepatic zone: periportal (PP)",
}
ADJACENT_SECTION_EVIDENCE = (
    "Adjacent cryosections were used to generate matching reduced "
    "representation bisulfite sequencing"
)


class GSE105127AuditError(RuntimeError):
    """Raised when frozen source evidence does not meet the inclusion requirements."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_series_fields(payload: bytes) -> dict[str, list[str]]:
    try:
        text = gzip.decompress(payload).decode("utf-8", errors="strict")
    except (OSError, EOFError, UnicodeDecodeError) as error:
        raise GSE105127AuditError("GEO SOFT is not strict UTF-8 gzip text") from error
    fields: dict[str, list[str]] = {}
    for line in text.splitlines():
        if line.startswith("!Series_") and " = " in line:
            key, value = line.split(" = ", 1)
            fields.setdefault(key, []).append(value.strip())
    return fields


def parse_group_mapping(overall_design: str) -> dict[str, str]:
    by_participant: dict[str, str] = {}
    observed_counts: Counter[str] = Counter()
    for match in GROUP_PATTERN.finditer(overall_design):
        label = match["label"]
        for participant in (value.strip() for value in match["ids"].split(",")):
            if not participant or not participant.isdigit():
                raise GSE105127AuditError("invalid participant in overall design")
            if participant in by_participant:
                raise GSE105127AuditError("participant appears in multiple phenotype groups")
            by_participant[participant] = label
            observed_counts[label] += 1
    if dict(observed_counts) != EXPECTED_GROUP_COUNTS or len(by_participant) != 19:
        raise GSE105127AuditError("overall-design phenotype groups differ")
    return by_participant


def read_sample_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 114:
        raise GSE105127AuditError("sample table no longer contains 114 records")
    if len({row["accession"] for row in rows}) != 114:
        raise GSE105127AuditError("sample accessions are not unique")
    return rows


def validate_supplementary_files(assay: str, files: list[str]) -> str:
    if assay == "RRBS":
        endings = (".CG.bed.gz", ".CG.bw", ".CG.ct_coverage.bw")
        if len(files) != 3 or not all(
            any(value.endswith(ending) for value in files) for ending in endings
        ):
            raise GSE105127AuditError("RRBS processed-file set differs")
        return next(value for value in files if value.endswith(".CG.bed.gz"))
    if len(files) != 1 or not files[0].endswith("_RNA.coverage.bw"):
        raise GSE105127AuditError("RNA processed-file set differs")
    return files[0]


def fetch_article(url: str = ARTICLE_URL) -> tuple[bytes, dict[str, str]]:
    request = urllib.request.Request(
        url, headers={"User-Agent": "masld-bench-source-audit/1.0"}
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read(MAX_ARTICLE_BYTES + 1)
        if len(payload) > MAX_ARTICLE_BYTES:
            raise GSE105127AuditError("article XML exceeds size limit")
        declared = response.headers.get("Content-Length")
        if declared is not None and int(declared) != len(payload):
            raise GSE105127AuditError("article XML Content-Length differs")
        receipt = {
            key.lower(): value
            for key in ("Content-Length", "Last-Modified", "ETag")
            if (value := response.headers.get(key)) is not None
        }
    return payload, receipt


def normalize_text(value: str) -> str:
    return " ".join(value.split())


def audit_article(payload: bytes) -> dict[str, object]:
    try:
        root = ET.fromstring(payload)
    except ET.ParseError as error:
        raise GSE105127AuditError("article source is not valid XML") from error
    article_text = normalize_text(" ".join(root.itertext()))
    doi_values = {
        normalize_text("".join(node.itertext())).lower()
        for node in root.findall(".//article-id")
        if node.attrib.get("pub-id-type") == "doi"
    }
    if ARTICLE_DOI not in doi_values:
        raise GSE105127AuditError("article DOI differs")
    if ADJACENT_SECTION_EVIDENCE not in article_text:
        raise GSE105127AuditError("adjacent-section evidence is absent")
    license_links = sorted(
        {
            value
            for node in root.findall(".//license-ref")
            for key, value in node.attrib.items()
            if key.endswith("href")
        }
    )
    return {
        "doi": ARTICLE_DOI,
        "pmcid": "PMC6175862",
        "adjacent_section_evidence_present": True,
        "article_license_links": license_links,
        "article_license_does_not_assign_data_or_weight_rights": True,
    }


def fold_mapping(group_mapping: dict[str, str]) -> dict[str, int]:
    folds: dict[str, int] = {}
    for label in ("NC", "HO", "STEA", "EARLY"):
        participants = sorted(
            participant
            for participant, observed_label in group_mapping.items()
            if observed_label == label
        )
        for index, participant in enumerate(participants):
            folds[participant] = index
    return folds


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument(
        "--metadata-artifacts-sha256",
        default=EXPECTED_METADATA_ARTIFACTS_SHA256,
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--article-url", default=ARTICLE_URL)
    args = parser.parse_args()
    if sha256_file(args.metadata / "ARTIFACTS.json") != args.metadata_artifacts_sha256:
        raise GSE105127AuditError("metadata ARTIFACTS SHA-256 differs")

    soft_path = args.metadata / "raw" / "GSE105127_family.soft.gz"
    sample_path = args.metadata / "samples" / "GSE105127.tsv"
    fields = parse_series_fields(soft_path.read_bytes())
    overall_designs = fields.get("!Series_overall_design", [])
    if len(overall_designs) != 1:
        raise GSE105127AuditError("expected exactly one overall-design field")
    if fields.get("!Series_pubmed_id") != ["30297808"]:
        raise GSE105127AuditError("GEO PubMed identifier differs")
    group_mapping = parse_group_mapping(overall_designs[0])
    rows = read_sample_rows(sample_path)
    article_payload, article_receipt = fetch_article(args.article_url)
    article = audit_article(article_payload)
    folds = fold_mapping(group_mapping)

    sample_output: list[dict[str, object]] = []
    combinations: Counter[tuple[str, str, str]] = Counter()
    for row in rows:
        match = TITLE_PATTERN.fullmatch(row["title"])
        if match is None:
            raise GSE105127AuditError("sample title violates donor-zone-assay contract")
        participant = match["participant"]
        zone = match["zone"]
        assay = match["assay"]
        if row["candidate_participant_id"] != participant:
            raise GSE105127AuditError("sample table participant join differs")
        if participant not in group_mapping:
            raise GSE105127AuditError("sample participant absent from overall design")
        characteristics = json.loads(row["characteristics_json"])
        if ZONE_CHARACTERISTICS[zone] not in characteristics:
            raise GSE105127AuditError("sample zone characteristic differs")
        if "isolation: laser capture microdissection" not in characteristics:
            raise GSE105127AuditError("sample isolation characteristic differs")
        files = json.loads(row["supplementary_files_json"])
        primary_processed_file = validate_supplementary_files(assay, files)
        combinations[(participant, zone, assay)] += 1
        sample_output.append(
            {
                "participant_id": participant,
                "phenotype": group_mapping[participant],
                "phenotype_name": GROUP_NAMES[group_mapping[participant]],
                "outer_fold": folds[participant],
                "sample_accession": row["accession"],
                "zone": zone,
                "assay": assay,
                "pairing": "adjacent_section",
                "modality_status": "observed",
                "native_genome_build": "1000_Genomes_GRCh37_hg19",
                "primary_processed_file": primary_processed_file,
            }
        )
    expected = {
        (participant, zone, assay)
        for participant in group_mapping
        for zone in ("CV", "IZ", "PP")
        for assay in ("RNA", "RRBS")
    }
    if set(combinations) != expected or any(value != 1 for value in combinations.values()):
        raise GSE105127AuditError("donor-zone-assay cross-product differs")

    args.output.mkdir(parents=True, exist_ok=False)
    raw = args.output / "raw"
    raw.mkdir()
    article_path = raw / "PMC6175862.fullTextXML"
    article_path.write_bytes(article_payload)
    participant_output = [
        {
            "participant_id": participant,
            "phenotype": group_mapping[participant],
            "phenotype_name": GROUP_NAMES[group_mapping[participant]],
            "outer_fold": folds[participant],
            "biological_unit": "participant",
            "assay_records": 6,
        }
        for participant in sorted(group_mapping)
    ]
    for name, output_rows in (
        ("participant_join.tsv", participant_output),
        ("sample_join.tsv", sample_output),
    ):
        with (args.output / name).open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=tuple(output_rows[0]),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(output_rows)

    join = {
        "schema_version": "masld-bench-gse105127-authoritative-join-v1",
        "status": "pass",
        "biological_unit": "participant",
        "participants": 19,
        "assay_records": 114,
        "phenotype_counts": dict(sorted(Counter(group_mapping.values()).items())),
        "outer_fold_participant_counts": dict(
            sorted(Counter(folds.values()).items())
        ),
        "zones": ["CV", "IZ", "PP"],
        "assays": ["RNA", "RRBS"],
        "pairing": "adjacent_section",
        "same_cell_or_same_section_pairing": False,
        "cells_sections_libraries_as_independent_replicates": False,
        "source_evidence": {
            "geo_overall_design": True,
            "geo_sample_titles_and_characteristics": True,
            "article": article,
        },
        "model_training_activated": False,
        "remaining_blockers": [
            "download_and_assay_native_qc_of_all_57_RRBS_BEDs",
            "RRBS_coordinate_and_coverage_semantics",
            "failure_aware_hg19_to_GRCh38p14_CpG_crosswalk",
            "RNA_count_or_read_level_processing_contract",
            "task_specific_training_activation_contract",
        ],
    }
    rights = {
        "schema_version": "masld-bench-gse105127-rights-v1",
        "access_tier": "public",
        "new_dua_or_controlled_access_required": False,
        "internal_nonclinical_research_use": "allowed_by_project_public_data_definition",
        "source_terms_apply": True,
        "explicit_data_license_detected": False,
        "article_license_links": article["article_license_links"],
        "article_license_does_not_assign_data_or_weight_rights": True,
        "raw_or_processed_data_redistribution": "prohibited_by_project_policy",
        "released_weight_redistribution": "requires_model_specific_terms_and_legal_review",
        "clinical_use": False,
        "open_champion_eligibility": "blocked_until_derivative_weight_review_if_this_source_materially_trains_the_champion",
    }
    receipt = {
        "article_url": args.article_url,
        "article_sha256": sha256(article_payload).hexdigest(),
        "article_bytes": len(article_payload),
        "http_headers": article_receipt,
        "metadata_artifacts_sha256": args.metadata_artifacts_sha256,
        "soft_sha256": sha256_file(soft_path),
        "sample_table_sha256": sha256_file(sample_path),
    }
    for name, value in (
        ("join_audit.json", join),
        ("rights_contract.json", rights),
        ("source_receipt.json", receipt),
    ):
        (args.output / name).write_text(
            json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(json.dumps({"join": join, "rights": rights}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
