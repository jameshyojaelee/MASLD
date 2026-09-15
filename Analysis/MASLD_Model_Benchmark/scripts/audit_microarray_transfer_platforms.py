#!/usr/bin/env python3
"""Audit platform axes, CEL formats, rights, and cohort-family overlap."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import html
import json
from pathlib import Path
import re
import urllib.request

from scripts.acquire_audit_microarray_transfer_raw import identify_cel_format


RIGHTS_URL = "https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html"
RIGHTS_MAX_BYTES = 2 * 1024 * 1024
PLATFORMS = {
    "GSE49541": {
        "platform_id": "GPL570",
        "title": "[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array",
        "technology": "in situ oligonucleotide",
        "organism": "Homo sapiens",
        "rows": 54_675,
        "selected_columns": (
            "ID",
            "GB_ACC",
            "SPOT_ID",
            "Gene Symbol",
            "ENTREZ_GENE_ID",
            "RefSeq Transcript ID",
        ),
    },
    "GSE83452": {
        "platform_id": "GPL16686",
        "title": (
            "[HuGene-2_0-st] Affymetrix Human Gene 2.0 ST Array "
            "[transcript (gene) version]"
        ),
        "technology": "in situ oligonucleotide",
        "organism": "Homo sapiens",
        "rows": 53_981,
        "selected_columns": (
            "ID",
            "RANGE_STRAND",
            "RANGE_START",
            "RANGE_END",
            "total_probes",
            "GB_ACC",
            "SPOT_ID",
            "RANGE_GB",
        ),
    },
}


class MicroarrayPlatformAuditError(RuntimeError):
    """Raised when a platform, rights, or overlap requirement differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_artifact_hash(path: Path, expected: str, label: str) -> None:
    observed = sha256_file(path / "ARTIFACTS.json")
    if observed != expected:
        raise MicroarrayPlatformAuditError(f"{label} ARTIFACTS SHA-256 differs")


def parse_platform_soft(
    path: Path, expected_platform: str
) -> tuple[dict[str, list[str]], tuple[str, ...], list[dict[str, str]]]:
    metadata: dict[str, list[str]] = {}
    header: tuple[str, ...] = ()
    rows: list[dict[str, str]] = []
    active = False
    in_table = False
    try:
        handle = gzip.open(path, mode="rt", encoding="utf-8", errors="strict", newline="")
        with handle:
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")
                if line.startswith("^PLATFORM = "):
                    platform = line.split(" = ", 1)[1].strip()
                    active = platform == expected_platform
                    in_table = False
                    continue
                if not active:
                    continue
                if line.lower() == "!platform_table_begin":
                    in_table = True
                    continue
                if line.lower() == "!platform_table_end":
                    in_table = False
                    active = False
                    continue
                if in_table:
                    values = next(csv.reader([line], delimiter="\t"))
                    if not header:
                        header = tuple(values)
                        if not header or len(set(header)) != len(header):
                            raise MicroarrayPlatformAuditError(
                                f"{expected_platform} table header differs"
                            )
                    else:
                        if len(values) != len(header):
                            raise MicroarrayPlatformAuditError(
                                f"{expected_platform} table row width differs"
                            )
                        rows.append(dict(zip(header, values, strict=True)))
                    continue
                if line.startswith("!Platform_") and " = " in line:
                    key, value = line.split(" = ", 1)
                    metadata.setdefault(key, []).append(value.strip())
    except (OSError, EOFError, UnicodeDecodeError, csv.Error) as error:
        raise MicroarrayPlatformAuditError(
            f"{expected_platform} SOFT parsing failed"
        ) from error
    if not header or not rows:
        raise MicroarrayPlatformAuditError(f"{expected_platform} table is absent")
    return metadata, header, rows


def validate_platform(
    series: str,
    metadata: dict[str, list[str]],
    header: tuple[str, ...],
    rows: list[dict[str, str]],
) -> dict[str, object]:
    contract = PLATFORMS[series]
    expected_scalar = {
        "!Platform_title": str(contract["title"]),
        "!Platform_technology": str(contract["technology"]),
        "!Platform_organism": str(contract["organism"]),
        "!Platform_data_row_count": str(contract["rows"]),
    }
    for key, expected in expected_scalar.items():
        if metadata.get(key) != [expected]:
            raise MicroarrayPlatformAuditError(
                f"{contract['platform_id']} {key} differs"
            )
    selected = tuple(str(value) for value in contract["selected_columns"])
    if not set(selected).issubset(header):
        raise MicroarrayPlatformAuditError(
            f"{contract['platform_id']} selected annotation columns differ"
        )
    if len(rows) != int(contract["rows"]):
        raise MicroarrayPlatformAuditError(
            f"{contract['platform_id']} row count differs"
        )
    feature_ids = [row["ID"] for row in rows]
    if any(not value for value in feature_ids) or len(set(feature_ids)) != len(feature_ids):
        raise MicroarrayPlatformAuditError(
            f"{contract['platform_id']} feature IDs are absent or duplicated"
        )
    result: dict[str, object] = {
        "platform_id": contract["platform_id"],
        "title": contract["title"],
        "technology": contract["technology"],
        "organism": contract["organism"],
        "feature_rows": len(rows),
        "feature_ids_unique": True,
        "table_columns": list(header),
        "selected_annotation_columns": list(selected),
    }
    if series == "GSE49541":
        mapping = Counter()
        for row in rows:
            values = {
                value.strip()
                for value in row["ENTREZ_GENE_ID"].split("///")
                if value.strip() and value.strip() != "---"
            }
            if not values:
                mapping["unmapped_or_control"] += 1
            elif len(values) == 1:
                mapping["one_entrez_id"] += 1
            else:
                mapping["multiple_entrez_ids"] += 1
        result["geo_entrez_mapping_states"] = dict(sorted(mapping.items()))
    else:
        result["geo_entrez_mapping_states"] = {
            "not_provided_by_GPL16686_table": len(rows)
        }
    return result


def write_selected_platform_rows(
    path: Path,
    series: str,
    rows: list[dict[str, str]],
) -> None:
    contract = PLATFORMS[series]
    selected = tuple(str(value) for value in contract["selected_columns"])
    fields = ("series", "platform_id", *selected)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "series": series,
                    "platform_id": contract["platform_id"],
                    **{field: row[field] for field in selected},
                }
            )


def audit_existing_cel_headers(raw: Path) -> dict[str, object]:
    result: dict[str, object] = {}
    all_hashes: list[str] = []
    all_accessions: list[str] = []
    for series in PLATFORMS:
        path = raw / f"{series}_cel_manifest.tsv"
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "series",
                "sample_accession",
                "decompressed_sha256",
                "header_hex_32",
            }
            if not required.issubset(reader.fieldnames or ()):
                raise MicroarrayPlatformAuditError("raw CEL manifest fields differ")
            rows = [dict(row) for row in reader]
        expected = 72 if series == "GSE49541" else 231
        if len(rows) != expected or any(row["series"] != series for row in rows):
            raise MicroarrayPlatformAuditError(f"{series} raw CEL manifest differs")
        formats = Counter(
            identify_cel_format(bytes.fromhex(row["header_hex_32"])) for row in rows
        )
        result[series] = {
            "cel_records": len(rows),
            "cel_formats": dict(sorted(formats.items())),
            "supported_CEL_signature_exact": True,
        }
        all_hashes.extend(row["decompressed_sha256"] for row in rows)
        all_accessions.extend(row["sample_accession"] for row in rows)
    if len(all_hashes) != len(set(all_hashes)):
        raise MicroarrayPlatformAuditError("cross-cohort decompressed CEL duplicate found")
    if len(all_accessions) != len(set(all_accessions)):
        raise MicroarrayPlatformAuditError("cross-cohort GSM overlap found")
    result["cross_cohort"] = {
        "decompressed_CEL_hash_duplicates": 0,
        "GSM_accession_overlap": 0,
    }
    return result


def fetch_rights() -> dict[str, object]:
    request = urllib.request.Request(
        RIGHTS_URL, headers={"User-Agent": "masld-bench-microarray-rights-audit/1.0"}
    )
    with urllib.request.urlopen(request, timeout=300) as response:
        payload = response.read(RIGHTS_MAX_BYTES + 1)
        headers = {
            key.lower(): value
            for key in ("Content-Length", "Last-Modified", "ETag")
            if (value := response.headers.get(key)) is not None
        }
    if not payload or len(payload) > RIGHTS_MAX_BYTES:
        raise MicroarrayPlatformAuditError("GEO rights authority size differs")
    text = html.unescape(re.sub(r"<[^>]+>", " ", payload.decode("utf-8", "strict")))
    normalized = " ".join(text.split())
    required = (
        "NCBI places no restrictions on the use or distribution of the GEO data",
        "some submitters may claim patent, copyright, or other intellectual property rights",
        "cannot provide comment or unrestricted permission",
    )
    if not all(value in normalized for value in required):
        raise MicroarrayPlatformAuditError("GEO rights authority text differs")
    return {
        "url": RIGHTS_URL,
        "bytes": len(payload),
        "sha256": sha256(payload).hexdigest(),
        "http_headers": headers,
        "authority_statements_verified": True,
        "new_DUA_or_controlled_access_required": False,
        "internal_nonclinical_research_use_allowed": True,
        "raw_CEL_redistribution_by_project": "prohibited_conservative_default",
        "derived_matrix_and_weight_redistribution": "model_and_source_specific_review_required",
    }


def audit_overlap(source: Path) -> dict[str, object]:
    audit = json.loads((source / "source_audit.json").read_text(encoding="utf-8"))
    gse49541 = audit["cohorts"]["GSE49541"]
    gse83452 = audit["cohorts"]["GSE83452"]
    if "SubSeries of: GSE31803" not in gse49541["relations"]:
        raise MicroarrayPlatformAuditError("GSE49541 parent SuperSeries alias differs")
    if any("SubSeries of:" in value for value in gse83452["relations"]):
        raise MicroarrayPlatformAuditError("GSE83452 unexpected SuperSeries relation")
    return {
        "GSE49541": {
            "cohort_family_id": "gse31803_gse49541_fibrosis_array",
            "accession_aliases": ["GSE31803", "GSE49541"],
            "relation": "GSE49541_is_SubSeries_of_GSE31803",
            "same_family_split_required": True,
            "independent_cohort_double_count_allowed": False,
        },
        "GSE83452": {
            "cohort_family_id": "antwerp_inserm_shared",
            "accession_aliases": ["GSE106737", "GSE83452"],
            "all_repeated_records_same_participant_group_required": True,
            "cross_accession_overlap_audit_required_before_split": True,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--source-artifacts-sha256", required=True)
    parser.add_argument("--participants", type=Path, required=True)
    parser.add_argument("--participants-artifacts-sha256", required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--raw-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    require_artifact_hash(arguments.source, arguments.source_artifacts_sha256, "source")
    require_artifact_hash(
        arguments.participants,
        arguments.participants_artifacts_sha256,
        "participant",
    )
    require_artifact_hash(arguments.raw, arguments.raw_artifacts_sha256, "raw")
    arguments.output.mkdir(parents=True, exist_ok=False)
    platform_summary: dict[str, object] = {}
    for series, contract in PLATFORMS.items():
        metadata, header, rows = parse_platform_soft(
            arguments.source / f"raw/{series}_family.soft.gz",
            str(contract["platform_id"]),
        )
        platform_summary[series] = validate_platform(series, metadata, header, rows)
        write_selected_platform_rows(
            arguments.output / f"{contract['platform_id']}_selected_annotation.tsv",
            series,
            rows,
        )
    rights = fetch_rights()
    (arguments.output / "rights_receipt.json").write_text(
        json.dumps(rights, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary = {
        "schema_version": "masld-bench-microarray-platform-rights-overlap-audit-v1",
        "status": "pass_platform_axes_CEL_signatures_rights_and_overlap",
        "platforms": platform_summary,
        "cel_headers": audit_existing_cel_headers(arguments.raw),
        "cohort_overlap": audit_overlap(arguments.source),
        "rights": rights,
        "official_platform_axis_complete": True,
        "raw_CEL_format_complete": True,
        "participant_grouping_complete": True,
        "rights_audit_complete": True,
        "platform_to_entrez_mapping_complete": False,
        "entrez_to_gencode_v49_mapping_complete": False,
        "normalization_run": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
        "next_step": (
            "freeze_existing_single_array_runtime_and_authoritative_annotation_then_"
            "apply_only_outer_training_frozen_transforms_to_each_held_array"
        ),
    }
    (arguments.output / "platform_rights_overlap_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
