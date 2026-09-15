#!/usr/bin/env python3
"""Freeze CEL reader fixtures and GPL16686 annotation availability."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import html
import json
from pathlib import Path
import re
import tarfile
import urllib.error
import urllib.request

from scripts.acquire_audit_microarray_transfer_raw import identify_cel_format


THERMO_PAGE = (
    "https://www.thermofisher.com/us/en/home/life-science/microarray-analysis/"
    "microarray-data-analysis/genechip-array-annotation-files.html"
)
THERMO_ASSET_NAME = "HuGene-2_0-st-v1-na36-hg19-transcript-csv.zip"
NCBI_ANNOT = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL16nnn/GPL16686/annot/"
    "GPL16686.annot.gz"
)


class MicroarrayRuntimePreparationError(RuntimeError):
    """Raised when a frozen reader or annotation fixture differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_artifact_hash(path: Path, expected: str, label: str) -> None:
    if sha256_file(path / "ARTIFACTS.json") != expected:
        raise MicroarrayRuntimePreparationError(f"{label} ARTIFACTS SHA-256 differs")


def read_first_manifest_row(path: Path, expected_series: str) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "series",
            "sample_accession",
            "member_name",
            "decompressed_bytes",
            "decompressed_sha256",
            "header_hex_32",
        }
        if not required.issubset(reader.fieldnames or ()):
            raise MicroarrayRuntimePreparationError("raw CEL manifest fields differ")
        row = next(reader, None)
    if row is None or row["series"] != expected_series:
        raise MicroarrayRuntimePreparationError("raw CEL fixture row differs")
    return dict(row)


def extract_fixture(raw: Path, output: Path, series: str) -> dict[str, object]:
    row = read_first_manifest_row(raw / f"{series}_cel_manifest.tsv", series)
    tar_path = raw / f"raw/{series}_RAW.tar"
    target = output / "fixtures" / f"{series}_{row['sample_accession']}.CEL"
    target.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, mode="r:") as archive:
        matches = [member for member in archive if member.name == row["member_name"]]
        if len(matches) != 1 or not matches[0].isfile():
            raise MicroarrayRuntimePreparationError("exact CEL fixture member is absent")
        source = archive.extractfile(matches[0])
        if source is None:
            raise MicroarrayRuntimePreparationError("CEL fixture member is unreadable")
        decompressor = gzip.GzipFile(fileobj=source, mode="rb")
        digest = sha256()
        total = 0
        header = b""
        with target.open("xb") as handle:
            for block in iter(lambda: decompressor.read(8 * 1024 * 1024), b""):
                if len(header) < 32:
                    header += block[: 32 - len(header)]
                digest.update(block)
                total += len(block)
                handle.write(block)
    if (
        total != int(row["decompressed_bytes"])
        or digest.hexdigest() != row["decompressed_sha256"]
        or header.hex() != row["header_hex_32"]
        or identify_cel_format(header) != "calvin_v1"
    ):
        raise MicroarrayRuntimePreparationError("extracted CEL fixture integrity differs")
    return {
        "series": series,
        "sample_accession": row["sample_accession"],
        "source_tar_member": row["member_name"],
        "path": str(target.relative_to(output)),
        "bytes": total,
        "sha256": digest.hexdigest(),
        "format": "calvin_v1",
        "expression_values_read": False,
    }


def audit_gpl16686_geo_table(platform: Path) -> dict[str, object]:
    path = platform / "GPL16686_selected_annotation.tsv"
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        rows = sum(1 for _ in reader)
    expected = (
        "series",
        "platform_id",
        "ID",
        "RANGE_STRAND",
        "RANGE_START",
        "RANGE_END",
        "total_probes",
        "GB_ACC",
        "SPOT_ID",
        "RANGE_GB",
    )
    if fields != expected or rows != 53_981:
        raise MicroarrayRuntimePreparationError("GPL16686 official GEO table differs")
    return {
        "authority": "NCBI_GEO_GPL16686_platform_table",
        "rows": rows,
        "columns": list(fields),
        "entrez_or_ensembl_gene_column_present": False,
        "coordinates_or_GB_ACC_used_to_infer_gene_mapping": False,
        "mapping_state_for_all_features": "join_unresolved",
    }


def fetch_thermo_authority() -> dict[str, object]:
    request = urllib.request.Request(
        THERMO_PAGE,
        headers={"User-Agent": "masld-bench-microarray-annotation-audit/1.0"},
    )
    with urllib.request.urlopen(request, timeout=300) as response:
        payload = response.read(8 * 1024 * 1024 + 1)
        page_final_url = response.geturl()
    if not payload or len(payload) > 8 * 1024 * 1024:
        raise MicroarrayRuntimePreparationError("Thermo annotation page size differs")
    text = html.unescape(payload.decode("utf-8", "strict"))
    matches = sorted(set(re.findall(r'https?[^"\'<> ]+' + re.escape(THERMO_ASSET_NAME), text)))
    if len(matches) != 1 or "Release 36" not in text or "HuGene-2_0-st-v1" not in text:
        raise MicroarrayRuntimePreparationError("Thermo GPL16686 authority entry differs")
    asset_url = matches[0].replace("&amp;", "&")
    request = urllib.request.Request(
        asset_url,
        method="HEAD",
        headers={"User-Agent": "masld-bench-microarray-annotation-audit/1.0"},
    )
    asset_result: dict[str, object]
    try:
        with urllib.request.urlopen(request, timeout=300) as response:
            asset_result = {
                "http_status": response.status,
                "final_url": response.geturl(),
                "content_type": response.headers.get("Content-Type", ""),
                "content_length": response.headers.get("Content-Length", ""),
            }
    except urllib.error.HTTPError as error:
        asset_result = {
            "http_status": error.code,
            "final_url": error.geturl(),
            "content_type": error.headers.get("Content-Type", ""),
            "content_length": error.headers.get("Content-Length", ""),
        }
    final_url = str(asset_result["final_url"])
    downloadable = (
        int(asset_result["http_status"]) == 200
        and "auth/login" not in final_url
        and "zip" in str(asset_result["content_type"]).lower()
    )
    return {
        "page_url": THERMO_PAGE,
        "page_final_url": page_final_url,
        "page_sha256": sha256(payload).hexdigest(),
        "official_entry": "HuGene-2_0-st-v1 Transcript Cluster Annotations, CSV, Release 36",
        "official_asset_url": asset_url,
        "official_asset_declared_build": "hg19",
        "asset_head": asset_result,
        "openly_downloadable_without_login": downloadable,
    }


def audit_ncbi_annotation_endpoint() -> dict[str, object]:
    request = urllib.request.Request(
        NCBI_ANNOT,
        method="HEAD",
        headers={"User-Agent": "masld-bench-microarray-annotation-audit/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=300) as response:
            status = response.status
            final_url = response.geturl()
    except urllib.error.HTTPError as error:
        status = error.code
        final_url = error.geturl()
    return {"url": NCBI_ANNOT, "http_status": status, "final_url": final_url}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--raw-artifacts-sha256", required=True)
    parser.add_argument("--platform", type=Path, required=True)
    parser.add_argument("--platform-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    require_artifact_hash(arguments.raw, arguments.raw_artifacts_sha256, "raw")
    require_artifact_hash(arguments.platform, arguments.platform_artifacts_sha256, "platform")
    arguments.output.mkdir(parents=True, exist_ok=False)
    fixtures = [
        extract_fixture(arguments.raw, arguments.output, series)
        for series in ("GSE49541", "GSE83452")
    ]
    geo = audit_gpl16686_geo_table(arguments.platform)
    thermo = fetch_thermo_authority()
    ncbi = audit_ncbi_annotation_endpoint()
    if thermo["openly_downloadable_without_login"]:
        mapping_state = "authoritative_asset_accessible_crosswalk_not_yet_parsed"
        availability = "derivable_not_processed"
    else:
        mapping_state = "authoritative_release_identified_but_asset_requires_login"
        availability = "unavailable_permission"
    summary = {
        "schema_version": "masld-bench-microarray-runtime-preparation-v1",
        "status": "pass_reader_fixtures_and_authoritative_annotation_disposition",
        "fixtures": fixtures,
        "gpl16686_geo_platform_table": geo,
        "gpl16686_thermo_authority": thermo,
        "gpl16686_ncbi_annot_endpoint": ncbi,
        "gpl16686_authoritative_mapping_state": mapping_state,
        "gpl16686_mapping_availability": availability,
        "gene_mapping_inferred_from_GB_ACC_or_coordinates": False,
        "gpl16686_to_gencode_v49_crosswalk_complete": False,
        "CEL_expression_values_read": False,
        "normalization_run": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
    }
    (arguments.output / "runtime_preparation_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
