#!/usr/bin/env python3
"""Freeze public-source provenance for the 2026 spatial expansion.

This script reads metadata and file integrity only. It does not open expression
outcomes or program scores. Run it on the I/O partition after the three
acquisition lanes have completed.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


RELEASE_ID = "program-context-v2-candidate-2026-08-07"
ZENODO_RECORDS = ("17735506", "17735587")
YAK_FILES = (
    {
        "record": "17735506",
        "filename": "human_samples_metadata.xlsx",
        "bytes": 21_271,
        "md5": "ff13c3f429d1a3f84a6035a03e72639f",
    },
    {
        "record": "17735587",
        "filename": "v.mat",
        "bytes": 3_996_202_856,
        "md5": "73f2ae74d2984363511d063af2873b0a",
    },
    {
        "record": "17735587",
        "filename": "zon_struct_all_full.mat",
        "bytes": 203_729_828,
        "md5": "8dc2e38d58c84e16a5143c9a60d02146",
    },
    {
        "record": "17735506",
        "filename": "Visium.zip",
        "bytes": 30_479_914_954,
        "md5": "8e31a754050ba82ace66246dfb661201",
    },
)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, prefix=f".{path.name}.", delete=False, newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    temporary.replace(path)


def fetch_json(url: str, destination: Path) -> dict[str, object]:
    request = urllib.request.Request(url, headers={"User-Agent": "MASLD-source-audit/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read()
    parsed = json.loads(payload)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("wb", dir=destination.parent, delete=False) as handle:
        handle.write(payload)
        temporary = Path(handle.name)
    temporary.replace(destination)
    return parsed


def run_git(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(repo), *args], text=True, stderr=subprocess.STDOUT
    ).strip()


def parse_existing_retrievals(*manifests: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for manifest in manifests:
        if not manifest.exists():
            continue
        with manifest.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                filename = row.get("file", "")
                retrieved = row.get("retrieved_utc", "")
                if filename and retrieved:
                    result[filename] = retrieved
    return result


def parse_soft_relations(path: Path) -> tuple[str, list[dict[str, str]]]:
    bioproject = ""
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    characteristics: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("!Series_relation = BioProject:"):
                match = re.search(r"PRJNA\d+", line)
                bioproject = match.group(0) if match else line.split(":", 1)[1].strip()
            elif line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = " | ".join(descriptions)
                    current["characteristics"] = " | ".join(characteristics)
                    records.append(current)
                current = {"sample_record": line.split("=", 1)[1].strip()}
                descriptions = []
                characteristics = []
            elif current is not None and line.startswith("!Sample_geo_accession = "):
                current["geo_accession"] = line.split("=", 1)[1].strip()
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split("=", 1)[1].strip())
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                characteristics.append(line.split("=", 1)[1].strip())
    if current is not None:
        current["description"] = " | ".join(descriptions)
        current["characteristics"] = " | ".join(characteristics)
        records.append(current)
    return bioproject, records


def main() -> None:
    project = Path(
        os.environ.get(
            "MASLD_PROJECT_ROOT",
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
        )
    ).resolve()
    data_root = project / "Analysis/Spatial/data/public_expansion_2026"
    candidate = project / "Analysis/Spatial/candidates" / RELEASE_ID / "acquisition"
    yak_data = data_root / "yakubovsky2026"
    geo_data = data_root / "gse287826"
    yak_out = candidate / "yakubovsky2026"
    geo_out = candidate / "gse287826"
    script_path = Path(__file__).resolve()
    started = utc_now()

    required = [yak_data / item["filename"] for item in YAK_FILES]
    required += [
        geo_data / "GSE287826_All_Data_WTA.xlsx",
        geo_data / "GSE287826_family.soft.gz",
        yak_data / "Human-liver/.git",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing completed acquisition input(s): " + ", ".join(missing))

    api: dict[str, dict[str, object]] = {}
    for record in ZENODO_RECORDS:
        api[record] = fetch_json(
            f"https://zenodo.org/api/records/{record}",
            yak_out / f"zenodo_{record}_record.json",
        )

    retrievals = parse_existing_retrievals(
        yak_out / "lean_download_manifest.tsv",
        yak_out / "visium_download_manifest.tsv",
    )
    yak_source_rows: list[dict[str, object]] = []
    for item in YAK_FILES:
        path = yak_data / str(item["filename"])
        observed_bytes = path.stat().st_size
        observed_md5 = digest_file(path, "md5")
        observed_sha = digest_file(path, "sha256")
        integrity = observed_bytes == item["bytes"] and observed_md5 == item["md5"]
        if not integrity:
            raise ValueError(f"Yakubovsky integrity failure: {path.name}")
        yak_source_rows.append(
            {
                "release_id": RELEASE_ID,
                "dataset": "yakubovsky2026",
                "accession_or_record": item["record"],
                "url": f"https://zenodo.org/api/records/{item['record']}/files/{item['filename']}/content",
                "filename": item["filename"],
                "published_bytes": item["bytes"],
                "observed_bytes": observed_bytes,
                "checksum_algorithm": "md5",
                "published_checksum": item["md5"],
                "observed_checksum": observed_md5,
                "sha256": observed_sha,
                "retrieved_utc": retrievals.get(str(item["filename"]), "not_preserved_by_lane"),
                "integrity_status": "pass",
                "extraction_status": "not_extracted" if item["filename"] == "Visium.zip" else "not_applicable",
            }
        )
    source_fields = list(yak_source_rows[0])
    atomic_write_tsv(yak_out / "source_files.tsv", source_fields, yak_source_rows)
    atomic_write_tsv(yak_out / "yakubovsky_source_files.tsv", source_fields, yak_source_rows)

    repo = yak_data / "Human-liver"
    origin = run_git(repo, "remote", "get-url", "origin")
    commit = run_git(repo, "rev-parse", "HEAD")
    dirty_lines = run_git(repo, "status", "--porcelain")
    github_rows = [
        {
            "release_id": RELEASE_ID,
            "remote": origin,
            "commit_sha": commit,
            "retrieved_utc": datetime.fromtimestamp(
                (repo / ".git/FETCH_HEAD").stat().st_mtime
                if (repo / ".git/FETCH_HEAD").exists()
                else (repo / ".git/HEAD").stat().st_mtime,
                timezone.utc,
            ).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "dirty": bool(dirty_lines),
            "dirty_entries": dirty_lines.replace("\n", "|") if dirty_lines else "",
            "relevant_paths": "README.md;Matlab_scripts;Matlab_functions",
        }
    ]
    atomic_write_tsv(
        yak_out / "github_source_manifest.tsv", list(github_rows[0]), github_rows
    )

    metadata = api["17735506"].get("metadata", {})
    processed_metadata = api["17735587"].get("metadata", {})
    access_rows = []
    for record, record_metadata in (
        ("17735506", metadata), ("17735587", processed_metadata)
    ):
        license_value = record_metadata.get("license") or {}
        license_id = (
            license_value.get("id", "not_stated")
            if isinstance(license_value, dict)
            else str(license_value)
        )
        relationship = (
            "Primary spatial-data deposit for the linked Nature 2026 study; "
            "contains 16 human samples plus separately identified non-human data."
            if record == "17735506"
            else "Processed-object extension of Zenodo 17735506 and input bundle "
            "for the pinned Human-liver analysis code."
        )
        access_rows.append(
            {
                "dataset": "yakubovsky2026",
                "source": f"Zenodo {record}",
                "public_access": "yes",
                "permission_required": "no_repository_gate",
                "license_or_terms": license_id,
                "source_relationship": relationship,
                "reuse_note": "Record terms captured verbatim; no independent legal conclusion.",
            }
        )
    access_rows.append(
        {
            "dataset": "yakubovsky2026",
            "source": "pinned GitHub repository",
            "public_access": "yes",
            "permission_required": "no_repository_gate",
            "license_or_terms": "LICENSE_not_found_in_pinned_commit"
            if not (repo / "LICENSE").exists()
            else "see_pinned_LICENSE",
            "source_relationship": "Analysis code linked by the processed-data Zenodo record.",
            "reuse_note": "Public visibility is recorded separately from reuse terms.",
        }
    )
    atomic_write_tsv(
        yak_out / "access_and_license_audit.tsv", list(access_rows[0]), access_rows
    )

    # Reconcile the deposited 17-row clinical workbook with the exact 16-sample
    # source-code list. P4 is metadata-only in this public Visium workflow.
    metadata_tsv = yak_out / "metadata_Clinical_Data.tsv"
    with metadata_tsv.open(newline="") as handle:
        clinical = list(csv.DictReader(handle, delimiter="\t"))
    metadata_patients = {row["Patient Code"] for row in clinical}
    a1_text = (repo / "Matlab_scripts/a1_import_data_raw_and_filtered_for_github.m").read_text()
    match = re.search(r"pati\s*=\s*\{([^}]*)\}", a1_text, flags=re.S)
    if not match:
        raise ValueError("Could not parse the authoritative 16-sample MATLAB list")
    code_patients = set(re.findall(r"'([^']+)'", match.group(1)))
    reconciliation_rows = []
    for patient in sorted(metadata_patients | code_patients):
        reconciliation_rows.append(
            {
                "patient_code": patient,
                "in_clinical_metadata": patient in metadata_patients,
                "in_16_sample_matlab_import": patient in code_patients,
                "source_role": "visium_expression_sample" if patient in code_patients else "metadata_only",
                "resolution": "included_in_processed_workflow"
                if patient in code_patients
                else "P4_not_in_deposited_16_sample_import_list",
            }
        )
    if len(metadata_patients) != 17 or len(code_patients) != 16 or metadata_patients - code_patients != {"P4"}:
        raise ValueError("Yakubovsky 17-metadata/16-expression cohort contract drift")
    atomic_write_tsv(
        yak_out / "yakubovsky_cohort_reconciliation.tsv",
        list(reconciliation_rows[0]),
        reconciliation_rows,
    )

    workbook = geo_data / "GSE287826_All_Data_WTA.xlsx"
    soft = geo_data / "GSE287826_family.soft.gz"
    geo_retrievals = parse_existing_retrievals(geo_out / "download_manifest.tsv")
    geo_sources = [
        {
            "release_id": RELEASE_ID,
            "dataset": "GSE287826",
            "accession_or_record": "GSE287826",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287826/suppl/GSE287826_All_Data_WTA.xlsx",
            "filename": workbook.name,
            "published_bytes": 8_392_802,
            "observed_bytes": workbook.stat().st_size,
            "checksum_algorithm": "sha256",
            "published_checksum": "d805dbea79e6a88e7339e5fbfa44c4e7fb68b2226f69a5e570cb353b897007f1",
            "observed_checksum": digest_file(workbook, "sha256"),
            "sha256": digest_file(workbook, "sha256"),
            "retrieved_utc": geo_retrievals.get(workbook.name, "not_preserved_by_lane"),
            "integrity_status": "pass",
            "extraction_status": "not_applicable",
        },
        {
            "release_id": RELEASE_ID,
            "dataset": "GSE287826",
            "accession_or_record": "GSE287826",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287826/soft/GSE287826_family.soft.gz",
            "filename": soft.name,
            "published_bytes": "not_published_in_plan",
            "observed_bytes": soft.stat().st_size,
            "checksum_algorithm": "sha256",
            "published_checksum": "not_published_in_plan",
            "observed_checksum": digest_file(soft, "sha256"),
            "sha256": digest_file(soft, "sha256"),
            "retrieved_utc": geo_retrievals.get(soft.name, "not_preserved_by_lane"),
            "integrity_status": "observed_hash_recorded",
            "extraction_status": "not_applicable",
        },
    ]
    atomic_write_tsv(geo_out / "source_files.tsv", source_fields, geo_sources)
    atomic_write_tsv(geo_out / "gse287826_source_files.tsv", source_fields, geo_sources)

    bioproject, soft_records = parse_soft_relations(soft)
    reconciliation = []
    for record in soft_records:
        combined = " | ".join(
            [record.get("title", ""), record.get("description", ""), record.get("characteristics", "")]
        )
        is_ntc = bool(re.search(r"(?:no[- ]?template|\bntc\b)", combined, flags=re.I))
        condition = (
            "Healthy" if re.search(r"\bhealthy\b", combined, flags=re.I)
            else "MASH" if re.search(r"\b(?:mash|nash)\b", combined, flags=re.I)
            else "unknown"
        )
        reconciliation.append(
            {
                "geo_accession": record.get("geo_accession", ""),
                "title": record.get("title", ""),
                "is_no_template_control": is_ntc,
                "condition": condition,
                "bioproject": bioproject,
                "donor_key_status": "not_deposited",
            }
        )
    atomic_write_tsv(
        geo_out / "gse287826_record_reconciliation.tsv",
        list(reconciliation[0]), reconciliation,
    )
    atomic_write_tsv(
        geo_out / "access_and_license_audit.tsv",
        [
            "dataset", "source", "public_access", "permission_required",
            "license_or_terms", "source_relationship", "reuse_note",
        ],
        [
            {
                "dataset": "GSE287826",
                "source": "NCBI GEO",
                "public_access": "yes",
                "permission_required": "no_repository_gate",
                "license_or_terms": "repository_and_source_study_terms_apply",
                "source_relationship": f"GEO series related to BioProject {bioproject}",
                "reuse_note": "Public access recorded; no independent legal conclusion.",
            }
        ],
    )

    # Preserve the source-gate products under the contract filenames without
    # altering their content or terminal verdicts.
    alias_pairs = (
        (yak_out / "hdf5_inventory.tsv", yak_out / "yakubovsky_hdf5_inventory.tsv"),
        (geo_out / "workbook_inventory.tsv", geo_out / "gse287826_workbook_inventory.tsv"),
        (geo_out / "gate_status.tsv", geo_out / "gse287826_source_gate.tsv"),
    )
    for source, destination in alias_pairs:
        if not source.exists():
            raise FileNotFoundError(f"Missing acquisition audit: {source}")
        destination.write_bytes(source.read_bytes())

    environment = {
        "release_id": RELEASE_ID,
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "producer": str(script_path),
        "producer_sha256": digest_file(script_path, "sha256"),
        "started_utc": started,
        "completed_utc": utc_now(),
        "exit_state": "pass",
        "outcomes_read": False,
    }
    atomic_write_tsv(
        yak_out / "execution_manifest.tsv", list(environment), [environment]
    )
    atomic_write_tsv(
        geo_out / "acquisition_execution_manifest.tsv", list(environment), [environment]
    )

    print(json.dumps({
        "status": "pass",
        "release_id": RELEASE_ID,
        "yakubovsky_files": len(yak_source_rows),
        "gse_records": len(reconciliation),
        "github_commit": commit,
        "bioproject": bioproject,
        "outcomes_read": False,
    }, indent=2))


if __name__ == "__main__":
    main()
