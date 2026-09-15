#!/usr/bin/env python3
"""Resumably download and verify one participant-contained GSE268273 bundle."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
from typing import Any

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


class GSE268273DownloadError(RuntimeError):
    """Raised when a raw bundle cannot be included safely."""


def file_digests(path: Path) -> tuple[str, str]:
    md5_digest = hashlib.md5(usedforsecurity=False)
    sha256_digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            md5_digest.update(block)
            sha256_digest.update(block)
    return md5_digest.hexdigest(), sha256_digest.hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_bundle_rows(plan_root: Path, bundle_id: int) -> list[dict[str, str]]:
    verify_frozen_tree(plan_root)
    with (plan_root / "fastq_files.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "bundle_id",
            "row_id",
            "run_accession",
            "experiment_accession",
            "biosample_accession",
            "declared_ena_library_layout",
            "effective_library_layout",
            "layout_evidence",
            "ena_read_count",
            "ena_base_count",
            "file_part_index",
            "file_part_count",
            "read_role",
            "https_url",
            "relative_path",
            "fastq_bytes",
            "fastq_md5",
        }
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise GSE268273DownloadError("frozen FASTQ plan schema differs")
        forbidden = {"fibrosis", "outcome", "label", "disease", "sex", "nas", "bmi"}
        if forbidden & set(reader.fieldnames):
            raise GSE268273DownloadError("FASTQ plan violates the label firewall")
        rows = [row for row in reader if int(row["bundle_id"]) == bundle_id]
    if not rows:
        raise GSE268273DownloadError(f"bundle {bundle_id} is absent from the plan")
    return rows


def fastq_first_record(path: Path) -> tuple[str, str]:
    try:
        with gzip.open(path, "rt", encoding="ascii") as handle:
            header = handle.readline().rstrip("\n")
            sequence = handle.readline().rstrip("\n")
            separator = handle.readline().rstrip("\n")
            quality = handle.readline().rstrip("\n")
    except (OSError, UnicodeDecodeError) as error:
        raise GSE268273DownloadError(f"FASTQ first record is unreadable: {path}") from error
    if (
        not header.startswith("@")
        or not sequence
        or not separator.startswith("+")
        or len(sequence) != len(quality)
    ):
        raise GSE268273DownloadError(f"FASTQ first record differs: {path}")
    tokens = header.split()
    if len(tokens) < 2 or not tokens[-1].endswith(("/1", "/2")):
        raise GSE268273DownloadError(f"FASTQ first-read mate tag differs: {path}")
    return tokens[0], tokens[-1][-1]


def validate_member(
    path: Path,
    expected_bytes: int,
    expected_md5: str,
    expected_read_role: str | None = None,
) -> str:
    if path.stat().st_size != expected_bytes:
        raise GSE268273DownloadError(f"FASTQ byte size differs: {path}")
    observed_md5, observed_sha256 = file_digests(path)
    if observed_md5 != expected_md5:
        raise GSE268273DownloadError(f"FASTQ ENA MD5 differs: {path}")
    result = subprocess.run(
        ["gzip", "-t", str(path)],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise GSE268273DownloadError(
            f"FASTQ gzip integrity differs: {path}: {result.stderr.strip()}"
        )
    if expected_read_role is not None:
        _, observed_mate = fastq_first_record(path)
        expected_mate = {
            "single_end_read": "1",
            "paired_end_mate_1": "1",
            "paired_end_mate_2": "2",
        }.get(expected_read_role)
        if expected_mate is None or observed_mate != expected_mate:
            raise GSE268273DownloadError(
                f"FASTQ first-read role differs: {path}: {expected_read_role}"
            )
    return observed_sha256


def download_bundle(
    *, plan_root: Path, state_root: Path, bundle_id: int
) -> dict[str, Any]:
    rows = read_bundle_rows(plan_root, bundle_id)
    layout_by_participant: defaultdict[str, set[str]] = defaultdict(set)
    for row in rows:
        layout_by_participant[row["row_id"]].add(row["effective_library_layout"])
    if any(len(layouts) != 1 for layouts in layout_by_participant.values()):
        raise GSE268273DownloadError("one participant contains mixed read layouts")
    layout_participants = Counter(
        next(iter(layouts)) for layouts in layout_by_participant.values()
    )
    bundle = state_root / f"bundle_{bundle_id:02d}"
    if (bundle / "COMPLETE").exists():
        verify_frozen_tree(bundle)
        return json.loads((bundle / "bundle_receipt.json").read_text(encoding="utf-8"))
    bundle.mkdir(parents=True, exist_ok=True)
    receipts = bundle / "file_receipts"
    receipts.mkdir(exist_ok=True)
    for row in rows:
        relative = PurePosixPath(row["relative_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise GSE268273DownloadError("FASTQ relative path is unsafe")
        destination = bundle.joinpath(*relative.parts)
        destination.parent.mkdir(parents=True, exist_ok=True)
        expected_bytes = int(row["fastq_bytes"])
        expected_md5 = row["fastq_md5"]
        receipt_path = receipts / f"{row['run_accession']}.{row['file_part_index']}.json"
        if destination.exists():
            observed_sha256 = validate_member(
                destination, expected_bytes, expected_md5, row["read_role"]
            )
        else:
            partial = destination.with_name(destination.name + ".partial")
            result = subprocess.run(
                [
                    "curl",
                    "--fail",
                    "--location",
                    "--retry",
                    "8",
                    "--retry-all-errors",
                    "--connect-timeout",
                    "30",
                    "--continue-at",
                    "-",
                    "--output",
                    str(partial),
                    row["https_url"],
                ],
                check=False,
            )
            if result.returncode != 0:
                raise GSE268273DownloadError(
                    f"curl failed for {row['run_accession']} part {row['file_part_index']}"
                )
            try:
                observed_sha256 = validate_member(
                    partial, expected_bytes, expected_md5, row["read_role"]
                )
            except GSE268273DownloadError:
                attempt = os.environ.get("SLURM_JOB_ID", "manual")
                quarantine = destination.with_name(
                    destination.name + f".invalid-attempt-{attempt}"
                )
                if quarantine.exists():
                    raise GSE268273DownloadError(
                        f"invalid FASTQ quarantine already exists: {quarantine}"
                    )
                partial.rename(quarantine)
                raise
            if destination.exists():
                raise GSE268273DownloadError(
                    f"refusing to overwrite concurrently created FASTQ: {destination}"
                )
            partial.rename(destination)
        observed = {
            "schema_version": "masld-bench-gse268273-fastq-member-v1",
            "bundle_id": bundle_id,
            "row_id": row["row_id"],
            "run_accession": row["run_accession"],
            "file_part_index": int(row["file_part_index"]),
            "relative_path": row["relative_path"],
            "bytes": expected_bytes,
            "ena_md5": expected_md5,
            "sha256": observed_sha256,
            "gzip_integrity": "passed",
            "labels_accessed": False,
        }
        if receipt_path.exists():
            if json.loads(receipt_path.read_text(encoding="utf-8")) != observed:
                raise GSE268273DownloadError("existing FASTQ receipt differs")
        else:
            receipt_path.write_text(
                json.dumps(observed, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
    by_run: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_run[row["run_accession"]].append(row)
    for run, run_rows in by_run.items():
        if run_rows[0]["effective_library_layout"] != "paired_end":
            continue
        ordered = sorted(run_rows, key=lambda row: int(row["file_part_index"]))
        if len(ordered) != 2:
            raise GSE268273DownloadError("paired-end run lacks exactly two mates")
        identities = []
        for row in ordered:
            relative = PurePosixPath(row["relative_path"])
            identities.append(
                fastq_first_record(bundle.joinpath(*relative.parts))[0]
            )
        if identities[0] != identities[1]:
            raise GSE268273DownloadError(
                f"paired FASTQ first-read identities differ: {run}"
            )
    receipt = {
        "schema_version": "masld-bench-gse268273-raw-bundle-v1",
        "status": "complete_verified",
        "bundle_id": bundle_id,
        "participants": len({row["row_id"] for row in rows}),
        "effective_layout_participants": dict(sorted(layout_participants.items())),
        "technical_runs": len({row["run_accession"] for row in rows}),
        "fastq_files": len(rows),
        "fastq_bytes": sum(int(row["fastq_bytes"]) for row in rows),
        "ena_md5_verified_for_every_file": True,
        "gzip_integrity_verified_for_every_file": True,
        "fastq_first_record_layout_verified_for_every_file": True,
        "paired_first_record_identity_verified_for_every_run": True,
        "technical_files_not_biological_replicates": True,
        "labels_accessed": False,
        "fit_or_score_performed": False,
        "plan_artifacts_sha256": sha256_file(plan_root / "ARTIFACTS.json"),
    }
    (bundle / "bundle_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(
        bundle,
        {
            "artifact_class": "gse268273_verified_raw_fastq_bundle",
            "bundle_id": bundle_id,
            "participants": receipt["participants"],
            "fastq_files": receipt["fastq_files"],
            "fastq_bytes": receipt["fastq_bytes"],
            "status": "passed",
        },
    )
    verify_frozen_tree(bundle)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--state-root", type=Path, required=True)
    parser.add_argument("--bundle-id", type=int, required=True)
    args = parser.parse_args()
    receipt = download_bundle(
        plan_root=args.plan_root,
        state_root=args.state_root,
        bundle_id=args.bundle_id,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
