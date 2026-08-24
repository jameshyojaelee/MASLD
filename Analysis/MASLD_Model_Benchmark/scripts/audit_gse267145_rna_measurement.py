#!/usr/bin/env python3
"""Freeze the defensible measurement semantics of deposited GSE267145 RNA."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
import math
from pathlib import Path, PurePosixPath
import re
import tarfile
import urllib.request


class RNAMeasurementError(RuntimeError):
    """Raised when source lineage or deposited matrix semantics differ."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path, maximum_bytes: int = 50_000_000) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "masld-bench/1.0"})
    with urllib.request.urlopen(request, timeout=180) as response, destination.open(
        "xb"
    ) as handle:
        digest = sha256()
        size = 0
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            size += len(block)
            if size > maximum_bytes:
                raise RNAMeasurementError("PISCES source archive exceeds size ceiling")
            handle.write(block)
            digest.update(block)
        headers = {
            key.lower(): value
            for key, value in response.headers.items()
            if key.lower() in {"content-length", "last-modified", "etag"}
        }
    return {"url": url, "bytes": size, "sha256": digest.hexdigest(), "headers": headers}


def audit_source_archive(path: Path, expected_tag_root: str) -> dict:
    required = {
        f"{expected_tag_root}/LICENSE": None,
        f"{expected_tag_root}/pisces/R/make_expression_matrix.R": None,
    }
    member_count = 0
    with tarfile.open(path, "r:gz") as archive:
        for member in archive:
            member_count += 1
            posix = PurePosixPath(member.name)
            if posix.is_absolute() or ".." in posix.parts or member.issym() or member.islnk():
                raise RNAMeasurementError("unsafe PISCES archive member")
            if member.name in required:
                if not member.isfile() or member.size > 5_000_000:
                    raise RNAMeasurementError("required PISCES source member is invalid")
                handle = archive.extractfile(member)
                if handle is None:
                    raise RNAMeasurementError("required PISCES source member is unreadable")
                required[member.name] = handle.read()
    if member_count == 0 or any(value is None for value in required.values()):
        raise RNAMeasurementError("required PISCES source members are missing")
    license_text = required[f"{expected_tag_root}/LICENSE"].decode("utf-8")
    script = required[
        f"{expected_tag_root}/pisces/R/make_expression_matrix.R"
    ].decode("utf-8")
    tokens = (
        "library(tximport)",
        "summarizeToGene(txi",
        "txi.gene$counts",
        ".counts.txt",
        "DESeqDataSetFromTximport",
    )
    if "Apache License" not in license_text or any(token not in script for token in tokens):
        raise RNAMeasurementError("PISCES count-generation evidence differs")
    return {
        "archive_members": member_count,
        "license": "Apache-2.0",
        "make_expression_matrix_sha256": sha256(script.encode()).hexdigest(),
        "evidence_tokens": list(tokens),
    }


def audit_matrix(path: Path, expected_rows: int = 43_285, expected_samples: int = 262) -> dict:
    genes: set[str] = set()
    row_count = 0
    values = 0
    zeros = 0
    fractional = 0
    minimum = math.inf
    maximum = -math.inf
    column_sums: list[float] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) != expected_samples + 1 or header[0] != "ensembl_gene_id":
            raise RNAMeasurementError("RNA matrix header differs")
        if len(set(header[1:])) != expected_samples:
            raise RNAMeasurementError("RNA sample axis is duplicated")
        column_sums = [0.0] * expected_samples
        for row_count, row in enumerate(reader, start=1):
            if len(row) != len(header) or not re.fullmatch(r"ENSG\d+(?:\.\d+)?", row[0]):
                raise RNAMeasurementError(f"RNA row {row_count} differs")
            if row[0] in genes:
                raise RNAMeasurementError("RNA feature axis is duplicated")
            genes.add(row[0])
            for index, text in enumerate(row[1:]):
                value = float(text)
                if not math.isfinite(value) or value < 0:
                    raise RNAMeasurementError("RNA value is not finite nonnegative")
                values += 1
                zeros += int(value == 0)
                fractional += int(not value.is_integer())
                minimum = min(minimum, value)
                maximum = max(maximum, value)
                column_sums[index] += value
    if row_count != expected_rows or not fractional or any(total <= 0 for total in column_sums):
        raise RNAMeasurementError("RNA matrix dimensions or fractional semantics differ")
    return {
        "feature_rows": row_count,
        "samples": expected_samples,
        "numeric_values": values,
        "zero_values": zeros,
        "fractional_values": fractional,
        "fractional_fraction": fractional / values,
        "minimum": minimum,
        "maximum": maximum,
        "minimum_column_sum": min(column_sums),
        "maximum_column_sum": max(column_sums),
        "all_nonnegative": True,
        "all_integer": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix-audit", type=Path, required=True)
    parser.add_argument("--matrix-artifacts-sha256", required=True)
    parser.add_argument("--join-audit", type=Path, required=True)
    parser.add_argument("--join-artifacts-sha256", required=True)
    parser.add_argument("--pisces-url", required=True)
    parser.add_argument("--pisces-ref-url", required=True)
    parser.add_argument("--pisces-tag", required=True)
    parser.add_argument("--pisces-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for root, expected in (
        (args.matrix_audit, args.matrix_artifacts_sha256),
        (args.join_audit, args.join_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise RNAMeasurementError("input ARTIFACTS SHA-256 differs")
    if not re.fullmatch(r"[0-9a-f]{40}", args.pisces_commit):
        raise RNAMeasurementError("PISCES commit is malformed")
    args.output.mkdir(parents=True, exist_ok=False)
    ref_path = args.output / "pisces-tag-ref.json"
    ref_receipt = download(args.pisces_ref_url, ref_path, maximum_bytes=1_000_000)
    ref = json.loads(ref_path.read_text(encoding="utf-8"))
    if (
        ref.get("ref") != f"refs/tags/{args.pisces_tag}"
        or ref.get("object", {}).get("type") != "commit"
        or ref.get("object", {}).get("sha") != args.pisces_commit
    ):
        raise RNAMeasurementError("PISCES tag does not resolve to the frozen commit")
    archive = args.output / f"pisces-{args.pisces_tag}.tar.gz"
    receipt = download(args.pisces_url, archive)
    source = audit_source_archive(archive, f"pisces-{args.pisces_tag.lstrip('v')}")
    matrix = audit_matrix(args.matrix_audit / "raw" / "GSE269412_RNA.txt.gz")
    evidence = {
        "schema_version": "masld-bench-gse267145-rna-measurement-v1",
        "status": "pass_with_constrained_use",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "deposited_measurement": "nonnegative_fractional_gene_expression_estimates",
        "not_raw_integer_molecule_counts": True,
        "public_pisces_prepublication_tag": args.pisces_tag,
        "public_pisces_commit": args.pisces_commit,
        "public_code_is_exact_private_generator": False,
        "source_lineage_interpretation": "The last public PISCES tag before the deposit uses Salmon, tximport, summarizeToGene, and tximport-aware DESeq2. This explains why a file called counts can contain fractional abundance estimates, but the paper does not pin the exact PISCES/private exon-pipeline revision.",
        "allowed_model_input": "continuous_expression_only_after_outer_training_fold_fit",
        "allowed_transforms": ["log1p_library_scaled", "rank_or_quantile_mapper"],
        "forbidden_uses": [
            "raw_integer_count_likelihood",
            "rounding_or_flooring_to_manufacture_counts",
            "DESeq2_reproduction_without_tximport_lengths_and_offsets",
            "normalization_fit_on_held_participants",
        ],
        "matrix": matrix,
        "tag_ref_receipt": ref_receipt,
        "source_archive": {**receipt, **source},
        "matrix_artifacts_sha256": args.matrix_artifacts_sha256,
        "join_artifacts_sha256": args.join_artifacts_sha256,
    }
    (args.output / "measurement_evidence.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(evidence, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
