#!/usr/bin/env python3
"""Strictly admit one-member gzip inputs for GSE105127 revision v2.

The preflight authenticates compressed bytes and fully inflates each stream to
verify its gzip trailer.  It does not parse RNA, RRBS, or model measurements.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
from typing import Any
import zlib

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


PLAN_ARTIFACTS_SHA256 = "c4f91d2a3fedc7b3b1802a0d6cf3d1b1803ce1bea70cc343f6e18f744783b973"
REFERENCE_ARTIFACTS_SHA256 = "9575aa24dc39f855b961eecf53e29243b4a45fca6770c7337e05901bc7d1ec89"
TARGET_FASTA_SHA256 = "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215"
TARGET_GTF_SHA256 = "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
FORWARD_CHAIN_MD5 = "35887f73fe5e2231656504d1f6430900"
REVERSE_CHAIN_MD5 = "ff3031d93792f4cbb86af44055efd903"
FORBIDDEN_FIELDS = frozenset(
    {"phenotype", "label", "outcome", "disease", "fibrosis", "nas", "sex", "age", "bmi", "outer_fold"}
)


class GSE105127GzipPreflightError(ValueError):
    """Raised when a revision-v2 compressed input is not one exact gzip member."""


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE105127GzipPreflightError(f"TSV lacks a header: {path}")
        fields = tuple(reader.fieldnames)
        if FORBIDDEN_FIELDS & set(fields):
            raise GSE105127GzipPreflightError("source plan violates the outcome firewall")
        return fields, list(reader)


def digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def strict_single_member_gzip(
    path: Path,
    *,
    expected_sha256: str,
    expected_bytes: int | None = None,
    expected_md5: str | None = None,
) -> dict[str, Any]:
    if not path.is_file():
        raise GSE105127GzipPreflightError(f"compressed input is absent: {path}")
    if len(expected_sha256) != 64 or any(
        character not in "0123456789abcdef" for character in expected_sha256
    ):
        raise GSE105127GzipPreflightError("expected compressed SHA-256 is invalid")
    if expected_md5 is not None and (
        len(expected_md5) != 32
        or any(character not in "0123456789abcdef" for character in expected_md5)
    ):
        raise GSE105127GzipPreflightError("expected compressed MD5 is invalid")
    compressed_sha = hashlib.sha256()
    compressed_md5 = hashlib.md5(usedforsecurity=False)
    uncompressed_sha = hashlib.sha256()
    decompressor = zlib.decompressobj(wbits=31)
    compressed_bytes = 0
    uncompressed_bytes = 0
    trailing_bytes = 0
    trailing_prefix = bytearray()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            compressed_sha.update(block)
            compressed_md5.update(block)
            compressed_bytes += len(block)
            if decompressor.eof:
                trailing_bytes += len(block)
                trailing_prefix.extend(block[: max(0, 16 - len(trailing_prefix))])
                continue
            try:
                output = decompressor.decompress(block)
            except zlib.error as error:
                raise GSE105127GzipPreflightError(
                    "gzip member failed DEFLATE or trailer validation"
                ) from error
            uncompressed_sha.update(output)
            uncompressed_bytes += len(output)
            if decompressor.eof and decompressor.unused_data:
                unused = decompressor.unused_data
                trailing_bytes += len(unused)
                trailing_prefix.extend(unused[:16])
    try:
        output = decompressor.flush()
    except zlib.error as error:
        raise GSE105127GzipPreflightError(
            "gzip member failed final trailer validation"
        ) from error
    uncompressed_sha.update(output)
    uncompressed_bytes += len(output)
    if not decompressor.eof or decompressor.unconsumed_tail:
        raise GSE105127GzipPreflightError("gzip member is truncated or incomplete")
    if trailing_bytes:
        disposition = (
            "additional_gzip_member"
            if bytes(trailing_prefix[:2]) == b"\x1f\x8b"
            else "trailing_garbage"
        )
        raise GSE105127GzipPreflightError(
            f"gzip stream is not exactly one member: {disposition}"
        )
    observed_sha = compressed_sha.hexdigest()
    observed_md5 = compressed_md5.hexdigest()
    if observed_sha != expected_sha256:
        raise GSE105127GzipPreflightError("compressed SHA-256 differs")
    if expected_bytes is not None and compressed_bytes != expected_bytes:
        raise GSE105127GzipPreflightError("compressed byte count differs")
    if expected_md5 is not None and observed_md5 != expected_md5:
        raise GSE105127GzipPreflightError("compressed MD5 differs")
    result = subprocess.run(
        ["gzip", "-t", str(path)], capture_output=True, text=True, check=False
    )
    if result.returncode != 0 or result.stdout or result.stderr:
        raise GSE105127GzipPreflightError("gzip -t did not pass silently")
    return {
        "path": str(path.resolve()),
        "compressed_bytes": compressed_bytes,
        "compressed_sha256": observed_sha,
        "compressed_md5": observed_md5,
        "gzip_member_count": 1,
        "trailing_bytes": 0,
        "uncompressed_bytes": uncompressed_bytes,
        "uncompressed_sha256": uncompressed_sha.hexdigest(),
        "gzip_test_exit_status": 0,
        "gzip_test_passed": True,
    }


def require_frozen(root: Path, expected_sha256: str | None = None) -> dict[str, Any]:
    if expected_sha256 is not None and sha256_file(root / "ARTIFACTS.json") != expected_sha256:
        raise GSE105127GzipPreflightError(f"frozen artifact identity differs: {root}")
    return verify_frozen_tree(root)


def preflight_rsem_reference(target_fasta: Path, target_gtf: Path) -> list[dict[str, Any]]:
    return [
        strict_single_member_gzip(target_fasta, expected_sha256=TARGET_FASTA_SHA256),
        strict_single_member_gzip(target_gtf, expected_sha256=TARGET_GTF_SHA256),
    ]


def preflight_rrbs(plan_root: Path, bundle_id: int) -> list[dict[str, Any]]:
    require_frozen(plan_root, PLAN_ARTIFACTS_SHA256)
    _, rows = read_tsv(plan_root / "rrbs_rows.tsv")
    selected = [row for row in rows if int(row["bundle_id"]) == bundle_id]
    if not selected:
        raise GSE105127GzipPreflightError("RRBS bundle is empty")
    output = []
    for row in selected:
        result = strict_single_member_gzip(
            Path(row["source_bed_path"]),
            expected_sha256=row["source_bed_sha256"],
            expected_bytes=int(row["source_bed_bytes"]),
        )
        result.update({"row_id": row["row_id"], "input_role": "rrbs_source_bed"})
        output.append(result)
    return output


def preflight_rna(plan_root: Path, raw_root: Path, bundle_id: int) -> list[dict[str, Any]]:
    require_frozen(plan_root, PLAN_ARTIFACTS_SHA256)
    _, rows = read_tsv(plan_root / "rna_rows.tsv")
    selected = [row for row in rows if int(row["bundle_id"]) == bundle_id]
    if not selected:
        raise GSE105127GzipPreflightError("RNA bundle is empty")
    output = []
    for row in selected:
        participant = raw_root / "participants" / row["row_id"]
        require_frozen(participant)
        receipt = json.loads((participant / "receipt.json").read_text(encoding="utf-8"))
        if (
            receipt.get("status") != "complete_verified"
            or receipt.get("labels_accessed") is not False
            or receipt.get("fit_or_score_performed") is not False
            or receipt.get("row_id") != row["row_id"]
            or receipt.get("md5") != row["fastq_md5"]
        ):
            raise GSE105127GzipPreflightError("frozen RNA source receipt differs")
        result = strict_single_member_gzip(
            participant / "reads.fastq.gz",
            expected_sha256=str(receipt["sha256"]),
            expected_bytes=int(row["fastq_bytes"]),
            expected_md5=row["fastq_md5"],
        )
        result.update({"row_id": row["row_id"], "input_role": "rna_fastq"})
        output.append(result)
    return output


def preflight_crosswalk(
    plan_root: Path, reference_root: Path, collapsed_root: Path
) -> list[dict[str, Any]]:
    require_frozen(plan_root, PLAN_ARTIFACTS_SHA256)
    require_frozen(reference_root, REFERENCE_ARTIFACTS_SHA256)
    reference_receipt = json.loads((reference_root / "receipt.json").read_text(encoding="utf-8"))
    output = []
    for name, expected_md5, receipt_key in (
        ("hg19ToHg38.over.chain.gz", FORWARD_CHAIN_MD5, "forward_chain_sha256"),
        ("hg38ToHg19.over.chain.gz", REVERSE_CHAIN_MD5, "reverse_chain_sha256"),
    ):
        path = reference_root / name
        result = strict_single_member_gzip(
            path,
            expected_sha256=str(reference_receipt[receipt_key]),
            expected_md5=expected_md5,
        )
        result["input_role"] = "lift_chain"
        output.append(result)
    _, rows = read_tsv(plan_root / "rrbs_rows.tsv")
    if len(rows) != 57:
        raise GSE105127GzipPreflightError("crosswalk RRBS row census differs")
    for row in rows:
        participant = collapsed_root / "participants" / row["row_id"]
        require_frozen(participant)
        receipt = json.loads((participant / "receipt.json").read_text(encoding="utf-8"))
        if (
            receipt.get("status") != "passed"
            or receipt.get("row_id") != row["row_id"]
            or receipt.get("labels_accessed") is not False
            or receipt.get("fit_or_score_performed") is not False
        ):
            raise GSE105127GzipPreflightError("revised collapsed RRBS receipt differs")
        result = strict_single_member_gzip(
            participant / "cpg_counts.hg19.tsv.gz",
            expected_sha256=str(receipt["collapsed_sha256"]),
        )
        result.update({"row_id": row["row_id"], "input_role": "collapsed_rrbs"})
        output.append(result)
    return output


def run(
    *,
    stage: str,
    plan_root: Path | None,
    target_fasta: Path | None,
    target_gtf: Path | None,
    raw_root: Path | None,
    reference_root: Path | None,
    collapsed_root: Path | None,
    bundle_id: int | None,
    output: Path,
) -> None:
    if output.exists():
        raise GSE105127GzipPreflightError("refusing to overwrite gzip preflight")
    if stage == "rsem_reference":
        if target_fasta is None or target_gtf is None:
            raise GSE105127GzipPreflightError("RSEM preflight inputs are incomplete")
        files = preflight_rsem_reference(target_fasta, target_gtf)
    elif stage == "rrbs_collapse":
        if plan_root is None or bundle_id is None:
            raise GSE105127GzipPreflightError("RRBS preflight inputs are incomplete")
        files = preflight_rrbs(plan_root, bundle_id)
    elif stage == "rna_quantification":
        if plan_root is None or raw_root is None or bundle_id is None:
            raise GSE105127GzipPreflightError("RNA preflight inputs are incomplete")
        files = preflight_rna(plan_root, raw_root, bundle_id)
    elif stage == "cpg_crosswalk":
        if plan_root is None or reference_root is None or collapsed_root is None:
            raise GSE105127GzipPreflightError("crosswalk preflight inputs are incomplete")
        files = preflight_crosswalk(plan_root, reference_root, collapsed_root)
    else:
        raise GSE105127GzipPreflightError(f"unsupported preflight stage: {stage}")
    output.mkdir(mode=0o750)
    receipt = {
        "schema_version": "masld-bench-gse105127-gzip-preflight-v2",
        "status": "passed_exactly_one_gzip_member_per_input",
        "stage": stage,
        "bundle_id": bundle_id,
        "files": files,
        "file_count": len(files),
        "all_gzip_test_passed": True,
        "all_exactly_one_member": True,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "molecular_values_parsed": False,
        "labels_accessed": False,
        "outcomes_accessed": False,
        "fit_or_score_performed": False,
    }
    write_json_exclusive(output / "receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse105127_gzip_preflight_v2",
            "stage": stage,
            "bundle_id": bundle_id,
            "molecular_values_parsed": False,
            "status": "passed",
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        required=True,
        choices=("rsem_reference", "rrbs_collapse", "rna_quantification", "cpg_crosswalk"),
    )
    parser.add_argument("--plan-root", type=Path)
    parser.add_argument("--target-fasta", type=Path)
    parser.add_argument("--target-gtf", type=Path)
    parser.add_argument("--raw-root", type=Path)
    parser.add_argument("--reference-root", type=Path)
    parser.add_argument("--collapsed-root", type=Path)
    parser.add_argument("--bundle-id", type=int)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(**vars(args))
    print(json.dumps({"output": args.output.as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
