#!/usr/bin/env python3
"""Acquire and exhaustively validate the 57 processed GSE105127 RRBS BEDs."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import time
import urllib.parse
import urllib.request


EXPECTED_JOIN_ARTIFACTS_SHA256 = (
    "574eeee1119054a2d7b55066804accf29fb45f23b591dd26f9d61d907b27a75f"
)
MAX_FILE_BYTES = 512 * 1024 * 1024
MAX_TOTAL_BYTES = 32 * 1024 * 1024 * 1024
DOWNLOAD_WORKERS = 4
AUDIT_WORKERS = 4
AWK_AUDIT = r"""
BEGIN { OFS="\t" }
/^track / { track_rows++; next }
{
    rows++
    if (NF != 11) { bad_width++; next }
    if ($1 !~ /^chr[[:alnum:]_.-]+$/) bad_chrom++
    if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $7 !~ /^[0-9]+$/ || $8 !~ /^[0-9]+$/) bad_coordinate++
    if ($4 !~ /^[0-9]+([.][0-9]+)?$/) bad_percent++
    if ($5 !~ /^[0-9]+$/ || $11 !~ /^[0-9]+$/) bad_coverage++
    start=$2+0; end=$3+0; pct=$4+0; cov=$5+0; other_cov=$11+0
    if (start < 0 || end != start + 1 || $7+0 != start || $8+0 != end) bad_interval++
    if (pct < 0 || pct > 100) bad_percent_range++
    if (cov < 1 || other_cov < 0) bad_coverage_range++
    if ($6 != "+" && $6 != "-") bad_strand++
    if ($9 !~ /^[0-9]+,[0-9]+,[0-9]+$/) bad_rgb++
    if ($10 !~ /^[0-9]+$/) bad_detail++
    detail=$10+0
    implied=pct*cov/100.0
    rounded=int(implied+0.5)
    residual=implied-rounded
    if (residual < 0) residual=-residual
    tolerance=cov*0.00005+0.00001
    if (residual > tolerance) bad_implied_methylated_count++
    if (($1 in last_start) && start < last_start[$1]) bad_sort++
    last_start[$1]=start
    chrom[$1]++
    plus += ($6 == "+")
    minus += ($6 == "-")
    sum_cov += cov
    sum_other_cov += other_cov
    sum_detail += detail
    detail_nonzero += (detail > 0)
    if (rows == 1 || detail > max_detail) max_detail=detail
    sum_weighted_methylated += implied
    ge5 += (cov >= 5)
    ge10 += (cov >= 10)
    ge20 += (cov >= 20)
    ge50 += (cov >= 50)
    ge100 += (cov >= 100)
    if (rows == 1 || cov < min_cov) min_cov=cov
    if (rows == 1 || cov > max_cov) max_cov=cov
}
END {
    error_total=bad_width+bad_chrom+bad_coordinate+bad_percent+bad_coverage+bad_interval+bad_percent_range+bad_coverage_range+bad_strand+bad_rgb+bad_detail+bad_implied_methylated_count+bad_sort
    printf "rows\t%d\n", rows
    printf "track_rows\t%d\n", track_rows
    printf "chromosomes\t%d\n", length(chrom)
    printf "plus_rows\t%d\n", plus
    printf "minus_rows\t%d\n", minus
    printf "coverage_sum\t%.0f\n", sum_cov
    printf "other_strand_coverage_sum\t%.0f\n", sum_other_cov
    printf "weighted_methylated_sum\t%.6f\n", sum_weighted_methylated
    printf "detail_sum\t%.0f\n", sum_detail
    printf "detail_nonzero_rows\t%d\n", detail_nonzero
    printf "max_detail\t%d\n", max_detail
    printf "coverage_ge5\t%d\n", ge5
    printf "coverage_ge10\t%d\n", ge10
    printf "coverage_ge20\t%d\n", ge20
    printf "coverage_ge50\t%d\n", ge50
    printf "coverage_ge100\t%d\n", ge100
    printf "min_coverage\t%d\n", min_cov
    printf "max_coverage\t%d\n", max_cov
    printf "bad_width\t%d\n", bad_width
    printf "bad_chrom\t%d\n", bad_chrom
    printf "bad_coordinate\t%d\n", bad_coordinate
    printf "bad_percent\t%d\n", bad_percent
    printf "bad_coverage\t%d\n", bad_coverage
    printf "bad_interval\t%d\n", bad_interval
    printf "bad_percent_range\t%d\n", bad_percent_range
    printf "bad_coverage_range\t%d\n", bad_coverage_range
    printf "bad_strand\t%d\n", bad_strand
    printf "bad_rgb\t%d\n", bad_rgb
    printf "bad_detail\t%d\n", bad_detail
    printf "bad_implied_methylated_count\t%d\n", bad_implied_methylated_count
    printf "bad_sort\t%d\n", bad_sort
    printf "error_total\t%d\n", error_total
}
"""


class RRBSAuditError(RuntimeError):
    """Raised when an RRBS source or assay-native invariant differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rrbs_rows(join: Path) -> list[dict[str, str]]:
    with (join / "sample_join.tsv").open(encoding="utf-8", newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["assay"] == "RRBS"
        ]
    if len(rows) != 57 or len({row["sample_accession"] for row in rows}) != 57:
        raise RRBSAuditError("authoritative join does not contain 57 unique RRBS samples")
    if any(row["pairing"] != "adjacent_section" for row in rows):
        raise RRBSAuditError("RRBS join no longer has adjacent-section topology")
    return rows


def https_url(value: str) -> str:
    prefix = "ftp://ftp.ncbi.nlm.nih.gov/"
    if not value.startswith(prefix):
        raise RRBSAuditError("processed RRBS URL is not an NCBI FTP source")
    return "https://ftp.ncbi.nlm.nih.gov/" + value[len(prefix) :]


def download_one(row: dict[str, str], raw: Path) -> dict[str, object]:
    url = https_url(row["primary_processed_file"])
    source_name = Path(urllib.parse.urlparse(url).path).name
    target = raw / f"{row['sample_accession']}__{source_name}"
    temporary = target.with_suffix(target.suffix + ".part")
    if target.exists() or temporary.exists():
        raise RRBSAuditError(f"refusing to overwrite download: {target}")
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "masld-bench-rrbs-acquisition/1.0"}
            )
            digest = sha256()
            byte_count = 0
            with urllib.request.urlopen(request, timeout=300) as response:
                declared = response.headers.get("Content-Length")
                with temporary.open("xb") as handle:
                    while True:
                        block = response.read(8 * 1024 * 1024)
                        if not block:
                            break
                        byte_count += len(block)
                        if byte_count > MAX_FILE_BYTES:
                            raise RRBSAuditError(f"RRBS file exceeds size cap: {url}")
                        digest.update(block)
                        handle.write(block)
                if declared is not None and int(declared) != byte_count:
                    raise RRBSAuditError(
                        f"Content-Length differs for {url}: {declared} != {byte_count}"
                    )
            os.replace(temporary, target)
            return {
                **row,
                "source_url": url,
                "local_file": str(target.relative_to(raw.parent)),
                "bytes": byte_count,
                "sha256": digest.hexdigest(),
                "download_attempts": attempt,
            }
        except Exception as error:  # retry bounded network failures, then fail closed
            last_error = error
            temporary.unlink(missing_ok=True)
            if attempt < 3:
                time.sleep(attempt * 2)
    raise RRBSAuditError(f"download failed after three attempts: {url}") from last_error


def reuse_one(row: dict[str, str], source: Path, raw: Path) -> dict[str, object]:
    url = https_url(row["primary_processed_file"])
    source_name = Path(urllib.parse.urlparse(url).path).name
    candidates = list(source.glob(f"{row['sample_accession']}__{source_name}"))
    if len(candidates) != 1:
        raise RRBSAuditError(f"reusable raw file census differs: {row['sample_accession']}")
    origin = candidates[0].resolve(strict=True)
    target = raw / origin.name
    if target.exists():
        raise RRBSAuditError(f"refusing to overwrite reusable raw file: {target}")
    os.link(origin, target)
    byte_count = target.stat().st_size
    request = urllib.request.Request(
        url,
        method="HEAD",
        headers={"User-Agent": "masld-bench-rrbs-reuse-audit/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        declared = response.headers.get("Content-Length")
    if declared is not None and int(declared) != byte_count:
        raise RRBSAuditError(f"reused file size differs from NCBI source: {url}")
    return {
        **row,
        "source_url": url,
        "local_file": str(target.relative_to(raw.parent)),
        "bytes": byte_count,
        "sha256": sha256_file(target),
        "download_attempts": 0,
        "reused_from_failed_attempt": str(origin),
    }


def parse_awk_output(payload: str) -> dict[str, int | float]:
    observed: dict[str, int | float] = {}
    float_keys = {"weighted_methylated_sum"}
    for line in payload.splitlines():
        key, value = line.split("\t", 1)
        observed[key] = float(value) if key in float_keys else int(value)
    return observed


def audit_one(receipt: dict[str, object], root: Path) -> dict[str, object]:
    path = root / str(receipt["local_file"])
    decompressor = subprocess.Popen(
        ["gzip", "-dc", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env={**os.environ, "LC_ALL": "C"},
    )
    if decompressor.stdout is None:
        raise RRBSAuditError("gzip stdout pipe was not created")
    audit = subprocess.run(
        ["awk", AWK_AUDIT],
        stdin=decompressor.stdout,
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "LC_ALL": "C"},
    )
    decompressor.stdout.close()
    gzip_stderr = decompressor.stderr.read().decode("utf-8", errors="replace") if decompressor.stderr else ""
    gzip_status = decompressor.wait()
    if gzip_status != 0:
        raise RRBSAuditError(f"gzip integrity/read failure for {path}: {gzip_stderr}")
    if audit.returncode != 0:
        raise RRBSAuditError(f"awk assay audit failed for {path}: {audit.stderr}")
    metrics = parse_awk_output(audit.stdout)
    if metrics.get("track_rows") != 1:
        raise RRBSAuditError(f"expected one BED track header: {path}")
    if int(metrics.get("rows", 0)) <= 0 or metrics.get("error_total") != 0:
        raise RRBSAuditError(f"RRBS row invariants failed for {path}: {metrics}")
    return {**receipt, **metrics, "hard_qc_state": "observed"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--join", type=Path, required=True)
    parser.add_argument(
        "--join-artifacts-sha256", default=EXPECTED_JOIN_ARTIFACTS_SHA256
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reuse-raw", type=Path)
    args = parser.parse_args()
    if sha256_file(args.join / "ARTIFACTS.json") != args.join_artifacts_sha256:
        raise RRBSAuditError("join ARTIFACTS SHA-256 differs")
    rows = read_rrbs_rows(args.join)
    args.output.mkdir(parents=True, exist_ok=False)
    raw = args.output / "raw"
    raw.mkdir()

    receipts: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=DOWNLOAD_WORKERS) as executor:
        if args.reuse_raw is None:
            futures = {executor.submit(download_one, row, raw): row for row in rows}
        else:
            source = args.reuse_raw.resolve(strict=True)
            futures = {
                executor.submit(reuse_one, row, source, raw): row for row in rows
            }
        for future in as_completed(futures):
            receipts.append(future.result())
    receipts.sort(key=lambda row: str(row["sample_accession"]))
    total_bytes = sum(int(row["bytes"]) for row in receipts)
    if total_bytes > MAX_TOTAL_BYTES:
        raise RRBSAuditError("downloaded RRBS corpus exceeds total size cap")

    audits: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=AUDIT_WORKERS) as executor:
        futures = {
            executor.submit(audit_one, receipt, args.output): receipt
            for receipt in receipts
        }
        for future in as_completed(futures):
            audits.append(future.result())
    audits.sort(key=lambda row: str(row["sample_accession"]))

    fields = tuple(audits[0])
    with (args.output / "rrbs_file_qc.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(audits)
    summary = {
        "schema_version": "masld-bench-gse105127-rrbs-source-qc-v1",
        "status": "pass",
        "biological_unit": "participant",
        "participants": len({str(row["participant_id"]) for row in audits}),
        "rrbs_files": len(audits),
        "zones": dict(sorted(Counter(str(row["zone"]) for row in audits).items())),
        "phenotypes": dict(
            sorted(Counter(str(row["phenotype"]) for row in audits).items())
        ),
        "compressed_bytes": total_bytes,
        "cytosine_rows": sum(int(row["rows"]) for row in audits),
        "hard_qc_failures": sum(row["hard_qc_state"] != "observed" for row in audits),
        "native_build": "1000_Genomes_GRCh37_hg19",
        "bed_coordinate_evidence": "0_based_half_open_single_base_intervals_with_equal_thick_start_end",
        "row_schema": [
            "chrom",
            "start_0based",
            "end_0based_half_open",
            "methylation_percent",
            "strand_coverage",
            "strand",
            "thick_start",
            "thick_end",
            "item_rgb",
            "bed_detail_auxiliary_nonnegative_integer_semantics_unresolved",
            "other_strand_coverage",
        ],
        "methylation_likelihood_contract": "coverage_aware_binomial_or_beta_binomial_at_participant_level",
        "methylation_fraction_as_rna_count_forbidden": True,
        "strand_rows_require_cpg_pair_collapse_before_modeling": True,
        "coordinate_remap_complete": False,
        "model_training_activated": False,
        "remaining_blockers": [
            "freeze_exact_1000_Genomes_GRCh37_reference",
            "validate_reference_base_and_CpG_strand_pair_semantics",
            "resolve_or_exclude_bed_detail_auxiliary_integer_before_modeling",
            "freeze_failure_aware_hg19_to_GRCh38p14_CpG_crosswalk",
            "freeze_RNA_processing_and_task_specific_activation_contract",
        ],
    }
    (args.output / "qc_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (args.output / "download_receipt.json").write_text(
        json.dumps(receipts, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
