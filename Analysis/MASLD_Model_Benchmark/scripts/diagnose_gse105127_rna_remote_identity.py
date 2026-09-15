#!/usr/bin/env python3
"""Census GSE105127 RNA remote byte identities without downloading FASTQs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import time
import urllib.error
import urllib.request

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


class RemoteIdentityError(RuntimeError):
    """Raised when the label-free plan or completed local output file differs."""


def parse_content_range(value: str | None) -> int | None:
    if not value or "/" not in value:
        return None
    total = value.rsplit("/", 1)[1]
    return int(total) if total.isdigit() else None


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        forbidden = {"phenotype", "label", "outcome", "disease", "fibrosis", "nas", "sex", "age", "bmi", "outer_fold"}
        if forbidden & fields:
            raise RemoteIdentityError("RNA remote census violates the label firewall")
        rows = list(reader)
    if len(rows) != 57 or len({row["row_id"] for row in rows}) != 57:
        raise RemoteIdentityError("RNA plan row identity differs")
    return rows


def remote_total(url: str) -> dict[str, object]:
    errors = []
    for attempt in range(1, 4):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "masld-bench-rna-identity/1.0"},
                method="HEAD",
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                length = response.headers.get("Content-Length")
                return {
                    "status": "head_complete",
                    "http_status": response.status,
                    "content_length": int(length) if length and length.isdigit() else None,
                    "content_range_total": parse_content_range(response.headers.get("Content-Range")),
                    "etag": response.headers.get("ETag"),
                    "last_modified": response.headers.get("Last-Modified"),
                    "final_url": response.geturl(),
                    "attempts": attempt,
                }
        except (OSError, urllib.error.URLError, ValueError) as error:
            errors.append(f"{type(error).__name__}: {error}"[:500])
            if attempt < 3:
                time.sleep(attempt * 2)
    return {
        "status": "head_unavailable",
        "http_status": None,
        "content_length": None,
        "content_range_total": None,
        "etag": None,
        "last_modified": None,
        "final_url": None,
        "attempts": 3,
        "errors": errors,
    }


def run(*, plan_root: Path, local_root: Path, output: Path) -> dict[str, object]:
    if output.exists():
        raise RemoteIdentityError(f"refusing to overwrite remote census: {output}")
    verify_frozen_tree(plan_root)
    rows = read_rows(plan_root / "rna_rows.tsv")
    output.mkdir(parents=True)
    results = []
    for row in rows:
        local = local_root / "participants" / row["row_id"]
        local_complete = (local / "COMPLETE").is_file()
        local_receipt = None
        if local_complete:
            verify_frozen_tree(local)
            local_receipt = json.loads((local / "receipt.json").read_text(encoding="utf-8"))
            if (
                local_receipt.get("row_id") != row["row_id"]
                or local_receipt.get("bytes") != int(row["fastq_bytes"])
                or local_receipt.get("md5") != row["fastq_md5"].lower()
            ):
                raise RemoteIdentityError("completed RNA member differs from plan")
        remote = remote_total(row["fastq_url"])
        observed_total = remote.get("content_length") or remote.get("content_range_total")
        results.append({
            "bundle_id": int(row["bundle_id"]),
            "row_id": row["row_id"],
            "run_accession": row["run_accession"],
            "planned_bytes": int(row["fastq_bytes"]),
            "planned_md5": row["fastq_md5"].lower(),
            "local_complete": local_complete,
            "local_bytes": local_receipt.get("bytes") if local_receipt else None,
            "remote": remote,
            "remote_total_bytes": observed_total,
            "remote_bytes_match_plan": observed_total == int(row["fastq_bytes"]),
        })
    receipt = {
        "schema_version": "masld-bench-gse105127-rna-remote-identity-census-v1",
        "status": "complete_outcome_free_remote_identity_census",
        "planned_rows": len(results),
        "local_complete_rows": sum(row["local_complete"] for row in results),
        "remote_head_complete_rows": sum(row["remote"]["status"] == "head_complete" for row in results),
        "remote_bytes_match_rows": sum(row["remote_bytes_match_plan"] for row in results),
        "remote_bytes_differ_rows": sum(
            row["remote_total_bytes"] is not None and not row["remote_bytes_match_plan"]
            for row in results
        ),
        "remote_total_unavailable_rows": sum(row["remote_total_bytes"] is None for row in results),
        "labels_accessed": False,
        "outcomes_accessed": False,
        "fit_or_score_performed": False,
    }
    write_json_exclusive(output / "row_identity.json", {"rows": results})
    write_json_exclusive(output / "receipt.json", receipt)
    freeze_tree(output, {
        "artifact_class": "gse105127_rna_remote_identity_census",
        "labels_accessed": False,
        "outcomes_accessed": False,
        "fit_or_score_performed": False,
        "status": "complete",
    })
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--local-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(
        plan_root=args.plan_root,
        local_root=args.local_root,
        output=args.output,
    ), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
