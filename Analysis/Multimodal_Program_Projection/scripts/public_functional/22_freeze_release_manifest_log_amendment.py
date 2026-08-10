#!/usr/bin/env python3
"""Freeze the post-job manifest audit before excluding mutable SLURM logs."""

from __future__ import annotations

import datetime as dt
import json
import shutil

from public_functional_common import CANDIDATE_ROOT, atomic_write_text, sha256_file


EXPECTED_RELEASE_SHA256 = "22a2cf27b0009954f92e1e7a5599232da0681a9ddddb7f3ed88d821507f4393a"
ARCHIVE = CANDIDATE_ROOT / "superseded_bundles/2026-08-09_pre_log_exclusion_fix"


def main() -> None:
    validated_path = CANDIDATE_ROOT / "VALIDATED.json"
    release_path = CANDIDATE_ROOT / "release_manifest.tsv"
    validated = json.loads(validated_path.read_text())
    if validated.get("release_manifest_sha256") != EXPECTED_RELEASE_SHA256:
        raise RuntimeError("Unexpected terminal release identity")
    if sha256_file(release_path) != EXPECTED_RELEASE_SHA256:
        raise RuntimeError("release_manifest.tsv has drifted before amendment")
    if ARCHIVE.exists():
        raise RuntimeError(f"Archive already exists: {ARCHIVE}")

    files = [
        CANDIDATE_ROOT / "VALIDATED.json",
        CANDIDATE_ROOT / "release_manifest.tsv",
        CANDIDATE_ROOT / "validation_report.tsv",
        CANDIDATE_ROOT / "logs/final_19655697.out",
        CANDIDATE_ROOT / "logs/final_19655697.err",
    ]
    rows = []
    for source in files:
        if not source.is_file():
            raise RuntimeError(f"Missing pre-amendment artifact: {source}")
        destination = ARCHIVE / source.relative_to(CANDIDATE_ROOT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        rows.append(
            {
                "source": str(source.relative_to(CANDIDATE_ROOT)),
                "archive": str(destination.relative_to(CANDIDATE_ROOT)),
                "size_bytes": source.stat().st_size,
                "sha256": sha256_file(source),
            }
        )
    lines = ["source\tarchive\tsize_bytes\tsha256"]
    lines.extend("\t".join(str(row[key]) for key in ["source", "archive", "size_bytes", "sha256"]) for row in rows)
    manifest = ARCHIVE / "archive_manifest.tsv"
    atomic_write_text(manifest, "\n".join(lines) + "\n")

    amendment = {
        "amendment_id": "release-manifest-mutable-log-exclusion-v1",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "superseded_release_manifest_sha256": EXPECTED_RELEASE_SHA256,
        "archive_manifest_sha256": sha256_file(manifest),
        "observed_failure": "post-job independent audit found final_19655697.out changed after the in-job release manifest hashed it",
        "scientific_artifacts_affected": False,
        "required_correction": "exclude the entire mutable logs directory from the checksummed scientific payload and rerun terminal validation without rerunning scientific analyses",
        "allowed_changes": [
            "13_validate_final.py manifest exclusion rule",
            "new terminal validation freeze and wrapper",
            "release/validation documentation",
        ],
        "prohibited_changes": [
            "scientific effects, gates, contrasts, programs, evidence classes, or verdict logic",
            "mutation or deletion of prior bundles",
        ],
    }
    destination = CANDIDATE_ROOT / "RELEASE_MANIFEST_LOG_AMENDMENT_01.json"
    atomic_write_text(destination, json.dumps(amendment, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "frozen", "amendment_sha256": sha256_file(destination), "archived_files": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
