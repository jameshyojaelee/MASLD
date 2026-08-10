#!/usr/bin/env python3
"""Freeze an official array-annotation remediation after GEO lacks symbols."""

from __future__ import annotations

import datetime as dt
import json

from public_functional_common import CANDIDATE_ROOT, atomic_write_text, require_sealed, stable_json_sha256


def main() -> None:
    seal = require_sealed()
    destination = CANDIDATE_ROOT / "SOURCE_AMENDMENT_02.json"
    if destination.exists():
        raise RuntimeError(f"Amendment already exists: {destination}")
    error_log = CANDIDATE_ROOT / "logs/bulk_19655508_2.err"
    expected = "GPL16686 probe or gene-symbol column not found"
    if not error_log.is_file() or expected not in error_log.read_text(encoding="utf-8"):
        raise RuntimeError("The expected fail-closed GPL16686 mapping error was not reproduced")
    amendment = {
        "amendment_id": "SOURCE_AMENDMENT_02_GPL16686_BIOCONDUCTOR_ANNOTATION",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "parent_specification_sha256": seal["specification_sha256"],
        "reason": "The deposited GPL16686 table contains transcript-cluster IDs and genomic/accession fields but no gene-symbol column; the official Bioconductor transcript-cluster annotation resolves that exact platform mapping gate.",
        "package": "hugene20sttranscriptcluster.db",
        "version": "8.8.0",
        "url": "https://bioconductor.org/packages/release/data/annotation/src/contrib/hugene20sttranscriptcluster.db_8.8.0.tar.gz",
        "published_md5": "0b929a3a959662e8a7265f58b81b4e35",
        "scope": "probe-to-gene annotation only; no expression outcome, contrast, program, weight, class, multiplicity, or promotion rule changes",
        "mapping_rule": "discard probes mapping to zero or multiple gene symbols; among remaining probes select highest median expression per gene outcome-blindly; median-across-probes sensitivity remains required",
    }
    amendment["amendment_sha256"] = stable_json_sha256(amendment)
    atomic_write_text(destination, json.dumps(amendment, indent=2, sort_keys=True) + "\n")
    atomic_write_text(CANDIDATE_ROOT / "SOURCE_AMENDMENT_02_SHA256", amendment["amendment_sha256"] + "\n")
    print(json.dumps(amendment, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
