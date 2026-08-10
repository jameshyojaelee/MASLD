#!/usr/bin/env python3
"""Freeze a same-study source-gate remediation before inspecting the OS-HLO object."""

from __future__ import annotations

import datetime as dt
import json

from public_functional_common import CANDIDATE_ROOT, atomic_write_text, read_tsv, require_sealed, sha256_file, stable_json_sha256


def main() -> None:
    seal = require_sealed()
    destination = CANDIDATE_ROOT / "SOURCE_AMENDMENT_01.json"
    if destination.exists():
        raise RuntimeError(f"Amendment already exists: {destination}")
    gates = read_tsv(CANDIDATE_ROOT / "source_gate_status.tsv")
    row = [item for item in gates if item["dataset_id"] == "GSE207889"]
    if len(row) != 1 or row[0]["status"] != "rejected_qc":
        raise RuntimeError("Amendment is allowed only after the specified merged-H5AD source gate rejects")
    if "obs fields=barcode,condition,culture,doublet_score,sample" not in row[0]["detail"]:
        raise RuntimeError("Observed GSE207889 source failure differs from the reviewed failure")
    amendment = {
        "amendment_id": "SOURCE_AMENDMENT_01_GSE207889_OS_H5AD",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "parent_specification_sha256": seal["specification_sha256"],
        "failed_gate_table_sha256": sha256_file(CANDIDATE_ROOT / "source_gate_status.tsv"),
        "reason": "The preregistered merged H5AD contains source sample/condition and integer raw counts in X but no source cell-lineage annotation; the same GEO record deposits an OS-HLO-specific H5AD intended to resolve that exact lineage gate.",
        "dataset_id": "GSE207889",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207889/suppl/GSE207889_os.h5ad.gz",
        "local_name": "GSE207889_os.h5ad.gz",
        "scope": "source-schema remediation only; same study, treatments, replicates, hypotheses, contrasts, programs, weights, multiplicity, and promotion rules",
        "prohibitions": [
            "no program or evidence-class outcome may be read before the amended source gate",
            "no self-derived cell annotation may substitute for a missing deposited lineage field",
            "failure of the OS-HLO source gate terminates the scHLO branch",
        ],
    }
    amendment["amendment_sha256"] = stable_json_sha256(amendment)
    atomic_write_text(destination, json.dumps(amendment, indent=2, sort_keys=True) + "\n")
    atomic_write_text(CANDIDATE_ROOT / "SOURCE_AMENDMENT_01_SHA256", amendment["amendment_sha256"] + "\n")
    print(json.dumps(amendment, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
