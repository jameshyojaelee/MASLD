#!/usr/bin/env python3
"""Acquire exact allowlisted training-only gene-perturbation inputs."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any
import urllib.request


SOURCES: dict[str, list[dict[str, Any]]] = {
    "gse264667_hepg2": [
        {
            "name": "GSE264667_hepg2_raw_singlecell_01.h5ad",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_hepg2_raw_singlecell_01.h5ad",
            "expected_size_bytes": 5614460941,
        }
    ],
    "gse242934_perturbseq": [
        {
            "name": "GSM7775285_RNA_756_barcodes.tsv.gz",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7775nnn/GSM7775285/suppl/GSM7775285_RNA_756_barcodes.tsv.gz",
            "expected_size_bytes": 30410,
        },
        {
            "name": "GSM7775285_RNA_756_features.tsv.gz",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7775nnn/GSM7775285/suppl/GSM7775285_RNA_756_features.tsv.gz",
            "expected_size_bytes": 232328,
        },
        {
            "name": "GSM7775285_RNA_756_matrix.mtx.gz",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7775nnn/GSM7775285/suppl/GSM7775285_RNA_756_matrix.mtx.gz",
            "expected_size_bytes": 112209964,
        },
        {
            "name": "GSM7775286_BARCODE_10x_RNA-756_RNA-757.txt.gz",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7775nnn/GSM7775286/suppl/GSM7775286_BARCODE_10x_RNA-756_RNA-757.txt.gz",
            "expected_size_bytes": 159905,
        },
    ],
}


class AcquisitionError(RuntimeError):
    """Raised when one exact source file differs."""


def acquire(spec: dict[str, Any], target: Path) -> dict[str, Any]:
    request = urllib.request.Request(
        spec["url"],
        headers={"User-Agent": "MASLD-Model-Benchmark/1.0 training-input-acquisition"},
    )
    value = sha256()
    observed = 0
    with urllib.request.urlopen(request, timeout=180) as response, target.open("xb") as output:  # noqa: S310
        declared = response.headers.get("Content-Length")
        if declared is not None and int(declared) != spec["expected_size_bytes"]:
            raise AcquisitionError(f"source Content-Length differs for {spec['name']}")
        while block := response.read(8 * 1024 * 1024):
            output.write(block)
            value.update(block)
            observed += len(block)
    if observed != spec["expected_size_bytes"]:
        raise AcquisitionError(f"downloaded size differs for {spec['name']}")
    return {
        **spec,
        "observed_size_bytes": observed,
        "sha256": value.hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=tuple(SOURCES), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise AcquisitionError("refusing to overwrite training-input acquisition")
    args.output.mkdir(parents=True, exist_ok=False)
    raw = args.output / "raw"
    raw.mkdir()
    records = [acquire(spec, raw / spec["name"]) for spec in SOURCES[args.source]]
    result = {
        "schema_version": "masld-bench-gene-perturbation-training-input-acquisition-v1",
        "source": args.source,
        "records": records,
        "allowlist_exact": True,
        "outcome_model_fit": False,
        "biological_replication_claim": False,
        "gse313774_accessed": False,
    }
    (args.output / "acquisition.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
