#!/usr/bin/env python3
"""Freeze one GSE238219 SRA experiment/run ancestry record."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


EXPERIMENTS = {
    "SRX21157431": ("GSM7660623", "SAMN36706212"),
    "SRX21157432": ("GSM7660624", "SAMN36706211"),
    "SRX21157433": ("GSM7660625", "SAMN36706210"),
    "SRX21157434": ("GSM7660626", "SAMN36706209"),
    "SRX21157435": ("GSM7660627", "SAMN36706208"),
}


class SraAuditError(RuntimeError):
    """Raised when SRA metadata differs from the frozen sample mapping."""


def digest(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def fetch_xml(accession: str, delay_seconds: int) -> tuple[str, bytes]:
    if delay_seconds:
        time.sleep(delay_seconds)
    query = urllib.parse.urlencode(
        {
            "db": "sra",
            "id": accession,
            "rettype": "full",
            "retmode": "xml",
        }
    )
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "MASLD-Model-Benchmark/1.0 sra-ancestry-audit"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:  # noqa: S310
        payload = response.read(5 * 1024 * 1024 + 1)
    if len(payload) > 5 * 1024 * 1024:
        raise SraAuditError("SRA XML exceeds 5 MiB ceiling")
    return url, payload


def values(root: ET.Element, tag: str) -> list[str]:
    result = []
    for element in root.iter():
        if local_name(element.tag) == tag and element.text and element.text.strip():
            result.append(element.text.strip())
    return sorted(set(result))


def records(root: ET.Element, tag: str) -> list[dict[str, str]]:
    return [
        dict(sorted(element.attrib.items()))
        for element in root.iter()
        if local_name(element.tag) == tag
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", choices=tuple(EXPERIMENTS), required=True)
    parser.add_argument("--delay-seconds", type=int, choices=range(0, 10), default=0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SraAuditError("refusing to overwrite SRA ancestry audit")
    sample, expected_biosample = EXPERIMENTS[args.experiment]
    url, payload = fetch_xml(args.experiment, args.delay_seconds)
    root = ET.fromstring(payload)
    accessions = values(root, "PRIMARY_ID") + values(root, "EXTERNAL_ID")
    if args.experiment not in accessions or expected_biosample not in accessions:
        raise SraAuditError("SRA experiment or BioSample mapping differs")
    summary = {
        "schema_version": "masld-bench-gse238219-sra-ancestry-audit-v1",
        "geo_sample": sample,
        "sra_experiment": args.experiment,
        "expected_biosample": expected_biosample,
        "source_url": url,
        "source_sha256": digest(payload),
        "source_size_bytes": len(payload),
        "identifier_values": sorted(set(accessions)),
        "experiment_records": records(root, "EXPERIMENT"),
        "sample_records": records(root, "SAMPLE"),
        "run_records": records(root, "RUN"),
        "library_strategy": values(root, "LIBRARY_STRATEGY"),
        "library_source": values(root, "LIBRARY_SOURCE"),
        "library_selection": values(root, "LIBRARY_SELECTION"),
        "instrument_model": values(root, "INSTRUMENT_MODEL"),
        "titles": values(root, "TITLE"),
        "biological_unit_status": "SRA_library_and_run_topology_only_not_independent_culture_proof",
        "gse313774_accessed": False,
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "sra_record.xml").write_bytes(payload)
    (args.output / "sra_ancestry.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
