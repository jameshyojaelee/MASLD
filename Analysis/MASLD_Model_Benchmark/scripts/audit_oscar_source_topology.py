#!/usr/bin/env python3
"""Audit OSCAR GSA capture topology against the source-paper methods."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import urllib.request


FULLTEXT_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10617096/fullTextXML"


class TopologyError(RuntimeError):
    """Raised when the OSCAR source topology differs."""


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gsa-html", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.gsa_html.is_file() or args.output.exists():
        raise TopologyError("OSCAR topology request differs")
    html = args.gsa_html.read_text(encoding="utf-8")
    blocks = re.findall(
        r'(<tr class="experiment">.*?)(?=<tr class="experiment">|\Z)', html, flags=re.DOTALL
    )
    experiments = []
    for block in blocks:
        match = re.search(
            r'browse/CRA009621/(CRX\d+)">[^<]+</a>\s*</td>\s*'
            r'<td>([^<]+)</td>.*?biosample/browse/(SAMC\d+)',
            block,
            flags=re.DOTALL,
        )
        if match is None:
            raise TopologyError("OSCAR GSA experiment record differs")
        accession, name, biosample = match.groups()
        run_accessions = re.findall(r'browse/CRA009621/(CRR\d+)">', block)
        run_aliases = re.findall(r'<td colspan="2">(OSCAR_[^<]+)</td>', block)
        files = re.findall(r'<strong>File: </strong>(CRR\d+_[^<\s]+)', block)
        if len(run_accessions) != len(set(run_accessions)) or len(run_aliases) != len(run_accessions):
            raise TopologyError(f"OSCAR GSA run topology differs for {accession}")
        if len(files) != 2 * len(run_accessions):
            raise TopologyError(f"OSCAR GSA file topology differs for {accession}")
        for run in run_accessions:
            expected = {f"{run}_f1.fastq.gz", f"{run}_r2.fastq.gz"}
            if set(file for file in files if file.startswith(f"{run}_")) != expected:
                raise TopologyError(f"OSCAR paired FASTQ topology differs for {run}")
        name_match = re.fullmatch(r"Organoid_OSCAR_(EM|DM)_(\d+|pilot)", name)
        if name_match is None:
            raise TopologyError(f"OSCAR experiment name differs: {name}")
        context, lane = name_match.groups()
        experiments.append({
            "experiment": accession,
            "name": name,
            "biosample": biosample,
            "context": context,
            "capture_lane": lane,
            "production": lane != "pilot",
            "run_count": len(run_accessions),
            "run_accessions": run_accessions,
            "run_aliases": run_aliases,
            "fastq_file_count": len(files),
        })
    production = [record for record in experiments if record["production"]]
    pilot = [record for record in experiments if not record["production"]]
    production_by_context = {
        context: sum(record["context"] == context for record in production)
        for context in ("DM", "EM")
    }
    if (
        len(experiments) != 11
        or len(production) != 9
        or len(pilot) != 2
        or production_by_context != {"DM": 6, "EM": 3}
        or len({record["biosample"] for record in experiments}) != 11
    ):
        raise TopologyError("OSCAR GSA capture-lane census differs")
    args.output.mkdir(parents=True, exist_ok=False)
    fulltext_path = args.output / "PMC10617096_fulltext.xml"
    request = urllib.request.Request(FULLTEXT_URL, headers={"User-Agent": "MASLD-model-benchmark/1"})
    with urllib.request.urlopen(request, timeout=120) as response:
        if response.status != 200 or response.geturl() != FULLTEXT_URL:
            raise TopologyError("OSCAR full-text response differs")
        fulltext_path.write_bytes(response.read())
    xml_text = fulltext_path.read_text(encoding="utf-8")
    fulltext = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", xml_text))
    required_methods = (
        "500,000 cell clusters were infected",
        "30,000 BECs from EM condition",
        "60,000 differentiated hepatocyte-like cells from DM condition",
        "10,000 input cells for each lane",
        "Cells with unique sgRNA were retained for downstream analysis",
    )
    if any(phrase not in fulltext for phrase in required_methods):
        raise TopologyError("OSCAR source-paper method evidence differs")
    if "creativecommons.org/licenses/by/4.0" not in xml_text:
        raise TopologyError("OSCAR article license differs")
    result = {
        "schema_version": "masld-bench-oscar-source-topology-v1",
        "dataset_id": "cra009621_oscar",
        "gsa_html_sha256": digest(args.gsa_html),
        "fulltext_url": FULLTEXT_URL,
        "fulltext_sha256": digest(fulltext_path),
        "article_license": "CC_BY_4.0",
        "experiment_count": len(experiments),
        "production_capture_lane_count": len(production),
        "pilot_capture_lane_count": len(pilot),
        "production_capture_lanes_by_context": production_by_context,
        "registered_run_count": sum(record["run_count"] for record in experiments),
        "registered_fastq_file_count": sum(record["fastq_file_count"] for record in experiments),
        "distinct_biosample_count": len({record["biosample"] for record in experiments}),
        "experiments": experiments,
        "method_evidence": {
            "upstream_screen": "The paper describes one 500,000-cluster infection workflow followed by FACS, expansion, passage, and context exposure; it does not report independent screen cultures or mouse-of-origin ancestry.",
            "capture": "The paper describes approximately 30,000 EM and 60,000 DM cells captured with 10,000 input cells per Chromium lane, matching the three EM and six DM production experiments in GSA.",
            "filter": "The deposited author metadata contains cells retained after unique-sgRNA assignment; this is a source preprocessing gate, not a biological replicate definition."
        },
        "admission_conclusion": "The nine production experiments and BioSamples are capture-lane records, not nine biological replicates. Independent mouse-of-origin and culture ancestry remain unresolved; at most one shared screen pool is admitted for descriptive use.",
        "outcome_model_fit": False,
        "gse313774_accessed": False,
    }
    (args.output / "source_topology.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
