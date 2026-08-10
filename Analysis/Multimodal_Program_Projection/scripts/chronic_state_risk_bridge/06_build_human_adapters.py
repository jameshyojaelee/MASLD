#!/usr/bin/env python3
"""Build outcome-blind human expression adapters from deposited metadata."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

from bridge_common import CANDIDATE_ROOT, PROJECT_ROOT, require_validated_seal, write_tsv
from importlib.machinery import SourceFileLoader


PLAN41_ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09"
PAIR_GATE = SourceFileLoader(
    "pair_gate",
    str(Path(__file__).with_name("03_gate_human_pairs.py")),
).load_module()


def sample_rows(dataset: str, matrix_name: str) -> list[dict[str, object]]:
    sample, fields, _ = PAIR_GATE.matrix_metadata(CANDIDATE_ROOT / f"sources/{dataset}/{matrix_name}")
    accessions = sample["!Sample_geo_accession"]
    titles = sample["!Sample_title"]
    rows: list[dict[str, object]] = []
    for index, accession in enumerate(accessions):
        if dataset == "GSE83452":
            timepoint = fields["time"][index]
            intervention = fields["type of intervention"][index]
            histology = fields["liver status"][index]
            nas = ""
            bmi = ""
        else:
            timepoint = fields["bariatric surgery"][index]
            intervention = "bariatric_surgery"
            histology = fields.get("diagnosis", [""] * len(accessions))[index]
            nas = fields["nas"][index]
            bmi = fields["bmi"][index]
        rows.append({
            "dataset_id": dataset,
            "sample_id": accession,
            "sample_title": titles[index],
            "timepoint": timepoint,
            "intervention": intervention,
            "histology": histology,
            "nas": nas,
            "bmi": bmi,
            "reference_eligible": str(timepoint in {"baseline", "before surgery"}).lower(),
            "source_basis": "deposited GEO series-matrix sample metadata",
        })
    return rows


def gse106737_rows() -> list[dict[str, object]]:
    path = PLAN41_ROOT / "analyses/GSE106737/analysis_sample_manifest.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        source = list(csv.DictReader(handle, delimiter="\t"))
    return [{
        "dataset_id": "GSE106737",
        "sample_id": row["sample_id"],
        "sample_title": row["sample_id"],
        "timepoint": row["timepoint"],
        "intervention": row["treatment"],
        "histology": row["condition"],
        "nas": "",
        "bmi": "",
        "reference_eligible": str(row["timepoint"] == "baseline").lower(),
        "source_basis": "terminal validated Plan 41 sample manifest",
    } for row in source]


def gpl11532_mapping() -> list[dict[str, object]]:
    path = CANDIDATE_ROOT / "sources/GSE48452/GSE48452_family.soft.gz"
    rows: list[dict[str, object]] = []
    active = False
    header: list[str] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                active = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not active:
                continue
            values = line.rstrip("\n").split("\t")
            if header is None:
                header = values
                continue
            record = dict(zip(header, values))
            symbols = set()
            for assignment in record.get("gene_assignment", "").split(" /// "):
                fields = [part.strip() for part in assignment.split(" // ")]
                if len(fields) > 1 and fields[1] not in {"", "---"}:
                    symbols.add(fields[1])
            rows.append({
                "probe_id": record["ID"],
                "symbol": next(iter(symbols)) if len(symbols) == 1 else "",
                "mapping_status": "unambiguous_single_symbol" if len(symbols) == 1 else "ambiguous_or_unmapped",
                "n_symbols": len(symbols),
            })
    if len(rows) != 33297:
        raise RuntimeError(f"Unexpected GPL11532 row count: {len(rows)}")
    return rows


def main() -> None:
    require_validated_seal()
    human_rows = (
        gse106737_rows()
        + sample_rows("GSE83452", "GSE83452_series_matrix.txt.gz")
        + sample_rows("GSE48452", "GSE48452_series_matrix.txt.gz")
    )
    expected = {"GSE106737": 111, "GSE83452": 231, "GSE48452": 73}
    observed = {dataset: sum(row["dataset_id"] == dataset for row in human_rows) for dataset in expected}
    if observed != expected:
        raise RuntimeError(f"Human adapter sample-count drift: {observed}")
    write_tsv(
        CANDIDATE_ROOT / "human_reversal/human_sample_manifest.tsv",
        human_rows,
        list(human_rows[0]),
    )
    mapping = gpl11532_mapping()
    write_tsv(
        CANDIDATE_ROOT / "human_reversal/GPL11532_probe_symbol.tsv",
        mapping,
        list(mapping[0]),
    )
    print("HUMAN_ADAPTERS_COMPLETE")


if __name__ == "__main__":
    main()
