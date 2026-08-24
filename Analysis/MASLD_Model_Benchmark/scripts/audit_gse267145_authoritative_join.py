#!/usr/bin/env python3
"""Triangulate and freeze the GSE267145 RNA-H3K27ac participant join."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re
import shlex
import urllib.request
import xml.etree.ElementTree as ET


SAMPLE_FIELDS = (
    "series",
    "accession",
    "title",
    "source_name",
    "candidate_participant_id",
    "characteristics_json",
    "data_processing_json",
    "supplementary_files_json",
)
STAGES = ("NOR", "NAFL", "NASH_F0", "NASH_F1", "NASH_F23")


class GSE267145JoinError(RuntimeError):
    """Raised when one line of participant-join evidence is inconsistent."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_samples(path: Path, expected_series: str) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != SAMPLE_FIELDS:
            raise GSE267145JoinError(f"sample fields differ: {path}")
        rows = [dict(row) for row in reader]
    if not rows or any(row["series"] != expected_series for row in rows):
        raise GSE267145JoinError(f"sample series differs: {path}")
    return rows


def parse_characteristics(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in json.loads(value):
        if ":" not in item:
            raise GSE267145JoinError(f"malformed characteristic: {item}")
        key, observed = item.split(":", 1)
        key, observed = key.strip(), observed.strip()
        if key in result:
            raise GSE267145JoinError(f"duplicate characteristic: {key}")
        result[key] = observed
    return result


def participant_fold(participant_id: str, stage: str, stage_rank: int) -> int:
    """Return a deterministic stage-balanced five-fold assignment."""

    if stage not in STAGES or stage_rank < 0:
        raise GSE267145JoinError("invalid participant fold input")
    return stage_rank % 5


def build_join_rows(
    h3_rows: list[dict[str, str]], rna_rows: list[dict[str, str]]
) -> list[dict[str, str]]:
    rna_by_id = {row["candidate_participant_id"]: row for row in rna_rows}
    h3_by_id = {row["candidate_participant_id"]: row for row in h3_rows}
    if len(rna_by_id) != len(rna_rows) or len(h3_by_id) != len(h3_rows):
        raise GSE267145JoinError("deposited participant IDs are duplicated")
    if len(h3_by_id) != 99 or not set(h3_by_id) <= set(rna_by_id):
        raise GSE267145JoinError("99-person H3K27ac axis is not a subset of RNA")
    staged: dict[str, list[str]] = defaultdict(list)
    characteristics: dict[str, dict[str, str]] = {}
    for participant, row in h3_by_id.items():
        expected_prefix = row["title"].split("_", 1)[0]
        if expected_prefix != participant or rna_by_id[participant]["title"] != participant:
            raise GSE267145JoinError(f"deposited titles do not identify {participant}")
        values = parse_characteristics(row["characteristics_json"])
        required = {
            "tissue",
            "steatosis",
            "balloning",
            "lobular inflammation",
            "lobular necrosis",
            "fibrosis nas_tidy",
            "Stage",
            "Sex",
            "chip antibody",
        }
        if set(values) != required or values["Stage"] not in STAGES:
            raise GSE267145JoinError(f"phenotype fields differ for {participant}")
        if values["tissue"] != "Liver" or values["chip antibody"] != "H3K27ac":
            raise GSE267145JoinError(f"assay identity differs for {participant}")
        characteristics[participant] = values
        staged[values["Stage"]].append(participant)
    ranked: dict[str, int] = {}
    for stage in STAGES:
        ordered = sorted(
            staged[stage],
            key=lambda participant: (
                sha256(f"20260824\0{stage}\0{participant}".encode()).digest(),
                participant,
            ),
        )
        ranked.update({participant: rank for rank, participant in enumerate(ordered)})
    result: list[dict[str, str]] = []
    for participant in sorted(h3_by_id):
        h3 = h3_by_id[participant]
        rna = rna_by_id[participant]
        values = characteristics[participant]
        result.append(
            {
                "cohort_family_id": "gse267145_znf469_human_liver",
                "participant_id": participant,
                "rna_gsm": rna["accession"],
                "h3k27ac_gsm": h3["accession"],
                "stage": values["Stage"],
                "steatosis": values["steatosis"],
                "ballooning": values["balloning"],
                "lobular_inflammation": values["lobular inflammation"],
                "lobular_necrosis": values["lobular necrosis"],
                "fibrosis": values["fibrosis nas_tidy"],
                "sex": values["Sex"],
                "pairing": "same_sample_different_aliquot",
                "outer_fold": str(
                    participant_fold(
                        participant, values["Stage"], ranked[participant]
                    )
                ),
                "rna_measurement_state": "observed_unit_semantics_unresolved",
                "h3k27ac_measurement_state": "observed_integer_counts",
            }
        )
    return result


def _download(url: str, path: Path) -> dict[str, object]:
    request = urllib.request.Request(url, headers={"User-Agent": "masld-bench/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response, path.open("xb") as handle:
        digest = sha256()
        size = 0
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            handle.write(block)
            digest.update(block)
            size += len(block)
        headers = {
            key.lower(): value
            for key, value in response.headers.items()
            if key.lower() in {"content-length", "last-modified", "etag"}
        }
    return {"url": url, "bytes": size, "sha256": digest.hexdigest(), "headers": headers}


def verify_super_series(path: Path) -> None:
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        text = handle.read()
    if (
        "!Series_geo_accession = GSE267145" not in text
        or "GSE267119" not in text
        or "GSE269412" not in text
        or "39998893" not in text
    ):
        raise GSE267145JoinError("SuperSeries/subseries/citation evidence differs")


def verify_article(path: Path) -> None:
    root = ET.parse(path).getroot()
    text = re.sub(r"\s+", " ", " ".join(root.itertext())).lower()
    if not re.search(
        r"matched transcriptomes and epigenomes of 108 human liver biopsies", text
    ):
        raise GSE267145JoinError("published matched-biopsy statement is missing")
    if "pmc11071482" not in text and "39998893" not in text:
        raise GSE267145JoinError("article identity differs")


def verify_matrix_axes(matrix_root: Path, participants: set[str]) -> dict[str, int]:
    h3 = matrix_root / "raw" / "GSE267119_H3K27ac.txt.gz"
    rna = matrix_root / "raw" / "GSE269412_RNA.txt.gz"
    with gzip.open(h3, "rt", encoding="utf-8") as handle:
        h3_titles = shlex.split(next(handle))
    with gzip.open(rna, "rt", encoding="utf-8") as handle:
        rna_header = next(csv.reader(handle, delimiter="\t"))
    h3_ids = {value.split("_", 1)[0] for value in h3_titles}
    rna_ids = set(rna_header[1:])
    if h3_ids != participants or not participants <= rna_ids:
        raise GSE267145JoinError("matrix and participant join axes differ")
    return {"h3k27ac_columns": len(h3_titles), "rna_columns": len(rna_header) - 1}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-audit", type=Path, required=True)
    parser.add_argument("--metadata-artifacts-sha256", required=True)
    parser.add_argument("--matrix-audit", type=Path, required=True)
    parser.add_argument("--matrix-artifacts-sha256", required=True)
    parser.add_argument("--super-series-url", required=True)
    parser.add_argument("--article-xml-url", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (
        sha256_file(args.metadata_audit / "ARTIFACTS.json")
        != args.metadata_artifacts_sha256
        or sha256_file(args.matrix_audit / "ARTIFACTS.json")
        != args.matrix_artifacts_sha256
    ):
        raise GSE267145JoinError("input audit ARTIFACTS SHA-256 differs")
    args.output.mkdir(parents=True, exist_ok=False)
    raw = args.output / "raw"
    raw.mkdir()
    receipts = {
        "super_series": _download(
            args.super_series_url, raw / "GSE267145_family.soft.gz"
        ),
        "article_xml": _download(args.article_xml_url, raw / "PMC11071482.xml"),
    }
    verify_super_series(raw / "GSE267145_family.soft.gz")
    verify_article(raw / "PMC11071482.xml")
    h3_rows = read_samples(
        args.metadata_audit / "samples" / "GSE267119.tsv", "GSE267119"
    )
    rna_rows = read_samples(
        args.metadata_audit / "samples" / "GSE269412.tsv", "GSE269412"
    )
    joins = build_join_rows(h3_rows, rna_rows)
    matrix_axes = verify_matrix_axes(
        args.matrix_audit, {row["participant_id"] for row in joins}
    )
    fields = tuple(joins[0])
    with (args.output / "participant_join.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(joins)
    stage_counts = dict(sorted(Counter(row["stage"] for row in joins).items()))
    fold_stage_counts = {
        fold: dict(
            sorted(
                Counter(
                    row["stage"] for row in joins if row["outer_fold"] == fold
                ).items()
            )
        )
        for fold in map(str, range(5))
    }
    evidence = {
        "schema_version": "masld-bench-gse267145-authoritative-join-v1",
        "status": "pass",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "participants_joined": len(joins),
        "biological_unit": "participant",
        "pairing": "same_sample_different_aliquot",
        "stage_counts": stage_counts,
        "fold_stage_counts": fold_stage_counts,
        "matrix_axes": matrix_axes,
        "evidence_lines": [
            "official_GEO_SuperSeries_contains_both_subseries_and_PubMed_39998893",
            "published_article_states_matched_transcriptomes_and_epigenomes_of_108_liver_biopsies",
            "99_source_deposited_H3K27ac_participant_IDs_match_exactly_one_RNA_title_and_both_matrix_axes",
            "H3K27ac_records_supply_source_native_histology_and_sex",
        ],
        "participant_join_authoritative": True,
        "title_match_alone_used_as_authority": False,
        "activates_topology_only": True,
        "activates_model_training": False,
        "remaining_blockers": [
            "PISCES_RNA_fractional_unit_and_error_model",
            "Ensembl98_to_GENCODE49_gene_crosswalk",
            "H3K27ac_hg38_to_GRCh38p14_coordinate_and_contig_crosswalk",
            "assay_native_QC_and_rights_release_contract",
        ],
        "source_receipts": receipts,
        "metadata_artifacts_sha256": args.metadata_artifacts_sha256,
        "matrix_artifacts_sha256": args.matrix_artifacts_sha256,
    }
    (args.output / "join_evidence.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(evidence, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
