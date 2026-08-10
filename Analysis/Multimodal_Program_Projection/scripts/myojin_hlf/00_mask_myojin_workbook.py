#!/usr/bin/env python3
"""Phase-A custodian: authenticate S1 and export only outcome-masked fields.

This process is allowed to parse the raw workbook.  It never writes or prints a
MAGeCK effect, score, p value, FDR, rank, goodsgrna value, or derived hit call.
The independent freezer consumes only the allow-listed masked table.
"""

from __future__ import annotations

import argparse
import collections
import datetime as dt
import re
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from myojin_firewall_common import (
    CANDIDATE_ROOT,
    RELEASE_ID,
    atomic_write_text,
    require_within,
    sha256_file,
    stable_bundle_sha256,
    write_tsv,
)


NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS = {"m": NS_MAIN, "r": NS_REL}
EXPECTED_SHEET = "Supplementary Table1"
EXPECTED_DIMENSION = "A1:Q18346"
EXPECTED_RECORDS = 18343
EXPECTED_BLOCKS = {
    "B2": "Palmitic acid (D21 vs D0)",
    "J2": "Palmitic acid  vs Vehicle",
}
EXPECTED_HEADERS = {
    "A3": "id",
    "B3": "num",
    "C3": "pos|score",
    "D3": "pos|p-value",
    "E3": "pos|fdr",
    "F3": "pos|rank",
    "G3": "pos|goodsgrna",
    "H3": "pos|lfc",
    "I3": "-Log10pvalue_Pos",
    "J3": "num",
    "K3": "pos|score",
    "L3": "pos|p-value",
    "M3": "pos|fdr",
    "N3": "pos|rank",
    "O3": "pos|goodsgrna",
    "P3": "pos|lfc",
    "Q3": "-Log10pvalue_Pos",
}
MASKED_FIELDS = [
    "source_row",
    "source_gene_id",
    "guide_count_d21_d0",
    "guide_count_pa_vehicle",
    "d21_d0_score_present",
    "d21_d0_pvalue_present",
    "d21_d0_fdr_present",
    "d21_d0_lfc_present",
    "pa_vehicle_score_present",
    "pa_vehicle_pvalue_present",
    "pa_vehicle_fdr_present",
    "pa_vehicle_lfc_present",
    "pa_vehicle_primary_fields_complete",
    "source_gene_id_duplicate",
    "outcome_values_exported",
]
FORBIDDEN_EXPORT_NAMES = {
    "score",
    "pvalue",
    "fdr",
    "rank",
    "goodsgrna",
    "lfc",
    "protective_hit",
}


def column_letters(cell_reference: str) -> str:
    match = re.match(r"^([A-Z]+)", cell_reference)
    if not match:
        raise ValueError(f"Invalid cell reference: {cell_reference}")
    return match.group(1)


def shared_strings(archive: ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return [
        "".join(node.text or "" for node in item.iterfind(".//m:t", NS))
        for item in root.findall("m:si", NS)
    ]


def cell_value(cell: ET.Element, strings: list[str]) -> str:
    cell_type = cell.attrib.get("t", "")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.iterfind(".//m:t", NS))
    value = cell.find("m:v", NS)
    raw = "" if value is None else (value.text or "")
    if cell_type == "s" and raw:
        return strings[int(raw)]
    if cell_type == "b":
        return "TRUE" if raw == "1" else "FALSE"
    return raw


def parse_integer(value: str, cell: str) -> int | str:
    if value == "":
        return ""
    number = float(value)
    if not number.is_integer() or number < 0:
        raise ValueError(f"Expected nonnegative integer in {cell}; observed {value!r}")
    return int(number)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=CANDIDATE_ROOT)
    parser.add_argument("--source-url", default="https://links.lww.com/HC9/C218")
    parser.add_argument("--retrieved-at-utc", default="")
    args = parser.parse_args()

    outdir = require_within(args.outdir, CANDIDATE_ROOT)
    workbook = require_within(args.workbook, outdir)
    source_dir = outdir / "source"
    source_dir.mkdir(parents=True, exist_ok=True)
    retrieved = args.retrieved_at_utc or dt.datetime.now(dt.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )

    with ZipFile(workbook) as archive:
        names = set(archive.namelist())
        required_parts = {
            "xl/workbook.xml",
            "xl/worksheets/sheet1.xml",
            "xl/_rels/workbook.xml.rels",
        }
        missing_parts = sorted(required_parts - names)
        if missing_parts:
            raise ValueError(f"Workbook missing required parts: {missing_parts}")
        strings = shared_strings(archive)
        workbook_xml = ET.fromstring(archive.read("xl/workbook.xml"))
        sheets = workbook_xml.find("m:sheets", NS)
        sheet_names = [node.attrib["name"] for node in sheets] if sheets is not None else []
        if sheet_names != [EXPECTED_SHEET]:
            raise ValueError(f"Expected one sheet {EXPECTED_SHEET!r}; observed {sheet_names}")
        sheet_xml = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))

    dimension_node = sheet_xml.find("m:dimension", NS)
    dimension = "" if dimension_node is None else dimension_node.attrib.get("ref", "")
    if dimension != EXPECTED_DIMENSION:
        raise ValueError(f"Expected used range {EXPECTED_DIMENSION}; observed {dimension}")

    rows_by_number: dict[int, dict[str, str]] = {}
    formula_counts = collections.Counter()
    nonmissing_counts = collections.Counter()
    for row in sheet_xml.findall(".//m:row", NS):
        row_number = int(row.attrib["r"])
        parsed: dict[str, str] = {}
        for cell in row.findall("m:c", NS):
            reference = cell.attrib["r"]
            column = column_letters(reference)
            parsed[column] = cell_value(cell, strings)
            if cell.find("m:f", NS) is not None:
                formula_counts[column] += 1
            if parsed[column] != "":
                nonmissing_counts[column] += 1
        rows_by_number[row_number] = parsed

    for reference, expected in {**EXPECTED_BLOCKS, **EXPECTED_HEADERS}.items():
        row_number = int(re.sub(r"^[A-Z]+", "", reference))
        observed = rows_by_number.get(row_number, {}).get(column_letters(reference), "")
        if observed != expected:
            raise ValueError(
                f"Header mismatch at {reference}: expected {expected!r}; observed {observed!r}"
            )

    data_rows = [rows_by_number.get(index, {}) for index in range(4, 18347)]
    if len(data_rows) != EXPECTED_RECORDS:
        raise AssertionError("Internal record-range error")
    gene_counts = collections.Counter(row.get("A", "").strip() for row in data_rows)
    if "" in gene_counts:
        raise ValueError(f"Blank source gene identifiers: {gene_counts['']}")

    masked_rows = []
    presence_columns = {
        "d21_d0_score_present": "C",
        "d21_d0_pvalue_present": "D",
        "d21_d0_fdr_present": "E",
        "d21_d0_lfc_present": "H",
        "pa_vehicle_score_present": "K",
        "pa_vehicle_pvalue_present": "L",
        "pa_vehicle_fdr_present": "M",
        "pa_vehicle_lfc_present": "P",
    }
    for source_row, row in enumerate(data_rows, start=4):
        gene_id = row.get("A", "").strip()
        record: dict[str, object] = {
            "source_row": source_row,
            "source_gene_id": gene_id,
            "guide_count_d21_d0": parse_integer(row.get("B", ""), f"B{source_row}"),
            "guide_count_pa_vehicle": parse_integer(row.get("J", ""), f"J{source_row}"),
            "source_gene_id_duplicate": str(gene_counts[gene_id] > 1).upper(),
            "outcome_values_exported": "FALSE",
        }
        for field, column in presence_columns.items():
            record[field] = str(row.get(column, "") != "").upper()
        record["pa_vehicle_primary_fields_complete"] = str(
            all(bool(row.get(column, "")) for column in ("K", "L", "M", "P"))
        ).upper()
        masked_rows.append(record)

    # Defense in depth: the masked table has a fixed allow-list, and no bare
    # outcome value field can be added without breaking this assertion.
    for field in MASKED_FIELDS:
        normalized = re.sub(r"[^a-z0-9]", "", field.lower())
        if normalized in FORBIDDEN_EXPORT_NAMES:
            raise AssertionError(f"Forbidden result value field in masked output: {field}")
    if any(row["outcome_values_exported"] != "FALSE" for row in masked_rows):
        raise AssertionError("Outcome export marker is not uniformly FALSE")

    masked_path = outdir / "masked_screen_schema.tsv"
    write_tsv(masked_path, masked_rows, MASKED_FIELDS)

    header_map = {column_letters(ref): value for ref, value in EXPECTED_HEADERS.items()}
    schema_rows = []
    for column in list("ABCDEFGHIJKLMNOPQ"):
        block = (
            "identifier"
            if column == "A"
            else "pa_d21_vs_d0"
            if column in "BCDEFGHI"
            else "pa_vs_vehicle"
        )
        source_field = header_map[column]
        export_policy = (
            "source_gene_identifier"
            if column == "A"
            else "neutral_guide_count"
            if column in {"B", "J"}
            else "missingness_only_no_values"
        )
        schema_rows.append(
            {
                "column": column,
                "block": block,
                "source_field": source_field,
                "n_nonmissing_including_headers": nonmissing_counts[column],
                "n_formulas": formula_counts[column],
                "blind_phase_export_policy": export_policy,
            }
        )
    schema_path = outdir / "source_schema_audit.tsv"
    write_tsv(
        schema_path,
        schema_rows,
        [
            "column",
            "block",
            "source_field",
            "n_nonmissing_including_headers",
            "n_formulas",
            "blind_phase_export_policy",
        ],
    )

    source_manifest_path = outdir / "source_manifest.tsv"
    write_tsv(
        source_manifest_path,
        [
            {
                "source_id": "myojin_2026_supplemental_table_s1",
                "source_url": args.source_url,
                "source_paper_doi": "10.1097/HC9.0000000000000884",
                "source_paper_pmid": "41564380",
                "source_paper_pmcid": "PMC12826143",
                "source_paper_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12826143/",
                "supplement_relationship": "Supplemental_Table_S1",
                "retrieved_at_utc": retrieved,
                "relative_path": str(workbook.relative_to(outdir)),
                "bytes": workbook.stat().st_size,
                "sha256": sha256_file(workbook),
                "sheet_count": 1,
                "sheet_names": EXPECTED_SHEET,
                "used_range": dimension,
                "gene_records": EXPECTED_RECORDS,
                "license_status": "article_CC_BY_4.0_supplement_publisher_hosted",
                "license_url": "https://creativecommons.org/licenses/by/4.0/",
                "raw_t0_counts_status": "not_recovered_public_no_permission_search_2026-08-07",
                "outcome_values_exported_to_blind_phase": "FALSE",
                "release_id": RELEASE_ID,
            }
        ],
        [
            "source_id",
            "source_url",
            "source_paper_doi",
            "source_paper_pmid",
            "source_paper_pmcid",
            "source_paper_url",
            "supplement_relationship",
            "retrieved_at_utc",
            "relative_path",
            "bytes",
            "sha256",
            "sheet_count",
            "sheet_names",
            "used_range",
            "gene_records",
            "license_status",
            "license_url",
            "raw_t0_counts_status",
            "outcome_values_exported_to_blind_phase",
            "release_id",
        ],
    )

    bundle_hash, bundle_rows = stable_bundle_sha256(
        [masked_path, schema_path, source_manifest_path], outdir
    )
    phase_a_manifest = outdir / "phase_a_manifest.tsv"
    write_tsv(
        phase_a_manifest,
        [dict(row, release_id=RELEASE_ID) for row in bundle_rows],
        ["relative_path", "bytes", "sha256", "release_id"],
    )
    atomic_write_text(
        outdir / "PHASE_A_READY",
        "release_id\tstatus\tbundle_sha256\toutcome_values_exported\n"
        f"{RELEASE_ID}\tmasked_source_ready\t{bundle_hash}\tFALSE\n",
    )
    print(
        "PHASE_A_READY "
        f"records={EXPECTED_RECORDS} duplicate_ids={sum(v > 1 for v in gene_counts.values())} "
        f"bundle_sha256={bundle_hash} outcome_values_exported=FALSE"
    )


if __name__ == "__main__":
    main()
