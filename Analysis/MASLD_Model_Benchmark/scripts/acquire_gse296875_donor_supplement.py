#!/usr/bin/env python3
"""Verify and inventory the authoritative GSE296875 donor-phenotype supplement.

The frozen TaskSpec pins one publisher workbook by SHA-256 and byte count.  A
workbook that differs in either is refused outright: a phenotype fixture may
never be built from a substitute file.  Nothing here recodes, imputes, or
infers a stage; the sheet is dumped losslessly so that the endpoint freeze runs
against the observed schema rather than a guessed one.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import tomllib
from typing import Any, Iterator
import xml.etree.ElementTree as ET
import zipfile

from masld_bench.artifacts import freeze_tree, write_json_exclusive, write_text_exclusive
from masld_bench.hashing import sha256_file


EXPECTED_SHEET = "S1 donorInfo"
LOW_CARDINALITY_LIMIT = 60
# The deposited sheet writes an unobserved donor field as the literal string
# "NA", not as an empty cell.  Treating a sentinel as an observed value would
# silently turn missing pathology into a score, so both are counted apart.
MISSING_SENTINELS = ("NA",)

_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_DOCUMENT_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_PACKAGE_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"


class SupplementAcquisitionError(ValueError):
    """Raised when the publisher supplement or its worksheet differs."""


def _column_index(cell_reference: str) -> int:
    match = re.fullmatch(r"([A-Z]+)[1-9][0-9]*", cell_reference)
    if match is None:
        raise SupplementAcquisitionError(f"XLSX cell reference differs: {cell_reference!r}")
    index = 0
    for letter in match.group(1):
        index = index * 26 + ord(letter) - ord("A") + 1
    return index - 1


def sheet_names(workbook_path: Path) -> list[str]:
    """Return worksheet names in deposited order."""

    with zipfile.ZipFile(workbook_path) as archive:
        try:
            workbook = ET.parse(archive.open("xl/workbook.xml")).getroot()
        except (KeyError, ET.ParseError) as error:
            raise SupplementAcquisitionError("supplement workbook structure differs") from error
        return [
            str(sheet.get("name"))
            for sheet in workbook.findall(f".//{{{_MAIN_NS}}}sheet")
        ]


def _shared_strings(archive: zipfile.ZipFile) -> list[str]:
    try:
        handle = archive.open("xl/sharedStrings.xml")
    except KeyError:
        return []
    with handle:
        try:
            root = ET.parse(handle).getroot()
        except ET.ParseError as error:
            raise SupplementAcquisitionError("supplement shared strings differ") from error
    return [
        "".join(value.text or "" for value in item.iter(f"{{{_MAIN_NS}}}t"))
        for item in root.findall(f"{{{_MAIN_NS}}}si")
    ]


def iter_sheet_cells(workbook_path: Path, sheet_name: str) -> Iterator[list[dict[str, Any]]]:
    """Yield one lossless cell record list per worksheet row.

    Each record carries the deposited reference, the declared cell type, the
    raw ``<v>`` text, and the rendered value.  Numeric text is preserved so a
    float repr can never silently alter a deposited pathology score.
    """

    with zipfile.ZipFile(workbook_path) as archive:
        try:
            workbook = ET.parse(archive.open("xl/workbook.xml")).getroot()
            relationships = ET.parse(archive.open("xl/_rels/workbook.xml.rels")).getroot()
        except (KeyError, ET.ParseError) as error:
            raise SupplementAcquisitionError("supplement workbook structure differs") from error
        matches = [
            sheet
            for sheet in workbook.findall(f".//{{{_MAIN_NS}}}sheet")
            if sheet.get("name") == sheet_name
        ]
        if len(matches) != 1:
            raise SupplementAcquisitionError(
                f"supplement worksheet identity differs: {sheet_name!r}"
            )
        targets = {
            item.get("Id"): str(item.get("Target"))
            for item in relationships.findall(f".//{{{_PACKAGE_REL_NS}}}Relationship")
        }
        target = targets.get(matches[0].get(f"{{{_DOCUMENT_REL_NS}}}id"))
        if target is None:
            raise SupplementAcquisitionError("supplement worksheet relationship is absent")
        member = f"xl/{target.lstrip('/')}" if not target.startswith("xl/") else target
        shared = _shared_strings(archive)
        try:
            sheet_handle = archive.open(member)
        except KeyError as error:
            raise SupplementAcquisitionError(f"supplement worksheet is absent: {member}") from error
        with sheet_handle:
            try:
                for _, element in ET.iterparse(sheet_handle, events=("end",)):
                    if element.tag != f"{{{_MAIN_NS}}}row":
                        continue
                    row: list[dict[str, Any]] = []
                    for cell in element.findall(f"{{{_MAIN_NS}}}c"):
                        reference = cell.get("r", "")
                        column = _column_index(reference)
                        if column < len(row):
                            raise SupplementAcquisitionError("supplement cell order differs")
                        while len(row) < column:
                            row.append(
                                {"ref": None, "type": "blank", "raw": None, "value": None}
                            )
                        row.append(_render_cell(cell, reference, shared))
                    yield row
                    element.clear()
            except ET.ParseError as error:
                raise SupplementAcquisitionError("supplement worksheet XML differs") from error


def _render_cell(cell: ET.Element, reference: str, shared: list[str]) -> dict[str, Any]:
    cell_type = cell.get("t")
    if cell_type == "inlineStr":
        inline = cell.find(f"{{{_MAIN_NS}}}is")
        text = (
            "".join(value.text or "" for value in inline.iter(f"{{{_MAIN_NS}}}t"))
            if inline is not None
            else None
        )
        return {"ref": reference, "type": "inlineStr", "raw": text, "value": text}
    value = cell.find(f"{{{_MAIN_NS}}}v")
    raw = None if value is None else value.text
    if raw is None:
        return {"ref": reference, "type": "blank", "raw": None, "value": None}
    if cell_type == "s":
        try:
            text = shared[int(raw)]
        except (ValueError, IndexError) as error:
            raise SupplementAcquisitionError("supplement shared-string index differs") from error
        return {"ref": reference, "type": "s", "raw": raw, "value": text}
    if cell_type in {"str", "e"}:
        return {"ref": reference, "type": cell_type, "raw": raw, "value": raw}
    if cell_type == "b":
        return {"ref": reference, "type": "b", "raw": raw, "value": raw == "1"}
    if cell_type is None or cell_type == "n":
        try:
            float(raw)
        except ValueError as error:
            raise SupplementAcquisitionError("supplement numeric cell differs") from error
        return {"ref": reference, "type": "n", "raw": raw, "value": raw}
    raise SupplementAcquisitionError(f"supplement cell type differs: {cell_type!r}")


def _is_blank(record: dict[str, Any]) -> bool:
    value = record["value"]
    return value is None or (isinstance(value, str) and not value.strip())


def summarize_sheet(rows: list[list[dict[str, Any]]]) -> dict[str, Any]:
    """Describe the observed worksheet schema without recoding any value."""

    if not rows:
        raise SupplementAcquisitionError("supplement worksheet is empty")
    # The first populated row of this sheet is a prose caption occupying one
    # cell.  A header must name at least two columns, which skips the caption
    # without hard-coding its text.
    header_index = next(
        (
            index
            for index, row in enumerate(rows)
            if sum(1 for cell in row if not _is_blank(cell)) >= 2
        ),
        None,
    )
    if header_index is None:
        raise SupplementAcquisitionError("supplement worksheet has no header row")
    header = [
        str(cell["value"]).strip() if not _is_blank(cell) else ""
        for cell in rows[header_index]
    ]
    while header and header[-1] == "":
        header.pop()
    if not header:
        raise SupplementAcquisitionError("supplement header row is blank")
    if len(set(header)) != len(header):
        raise SupplementAcquisitionError(f"supplement header is not unique: {header}")
    body = rows[header_index + 1 :]
    body = [row for row in body if any(not _is_blank(cell) for cell in row)]
    columns = []
    for index, name in enumerate(header):
        non_blank = []
        for row in body:
            record = row[index] if index < len(row) else {"value": None}
            if not _is_blank(record):
                non_blank.append(str(record["value"]).strip())
        observed = [value for value in non_blank if value not in MISSING_SENTINELS]
        distinct = sorted(set(non_blank))
        columns.append(
            {
                "column_index": index,
                "name": name,
                "non_blank_rows": len(non_blank),
                "blank_rows": len(body) - len(non_blank),
                "sentinel_missing_rows": len(non_blank) - len(observed),
                "observed_rows": len(observed),
                "distinct_count": len(distinct),
                "cell_types": sorted(
                    {
                        (row[index]["type"] if index < len(row) else "blank")
                        for row in body
                    }
                ),
                "distinct_values": distinct if len(distinct) <= LOW_CARDINALITY_LIMIT else None,
            }
        )
    return {
        "header_row_index": header_index,
        "header": header,
        "caption_rows": [
            [cell["value"] for cell in row if not _is_blank(cell)]
            for row in rows[:header_index]
        ],
        "missing_sentinels": list(MISSING_SENTINELS),
        "data_rows": len(body),
        "columns": columns,
    }


def _tsv(header: list[str], rows: list[list[dict[str, Any]]], header_index: int) -> str:
    lines = ["\t".join(header)]
    for row in rows[header_index + 1 :]:
        if all(_is_blank(cell) for cell in row):
            continue
        fields = []
        for index in range(len(header)):
            record = row[index] if index < len(row) else {"value": None}
            if _is_blank(record):
                fields.append("")
                continue
            text = str(record["value"]).strip()
            if any(character in text for character in "\t\r\n"):
                raise SupplementAcquisitionError(f"supplement cell holds a delimiter: {text!r}")
            fields.append(text)
        lines.append("\t".join(fields))
    return "\n".join(lines) + "\n"


def run(*, workbook: Path, taskspec: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise SupplementAcquisitionError("refusing to overwrite supplement artifact")
    with taskspec.open("rb") as handle:
        spec = tomllib.load(handle)
    digest = sha256_file(workbook)
    size = workbook.stat().st_size
    if digest != spec["source_supplement_sha256"]:
        raise SupplementAcquisitionError(
            f"supplement SHA-256 differs: observed {digest}, "
            f"frozen {spec['source_supplement_sha256']}"
        )
    if size != spec["source_supplement_size_bytes"]:
        raise SupplementAcquisitionError(
            f"supplement byte count differs: observed {size}, "
            f"frozen {spec['source_supplement_size_bytes']}"
        )
    if spec["source_sheet"] != EXPECTED_SHEET:
        raise SupplementAcquisitionError("frozen source sheet name differs")
    names = sheet_names(workbook)
    if EXPECTED_SHEET not in names:
        raise SupplementAcquisitionError(
            f"authoritative sheet is absent; observed sheets {names}"
        )
    rows = list(iter_sheet_cells(workbook, EXPECTED_SHEET))
    summary = summarize_sheet(rows)
    declared = list(spec["source_fields"])
    missing_fields = [field for field in declared if field not in summary["header"]]

    output.mkdir(mode=0o750, parents=True)
    (output / "supplement").mkdir(mode=0o750)
    (output / "inventory").mkdir(mode=0o750)
    shutil.copyfile(workbook, output / "supplement" / workbook.name)
    (output / "supplement" / workbook.name).chmod(0o440)
    write_json_exclusive(
        output / "inventory" / "workbook_sheets.json",
        {"sheet_names": names, "authoritative_sheet": EXPECTED_SHEET},
    )
    write_json_exclusive(output / "inventory" / "s1_donorinfo_schema.json", summary)
    write_json_exclusive(
        output / "inventory" / "s1_donorinfo_cells.json",
        [
            [
                {key: cell[key] for key in ("ref", "type", "raw", "value")}
                for cell in row
            ]
            for row in rows
        ],
    )
    write_text_exclusive(
        output / "inventory" / "s1_donorinfo.tsv",
        _tsv(summary["header"], rows, summary["header_row_index"]),
    )
    receipt = {
        "source_url": spec["source_supplement"],
        "source_doi": spec["source_doi"],
        "source_sheet": EXPECTED_SHEET,
        "supplement_sha256": digest,
        "supplement_size_bytes": size,
        "sheet_names": names,
        "declared_source_fields": declared,
        "declared_source_fields_absent_from_header": missing_fields,
        "data_rows": summary["data_rows"],
        "header_row_index": summary["header_row_index"],
        "header": summary["header"],
        "missing_sentinels": list(MISSING_SENTINELS),
        "analyzed_donors_declared": spec["analyzed_donors"],
        "recoding": "none",
        "missingness": "explicit_endpoint_mask_never_absence",
    }
    write_json_exclusive(output / "inventory" / "acquisition_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse296875_donor_supplement_acquisition",
            "dataset_id": "gse296875",
            "role": "secondary_development_only",
            "external_or_sealed": False,
            "champion_eligible": False,
            "supplement_sha256": digest,
            "supplement_size_bytes": size,
            "source_sheet": EXPECTED_SHEET,
            "automatic_download": False,
            "download_authority": "explicit_phenotype_lane_instruction_task_0",
            "redistribution_class": "source_reference_only",
        },
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", required=True, type=Path)
    parser.add_argument("--taskspec", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = run(
        workbook=arguments.workbook,
        taskspec=arguments.taskspec,
        output=arguments.output,
    )
    print(json.dumps({"output": arguments.output.as_posix(), "receipt": receipt}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
