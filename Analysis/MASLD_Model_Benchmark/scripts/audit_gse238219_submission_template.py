#!/usr/bin/env python3
"""Audit the author-supplied GEO sample/replicate submission template."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import io
import json
from pathlib import Path, PurePosixPath
import re
import tarfile


class SubmissionTemplateError(RuntimeError):
    """Raised when the author-submitted replicate template differs."""


def safe_name(value: str) -> None:
    pure = PurePosixPath(value)
    if pure.is_absolute() or ".." in pure.parts:
        raise SubmissionTemplateError(f"unsafe archive member: {value}")


def extract_template(source: Path) -> bytes:
    outer_name = "GSM7660624_P_Seq2_crispr_analysis.tar.gz"
    inner_suffix = "/PerturbSeq_seq_template_PSG.csv"
    with tarfile.open(source, mode="r:") as outer:
        matches = [member for member in outer if member.isfile() and member.name == outer_name]
        if len(matches) != 1:
            raise SubmissionTemplateError("submission-template outer archive differs")
        stream = outer.extractfile(matches[0])
        if stream is None:
            raise SubmissionTemplateError("cannot read submission-template archive")
        with tarfile.open(fileobj=stream, mode="r|gz") as nested:
            for member in nested:
                safe_name(member.name)
                if member.isfile() and member.name.endswith(inner_suffix):
                    handle = nested.extractfile(member)
                    if handle is None:
                        raise SubmissionTemplateError("cannot read submission template")
                    return handle.read()
    raise SubmissionTemplateError("submission template is absent")


def nonempty(row: list[str]) -> list[str]:
    return [value.strip() for value in row if value.strip()]


def parse(payload: bytes) -> dict[str, object]:
    rows = list(csv.reader(io.StringIO(payload.decode("utf-8-sig"), newline="")))
    biological_definition = next(
        (values[0] for row in rows if (values := nonempty(row)) and "Biological replicates of the same sample" in values[0]),
        None,
    )
    technical_definition = next(
        (values[0] for row in rows if (values := nonempty(row)) and "Technical replicates" in values[0]),
        None,
    )
    if biological_definition is None or technical_definition is None:
        raise SubmissionTemplateError("replicate definitions are absent")
    sample_header_index = next(
        index for index, row in enumerate(rows) if row and row[0].strip() == "*library name"
    )
    sample_header = [value.strip() for value in rows[sample_header_index]]
    sample_rows = []
    for row in rows[sample_header_index + 1 :]:
        if not row or not row[0].strip():
            break
        padded = row + [""] * (len(sample_header) - len(row))
        record = dict(zip(sample_header, (value.strip() for value in padded), strict=False))
        sample_rows.append(
            {
                "library_name": record["*library name"],
                "title": record["*title"],
                "organism": record["*organism"],
                "cell_line": record["**cell line"],
                "cell_type": record["**cell type"],
                "genotype": record["genotype"],
                "treatment": record["treatment"],
                "instrument_model_as_submitted": record["*instrument model"],
            }
        )
    expected_titles = [f"Perturb-Seq replicate {index}" for index in range(1, 6)]
    if [record["title"] for record in sample_rows] != expected_titles:
        raise SubmissionTemplateError("five replicate titles differ")
    lane_pattern = re.compile(r"_L(\d{3})_")
    lane_sets: dict[str, list[str]] = {}
    for library in (record["library_name"] for record in sample_rows):
        lanes = sorted(
            {
                match.group(1)
                for row in rows
                for value in row
                if library in value
                for match in [lane_pattern.search(value)]
                if match is not None
            }
        )
        lane_sets[library] = lanes
    return {
        "biological_replicate_definition": biological_definition,
        "technical_replicate_definition": technical_definition,
        "sample_rows": sample_rows,
        "sequencing_lane_ids_by_library": lane_sets,
        "author_submission_designation": "five separately titled replicate sample rows",
        "evidence_interpretation": "The author-supplied template distinguishes biological replicate rows from technical sequencing lanes and places five titled replicates in separate sample rows. It proves five submitted samples/captures and author replicate designation, but does not itself document independent culture or transduction ancestry.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.source.is_file() or args.output.exists():
        raise SubmissionTemplateError("submission-template request differs")
    payload = extract_template(args.source)
    record = {
        "schema_version": "masld-bench-gse238219-submission-template-audit-v1",
        "source_path": args.source.as_posix(),
        "template_sha256": sha256(payload).hexdigest(),
        "template_size_bytes": len(payload),
        **parse(payload),
        "source_perturbation_effect_tables_read": False,
        "gse313774_accessed": False,
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "submission_template_audit.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
