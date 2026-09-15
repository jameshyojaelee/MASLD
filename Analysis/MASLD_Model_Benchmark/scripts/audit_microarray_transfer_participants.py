#!/usr/bin/env python3
"""Freeze participant and repeated-biopsy joins for GSE49541 and GSE83452."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path
import re


SOURCE_FIELDS = (
    "series",
    "accession",
    "title",
    "source_name",
    "organism",
    "platform_id",
    "library_strategy",
    "description_json",
    "characteristics_json",
    "data_processing_json",
    "relations_json",
    "supplementary_files_json",
)
GSE49541_TITLE = re.compile(r"^NAFLD liver biopsy tissue (\d+)$")
GSE83452_TITLE = re.compile(r"^liver biopsy (\d+)(?: \((\d+)\))?$")
GSE83452_SAMPLE = re.compile(r"^patient (\d+)(?: \((\d+)\))?$")


class MicroarrayParticipantAuditError(RuntimeError):
    """Raised when deposited participant or repeat evidence is inconsistent."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_samples(path: Path, expected_series: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != SOURCE_FIELDS:
            raise MicroarrayParticipantAuditError("source sample fields differ")
        rows = [dict(row) for row in reader]
    if not rows or any(row["series"] != expected_series for row in rows):
        raise MicroarrayParticipantAuditError("source sample series differs")
    return rows


def parse_characteristics(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in json.loads(value):
        if ":" not in item:
            raise MicroarrayParticipantAuditError(f"unkeyed characteristic: {item}")
        key, observed = item.split(":", 1)
        key, observed = key.strip(), observed.strip()
        if key in result or not observed:
            raise MicroarrayParticipantAuditError(f"duplicate or empty characteristic: {key}")
        result[key] = observed
    return result


def build_gse49541(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for row in rows:
        match = GSE49541_TITLE.fullmatch(row["title"])
        values = parse_characteristics(row["characteristics_json"])
        if (
            match is None
            or set(values) != {"Stage", "tissue"}
            or values["tissue"] != "liver"
            or row["platform_id"] != "GPL570"
        ):
            raise MicroarrayParticipantAuditError("GSE49541 identity or assay differs")
        stage = values["Stage"]
        mapping = {
            "mild (fibrosis stage 0-1)": "mild_f0_f1",
            "advanced (fibrosis stage 3-4)": "advanced_f3_f4",
        }
        if stage not in mapping:
            raise MicroarrayParticipantAuditError("GSE49541 fibrosis group differs")
        source_key = match.group(1)
        result.append(
            {
                "cohort_family_id": "gse31803_gse49541_fibrosis_array",
                "series": "GSE49541",
                "participant_id": f"gse49541::{source_key}",
                "sample_accession": row["accession"],
                "source_sample_key": source_key,
                "timepoint": "baseline_cross_sectional",
                "repeat_topology": "one_record_per_participant",
                "fibrosis_stage_group": mapping[stage],
                "exact_fibrosis_stage": "structurally_missing",
                "nash_status": "structurally_missing",
                "age": "structurally_missing",
                "sex": "structurally_missing",
                "intervention": "not_applicable",
                "primary_transfer_evaluable": True,
                "paired_expression_stress_evaluable": False,
                "paired_nash_transition_evaluable": False,
            }
        )
    counts = Counter(str(row["fibrosis_stage_group"]) for row in result)
    if (
        len(result) != 72
        or len({str(row["participant_id"]) for row in result}) != 72
        or counts != Counter({"mild_f0_f1": 40, "advanced_f3_f4": 32})
    ):
        raise MicroarrayParticipantAuditError("GSE49541 participant census differs")
    return result


def build_gse83452(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    parsed: list[tuple[dict[str, str], dict[str, str], re.Match[str]]] = []
    baseline_by_number: dict[str, tuple[dict[str, str], dict[str, str]]] = {}
    required = {
        "sample name",
        "liver status",
        "type of intervention",
        "time",
        "age",
        "gender",
        "scan date",
        "tissue",
    }
    for row in rows:
        values = parse_characteristics(row["characteristics_json"])
        title = GSE83452_TITLE.fullmatch(row["title"])
        sample = GSE83452_SAMPLE.fullmatch(values.get("sample name", ""))
        if (
            title is None
            or sample is None
            or title.groups() != sample.groups()
            or set(values) != required
            or values["tissue"] != "liver biopsy"
            or values["time"] not in {"baseline", "follow-up"}
            or values["liver status"] not in {"NASH", "no NASH", "undefined"}
            or values["type of intervention"] not in {"BS", "Diet"}
            or values["gender"] not in {"female", "male"}
            or row["platform_id"] != "GPL16686"
        ):
            raise MicroarrayParticipantAuditError("GSE83452 identity or phenotype differs")
        try:
            age = int(values["age"])
        except ValueError as error:
            raise MicroarrayParticipantAuditError("GSE83452 age is not integral") from error
        if not 18 <= age <= 90:
            raise MicroarrayParticipantAuditError("GSE83452 age is outside admission bounds")
        if values["time"] == "baseline":
            if sample.group(2) is not None or sample.group(1) in baseline_by_number:
                raise MicroarrayParticipantAuditError("GSE83452 baseline key differs")
            baseline_by_number[sample.group(1)] = (row, values)
        parsed.append((row, values, sample))
    result: list[dict[str, object]] = []
    for row, values, sample in parsed:
        record_number, baseline_number = sample.groups()
        is_followup = values["time"] == "follow-up"
        if not is_followup and baseline_number is not None:
            raise MicroarrayParticipantAuditError("baseline record has a repeat pointer")
        if baseline_number is not None:
            if baseline_number not in baseline_by_number:
                raise MicroarrayParticipantAuditError("follow-up baseline pointer is absent")
            baseline_row, baseline_values = baseline_by_number[baseline_number]
            if (
                baseline_values["gender"] != values["gender"]
                or baseline_values["type of intervention"] != values["type of intervention"]
                or int(values["age"]) - int(baseline_values["age"]) not in {1, 2}
            ):
                raise MicroarrayParticipantAuditError("paired participant metadata differs")
            participant_id = f"gse83452::{baseline_row['accession']}"
            baseline_gsm = baseline_row["accession"]
            repeat_topology = "paired_baseline_one_year_followup"
        elif is_followup:
            participant_id = f"gse83452::{row['accession']}"
            baseline_gsm = "structurally_missing"
            repeat_topology = "followup_only_unpaired"
        else:
            participant_id = f"gse83452::{row['accession']}"
            baseline_gsm = row["accession"]
            repeat_topology = "baseline_only_or_pair_anchor"
        status = values["liver status"]
        result.append(
            {
                "cohort_family_id": "antwerp_inserm_shared",
                "series": "GSE83452",
                "participant_id": participant_id,
                "sample_accession": row["accession"],
                "source_sample_number": record_number,
                "paired_baseline_source_number": baseline_number or "not_applicable",
                "paired_baseline_accession": baseline_gsm,
                "timepoint": values["time"],
                "repeat_topology": repeat_topology,
                "nash_status": status.lower().replace(" ", "_"),
                "age": int(values["age"]),
                "sex": values["gender"],
                "intervention": values["type of intervention"],
                "scan_date": values["scan date"],
                "primary_transfer_evaluable": values["time"] == "baseline"
                and status != "undefined",
                "paired_expression_stress_evaluable": baseline_number is not None,
                "paired_nash_transition_evaluable": baseline_number is not None
                and status != "undefined"
                and baseline_by_number[baseline_number][1]["liver status"] != "undefined",
            }
        )
    time_counts = Counter(str(row["timepoint"]) for row in result)
    source_counts = Counter(
        (str(row["timepoint"]), str(row["nash_status"])) for row in result
    )
    if (
        len(result) != 231
        or len({str(row["participant_id"]) for row in result}) != 171
        or time_counts != Counter({"baseline": 152, "follow-up": 79})
        or sum(row["repeat_topology"] == "paired_baseline_one_year_followup" for row in result) != 60
        or sum(row["repeat_topology"] == "followup_only_unpaired" for row in result) != 19
        or sum(bool(row["primary_transfer_evaluable"]) for row in result) != 148
        or sum(bool(row["paired_expression_stress_evaluable"]) for row in result)
        != 60
        or sum(bool(row["paired_nash_transition_evaluable"]) for row in result)
        != 54
        or source_counts
        != Counter(
            {
                ("baseline", "nash"): 104,
                ("baseline", "no_nash"): 44,
                ("baseline", "undefined"): 4,
                ("follow-up", "nash"): 22,
                ("follow-up", "no_nash"): 54,
                ("follow-up", "undefined"): 3,
            }
        )
    ):
        raise MicroarrayParticipantAuditError("GSE83452 participant census differs")
    return result


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--source-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if sha256_file(arguments.source / "ARTIFACTS.json") != arguments.source_artifacts_sha256:
        raise MicroarrayParticipantAuditError("source ARTIFACTS SHA-256 differs")
    source_manifest = json.loads(
        (arguments.source / "ARTIFACTS.json").read_text(encoding="utf-8")
    )
    if (
        source_manifest["metadata"].get("official_geo_sample_records") != 303
        or source_manifest["metadata"].get("model_training_activated") is not False
    ):
        raise MicroarrayParticipantAuditError("source artifact role differs")
    arguments.output.mkdir(parents=True, exist_ok=False)
    gse49541 = build_gse49541(
        read_samples(arguments.source / "samples/GSE49541.tsv", "GSE49541")
    )
    gse83452 = build_gse83452(
        read_samples(arguments.source / "samples/GSE83452.tsv", "GSE83452")
    )
    write_rows(arguments.output / "gse49541_participants.tsv", gse49541)
    write_rows(arguments.output / "gse83452_records.tsv", gse83452)
    summary = {
        "schema_version": "masld-bench-microarray-transfer-participant-audit-v2",
        "status": "pass_participant_topology_and_source_phenotypes",
        "gse49541": {
            "cohort_family_id": "gse31803_gse49541_fibrosis_array",
            "participants": 72,
            "records": 72,
            "fibrosis_groups": {"mild_f0_f1": 40, "advanced_f3_f4": 32},
            "exact_fibrosis_stage_available": False,
            "age_available": False,
            "sex_available": False,
        },
        "gse83452": {
            "cohort_family_id": "antwerp_inserm_shared",
            "accession_aliases": ["GSE106737", "GSE83452"],
            "participants": 171,
            "records": 231,
            "baseline_records": 152,
            "followup_records": 79,
            "paired_participants": 60,
            "unpaired_followup_participants": 19,
            "baseline_nash_transfer_evaluable": 148,
            "paired_expression_stress_evaluable": 60,
            "paired_nash_transition_evaluable": 54,
            "age_available": True,
            "sex_available": True,
            "intervention_available": True,
        },
        "participant_join_authoritative_from_deposited_geo_fields": True,
        "raw_array_acquisition_complete": False,
        "platform_probe_mapping_complete": False,
        "rights_audit_complete": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
        "claims": [
            "GSE49541 supports external-development transfer for a source-deposited F0-F1 versus F3-F4 fibrosis group endpoint, not exact fibrosis-stage regression.",
            "GSE83452 supports baseline NASH-status transfer and a separately named paired one-year intervention stress test.",
            "Undefined NASH records remain explicit and are not mapped to a negative class.",
            "All 60 deposited baseline-follow-up pairs support outcome-blind expression-change stress evaluation; only the 54 pairs with defined NASH status at both visits support NASH-transition summaries.",
            "All repeated records from one participant must remain in one split.",
            "GSE106737 and GSE83452 are one project-exposed cohort family; cross-accession aliases are frozen by the separate label-blind fingerprint audit.",
        ],
    }
    (arguments.output / "participant_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
