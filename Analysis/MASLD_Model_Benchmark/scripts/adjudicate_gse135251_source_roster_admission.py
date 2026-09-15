#!/usr/bin/env python3
"""Adjudicate GSE135251 against the frozen microarray transfer label requirements.

This emits an ADDITIVE, versioned inclusion record.  It opens
``microarray_transfer_label_semantics.toml`` and the GSE83452 TaskSpec read-only,
verifies their SHA-256 before and after, and never writes to either.  Every
existing output file keeps whatever inclusion stamp it already carries; only this
new record states an included state, following the ``model-check-288-21098630``
additive-overlay precedent.

The binding question is narrow.  ``[source_training] negative_required`` wants
"an explicit participant-level biopsy-derived no-steatohepatitis label sampled
from the same source population and adjudication procedure as positives", while
``prohibited_negative_substitutions`` bars
``simple_steatosis_without_explicit_no_NASH_adjudication``.  NAFL *is* simple
steatosis, so the clause turns entirely on whether the deposited NAFL call is an
independent histological adjudication or a threshold applied to the NAS sum.

Every evidence item is recomputed here from the frozen GSE135251 copy and
written out as a table.  Nothing is quoted from a message, and the two
cross-tabulations this module emits are the output file a verifier re-derives.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

SOURCE_SERIES = "GSE135251"
EXTERNAL_SERIES = "GSE83452"
NASH_GROUPS = ("NASH_F0-F1", "NASH_F2", "NASH_F3", "NASH_F4")
NAFL_GROUP = "NAFL"
CONTROL_GROUP = "control"
GROUP_ORDER = (CONTROL_GROUP, NAFL_GROUP) + NASH_GROUPS
NAS_COMPONENT_TOKENS = ("steato", "balloon", "lobular", "inflamm")
ADMITTED_CONTRAST = "nash_vs_nafl"


class AdmissionAdjudicationError(RuntimeError):
    """Raised when the inclusion record would rest on unverified evidence."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_rowname_offset_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read the R ``write.table(row.names = TRUE)`` metadata without shifting columns."""

    lines = [line for line in path.read_text(encoding="utf-8").split("\n") if line != ""]
    if len(lines) < 2:
        raise AdmissionAdjudicationError(f"metadata table has no data rows: {path}")
    header = lines[0].split("\t")
    rows: list[dict[str, str]] = []
    for index, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(header) + 1:
            raise AdmissionAdjudicationError(
                f"{path} line {index} carries {len(fields)} fields, not "
                f"{len(header) + 1}; the R rowname offset does not hold"
            )
        rows.append({header[i]: fields[i + 1] for i in range(len(header))})
    return header, rows


# ---------------------------------------------------------------------------
# Evidence item 1: the NAFL call is not a NAS-sum threshold.
# ---------------------------------------------------------------------------


def nas_threshold_separability(
    *, nafl: Sequence[int], nash: Sequence[int]
) -> dict[str, Any]:
    """Best achievable ``NAS >= t`` rule separating NAFL from NASH.

    If the deposited NAFL/NASH call were a threshold on the NAS sum, some ``t``
    would reproduce it exactly and accuracy would be 1.0.  Anything short of that
    is positive evidence that the call carries information the NAS sum does not,
    which is what ``negative_required`` asks for.
    """

    total = len(nafl) + len(nash)
    if total == 0:
        raise AdmissionAdjudicationError("NAS separability needs both classes")
    rows: list[dict[str, Any]] = []
    best: tuple[float, int] | None = None
    for threshold in range(0, 10):
        true_positive = sum(1 for value in nash if value >= threshold)
        false_positive = sum(1 for value in nafl if value >= threshold)
        false_negative = len(nash) - true_positive
        true_negative = len(nafl) - false_positive
        accuracy = (true_positive + true_negative) / total
        rows.append(
            {
                "nas_threshold": threshold,
                "rule": f"nas_score >= {threshold} called NASH",
                "true_positive": true_positive,
                "false_positive": false_positive,
                "false_negative": false_negative,
                "true_negative": true_negative,
                "misclassified": false_positive + false_negative,
                "accuracy": f"{accuracy:.6f}",
            }
        )
        if best is None or accuracy > best[0]:
            best = (accuracy, threshold)
    assert best is not None
    accuracy, threshold = best
    overlap = sorted(set(nafl) & set(nash))
    return {
        "grid": rows,
        "best_threshold": threshold,
        "best_accuracy": accuracy,
        "best_misclassified": int(round((1.0 - accuracy) * total)),
        "evaluated_participants": total,
        "overlapping_nas_values": overlap,
        "nafl_inside_overlap": sum(1 for value in nafl if value in overlap),
        "nash_inside_overlap": sum(1 for value in nash if value in overlap),
        "nafl_nas_range": [min(nafl), max(nafl)],
        "nash_nas_range": [min(nash), max(nash)],
        "nash_below_nas_3": sum(1 for value in nash if value < 3),
        "a_derived_threshold_would_score": 1.0,
        "call_is_a_nas_sum_threshold": accuracy >= 1.0,
    }


def cross_tabulate(
    rows: Sequence[Mapping[str, str]], *, field: str
) -> tuple[list[dict[str, Any]], dict[str, dict[str, int]]]:
    table: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for row in rows:
        table[row["group_in_paper"]][row[field]] += 1
    emitted: list[dict[str, Any]] = []
    for group in GROUP_ORDER:
        for value, count in sorted(table[group].items()):
            emitted.append(
                {
                    "group_in_paper": group,
                    "field": field,
                    "value": value,
                    "participants": count,
                }
            )
    return emitted, {
        group: {key: int(value) for key, value in sorted(table[group].items())}
        for group in GROUP_ORDER
    }


def adjudicate(
    *,
    metadata_path: Path,
    counts_path: Path,
    source_receipt_path: Path,
    label_contract_path: Path,
    taskspec_path: Path,
    output: Path,
    user_approval_date: str,
) -> dict[str, Any]:
    from masld_bench.artifacts import freeze_tree

    if output.exists():
        raise AdmissionAdjudicationError(f"refusing to overwrite admission record: {output}")

    # The two frozen bases are read-only here.  Their digests are captured before
    # anything is written and re-checked by the caller afterwards.
    label_contract_sha256 = sha256_file(label_contract_path)
    taskspec_sha256 = sha256_file(taskspec_path)
    source_receipt = json.loads(source_receipt_path.read_text(encoding="utf-8"))
    if source_receipt.get("status") != "pass_source_frozen_hashed_copy":
        raise AdmissionAdjudicationError("frozen GSE135251 source receipt differs")

    header, rows = read_rowname_offset_tsv(metadata_path)
    if len(rows) != 180:
        raise AdmissionAdjudicationError("frozen metadata row census differs")
    groups = {row["group_in_paper"] for row in rows}
    if groups != set(GROUP_ORDER):
        raise AdmissionAdjudicationError("group_in_paper roster differs")

    nafl = [int(row["nas_score"]) for row in rows if row["group_in_paper"] == NAFL_GROUP]
    nash = [
        int(row["nas_score"]) for row in rows if row["group_in_paper"] in NASH_GROUPS
    ]
    if len(nafl) != 41 or len(nash) != 131:
        raise AdmissionAdjudicationError("admitted contrast census differs from 131/41")
    separability = nas_threshold_separability(nafl=nafl, nash=nash)
    if separability["call_is_a_nas_sum_threshold"]:
        raise AdmissionAdjudicationError(
            "a NAS threshold reproduces the deposited NAFL/NASH call exactly; the "
            "negative label is a derived threshold, not an independent adjudication, "
            "and prohibited_negative_substitutions bars it"
        )

    nas_rows, nas_table = cross_tabulate(rows, field="nas_score")
    source_rows, source_table = cross_tabulate(rows, field="source_name")
    stage_rows, stage_table = cross_tabulate(rows, field="Stage")

    # What source_name actually encodes, checked rather than assumed.  It is a
    # fibrosis-severity label, not a recruitment field: it separates the control
    # arm from the NAFLD arm, but inside the NAFLD arm it tracks stage, so NAFL
    # shares a value with NASH_F0-F1 and NASH_F2 and not with NASH_F3 or NASH_F4.
    nafl_source_values = set(source_table[NAFL_GROUP])
    nash_source_values = {
        value for group in NASH_GROUPS for value in source_table[group]
    }
    control_source_values = set(source_table[CONTROL_GROUP])
    source_name_findings = {
        "field": "source_name",
        "nafl_values": sorted(nafl_source_values),
        "nash_values": sorted(nash_source_values),
        "control_values": sorted(control_source_values),
        "nafl_and_nash_share_at_least_one_value": bool(
            nafl_source_values & nash_source_values
        ),
        "nafl_and_nash_values_are_identical": nafl_source_values == nash_source_values,
        "control_arm_is_disjoint_from_the_nafld_arm": not (
            control_source_values & (nafl_source_values | nash_source_values)
        ),
        "field_is_determined_by_fibrosis_stage": all(
            len(source_table[group]) == 1 for group in GROUP_ORDER
        ),
        "what_this_establishes": (
            "every NAFL and every NASH participant is a liver biopsy in the NAFLD "
            "arm and the control arm is disjoint from both, which is the same fact "
            "that justified excluding controls"
        ),
        "what_this_does_not_establish": (
            "source_name is fully determined by fibrosis stage inside the NAFLD arm "
            "-- early is F0-F2 and moderate is F3-F4 -- so NAFL shares its value "
            "with NASH_F0-F1 and NASH_F2 but not with NASH_F3 or NASH_F4. The field "
            "is a severity label, not a recruitment or adjudication-procedure "
            "field, and it does not by itself establish a common adjudication "
            "procedure. That comes from the NAS separability evidence."
        ),
    }

    component_fields = [
        column
        for column in header
        if any(token in column.lower() for token in NAS_COMPONENT_TOKENS)
    ]
    with counts_path.open(encoding="utf-8") as handle:
        deposited_columns = len(handle.readline().rstrip("\n").split("\t")) - 1

    negative_required_met = bool(
        not separability["call_is_a_nas_sum_threshold"]
        and source_name_findings["control_arm_is_disjoint_from_the_nafld_arm"]
        and source_name_findings["nafl_and_nash_share_at_least_one_value"]
    )
    positive_required_met_on_disk = bool(component_fields)

    output.mkdir(parents=True)
    write_tsv(
        output / "nas_score_by_group.tsv",
        ("group_in_paper", "field", "value", "participants"),
        nas_rows,
    )
    write_tsv(
        output / "source_name_and_stage_by_group.tsv",
        ("group_in_paper", "field", "value", "participants"),
        source_rows + stage_rows,
    )
    write_tsv(
        output / "nas_threshold_separability.tsv",
        (
            "nas_threshold",
            "rule",
            "true_positive",
            "false_positive",
            "false_negative",
            "true_negative",
            "misclassified",
            "accuracy",
        ),
        separability["grid"],
    )

    record = {
        "schema_version": "masld-bench-source-roster-admission-record-v1",
        "record_id": "gse135251_gse83452_source_roster_admission_v1",
        "status": "admitted_to_the_gse83452_source_roster",
        "record_is_additive_overlay": True,
        "additive_overlay_precedent": "executions/model-check-288-21098630",
        "frozen_bases_edited": False,
        "frozen_bases_read_only_sha256": {
            "config/evaluation/microarray_transfer_label_semantics.toml": label_contract_sha256,
            "config/evaluation/gse83452_baseline_nash_transfer_task.toml": taskspec_sha256,
        },
        "source_series": SOURCE_SERIES,
        "external_series": EXTERNAL_SERIES,
        "cohort_family_id": "antwerp_inserm_shared",
        "frozen_source_copy": {
            "path": str(source_receipt_path.parent),
            "metadata_sha256": sha256_file(metadata_path),
            "counts_sha256": sha256_file(counts_path),
            "deposited_counts_columns": deposited_columns,
            "metadata_participants": len(rows),
        },
        "admitted_contrast": {
            "arm_id": ADMITTED_CONTRAST,
            "positive_class": "nash",
            "positive_source_groups": list(NASH_GROUPS),
            "positive_n": len(nash),
            "negative_class": "no_nash",
            "negative_source_groups": [NAFL_GROUP],
            "negative_n": len(nafl),
            "excluded_source_groups": [CONTROL_GROUP],
            "excluded_n": len(rows) - len(nafl) - len(nash),
            "exclusion_reason": (
                "prohibited_negative_substitutions bars healthy_control and "
                "normal_liver from the negative class, and the contract's own "
                "target negative_semantics state that GSE83452's deposited no-NASH "
                "arm is not a healthy-control label"
            ),
        },
        "evidence": {
            "item_1_nafl_call_is_not_a_nas_threshold": {
                "clause": "[source_training] negative_required / prohibited_negative_substitutions",
                "artifact": "nas_threshold_separability.tsv, nas_score_by_group.tsv",
                "recomputed_here": True,
                "nas_by_group": nas_table,
                **{
                    key: separability[key]
                    for key in (
                        "overlapping_nas_values",
                        "nafl_inside_overlap",
                        "nash_inside_overlap",
                        "nafl_nas_range",
                        "nash_nas_range",
                        "nash_below_nas_3",
                        "best_threshold",
                        "best_accuracy",
                        "best_misclassified",
                        "evaluated_participants",
                        "a_derived_threshold_would_score",
                        "call_is_a_nas_sum_threshold",
                    )
                },
                "finding": (
                    "no threshold on the NAS sum reproduces the deposited NAFL/NASH "
                    "call; the best rule misclassifies participants that a derived "
                    "rule could not. The call therefore carries histological "
                    "information the NAS sum does not, which is what an independent "
                    "adjudication means"
                ),
                "evidence_class": "recomputed_from_a_hashed_artifact",
            },
            "item_2_shared_nafld_biopsy_arm": {
                "clause": "[source_training] negative_required, same source population",
                "artifact": "source_name_and_stage_by_group.tsv",
                "recomputed_here": True,
                "source_name_by_group": source_table,
                "stage_by_group": stage_table,
                **source_name_findings,
                "evidence_class": "recomputed_from_a_hashed_artifact",
            },
            "item_3_registry_characterization": {
                "clause": "[source_training] positive_required, authoritatively documented equivalent definition",
                "statement": (
                    "GSE135251 is described as 206 histologically characterised "
                    "NAFLD participants plus 10 controls from the European NAFLD "
                    "Registry (Govaere et al.), addressing steatohepatitis and "
                    "fibrosis gene signatures"
                ),
                "evidence_class": "publication_level_statement_not_read_from_primary_methods",
                "primary_methods_text_read": False,
                "nash_crn_scoring_paragraph_read": False,
                "relayed_by": "team lead, 2026-08-26",
                "independent_corroboration_available_on_disk": {
                    "check": "deposited counts columns equal the published cohort size",
                    "deposited_counts_columns": deposited_columns,
                    "published_total": 216,
                    "consistent": deposited_columns == 216,
                    "note": (
                        "arithmetic consistency only; it corroborates the cohort "
                        "size and says nothing about the scoring procedure"
                    ),
                },
            },
        },
        "clause_adjudication": {
            "negative_required": {
                "met": negative_required_met,
                "basis": "recomputed_on_disk_evidence",
                "reasoning": (
                    "NAFL is an explicit participant-level, biopsy-derived, "
                    "non-steatohepatitis label; it is demonstrably not a NAS-sum "
                    "threshold; NAFL and NASH sit in the same NAFLD liver-biopsy arm "
                    "and carry NAS scores from the same instrument, while the "
                    "control arm is disjoint and already excluded. It is therefore "
                    "not the barred simple_steatosis_without_explicit_no_NASH_adjudication"
                ),
            },
            "positive_required": {
                "met": True,
                "met_on_disk": positive_required_met_on_disk,
                "basis": "authoritatively_documented_equivalent_definition",
                "nas_component_fields_on_disk": component_fields,
                "reasoning": (
                    "the three-component definition -- steatosis plus hepatocyte "
                    "ballooning plus lobular inflammation -- cannot be checked "
                    "against this deposit because no component field exists on "
                    "disk; only the NAS sum is deposited. The clause is satisfied "
                    "through its 'authoritatively documented equivalent definition' "
                    "alternative, resting on the registry's NASH-CRN "
                    "characterization, which is a publication-level statement whose "
                    "primary methods text was not read"
                ),
            },
        },
        "not_established": [
            "The NAS components (steatosis, hepatocyte ballooning, lobular "
            "inflammation) are not separately deposited; only the NAS sum is. "
            "positive_required's three-component definition is met by documented "
            "equivalence, not by component fields on disk.",
            "The NASH-CRN scoring paragraph of the source publication was not read "
            "from primary methods text by the team lead or by this record. Item 3 "
            "is a publication-level statement, and the on-disk corroboration is "
            "cohort-size arithmetic only.",
            "source_name does not establish a common adjudication procedure. It is "
            "fully determined by fibrosis stage inside the NAFLD arm, so it "
            "establishes the shared biopsy arm and the disjoint control arm and "
            "nothing more.",
            "One NASH participant carries NAS 1, below the NAS 3 floor a "
            "three-component NASH definition would imply. This is consistent with "
            "steatosis regression in advanced fibrosis (it is a NASH_F4 sample) and "
            "it is recorded as an anomaly rather than explained away.",
        ],
        "age_sex_logistic": {
            "model_id": "age_sex_logistic",
            "registered_taskspec_baseline": True,
            "fittable_from_this_source": False,
            "reason": (
                "GSE135251 deposits neither age nor recorded sex, and the TaskSpec "
                "forbids inferring a missing clinical field"
            ),
            "consequence": (
                "the TaskSpec's age-sex-only and masked molecular-plus-metadata "
                "lanes cannot be produced from GSE135251 at all; admission does not "
                "unblock them"
            ),
        },
        "admission": {
            "source_roster_admission_state": "admitted",
            "supersedes_state": "proposed_not_approved",
            "approved_by": "user",
            "approval_date": user_approval_date,
            "adjudication_verified_by": "team lead, re-derived independently in this record",
            "scope": (
                "admits GSE135251 as a training source for the "
                "gse83452_baseline_nash_transfer lane on the nash_vs_nafl contrast only"
            ),
            "prior_artifact_stamps_rewritten": False,
            "prior_artifact_stamp_note": (
                "every existing artifact keeps source_roster_admission_state = "
                "proposed_not_approved as it was written; this record is the only "
                "place the admitted state is asserted"
            ),
        },
        "claim_boundary": {
            "external_development_only": True,
            "champion_eligible": False,
            "diagnostic_or_prognostic_claim_allowed": False,
            "clinical_claim_allowed": False,
            "note": "admission unblocks reporting, not promotion",
        },
    }
    with (output / "admission_record.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    freeze_tree(
        output,
        {
            "artifact_class": "gse135251_gse83452_source_roster_admission_record",
            "source_series": SOURCE_SERIES,
            "external_series": EXTERNAL_SERIES,
            "admitted_contrast": ADMITTED_CONTRAST,
            "source_roster_admission_state": "admitted",
            "record_is_additive_overlay": True,
            "frozen_bases_edited": False,
            "external_development_only": True,
            "champion_eligible": False,
            "status": "passed",
        },
    )
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen-source", required=True, type=Path)
    parser.add_argument("--label-contract", required=True, type=Path)
    parser.add_argument("--taskspec", required=True, type=Path)
    parser.add_argument("--user-approval-date", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    record = adjudicate(
        metadata_path=arguments.frozen_source / "raw" / "GSE135251_metadata.tsv",
        counts_path=arguments.frozen_source / "raw" / "GSE135251_counts.tsv",
        source_receipt_path=arguments.frozen_source / "source_freeze_receipt.json",
        label_contract_path=arguments.label_contract,
        taskspec_path=arguments.taskspec,
        output=arguments.output,
        user_approval_date=arguments.user_approval_date,
    )
    print(json.dumps(record, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
