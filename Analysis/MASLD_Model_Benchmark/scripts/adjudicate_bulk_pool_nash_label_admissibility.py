#!/usr/bin/env python3
"""Adjudicate GSE162694, GSE240729 and GSE213621 against the frozen label requirements.

Same operation as the GSE135251 inclusion record and the same additive-overlay
discipline: ``microarray_transfer_label_semantics.toml`` and the GSE83452
TaskSpec are opened read-only, their SHA-256 captured before and after, and
neither is ever written.

SCOPE, which matters more here than it did for GSE135251.  The only
label-semantics authority in this repository is
``contract_id = "gse83452_source_native_nash_semantics_v1"``.  It adjudicates a
NASH / no-NASH source label.  It is NOT a fibrosis-label requirement, and no
fibrosis-label requirements exists.  So this record answers exactly one question per
cohort - may this cohort serve as a NASH training source for the
gse83452_baseline_nash_transfer lane - and it deliberately does not answer
whether a fibrosis stage is an eligible label for some future endpoint,
because there is no frozen text to adjudicate that against.

A cohort that fails is recorded with the clause it fails and the literal
prohibited substitution it matches.  The module raises rather than emit an
included state for a cohort that did not pass, so a FAIL cannot be upgraded by
editing a field.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

# Read verbatim from the frozen requirements' [source_training] block.  Kept here as
# literals so a drift in the requirement shows up as a mismatch, not as silence.
PROHIBITED_POSITIVE = (
    "generic_MASLD",
    "generic_NAFLD",
    "fibrosis_stage_only",
    "NAS_threshold_without_source_NASH_adjudication",
    "clinical_code_only",
    "transcriptome_inferred_state",
)
PROHIBITED_NEGATIVE = (
    "healthy_control",
    "normal_liver",
    "simple_steatosis_without_explicit_no_NASH_adjudication",
    "generic_non_MASLD",
    "undefined",
    "missing",
)

COHORTS = ("GSE162694", "GSE240729", "GSE213621")
LABEL_FIELD = {
    "GSE162694": "condition",
    "GSE240729": "fibrosisscore",
    "GSE213621": "fibrotic_stage",
}
EXPECTED_CENSUS = {
    "GSE162694": {
        "Control": 31,
        "NASH_F0": 35,
        "NASH_F1": 30,
        "NASH_F2": 27,
        "NASH_F3": 8,
        "NASH_F4": 12,
    },
    "GSE240729": {"F0": 9, "F1": 16, "F2": 25, "F3": 10, "F4": 6},
    "GSE213621": {"Control": 68, "F0F1": 97, "F2": 107, "F3F4": 95},
}
EXPECTED_ROWS = {"GSE162694": 143, "GSE240729": 66, "GSE213621": 367}


class LabelAdmissibilityError(RuntimeError):
    """Raised when the record would misstate a verdict or leak an included state."""


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
    """Read R ``write.table(row.names = TRUE)`` output, failing closed on the offset."""

    text = path.read_text(encoding="utf-8", errors="replace")
    lines = [line for line in text.split("\n") if line != ""]
    if len(lines) < 2:
        raise LabelAdmissibilityError(f"metadata table has no data rows: {path}")
    header = lines[0].split("\t")
    rows: list[dict[str, str]] = []
    for index, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(header) + 1:
            raise LabelAdmissibilityError(
                f"{path} line {index} carries {len(fields)} fields, not "
                f"{len(header) + 1}; the R rowname offset does not hold"
            )
        rows.append({header[i]: fields[i + 1] for i in range(len(header))})
    return header, rows


def classify_label_values(
    *, cohort: str, values: Mapping[str, int]
) -> dict[str, Any]:
    """Sort a cohort's deposited label values into NASH, no-NASH and control arms.

    The test is deliberately literal.  A value is an eligible positive only if
    it names steatohepatitis.  A value is an eligible negative only if it names
    a biopsy-derived non-steatohepatitis state that is not a healthy control.  A
    fibrosis stage names neither, which is exactly why the requirements list
    ``fibrosis_stage_only`` as a prohibited positive substitution.
    """

    nash: dict[str, int] = {}
    control: dict[str, int] = {}
    fibrosis_only: dict[str, int] = {}
    steatosis_no_nash: dict[str, int] = {}
    for value, count in sorted(values.items()):
        token = value.strip().lower()
        if token.startswith("nash") or token.startswith("mash"):
            nash[value] = count
        elif token in {"control", "normal", "healthy", "healthy_control"}:
            control[value] = count
        elif token in {"nafl", "steatosis", "simple_steatosis", "no_nash", "no nash"}:
            steatosis_no_nash[value] = count
        else:
            fibrosis_only[value] = count
    return {
        "cohort": cohort,
        "nash_values": nash,
        "healthy_control_values": control,
        "explicit_no_nash_values": steatosis_no_nash,
        "fibrosis_stage_or_other_values": fibrosis_only,
        "has_explicit_nash_positive": bool(nash),
        "has_admissible_negative": bool(steatosis_no_nash),
        "only_negative_is_a_healthy_control_arm": bool(control) and not steatosis_no_nash,
        "label_is_fibrosis_stage_only": bool(fibrosis_only) and not nash,
    }


def adjudicate_cohort(*, cohort: str, arms: Mapping[str, Any]) -> dict[str, Any]:
    """Return a per-cohort verdict citing the clause and the literal substitution."""

    failures: list[dict[str, str]] = []

    if arms["label_is_fibrosis_stage_only"]:
        failures.append(
            {
                "clause": "[source_training] positive_required",
                "prohibited_substitution": "fibrosis_stage_only",
                "detail": (
                    "the deposit carries a fibrosis stage and no steatohepatitis "
                    "label, so there is no NASH positive to adjudicate"
                ),
            }
        )
    elif not arms["has_explicit_nash_positive"]:
        failures.append(
            {
                "clause": "[source_training] positive_required",
                "prohibited_substitution": "generic_MASLD",
                "detail": "no explicit participant-level steatohepatitis label is deposited",
            }
        )

    if not arms["has_admissible_negative"]:
        if arms["only_negative_is_a_healthy_control_arm"]:
            failures.append(
                {
                    "clause": "[source_training] negative_required",
                    "prohibited_substitution": "healthy_control",
                    "detail": (
                        "the only non-NASH arm deposited is a healthy control arm; "
                        "no biopsy-derived no-steatohepatitis arm from the same "
                        "source population exists"
                    ),
                }
            )
        else:
            failures.append(
                {
                    "clause": "[source_training] negative_required",
                    "prohibited_substitution": "missing",
                    "detail": "no admissible non-steatohepatitis negative arm is deposited",
                }
            )

    for failure in failures:
        substitution = failure["prohibited_substitution"]
        if substitution not in set(PROHIBITED_POSITIVE) | set(PROHIBITED_NEGATIVE):
            raise LabelAdmissibilityError(
                f"{cohort} cites {substitution}, which the contract does not name"
            )
    verdict = "does_not_pass" if failures else "passes"
    return {
        "cohort": cohort,
        "verdict": verdict,
        "admissible_as_a_nash_training_source": not failures,
        "failed_clauses": failures,
        "clauses_failed": len(failures),
        **{key: arms[key] for key in sorted(arms) if key != "cohort"},
    }


def adjudicate(
    *,
    metadata_root: Path,
    label_contract_path: Path,
    taskspec_path: Path,
    output: Path,
) -> dict[str, Any]:
    from masld_bench.artifacts import freeze_tree

    if output.exists():
        raise LabelAdmissibilityError(f"refusing to overwrite record: {output}")
    label_contract_sha256 = sha256_file(label_contract_path)
    taskspec_sha256 = sha256_file(taskspec_path)
    contract_text = label_contract_path.read_text(encoding="utf-8")
    for substitution in PROHIBITED_POSITIVE + PROHIBITED_NEGATIVE:
        if f'"{substitution}"' not in contract_text:
            raise LabelAdmissibilityError(
                f"the frozen contract no longer names {substitution}; this "
                f"record's literals have drifted from the authority"
            )
    if 'contract_id = "gse83452_source_native_nash_semantics_v1"' not in contract_text:
        raise LabelAdmissibilityError("the frozen label contract identity differs")

    census_rows: list[dict[str, Any]] = []
    substrate_rows: list[dict[str, Any]] = []
    verdicts: list[dict[str, Any]] = []
    source_sha256: dict[str, str] = {}

    for cohort in COHORTS:
        path = metadata_root / cohort / f"{cohort}_metadata.tsv"
        source_sha256[path.name] = sha256_file(path)
        header, rows = read_rowname_offset_tsv(path)
        if len(rows) != EXPECTED_ROWS[cohort]:
            raise LabelAdmissibilityError(f"{cohort} row census differs")
        field = LABEL_FIELD[cohort]
        if field not in header:
            raise LabelAdmissibilityError(f"{cohort} lacks its label field {field}")
        census = collections.Counter(row[field] for row in rows)
        if dict(census) != EXPECTED_CENSUS[cohort]:
            raise LabelAdmissibilityError(
                f"{cohort} label census differs: {dict(sorted(census.items()))}"
            )
        for value, count in sorted(census.items()):
            census_rows.append(
                {
                    "cohort": cohort,
                    "label_field": field,
                    "label_value": value,
                    "samples": count,
                }
            )
        arms = classify_label_values(cohort=cohort, values=dict(census))
        verdicts.append(adjudicate_cohort(cohort=cohort, arms=arms))

        # Substrate facts recorded alongside the verdict: the tissue field, and
        # whether any sample identifier repeats.  Neither changes the verdict; both
        # are things a later endpoint would otherwise have to rediscover.
        tissue_fields = [
            column
            for column in header
            if column.lower() in {"cell_type", "tissue", "source_name"}
        ]
        tissue_values: dict[str, list[str]] = {}
        for column in tissue_fields:
            tissue_values[column] = sorted({row[column] for row in rows})
        identifier = next(
            (
                column
                for column in ("BioSample", "Sample Name", "sample_id", "Run", "SRX")
                if column in header
            ),
            None,
        )
        repeats = 0
        distinct = None
        if identifier is not None:
            counts = collections.Counter(row[identifier] for row in rows)
            distinct = len(counts)
            repeats = sum(1 for value in counts.values() if value > 1)
        substrate_rows.append(
            {
                "cohort": cohort,
                "metadata_rows": len(rows),
                "metadata_header_fields": len(header),
                "rowname_offset_invariant_holds": True,
                "label_field": field,
                "tissue_like_fields": json.dumps(tissue_values, sort_keys=True),
                "identifier_field": identifier or "",
                "distinct_identifiers": "" if distinct is None else distinct,
                "identifiers_appearing_more_than_once": repeats,
            }
        )

    admitted = [row["cohort"] for row in verdicts if row["admissible_as_a_nash_training_source"]]
    rejected = [row["cohort"] for row in verdicts if not row["admissible_as_a_nash_training_source"]]
    if admitted:
        raise LabelAdmissibilityError(
            "a cohort passed unexpectedly; an admission record must be written "
            "with its own evidence rather than emitted from this rejection sweep: "
            f"{admitted}"
        )

    output.mkdir(parents=True)
    write_tsv(
        output / "label_value_census.tsv",
        ("cohort", "label_field", "label_value", "samples"),
        census_rows,
    )
    write_tsv(
        output / "substrate_facts.tsv",
        (
            "cohort",
            "metadata_rows",
            "metadata_header_fields",
            "rowname_offset_invariant_holds",
            "label_field",
            "tissue_like_fields",
            "identifier_field",
            "distinct_identifiers",
            "identifiers_appearing_more_than_once",
        ),
        substrate_rows,
    )
    write_tsv(
        output / "clause_failures.tsv",
        ("cohort", "clause", "prohibited_substitution", "detail"),
        [
            {"cohort": row["cohort"], **failure}
            for row in verdicts
            for failure in row["failed_clauses"]
        ],
    )

    record = {
        "schema_version": "masld-bench-source-roster-admission-record-v1",
        "record_id": "bulk_existing_source_pool_nash_label_admissibility_v1",
        "status": "no_cohort_admitted",
        "record_is_additive_overlay": True,
        "additive_overlay_precedent": "executions/model-check-288-21098630",
        "frozen_bases_edited": False,
        "frozen_bases_read_only_sha256": {
            "config/evaluation/microarray_transfer_label_semantics.toml": label_contract_sha256,
            "config/evaluation/gse83452_baseline_nash_transfer_task.toml": taskspec_sha256,
        },
        "cohort_family_id": "bulk_existing_source_pool",
        "cohorts_adjudicated": list(COHORTS),
        "admitted": admitted,
        "does_not_pass": rejected,
        "source_metadata_sha256": source_sha256,
        "adjudicated_question": (
            "may this cohort serve as a NASH training source for the "
            "gse83452_baseline_nash_transfer lane, under "
            "contract_id gse83452_source_native_nash_semantics_v1"
        ),
        "scope_limits": {
            "this_is_a_nash_label_contract_not_a_fibrosis_label_contract": True,
            "a_frozen_fibrosis_label_semantics_contract_exists": False,
            "fibrosis_label_admissibility_adjudicated_here": False,
            "note": (
                "microarray_transfer_label_semantics.toml is the only "
                "label-semantics authority in this repository and it adjudicates a "
                "NASH / no-NASH source label. fibrosis_stage_only appears in "
                "prohibited_positive_substitutions as a rule about NASH positives. "
                "It does NOT establish that a fibrosis stage is an inadmissible "
                "label for a fibrosis endpoint. These cohorts failing here means "
                "they cannot be NASH training sources; it decides nothing about a "
                "composition-versus-fibrosis endpoint, which has no frozen label "
                "contract to be adjudicated against."
            ),
        },
        "verdicts": verdicts,
        "admission": {
            "source_roster_admission_state": "proposed_not_approved",
            "state_unchanged_by_this_record": True,
            "approved_by": None,
            "approval_date": None,
            "reason_no_approval_recorded": (
                "no cohort passed, so there is nothing to approve; this record "
                "asserts no admitted state for any cohort"
            ),
            "prior_artifact_stamps_rewritten": False,
        },
        "consequence": (
            "GSE135251 remains the only bulk_existing_source_pool cohort that has "
            "passed the frozen label contract. The contract's "
            "eligible_source_roster = [] is therefore not an oversight; on the "
            "evidence deposited for these three cohorts it is correct"
        ),
        "claim_boundary": {
            "external_development_only": True,
            "champion_eligible": False,
            "diagnostic_or_prognostic_claim_allowed": False,
            "clinical_claim_allowed": False,
        },
    }
    with (output / "admissibility_record.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    freeze_tree(
        output,
        {
            "artifact_class": "bulk_existing_source_pool_nash_label_admissibility_record",
            "cohort_family_id": "bulk_existing_source_pool",
            "cohorts_adjudicated": len(COHORTS),
            "admitted": len(admitted),
            "does_not_pass": len(rejected),
            "record_is_additive_overlay": True,
            "frozen_bases_edited": False,
            "status": "passed",
        },
    )
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-root", required=True, type=Path)
    parser.add_argument("--label-contract", required=True, type=Path)
    parser.add_argument("--taskspec", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    record = adjudicate(
        metadata_root=arguments.metadata_root,
        label_contract_path=arguments.label_contract,
        taskspec_path=arguments.taskspec,
        output=arguments.output,
    )
    print(json.dumps(record, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
