#!/usr/bin/env python3
"""PASS-01 fixtures plus PASS-06 MASLD Gene Catalog release validation.

Legacy ``passport_*`` identifiers remain part of the sealed v1 file contract.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import stat
from collections import Counter, defaultdict, deque
from pathlib import Path, PurePosixPath
from typing import Mapping

import pandas as pd
import pyarrow.parquet as pq

import build_evidence_passport_ui as passport_ui
import generate_evidence_passports as contract
import render_evidence_passport_ui as passport_ui_renderer


SCRIPT_PATH = Path(__file__).resolve()
LOGICAL_PRODUCER_ID = "scripts/portal/validate_evidence_passports.py"

class PassportValidationError(contract.PassportContractError):
    pass


REQUIRED_BUNDLE_FILES = {
    *contract.PARQUET_SCHEMAS,
    *contract.TSV_SCHEMAS,
    "passport_input_selection.tsv",
    "passport_input_selection.signature.json",
    "fixture_source_evidence.tsv",
}
PRODUCTION_PRETERMINAL_FILES = {
    *contract.PRODUCTION_PARQUET_SCHEMAS,
    *contract.PRODUCTION_TSV_SCHEMAS,
    *contract.PORTAL_PARQUET_SCHEMAS,
    "passport_input_selection.tsv",
    "passport_input_selection.signature.json",
    "portal_export/passport_dictionary.json",
    "portal_export/passport_source_graph.json",
    "portal_export/passport_domain_coverage.json",
    "portal_export/passport_bundle_status.json",
    "portal_candidate/index.html",
    "portal_candidate/ui_contract.json",
    "portal_candidate/ui_source_manifest.tsv",
    "portal_candidate/review/overview.png",
    "portal_candidate/review/THRB.png",
    "portal_candidate/review/HKDC1.png",
    "portal_candidate/review/GLP1R.png",
    "portal_candidate/review/MTARC1.png",
    "portal_candidate/review/review_manifest.tsv",
}
PRODUCTION_PRETERMINAL_FILES -= {
    "passport_validation_report.tsv",
    "passport_plan60_handoff.tsv",
}
FORBIDDEN_COLUMN_PATTERN = re.compile(
    r"(^|_)(score|rank|ranking|priority|tier|posterior|probability|"
    r"n_modalities|modality_count|combined_evidence|total_evidence)($|_)",
    re.IGNORECASE,
)
NUMERIC_EVIDENCE_COLUMNS = [
    "estimate",
    "standard_error",
    "ci_lower",
    "ci_upper",
    "p_value",
    "q_value",
]
ENSG_PATTERN = re.compile(r"^ENSG\d{11}$")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
TERMINAL_PROVENANCE_FILENAME = "passport_terminal_provenance.tsv"
TERMINAL_PROVENANCE_CHECK_ID = "PASS06_TERMINAL_PROVENANCE"
TERMINAL_CHAIN_PATHS = frozenset(
    {
        "PASS06_VALIDATED",
        "FIXTURE_PASS06_VALIDATED",
        "passport_input_selection.tsv",
        "passport_release_manifest.tsv",
        "passport_validation_report.tsv",
        "passport_plan60_handoff.tsv",
        TERMINAL_PROVENANCE_FILENAME,
    }
)
POST_SEAL_EXCLUDED_TREES = frozenset({"fixtures", "superseded_bundles"})


def release_bundle_uri(analysis_release_id: str) -> str:
    """Return the stable release-relative URI used inside a movable bundle."""

    return f"release://{analysis_release_id}/"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PassportValidationError("EMPTY_TABLE", str(path))
        return list(reader)


def table_columns(path: Path) -> list[str]:
    if path.suffix == ".parquet":
        return pq.read_schema(path).names
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            return next(reader)
        except StopIteration as exc:
            raise PassportValidationError("EMPTY_TABLE", str(path)) from exc


def check_inventory(bundle: Path) -> None:
    missing = sorted(
        name for name in REQUIRED_BUNDLE_FILES if not (bundle / name).is_file()
    )
    if missing:
        raise PassportValidationError("MISSING_BUNDLE_FILES", ", ".join(missing))


def check_production_inventory(bundle: Path) -> None:
    missing = sorted(
        name for name in PRODUCTION_PRETERMINAL_FILES if not (bundle / name).is_file()
    )
    if missing:
        raise PassportValidationError("MISSING_PRODUCTION_FILES", ", ".join(missing))
    if (bundle / "PASS06_VALIDATED").exists():
        raise PassportValidationError(
            "STALE_PASS06_SEAL",
            "terminal seal must not predate the current validation run",
        )


def check_production_schemas(bundle: Path) -> None:
    for filename, expected in {
        **contract.PRODUCTION_PARQUET_SCHEMAS,
        **contract.PORTAL_PARQUET_SCHEMAS,
    }.items():
        observed = pq.read_schema(bundle / filename)
        if observed.names != list(expected):
            raise PassportValidationError(
                "PARQUET_SCHEMA_COLUMNS",
                f"{filename}: expected {list(expected)}, observed {observed.names}",
            )
        for field in observed:
            if field.type != expected[field.name]:
                raise PassportValidationError(
                    "PARQUET_SCHEMA_TYPE",
                    f"{filename}:{field.name}: expected {expected[field.name]}, observed {field.type}",
                )
    for filename, expected in contract.PRODUCTION_TSV_SCHEMAS.items():
        if filename in {
            "passport_validation_report.tsv",
            "passport_plan60_handoff.tsv",
        }:
            continue
        observed = table_columns(bundle / filename)
        if observed != expected:
            raise PassportValidationError(
                "TSV_SCHEMA_COLUMNS",
                f"{filename}: expected {expected}, observed {observed}",
            )


def check_forbidden_columns(bundle: Path) -> None:
    offenders = []
    for name in sorted(contract.PARQUET_SCHEMAS):
        for column in table_columns(bundle / name):
            if FORBIDDEN_COLUMN_PATTERN.search(column):
                offenders.append(f"{name}:{column}")
    for name in sorted(contract.TSV_SCHEMAS):
        if name == "passport_data_dictionary.tsv":
            continue
        for column in table_columns(bundle / name):
            if FORBIDDEN_COLUMN_PATTERN.search(column):
                offenders.append(f"{name}:{column}")
    if offenders:
        raise PassportValidationError("FORBIDDEN_SCORE_OR_RANK", ", ".join(offenders))


def check_production_forbidden_columns(bundle: Path) -> None:
    offenders = []
    parquet_files = {
        **contract.PRODUCTION_PARQUET_SCHEMAS,
        **contract.PORTAL_PARQUET_SCHEMAS,
    }
    for name in sorted(parquet_files):
        for column in table_columns(bundle / name):
            if FORBIDDEN_COLUMN_PATTERN.search(column):
                offenders.append(f"{name}:{column}")
    for name in sorted(contract.PRODUCTION_TSV_SCHEMAS):
        if name in {
            "passport_data_dictionary.tsv",
            "passport_validation_report.tsv",
            "passport_plan60_handoff.tsv",
        }:
            continue
        for column in table_columns(bundle / name):
            if FORBIDDEN_COLUMN_PATTERN.search(column):
                offenders.append(f"{name}:{column}")
    if offenders:
        raise PassportValidationError("FORBIDDEN_SCORE_OR_RANK", ", ".join(offenders))


def check_schemas(bundle: Path) -> None:
    for filename, expected in contract.PARQUET_SCHEMAS.items():
        observed = pq.read_schema(bundle / filename)
        if observed.names != list(expected):
            raise PassportValidationError(
                "PARQUET_SCHEMA_COLUMNS",
                f"{filename}: expected {list(expected)}, observed {observed.names}",
            )
        for field in observed:
            if field.type != expected[field.name]:
                raise PassportValidationError(
                    "PARQUET_SCHEMA_TYPE",
                    f"{filename}:{field.name}: expected {expected[field.name]}, observed {field.type}",
                )
    for filename, expected in contract.TSV_SCHEMAS.items():
        observed = table_columns(bundle / filename)
        if observed != expected:
            raise PassportValidationError(
                "TSV_SCHEMA_COLUMNS",
                f"{filename}: expected {expected}, observed {observed}",
            )


def check_identity(bundle: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genes = pd.read_parquet(bundle / "passport_gene_index.parquet", engine="pyarrow")
    evidence = pd.read_parquet(
        bundle / "passport_evidence_long.parquet", engine="pyarrow"
    )
    coverage = pd.read_parquet(
        bundle / "passport_coverage_long.parquet", engine="pyarrow"
    )
    if genes["passport_id"].duplicated().any():
        values = genes.loc[
            genes["passport_id"].duplicated(False), "passport_id"
        ].tolist()
        raise PassportValidationError("DUPLICATE_PASSPORT_ID", str(values))
    if genes["ensembl_id"].duplicated().any():
        raise PassportValidationError(
            "DUPLICATE_ENSEMBL_ID", "gene index must be one row per stable ID"
        )
    if genes["passport_id"].astype(str).tolist() != sorted(
        genes["passport_id"].astype(str)
    ):
        raise PassportValidationError(
            "NONDETERMINISTIC_KEY_ORDER",
            "gene index must use stable Catalog-entry ID order, never an effect-derived order",
        )
    invalid = [
        value
        for value in genes["ensembl_id"].astype(str)
        if not ENSG_PATTERN.fullmatch(value)
    ]
    if invalid:
        raise PassportValidationError("INVALID_ENSEMBL_ID", str(invalid[:5]))
    expected_passports = (
        genes["analysis_release_id"].astype(str) + ":" + genes["ensembl_id"].astype(str)
    )
    if not expected_passports.equals(genes["passport_id"].astype(str)):
        raise PassportValidationError(
            "PASSPORT_ID_DERIVATION", "Catalog entry ID must be release:Ensembl"
        )
    actual_collision = genes["symbol"].duplicated(keep=False)
    if not actual_collision.equals(genes["symbol_collision"].astype(bool)):
        raise PassportValidationError(
            "SYMBOL_COLLISION_FLAG",
            "duplicate symbols must remain separate Ensembl rows and carry symbol_collision=true",
        )
    if evidence["evidence_result_id"].duplicated().any():
        raise PassportValidationError(
            "DUPLICATE_EVIDENCE_RESULT", "evidence_result_id is not unique"
        )
    evidence_order = list(
        zip(
            evidence["passport_id"].astype(str),
            evidence["evidence_result_id"].astype(str),
        )
    )
    if evidence_order != sorted(evidence_order):
        raise PassportValidationError(
            "NONDETERMINISTIC_KEY_ORDER", "evidence rows must use stable identity order"
        )
    composite = [
        "passport_id",
        "assay",
        "dataset_id",
        "phenotype",
        "context",
        "contrast_or_exposure",
        "gate_id",
    ]
    if evidence.duplicated(composite).any():
        raise PassportValidationError("DUPLICATE_INFERENTIAL_RESULT", str(composite))
    if coverage["coverage_result_id"].duplicated().any():
        raise PassportValidationError(
            "DUPLICATE_COVERAGE_RESULT", "coverage_result_id is not unique"
        )
    coverage_order = list(
        zip(
            coverage["passport_id"].astype(str),
            coverage["coverage_result_id"].astype(str),
        )
    )
    if coverage_order != sorted(coverage_order):
        raise PassportValidationError(
            "NONDETERMINISTIC_KEY_ORDER", "coverage rows must use stable identity order"
        )
    return genes, evidence, coverage


def _normalized_records(
    rows: list[dict[str, object]], columns: list[str]
) -> list[dict[str, str]]:
    return [
        {
            column: ""
            if row.get(column) is None or pd.isna(row.get(column))
            else str(row.get(column))
            for column in columns
        }
        for row in rows
    ]


def _compare_parquet_to_rows(
    bundle: Path, filename: str, rows: list[dict[str, object]], schema
) -> None:
    expected = contract.coerce_dataframe(rows, schema)
    sort_columns = [
        column
        for column in (
            "passport_id",
            "evidence_result_id",
            "coverage_result_id",
            "program_uid",
            "ensembl_id",
            "dataset_id",
            "program_context_id",
            "source_input_row_id",
        )
        if column in expected.columns
    ]
    if sort_columns:
        expected = expected.sort_values(sort_columns, kind="mergesort").reset_index(
            drop=True
        )
    observed = pd.read_parquet(bundle / filename, engine="pyarrow")
    if sort_columns:
        observed = observed.sort_values(sort_columns, kind="mergesort").reset_index(
            drop=True
        )
    try:
        pd.testing.assert_frame_equal(
            observed, expected, check_dtype=True, check_like=False
        )
    except AssertionError as exc:
        raise PassportValidationError(
            "ADAPTER_REDERIVATION_MISMATCH", f"{filename}: {str(exc)[:500]}"
        ) from exc


def rederive_production_bundle(bundle: Path, fixture_mode: bool):
    import build_evidence_passport_bundle as builder

    selection_path = bundle / "passport_input_selection.tsv"
    selected, signature = builder._resolve_selected_inputs(  # noqa: SLF001
        selection_path, allow_fixture=fixture_mode
    )
    assembly, _ = builder.adapt_selection(
        selected, signature, allow_fixture=fixture_mode
    )
    genes = sorted(assembly.genes.values(), key=lambda row: row["passport_id"])
    evidence = sorted(
        assembly.evidence,
        key=lambda row: (row["passport_id"], row["evidence_result_id"]),
    )
    coverage = sorted(
        assembly.coverage,
        key=lambda row: (row["passport_id"], row["coverage_result_id"]),
    )
    programs = sorted(assembly.programs.values(), key=lambda row: row["program_uid"])
    memberships = sorted(
        assembly.memberships,
        key=lambda row: (
            row["program_uid"],
            row["ensembl_id"],
            row["source_input_row_id"],
        ),
    )
    contexts = sorted(
        assembly.program_context,
        key=lambda row: (
            row["program_uid"],
            row["dataset_id"],
            row["program_context_id"],
        ),
    )
    for filename, rows in (
        ("passport_gene_index.parquet", genes),
        ("passport_evidence_long.parquet", evidence),
        ("passport_coverage_long.parquet", coverage),
        ("passport_program_index.parquet", programs),
        ("passport_program_membership.parquet", memberships),
        ("passport_program_context.parquet", contexts),
    ):
        _compare_parquet_to_rows(
            bundle, filename, rows, contract.PRODUCTION_PARQUET_SCHEMAS[filename]
        )
    expected_nodes, expected_edges = builder._source_graph(selected, assembly)  # noqa: SLF001
    for filename, expected_rows in (
        ("passport_source_nodes.tsv", expected_nodes),
        ("passport_source_edges.tsv", expected_edges),
        ("passport_gen_identity_status.tsv", assembly.gen_identity_status),
        ("passport_domain_coverage.tsv", builder._domain_coverage(selected, assembly)),  # noqa: SLF001
    ):
        columns = contract.PRODUCTION_TSV_SCHEMAS[filename]
        if read_tsv(bundle / filename) != _normalized_records(expected_rows, columns):
            raise PassportValidationError(
                "SOURCE_GRAPH_REDERIVATION_MISMATCH", filename
            )
    return assembly, selected, signature


def check_program_boundaries(
    bundle: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    programs = pd.read_parquet(
        bundle / "passport_program_index.parquet", engine="pyarrow"
    )
    membership = pd.read_parquet(
        bundle / "passport_program_membership.parquet", engine="pyarrow"
    )
    context = pd.read_parquet(
        bundle / "passport_program_context.parquet", engine="pyarrow"
    )
    if programs["program_uid"].duplicated().any():
        raise PassportValidationError("DUPLICATE_PROGRAM_UID", "program index")
    if programs["member_gene_call_expansion_authorized"].fillna(True).any():
        raise PassportValidationError(
            "PROGRAM_GENE_EXPANSION", "program index authorizes member-gene calls"
        )
    if membership["gene_call_expansion_authorized"].fillna(True).any():
        raise PassportValidationError(
            "PROGRAM_GENE_EXPANSION", "membership authorizes member-gene calls"
        )
    if (~membership["membership_only"].fillna(False)).any():
        raise PassportValidationError(
            "PROGRAM_MEMBERSHIP_SCOPE", "membership_only must be true"
        )
    if context["gene_call_expansion_authorized"].fillna(True).any():
        raise PassportValidationError(
            "PROGRAM_GENE_EXPANSION", "program context authorizes gene calls"
        )
    program_ids = set(programs["program_uid"].astype(str))
    for label, frame in (("membership", membership), ("context", context)):
        unresolved = set(frame["program_uid"].astype(str)) - program_ids
        if unresolved:
            raise PassportValidationError(
                "PROGRAM_REFERENCE", f"{label}:{sorted(unresolved)}"
            )
    if membership.duplicated(
        ["program_uid", "ensembl_id", "source_input_row_id"]
    ).any():
        raise PassportValidationError(
            "DUPLICATE_PROGRAM_MEMBERSHIP", "program/gene/source row"
        )
    if context["program_context_id"].duplicated().any():
        raise PassportValidationError("DUPLICATE_PROGRAM_CONTEXT", "program_context_id")
    check_semantics(context)
    skipped = {
        row["dataset_id"]
        for row in read_tsv(bundle / "passport_dataset_status.tsv")
        if row["status"] == "skipped_source_gate"
    }
    leaked = skipped & set(context["dataset_id"].astype(str))
    if leaked:
        raise PassportValidationError(
            "SKIPPED_DATASET_PROGRAM_LEAK",
            f"dataset-level skips generated program rows: {sorted(leaked)}",
        )
    assay_status = read_tsv(bundle / "passport_assay_status.tsv")
    if any(
        row["result_grain"] == "gene"
        or row["gene_expansion_authorized"].lower() != "false"
        for row in assay_status
    ):
        raise PassportValidationError(
            "MYOJIN_GENE_NEGATIVE_EXPANSION",
            "Myojin status is restricted to assay/class/program non-support",
        )
    return programs, membership, context


def selected_source_rows(
    bundle: Path, fixture_mode: bool
) -> tuple[list[dict[str, str]], str]:
    selection_path = bundle / "passport_input_selection.tsv"
    selected = contract.validate_signed_selection(
        selection_path, allow_fixture=fixture_mode
    )
    if len(selected) != 1:
        raise PassportValidationError(
            "FIXTURE_SELECTION_CARDINALITY", "fixture requires one source artifact"
        )
    row = selected[0]
    artifact = Path(row["artifact_path"])
    if not artifact.is_absolute():
        artifact = bundle / artifact
    return read_tsv(artifact), row["artifact_sha256"]


def check_source_gene_universe(
    genes: pd.DataFrame, source_rows: list[dict[str, str]], fixture_mode: bool
) -> None:
    if not fixture_mode:
        return
    ambiguous = [
        row.get("source_input_row_id", "")
        for row in source_rows
        if any(separator in row.get("ensembl_id", "") for separator in (";", ",", "|"))
    ]
    if ambiguous:
        raise PassportValidationError(
            "MULTI_ENSEMBL_AMBIGUITY",
            f"one source row cannot be split across multiple stable identities: {ambiguous}",
        )
    expected = {row["ensembl_id"].split(".", 1)[0] for row in source_rows}
    observed = set(genes["ensembl_id"].astype(str))
    if expected != observed:
        raise PassportValidationError(
            "SOURCE_GENE_LOSS_OR_GAIN",
            f"missing={sorted(expected - observed)}; unexpected={sorted(observed - expected)}",
        )


def _assert_vocab(
    frame: pd.DataFrame, column: str, values: Mapping[str, str], code: str
) -> None:
    observed = set(frame[column].dropna().astype(str))
    unknown = sorted(observed - set(values))
    if unknown:
        raise PassportValidationError(code, f"{column}: {unknown}")


def check_semantics(evidence: pd.DataFrame) -> None:
    _assert_vocab(evidence, "call_state", contract.CALL_STATES, "UNKNOWN_CALL_STATE")
    _assert_vocab(
        evidence,
        "testability_state",
        contract.TESTABILITY_STATES,
        "UNKNOWN_TESTABILITY",
    )
    _assert_vocab(
        evidence,
        "testability_reason",
        contract.TESTABILITY_STATES,
        "UNKNOWN_TESTABILITY_REASON",
    )
    if "direction" in evidence:
        _assert_vocab(evidence, "direction", contract.DIRECTIONS, "UNKNOWN_DIRECTION")
    if "evidence_domain" in evidence:
        _assert_vocab(
            evidence, "evidence_domain", contract.EVIDENCE_DOMAINS, "UNKNOWN_DOMAIN"
        )

    tested_negative = evidence["call_state"] == "tested_negative"
    if (evidence.loc[tested_negative, "testability_state"] != "testable").any():
        raise PassportValidationError(
            "TESTED_NEGATIVE_REQUIRES_TESTABLE",
            "tested_negative is authorized only inside an explicit adequate test",
        )
    negative_rows = evidence.loc[tested_negative]
    if not negative_rows.empty:
        missing_rule = (
            negative_rows["negative_call_rule_id"].isna()
            | (
                negative_rows["negative_call_rule_id"]
                .astype("string")
                .fillna("")
                .str.strip()
                == ""
            )
            | negative_rows["negative_decision_boundary"].isna()
            | negative_rows["negative_margin"].isna()
            | (negative_rows["negative_margin"] < 0)
            | (negative_rows["negative_call_passed"].fillna(False) != True)  # noqa: E712
        )
        if missing_rule.any():
            raise PassportValidationError(
                "TESTED_NEGATIVE_RULE_REQUIRED",
                "tested_negative requires an explicit rule ID, finite decision boundary, "
                "nonnegative margin, and negative_call_passed=true",
            )
    nonnegative_rows = evidence.loc[~tested_negative]
    if (
        nonnegative_rows["negative_call_rule_id"].notna().any()
        or nonnegative_rows["negative_decision_boundary"].notna().any()
        or nonnegative_rows["negative_margin"].notna().any()
        or nonnegative_rows["negative_call_passed"].fillna(False).any()
    ):
        raise PassportValidationError(
            "NEGATIVE_METADATA_WITHOUT_CALL",
            "informative-negative metadata is reserved for tested_negative rows",
        )
    supported = evidence["call_state"].isin(["supported", "discordant"])
    if (evidence.loc[supported, "testability_state"] != "testable").any():
        raise PassportValidationError(
            "SUPPORTED_REQUIRES_TESTABLE",
            "supported/discordant rows must be explicitly testable",
        )
    no_numeric = evidence["call_state"].isin(["untestable", "not_applicable"])
    if evidence.loc[no_numeric, NUMERIC_EVIDENCE_COLUMNS].notna().any(axis=None):
        raise PassportValidationError(
            "FABRICATED_UNTESTABLE_STATISTIC",
            "untestable/not_applicable rows must not carry effects or p/q",
        )
    not_applicable = evidence["call_state"] == "not_applicable"
    if (
        evidence.loc[not_applicable, "testability_state"] != "assay_out_of_scope"
    ).any():
        raise PassportValidationError(
            "NOT_APPLICABLE_SCOPE", "not_applicable requires assay_out_of_scope"
        )

    if (
        evidence["provenance_state"].isna().any()
        or (evidence["provenance_state"].astype("string").fillna("") == "").any()
    ):
        raise PassportValidationError(
            "UNKNOWN_PROVENANCE", "accepted rows require resolved provenance"
        )
    _assert_vocab(
        evidence, "provenance_state", contract.PROVENANCE_STATES, "UNKNOWN_PROVENANCE"
    )
    expected_dependence = evidence["provenance_state"] != "independent"
    observed_dependence = evidence["source_dependent"].astype(bool)
    if (expected_dependence.to_numpy() != observed_dependence.to_numpy()).any():
        raise PassportValidationError(
            "SOURCE_DEPENDENT_DERIVATION",
            "source_dependent must equal provenance_state != independent",
        )


def check_source_call_preservation(
    evidence: pd.DataFrame, source_rows: list[dict[str, str]], source_artifact_sha: str
) -> None:
    source = {row["source_input_row_id"]: row for row in source_rows}
    if len(source) != len(source_rows):
        raise PassportValidationError(
            "DUPLICATE_SOURCE_ROW_ID", "source input IDs are not unique"
        )
    for record in evidence.to_dict("records"):
        row_id = record["source_input_row_id"]
        if row_id not in source:
            raise PassportValidationError("SOURCE_ROW_UNRESOLVED", str(row_id))
        frozen = source[row_id]
        if record["source_artifact_sha256"] != source_artifact_sha:
            raise PassportValidationError("SOURCE_ARTIFACT_TRACE", str(row_id))
        if (
            record["call_state"] != frozen["source_call_state"]
            or record["source_call_state"] != frozen["source_call_state"]
        ):
            raise PassportValidationError("SOURCE_CALL_MUTATION", str(row_id))
        if (
            record["testability_state"] != frozen["source_testability_state"]
            or record["source_testability_state"] != frozen["source_testability_state"]
        ):
            raise PassportValidationError("SOURCE_TESTABILITY_MUTATION", str(row_id))
        if (
            record["provenance_state"] != frozen["source_provenance_state"]
            or record["source_provenance_state"] != frozen["source_provenance_state"]
        ):
            raise PassportValidationError("SOURCE_PROVENANCE_MUTATION", str(row_id))
        negative_pairs = (
            ("negative_call_rule_id", "source_negative_call_rule_id"),
            ("negative_decision_boundary", "source_negative_decision_boundary"),
            ("negative_margin", "source_negative_margin"),
            ("negative_call_passed", "source_negative_call_passed"),
        )
        for source_field, output_field in negative_pairs:
            frozen_value = frozen[source_field]
            values = (record[source_field], record[output_field])
            for value in values:
                if frozen_value == "" and pd.isna(value):
                    continue
                if source_field in {"negative_decision_boundary", "negative_margin"}:
                    equal = frozen_value != "" and float(value) == float(frozen_value)
                elif source_field == "negative_call_passed":
                    expected_bool = {"True": True, "False": False}.get(frozen_value)
                    equal = expected_bool is not None and bool(value) == expected_bool
                else:
                    equal = str(value) == frozen_value
                if not equal:
                    raise PassportValidationError(
                        "SOURCE_NEGATIVE_RULE_MUTATION", f"{row_id}:{source_field}"
                    )
        expected_hash = contract.source_row_hash(frozen)
        if (
            frozen["source_row_sha256"] != expected_hash
            or record["source_row_sha256"] != expected_hash
        ):
            raise PassportValidationError("SOURCE_ROW_HASH", str(row_id))
        if (
            record["allowed_wording"] != frozen["allowed_wording"]
            or record["limitation"] != frozen["limitation"]
        ):
            raise PassportValidationError("SOURCE_WORDING_MUTATION", str(row_id))


def check_referential_integrity(
    genes: pd.DataFrame, evidence: pd.DataFrame, coverage: pd.DataFrame
) -> None:
    passports = set(genes["passport_id"].astype(str))
    for label, frame in (("evidence", evidence), ("coverage", coverage)):
        unresolved = set(frame["passport_id"].astype(str)) - passports
        if unresolved:
            raise PassportValidationError(
                "GENE_REFERENCE", f"{label}: {sorted(unresolved)}"
            )
        mapping = genes.set_index("passport_id")["ensembl_id"].astype(str).to_dict()
        bad = frame[
            frame.apply(
                lambda row: mapping[str(row["passport_id"])] != str(row["ensembl_id"]),
                axis=1,
            )
        ]
        if not bad.empty:
            raise PassportValidationError("GENE_IDENTITY_MUTATION", label)
    coverage_key = coverage.set_index(["passport_id", "assay", "dataset_id"])
    for row in evidence.to_dict("records"):
        key = (row["passport_id"], row["assay"], row["dataset_id"])
        if key not in coverage_key.index:
            raise PassportValidationError("COVERAGE_ROW_MISSING", str(key))
        cov = coverage_key.loc[key]
        if isinstance(cov, pd.DataFrame):
            raise PassportValidationError("DUPLICATE_COVERAGE_KEY", str(key))
        for column in ("call_state", "testability_state", "testability_reason"):
            if cov[column] != row[column]:
                raise PassportValidationError(
                    "COVERAGE_CALL_MUTATION", f"{key}:{column}"
                )


def _derived_reachable(
    edges: list[dict[str, str]], start: str, node_types: dict[str, str]
) -> set[str]:
    adjacency: dict[str, list[str]] = defaultdict(list)
    for edge in edges:
        if edge["edge_type"] == "derived_from":
            adjacency[edge["from_node_id"]].append(edge["to_node_id"])
    reached: set[str] = set()
    queue = deque([start])
    while queue:
        current = queue.popleft()
        for nxt in adjacency[current]:
            if nxt not in reached:
                reached.add(nxt)
                queue.append(nxt)
    return {node_types[node] for node in reached if node in node_types}


def check_source_graph(bundle: Path, evidence: pd.DataFrame) -> None:
    nodes = read_tsv(bundle / "passport_source_nodes.tsv")
    edges = read_tsv(bundle / "passport_source_edges.tsv")
    node_ids = [row["source_node_id"] for row in nodes]
    if len(node_ids) != len(set(node_ids)):
        raise PassportValidationError("DUPLICATE_SOURCE_NODE", "source_node_id")
    node_types = {row["source_node_id"]: row["node_type"] for row in nodes}
    unknown_types = sorted(set(node_types.values()) - set(contract.NODE_TYPES))
    if unknown_types:
        raise PassportValidationError("UNKNOWN_SOURCE_NODE_TYPE", str(unknown_types))
    edge_ids = [row["source_edge_id"] for row in edges]
    if len(edge_ids) != len(set(edge_ids)):
        raise PassportValidationError("DUPLICATE_SOURCE_EDGE", "source_edge_id")
    for edge in edges:
        if edge["from_node_id"] == edge["to_node_id"]:
            raise PassportValidationError("SOURCE_SELF_EDGE", edge["source_edge_id"])
        if (
            edge["from_node_id"] not in node_types
            or edge["to_node_id"] not in node_types
        ):
            raise PassportValidationError(
                "SOURCE_EDGE_REFERENCE", edge["source_edge_id"]
            )
        if edge["edge_type"] not in contract.EDGE_TYPES:
            raise PassportValidationError("UNKNOWN_SOURCE_EDGE_TYPE", edge["edge_type"])
    for source_node in set(evidence["source_node_id"].astype(str)):
        if source_node not in node_types:
            raise PassportValidationError("EVIDENCE_SOURCE_NODE", source_node)
        reached_types = _derived_reachable(edges, source_node, node_types)
        if not {"dataset", "frozen_release"}.issubset(reached_types):
            raise PassportValidationError(
                "SOURCE_TRACE_INCOMPLETE",
                f"{source_node} reaches {sorted(reached_types)}",
            )
    linked_results = {
        edge["evidence_result_id"] for edge in edges if edge["evidence_result_id"]
    }
    missing = set(evidence["evidence_result_id"].astype(str)) - linked_results
    if missing:
        raise PassportValidationError(
            "SOURCE_RESULT_EDGE_MISSING", str(sorted(missing))
        )


def check_rulebook_and_experiments(bundle: Path, genes: pd.DataFrame) -> None:
    rules = read_tsv(bundle / "passport_rulebook.tsv")
    rule_ids = [row["rule_id"] for row in rules]
    if len(rule_ids) != len(set(rule_ids)):
        raise PassportValidationError("DUPLICATE_RULE_ID", "rulebook rule_id")
    role_rules = {
        row["rule_id"]: row for row in rules if row["rule_type"] == "role_hypothesis"
    }
    experiment_rules = {
        row["rule_id"]: row for row in rules if row["rule_type"] == "next_experiment"
    }
    _assert_vocab(
        genes,
        "primary_evidence_class",
        contract.PRIMARY_EVIDENCE_CLASSES,
        "UNKNOWN_PRIMARY_CLASS",
    )
    _assert_vocab(genes, "role_hypothesis", contract.ROLE_HYPOTHESES, "UNKNOWN_ROLE")
    for row in genes.to_dict("records"):
        expected_role, expected_role_rule, expected_experiment = (
            contract.route_hypothesis(row["primary_evidence_class"])
        )
        if (
            row["role_hypothesis"] != expected_role
            or row["role_rule_id"] != expected_role_rule
            or row["next_experiment_rule_id"] != expected_experiment
        ):
            raise PassportValidationError(
                "ROLE_ROUTE_MUTATION", str(row["passport_id"])
            )
        if (
            expected_role_rule not in role_rules
            or expected_experiment not in experiment_rules
        ):
            raise PassportValidationError("RULE_REFERENCE", str(row["passport_id"]))

    experiments = read_tsv(bundle / "passport_next_experiment.tsv")
    if len(experiments) != len(genes):
        raise PassportValidationError(
            "EXPERIMENT_CARDINALITY", "exactly one experiment per gene"
        )
    counts = Counter(row["passport_id"] for row in experiments)
    if any(value != 1 for value in counts.values()) or set(counts) != set(
        genes["passport_id"]
    ):
        raise PassportValidationError(
            "EXPERIMENT_GENE_REFERENCE", "one experiment per Gene Catalog entry"
        )
    gene_map = genes.set_index("passport_id").to_dict("index")
    for row in experiments:
        gene = gene_map[row["passport_id"]]
        if row["experiment_rule_id"] != gene["next_experiment_rule_id"]:
            raise PassportValidationError(
                "EXPERIMENT_ROUTE_MUTATION", row["passport_id"]
            )
        expected = contract.EXPERIMENT_RULES[row["experiment_rule_id"]]
        for field, value in expected.items():
            if row[field] != value:
                raise PassportValidationError(
                    "EXPERIMENT_TEXT_MUTATION", f"{row['passport_id']}:{field}"
                )
        if row["disclaimer"] != contract.DISCLAIMER:
            raise PassportValidationError("EXPERIMENT_DISCLAIMER", row["passport_id"])


def check_data_dictionary(bundle: Path) -> None:
    dictionary = read_tsv(bundle / "passport_data_dictionary.tsv")
    keys = [(row["table_name"], row["column_name"]) for row in dictionary]
    if len(keys) != len(set(keys)):
        raise PassportValidationError(
            "DUPLICATE_DICTIONARY_KEY", "table_name + column_name"
        )
    dictionary_keys = set(keys)
    for filename in [*contract.PARQUET_SCHEMAS, *contract.TSV_SCHEMAS]:
        if filename == "passport_data_dictionary.tsv":
            continue
        for column in table_columns(bundle / filename):
            if (filename, column) not in dictionary_keys:
                raise PassportValidationError(
                    "DICTIONARY_COLUMN_MISSING", f"{filename}:{column}"
                )
    vocab_rows = read_tsv(bundle / "passport_controlled_vocabularies.tsv")
    observed = defaultdict(set)
    for row in vocab_rows:
        observed[row["vocabulary_name"]].add(row["value"])
    for name, expected in contract.VOCABULARIES.items():
        if observed[name] != set(expected):
            raise PassportValidationError(
                "VOCABULARY_DEFINITION_DRIFT", f"{name}: {sorted(observed[name])}"
            )


def check_production_data_dictionary(bundle: Path) -> None:
    dictionary = read_tsv(bundle / "passport_data_dictionary.tsv")
    keys = [(row["table_name"], row["column_name"]) for row in dictionary]
    if len(keys) != len(set(keys)):
        raise PassportValidationError(
            "DUPLICATE_DICTIONARY_KEY", "table_name + column_name"
        )
    dictionary_keys = set(keys)
    schemas = {
        **contract.PRODUCTION_PARQUET_SCHEMAS,
        **contract.PORTAL_PARQUET_SCHEMAS,
        **contract.PRODUCTION_TSV_SCHEMAS,
    }
    for filename, schema in schemas.items():
        if filename == "passport_data_dictionary.tsv":
            continue
        columns = list(schema)
        for column in columns:
            if (filename, column) not in dictionary_keys:
                raise PassportValidationError(
                    "DICTIONARY_COLUMN_MISSING", f"{filename}:{column}"
                )
    vocab_rows = read_tsv(bundle / "passport_controlled_vocabularies.tsv")
    observed = defaultdict(set)
    for row in vocab_rows:
        observed[row["vocabulary_name"]].add(row["value"])
    for name, expected in contract.VOCABULARIES.items():
        if observed[name] != set(expected):
            raise PassportValidationError(
                "VOCABULARY_DEFINITION_DRIFT", f"{name}: {sorted(observed[name])}"
            )


def check_portal_exports(bundle: Path) -> None:
    genes = pd.read_parquet(bundle / "passport_gene_index.parquet", engine="pyarrow")
    summary = pd.read_parquet(
        bundle / "portal_export/passport_summary.parquet", engine="pyarrow"
    )
    expected_summary = genes.copy()
    expected_summary["alphabetical_key"] = (
        expected_summary["symbol"].astype(str).str.upper()
        + "|"
        + expected_summary["ensembl_id"].astype(str)
    ).astype("string")
    try:
        pd.testing.assert_frame_equal(summary, expected_summary, check_dtype=True)
    except AssertionError as exc:
        raise PassportValidationError(
            "PORTAL_SUMMARY_MUTATION", str(exc)[:500]
        ) from exc
    exports = (
        ("passport_evidence_long.parquet", "portal_export/passport_evidence.parquet"),
        (
            "passport_program_index.parquet",
            "portal_export/passport_program_index.parquet",
        ),
        (
            "passport_program_membership.parquet",
            "portal_export/passport_program_membership.parquet",
        ),
        (
            "passport_program_context.parquet",
            "portal_export/passport_program_context.parquet",
        ),
    )
    for full_name, export_name in exports:
        full = pd.read_parquet(bundle / full_name, engine="pyarrow")
        export = pd.read_parquet(bundle / export_name, engine="pyarrow")
        try:
            pd.testing.assert_frame_equal(full, export, check_dtype=True)
        except AssertionError as exc:
            raise PassportValidationError(
                "PORTAL_EXPORT_MUTATION", export_name
            ) from exc
    status = json.loads(
        (bundle / "portal_export/passport_bundle_status.json").read_text(
            encoding="utf-8"
        )
    )
    if status.get("combined_score_constructed") is not False:
        raise PassportValidationError("FORBIDDEN_SCORE_OR_RANK", "portal status")
    if status.get("scientific_call_recomputed") is not False:
        raise PassportValidationError("PORTAL_RECOMPUTED_CALL", "portal status")
    if status.get("default_order") != "alphabetical":
        raise PassportValidationError(
            "PORTAL_DEFAULT_ORDER", str(status.get("default_order"))
        )
    if status.get("gene_level_domain_scope") != ["genetics", "transcriptomics"]:
        raise PassportValidationError(
            "PORTAL_GENE_DOMAIN_SCOPE", str(status.get("gene_level_domain_scope"))
        )
    if status.get("cross_species_populated") is not False:
        raise PassportValidationError(
            "CROSS_SPECIES_CANDIDATE_BOUNDARY", "portal status"
        )
    graph = json.loads(
        (bundle / "portal_export/passport_source_graph.json").read_text(
            encoding="utf-8"
        )
    )
    if graph.get("nodes") != read_tsv(bundle / "passport_source_nodes.tsv"):
        raise PassportValidationError("PORTAL_SOURCE_GRAPH_MUTATION", "nodes")
    if graph.get("edges") != read_tsv(bundle / "passport_source_edges.tsv"):
        raise PassportValidationError("PORTAL_SOURCE_GRAPH_MUTATION", "edges")
    domain_rows = read_tsv(bundle / "passport_domain_coverage.tsv")
    domain_export = json.loads(
        (bundle / "portal_export/passport_domain_coverage.json").read_text(
            encoding="utf-8"
        )
    ).get("records")
    if not isinstance(domain_export, list):
        raise PassportValidationError(
            "PORTAL_DOMAIN_COVERAGE_SCHEMA", "records must be a JSON array"
        )
    normalized_export = [
        {
            key: ""
            if value is None
            else str(value).lower()
            if isinstance(value, bool)
            else str(value)
            for key, value in row.items()
        }
        for row in domain_export
    ]
    if normalized_export != domain_rows:
        raise PassportValidationError(
            "PORTAL_DOMAIN_COVERAGE_MUTATION",
            "coverage UI must reproduce the frozen TSV",
        )


def _ui_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "t"}:
        return True
    if text in {"false", "0", "no", "f", ""}:
        return False
    raise PassportValidationError("UI_EMBEDDED_BOOLEAN", repr(value))


def _ui_number(value: object) -> int | float | None:
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    number = float(value)
    return int(number) if number.is_integer() else number


def _ui_string(value: object) -> str:
    return "" if value is None or pd.isna(value) else str(value)


def _expected_ui_evidence(row: Mapping[str, object]) -> dict[str, object]:
    return {
        "id": _ui_string(row["evidence_result_id"]),
        "domain": _ui_string(row["evidence_domain"]),
        "assay": _ui_string(row["assay"]),
        "dataset": _ui_string(row["dataset_id"]),
        "phenotype": _ui_string(row["phenotype"]),
        "context": _ui_string(row["context"]),
        "grain": _ui_string(row["biological_unit"]),
        "contrast": _ui_string(row["contrast_or_exposure"]),
        "unit": _ui_string(row["effect_unit"]),
        "estimate": _ui_number(row["estimate"]),
        "q": _ui_number(row["q_value"]),
        "direction": _ui_string(row["direction"]),
        "testability": _ui_string(row["testability_state"]),
        "testabilityReason": _ui_string(row["testability_reason"]),
        "call": _ui_string(row["call_state"]),
        "provenance": _ui_string(row["provenance_state"]),
        "sourceDependent": _ui_bool(row["source_dependent"]),
        "sourceNode": _ui_string(row["source_node_id"]),
        "sourceRelease": _ui_string(row["source_release_id"]),
        "sourceRow": _ui_string(row["source_input_row_id"]),
        "wording": _ui_string(row["allowed_wording"]),
        "limitation": _ui_string(row["limitation"]),
    }


def _expected_ui_program_context(row: Mapping[str, object]) -> dict[str, object]:
    return {
        "id": _ui_string(row["program_context_id"]),
        "programId": _ui_string(row["program_uid"]),
        "program": _ui_string(row["program_label"]),
        "cellType": _ui_string(row["cell_type"]),
        "dataset": _ui_string(row["dataset_id"]),
        "assay": _ui_string(row["assay"]),
        "grain": _ui_string(row["biological_unit"]),
        "contrast": _ui_string(row["contrast_or_exposure"]),
        "unit": _ui_string(row["effect_unit"]),
        "estimate": _ui_number(row["estimate"]),
        "q": _ui_number(row["q_value"]),
        "direction": _ui_string(row["direction"]),
        "testability": _ui_string(row["testability_state"]),
        "testabilityReason": _ui_string(row["testability_reason"]),
        "call": _ui_string(row["call_state"]),
        "provenance": _ui_string(row["provenance_state"]),
        "sourceDependent": _ui_bool(row["source_dependent"]),
        "sourceRelease": _ui_string(row["source_release_id"]),
        "wording": _ui_string(row["allowed_wording"]),
        "limitation": _ui_string(row["limitation"]),
        "geneExpansionAuthorized": _ui_bool(
            row["gene_call_expansion_authorized"]
        ),
    }


def _expected_ui_assay_record(row: Mapping[str, object]) -> dict[str, object]:
    return {
        "id": _ui_string(row["assay_status_id"]),
        "grain": _ui_string(row["result_grain"]),
        "entity": _ui_string(row["entity_id"]),
        "assay": _ui_string(row["assay"]),
        "dataset": _ui_string(row["dataset_id"]),
        "status": _ui_string(row["status"]),
        "call": _ui_string(row["call_state"]),
        "testability": _ui_string(row["testability_state"]),
        "estimate": _ui_number(row["estimate"]),
        "q": _ui_number(row["q_value"]),
        "unit": _ui_string(row["effect_unit"]),
        "provenance": _ui_string(row["provenance_state"]),
        "sourceDependent": _ui_bool(row["source_dependent"]),
        "sourceRelease": _ui_string(row["source_release_id"]),
        "wording": _ui_string(row["allowed_wording"]),
        "limitation": _ui_string(row["limitation"]),
        "geneExpansionAuthorized": _ui_bool(row["gene_expansion_authorized"]),
    }


def check_candidate_ui(
    bundle: Path, *, fixture_mode: bool, sealed: bool = False
) -> None:
    """Verify the local UI against its exact frozen source tables and renders."""

    ui_root = bundle / passport_ui.UI_RELATIVE_ROOT
    paths = {
        "index": ui_root / "index.html",
        "contract": ui_root / "ui_contract.json",
        "sources": ui_root / "ui_source_manifest.tsv",
        "review": ui_root / "review/review_manifest.tsv",
    }
    for label, path in paths.items():
        if not path.is_file():
            raise PassportValidationError("UI_FILE_MISSING", f"{label}:{path}")
        if path.is_symlink():
            raise PassportValidationError("UI_FILE_SYMLINK", f"{label}:{path}")

    source_rows = read_tsv(paths["sources"])
    if table_columns(paths["sources"]) != list(passport_ui.SOURCE_MANIFEST_COLUMNS):
        raise PassportValidationError("UI_SOURCE_MANIFEST_SCHEMA", str(paths["sources"]))
    expected_sources = {
        relative: role for role, relative in passport_ui.SOURCE_TABLES
    }
    indexed_sources = {row["relative_path"]: row for row in source_rows}
    if len(indexed_sources) != len(source_rows) or set(indexed_sources) != set(
        expected_sources
    ):
        raise PassportValidationError(
            "UI_SOURCE_MANIFEST_INVENTORY",
            f"expected={sorted(expected_sources)};observed={sorted(indexed_sources)}",
        )
    for relative, role in expected_sources.items():
        row = indexed_sources[relative]
        if row["source_role"] != role:
            raise PassportValidationError("UI_SOURCE_ROLE", relative)
        clean = PurePosixPath(relative)
        if clean.is_absolute() or ".." in clean.parts:
            raise PassportValidationError("UI_SOURCE_PATH", relative)
        source = bundle / relative
        if not source.is_file() or source.is_symlink():
            raise PassportValidationError("UI_SOURCE_FILE", relative)
        if row["sha256"] != contract.sha256_file(source):
            raise PassportValidationError("UI_SOURCE_HASH", relative)
        if int(row["bytes"]) != source.stat().st_size:
            raise PassportValidationError("UI_SOURCE_BYTES", relative)

    try:
        ui_contract = json.loads(paths["contract"].read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise PassportValidationError("UI_CONTRACT_JSON", str(exc)) from exc
    expected_contract_keys = {
        "contract_version", "analysis_release_id", "alphabetical_default",
        "embedded_payload_version", "embedded_payload_sha256", "fixture_mode",
        "filters", "forbidden_gene_ordering", "hero_symbols", "required_views",
        "source_manifest_sha256", "generator", "generator_sha256", "renderer",
        "renderer_sha256", "render_width_px", "render_height_px",
        "canonical_promotion_authorized", "scientific_call_recomputed",
    }
    if set(ui_contract) != expected_contract_keys:
        raise PassportValidationError(
            "UI_CONTRACT_SCHEMA",
            f"missing={sorted(expected_contract_keys - set(ui_contract))};"
            f"extra={sorted(set(ui_contract) - expected_contract_keys)}",
        )
    genes = pd.read_parquet(bundle / "passport_gene_index.parquet", engine="pyarrow")
    release_ids = set(genes["analysis_release_id"].astype(str))
    if len(release_ids) != 1:
        raise PassportValidationError("UI_RELEASE_CARDINALITY", str(sorted(release_ids)))
    release_id = next(iter(release_ids))
    expected_contract_values = {
        "contract_version": passport_ui.UI_CONTRACT_VERSION,
        "analysis_release_id": release_id,
        "alphabetical_default": True,
        "embedded_payload_version": passport_ui.EMBEDDED_PAYLOAD_VERSION,
        "fixture_mode": fixture_mode,
        "filters": ["search", "call", "testability", "provenance"],
        "forbidden_gene_ordering": ["effect_order", "composite_order"],
        "hero_symbols": list(passport_ui.HERO_SYMBOLS),
        "required_views": ["overview", *passport_ui.HERO_SYMBOLS],
        "source_manifest_sha256": contract.sha256_file(paths["sources"]),
        "generator": "scripts/portal/build_evidence_passport_ui.py",
        "renderer": "scripts/portal/render_evidence_passport_ui.py",
        "render_width_px": 1440,
        "render_height_px": 1200,
        "canonical_promotion_authorized": False,
        "scientific_call_recomputed": False,
    }
    for field, expected in expected_contract_values.items():
        if ui_contract[field] != expected:
            raise PassportValidationError(
                "UI_CONTRACT_SEMANTICS", f"{field}={ui_contract[field]!r}"
            )

    release_manifest = {
        row["relative_path"]: row
        for row in read_tsv(bundle / "passport_release_manifest.tsv")
    }
    generator_hash = str(ui_contract["generator_sha256"])
    renderer_hash = str(ui_contract["renderer_sha256"])
    if sealed:
        for relative in (
            "portal_candidate/index.html",
            "portal_candidate/ui_contract.json",
            "portal_candidate/ui_source_manifest.tsv",
        ):
            row = release_manifest.get(relative)
            if (
                row is None
                or row["producer"] != ui_contract["generator"]
                or row["producer_sha256"] != generator_hash
            ):
                raise PassportValidationError("UI_GENERATOR_PROVENANCE", relative)
        renderer_paths = [
            "portal_candidate/review/overview.png",
            *(
                f"portal_candidate/review/{symbol}.png"
                for symbol in passport_ui.HERO_SYMBOLS
            ),
            "portal_candidate/review/review_manifest.tsv",
        ]
        for relative in renderer_paths:
            row = release_manifest.get(relative)
            if (
                row is None
                or row["producer"] != ui_contract["renderer"]
                or row["producer_sha256"] != renderer_hash
            ):
                raise PassportValidationError("UI_RENDERER_PROVENANCE", relative)
    else:
        if generator_hash != contract.sha256_file(passport_ui.SCRIPT_PATH):
            raise PassportValidationError("UI_GENERATOR_PROVENANCE", generator_hash)
        if renderer_hash != contract.sha256_file(passport_ui_renderer.SCRIPT_PATH):
            raise PassportValidationError("UI_RENDERER_PROVENANCE", renderer_hash)

    index_text = paths["index"].read_text(encoding="utf-8")
    if len(re.findall(r'<script id="passport-data" type="application/json">', index_text)) != 1:
        raise PassportValidationError("UI_EMBEDDED_PAYLOAD_CARDINALITY", str(paths["index"]))
    match = re.search(
        r'<script id="passport-data" type="application/json">(.*?)</script>',
        index_text,
        flags=re.DOTALL,
    )
    if match is None:
        raise PassportValidationError("UI_EMBEDDED_PAYLOAD_MISSING", str(paths["index"]))
    embedded_text = match.group(1)
    embedded_hash = hashlib.sha256(embedded_text.encode("utf-8")).hexdigest()
    if embedded_hash != ui_contract["embedded_payload_sha256"]:
        raise PassportValidationError("UI_EMBEDDED_PAYLOAD_HASH", embedded_hash)
    try:
        payload = json.loads(embedded_text)
    except json.JSONDecodeError as exc:
        raise PassportValidationError("UI_EMBEDDED_PAYLOAD_JSON", str(exc)) from exc
    if re.search(r"<(?:script|img)\b[^>]*\bsrc\s*=", index_text, re.IGNORECASE):
        raise PassportValidationError("UI_EXTERNAL_ASSET", "script/img src")
    if re.search(r"<link\b[^>]*\bhref\s*=", index_text, re.IGNORECASE):
        raise PassportValidationError("UI_EXTERNAL_ASSET", "link href")
    required_visible_tokens = (
        "Source-authorized claim", "Limitation", "Testability reason",
        "Provenance", "Domain boundary", "Next discriminating experiment",
        "Falsifying outcome", "gene grain", "program grain", "assay",
        'data-passport=', 'searchParams.set("passport"',
        'searchParams.set("view",view)', 'p.get("view")',
        'view==="program-context"', 'view==="assay-records"',
        'id="program-context-view"', 'id="assay-record-view"',
        f'data-ui-contract="{passport_ui.UI_CONTRACT_VERSION}"',
    )
    missing_tokens = [token for token in required_visible_tokens if token not in index_text]
    if missing_tokens:
        raise PassportValidationError("UI_VISIBLE_SEMANTICS", str(missing_tokens))
    if "data-symbol=" in index_text or "find(g=>g.symbol===" in index_text:
        raise PassportValidationError(
            "UI_IDENTITY_NAVIGATION", "general navigation must use Catalog entry IDs"
        )
    if re.search(r"https?://", index_text, re.IGNORECASE):
        raise PassportValidationError(
            "UI_EXTERNAL_ASSET", "self-contained candidate HTML must be offline"
        )

    expected_payload_keys = {
        "payloadVersion", "analysisReleaseId", "alphabeticalDefault",
        "fixtureMode", "heroPassports", "grainBoundary", "domains",
        "programContextRecords", "assayGrainRecords", "genes",
    }
    if set(payload) != expected_payload_keys:
        raise PassportValidationError(
            "UI_EMBEDDED_PAYLOAD_SCHEMA",
            f"missing={sorted(expected_payload_keys - set(payload))};"
            f"extra={sorted(set(payload) - expected_payload_keys)}",
        )
    if (
        payload.get("payloadVersion") != passport_ui.EMBEDDED_PAYLOAD_VERSION
        or payload.get("analysisReleaseId") != release_id
        or payload.get("alphabeticalDefault") is not True
        or payload.get("fixtureMode") is not fixture_mode
    ):
        raise PassportValidationError("UI_EMBEDDED_PAYLOAD_SEMANTICS", "header")
    payload_genes = payload.get("genes")
    if not isinstance(payload_genes, list) or len(payload_genes) != len(genes):
        raise PassportValidationError("UI_GENE_UNIVERSE", str(type(payload_genes)))
    ordered_source = genes.assign(
        _symbol_key=genes["symbol"].astype(str).str.casefold(),
        _ensembl_key=genes["ensembl_id"].astype(str),
    ).sort_values(["_symbol_key", "_ensembl_key"], kind="stable")
    observed_ids = [str(row.get("passportId", "")) for row in payload_genes]
    expected_ids = ordered_source["passport_id"].astype(str).tolist()
    if observed_ids != expected_ids or len(observed_ids) != len(set(observed_ids)):
        raise PassportValidationError("UI_GENE_ALPHABETICAL_ORDER", "Catalog entry IDs")

    evidence = pd.read_parquet(
        bundle / "passport_evidence_long.parquet", engine="pyarrow"
    )
    evidence_by_passport: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in evidence.sort_values(
        ["passport_id", "evidence_domain", "assay", "evidence_result_id"],
        kind="stable",
    ).to_dict("records"):
        evidence_by_passport[str(row["passport_id"])].append(
            _expected_ui_evidence(row)
        )
    experiments = {
        row["passport_id"]: row
        for row in read_tsv(bundle / "passport_next_experiment.tsv")
    }
    source_gene_by_id = {
        str(row["passport_id"]): row for row in genes.to_dict("records")
    }
    for row in payload_genes:
        passport_id = str(row.get("passportId", ""))
        source = source_gene_by_id[passport_id]
        expected_evidence = evidence_by_passport.get(passport_id, [])
        expected_gene = {
            "passportId": passport_id,
            "symbol": _ui_string(source["symbol"]),
            "ensembl": _ui_string(source["ensembl_id"]),
            "primaryClass": _ui_string(source["primary_evidence_class"]),
            "role": _ui_string(source["role_hypothesis"]),
            "claim": _ui_string(source["allowed_wording"]),
            "limitation": _ui_string(source["limitation"]),
            "calls": sorted({item["call"] for item in expected_evidence}),
            "testability": sorted(
                {item["testability"] for item in expected_evidence}
            ),
            "provenance": sorted(
                {item["provenance"] for item in expected_evidence}
            ),
            "sourceDependence": sorted(
                {
                    "source_dependent" if item["sourceDependent"] else "independent"
                    for item in expected_evidence
                }
            ),
            "domains": sorted({item["domain"] for item in expected_evidence}),
            "evidence": expected_evidence,
        }
        experiment = experiments[passport_id]
        expected_gene["experiment"] = {
            "model": experiment["biological_model"],
            "context": experiment["context"],
            "perturbation": experiment["perturbation"],
            "readout": experiment["primary_readout"],
            "falsifier": experiment["falsifying_outcome"],
            "disclaimer": experiment["disclaimer"],
        }
        if row != expected_gene:
            raise PassportValidationError("UI_GENE_SEMANTIC_DRIFT", passport_id)

    hero_passports = payload.get("heroPassports")
    if not isinstance(hero_passports, dict) or set(hero_passports) != set(
        passport_ui.HERO_SYMBOLS
    ):
        raise PassportValidationError("UI_HERO_SET", repr(hero_passports))
    for hero_index, symbol in enumerate(passport_ui.HERO_SYMBOLS):
        matches = genes.loc[genes["symbol"].astype(str) == symbol, "passport_id"]
        expected = (
            expected_ids[hero_index]
            if fixture_mode and len(matches) != 1
            else (str(matches.iloc[0]) if len(matches) == 1 else None)
        )
        if (not fixture_mode and len(matches) != 1) or hero_passports[symbol] != expected:
            raise PassportValidationError("UI_HERO_IDENTITY", symbol)

    expected_domains = []
    for row in sorted(
        read_tsv(bundle / "passport_domain_coverage.tsv"),
        key=lambda value: value["evidence_domain"],
    ):
        expected_domains.append(
            {
                "domain": row["evidence_domain"],
                "geneStatus": row["gene_evidence_status"],
                "programStatus": row["program_context_status"],
                "assayStatus": row["assay_status"],
                "boundary": row["coverage_boundary"],
                "wording": row["allowed_wording"],
                "limitation": row["limitation"],
            }
        )
    if payload.get("domains") != expected_domains:
        raise PassportValidationError("UI_DOMAIN_SEMANTIC_DRIFT", "domains")
    program_context = pd.read_parquet(
        bundle / "passport_program_context.parquet", engine="pyarrow"
    )
    expected_program_context = [
        _expected_ui_program_context(row)
        for row in program_context.sort_values(
            ["program_label", "dataset_id", "assay", "program_context_id"],
            kind="stable",
        ).to_dict("records")
    ]
    if payload.get("programContextRecords") != expected_program_context:
        raise PassportValidationError(
            "UI_PROGRAM_CONTEXT_SEMANTIC_DRIFT", "programContextRecords"
        )
    assay_status = pd.read_csv(
        bundle / "passport_assay_status.tsv",
        sep="\t", dtype="string", keep_default_na=False,
    )
    expected_assay_records = [
        _expected_ui_assay_record(row)
        for row in assay_status.sort_values(
            ["result_grain", "entity_id", "assay_status_id"], kind="stable"
        ).to_dict("records")
    ]
    if payload.get("assayGrainRecords") != expected_assay_records:
        raise PassportValidationError(
            "UI_ASSAY_GRAIN_SEMANTIC_DRIFT", "assayGrainRecords"
        )
    grain = payload.get("grainBoundary")
    expected_grain = {
        "gene": "Gene pages show only frozen gene-grain evidence and coverage.",
        "program": (
            "Cell programs and spatial contexts remain at program grain; membership "
            "does not create a gene call."
        ),
        "assay": (
            "Functional challenge records remain at assay, class, or program grain; "
            "they do not create a gene call."
        ),
        "programs": pq.read_table(
            bundle / "passport_program_index.parquet"
        ).num_rows,
        "programContexts": pq.read_table(
            bundle / "passport_program_context.parquet"
        ).num_rows,
        "assayRecords": len(read_tsv(bundle / "passport_assay_status.tsv")),
    }
    if grain != expected_grain:
        raise PassportValidationError("UI_GRAIN_BOUNDARY", repr(grain))

    review_rows = read_tsv(paths["review"])
    if table_columns(paths["review"]) != list(
        passport_ui_renderer.REVIEW_MANIFEST_COLUMNS
    ):
        raise PassportValidationError("UI_REVIEW_MANIFEST_SCHEMA", str(paths["review"]))
    by_review = {row["review_id"]: row for row in review_rows}
    expected_review_ids = {
        "passport_ui_review:overview",
        *(f"passport_ui_review:{symbol}" for symbol in passport_ui.HERO_SYMBOLS),
    }
    if len(by_review) != len(review_rows) or set(by_review) != expected_review_ids:
        raise PassportValidationError("UI_REVIEW_MANIFEST_INVENTORY", str(sorted(by_review)))
    observed_png_hashes = set()
    browser_versions = set()
    browser_hashes = set()
    for view in ("overview", *passport_ui.HERO_SYMBOLS):
        row = by_review[f"passport_ui_review:{view}"]
        gene_symbol = "" if view == "overview" else view
        relative = f"portal_candidate/review/{view}.png"
        expected_values = {
            "relative_path": relative,
            "view": "overview" if view == "overview" else "gene",
            "gene_symbol": gene_symbol,
            "width_px": "1440",
            "height_px": "1200",
            "renderer": ui_contract["renderer"],
            "renderer_sha256": renderer_hash,
            "browser_label": "pinned_puppeteer_chrome",
            "mesa_module": passport_ui_renderer.EXPECTED_MESA_MODULE,
            "source_index_sha256": contract.sha256_file(paths["index"]),
            "ui_contract_sha256": contract.sha256_file(paths["contract"]),
            "render_query": (
                "?review=1" if view == "overview" else f"?gene={view}&review=1"
            ),
            "visible_fields": (
                "gene_grain_boundary;native_grain_navigation"
                if view == "overview"
                else (
                    "claim;limitation;call;testability;testability_reason;"
                    "provenance;source_dependence;grain_boundary;"
                    "next_experiment;falsifier"
                )
            ),
            "viewport_check": "pass",
        }
        for field, expected in expected_values.items():
            if row[field] != expected:
                raise PassportValidationError(
                    "UI_REVIEW_MANIFEST_SEMANTICS", f"{view}:{field}"
                )
        if not row["browser_version"].strip() or not SHA256_PATTERN.fullmatch(
            row["browser_sha256"]
        ):
            raise PassportValidationError("UI_RENDERER_PROVENANCE", view)
        try:
            max_bottom = float(row["max_required_bottom_px"])
        except ValueError as exc:
            raise PassportValidationError(
                "UI_REVIEW_VIEWPORT", f"{view}:nonnumeric"
            ) from exc
        if not 0 < max_bottom <= 1200:
            raise PassportValidationError(
                "UI_REVIEW_VIEWPORT", f"{view}:{max_bottom}"
            )
        image = bundle / relative
        if not image.is_file() or image.is_symlink():
            raise PassportValidationError("UI_REVIEW_IMAGE_MISSING", relative)
        if row["sha256"] != contract.sha256_file(image):
            raise PassportValidationError("UI_REVIEW_IMAGE_HASH", relative)
        if int(row["bytes"]) != image.stat().st_size or image.stat().st_size < 8000:
            raise PassportValidationError("UI_REVIEW_IMAGE_BYTES", relative)
        if passport_ui_renderer.png_dimensions(image) != (1440, 1200):
            raise PassportValidationError("UI_REVIEW_IMAGE_DIMENSIONS", relative)
        observed_png_hashes.add(row["sha256"])
        browser_versions.add(row["browser_version"])
        browser_hashes.add(row["browser_sha256"])
    if len(observed_png_hashes) != 5:
        raise PassportValidationError(
            "UI_REVIEW_IMAGE_DUPLICATION", "each requested view must render distinctly"
        )
    if len(browser_versions) != 1 or len(browser_hashes) != 1:
        raise PassportValidationError(
            "UI_RENDERER_PROVENANCE", "review views used inconsistent browser provenance"
        )
    if sealed:
        import build_evidence_passport_bundle as builder

        expected_environment = builder._review_environment(bundle)  # noqa: SLF001
        for relative in (
            "portal_candidate/review/overview.png",
            *(f"portal_candidate/review/{symbol}.png" for symbol in passport_ui.HERO_SYMBOLS),
            "portal_candidate/review/review_manifest.tsv",
        ):
            row = release_manifest.get(relative)
            if row is None or row["environment"] != expected_environment:
                raise PassportValidationError(
                    "UI_RENDERER_ENVIRONMENT", relative
                )


def check_domain_coverage(bundle: Path, fixture_mode: bool) -> None:
    rows = read_tsv(bundle / "passport_domain_coverage.tsv")
    by_domain = {row["evidence_domain"]: row for row in rows}
    if len(rows) != len(by_domain) or set(by_domain) != set(contract.EVIDENCE_DOMAINS):
        raise PassportValidationError(
            "DOMAIN_COVERAGE_FAMILY",
            f"expected={list(contract.EVIDENCE_DOMAINS)};observed={sorted(by_domain)}",
        )
    count_fields = (
        "n_gene_evidence_rows",
        "n_programs",
        "n_program_membership_rows",
        "n_program_context_rows",
        "n_assay_status_rows",
    )
    for domain, row in by_domain.items():
        try:
            counts = {field: int(row[field]) for field in count_fields}
        except ValueError as exc:
            raise PassportValidationError(
                "DOMAIN_COVERAGE_COUNT", f"{domain}: non-integer count"
            ) from exc
        if any(value < 0 for value in counts.values()):
            raise PassportValidationError(
                "DOMAIN_COVERAGE_COUNT", f"{domain}: negative count"
            )
        gene_populated = counts["n_gene_evidence_rows"] > 0
        if row["gene_evidence_status"] != (
            "populated" if gene_populated else "not_populated"
        ):
            raise PassportValidationError("DOMAIN_GENE_STATUS", domain)
        if row["gene_evidence_authorized"] != str(gene_populated).lower():
            raise PassportValidationError("DOMAIN_GENE_AUTHORIZATION", domain)
        if row["member_gene_expansion_authorized"] != "false":
            raise PassportValidationError("PROGRAM_GENE_EXPANSION", domain)
        if gene_populated and domain not in {"genetics", "transcriptomics"}:
            raise PassportValidationError("GENE_DOMAIN_OUT_OF_CANDIDATE_SCOPE", domain)
    for domain in (
        "spatial",
        "proteomics",
        "chromatin",
        "functional",
        "cross_species",
        "cell_context",
    ):
        if by_domain[domain]["gene_evidence_status"] != "not_populated":
            raise PassportValidationError("GENE_DOMAIN_OUT_OF_CANDIDATE_SCOPE", domain)
    cross_species = by_domain["cross_species"]
    if any(int(cross_species[field]) for field in count_fields):
        raise PassportValidationError(
            "CROSS_SPECIES_CANDIDATE_BOUNDARY",
            "cross-species evidence is not populated in this candidate",
        )
    identity_rows = read_tsv(bundle / "passport_gen_identity_status.tsv")
    if fixture_mode:
        if len(identity_rows) > 1:
            raise PassportValidationError(
                "GEN_IDENTITY_STATUS_CARDINALITY", str(len(identity_rows))
            )
    elif len(identity_rows) != 1:
        raise PassportValidationError(
            "GEN_IDENTITY_STATUS_CARDINALITY",
            "production requires one adjudication terminal",
        )
    if identity_rows:
        identity = identity_rows[0]
        if (
            identity["outcome_fields_used"] != "false"
            or identity["merge_discordant_calls"] != "false"
            or identity["one_row_per_ensembl"] != "true"
            or identity["canonical_promotion_authorized"] != "false"
        ):
            raise PassportValidationError(
                "GEN_IDENTITY_STATUS_SEMANTICS",
                "outcome-blind identity boundary failed",
            )


def check_gate_status(bundle: Path) -> None:
    rows = {
        row["gate_id"]: row for row in read_tsv(bundle / "passport_gate_status.tsv")
    }
    if rows.get("PASS01_SCHEMA_FIXTURE", {}).get("status") != "pass":
        raise PassportValidationError("PASS01_GATE", "schema fixture must pass")
    if any(row["promotion_allowed"].lower() != "false" for row in rows.values()):
        raise PassportValidationError(
            "PREMATURE_PROMOTION", "PASS-01 fixture cannot be promoted"
        )
    if rows.get("PORTAL_PRODUCTION_READY", {}).get("status") != "blocked_pass01_only":
        raise PassportValidationError(
            "PRODUCTION_GATE", "real portal population must remain blocked"
        )


def check_release_manifest(bundle: Path) -> None:
    manifest_path = bundle / "passport_release_manifest.tsv"
    rows = read_tsv(manifest_path)
    indexed = {row["relative_path"]: row for row in rows}
    if len(indexed) != len(rows):
        raise PassportValidationError("DUPLICATE_MANIFEST_PATH", "release manifest")
    actual = {
        path.name
        for path in bundle.iterdir()
        if path.is_file() and path != manifest_path
    }
    if actual != set(indexed):
        raise PassportValidationError(
            "MANIFEST_INVENTORY",
            f"missing={sorted(actual - set(indexed))}; stale={sorted(set(indexed) - actual)}",
        )
    for name, row in indexed.items():
        path = bundle / name
        if contract.sha256_file(path) != row["sha256"]:
            raise PassportValidationError("MANIFEST_HASH", name)
        if path.stat().st_size != int(row["bytes"]):
            raise PassportValidationError("MANIFEST_BYTES", name)
        if row["producer_sha256"] != contract.sha256_file(contract.SCRIPT_PATH):
            raise PassportValidationError("MANIFEST_PRODUCER_HASH", name)
        if row["release_status"] != "synthetic_fixture_not_for_promotion":
            raise PassportValidationError("MANIFEST_RELEASE_STATUS", name)


def check_production_manifest(bundle: Path, allowed_status: str) -> None:
    import build_evidence_passport_bundle as builder

    manifest_path = bundle / "passport_release_manifest.tsv"
    temp_residue = sorted(
        path.relative_to(bundle).as_posix()
        for path in bundle.rglob("*")
        if path.is_file()
        and path.name.startswith(".")
        and path.name.endswith(".tmp")
        and "fixtures" not in path.relative_to(bundle).parts
        and "superseded_bundles" not in path.relative_to(bundle).parts
    )
    if temp_residue:
        raise PassportValidationError("MANIFEST_TEMP_RESIDUE", ", ".join(temp_residue))
    rows = read_tsv(manifest_path)
    indexed = {row["relative_path"]: row for row in rows}
    if len(indexed) != len(rows):
        raise PassportValidationError("DUPLICATE_MANIFEST_PATH", "release manifest")
    exclusions = {
        "passport_release_manifest.tsv",
        "passport_validation_report.tsv",
        "passport_plan60_handoff.tsv",
        TERMINAL_PROVENANCE_FILENAME,
        "PASS06_VALIDATED",
        "FIXTURE_PASS06_VALIDATED",
        "plan60_terminal_artifacts.tsv",
        "plan60_terminal_artifacts.signature.json",
    }
    actual = {
        path.relative_to(bundle).as_posix()
        for path in bundle.rglob("*")
        if path.is_file()
        and path.relative_to(bundle).as_posix() not in exclusions
        and "fixtures" not in path.relative_to(bundle).parts
        and "superseded_bundles" not in path.relative_to(bundle).parts
    }
    if actual != set(indexed):
        raise PassportValidationError(
            "MANIFEST_INVENTORY",
            f"missing={sorted(actual - set(indexed))}; stale={sorted(set(indexed) - actual)}",
        )
    producer_paths = builder.manifest_producer_paths()
    for name, row in indexed.items():
        path = bundle / name
        if contract.sha256_file(path) != row["sha256"]:
            raise PassportValidationError("MANIFEST_HASH", name)
        if path.stat().st_size != int(row["bytes"]):
            raise PassportValidationError("MANIFEST_BYTES", name)
        producer = producer_paths.get(row["producer"])
        if producer is None:
            raise PassportValidationError(
                "MANIFEST_PRODUCER_NOT_ALLOWED", f"{name}:{row['producer']}"
            )
        if not producer.is_file() or producer.is_symlink():
            raise PassportValidationError(
                "MANIFEST_PRODUCER_MISSING", f"{name}:{producer}"
            )
        if row["producer_sha256"] != contract.sha256_file(producer):
            raise PassportValidationError(
                "MANIFEST_PRODUCER_HASH", f"{name}:{row['producer']}"
            )
        if row["release_status"] != allowed_status:
            raise PassportValidationError("MANIFEST_RELEASE_STATUS", name)


def _require_release_relative_path(value: str, error_code: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if (
        not value
        or path.is_absolute()
        or value.startswith(("/", "~"))
        or ".." in path.parts
        or str(path) in {"", "."}
    ):
        raise PassportValidationError(error_code, value)
    return path


def _require_unsymlinked_bundle_root(bundle: Path) -> Path:
    """Return a resolved directory after rejecting symlinks in its lexical path."""

    lexical = Path(os.path.abspath(os.fspath(bundle)))
    current = Path(lexical.anchor)
    for part in lexical.parts[1:]:
        current /= part
        try:
            mode = current.lstat().st_mode
        except FileNotFoundError as exc:
            raise PassportValidationError(
                "TERMINAL_BUNDLE_MISSING", str(current)
            ) from exc
        if stat.S_ISLNK(mode):
            raise PassportValidationError("TERMINAL_BUNDLE_ROOT_SYMLINK", str(current))
    if not stat.S_ISDIR(lexical.lstat().st_mode):
        raise PassportValidationError("TERMINAL_BUNDLE_NOT_DIRECTORY", str(lexical))
    return lexical.resolve(strict=True)


def _require_contained_regular_file(
    bundle: Path,
    relative_path: str,
    *,
    symlink_code: str,
) -> Path:
    """Require a regular, non-symlinked artifact contained by ``bundle``."""

    relative = _require_release_relative_path(
        relative_path, "TERMINAL_ARTIFACT_UNSAFE_PATH"
    )
    artifact = bundle.joinpath(*relative.parts)
    current = bundle
    for part in relative.parts:
        current /= part
        try:
            mode = current.lstat().st_mode
        except FileNotFoundError as exc:
            raise PassportValidationError(
                "TERMINAL_ARTIFACT_MISSING", relative_path
            ) from exc
        if stat.S_ISLNK(mode):
            raise PassportValidationError(symlink_code, relative_path)
    if not stat.S_ISREG(artifact.lstat().st_mode):
        raise PassportValidationError("TERMINAL_ARTIFACT_NOT_REGULAR", relative_path)
    resolved_artifact = artifact.resolve(strict=True)
    try:
        resolved_artifact.relative_to(bundle)
    except ValueError as exc:
        raise PassportValidationError(
            "TERMINAL_ARTIFACT_OUTSIDE_BUNDLE", relative_path
        ) from exc
    return resolved_artifact


def _reject_active_bundle_symlinks(bundle: Path, payload_paths: set[str]) -> None:
    """Reject symlinks anywhere in the active sealed-bundle namespace."""

    stack = [bundle]
    while stack:
        directory = stack.pop()
        with os.scandir(directory) as entries:
            for entry in entries:
                path = Path(entry.path)
                relative = path.relative_to(bundle).as_posix()
                if entry.is_symlink():
                    if relative in TERMINAL_CHAIN_PATHS:
                        code = "TERMINAL_ARTIFACT_SYMLINK"
                    elif relative in payload_paths:
                        code = "TERMINAL_PAYLOAD_SYMLINK"
                    else:
                        code = "TERMINAL_INTERNAL_SYMLINK"
                    raise PassportValidationError(code, relative)
                if entry.is_dir(follow_symlinks=False):
                    if path.parent == bundle and path.name in POST_SEAL_EXCLUDED_TREES:
                        continue
                    stack.append(path)


def _verify_relocatable_payload_manifest(
    bundle: Path, manifest_path: Path, release_status: str, selection_sha256: str
) -> None:
    rows = read_tsv(manifest_path)
    if (
        table_columns(manifest_path)
        != contract.TSV_SCHEMAS["passport_release_manifest.tsv"]
    ):
        raise PassportValidationError("TERMINAL_MANIFEST_SCHEMA", str(manifest_path))
    indexed = {row["relative_path"]: row for row in rows}
    if len(indexed) != len(rows):
        raise PassportValidationError(
            "TERMINAL_MANIFEST_DUPLICATE_PATH", str(manifest_path)
        )
    exclusions = {
        "passport_release_manifest.tsv",
        "passport_validation_report.tsv",
        "passport_plan60_handoff.tsv",
        TERMINAL_PROVENANCE_FILENAME,
        "PASS06_VALIDATED",
        "FIXTURE_PASS06_VALIDATED",
        "plan60_terminal_artifacts.tsv",
        "plan60_terminal_artifacts.signature.json",
    }
    actual = {
        path.relative_to(bundle).as_posix()
        for path in bundle.rglob("*")
        if path.is_file()
        and path.relative_to(bundle).as_posix() not in exclusions
        and "fixtures" not in path.relative_to(bundle).parts
        and "superseded_bundles" not in path.relative_to(bundle).parts
    }
    if actual != set(indexed):
        raise PassportValidationError(
            "TERMINAL_MANIFEST_INVENTORY",
            f"missing={sorted(actual - set(indexed))}; stale={sorted(set(indexed) - actual)}",
        )
    for relative_path, row in indexed.items():
        _require_release_relative_path(relative_path, "TERMINAL_MANIFEST_ABSOLUTE_PATH")
        artifact = _require_contained_regular_file(
            bundle, relative_path, symlink_code="TERMINAL_PAYLOAD_SYMLINK"
        )
        if not SHA256_PATTERN.fullmatch(row["sha256"]):
            raise PassportValidationError(
                "TERMINAL_MANIFEST_INVALID_HASH", relative_path
            )
        if contract.sha256_file(artifact) != row["sha256"]:
            raise PassportValidationError("TERMINAL_MANIFEST_HASH", relative_path)
        if artifact.stat().st_size != int(row["bytes"]):
            raise PassportValidationError("TERMINAL_MANIFEST_BYTES", relative_path)
        _require_release_relative_path(
            row["producer"], "TERMINAL_MANIFEST_ABSOLUTE_PRODUCER"
        )
        if not SHA256_PATTERN.fullmatch(row["producer_sha256"]):
            raise PassportValidationError(
                "TERMINAL_MANIFEST_INVALID_PRODUCER_HASH", relative_path
            )
        if row["upstream_artifact_sha256"] != selection_sha256:
            raise PassportValidationError(
                "TERMINAL_MANIFEST_SELECTION_HASH", relative_path
            )
        if row["release_status"] != release_status:
            raise PassportValidationError(
                "TERMINAL_MANIFEST_RELEASE_STATUS", relative_path
            )


def verify_sealed_bundle(bundle: Path) -> list[dict[str, str]]:
    """Verify a sealed PASS06 bundle without writing or resolving live inputs."""

    bundle = _require_unsymlinked_bundle_root(bundle)
    for relative_path in TERMINAL_CHAIN_PATHS:
        artifact = bundle / relative_path
        if artifact.is_symlink():
            raise PassportValidationError("TERMINAL_ARTIFACT_SYMLINK", relative_path)
    manifest = _require_contained_regular_file(
        bundle,
        "passport_release_manifest.tsv",
        symlink_code="TERMINAL_ARTIFACT_SYMLINK",
    )
    manifest_rows = read_tsv(manifest)
    payload_paths = {row.get("relative_path", "") for row in manifest_rows}
    _reject_active_bundle_symlinks(bundle, payload_paths)
    temp_residue = sorted(
        path.relative_to(bundle).as_posix()
        for path in bundle.rglob("*")
        if path.is_file()
        and path.name.startswith(".")
        and path.name.endswith(".tmp")
        and "fixtures" not in path.relative_to(bundle).parts
        and "superseded_bundles" not in path.relative_to(bundle).parts
    )
    if temp_residue:
        raise PassportValidationError("TERMINAL_TEMP_RESIDUE", ", ".join(temp_residue))
    seal_candidates = [
        name
        for name in ("PASS06_VALIDATED", "FIXTURE_PASS06_VALIDATED")
        if (bundle / name).is_file()
    ]
    if len(seal_candidates) != 1:
        raise PassportValidationError(
            "TERMINAL_SEAL_CARDINALITY", f"observed={seal_candidates}"
        )
    seal_name = seal_candidates[0]
    seal_path = _require_contained_regular_file(
        bundle, seal_name, symlink_code="TERMINAL_ARTIFACT_SYMLINK"
    )
    if table_columns(seal_path) != contract.PASS06_TERMINAL_COLUMNS:
        raise PassportValidationError("TERMINAL_SEAL_SCHEMA", str(seal_path))
    seal_rows = read_tsv(seal_path)
    if len(seal_rows) != 1:
        raise PassportValidationError("TERMINAL_SEAL_CARDINALITY", str(len(seal_rows)))
    seal = seal_rows[0]
    expected_status = (
        "validated_synthetic_pass02_06_fixture"
        if seal_name == "FIXTURE_PASS06_VALIDATED"
        else "passport_release_validated"
    )
    if seal["status"] != expected_status:
        raise PassportValidationError("TERMINAL_SEAL_STATUS", seal["status"])
    expected_booleans = {
        "automated_validation": "pass",
        "canonical_promotion_authorized": "false",
        "scientific_call_recomputed": "false",
    }
    for field, expected in expected_booleans.items():
        if seal[field] != expected:
            raise PassportValidationError(
                "TERMINAL_SEAL_SEMANTICS", f"{field}={seal[field]}"
            )
    for field in ("selection_sha256", "manifest_sha256", "validation_report_sha256"):
        if not SHA256_PATTERN.fullmatch(seal[field]):
            raise PassportValidationError("TERMINAL_SEAL_INVALID_HASH", field)

    selection = _require_contained_regular_file(
        bundle,
        "passport_input_selection.tsv",
        symlink_code="TERMINAL_ARTIFACT_SYMLINK",
    )
    validation_report = _require_contained_regular_file(
        bundle,
        "passport_validation_report.tsv",
        symlink_code="TERMINAL_ARTIFACT_SYMLINK",
    )
    handoff = _require_contained_regular_file(
        bundle,
        "passport_plan60_handoff.tsv",
        symlink_code="TERMINAL_ARTIFACT_SYMLINK",
    )
    terminal_provenance = _require_contained_regular_file(
        bundle,
        TERMINAL_PROVENANCE_FILENAME,
        symlink_code="TERMINAL_ARTIFACT_SYMLINK",
    )
    expected_hashes = {
        selection: seal["selection_sha256"],
        manifest: seal["manifest_sha256"],
        validation_report: seal["validation_report_sha256"],
    }
    for path, expected in expected_hashes.items():
        if contract.sha256_file(path) != expected:
            raise PassportValidationError(
                "TERMINAL_SEAL_HASH", f"{path.name}:{expected}"
            )

    if (
        table_columns(validation_report)
        != contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_validation_report.tsv"]
    ):
        raise PassportValidationError(
            "TERMINAL_VALIDATION_REPORT_SCHEMA", str(validation_report)
        )
    report_rows = read_tsv(validation_report)
    terminal_report_rows = [
        row for row in report_rows if row["check_id"] == TERMINAL_PROVENANCE_CHECK_ID
    ]
    if len(terminal_report_rows) != 1 or terminal_report_rows[0]["status"] != "pass":
        raise PassportValidationError(
            "TERMINAL_PROVENANCE_REPORT_CARDINALITY", str(terminal_report_rows)
        )
    try:
        terminal_detail = json.loads(terminal_report_rows[0]["detail"])
    except json.JSONDecodeError as exc:
        raise PassportValidationError(
            "TERMINAL_PROVENANCE_REPORT_JSON", str(exc)
        ) from exc
    expected_detail_keys = {
        "bundle_uri",
        "contract_version",
        "terminal_provenance_path",
        "terminal_provenance_sha256",
    }
    if set(terminal_detail) != expected_detail_keys:
        raise PassportValidationError(
            "TERMINAL_PROVENANCE_REPORT_SCHEMA", str(sorted(terminal_detail))
        )
    bundle_uri = release_bundle_uri(seal["analysis_release_id"])
    if terminal_detail["bundle_uri"] != bundle_uri:
        raise PassportValidationError(
            "TERMINAL_BUNDLE_URI", str(terminal_detail["bundle_uri"])
        )
    if terminal_detail["contract_version"] != "passport_terminal_provenance_v1":
        raise PassportValidationError(
            "TERMINAL_PROVENANCE_VERSION", str(terminal_detail["contract_version"])
        )
    if terminal_detail["terminal_provenance_path"] != TERMINAL_PROVENANCE_FILENAME:
        raise PassportValidationError(
            "TERMINAL_PROVENANCE_PATH", str(terminal_detail["terminal_provenance_path"])
        )
    if (
        not SHA256_PATTERN.fullmatch(str(terminal_detail["terminal_provenance_sha256"]))
        or contract.sha256_file(terminal_provenance)
        != terminal_detail["terminal_provenance_sha256"]
    ):
        raise PassportValidationError(
            "TERMINAL_PROVENANCE_HASH", TERMINAL_PROVENANCE_FILENAME
        )

    if (
        table_columns(terminal_provenance)
        != contract.PASS06_TERMINAL_PROVENANCE_COLUMNS
    ):
        raise PassportValidationError(
            "TERMINAL_PROVENANCE_SCHEMA", str(terminal_provenance)
        )
    provenance_rows = read_tsv(terminal_provenance)
    indexed_provenance = {row["artifact_role"]: row for row in provenance_rows}
    if set(indexed_provenance) != {"payload_manifest", "plan60_handoff"}:
        raise PassportValidationError(
            "TERMINAL_PROVENANCE_ROLES", str(sorted(indexed_provenance))
        )
    expected_paths = {
        "payload_manifest": "passport_release_manifest.tsv",
        "plan60_handoff": "passport_plan60_handoff.tsv",
    }
    for role, relative_path in expected_paths.items():
        row = indexed_provenance[role]
        if row["analysis_release_id"] != seal["analysis_release_id"]:
            raise PassportValidationError("TERMINAL_PROVENANCE_RELEASE", role)
        if row["bundle_uri"] != bundle_uri or row["relative_path"] != relative_path:
            raise PassportValidationError("TERMINAL_PROVENANCE_LOCATION", role)
        artifact = _require_contained_regular_file(
            bundle, relative_path, symlink_code="TERMINAL_ARTIFACT_SYMLINK"
        )
        if row["sha256"] != contract.sha256_file(artifact):
            raise PassportValidationError("TERMINAL_PROVENANCE_ARTIFACT_HASH", role)
        if int(row["bytes"]) != artifact.stat().st_size:
            raise PassportValidationError("TERMINAL_PROVENANCE_ARTIFACT_BYTES", role)
        _require_release_relative_path(
            row["producer"], "TERMINAL_PROVENANCE_ABSOLUTE_PRODUCER"
        )
        if not SHA256_PATTERN.fullmatch(row["producer_sha256"]):
            raise PassportValidationError("TERMINAL_PROVENANCE_PRODUCER_HASH", role)

    handoff_rows = read_tsv(handoff)
    if (
        table_columns(handoff)
        != contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_plan60_handoff.tsv"]
        or len(handoff_rows) != 1
    ):
        raise PassportValidationError("TERMINAL_HANDOFF_SCHEMA", str(handoff))
    handoff_row = handoff_rows[0]
    expected_handoff = {
        "analysis_release_id": seal["analysis_release_id"],
        "bundle_path": bundle_uri,
        "manifest_sha256": seal["manifest_sha256"],
        "selection_sha256": seal["selection_sha256"],
        "automated_validation": "pass",
        "manual_acceptance": seal["manual_acceptance"],
        "promotion_allowed": "false",
        "scientific_call_recomputed": "false",
    }
    for field, expected in expected_handoff.items():
        if handoff_row[field] != expected:
            raise PassportValidationError(
                "TERMINAL_HANDOFF_SEMANTICS", f"{field}={handoff_row[field]}"
            )

    release_status = (
        "synthetic_fixture_not_for_promotion"
        if seal_name == "FIXTURE_PASS06_VALIDATED"
        else "validated_candidate_handoff"
    )
    _verify_relocatable_payload_manifest(
        bundle, manifest, release_status, seal["selection_sha256"]
    )
    check_candidate_ui(
        bundle,
        fixture_mode=seal_name == "FIXTURE_PASS06_VALIDATED",
        sealed=True,
    )
    if seal_name == "PASS06_VALIDATED":
        validate_manual_acceptance(
            bundle / "passport_manual_acceptance.tsv",
            seal["analysis_release_id"],
            seal["selection_sha256"],
            bundle,
        )

    genes = pq.read_table(bundle / "passport_gene_index.parquet").num_rows
    evidence = pq.read_table(bundle / "passport_evidence_long.parquet").num_rows
    programs = pq.read_table(bundle / "passport_program_index.parquet").num_rows
    contexts = pq.read_table(bundle / "passport_program_context.parquet").num_rows
    observed_counts = {
        "n_genes": genes,
        "n_evidence_rows": evidence,
        "n_programs": programs,
        "n_program_context_rows": contexts,
    }
    for field, observed in observed_counts.items():
        if int(seal[field]) != observed:
            raise PassportValidationError(
                "TERMINAL_SEAL_COUNT", f"{field}:{seal[field]}!={observed}"
            )
    return [
        {
            "check_id": "PASS06_POST_SEAL_READ_ONLY",
            "status": "pass",
            "detail": (
                f"release={seal['analysis_release_id']};bundle_uri={bundle_uri};"
                f"payload_artifacts={len(read_tsv(manifest))};terminal_artifacts=2"
            ),
        }
    ]


def validate_bundle(bundle: Path, fixture_mode: bool = True) -> list[dict[str, str]]:
    bundle = bundle.resolve()
    if fixture_mode:
        contract.assert_fixture_path(bundle, contract.DEFAULT_CANDIDATE_ROOT.resolve())
    check_inventory(bundle)
    check_forbidden_columns(bundle)
    check_schemas(bundle)
    genes, evidence, coverage = check_identity(bundle)
    source_rows, source_artifact_sha = selected_source_rows(bundle, fixture_mode)
    check_source_gene_universe(genes, source_rows, fixture_mode)
    check_semantics(evidence)
    check_source_call_preservation(evidence, source_rows, source_artifact_sha)
    check_referential_integrity(genes, evidence, coverage)
    check_source_graph(bundle, evidence)
    check_rulebook_and_experiments(bundle, genes)
    check_data_dictionary(bundle)
    check_gate_status(bundle)
    check_release_manifest(bundle)
    return [
        {"check": "bundle_validation", "status": "pass", "detail": bundle.name},
        {"check": "gene_rows", "status": "pass", "detail": str(len(genes))},
        {"check": "evidence_rows", "status": "pass", "detail": str(len(evidence))},
        {"check": "coverage_rows", "status": "pass", "detail": str(len(coverage))},
        {
            "check": "score_rank_columns",
            "status": "pass",
            "detail": "no combined score, modality count, posterior, hidden tier, priority, or rank",
        },
    ]


def validate_production_bundle(
    bundle: Path, fixture_mode: bool = False
) -> tuple[list[dict[str, str]], object, object, object]:
    bundle = bundle.resolve()
    if fixture_mode:
        contract.assert_fixture_path(bundle, contract.DEFAULT_CANDIDATE_ROOT.resolve())
    elif bundle != contract.DEFAULT_CANDIDATE_ROOT.resolve():
        raise PassportValidationError(
            "OUTPUT_ROOT",
            f"production bundle must be {contract.DEFAULT_CANDIDATE_ROOT.resolve()}",
        )
    check_production_inventory(bundle)
    check_production_forbidden_columns(bundle)
    check_production_schemas(bundle)
    genes, evidence, coverage = check_identity(bundle)
    check_semantics(evidence)
    check_referential_integrity(genes, evidence, coverage)
    programs, membership, context = check_program_boundaries(bundle)
    assembly, selected, signature = rederive_production_bundle(bundle, fixture_mode)
    check_source_graph(bundle, evidence)
    linked_results = {
        edge["evidence_result_id"]
        for edge in read_tsv(bundle / "passport_source_edges.tsv")
        if edge["evidence_result_id"]
    }
    missing_program_edges = (
        set(context["program_context_id"].astype(str)) - linked_results
    )
    if missing_program_edges:
        raise PassportValidationError(
            "SOURCE_RESULT_EDGE_MISSING",
            f"program_context={sorted(missing_program_edges)}",
        )
    check_rulebook_and_experiments(bundle, genes)
    check_production_data_dictionary(bundle)
    check_domain_coverage(bundle, fixture_mode)
    check_portal_exports(bundle)
    check_candidate_ui(bundle, fixture_mode=fixture_mode)
    check_production_manifest(
        bundle,
        "synthetic_fixture_not_for_promotion"
        if fixture_mode
        else "assembled_candidate_prevalidation",
    )
    selection_sha = contract.sha256_file(bundle / "passport_input_selection.tsv")
    audits = read_tsv(bundle / "accepted_input_audit.tsv")
    if any(row["selection_sha256"] != selection_sha for row in audits):
        raise PassportValidationError(
            "ACCEPTED_INPUT_SELECTION_HASH", "accepted input audit"
        )
    quarantine_rows = read_tsv(bundle / "passport_identity_quarantine.tsv")
    permitted_audit_only_quarantine = {
        "PROGRAM_MEMBERSHIP_UNMAPPED_SOURCE",
        "GEN_IRREDUCIBLE_MULTI_ENSEMBL",
        "GEN_ENSEMBL_ABSENT_V49",
        "GEN_IDENTITY_CONFLICT",
        "GEN_ENSEMBL_UNRESOLVED",
        "GEN_INVALID_ENSEMBL_TOKEN",
    }
    fatal_quarantine = [
        row
        for row in quarantine_rows
        if row["reason_code"] not in permitted_audit_only_quarantine
    ]
    if fatal_quarantine:
        raise PassportValidationError(
            "IDENTITY_QUARANTINE_NONEMPTY",
            "selected gene evidence contains unresolved identity rows outside the explicit source-exclusion audit",
        )
    return (
        [
            {
                "check_id": "PASS00_SIGNED_SELECTION",
                "status": "pass",
                "detail": f"{len(selected)} accepted immutable inputs",
            },
            {
                "check_id": "PASS02_SOURCE_GRAPH",
                "status": "pass",
                "detail": "source graph rederived and result-linked",
            },
            {
                "check_id": "PASS03_GENE_TABLES",
                "status": "pass",
                "detail": f"genes={len(genes)};evidence={len(evidence)};coverage={len(coverage)}",
            },
            {
                "check_id": "PASS03_PROGRAM_BOUNDARY",
                "status": "pass",
                "detail": f"programs={len(programs)};membership={len(membership)};context={len(context)};gene_expansion=false",
            },
            {
                "check_id": "PASS03_IDENTITY_QUARANTINE_BOUNDARY",
                "status": "pass",
                "detail": f"audit_only_quarantine_rows={len(quarantine_rows)};fatal_rows=0",
            },
            {
                "check_id": "PASS03_DOMAIN_COVERAGE_BOUNDARY",
                "status": "pass",
                "detail": "gene evidence restricted to genetics/transcriptomics; program/assay grains remain non-expanding; cross-species absent",
            },
            {
                "check_id": "PASS04_DETERMINISTIC_EXPERIMENTS",
                "status": "pass",
                "detail": "one source-independent deterministic route per gene",
            },
            {
                "check_id": "PASS05_PORTAL_EXPORT",
                "status": "pass",
                "detail": (
                    "parquet/JSON exports and the offline review UI exactly reproduce "
                    "candidate tables; five review images verified"
                ),
            },
            {
                "check_id": "TESTED_NEGATIVE_SEMANTICS",
                "status": "pass",
                "detail": f"explicit_tested_negative_rows={(evidence['call_state'] == 'tested_negative').sum() + (context['call_state'] == 'tested_negative').sum()}; zero is permitted",
            },
            {
                "check_id": "NO_SCORE_OR_RANK",
                "status": "pass",
                "detail": "no combined score, modality count, tier, priority, posterior, or rank",
            },
            {
                "check_id": "ADAPTER_REDERIVATION",
                "status": "pass",
                "detail": "all evidence and program tables rederived from signed selected artifacts",
            },
        ],
        assembly,
        selected,
        signature,
    )


def validate_manual_acceptance(
    path: Path, analysis_release_id: str, selection_sha256: str, bundle: Path
) -> tuple[Path, Path]:
    if not path.is_file():
        raise PassportValidationError("MANUAL_ACCEPTANCE_MISSING", str(path))
    signature_path = path.with_name("passport_manual_acceptance.signature.json")
    if not signature_path.is_file():
        raise PassportValidationError(
            "MANUAL_ACCEPTANCE_SIGNATURE_MISSING", str(signature_path)
        )
    rows = read_tsv(path)
    if table_columns(path) != contract.MANUAL_ACCEPTANCE_COLUMNS or len(rows) != 1:
        raise PassportValidationError(
            "MANUAL_ACCEPTANCE_SCHEMA",
            "exactly one row with the frozen PASS06 schema is required",
        )
    row = rows[0]
    if (
        row["analysis_release_id"] != analysis_release_id
        or row["selection_sha256"] != selection_sha256
    ):
        raise PassportValidationError(
            "MANUAL_ACCEPTANCE_RELEASE", "release/selection mismatch"
        )
    booleans = (
        "call_states_reviewed",
        "provenance_states_reviewed",
        "visible_boundary_pass",
        "visible_testability_pass",
        "visible_source_dependence_pass",
        "visible_falsifier_pass",
    )
    if any(row[column].lower() != "true" for column in booleans):
        raise PassportValidationError(
            "MANUAL_ACCEPTANCE_FAILED", "all manual gates must pass"
        )
    if set(
        part.strip() for part in row["hero_genes_reviewed"].split(";") if part.strip()
    ) != {"THRB", "HKDC1", "GLP1R", "MTARC1"}:
        raise PassportValidationError(
            "MANUAL_HERO_GENE_SET", "manual review must cover THRB;HKDC1;GLP1R;MTARC1"
        )
    if row["decision"] != "accepted":
        raise PassportValidationError("MANUAL_ACCEPTANCE_FAILED", row["decision"])
    genes = pd.read_parquet(bundle / "passport_gene_index.parquet", engine="pyarrow")
    evidence = pd.read_parquet(
        bundle / "passport_evidence_long.parquet", engine="pyarrow"
    )
    ui_contract_path = bundle / "portal_candidate/ui_contract.json"
    try:
        ui_contract = json.loads(ui_contract_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PassportValidationError(
            "MANUAL_HERO_STATE_COVERAGE", "invalid UI contract"
        ) from exc
    if ui_contract.get("fixture_mode") is True:
        ordered = genes.assign(
            _symbol_key=genes["symbol"].astype(str).str.casefold(),
            _ensembl_key=genes["ensembl_id"].astype(str),
        ).sort_values(["_symbol_key", "_ensembl_key"], kind="stable")
        hero_ids = set(ordered["passport_id"].astype(str).head(4))
    else:
        hero_rows = genes.loc[
            genes["symbol"].astype(str).isin(passport_ui.HERO_SYMBOLS),
            ["symbol", "passport_id"],
        ]
        counts = Counter(hero_rows["symbol"].astype(str))
        if set(counts) != set(passport_ui.HERO_SYMBOLS) or any(
            count != 1 for count in counts.values()
        ):
            raise PassportValidationError(
                "MANUAL_HERO_STATE_COVERAGE", "hero identity is not unique"
            )
        hero_ids = set(hero_rows["passport_id"].astype(str))
    all_calls = {str(value) for value in evidence["call_state"] if str(value)}
    all_provenance = {
        str(value) for value in evidence["provenance_state"] if str(value)
    }
    hero_evidence = evidence.loc[evidence["passport_id"].astype(str).isin(hero_ids)]
    hero_calls = {str(value) for value in hero_evidence["call_state"] if str(value)}
    hero_provenance = {
        str(value) for value in hero_evidence["provenance_state"] if str(value)
    }
    if not all_calls.issubset(hero_calls) or not all_provenance.issubset(
        hero_provenance
    ):
        raise PassportValidationError(
            "MANUAL_HERO_STATE_COVERAGE",
            f"missing_calls={sorted(all_calls - hero_calls)};"
            f"missing_provenance={sorted(all_provenance - hero_provenance)}",
        )
    ui_bindings = {
        "ui_index_sha256": bundle / "portal_candidate/index.html",
        "ui_contract_sha256": bundle / "portal_candidate/ui_contract.json",
        "ui_source_manifest_sha256": (
            bundle / "portal_candidate/ui_source_manifest.tsv"
        ),
        "ui_review_manifest_sha256": (
            bundle / "portal_candidate/review/review_manifest.tsv"
        ),
    }
    for field, artifact in ui_bindings.items():
        if not artifact.is_file() or artifact.is_symlink():
            raise PassportValidationError(
                "MANUAL_ACCEPTANCE_UI_MISSING", artifact.as_posix()
            )
        if (
            not SHA256_PATTERN.fullmatch(row[field])
            or row[field] != contract.sha256_file(artifact)
        ):
            raise PassportValidationError(
                "MANUAL_ACCEPTANCE_UI_HASH", f"{field}:{artifact}"
            )
    signature = json.loads(signature_path.read_text(encoding="utf-8"))
    required = {
        "attestation_version",
        "acceptance_sha256",
        *ui_bindings,
        "signed_by",
        "signed_at_utc",
        "decision_register_id",
        "fixture_only",
    }
    if set(signature) != required:
        raise PassportValidationError(
            "MANUAL_ACCEPTANCE_SIGNATURE_SCHEMA",
            f"missing={sorted(required - set(signature))};"
            f"extra={sorted(set(signature) - required)}",
        )
    if (
        signature["attestation_version"]
        != contract.MANUAL_ACCEPTANCE_ATTESTATION_VERSION
    ):
        raise PassportValidationError(
            "MANUAL_ACCEPTANCE_SIGNATURE_VERSION", str(signature["attestation_version"])
        )
    if signature["acceptance_sha256"] != contract.sha256_file(path):
        raise PassportValidationError("MANUAL_ACCEPTANCE_SIGNATURE_HASH", str(path))
    if type(signature["fixture_only"]) is not bool or signature["fixture_only"]:
        raise PassportValidationError(
            "MANUAL_ACCEPTANCE_FIXTURE", "real manual review is required"
        )
    if (
        not str(signature["signed_by"]).strip()
        or not str(signature["decision_register_id"]).strip()
    ):
        raise PassportValidationError(
            "MANUAL_ACCEPTANCE_IDENTITY", "signer and decision ID are required"
        )
    if (
        signature["signed_by"] != row["reviewer"]
        or signature["signed_at_utc"] != row["reviewed_at_utc"]
    ):
        raise PassportValidationError(
            "MANUAL_ACCEPTANCE_IDENTITY", "reviewer/signature identity mismatch"
        )
    for field in ui_bindings:
        if signature[field] != row[field]:
            raise PassportValidationError(
                "MANUAL_ACCEPTANCE_SIGNATURE_UI_HASH", field
            )
    return path, signature_path


def finalize_pass06(
    bundle: Path,
    report: list[dict[str, str]],
    assembly,
    signature: Mapping[str, object],
    fixture_mode: bool,
    manual_acceptance: Path | None,
) -> Path:
    import build_evidence_passport_bundle as builder

    selection_sha = contract.sha256_file(bundle / "passport_input_selection.tsv")
    if fixture_mode:
        manual_status = "not_applicable_synthetic_fixture"
    else:
        if manual_acceptance is None:
            raise PassportValidationError(
                "MANUAL_ACCEPTANCE_MISSING",
                "--manual-acceptance is required for PASS06 terminal sealing",
            )
        manual_path, manual_signature = validate_manual_acceptance(
            manual_acceptance.resolve(),
            str(signature["analysis_release_id"]),
            selection_sha,
            bundle,
        )
        manual_destination = bundle / "passport_manual_acceptance.tsv"
        signature_destination = bundle / "passport_manual_acceptance.signature.json"
        for source, destination in (
            (manual_path, manual_destination),
            (manual_signature, signature_destination),
        ):
            if source != destination and (destination.exists() or destination.is_symlink()):
                raise PassportValidationError(
                    "MANUAL_ACCEPTANCE_DESTINATION_EXISTS", str(destination)
                )
        if manual_path != manual_destination:
            shutil.copy2(manual_path, manual_destination)
        if manual_signature != signature_destination:
            shutil.copy2(manual_signature, signature_destination)
        manual_status = "passed_signed_manual_acceptance"
        report.append(
            {
                "check_id": "PASS06_MANUAL_ACCEPTANCE",
                "status": "pass",
                "detail": "signed manual UI review passed",
            }
        )
    gate_rows = read_tsv(bundle / "passport_gate_status.tsv")
    for row in gate_rows:
        if row["gate_id"] == "PASS06_TERMINAL":
            row.update(
                status="fixture_validated" if fixture_mode else "pass",
                promotion_allowed="false",
                reason="Synthetic fixture only."
                if fixture_mode
                else "Automated and signed manual acceptance passed; Plan60 handoff allowed.",
            )
    contract.atomic_write_tsv(
        bundle / "passport_gate_status.tsv",
        gate_rows,
        contract.TSV_SCHEMAS["passport_gate_status.tsv"],
    )
    build_rows = read_tsv(bundle / "passport_build_status.tsv")
    for row in build_rows:
        if row["stage_id"] == "PASS06":
            row.update(
                status="fixture_validated" if fixture_mode else "pass",
                promotion_allowed="false",
                reason="terminal semantic validation complete",
            )
    contract.atomic_write_tsv(
        bundle / "passport_build_status.tsv",
        build_rows,
        contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_build_status.tsv"],
    )
    manifest = builder._write_manifest(  # noqa: SLF001
        bundle,
        selection_sha,
        "synthetic_fixture_not_for_promotion"
        if fixture_mode
        else "validated_candidate_handoff",
        finalized=True,
    )
    manifest_sha = contract.sha256_file(manifest)
    analysis_release_id = str(signature["analysis_release_id"])
    bundle_uri = release_bundle_uri(analysis_release_id)
    handoff = bundle / "passport_plan60_handoff.tsv"
    contract.atomic_write_tsv(
        handoff,
        [
            {
                "analysis_release_id": analysis_release_id,
                "bundle_path": bundle_uri,
                "manifest_sha256": manifest_sha,
                "selection_sha256": selection_sha,
                "automated_validation": "pass",
                "manual_acceptance": manual_status,
                "promotion_allowed": "false",
                "build_command": (
                    "python scripts/portal/build_evidence_passport_bundle.py "
                    f"--selection {contract.CANDIDATE_REL.as_posix()}/passport_input_selection.tsv"
                ),
                "environment": f"spatial:pandas-{pd.__version__}",
                "scientific_call_recomputed": "false",
                "omitted_input_inventory": (
                    f"dataset_skips={sum(row['status'] == 'skipped_source_gate' for row in read_tsv(bundle / 'passport_dataset_status.tsv'))};"
                    f"identity_quarantine={len(read_tsv(bundle / 'passport_identity_quarantine.tsv'))};"
                    "gene_domains_not_populated=spatial;proteomics;chromatin;functional;cross_species;cell_context"
                ),
            }
        ],
        contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_plan60_handoff.tsv"],
    )
    terminal_provenance = bundle / TERMINAL_PROVENANCE_FILENAME
    terminal_rows = []
    for artifact_role, artifact, producer_id, producer_path in (
        (
            "payload_manifest", manifest,
            builder.BUNDLE_LOGICAL_PRODUCER_ID, builder.SCRIPT_PATH,
        ),
        ("plan60_handoff", handoff, LOGICAL_PRODUCER_ID, SCRIPT_PATH),
    ):
        terminal_rows.append(
            {
                "analysis_release_id": analysis_release_id,
                "bundle_uri": bundle_uri,
                "artifact_role": artifact_role,
                "relative_path": artifact.relative_to(bundle).as_posix(),
                "sha256": contract.sha256_file(artifact),
                "bytes": artifact.stat().st_size,
                "producer": producer_id,
                "producer_sha256": contract.sha256_file(producer_path),
            }
        )
    contract.atomic_write_tsv(
        terminal_provenance,
        terminal_rows,
        contract.PASS06_TERMINAL_PROVENANCE_COLUMNS,
    )
    report.append(
        {
            "check_id": TERMINAL_PROVENANCE_CHECK_ID,
            "status": "pass",
            "detail": json.dumps(
                {
                    "bundle_uri": bundle_uri,
                    "contract_version": "passport_terminal_provenance_v1",
                    "terminal_provenance_path": TERMINAL_PROVENANCE_FILENAME,
                    "terminal_provenance_sha256": contract.sha256_file(
                        terminal_provenance
                    ),
                },
                sort_keys=True,
                separators=(",", ":"),
            ),
        }
    )
    contract.atomic_write_tsv(
        bundle / "passport_validation_report.tsv",
        report,
        contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_validation_report.tsv"],
    )
    validated_at = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    seal_name = "FIXTURE_PASS06_VALIDATED" if fixture_mode else "PASS06_VALIDATED"
    seal_path = bundle / seal_name
    contract.atomic_write_tsv(
        seal_path,
        [
            {
                "analysis_release_id": analysis_release_id,
                "status": "validated_synthetic_pass02_06_fixture"
                if fixture_mode
                else "passport_release_validated",
                "selection_sha256": selection_sha,
                "manifest_sha256": manifest_sha,
                "validation_report_sha256": contract.sha256_file(
                    bundle / "passport_validation_report.tsv"
                ),
                "n_genes": len(assembly.genes),
                "n_evidence_rows": len(assembly.evidence),
                "n_programs": len(assembly.programs),
                "n_program_context_rows": len(assembly.program_context),
                "automated_validation": "pass",
                "manual_acceptance": manual_status,
                "handoff_allowed": "false" if fixture_mode else "true",
                "canonical_promotion_authorized": "false",
                "scientific_call_recomputed": "false",
                "validated_at_utc": validated_at,
            }
        ],
        contract.PASS06_TERMINAL_COLUMNS,
    )
    verify_sealed_bundle(bundle)
    return seal_path


def run_self_test(candidate_root: Path) -> Path:
    candidate = contract.assert_candidate_root(candidate_root)
    fixture_root = contract.build_fixture_suite(candidate)
    report: list[dict[str, str]] = []
    differences = contract.compare_fixture_determinism(
        fixture_root / "positive", fixture_root / "positive_repeat"
    )
    if differences:
        raise PassportValidationError(
            "NONDETERMINISTIC_FIXTURE", ", ".join(differences)
        )
    report.append(
        {
            "check": "deterministic_rebuild",
            "status": "pass",
            "detail": "byte-identical fixture bundles",
        }
    )
    report.extend(validate_bundle(fixture_root / "positive", fixture_mode=True))
    zero_negative_results = validate_bundle(
        fixture_root / "positive_zero_tested_negative", fixture_mode=True
    )
    report.append(
        {
            "check": "zero_tested_negative_allowed",
            "status": "pass",
            "detail": f"validated {len(zero_negative_results)} bundle checks with zero tested-negative rows",
        }
    )

    expectations = read_tsv(fixture_root / "negative_fixture_expectations.tsv")
    for expectation in expectations:
        fixture = fixture_root / expectation["fixture"]
        expected_code = expectation["expected_error_code"]
        try:
            validate_bundle(fixture, fixture_mode=True)
        except PassportValidationError as exc:
            if exc.code != expected_code:
                raise PassportValidationError(
                    "NEGATIVE_FIXTURE_WRONG_FAILURE",
                    f"{fixture.name}: expected {expected_code}, observed {exc.code}: {exc.message}",
                ) from exc
            report.append(
                {
                    "check": f"negative_fixture:{fixture.name}",
                    "status": "pass",
                    "detail": f"rejected with {exc.code}",
                }
            )
        else:
            raise PassportValidationError(
                "NEGATIVE_FIXTURE_FALSE_PASS", f"{fixture.name} unexpectedly passed"
            )

    try:
        contract.validate_signed_selection(
            fixture_root / "deliberately_missing" / "passport_input_selection.tsv",
            allow_fixture=False,
        )
    except contract.PassportContractError as exc:
        if exc.code != "SIGNED_SELECTION_MISSING":
            raise
        report.append(
            {
                "check": "production_missing_selection",
                "status": "pass",
                "detail": "rejected with SIGNED_SELECTION_MISSING",
            }
        )
    else:
        raise PassportValidationError(
            "UNSIGNED_PRODUCTION_FALSE_PASS", "missing selection passed"
        )

    try:
        contract.validate_signed_selection(
            fixture_root / "positive" / "passport_input_selection.tsv",
            allow_fixture=False,
        )
    except contract.PassportContractError as exc:
        if exc.code != "FIXTURE_SELECTION_PRODUCTION":
            raise
        report.append(
            {
                "check": "fixture_attestation_production",
                "status": "pass",
                "detail": "rejected with FIXTURE_SELECTION_PRODUCTION",
            }
        )
    else:
        raise PassportValidationError(
            "FIXTURE_PRODUCTION_FALSE_PASS", "fixture unlocked production"
        )

    report_path = fixture_root / "fixture_self_test_report.tsv"
    contract.atomic_write_tsv(report_path, report, ["check", "status", "detail"])
    return report_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--fixture-mode", action="store_true")
    parser.add_argument("--pass02-06-fixture", action="store_true")
    parser.add_argument("--finalize-pass06", action="store_true")
    parser.add_argument("--verify-sealed", action="store_true")
    parser.add_argument("--manual-acceptance", type=Path)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument(
        "--candidate-root", type=Path, default=contract.DEFAULT_CANDIDATE_ROOT
    )
    parser.add_argument("--selection", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.verify_sealed:
        if args.bundle is None:
            raise PassportValidationError(
                "TERMINAL_BUNDLE_MISSING", "--bundle is required with --verify-sealed"
            )
        results = verify_sealed_bundle(args.bundle)
        print(json.dumps(results, sort_keys=True, separators=(",", ":")))
        return
    candidate = contract.assert_candidate_root(args.candidate_root)
    if args.self_test:
        report = run_self_test(candidate)
        print(
            f"PASS: PASS-01 positive and negative fixtures validated; report={report}"
        )
        return
    if args.fixture_mode:
        bundle = args.bundle or candidate / "fixtures" / "positive"
        results = validate_bundle(bundle, fixture_mode=True)
        print(json.dumps(results, sort_keys=True, separators=(",", ":")))
        return

    if args.pass02_06_fixture:
        if args.bundle is None:
            raise PassportValidationError(
                "FIXTURE_BUNDLE_MISSING",
                "--bundle is required for PASS02-06 fixture validation",
            )
        report, assembly, _, signature = validate_production_bundle(
            args.bundle, fixture_mode=True
        )
        if args.finalize_pass06:
            seal = finalize_pass06(
                args.bundle.resolve(),
                report,
                assembly,
                signature,
                fixture_mode=True,
                manual_acceptance=None,
            )
            print(f"PASS: PASS-02--PASS-06 fixture validated; seal={seal}")
        else:
            print(json.dumps(report, sort_keys=True, separators=(",", ":")))
        return

    selection = args.selection or candidate / "passport_input_selection.tsv"
    contract.validate_signed_selection(selection, allow_fixture=False)
    if args.bundle is None:
        raise PassportValidationError(
            "PRODUCTION_BUNDLE_MISSING", "--bundle is required for PASS-06 validation"
        )
    report, assembly, _, signature = validate_production_bundle(
        args.bundle, fixture_mode=False
    )
    if args.finalize_pass06:
        seal = finalize_pass06(
            args.bundle.resolve(),
            report,
            assembly,
            signature,
            fixture_mode=False,
            manual_acceptance=args.manual_acceptance,
        )
        print(f"PASS: PASS-06 validated candidate handoff; seal={seal}")
    else:
        print(json.dumps(report, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    try:
        main()
    except contract.PassportContractError as exc:
        raise SystemExit(str(exc)) from exc
