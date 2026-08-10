#!/usr/bin/env python3
"""Shared, standard-library-only contract checks for spatial-context v2.

This module deliberately does not read assay outcomes.  It validates the
interface between already-sealed program registries and dataset-native adapter
bundles.  Dataset-native producers remain responsible for their own null
models and effect estimates.
"""

from __future__ import annotations

import csv
import hashlib
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence


RELEASE_ID = "program-context-v2-candidate-2026-08-07"
SEMANTIC_CONTRACT_ID = "evidence-state-v2-2026-08-08"

ALLOWED_EVIDENCE_STATES = {
    "robust",
    "indeterminate",
    "tested_negative",
    "untestable",
    "not_applicable",
    "skipped",
    "source_dependent",
}
ALLOWED_INTERVAL_TYPES = {
    "confidence_interval",
    "donor_effect_range",
    "none",
    "not_applicable",
}
ALLOWED_SOURCE_DEPENDENCE = {"independent", "source_dependent", "descriptive"}
ALLOWED_GATE_STATES = {
    "pass",
    "indeterminate",
    "valid_null",
    "skipped",
    "untestable",
    "not_applicable",
}
ALLOWED_RESOLUTION = {"resolved", "unresolved", "not_applicable"}
PROHIBITED_COLUMNS = {
    "universal_score",
    "combined_score",
    "cross_assay_score",
    "cross_assay_rank",
    "overall_score",
    "overall_rank",
    "modality_count",
    "evidence_count",
}
PROHIBITED_EFFECT_UNITS = {
    "combined_score",
    "universal_score",
    "standardized_cross_assay_score",
    "modality_count",
    "evidence_count",
    "overall_rank",
}


SCHEMAS: dict[str, tuple[str, ...]] = {
    "sample_manifest.tsv": (
        "release_id",
        "dataset",
        "analysis_set_id",
        "biological_id",
        "technical_id",
        "biological_unit",
        "technical_unit",
        "include_primary",
        "gate_state",
    ),
    "gene_mapping_audit.tsv": (
        "release_id",
        "dataset",
        "analysis_set_id",
        "program_uid",
        "membership_sha256",
        "n_source_genes",
        "n_genes_measured",
        "retained_l1_weight",
        "mapping_status",
    ),
    "design_audit.tsv": (
        "release_id",
        "dataset",
        "analysis_set_id",
        "contrast_or_exposure",
        "biological_unit",
        "biological_unit_resolution",
        "n_biological",
        "technical_unit",
        "n_technical",
        "design_status",
        "biological_ids_sha256",
        "technical_ids_sha256",
    ),
    "program_testability.tsv": (
        "release_id",
        "registry_sha256",
        "dataset",
        "analysis_set_id",
        "program_uid",
        "membership_sha256",
        "testable",
        "n_genes_measured",
        "retained_l1_weight",
        "testability_reason",
        "evidence_state",
    ),
    "per_sample_program_scores.tsv": (
        "release_id",
        "registry_sha256",
        "dataset",
        "analysis_set_id",
        "program_uid",
        "membership_sha256",
        "biological_id",
        "technical_unit_count",
        "program_score",
        "score_unit",
    ),
    "program_effects.tsv": (
        "release_id",
        "registry_sha256",
        "membership_sha256",
        "program_uid",
        "legacy_program_id",
        "program_label",
        "cell_type",
        "dataset",
        "assay",
        "source_publication",
        "source_dependence",
        "analysis_set_id",
        "biological_unit",
        "biological_unit_resolution",
        "n_biological",
        "technical_unit",
        "n_technical",
        "contrast_or_exposure",
        "effect_unit",
        "estimate",
        "std_error",
        "matched_null_sd",
        "interval_low",
        "interval_high",
        "interval_type",
        "pvalue",
        "padj",
        "pvalue_method",
        "multiplicity_family",
        "n_genes_measured",
        "retained_l1_weight",
        "testable",
        "testability_reason",
        "direction_expected",
        "direction_observed",
        "descriptive_effect_direction",
        "inferential_test_direction",
        "direction_agreement",
        "heterogeneity_statistic",
        "heterogeneity_df",
        "heterogeneity_pvalue",
        "sensitivity_sign_agree",
        "robustness_pass",
        "negative_call_rule_id",
        "cross_assay_comparable",
        "evidence_state",
        "producer",
        "producer_sha256",
        "source_manifest_sha256",
    ),
    "sensitivity.tsv": (
        "release_id",
        "registry_sha256",
        "dataset",
        "analysis_set_id",
        "program_uid",
        "membership_sha256",
        "sensitivity_id",
        "estimate",
        "effect_unit",
        "direction_observed",
        "sign_agree",
        "status",
    ),
    "gate_status.tsv": (
        "release_id",
        "dataset",
        "analysis_set_id",
        "gate_status",
        "gate_reason",
        "outcomes_tested",
        "registry_sha256",
        "source_manifest_sha256",
        "execution_manifest_sha256",
    ),
    "source_manifest.tsv": (
        "path_scope",
        "relative_path",
        "bytes",
        "sha256",
        "source_role",
        "public_access",
        "retrieved_utc",
    ),
    "execution_manifest.tsv": (
        "role",
        "relative_path",
        "bytes",
        "sha256",
        "created_utc",
    ),
}


ADAPTER_REGISTRY_COLUMNS = (
    "adapter_id",
    "adapter_root",
    "dataset",
    "assay",
    "source_publication",
    "source_dependence",
    "biological_unit",
    "biological_unit_resolution",
    "technical_unit",
    "analysis_set_id",
    "contrast_or_exposure",
    "effect_unit",
    "cross_assay_comparable",
    "program_universe",
    "min_genes_testable",
    "min_retained_l1_weight",
    "gate_expectation",
    "registry_sha256",
    "ready_sha256",
    "status",
)


class ContractError(RuntimeError):
    """Raised when a candidate adapter violates the frozen contract."""


@dataclass(frozen=True)
class Seal:
    release_id: str
    registry_sha256: str
    membership_table_sha256: str
    ready_sha256: str
    registry_rows: tuple[dict[str, str], ...]
    registry_by_uid: Mapping[str, dict[str, str]]


@dataclass(frozen=True)
class Check:
    check_id: str
    status: str
    detail: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def canonical_id_hash(values: Iterable[str]) -> str:
    values_sorted = sorted(set(values))
    return sha256_text("\n".join(values_sorted))


def read_tsv(path: Path, required: Sequence[str] | None = None) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ContractError(f"missing required TSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ContractError(f"missing TSV header: {path}")
        header = list(reader.fieldnames)
        if len(header) != len(set(header)):
            raise ContractError(f"duplicate columns in {path}")
        prohibited = sorted(set(header) & PROHIBITED_COLUMNS)
        if prohibited:
            raise ContractError(f"prohibited universal/composite columns in {path}: {prohibited}")
        if required is not None:
            missing = [column for column in required if column not in header]
            if missing:
                raise ContractError(f"missing columns in {path}: {missing}")
        rows = [{key: (value if value is not None else "") for key, value in row.items()} for row in reader]
    return header, rows


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: _format_value(row.get(column, "")) for column in columns})


def _format_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return str(value)


def parse_bool(value: str, label: str) -> bool:
    normalized = value.strip().upper()
    if normalized == "TRUE":
        return True
    if normalized == "FALSE":
        return False
    raise ContractError(f"{label} must be TRUE or FALSE, got {value!r}")


def parse_int(value: str, label: str, *, allow_empty: bool = False) -> int | None:
    if value.strip() == "" and allow_empty:
        return None
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ContractError(f"{label} must be an integer, got {value!r}") from exc
    return parsed


def parse_float(value: str, label: str, *, allow_empty: bool = False) -> float | None:
    if value.strip() == "" and allow_empty:
        return None
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ContractError(f"{label} must be numeric, got {value!r}") from exc
    if not math.isfinite(parsed):
        raise ContractError(f"{label} must be finite, got {value!r}")
    return parsed


def _one_row(path: Path, required: Sequence[str]) -> dict[str, str]:
    _, rows = read_tsv(path, required)
    if len(rows) != 1:
        raise ContractError(f"expected exactly one row in {path}, found {len(rows)}")
    return rows[0]


def validate_hotspot_seal(hotspot_root: Path) -> Seal:
    ready_path = hotspot_root / "READY"
    validation_path = hotspot_root / "validation_status.tsv"
    manifest_path = hotspot_root / "release_manifest.tsv"
    registry_path = hotspot_root / "program_registry_v2.tsv"
    membership_path = hotspot_root / "program_membership_v2.tsv"

    ready = _one_row(
        ready_path,
        (
            "release_id",
            "status",
            "validation_status_sha256",
            "release_manifest_sha256",
            "registry_sha256",
            "membership_table_sha256",
            "external_outcomes_read",
        ),
    )
    if ready["release_id"] != RELEASE_ID:
        raise ContractError("Hotspot READY release_id does not match the candidate")
    if ready["status"] != "ready_for_external_testing":
        raise ContractError(f"Hotspot READY status is not ready_for_external_testing: {ready['status']}")
    if parse_bool(ready["external_outcomes_read"], "READY.external_outcomes_read"):
        raise ContractError("Hotspot READY says external outcomes were read")

    expected_hashes = {
        validation_path: ready["validation_status_sha256"],
        manifest_path: ready["release_manifest_sha256"],
        registry_path: ready["registry_sha256"],
        membership_path: ready["membership_table_sha256"],
    }
    for path, expected in expected_hashes.items():
        observed = sha256_file(path)
        if observed != expected:
            raise ContractError(f"Hotspot READY hash mismatch for {path}: {observed} != {expected}")

    validation = _one_row(
        validation_path,
        (
            "release_id",
            "status",
            "registry_sha256",
            "membership_table_sha256",
            "release_manifest_sha256",
            "external_outcomes_read",
        ),
    )
    if validation["release_id"] != RELEASE_ID or validation["status"] != "passed_independent_validation":
        raise ContractError("Hotspot validation_status is not the sealed, independently validated release")
    if parse_bool(validation["external_outcomes_read"], "validation_status.external_outcomes_read"):
        raise ContractError("Hotspot validation_status says external outcomes were read")
    for column in ("registry_sha256", "membership_table_sha256", "release_manifest_sha256"):
        if validation[column] != ready[column]:
            raise ContractError(f"Hotspot READY/validation_status mismatch for {column}")

    _, registry_rows = read_tsv(
        registry_path,
        (
            "release_id",
            "program_uid",
            "membership_sha256",
            "cell_type",
            "module",
            "module_name",
            "robust_display",
            "external_test_eligible",
        ),
    )
    if len(registry_rows) != 117:
        raise ContractError(f"Hotspot registry must contain 117 rows, found {len(registry_rows)}")
    registry_by_uid: dict[str, dict[str, str]] = {}
    for row in registry_rows:
        if row["release_id"] != RELEASE_ID:
            raise ContractError("registry release_id drift")
        uid = row["program_uid"]
        if not uid or uid in registry_by_uid:
            raise ContractError(f"duplicate or empty program_uid in registry: {uid!r}")
        parse_bool(row["robust_display"], f"registry.robust_display[{uid}]")
        parse_bool(row["external_test_eligible"], f"registry.external_test_eligible[{uid}]")
        registry_by_uid[uid] = row

    _, membership_rows = read_tsv(membership_path, ("program_uid", "membership_sha256"))
    membership_pairs = {(row["program_uid"], row["membership_sha256"]) for row in membership_rows}
    for uid, row in registry_by_uid.items():
        if (uid, row["membership_sha256"]) not in membership_pairs:
            raise ContractError(f"registry membership hash absent from membership table for {uid}")

    return Seal(
        release_id=RELEASE_ID,
        registry_sha256=ready["registry_sha256"],
        membership_table_sha256=ready["membership_table_sha256"],
        ready_sha256=sha256_file(ready_path),
        registry_rows=tuple(registry_rows),
        registry_by_uid=registry_by_uid,
    )


def program_universe(seal: Seal, selector: str) -> dict[str, dict[str, str]]:
    if selector == "all_registry":
        return dict(seal.registry_by_uid)
    if selector not in {"robust_display", "external_test_eligible"}:
        raise ContractError(f"unsupported program_universe selector: {selector}")
    return {
        uid: row
        for uid, row in seal.registry_by_uid.items()
        if parse_bool(row[selector], f"registry.{selector}[{uid}]")
    }


def _assert_same(value: str, expected: str, label: str) -> None:
    if value != expected:
        raise ContractError(f"{label}: {value!r} != {expected!r}")


def _rows_by_uid(rows: Sequence[dict[str, str]], label: str) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        uid = row.get("program_uid", "")
        if not uid or uid in result:
            raise ContractError(f"{label} has duplicate or empty program_uid: {uid!r}")
        result[uid] = row
    return result


def _resolve_manifest_path(
    adapter_root: Path,
    project_root: Path,
    filename: str,
    row: Mapping[str, str],
) -> Path:
    relative = Path(row["relative_path"])
    if relative.is_absolute():
        raise ContractError(f"{filename} paths must be relative: {relative}")
    if filename == "source_manifest.tsv":
        scope = row["path_scope"]
        if scope == "adapter_relative":
            base = adapter_root
        elif scope == "project_relative":
            base = project_root
        else:
            raise ContractError(f"invalid source-manifest path_scope: {scope!r}")
    else:
        base = adapter_root
    base_resolved = base.resolve()
    path = (base_resolved / relative).resolve()
    if not path.is_relative_to(base_resolved):
        raise ContractError(f"{filename} path escapes its declared scope: {relative}")
    return path


def _validate_manifest(
    adapter_root: Path,
    project_root: Path,
    filename: str,
    rows: Sequence[dict[str, str]],
) -> None:
    for row in rows:
        relative = row["relative_path"]
        path = _resolve_manifest_path(adapter_root, project_root, filename, row)
        if not path.is_file():
            raise ContractError(f"{filename} points to missing file: {relative}")
        observed_bytes = path.stat().st_size
        observed_hash = sha256_file(path)
        if parse_int(row["bytes"], f"{filename}.bytes[{relative}]") != observed_bytes:
            raise ContractError(f"{filename} byte mismatch for {relative}")
        if row["sha256"] != observed_hash:
            raise ContractError(f"{filename} hash mismatch for {relative}")


def validate_adapter(
    candidate_root: Path,
    adapter: Mapping[str, str],
    seal: Seal,
    project_root: Path,
) -> list[Check]:
    checks: list[Check] = []
    adapter_id = adapter["adapter_id"]
    adapter_root = candidate_root / adapter["adapter_root"]
    if not adapter_root.is_dir():
        raise ContractError(f"adapter_root does not exist for {adapter_id}: {adapter_root}")

    _assert_same(adapter["registry_sha256"], seal.registry_sha256, f"{adapter_id}.registry_sha256")
    _assert_same(adapter["ready_sha256"], seal.ready_sha256, f"{adapter_id}.ready_sha256")
    if adapter["source_dependence"] not in ALLOWED_SOURCE_DEPENDENCE:
        raise ContractError(f"invalid source_dependence for {adapter_id}")
    if adapter["biological_unit_resolution"] not in ALLOWED_RESOLUTION:
        raise ContractError(f"invalid biological_unit_resolution for {adapter_id}")
    if parse_bool(adapter["cross_assay_comparable"], f"{adapter_id}.cross_assay_comparable"):
        raise ContractError(f"{adapter_id} claims cross-assay effect comparability")
    if adapter["effect_unit"] in PROHIBITED_EFFECT_UNITS or not adapter["effect_unit"].strip():
        raise ContractError(f"{adapter_id} uses prohibited or empty effect unit")

    tables: dict[str, list[dict[str, str]]] = {}
    for filename, schema in SCHEMAS.items():
        _, rows = read_tsv(adapter_root / filename, schema)
        tables[filename] = rows
    checks.append(Check("required_artifacts", "pass", f"{adapter_id}: {len(SCHEMAS)} required artifacts"))

    source_manifest_hash = sha256_file(adapter_root / "source_manifest.tsv")
    execution_manifest_hash = sha256_file(adapter_root / "execution_manifest.tsv")
    _validate_manifest(adapter_root, project_root, "source_manifest.tsv", tables["source_manifest.tsv"])
    _validate_manifest(adapter_root, project_root, "execution_manifest.tsv", tables["execution_manifest.tsv"])
    checks.append(Check("manifest_integrity", "pass", f"{adapter_id}: source and execution hashes rederived"))

    gates = tables["gate_status.tsv"]
    if len(gates) != 1:
        raise ContractError(f"{adapter_id} must have exactly one gate row")
    gate = gates[0]
    for column, expected in (
        ("release_id", seal.release_id),
        ("dataset", adapter["dataset"]),
        ("analysis_set_id", adapter["analysis_set_id"]),
        ("gate_status", adapter["gate_expectation"]),
        ("registry_sha256", seal.registry_sha256),
        ("source_manifest_sha256", source_manifest_hash),
        ("execution_manifest_sha256", execution_manifest_hash),
    ):
        _assert_same(gate[column], expected, f"{adapter_id}.gate.{column}")
    if gate["gate_status"] not in ALLOWED_GATE_STATES:
        raise ContractError(f"invalid gate state for {adapter_id}: {gate['gate_status']}")
    outcomes_tested = parse_bool(gate["outcomes_tested"], f"{adapter_id}.outcomes_tested")
    if gate["gate_status"] in {"skipped", "untestable", "not_applicable"} and outcomes_tested:
        raise ContractError(f"{adapter_id} terminal gate cannot have outcomes_tested=TRUE")
    if gate["gate_status"] in {"pass", "indeterminate", "valid_null"} and not outcomes_tested:
        raise ContractError(f"{adapter_id} passing gate must have outcomes_tested=TRUE")

    universe = program_universe(seal, adapter["program_universe"])
    expected_uids = set(universe)
    if not expected_uids:
        raise ContractError(f"{adapter_id} selected an empty registry universe")
    mapping = _rows_by_uid(tables["gene_mapping_audit.tsv"], f"{adapter_id}.gene_mapping")
    testability = _rows_by_uid(tables["program_testability.tsv"], f"{adapter_id}.testability")
    effects = _rows_by_uid(tables["program_effects.tsv"], f"{adapter_id}.effects")
    if set(mapping) != expected_uids or set(testability) != expected_uids or set(effects) != expected_uids:
        raise ContractError(f"{adapter_id} does not emit exactly one mapping/testability/effect row per prespecified program")
    checks.append(Check("complete_program_universe", "pass", f"{adapter_id}: {len(expected_uids)} programs"))

    samples = tables["sample_manifest.tsv"]
    included_samples = []
    technical_to_biological: dict[str, str] = {}
    for row in samples:
        for column, expected in (
            ("release_id", seal.release_id),
            ("dataset", adapter["dataset"]),
            ("analysis_set_id", adapter["analysis_set_id"]),
            ("biological_unit", adapter["biological_unit"]),
            ("technical_unit", adapter["technical_unit"]),
        ):
            _assert_same(row[column], expected, f"{adapter_id}.sample.{column}")
        technical_id = row["technical_id"]
        if not technical_id or technical_id in technical_to_biological:
            raise ContractError(f"{adapter_id} has duplicate or empty technical_id: {technical_id!r}")
        technical_to_biological[technical_id] = row["biological_id"]
        if parse_bool(row["include_primary"], f"{adapter_id}.include_primary[{technical_id}]"):
            if adapter["biological_unit_resolution"] == "resolved" and not row["biological_id"]:
                raise ContractError(f"{adapter_id} included technical unit lacks biological_id")
            included_samples.append(row)

    observed_biological_ids = sorted(
        {row["biological_id"] for row in included_samples if row["biological_id"]}
    )
    observed_technical_ids = sorted(row["technical_id"] for row in included_samples)
    design_rows = tables["design_audit.tsv"]
    if len(design_rows) != 1:
        raise ContractError(f"{adapter_id} must have exactly one design row per fixture adapter")
    design = design_rows[0]
    for column, expected in (
        ("release_id", seal.release_id),
        ("dataset", adapter["dataset"]),
        ("analysis_set_id", adapter["analysis_set_id"]),
        ("contrast_or_exposure", adapter["contrast_or_exposure"]),
        ("biological_unit", adapter["biological_unit"]),
        ("biological_unit_resolution", adapter["biological_unit_resolution"]),
        ("technical_unit", adapter["technical_unit"]),
    ):
        _assert_same(design[column], expected, f"{adapter_id}.design.{column}")
    n_biological = parse_int(
        design["n_biological"],
        f"{adapter_id}.design.n_biological",
        allow_empty=adapter["biological_unit_resolution"] == "unresolved",
    )
    n_technical = parse_int(design["n_technical"], f"{adapter_id}.design.n_technical")
    if adapter["biological_unit_resolution"] == "unresolved":
        if n_biological is not None or observed_biological_ids:
            raise ContractError(
                f"{adapter_id} unresolved arrays/units must not be reported as biological replicates"
            )
    elif n_biological != len(observed_biological_ids):
        raise ContractError(f"{adapter_id} resolved biological count does not match the sample manifest")
    if n_technical != len(observed_technical_ids):
        raise ContractError(f"{adapter_id} design counts do not match the sample manifest")
    if design["biological_ids_sha256"] != canonical_id_hash(observed_biological_ids):
        raise ContractError(f"{adapter_id} biological ID-set hash mismatch")
    if design["technical_ids_sha256"] != canonical_id_hash(observed_technical_ids):
        raise ContractError(f"{adapter_id} technical ID-set hash mismatch")
    if gate["gate_status"] in {"pass", "indeterminate", "valid_null"}:
        if adapter["source_dependence"] == "independent" and (
            adapter["biological_unit_resolution"] != "resolved"
            or n_biological is None
            or n_biological <= 0
        ):
            raise ContractError(f"{adapter_id} independent inference lacks resolved biological units")
        if n_technical <= 0:
            raise ContractError(f"{adapter_id} inference lacks technical units")
    checks.append(
        Check(
            "biological_units",
            "pass",
            f"{adapter_id}: {n_biological if n_biological is not None else 'unknown'} "
            f"{adapter['biological_unit']}; {n_technical} {adapter['technical_unit']}",
        )
    )

    min_genes = parse_int(adapter["min_genes_testable"], f"{adapter_id}.min_genes_testable")
    min_weight = parse_float(adapter["min_retained_l1_weight"], f"{adapter_id}.min_retained_l1_weight")
    assert min_genes is not None and min_weight is not None

    technical_counts: dict[str, int] = {}
    for row in included_samples:
        technical_counts[row["biological_id"]] = technical_counts.get(row["biological_id"], 0) + 1

    score_rows = tables["per_sample_program_scores.tsv"]
    score_keys: set[tuple[str, str]] = set()
    for row in score_rows:
        uid = row["program_uid"]
        biological_id = row["biological_id"]
        key = (uid, biological_id)
        if key in score_keys:
            raise ContractError(f"{adapter_id} duplicate donor/program score: {key}")
        score_keys.add(key)
        if uid not in expected_uids or biological_id not in technical_counts:
            raise ContractError(f"{adapter_id} score references an out-of-contract program or biological unit")
        if parse_int(row["technical_unit_count"], f"{adapter_id}.technical_unit_count[{key}]") != technical_counts[biological_id]:
            raise ContractError(f"{adapter_id} donor score technical-unit count mismatch for {key}")
        parse_float(row["program_score"], f"{adapter_id}.program_score[{key}]")
        for column, expected in (
            ("release_id", seal.release_id),
            ("registry_sha256", seal.registry_sha256),
            ("dataset", adapter["dataset"]),
            ("analysis_set_id", adapter["analysis_set_id"]),
            ("membership_sha256", universe[uid]["membership_sha256"]),
        ):
            _assert_same(row[column], expected, f"{adapter_id}.score.{column}")

    for uid in sorted(expected_uids):
        reg = universe[uid]
        map_row = mapping[uid]
        test_row = testability[uid]
        effect = effects[uid]
        for row_label, row in (("mapping", map_row), ("testability", test_row), ("effect", effect)):
            _assert_same(row["release_id"], seal.release_id, f"{adapter_id}.{row_label}.release_id[{uid}]")
            _assert_same(row["dataset"], adapter["dataset"], f"{adapter_id}.{row_label}.dataset[{uid}]")
            _assert_same(row["analysis_set_id"], adapter["analysis_set_id"], f"{adapter_id}.{row_label}.analysis_set_id[{uid}]")
            _assert_same(row["membership_sha256"], reg["membership_sha256"], f"{adapter_id}.{row_label}.membership[{uid}]")

        _assert_same(test_row["registry_sha256"], seal.registry_sha256, f"{adapter_id}.testability.registry[{uid}]")
        _assert_same(effect["registry_sha256"], seal.registry_sha256, f"{adapter_id}.effect.registry[{uid}]")
        _assert_same(effect["legacy_program_id"], f"{reg['cell_type']}:{reg['module']}", f"{adapter_id}.legacy_program_id[{uid}]")
        _assert_same(effect["program_label"], reg["module_name"], f"{adapter_id}.program_label[{uid}]")
        _assert_same(effect["cell_type"], reg["cell_type"], f"{adapter_id}.cell_type[{uid}]")
        for column, expected in (
            ("assay", adapter["assay"]),
            ("source_publication", adapter["source_publication"]),
            ("source_dependence", adapter["source_dependence"]),
            ("biological_unit", adapter["biological_unit"]),
            ("biological_unit_resolution", adapter["biological_unit_resolution"]),
            ("technical_unit", adapter["technical_unit"]),
            ("contrast_or_exposure", adapter["contrast_or_exposure"]),
            ("effect_unit", adapter["effect_unit"]),
            ("source_manifest_sha256", source_manifest_hash),
        ):
            _assert_same(effect[column], expected, f"{adapter_id}.effect.{column}[{uid}]")
        if parse_bool(effect["cross_assay_comparable"], f"{adapter_id}.effect.cross_assay_comparable[{uid}]"):
            raise ContractError(f"{adapter_id} effect row claims cross-assay comparability")
        if effect["evidence_state"] not in ALLOWED_EVIDENCE_STATES:
            raise ContractError(f"{adapter_id} invalid evidence state for {uid}: {effect['evidence_state']}")
        if test_row["evidence_state"] != effect["evidence_state"]:
            raise ContractError(f"{adapter_id} testability/effect state mismatch for {uid}")

        n_genes = parse_int(effect["n_genes_measured"], f"{adapter_id}.n_genes_measured[{uid}]")
        retained = parse_float(effect["retained_l1_weight"], f"{adapter_id}.retained_l1_weight[{uid}]")
        assert n_genes is not None and retained is not None
        if not 0.0 <= retained <= 1.0:
            raise ContractError(f"{adapter_id} retained L1 weight outside [0,1] for {uid}")
        for row_label, row in (("mapping", map_row), ("testability", test_row)):
            if parse_int(row["n_genes_measured"], f"{adapter_id}.{row_label}.n_genes[{uid}]") != n_genes:
                raise ContractError(f"{adapter_id} measured-gene mismatch for {uid}")
            if not math.isclose(
                parse_float(row["retained_l1_weight"], f"{adapter_id}.{row_label}.weight[{uid}]") or 0.0,
                retained,
                rel_tol=0.0,
                abs_tol=1e-12,
            ):
                raise ContractError(f"{adapter_id} retained-weight mismatch for {uid}")

        is_testable = parse_bool(effect["testable"], f"{adapter_id}.effect.testable[{uid}]")
        if parse_bool(test_row["testable"], f"{adapter_id}.testability.testable[{uid}]") != is_testable:
            raise ContractError(f"{adapter_id} testable flag mismatch for {uid}")
        if is_testable and (n_genes < min_genes or retained < min_weight):
            raise ContractError(f"{adapter_id} marks under-covered program testable: {uid}")
        expected_score_keys = {(uid, biological_id) for biological_id in observed_biological_ids} if is_testable else set()
        observed_score_keys = {key for key in score_keys if key[0] == uid}
        if gate["gate_status"] in {"pass", "indeterminate", "valid_null"} and observed_score_keys != expected_score_keys:
            raise ContractError(f"{adapter_id} donor-level scores are incomplete for {uid}")

        if parse_int(
            effect["n_biological"],
            f"{adapter_id}.effect.n_biological[{uid}]",
            allow_empty=n_biological is None,
        ) != n_biological:
            raise ContractError(f"{adapter_id} effect biological-unit count mismatch for {uid}")
        if parse_int(effect["n_technical"], f"{adapter_id}.effect.n_technical[{uid}]") != n_technical:
            raise ContractError(f"{adapter_id} effect technical-unit count mismatch for {uid}")

        producer_path = project_root / effect["producer"]
        if not producer_path.is_file() or sha256_file(producer_path) != effect["producer_sha256"]:
            raise ContractError(f"{adapter_id} producer provenance mismatch for {uid}")

        numeric_columns = (
            "estimate",
            "std_error",
            "matched_null_sd",
            "interval_low",
            "interval_high",
            "pvalue",
            "padj",
            "heterogeneity_statistic",
            "heterogeneity_df",
            "heterogeneity_pvalue",
        )
        interval_type = effect["interval_type"]
        if interval_type not in ALLOWED_INTERVAL_TYPES:
            raise ContractError(f"{adapter_id} invalid interval_type for {uid}: {interval_type!r}")
        interval_low_text = effect["interval_low"].strip()
        interval_high_text = effect["interval_high"].strip()
        if bool(interval_low_text) != bool(interval_high_text):
            raise ContractError(f"{adapter_id} interval bounds must be both present or both absent for {uid}")
        has_interval = bool(interval_low_text)
        if has_interval and interval_type not in {"confidence_interval", "donor_effect_range"}:
            raise ContractError(f"{adapter_id} paired interval has incompatible interval_type for {uid}")
        if not has_interval and interval_type in {"confidence_interval", "donor_effect_range"}:
            raise ContractError(f"{adapter_id} declared interval_type lacks paired bounds for {uid}")
        interval_low = interval_high = None
        if has_interval:
            interval_low = parse_float(effect["interval_low"], f"{adapter_id}.interval_low[{uid}]")
            interval_high = parse_float(effect["interval_high"], f"{adapter_id}.interval_high[{uid}]")
            assert interval_low is not None and interval_high is not None
            if interval_low > interval_high:
                raise ContractError(f"{adapter_id} interval_low exceeds interval_high for {uid}")

        if effect["evidence_state"] in {"untestable", "not_applicable", "skipped"}:
            if (
                is_testable
                or any(effect[column].strip() for column in numeric_columns)
                or interval_type != "not_applicable"
                or effect["pvalue_method"] != "not_applicable"
            ):
                raise ContractError(f"{adapter_id} terminal state carries a test or numeric effect for {uid}")
        elif effect["evidence_state"] in {"robust", "indeterminate", "tested_negative"}:
            if not is_testable:
                raise ContractError(f"{adapter_id} tested state marked untestable for {uid}")
            estimate = parse_float(effect["estimate"], f"{adapter_id}.estimate[{uid}]")
            pvalue = parse_float(effect["pvalue"], f"{adapter_id}.pvalue[{uid}]")
            padj = parse_float(effect["padj"], f"{adapter_id}.padj[{uid}]")
            assert estimate is not None and pvalue is not None and padj is not None
            std_error = parse_float(
                effect["std_error"],
                f"{adapter_id}.std_error[{uid}]",
                allow_empty=True,
            )
            matched_null_sd = parse_float(
                effect["matched_null_sd"],
                f"{adapter_id}.matched_null_sd[{uid}]",
                allow_empty=True,
            )
            if std_error is not None and std_error <= 0:
                raise ContractError(f"{adapter_id} std_error must be positive for {uid}")
            if matched_null_sd is not None and matched_null_sd <= 0:
                raise ContractError(f"{adapter_id} matched_null_sd must be positive for {uid}")
            if std_error is not None and matched_null_sd is not None:
                raise ContractError(f"{adapter_id} cannot label one dispersion as both SE and matched-null SD for {uid}")
            if std_error is None and matched_null_sd is None and not has_interval:
                raise ContractError(
                    f"{adapter_id} tested effect lacks assay-native uncertainty or matched-null dispersion for {uid}"
                )
            if not (0.0 <= pvalue <= 1.0 and 0.0 <= padj <= 1.0):
                raise ContractError(f"{adapter_id} p/q outside [0,1] for {uid}")
            if not effect["pvalue_method"].strip() or effect["pvalue_method"] == "not_applicable":
                raise ContractError(f"{adapter_id} tested effect lacks an explicit pvalue_method for {uid}")
            if interval_type == "not_applicable":
                raise ContractError(f"{adapter_id} tested effect cannot use interval_type=not_applicable for {uid}")
            if interval_type == "donor_effect_range":
                assert interval_low is not None and interval_high is not None
                if not interval_low <= estimate <= interval_high:
                    raise ContractError(f"{adapter_id} donor_effect_range does not contain the combined estimate for {uid}")
            expected_direction = "positive" if estimate > 0 else "negative" if estimate < 0 else "zero"
            if effect["direction_observed"] != expected_direction:
                raise ContractError(f"{adapter_id} direction does not match estimate for {uid}")
            if effect["descriptive_effect_direction"] != expected_direction:
                raise ContractError(f"{adapter_id} descriptive direction does not match estimate for {uid}")
            inferential_direction = effect["inferential_test_direction"]
            if inferential_direction not in {"positive", "negative", "zero"}:
                raise ContractError(f"{adapter_id} lacks an inferential test direction for {uid}")
            direction_agreement = parse_bool(
                effect["direction_agreement"], f"{adapter_id}.direction_agreement[{uid}]"
            )
            if direction_agreement != (inferential_direction == expected_direction):
                raise ContractError(f"{adapter_id} direction-agreement flag is inconsistent for {uid}")
            heterogeneity_values = (
                effect["heterogeneity_statistic"].strip(),
                effect["heterogeneity_df"].strip(),
                effect["heterogeneity_pvalue"].strip(),
            )
            if any(heterogeneity_values) and not all(heterogeneity_values):
                raise ContractError(f"{adapter_id} heterogeneity fields must be all present or all blank for {uid}")
            if all(heterogeneity_values):
                h_stat = parse_float(effect["heterogeneity_statistic"], f"{adapter_id}.heterogeneity_statistic[{uid}]")
                h_df = parse_int(effect["heterogeneity_df"], f"{adapter_id}.heterogeneity_df[{uid}]")
                h_p = parse_float(effect["heterogeneity_pvalue"], f"{adapter_id}.heterogeneity_pvalue[{uid}]")
                assert h_stat is not None and h_df is not None and h_p is not None
                if h_stat < 0 or h_df <= 0 or not 0.0 <= h_p <= 1.0:
                    raise ContractError(f"{adapter_id} invalid heterogeneity summary for {uid}")
            robustness = parse_bool(effect["robustness_pass"], f"{adapter_id}.robustness_pass[{uid}]")
            sign_agree = parse_bool(effect["sensitivity_sign_agree"], f"{adapter_id}.sign_agree[{uid}]")
            if effect["evidence_state"] == "robust" and not (robustness and sign_agree):
                raise ContractError(f"{adapter_id} robust state lacks robustness/sign agreement for {uid}")
            if effect["evidence_state"] == "tested_negative" and robustness:
                raise ContractError(f"{adapter_id} tested_negative row has robustness_pass=TRUE for {uid}")
            negative_rule = effect["negative_call_rule_id"].strip()
            if effect["evidence_state"] == "tested_negative" and not negative_rule:
                raise ContractError(
                    f"{adapter_id} tested_negative requires an explicit adequate-negative rule for {uid}"
                )
            if effect["evidence_state"] != "tested_negative" and negative_rule:
                raise ContractError(f"{adapter_id} non-negative row carries negative_call_rule_id for {uid}")
            if effect["evidence_state"] == "indeterminate" and robustness:
                raise ContractError(f"{adapter_id} indeterminate row has robustness_pass=TRUE for {uid}")
        elif effect["evidence_state"] == "source_dependent":
            if adapter["source_dependence"] == "independent":
                raise ContractError(f"{adapter_id} source_dependent state conflicts with independent source label")

    sensitivity_rows = tables["sensitivity.tsv"]
    for row in sensitivity_rows:
        uid = row["program_uid"]
        if uid not in expected_uids:
            raise ContractError(f"{adapter_id} sensitivity references out-of-contract program {uid}")
        for column, expected in (
            ("release_id", seal.release_id),
            ("registry_sha256", seal.registry_sha256),
            ("dataset", adapter["dataset"]),
            ("analysis_set_id", adapter["analysis_set_id"]),
            ("membership_sha256", universe[uid]["membership_sha256"]),
            ("effect_unit", adapter["effect_unit"]),
        ):
            _assert_same(row[column], expected, f"{adapter_id}.sensitivity.{column}")
        if row["estimate"].strip():
            parse_float(row["estimate"], f"{adapter_id}.sensitivity.estimate[{uid}]")

    checks.append(Check("state_semantics", "pass", f"{adapter_id}: allowed states and missingness preserved"))
    checks.append(Check("assay_native_units", "pass", f"{adapter_id}: {adapter['effect_unit']}; cross-assay comparable=FALSE"))
    checks.append(Check("registry_ready_linkage", "pass", f"{adapter_id}: registry and READY hashes match"))
    return checks


def validate_candidate_contract(
    candidate_root: Path,
    adapter_registry_path: Path,
    hotspot_root: Path,
    project_root: Path,
) -> list[Check]:
    seal = validate_hotspot_seal(hotspot_root)
    _, adapters = read_tsv(adapter_registry_path, ADAPTER_REGISTRY_COLUMNS)
    if not adapters:
        raise ContractError("adapter_registry.tsv is empty")
    adapter_ids: set[str] = set()
    checks = [Check("hotspot_ready", "pass", f"sealed registry {seal.registry_sha256}; {len(seal.registry_rows)} rows")]
    effect_units_by_assay: dict[tuple[str, str], str] = {}
    for adapter in adapters:
        adapter_id = adapter["adapter_id"]
        if not adapter_id or adapter_id in adapter_ids:
            raise ContractError(f"duplicate or empty adapter_id: {adapter_id!r}")
        adapter_ids.add(adapter_id)
        key = (adapter["dataset"], adapter["assay"])
        prior = effect_units_by_assay.get(key)
        if prior is not None and prior != adapter["effect_unit"]:
            raise ContractError(f"dataset/assay {key} uses conflicting effect units")
        effect_units_by_assay[key] = adapter["effect_unit"]
        checks.extend(validate_adapter(candidate_root, adapter, seal, project_root))
    checks.append(Check("no_universal_score", "pass", "no composite/rank columns; all adapter effects remain assay-native"))
    return checks
