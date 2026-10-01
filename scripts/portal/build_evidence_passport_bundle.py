#!/usr/bin/env python3
"""Build PASS-02--PASS-05 MASLD Gene Catalog candidate products.

Legacy ``passport_*`` identifiers remain part of the sealed v1 file contract.

The builder is intentionally a pure adapter/reshaping layer.  It consumes only
artifacts named by a signed PASS-00 selection, preserves source calls verbatim,
and refuses to run in production until the coordinator attestation and terminal
GEN/Plan13 gates pass.  PASS-06 validation and its terminal seal are owned by
``validate_evidence_passports.py``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import tempfile
from collections import Counter, OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd
import pyarrow as pa

import adjudicate_gen_ensembl_identity as gen_identity
import build_evidence_passport_ui as passport_ui
import generate_evidence_passports as contract
import render_evidence_passport_ui as passport_ui_renderer


SCRIPT_PATH = Path(__file__).resolve()
PORTAL_SCRIPT_DIR = SCRIPT_PATH.parent
BUNDLE_LOGICAL_PRODUCER_ID = "scripts/portal/build_evidence_passport_bundle.py"
SELECTION_LOGICAL_PRODUCER_ID = "scripts/portal/prepare_passport_selection.py"
ATTESTATION_LOGICAL_PRODUCER_ID = "scripts/portal/prepare_passport_attestations.py"
VALIDATOR_LOGICAL_PRODUCER_ID = "scripts/portal/validate_evidence_passports.py"
GEN_IDENTITY_LOGICAL_PRODUCER_ID = "scripts/portal/adjudicate_gen_ensembl_identity.py"
UI_GENERATOR_LOGICAL_PRODUCER_ID = "scripts/portal/build_evidence_passport_ui.py"
UI_RENDERER_LOGICAL_PRODUCER_ID = "scripts/portal/render_evidence_passport_ui.py"
SELECTION_PRODUCER_PATH = PORTAL_SCRIPT_DIR / "prepare_passport_selection.py"
ATTESTATION_PRODUCER_PATH = PORTAL_SCRIPT_DIR / "prepare_passport_attestations.py"
VALIDATOR_PRODUCER_PATH = PORTAL_SCRIPT_DIR / "validate_evidence_passports.py"
GEN_IDENTITY_PRODUCER_PATH = PORTAL_SCRIPT_DIR / "adjudicate_gen_ensembl_identity.py"
UI_GENERATOR_PRODUCER_PATH = PORTAL_SCRIPT_DIR / "build_evidence_passport_ui.py"
UI_RENDERER_PRODUCER_PATH = PORTAL_SCRIPT_DIR / "render_evidence_passport_ui.py"
ENSG_PATTERN = re.compile(r"^ENSG\d{11}$")
TRUE_VALUES = {"true", "t", "1", "yes"}
FALSE_VALUES = {"false", "f", "0", "no"}
GEN_CLASS_MAP = {
    "convergent": "concordant",
    "genetic_only": "genetically_anchored",
    "disease_state_only": "established_state_associated",
    "neither": "unresolved",
    "indeterminate_not_jointly_testable": "untested",
    # 01_freeze_and_rederive.py --untestable-state: every Tier-1/2 SuSiE signal pair
    # failed coloc's shared-posterior check. Untested, never a negative.
    "genetic_untestable_shared_posterior": "untested",
}
# 01_freeze_and_rederive.py --state-rule canonical writes this rule with bulk_padj.
CANONICAL_STATE_RULE = "canonical_padj0.05_absLFC0.5"

GEN_REQUIRED_COLUMNS = {
    "gene_symbol", "ensembl_bulk", "ensembl_genetic", "bulk_tested",
    "primary_genetic_map_tested", "primary_genetic",
    "established_state_associated", "static_class", "bulk_logFC", "bulk_t",
    "bulk_AveExpr", "bulk_treat_fdr", "coloc_best_susie_pp4",
    "coloc_best_abf_pp4", "driving_gwas", "driving_trait",
}
# Method-specific COLOC drivers (rebuild_tier12.py, review item 5). The legacy
# driving_gwas/driving_trait name the SuSiE driver whenever a SuSiE PP4 exists.
GEN_ABF_DRIVER_COLUMNS = ("abf_driving_gwas", "abf_driving_trait")
GEN_ADJUDICATED_REQUIRED_COLUMNS = {
    *GEN_REQUIRED_COLUMNS,
    "adjudicated_ensembl_id", "adjudicated_symbol",
    "identity_resolution_rule", "source_row_count", "source_row_ids",
    "source_row_sha256s", "source_rows_json", "source_call_sha256",
    "identity_adjudication_status",
}
GENE_LEVEL_POPULATED_DOMAINS = {"genetics", "transcriptomics"}
HOTSPOT_SEMANTIC_COLUMNS = {
    "program_uid", "cell_type", "module_name", "membership_sha256",
    "adjudicated_state", "robust_display", "external_test_eligible",
    "tested_negative_authorized", "evidence_scope", "adjudication_reason",
}
HOTSPOT_MEMBERSHIP_COLUMNS = {
    "program_uid", "canonical_gene", "mapped_symbol", "mapped_symbol_status",
    "source_weight", "original_l1_weight", "membership_sha256",
}
PLAN13_REQUIRED_COLUMNS = {
    "release_id", "membership_sha256", "program_uid", "program_label",
    "cell_type", "dataset", "assay", "source_dependence", "biological_unit", "n_biological",
    "n_technical", "contrast_or_exposure", "effect_unit", "estimate",
    "std_error", "interval_low", "interval_high", "pvalue", "padj",
    "testable", "testability_reason", "direction_observed", "robustness_pass",
    "evidence_state", "source_effect_row_sha256",
}
ACCEPTED_GENE_COLUMNS = {
    "source_input_row_id", "ensembl_id", "symbol", "primary_evidence_class",
    "evidence_domain", "assay", "dataset_id", "phenotype", "context",
    "biological_unit", "contrast_or_exposure", "effect_unit", "estimate",
    "standard_error", "ci_lower", "ci_upper", "p_value", "q_value",
    "direction", "testability_state", "testability_reason", "call_state",
    "negative_call_rule_id", "negative_decision_boundary", "negative_margin",
    "negative_call_passed", "gate_id", "n_biological_units",
    "n_technical_units", "allowed_wording", "limitation",
}
PLAN13_NEGATIVE_COLUMNS = {
    "negative_call_rule_id", "negative_decision_boundary", "negative_margin",
    "negative_call_passed",
}


def _stable_id(prefix: str, *parts: str) -> str:
    digest = contract.stable_sha256({"parts": [str(value) for value in parts]})[:20]
    return f"{prefix}:{digest}"


def _clean(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _optional_float(value: Any, field_name: str) -> float | None:
    text = _clean(value)
    if not text:
        return None
    try:
        result = float(text)
    except ValueError as exc:
        raise contract.PassportContractError(
            "NONNUMERIC_SOURCE_VALUE", f"{field_name}={text!r}"
        ) from exc
    if not math.isfinite(result):
        raise contract.PassportContractError(
            "NONFINITE_SOURCE_VALUE", f"{field_name}={text!r}"
        )
    return result


def _optional_int(value: Any, field_name: str) -> int | None:
    parsed = _optional_float(value, field_name)
    if parsed is None:
        return None
    if parsed < 0 or int(parsed) != parsed:
        raise contract.PassportContractError(
            "INVALID_SOURCE_COUNT", f"{field_name}={parsed}"
        )
    return int(parsed)


def _bool(value: Any, field_name: str) -> bool:
    text = _clean(value).lower()
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    raise contract.PassportContractError("INVALID_SOURCE_BOOLEAN", f"{field_name}={value!r}")


def _optional_bool(value: Any, field_name: str) -> bool | None:
    if not _clean(value):
        return None
    return _bool(value, field_name)


def _split_values(value: str) -> list[str]:
    return sorted({part.strip() for part in str(value).split(";") if part.strip()})


def _bulk_state_rule(source: Mapping[str, Any]) -> tuple[str, str]:
    """Assay label and q-value column of the frozen established-state rule."""
    if _clean(source.get("established_state_rule")) == CANONICAL_STATE_RULE:
        return "limma_voom_qw_canonical", "bulk_padj"
    return "limma_voom_qw_TREAT", "bulk_treat_fdr"


def _genetic_call_state(untestable: bool, tested: bool, positive: bool) -> tuple[str, str]:
    """Call and testability state of a GEN genetics record."""
    if untestable:
        return "untestable", "insufficient_shared_posterior"
    if positive:
        return "supported", "testable"
    return ("indeterminate", "testable") if tested else ("untestable", "underpowered_source")


def _genetic_display(
    source: Mapping[str, Any], susie: float | None, abf: float | None
) -> tuple[float | None, str, str]:
    """PP4, study, and trait of the COLOC method a GEN genetics record reports.

    SuSiE is reported when it has a PP4, unless SuSiE is not positive and ABF is
    (ABF-only positive; GNMT: SuSiE 0.128 from UKBB_ALT, ABF 0.706 from
    2021_34841290_NAFLD_EUR). The ABF value is then paired with the ABF driver.
    Source tables without the ABF driver columns keep the legacy pairing: the
    SuSiE PP4 with the legacy (SuSiE) driver, or the ABF PP4 when SuSiE is absent,
    where the legacy driver is the ABF driver.
    """
    abf_only_positive = abf is not None and abf > 0.5 and (susie is None or susie <= 0.5)
    has_abf_driver = all(column in source for column in GEN_ABF_DRIVER_COLUMNS)
    if susie is not None and not (abf_only_positive and has_abf_driver):
        return susie, _clean(source["driving_gwas"]), _clean(source["driving_trait"])
    if not has_abf_driver:
        return abf, _clean(source["driving_gwas"]), _clean(source["driving_trait"])
    study = _clean(source["abf_driving_gwas"])
    if abf is not None and not study:
        raise contract.PassportContractError(
            "GEN_ABF_DRIVER_MISSING", f"{_clean(source.get('gene_symbol'))}: ABF PP4 without abf_driving_gwas"
        )
    return abf, study, _clean(source["abf_driving_trait"])


def _read_frame(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path, engine="pyarrow").astype("string").fillna("")
    return pd.read_csv(path, sep="\t", dtype="string", keep_default_na=False)


def _require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise contract.PassportContractError("SOURCE_SCHEMA", f"{label}: missing={missing}")


def _row_hash(row: Mapping[str, Any]) -> str:
    return contract.source_row_hash(row)


@dataclass(frozen=True)
class SelectedInput:
    row: Mapping[str, str]
    path: Path

    @property
    def input_id(self) -> str:
        return self.row["input_id"]

    @property
    def sha256(self) -> str:
        return self.row["artifact_sha256"]

    @property
    def source_node_id(self) -> str:
        return _stable_id("node:derived", self.input_id, self.sha256)


@dataclass
class IdentityReference:
    by_ensembl: dict[str, str]
    by_symbol: dict[str, list[str]]

    def canonical_symbol(self, ensembl_id: str) -> str | None:
        return self.by_ensembl.get(ensembl_id)


@dataclass
class Assembly:
    analysis_release_id: str
    genes: dict[str, dict[str, Any]] = field(default_factory=dict)
    evidence: list[dict[str, Any]] = field(default_factory=list)
    coverage: list[dict[str, Any]] = field(default_factory=list)
    programs: dict[str, dict[str, Any]] = field(default_factory=dict)
    memberships: list[dict[str, Any]] = field(default_factory=list)
    program_context: list[dict[str, Any]] = field(default_factory=list)
    dataset_status: list[dict[str, Any]] = field(default_factory=list)
    assay_status: list[dict[str, Any]] = field(default_factory=list)
    quarantine: list[dict[str, Any]] = field(default_factory=list)
    gen_identity_status: list[dict[str, Any]] = field(default_factory=list)
    adapter_audit: list[dict[str, Any]] = field(default_factory=list)


def _resolve_selected_inputs(
    selection_path: Path, allow_fixture: bool
) -> tuple[list[SelectedInput], Mapping[str, Any]]:
    rows = contract.validate_signed_selection(selection_path, allow_fixture=allow_fixture)
    signature = json.loads(
        selection_path.with_name("passport_input_selection.signature.json").read_text(
            encoding="utf-8"
        )
    )
    selected = []
    for row in rows:
        path = Path(row["artifact_path"])
        if not path.is_absolute():
            path = selection_path.parent / path
        selected.append(SelectedInput(row=row, path=path.resolve()))
    return selected, signature


def _load_identity(selected: Sequence[SelectedInput]) -> IdentityReference:
    candidates = [item for item in selected if item.row["adapter_id"] == "gencode_v49_identity_v1"]
    if len(candidates) != 1:
        raise contract.PassportContractError(
            "IDENTITY_REFERENCE_CARDINALITY", "exactly one GENCODE v49 identity artifact is required"
        )
    frame = _read_frame(candidates[0].path)
    id_column = "ensembl_base" if "ensembl_base" in frame else "gene_id"
    symbol_column = "gene_name" if "gene_name" in frame else "symbol"
    _require_columns(frame, {id_column, symbol_column}, "GENCODE v49 identity")
    pairs: list[tuple[str, str]] = []
    for record in frame.to_dict("records"):
        ensembl = _clean(record[id_column]).split(".", 1)[0]
        symbol = _clean(record[symbol_column])
        if not ENSG_PATTERN.fullmatch(ensembl) or not symbol:
            continue
        pairs.append((ensembl, symbol))
    by_ensembl: dict[str, str] = {}
    for ensembl, symbol in pairs:
        if ensembl in by_ensembl and by_ensembl[ensembl] != symbol:
            raise contract.PassportContractError(
                "IDENTITY_REFERENCE_CONFLICT", f"{ensembl}: {by_ensembl[ensembl]} vs {symbol}"
            )
        by_ensembl[ensembl] = symbol
    by_symbol: dict[str, list[str]] = {}
    for ensembl, symbol in sorted(by_ensembl.items()):
        by_symbol.setdefault(symbol, []).append(ensembl)
    if not by_ensembl:
        raise contract.PassportContractError("IDENTITY_REFERENCE_EMPTY", str(candidates[0].path))
    return IdentityReference(by_ensembl=by_ensembl, by_symbol=by_symbol)


def _quarantine(
    assembly: Assembly,
    item: SelectedInput,
    row_id: str,
    symbol: str,
    ensembl_text: str,
    candidates: Sequence[str],
    reason_code: str,
    reason: str,
    source_row_sha256: str,
) -> None:
    assembly.quarantine.append(
        {
            "quarantine_id": _stable_id("quarantine", item.input_id, row_id, reason_code),
            "input_id": item.input_id,
            "adapter_id": item.row["adapter_id"],
            "artifact_role": item.row["artifact_role"],
            "source_input_row_id": row_id,
            "source_symbol": symbol,
            "source_ensembl_text": ensembl_text,
            "normalized_ensembl_candidates": ";".join(sorted(set(candidates))),
            "reason_code": reason_code,
            "reason": reason,
            "source_artifact_sha256": item.sha256,
            "source_row_sha256": source_row_sha256,
        }
    )


def _resolve_ensembl(
    assembly: Assembly,
    identity: IdentityReference,
    item: SelectedInput,
    row_id: str,
    symbol: str,
    ensembl_values: Iterable[str],
    source_row_sha256: str,
) -> str | None:
    values: set[str] = set()
    raw_values = []
    for value in ensembl_values:
        raw = _clean(value)
        if not raw:
            continue
        raw_values.append(raw)
        if any(separator in raw for separator in (";", ",", "|")):
            split = re.split(r"[;,|]", raw)
        else:
            split = [raw]
        values.update(part.strip().split(".", 1)[0] for part in split if part.strip())
    if not values and symbol:
        values.update(identity.by_symbol.get(symbol, []))
    if len(values) != 1:
        _quarantine(
            assembly, item, row_id, symbol, ";".join(raw_values), sorted(values),
            "MULTI_ENSEMBL_AMBIGUITY" if len(values) > 1 else "ENSEMBL_UNRESOLVED",
            "Source row must resolve to exactly one GENCODE v49 stable Ensembl ID; rows are never split by symbol.",
            source_row_sha256,
        )
        return None
    ensembl = next(iter(values))
    if not ENSG_PATTERN.fullmatch(ensembl) or ensembl not in identity.by_ensembl:
        _quarantine(
            assembly, item, row_id, symbol, ";".join(raw_values), [ensembl],
            "ENSEMBL_NOT_IN_GENCODE_V49",
            "Stable ID is invalid or absent from the selected GENCODE v49 identity artifact.",
            source_row_sha256,
        )
        return None
    return ensembl


def _merge_text(left: str, right: str) -> str:
    return " | ".join(sorted({value for value in (left, right) if value}))


def _ensure_gene(
    assembly: Assembly,
    identity: IdentityReference,
    ensembl: str,
    primary_class: str,
    allowed_wording: str,
    limitation: str,
) -> dict[str, Any]:
    if primary_class not in contract.PRIMARY_EVIDENCE_CLASSES:
        raise contract.PassportContractError("UNKNOWN_PRIMARY_CLASS", primary_class)
    symbol = identity.canonical_symbol(ensembl)
    if symbol is None:
        raise contract.PassportContractError("ENSEMBL_NOT_IN_IDENTITY", ensembl)
    if ensembl in assembly.genes:
        gene = assembly.genes[ensembl]
        if gene["primary_evidence_class"] != primary_class:
            raise contract.PassportContractError(
                "CONFLICTING_PRIMARY_CLASS",
                f"{ensembl}: {gene['primary_evidence_class']} vs {primary_class}",
            )
        gene["allowed_wording"] = _merge_text(gene["allowed_wording"], allowed_wording)
        gene["limitation"] = _merge_text(gene["limitation"], limitation)
        return gene
    role, role_rule, experiment_rule = contract.route_hypothesis(primary_class)
    passport_id = f"{assembly.analysis_release_id}:{ensembl}"
    gene = {
        "analysis_release_id": assembly.analysis_release_id,
        "passport_id": passport_id,
        "ensembl_id": ensembl,
        "symbol": symbol,
        "symbol_collision": False,
        "primary_evidence_class": primary_class,
        "role_hypothesis": role,
        "role_rule_id": role_rule,
        "next_experiment_rule_id": experiment_rule,
        "allowed_wording": allowed_wording,
        "limitation": limitation,
    }
    assembly.genes[ensembl] = gene
    return gene


def _validate_negative_fields(record: Mapping[str, Any], call_state: str) -> tuple[Any, Any, Any, Any]:
    rule = _clean(record.get("negative_call_rule_id")) or None
    boundary = _optional_float(record.get("negative_decision_boundary"), "negative_decision_boundary")
    margin = _optional_float(record.get("negative_margin"), "negative_margin")
    passed = _optional_bool(record.get("negative_call_passed"), "negative_call_passed")
    if call_state == "tested_negative":
        if not rule or boundary is None or margin is None or margin < 0 or passed is not True:
            raise contract.PassportContractError(
                "TESTED_NEGATIVE_RULE_REQUIRED",
                "tested_negative requires explicit rule ID, decision boundary, nonnegative margin, and passing flag",
            )
    elif any(value is not None for value in (rule, boundary, margin, passed)):
        raise contract.PassportContractError(
            "NEGATIVE_METADATA_WITHOUT_CALL", "negative metadata is reserved for tested_negative"
        )
    return rule, boundary, margin, passed


def _add_gene_evidence(
    assembly: Assembly,
    item: SelectedInput,
    source_record: Mapping[str, Any],
    gene: Mapping[str, Any],
    values: Mapping[str, Any],
) -> None:
    call_state = str(values["call_state"])
    testability = str(values["testability_state"])
    if call_state not in contract.CALL_STATES or testability not in contract.TESTABILITY_STATES:
        raise contract.PassportContractError(
            "SOURCE_CALL_VOCABULARY", f"call={call_state}; testability={testability}"
        )
    negative = _validate_negative_fields(values, call_state)
    if call_state == "tested_negative" and testability != "testable":
        raise contract.PassportContractError(
            "TESTED_NEGATIVE_REQUIRES_TESTABLE", str(values.get("source_input_row_id", ""))
        )
    source_row_id = str(values["source_input_row_id"])
    source_hash = _row_hash(source_record)
    evidence_id = _stable_id("evidence", item.input_id, source_row_id, str(values["assay"]), str(values["dataset_id"]))
    provenance = item.row["provenance_state"]
    evidence = {
        "analysis_release_id": assembly.analysis_release_id,
        "evidence_result_id": evidence_id,
        "passport_id": gene["passport_id"],
        "ensembl_id": gene["ensembl_id"],
        "symbol": gene["symbol"],
        "evidence_domain": values["evidence_domain"],
        "assay": values["assay"],
        "dataset_id": values["dataset_id"],
        "source_node_id": item.source_node_id,
        "source_release_id": item.row["source_release_id"],
        "phenotype": values["phenotype"],
        "context": values["context"],
        "biological_unit": values["biological_unit"],
        "contrast_or_exposure": values["contrast_or_exposure"],
        "effect_unit": values["effect_unit"],
        "estimate": values.get("estimate"),
        "standard_error": values.get("standard_error"),
        "ci_lower": values.get("ci_lower"),
        "ci_upper": values.get("ci_upper"),
        "p_value": values.get("p_value"),
        "q_value": values.get("q_value"),
        "direction": values["direction"],
        "testability_state": testability,
        "testability_reason": values["testability_reason"],
        "call_state": call_state,
        "negative_call_rule_id": negative[0],
        "negative_decision_boundary": negative[1],
        "negative_margin": negative[2],
        "negative_call_passed": negative[3],
        "gate_id": values["gate_id"],
        "n_biological_units": values.get("n_biological_units"),
        "n_technical_units": values.get("n_technical_units"),
        "provenance_state": provenance,
        "source_dependent": provenance != "independent",
        "source_artifact_sha256": item.sha256,
        "source_input_row_id": source_row_id,
        "source_call_state": call_state,
        "source_testability_state": testability,
        "source_provenance_state": provenance,
        "source_negative_call_rule_id": negative[0],
        "source_negative_decision_boundary": negative[1],
        "source_negative_margin": negative[2],
        "source_negative_call_passed": negative[3],
        "source_row_sha256": source_hash,
        "allowed_wording": values["allowed_wording"],
        "limitation": values["limitation"],
    }
    assembly.evidence.append(evidence)
    assembly.coverage.append(
        {
            "analysis_release_id": assembly.analysis_release_id,
            "coverage_result_id": _stable_id("coverage", evidence_id),
            "passport_id": gene["passport_id"],
            "ensembl_id": gene["ensembl_id"],
            "symbol": gene["symbol"],
            "evidence_domain": values["evidence_domain"],
            "assay": values["assay"],
            "dataset_id": values["dataset_id"],
            "source_release_id": item.row["source_release_id"],
            "testability_state": testability,
            "testability_reason": values["testability_reason"],
            "call_state": call_state,
            "n_biological_units": values.get("n_biological_units"),
            "n_technical_units": values.get("n_technical_units"),
            "coverage_denominator": values["biological_unit"],
            "limitation": values["limitation"],
        }
    )


def _base_audit(item: SelectedInput, frame: pd.DataFrame) -> dict[str, Any]:
    return {
        "input_id": item.input_id,
        "adapter_id": item.row["adapter_id"],
        "artifact_role": item.row["artifact_role"],
        "artifact_grain": item.row["artifact_grain"],
        "source_rows": len(frame),
        "gene_rows": 0,
        "evidence_rows": 0,
        "coverage_rows": 0,
        "program_rows": 0,
        "membership_rows": 0,
        "dataset_status_rows": 0,
        "assay_status_rows": 0,
        "quarantine_rows": 0,
        "status": "pass",
        "reason": "source-faithful deterministic adapter completed",
    }


def adapt_gen_classes(
    assembly: Assembly,
    identity: IdentityReference,
    item: SelectedInput,
    adjudicated: bool = False,
) -> None:
    frame = _read_frame(item.path)
    _require_columns(
        frame,
        GEN_ADJUDICATED_REQUIRED_COLUMNS if adjudicated else GEN_REQUIRED_COLUMNS,
        item.row["artifact_role"],
    )
    audit = _base_audit(item, frame)
    before = (len(assembly.genes), len(assembly.evidence), len(assembly.coverage), len(assembly.quarantine))
    if adjudicated and frame["adjudicated_ensembl_id"].duplicated().any():
        raise contract.PassportContractError(
            "GEN_ADJUDICATED_DUPLICATE_ENSEMBL",
            "adjudicated GEN input must contain exactly one row per stable Ensembl ID",
        )
    for index, source in enumerate(frame.to_dict("records"), start=1):
        row_id = f"{item.input_id}:GEN:{index:06d}"
        row_sha = _row_hash(source)
        symbol = _clean(source["gene_symbol"])
        if adjudicated:
            ensembl = _clean(source["adjudicated_ensembl_id"])
            if (
                _clean(source["identity_adjudication_status"]) != "resolved_v49"
                or not ENSG_PATTERN.fullmatch(ensembl)
                or ensembl not in identity.by_ensembl
                or _clean(source["adjudicated_symbol"]) != identity.by_ensembl[ensembl]
            ):
                raise contract.PassportContractError(
                    "GEN_ADJUDICATED_IDENTITY_INVALID", f"row={row_id};ensembl={ensembl}"
                )
            try:
                source_rows = json.loads(_clean(source["source_rows_json"]))
            except json.JSONDecodeError as exc:
                raise contract.PassportContractError(
                    "GEN_ADJUDICATED_SOURCE_JSON", row_id
                ) from exc
            if (
                not isinstance(source_rows, list)
                or len(source_rows) != _optional_int(source["source_row_count"], "source_row_count")
                or not source_rows
            ):
                raise contract.PassportContractError(
                    "GEN_ADJUDICATED_SOURCE_CARDINALITY", row_id
                )
        else:
            ensembl = _resolve_ensembl(
                assembly, identity, item, row_id, symbol,
                [source["ensembl_bulk"], source["ensembl_genetic"]], row_sha,
            )
        if ensembl is None:
            continue
        static_class = _clean(source["static_class"])
        if static_class not in GEN_CLASS_MAP:
            raise contract.PassportContractError("GEN_STATIC_CLASS", static_class)
        primary_class = GEN_CLASS_MAP[static_class]
        genetic_untestable = static_class == "genetic_untestable_shared_posterior"
        gene = _ensure_gene(
            assembly, identity, ensembl, primary_class,
            "Every Tier-1/2 SuSiE-COLOC signal pair failed the shared-posterior check; genetic evidence is untested, not negative."
            if genetic_untestable else
            {
                "concordant": "Frozen genetics/state interface class is concordant.",
                "genetically_anchored": "Frozen static genetic evidence is present.",
                "established_state_associated": "Frozen established-state association is present.",
                "untested": "The gene was not jointly testable across the frozen interface.",
                "unresolved": "Neither frozen positive criterion passed; this is not a tested negative.",
            }[primary_class],
            "GEN terminal closure forbids context rescue, source-negative claims, and powered-null claims.",
        )
        bulk_tested = _bool(source["bulk_tested"], "bulk_tested")
        bulk_positive = _bool(source["established_state_associated"], "established_state_associated")
        bulk_call = "supported" if bulk_positive else "indeterminate" if bulk_tested else "untestable"
        bulk_assay, bulk_q_field = _bulk_state_rule(source)
        bulk_values = {
            "source_input_row_id": f"{row_id}:bulk",
            "evidence_domain": "transcriptomics", "assay": bulk_assay,
            "dataset_id": item.row["source_datasets"] or "canonical_human_bulk",
            "phenotype": "MASLD established disease state", "context": "pooled human liver",
            "biological_unit": item.row["biological_unit"],
            "contrast_or_exposure": "disease versus control",
            "effect_unit": "log2 fold change",
            "estimate": _optional_float(source["bulk_logFC"], "bulk_logFC"),
            "standard_error": None, "ci_lower": None, "ci_upper": None,
            "p_value": None, "q_value": _optional_float(source[bulk_q_field], bulk_q_field),
            "direction": (
                "positive" if (_optional_float(source["bulk_logFC"], "bulk_logFC") or 0) > 0
                else "negative" if (_optional_float(source["bulk_logFC"], "bulk_logFC") or 0) < 0
                else "unknown"
            ),
            "testability_state": "testable" if bulk_tested else "not_detected",
            "testability_reason": "testable" if bulk_tested else "not_detected",
            "call_state": bulk_call,
            "negative_call_rule_id": None, "negative_decision_boundary": None,
            "negative_margin": None, "negative_call_passed": None,
            "gate_id": "GEN_FROZEN_STATIC_CLASS",
            "n_biological_units": _optional_int(
                item.row["n_biological_units"], "selection.n_biological_units"
            ),
            "n_technical_units": _optional_int(
                item.row["n_technical_units"], "selection.n_technical_units"
            ),
            "allowed_wording": "Established-state associated." if bulk_positive else "No frozen positive state call; not an informative negative.",
            "limitation": "Cross-sectional pooled-liver association; source cohort denominators remain in the source release.",
        }
        _add_gene_evidence(assembly, item, source, gene, bulk_values)
        genetic_tested = _bool(source["primary_genetic_map_tested"], "primary_genetic_map_tested")
        genetic_positive = _bool(source["primary_genetic"], "primary_genetic")
        genetic_call, genetic_testability = _genetic_call_state(
            genetic_untestable, genetic_tested, genetic_positive
        )
        susie = _optional_float(source["coloc_best_susie_pp4"], "coloc_best_susie_pp4")
        abf = _optional_float(source["coloc_best_abf_pp4"], "coloc_best_abf_pp4")
        genetic_pp4, driving_gwas, driving_trait = _genetic_display(source, susie, abf)
        if genetic_call == "untestable":
            # An untestable record carries no statistic (validator rule
            # FABRICATED_UNTESTABLE_STATISTIC); the best PP4 of another study, or
            # the ABF value, would read as this record's test result.
            genetic_pp4, driving_gwas, driving_trait = None, "", ""
        genetic_values = {
            "source_input_row_id": f"{row_id}:genetics",
            "evidence_domain": "genetics", "assay": "SuSiE_COLOC_or_ABF_fallback",
            "dataset_id": driving_gwas or item.row["source_datasets"] or "GEN_frozen_interface",
            "phenotype": driving_trait or "source-defined genetic phenotype",
            "context": "bulk liver cis-eQTL interface",
            "biological_unit": item.row["biological_unit"],
            "contrast_or_exposure": "trait-locus colocalization",
            "effect_unit": "PP.H4", "estimate": genetic_pp4,
            "standard_error": None, "ci_lower": None, "ci_upper": None,
            "p_value": None, "q_value": None, "direction": "not_directional",
            "testability_state": genetic_testability,
            "testability_reason": genetic_testability,
            "call_state": genetic_call,
            "negative_call_rule_id": None, "negative_decision_boundary": None,
            "negative_margin": None, "negative_call_passed": None,
            "gate_id": "GEN_FROZEN_STATIC_CLASS", "n_biological_units": None,
            "n_technical_units": None,
            "allowed_wording": (
                "Genetically anchored." if genetic_positive
                else "SuSiE-COLOC untestable: no signal pair kept enough shared posterior; not a negative."
                if genetic_untestable
                else "No frozen positive genetic call; not a powered genetic null."
            ),
            "limitation": "Phenotype provenance and ancestry/eQTL-panel limitations remain source-specific.",
        }
        _add_gene_evidence(assembly, item, source, gene, genetic_values)
    after = (len(assembly.genes), len(assembly.evidence), len(assembly.coverage), len(assembly.quarantine))
    audit.update(
        gene_rows=after[0] - before[0], evidence_rows=after[1] - before[1],
        coverage_rows=after[2] - before[2], quarantine_rows=after[3] - before[3]
    )
    assembly.adapter_audit.append(audit)


def adapt_gen_identity_status(assembly: Assembly, item: SelectedInput) -> None:
    frame = _read_frame(item.path)
    _require_columns(frame, set(gen_identity.READY_COLUMNS), item.row["artifact_role"])
    if len(frame) != 1:
        raise contract.PassportContractError(
            "GEN_IDENTITY_READY_CARDINALITY", f"{item.input_id}:{len(frame)}"
        )
    if assembly.gen_identity_status:
        raise contract.PassportContractError(
            "GEN_IDENTITY_READY_CARDINALITY", "more than one selected identity READY"
        )
    assembly.gen_identity_status.append(
        {column: _clean(frame.iloc[0][column]) for column in gen_identity.READY_COLUMNS}
    )
    audit = _base_audit(item, frame)
    audit["reason"] = "outcome-blind GENCODE-v49 identity adjudication terminal imported"
    assembly.adapter_audit.append(audit)


def adapt_gen_identity_quarantine(assembly: Assembly, item: SelectedInput) -> None:
    frame = _read_frame(item.path)
    _require_columns(frame, set(gen_identity.QUARANTINE_COLUMNS), item.row["artifact_role"])
    audit = _base_audit(item, frame)
    allowed_codes = {
        "GEN_IRREDUCIBLE_MULTI_ENSEMBL",
        "GEN_ENSEMBL_ABSENT_V49",
        "GEN_IDENTITY_CONFLICT",
        "GEN_ENSEMBL_UNRESOLVED",
        "GEN_INVALID_ENSEMBL_TOKEN",
    }
    for source in frame.to_dict("records"):
        code = _clean(source["reason_code"])
        if code not in allowed_codes:
            raise contract.PassportContractError(
                "GEN_IDENTITY_QUARANTINE_REASON", code
            )
        row_id = _clean(source["source_row_id"])
        assembly.quarantine.append(
            {
                "quarantine_id": _stable_id(
                    "quarantine", item.input_id, row_id, code,
                    _clean(source["source_row_sha256"]),
                ),
                "input_id": item.input_id,
                "adapter_id": item.row["adapter_id"],
                "artifact_role": item.row["artifact_role"],
                "source_input_row_id": row_id,
                "source_symbol": _clean(source["source_symbol"]),
                "source_ensembl_text": ";".join(
                    value
                    for value in (
                        _clean(source["source_ensembl_bulk"]),
                        _clean(source["source_ensembl_genetic"]),
                    )
                    if value
                ),
                "normalized_ensembl_candidates": _clean(
                    source["normalized_ensembl_candidates"]
                ),
                "reason_code": code,
                "reason": _clean(source["reason"]),
                "source_artifact_sha256": item.sha256,
                "source_row_sha256": _clean(source["source_row_sha256"]),
            }
        )
    audit["quarantine_rows"] = len(frame)
    audit["reason"] = "upstream identity exclusions imported for visible audit; no gene calls created"
    assembly.adapter_audit.append(audit)


def adapt_accepted_gene_evidence(
    assembly: Assembly, identity: IdentityReference, item: SelectedInput
) -> None:
    frame = _read_frame(item.path)
    _require_columns(frame, ACCEPTED_GENE_COLUMNS, item.row["artifact_role"])
    audit = _base_audit(item, frame)
    before = (len(assembly.genes), len(assembly.evidence), len(assembly.coverage), len(assembly.quarantine))
    for source in frame.to_dict("records"):
        evidence_domain = _clean(source["evidence_domain"])
        if evidence_domain not in GENE_LEVEL_POPULATED_DOMAINS:
            raise contract.PassportContractError(
                "GENE_DOMAIN_OUT_OF_CANDIDATE_SCOPE",
                f"{item.input_id}:{evidence_domain}; only genetics and transcriptomics may populate gene-grain evidence",
            )
        row_id = _clean(source["source_input_row_id"])
        row_sha = _row_hash(source)
        ensembl = _resolve_ensembl(
            assembly, identity, item, row_id, _clean(source["symbol"]),
            [_clean(source["ensembl_id"])], row_sha,
        )
        if ensembl is None:
            continue
        gene = _ensure_gene(
            assembly, identity, ensembl, _clean(source["primary_evidence_class"]),
            _clean(source["allowed_wording"]), _clean(source["limitation"]),
        )
        values = {
            **source,
            "source_input_row_id": f"{item.input_id}:{row_id}",
            "estimate": _optional_float(source["estimate"], "estimate"),
            "standard_error": _optional_float(source["standard_error"], "standard_error"),
            "ci_lower": _optional_float(source["ci_lower"], "ci_lower"),
            "ci_upper": _optional_float(source["ci_upper"], "ci_upper"),
            "p_value": _optional_float(source["p_value"], "p_value"),
            "q_value": _optional_float(source["q_value"], "q_value"),
            "n_biological_units": _optional_int(source["n_biological_units"], "n_biological_units"),
            "n_technical_units": _optional_int(source["n_technical_units"], "n_technical_units"),
        }
        _add_gene_evidence(assembly, item, source, gene, values)
    after = (len(assembly.genes), len(assembly.evidence), len(assembly.coverage), len(assembly.quarantine))
    audit.update(
        gene_rows=after[0] - before[0], evidence_rows=after[1] - before[1],
        coverage_rows=after[2] - before[2], quarantine_rows=after[3] - before[3]
    )
    assembly.adapter_audit.append(audit)


def adapt_hotspot_semantics(assembly: Assembly, item: SelectedInput) -> None:
    frame = _read_frame(item.path)
    _require_columns(frame, HOTSPOT_SEMANTIC_COLUMNS, item.row["artifact_role"])
    audit = _base_audit(item, frame)
    for source in frame.to_dict("records"):
        uid = _clean(source["program_uid"])
        if uid in assembly.programs:
            raise contract.PassportContractError("DUPLICATE_PROGRAM_UID", uid)
        if _bool(source["tested_negative_authorized"], "tested_negative_authorized"):
            raise contract.PassportContractError(
                "HOTSPOT_GENE_NEGATIVE_EXPANSION",
                "Hotspot semantics authorize no tested-negative program or member-gene expansion",
            )
        assembly.programs[uid] = {
            "analysis_release_id": assembly.analysis_release_id,
            "program_uid": uid,
            "program_label": _clean(source["module_name"]),
            "cell_type": _clean(source["cell_type"]),
            "membership_sha256": _clean(source["membership_sha256"]),
            "source_release_id": item.row["source_release_id"],
            "program_call_scope": "program_membership_and_context_only",
            "member_gene_call_expansion_authorized": False,
        }
    audit["program_rows"] = len(frame)
    assembly.adapter_audit.append(audit)


def adapt_hotspot_membership(
    assembly: Assembly, identity: IdentityReference, item: SelectedInput
) -> None:
    frame = _read_frame(item.path)
    _require_columns(frame, HOTSPOT_MEMBERSHIP_COLUMNS, item.row["artifact_role"])
    audit = _base_audit(item, frame)
    before_q = len(assembly.quarantine)
    for index, source in enumerate(frame.to_dict("records"), start=1):
        row_id = f"{item.input_id}:MEM:{index:06d}"
        row_sha = _row_hash(source)
        uid = _clean(source["program_uid"])
        if uid not in assembly.programs:
            raise contract.PassportContractError("PROGRAM_MEMBERSHIP_ORPHAN", uid)
        mapping_status = _clean(source["mapped_symbol_status"])
        if mapping_status not in {
            "gencode_v49_unique_symbol_confirmed",
            "gencode_v49_unambiguous_ensembl_to_symbol",
        }:
            _quarantine(
                assembly, item, row_id, _clean(source["canonical_gene"]), "", [],
                "PROGRAM_MEMBERSHIP_UNMAPPED_SOURCE",
                "Frozen program membership is retained by membership hash, but this source member has no unique GENCODE v49 display link and cannot generate a gene call.",
                row_sha,
            )
            continue
        mapped_value = _clean(source["mapped_symbol"] or source["canonical_gene"])
        if mapping_status == "gencode_v49_unambiguous_ensembl_to_symbol":
            symbol = ""
            ensembl_values = [mapped_value, _clean(source["canonical_gene"])]
        else:
            symbol = mapped_value
            ensembl_values = []
        ensembl = _resolve_ensembl(
            assembly, identity, item, row_id, symbol, ensembl_values, row_sha
        )
        if ensembl is None:
            continue
        if _clean(source["membership_sha256"]) != assembly.programs[uid]["membership_sha256"]:
            raise contract.PassportContractError("PROGRAM_MEMBERSHIP_HASH", uid)
        assembly.memberships.append(
            {
                "analysis_release_id": assembly.analysis_release_id,
                "program_uid": uid,
                "ensembl_id": ensembl,
                "symbol": identity.canonical_symbol(ensembl),
                "source_weight": _optional_float(source["source_weight"], "source_weight"),
                "original_l1_weight": _optional_float(source["original_l1_weight"], "original_l1_weight"),
                "membership_sha256": _clean(source["membership_sha256"]),
                "membership_only": True,
                "gene_call_expansion_authorized": False,
                "source_input_row_id": row_id,
                "source_artifact_sha256": item.sha256,
                "source_row_sha256": row_sha,
            }
        )
    audit["membership_rows"] = len(assembly.memberships)
    audit["quarantine_rows"] = len(assembly.quarantine) - before_q
    assembly.adapter_audit.append(audit)


def _direction(value: str) -> str:
    value = value.lower()
    if value in {"positive", "up", "increased", "higher"}:
        return "positive"
    if value in {"negative", "down", "decreased", "lower"}:
        return "negative"
    if value in {"discordant"}:
        return "discordant"
    return "unknown"


def adapt_plan13_context(assembly: Assembly, item: SelectedInput) -> None:
    frame = _read_frame(item.path)
    _require_columns(frame, PLAN13_REQUIRED_COLUMNS, item.row["artifact_role"])
    if (frame["evidence_state"] == "tested_negative").any():
        missing_negative = sorted(PLAN13_NEGATIVE_COLUMNS - set(frame.columns))
        if missing_negative:
            raise contract.PassportContractError(
                "TESTED_NEGATIVE_RULE_REQUIRED",
                f"Plan13 tested_negative rows lack explicit fields: {missing_negative}",
            )
    audit = _base_audit(item, frame)
    for index, source in enumerate(frame.to_dict("records"), start=1):
        state = _clean(source["evidence_state"])
        if state == "skipped":
            continue
        uid = _clean(source["program_uid"])
        if uid not in assembly.programs:
            raise contract.PassportContractError("PLAN13_PROGRAM_ORPHAN", uid)
        if _clean(source["membership_sha256"]) != assembly.programs[uid]["membership_sha256"]:
            raise contract.PassportContractError("PLAN13_MEMBERSHIP_HASH", uid)
        call_map = {
            "robust": "supported", "source_dependent": "supported",
            "indeterminate": "indeterminate", "tested_negative": "tested_negative",
            "untestable": "untestable",
            "not_applicable": "not_applicable",
        }
        if state not in call_map:
            raise contract.PassportContractError("PLAN13_EVIDENCE_STATE", state)
        call = call_map[state]
        testable = _bool(source["testable"], "testable")
        testability = (
            "testable" if testable else
            "assay_out_of_scope" if call == "not_applicable" else
            "insufficient_gene_coverage"
        )
        negative = _validate_negative_fields(source, call)
        row_id = f"{item.input_id}:P13:{index:06d}"
        source_dependence = _clean(source["source_dependence"])
        if source_dependence == "independent":
            provenance = "independent"
        elif source_dependence == "source_dependent":
            provenance = item.row["provenance_state"]
            if provenance == "independent":
                raise contract.PassportContractError(
                    "PLAN13_PROVENANCE_UNDERSPECIFIED",
                    "selection must classify source_dependent Plan13 rows as partially_dependent or reused_source",
                )
        else:
            raise contract.PassportContractError(
                "PLAN13_SOURCE_DEPENDENCE", source_dependence
            )
        assembly.program_context.append(
            {
                "analysis_release_id": assembly.analysis_release_id,
                "program_context_id": _stable_id("program_context", item.input_id, row_id),
                "program_uid": uid,
                "program_label": assembly.programs[uid]["program_label"],
                "cell_type": assembly.programs[uid]["cell_type"],
                "dataset_id": _clean(source["dataset"]),
                "assay": _clean(source["assay"]),
                "source_node_id": item.source_node_id,
                "source_release_id": item.row["source_release_id"],
                "context": item.row["claim_scope"],
                "biological_unit": _clean(source["biological_unit"]),
                "contrast_or_exposure": _clean(source["contrast_or_exposure"]),
                "effect_unit": _clean(source["effect_unit"]),
                "estimate": _optional_float(source["estimate"], "estimate"),
                "standard_error": _optional_float(source["std_error"], "std_error"),
                "ci_lower": _optional_float(source["interval_low"], "interval_low"),
                "ci_upper": _optional_float(source["interval_high"], "interval_high"),
                "p_value": _optional_float(source["pvalue"], "pvalue"),
                "q_value": _optional_float(source["padj"], "padj"),
                "direction": _direction(_clean(source["direction_observed"])),
                "testability_state": testability,
                "testability_reason": testability,
                "call_state": call,
                "negative_call_rule_id": negative[0],
                "negative_decision_boundary": negative[1],
                "negative_margin": negative[2],
                "negative_call_passed": negative[3],
                "gate_id": "PLAN13_TERMINAL_INTEGRATION",
                "n_biological_units": _optional_int(source["n_biological"], "n_biological"),
                "n_technical_units": _optional_int(source["n_technical"], "n_technical"),
                "provenance_state": provenance,
                "source_dependent": provenance != "independent",
                "source_artifact_sha256": item.sha256,
                "source_input_row_id": row_id,
                "source_row_sha256": _clean(source["source_effect_row_sha256"]) or _row_hash(source),
                "allowed_wording": item.row["allowed_claim_wording"],
                "limitation": item.row["limitation"],
                "gene_call_expansion_authorized": False,
            }
        )
    audit["program_rows"] = len(assembly.program_context)
    assembly.adapter_audit.append(audit)


def adapt_dataset_status(assembly: Assembly, item: SelectedInput) -> None:
    frame = _read_frame(item.path)
    required = {
        "dataset", "assay", "dataset_gate", "figure4_verdict", "biological_unit",
        "n_biological", "n_technical", "interpretation", "source_dependence",
    }
    _require_columns(frame, required, item.row["artifact_role"])
    audit = _base_audit(item, frame)
    for index, source in enumerate(frame.to_dict("records"), start=1):
        row_id = f"{item.input_id}:DATASET:{index:04d}"
        gate = _clean(source["dataset_gate"])
        skipped = gate in {"skipped_source_gate", "skipped", "skipped_no_donor_key"}
        source_dependence = _clean(source["source_dependence"])
        if source_dependence == "independent":
            provenance = "independent"
        elif source_dependence == "source_dependent":
            provenance = item.row["provenance_state"]
            if provenance == "independent":
                raise contract.PassportContractError(
                    "PLAN13_PROVENANCE_UNDERSPECIFIED", _clean(source["dataset"])
                )
        else:
            raise contract.PassportContractError(
                "PLAN13_SOURCE_DEPENDENCE", source_dependence
            )
        assembly.dataset_status.append(
            {
                "dataset_status_id": _stable_id("dataset_status", item.input_id, row_id),
                "dataset_id": _clean(source["dataset"]),
                "assay": _clean(source["assay"]),
                "source_node_id": item.source_node_id,
                "source_release_id": item.row["source_release_id"],
                "dataset_gate": "skipped_source_gate" if skipped else gate,
                "status": "skipped_source_gate" if skipped else _clean(source["figure4_verdict"]),
                "reason": _clean(source["interpretation"]),
                "biological_unit": _clean(source["biological_unit"]),
                "n_biological_units": _optional_int(source["n_biological"], "n_biological"),
                "n_technical_units": _optional_int(source["n_technical"], "n_technical"),
                "provenance_state": provenance,
                "source_dependent": provenance != "independent",
                "source_artifact_sha256": item.sha256,
                "source_input_row_id": row_id,
                "source_row_sha256": _row_hash(source),
                "gene_expansion_authorized": False,
                "allowed_wording": item.row["allowed_claim_wording"],
                "limitation": item.row["limitation"],
            }
        )
    audit["dataset_status_rows"] = len(frame)
    assembly.adapter_audit.append(audit)


def _add_assay_status(
    assembly: Assembly, item: SelectedInput, source: Mapping[str, Any], row_id: str,
    grain: str, entity_id: str, status: str, call: str, testability: str,
    estimate_field: str = "", se_field: str = "", low_field: str = "",
    high_field: str = "", p_field: str = "", q_field: str = "",
    effect_unit: str = "source-native statistic",
) -> None:
    assembly.assay_status.append(
        {
            "assay_status_id": _stable_id("assay_status", item.input_id, row_id),
            "result_grain": grain, "entity_id": entity_id,
            "assay": "Myojin_HLF_palmitate_CRISPR", "dataset_id": item.row["source_datasets"],
            "source_node_id": item.source_node_id, "source_release_id": item.row["source_release_id"],
            "status": status, "call_state": call, "testability_state": testability,
            "estimate": _optional_float(source.get(estimate_field), estimate_field) if estimate_field else None,
            "standard_error": _optional_float(source.get(se_field), se_field) if se_field else None,
            "ci_lower": _optional_float(source.get(low_field), low_field) if low_field else None,
            "ci_upper": _optional_float(source.get(high_field), high_field) if high_field else None,
            "p_value": _optional_float(source.get(p_field), p_field) if p_field else None,
            "q_value": _optional_float(source.get(q_field), q_field) if q_field else None,
            "effect_unit": effect_unit, "provenance_state": item.row["provenance_state"],
            "source_dependent": item.row["provenance_state"] != "independent",
            "source_artifact_sha256": item.sha256, "source_input_row_id": row_id,
            "source_row_sha256": _row_hash(source), "gene_expansion_authorized": False,
            "allowed_wording": item.row["allowed_claim_wording"], "limitation": item.row["limitation"],
        }
    )


def adapt_myojin(assembly: Assembly, item: SelectedInput) -> None:
    frame = _read_frame(item.path)
    audit = _base_audit(item, frame)
    adapter = item.row["adapter_id"]
    if adapter == "myojin_assay_verdict_v1":
        _require_columns(frame, {"main_figure_eligible", "source_gate_pass", "verdict_reason_codes"}, item.row["artifact_role"])
        if len(frame) != 1 or _bool(frame.iloc[0]["main_figure_eligible"], "main_figure_eligible"):
            raise contract.PassportContractError("MYOJIN_NONSUPPORT_CONTRACT", "adapter accepts only terminal non-support")
        source = frame.iloc[0].to_dict()
        _add_assay_status(
            assembly, item, source, f"{item.input_id}:ASSAY", "assay", "MYOJIN_HLF",
            "nonconfirmatory", "indeterminate",
            "testable" if _bool(source["source_gate_pass"], "source_gate_pass") else "source_gate_failed",
        )
    elif adapter == "myojin_class_nonsupport_v1":
        _require_columns(frame, {"contrast", "estimate", "HC3_SE", "CI95_low", "CI95_high", "permutation_p", "BH_q"}, item.row["artifact_role"])
        for index, source in enumerate(frame.to_dict("records"), start=1):
            _add_assay_status(
                assembly, item, source, f"{item.input_id}:CLASS:{index:04d}", "class",
                _clean(source["contrast"]), "nonconfirmatory", "indeterminate", "testable",
                "estimate", "HC3_SE", "CI95_low", "CI95_high", "permutation_p", "BH_q",
                "palmitate-minus-vehicle class-model coefficient",
            )
    elif adapter == "myojin_program_nonsupport_v1":
        _require_columns(frame, {"program_uid", "testable", "testability_reason", "signed_effect", "empirical_p", "BH_q"}, item.row["artifact_role"])
        for index, source in enumerate(frame.to_dict("records"), start=1):
            is_testable = _bool(source["testable"], "testable")
            _add_assay_status(
                assembly, item, source, f"{item.input_id}:PROGRAM:{index:04d}", "program",
                _clean(source["program_uid"]), "nonconfirmatory" if is_testable else "untestable",
                "indeterminate" if is_testable else "untestable",
                "testable" if is_testable else "insufficient_gene_coverage",
                "signed_effect", "", "", "", "empirical_p", "BH_q", "matched-null signed program effect",
            )
    else:
        raise contract.PassportContractError("MYOJIN_ADAPTER", adapter)
    audit["assay_status_rows"] = len(assembly.assay_status)
    assembly.adapter_audit.append(audit)


def _validate_gen_identity_adjudication(
    selected: Sequence[SelectedInput],
    signature: Mapping[str, Any],
    allow_fixture: bool,
) -> None:
    raw_legacy = [
        item for item in selected if item.row["adapter_id"] == "gen_frozen_classes_v1"
    ]
    adjudicated = [
        item for item in selected if item.row["adapter_id"] == "gen_adjudicated_classes_v1"
    ]
    if raw_legacy and not allow_fixture:
        raise contract.PassportContractError(
            "RAW_GEN_IDENTITY_UNADJUDICATED",
            "production Gene Catalog entries require the outcome-blind GENCODE-v49 adjudicated GEN artifact",
        )
    if raw_legacy and adjudicated:
        raise contract.PassportContractError(
            "GEN_IDENTITY_INPUT_CONFLICT", "select either fixture-legacy or adjudicated GEN classes"
        )
    if not adjudicated:
        return
    adapter_groups = {
        adapter: [item for item in selected if item.row["adapter_id"] == adapter]
        for adapter in (
            "gen_raw_identity_source_v1",
            "gen_identity_adjudication_ready_v1",
            "gen_identity_adjudication_audit_v1",
            "gen_identity_adjudication_quarantine_v1",
            "gencode_v49_identity_v1",
        )
    }
    bad = {adapter: len(items) for adapter, items in adapter_groups.items() if len(items) != 1}
    if len(adjudicated) != 1:
        bad["gen_adjudicated_classes_v1"] = len(adjudicated)
    if bad:
        raise contract.PassportContractError(
            "GEN_IDENTITY_SELECTION_CARDINALITY", str(bad)
        )
    raw = adapter_groups["gen_raw_identity_source_v1"][0]
    ready_item = adapter_groups["gen_identity_adjudication_ready_v1"][0]
    audit_item = adapter_groups["gen_identity_adjudication_audit_v1"][0]
    quarantine_item = adapter_groups["gen_identity_adjudication_quarantine_v1"][0]
    identity_item = adapter_groups["gencode_v49_identity_v1"][0]
    classes_item = adjudicated[0]
    ready = _read_frame(ready_item.path)
    _require_columns(ready, set(gen_identity.READY_COLUMNS), "GEN identity READY")
    if len(ready) != 1:
        raise contract.PassportContractError(
            "GEN_IDENTITY_READY_CARDINALITY", str(len(ready))
        )
    record = ready.iloc[0]
    required_semantics = {
        "status": "identity_adjudication_complete_with_quarantine",
        "source_artifact_sha256": raw.sha256,
        "identity_artifact_sha256": identity_item.sha256,
        "adjudicated_artifact_sha256": classes_item.sha256,
        "quarantine_artifact_sha256": quarantine_item.sha256,
        "audit_artifact_sha256": audit_item.sha256,
        "producer_sha256": contract.sha256_file(gen_identity.SCRIPT_PATH),
        "outcome_fields_used": "false",
        "merge_discordant_calls": "false",
        "one_row_per_ensembl": "true",
        "canonical_promotion_authorized": "false",
        "fixture_only": str(bool(signature["fixture_only"])).lower(),
    }
    mismatches = {
        key: (_clean(record[key]), expected)
        for key, expected in required_semantics.items()
        if _clean(record[key]).lower() != expected.lower()
    }
    if _clean(record["release_id"]) != str(signature["analysis_release_id"]):
        mismatches["release_id"] = (
            _clean(record["release_id"]), str(signature["analysis_release_id"])
        )
    if mismatches:
        raise contract.PassportContractError(
            "GEN_IDENTITY_READY_STALE", str(mismatches)
        )

    raw_frame = _read_frame(raw.path)
    classes = _read_frame(classes_item.path)
    quarantine = _read_frame(quarantine_item.path)
    audit = _read_frame(audit_item.path)
    _require_columns(classes, GEN_ADJUDICATED_REQUIRED_COLUMNS, "GEN adjudicated classes")
    _require_columns(quarantine, set(gen_identity.QUARANTINE_COLUMNS), "GEN identity quarantine")
    _require_columns(audit, set(gen_identity.AUDIT_COLUMNS), "GEN identity audit")
    if len(audit) != 1:
        raise contract.PassportContractError(
            "GEN_IDENTITY_AUDIT_CARDINALITY", str(len(audit))
        )
    audit_record = audit.iloc[0]
    for field in (
        "release_id", "status", "source_artifact_sha256", "identity_artifact_sha256",
        "adjudicated_artifact_sha256", "quarantine_artifact_sha256", "producer_sha256",
        "n_source_rows", "n_adjudicated_genes", "n_quarantined_source_rows",
        "n_irreducible_multi_id", "n_absent_v49", "n_identity_conflict",
        "n_unresolved", "n_invalid_token", "outcome_fields_used",
        "merge_discordant_calls", "one_row_per_ensembl",
    ):
        if field in audit and _clean(audit_record[field]) != _clean(record[field]):
            raise contract.PassportContractError(
                "GEN_IDENTITY_AUDIT_READY_MISMATCH", field
            )
    count_expectations = {
        "n_source_rows": len(raw_frame),
        "n_adjudicated_genes": len(classes),
        "n_quarantined_source_rows": len(quarantine),
    }
    reason_counts = Counter(quarantine["reason_code"].astype(str))
    count_expectations.update(
        {
            "n_irreducible_multi_id": reason_counts["GEN_IRREDUCIBLE_MULTI_ENSEMBL"],
            "n_absent_v49": reason_counts["GEN_ENSEMBL_ABSENT_V49"],
            "n_identity_conflict": reason_counts["GEN_IDENTITY_CONFLICT"],
            "n_unresolved": reason_counts["GEN_ENSEMBL_UNRESOLVED"],
            "n_invalid_token": reason_counts["GEN_INVALID_ENSEMBL_TOKEN"],
        }
    )
    for field, expected in count_expectations.items():
        if _optional_int(record[field], field) != expected:
            raise contract.PassportContractError(
                "GEN_IDENTITY_COUNT_MISMATCH", f"{field}: expected={expected}; observed={record[field]}"
            )
    if classes["adjudicated_ensembl_id"].duplicated().any():
        raise contract.PassportContractError(
            "GEN_ADJUDICATED_DUPLICATE_ENSEMBL", "one-row-per-v49 invariant failed"
        )

    raw_records = {
        f"GEN_SOURCE:{index:06d}": {
            str(key): _clean(value) for key, value in source.items()
        }
        for index, source in enumerate(raw_frame.to_dict("records"), start=1)
    }
    accounted: list[str] = []
    for source in classes.to_dict("records"):
        row_ids = _split_values(_clean(source["source_row_ids"]))
        try:
            preserved = json.loads(_clean(source["source_rows_json"]))
        except json.JSONDecodeError as exc:
            raise contract.PassportContractError(
                "GEN_ADJUDICATED_SOURCE_JSON", _clean(source["adjudicated_ensembl_id"])
            ) from exc
        expected_rows = [raw_records[row_id] for row_id in row_ids if row_id in raw_records]
        if len(expected_rows) != len(row_ids) or preserved != expected_rows:
            raise contract.PassportContractError(
                "GEN_ADJUDICATED_SOURCE_PRESERVATION", _clean(source["adjudicated_ensembl_id"])
            )
        call_hashes = {
            gen_identity._row_hash(  # noqa: SLF001
                {key: value for key, value in row.items() if key not in gen_identity.IDENTITY_COLUMNS}
            )
            for row in expected_rows
        }
        if len(call_hashes) != 1 or _clean(source["source_call_sha256"]) not in call_hashes:
            raise contract.PassportContractError(
                "GEN_ADJUDICATED_CALL_CONFLICT", _clean(source["adjudicated_ensembl_id"])
            )
        for key, value in expected_rows[0].items():
            if _clean(source.get(key)) != value:
                raise contract.PassportContractError(
                    "GEN_ADJUDICATED_SOURCE_PRESERVATION", f"{row_ids[0]}:{key}"
                )
        accounted.extend(row_ids)
    for source in quarantine.to_dict("records"):
        row_id = _clean(source["source_row_id"])
        if row_id not in raw_records:
            raise contract.PassportContractError(
                "GEN_IDENTITY_QUARANTINE_ORPHAN", row_id
            )
        try:
            preserved = json.loads(_clean(source["source_row_json"]))
        except json.JSONDecodeError as exc:
            raise contract.PassportContractError(
                "GEN_IDENTITY_QUARANTINE_JSON", row_id
            ) from exc
        if preserved != raw_records[row_id]:
            raise contract.PassportContractError(
                "GEN_IDENTITY_QUARANTINE_PRESERVATION", row_id
            )
        accounted.append(row_id)
    if sorted(accounted) != sorted(raw_records) or len(accounted) != len(set(accounted)):
        raise contract.PassportContractError(
            "GEN_IDENTITY_SOURCE_ACCOUNTING",
            "every raw GEN row must occur exactly once in accepted source_rows_json or quarantine",
        )


def _validate_gate_artifacts(
    selected: Sequence[SelectedInput], signature: Mapping[str, Any], allow_fixture: bool
) -> None:
    _validate_gen_identity_adjudication(selected, signature, allow_fixture)
    by_id = {item.input_id: item for item in selected}
    closure = by_id.get(str(signature["gen_terminal_closure_input_id"]))
    gen_ready = by_id.get(str(signature["gen_terminal_ready_input_id"]))
    plan13_ready = by_id.get(str(signature["plan13_terminal_ready_input_id"]))
    if allow_fixture and not all((closure, gen_ready, plan13_ready)):
        raise contract.PassportContractError("FIXTURE_TERMINAL_REFERENCE", "fixture must exercise all terminal gates")
    if not all((closure, gen_ready, plan13_ready)):
        raise contract.PassportContractError("TERMINAL_GATE_REFERENCE", "signed GEN and Plan13 gates are required")
    closure_frame = _read_frame(closure.path)
    _require_columns(
        closure_frame,
        {"release_id", "terminal_status", "context_rescue_authorized", "negative_claim_authorized"},
        "GEN terminal closure",
    )
    if len(closure_frame) != 1:
        raise contract.PassportContractError("GEN_TERMINAL_CARDINALITY", str(len(closure_frame)))
    closure_row = closure_frame.iloc[0]
    if (
        _clean(closure_row["terminal_status"]) != "coverage_limited_terminal"
        or _bool(closure_row["context_rescue_authorized"], "context_rescue_authorized")
        or _bool(closure_row["negative_claim_authorized"], "negative_claim_authorized")
    ):
        raise contract.PassportContractError("GEN_TERMINAL_STALE", "unexpected GEN closure semantics")
    ready_frame = _read_frame(gen_ready.path)
    _require_columns(ready_frame, {"status", "closure_sha256", "context_rescue_authorized", "negative_claim_authorized"}, "GEN terminal READY")
    if len(ready_frame) != 1:
        raise contract.PassportContractError("GEN_READY_CARDINALITY", str(len(ready_frame)))
    ready = ready_frame.iloc[0]
    if (
        _clean(ready["status"]) != "coverage_limited_terminal_validated"
        or _clean(ready["closure_sha256"]) != closure.sha256
        or _bool(ready["context_rescue_authorized"], "context_rescue_authorized")
        or _bool(ready["negative_claim_authorized"], "negative_claim_authorized")
    ):
        raise contract.PassportContractError("GEN_TERMINAL_STALE", "GEN READY does not attest selected closure")
    semantic_outer = _read_frame(plan13_ready.path)
    if not allow_fixture:
        _require_columns(
            semantic_outer,
            {
                "release_id", "semantic_contract_id", "status", "final_ready_sha256",
                "n_robust", "n_indeterminate", "n_tested_negative",
                "historical_bundles_unchanged", "canonical_promotion_authorized",
            },
            "Plan13 semantic-v2 terminal",
        )
        if (
            len(semantic_outer) != 1
            or _clean(semantic_outer.iloc[0]["semantic_contract_id"])
            != "evidence-state-v2-2026-08-08"
            or _clean(semantic_outer.iloc[0]["status"])
            != "sealed_semantic_v2_candidate_no_canonical_promotion"
            or _optional_int(
                semantic_outer.iloc[0]["n_tested_negative"], "n_tested_negative"
            ) != 0
            or not _bool(
                semantic_outer.iloc[0]["historical_bundles_unchanged"],
                "historical_bundles_unchanged",
            )
            or _bool(
                semantic_outer.iloc[0]["canonical_promotion_authorized"],
                "canonical_promotion_authorized",
            )
        ):
            raise contract.PassportContractError(
                "PLAN13_SEMANTIC_TERMINAL_STALE",
                "selected Plan13 terminal is not the sealed semantic-v2 zero-negative candidate",
            )
        final_ready_inputs = [
            item for item in selected if item.row["adapter_id"] == "plan13_final_ready_v1"
        ]
        if len(final_ready_inputs) != 1:
            raise contract.PassportContractError(
                "PLAN13_FINAL_READY_CARDINALITY", str(len(final_ready_inputs))
            )
        final_ready = final_ready_inputs[0]
        if _clean(semantic_outer.iloc[0]["final_ready_sha256"]) != final_ready.sha256:
            raise contract.PassportContractError(
                "PLAN13_SEMANTIC_TERMINAL_STALE",
                "semantic-v2 outer seal does not hash the selected final READY",
            )
        for item in selected:
            if item.row["adapter_id"] in {"plan13_program_context_v1", "dataset_status_v1"}:
                if item.row["terminal_gate_input_id"] != final_ready.input_id:
                    raise contract.PassportContractError(
                        "PLAN13_FINAL_READY_LINK",
                        f"{item.input_id} must link to {final_ready.input_id}",
                    )
        p13 = _read_frame(final_ready.path)
    else:
        # The synthetic fixture uses one compact terminal row while exercising
        # the same complete/hash checks; this branch cannot unlock production.
        p13 = semantic_outer
    _require_columns(p13, {"release_id", "status", "mode", "plan13_complete"}, "Plan13 final READY")
    if len(p13) != 1 or not _bool(p13.iloc[0]["plan13_complete"], "plan13_complete"):
        raise contract.PassportContractError("PLAN13_INCOMPLETE", "plan13_complete must be TRUE")
    if not allow_fixture and _clean(p13.iloc[0]["mode"]).lower() != "real":
        raise contract.PassportContractError("PLAN13_FIXTURE_OR_PARTIAL", _clean(p13.iloc[0]["mode"]))
    if "fixture" in _clean(p13.iloc[0]["status"]).lower() and not allow_fixture:
        raise contract.PassportContractError("PLAN13_FIXTURE_OR_PARTIAL", _clean(p13.iloc[0]["status"]))
    p13_row = p13.iloc[0]
    for item in selected:
        if item.row["adapter_id"] == "plan13_program_context_v1":
            if "integrated_effects_sha256" not in p13 or _clean(p13_row["integrated_effects_sha256"]) != item.sha256:
                raise contract.PassportContractError(
                    "PLAN13_TERMINAL_STALE", "terminal READY does not hash the selected integrated effects"
                )
        if item.row["adapter_id"] == "dataset_status_v1":
            if "figure4_verdict_sha256" not in p13 or _clean(p13_row["figure4_verdict_sha256"]) != item.sha256:
                raise contract.PassportContractError(
                    "PLAN13_TERMINAL_STALE", "terminal READY does not hash the selected dataset verdict"
                )
    hotspot_inputs = [
        item for item in selected
        if item.row["adapter_id"] in {"hotspot_program_semantics_v1", "hotspot_program_membership_v1"}
    ]
    if hotspot_inputs:
        gates = [item for item in selected if item.row["adapter_id"] == "hotspot_semantic_ready_v1"]
        if len(gates) != 1:
            raise contract.PassportContractError("HOTSPOT_TERMINAL_CARDINALITY", str(len(gates)))
        hotspot_ready = _read_frame(gates[0].path)
        _require_columns(
            hotspot_ready,
            {"status", "adjudication_sha256", "n_tested_negative_authorized"},
            "Hotspot semantic READY",
        )
        if (
            len(hotspot_ready) != 1
            or _clean(hotspot_ready.iloc[0]["status"]) != "semantic_adjudication_validated"
            or _optional_int(
                hotspot_ready.iloc[0]["n_tested_negative_authorized"],
                "n_tested_negative_authorized",
            ) != 0
        ):
            raise contract.PassportContractError("HOTSPOT_TERMINAL_STALE", "semantic READY did not pass")
        semantic_inputs = [
            item for item in selected if item.row["adapter_id"] == "hotspot_program_semantics_v1"
        ]
        if len(semantic_inputs) != 1 or _clean(hotspot_ready.iloc[0]["adjudication_sha256"]) != semantic_inputs[0].sha256:
            raise contract.PassportContractError(
                "HOTSPOT_TERMINAL_STALE", "semantic READY does not hash selected adjudication"
            )
    myojin_inputs = [
        item for item in selected
        if item.row["adapter_id"] in {
            "myojin_assay_verdict_v1", "myojin_class_nonsupport_v1",
            "myojin_program_nonsupport_v1",
        }
    ]
    if myojin_inputs:
        gates = [item for item in selected if item.row["adapter_id"] == "myojin_terminal_ready_v1"]
        if len(gates) != 1:
            raise contract.PassportContractError("MYOJIN_TERMINAL_CARDINALITY", str(len(gates)))
        myojin_ready = _read_frame(gates[0].path)
        _require_columns(myojin_ready, {"status", "external_outcomes_read"}, "Myojin terminal READY")
        if (
            len(myojin_ready) != 1
            or _clean(myojin_ready.iloc[0]["status"]) != "phase_c_release_validated"
            or not _bool(myojin_ready.iloc[0]["external_outcomes_read"], "external_outcomes_read")
        ):
            raise contract.PassportContractError("MYOJIN_TERMINAL_STALE", "Phase C READY did not pass")


def _validate_terminal_links(selected: Sequence[SelectedInput]) -> None:
    by_id = {item.input_id: item for item in selected}
    for item in selected:
        adapter = item.row["adapter_id"]
        terminal_id = item.row["terminal_gate_input_id"]
        required_adapter: str | None = None
        if adapter == "gen_frozen_classes_v1":
            required_adapter = "gen_terminal_ready_v1"
        elif adapter == "gen_identity_adjudication_ready_v1":
            required_adapter = "gen_terminal_ready_v1"
        elif adapter in {
            "gen_adjudicated_classes_v1",
            "gen_identity_adjudication_audit_v1",
            "gen_identity_adjudication_quarantine_v1",
        }:
            required_adapter = "gen_identity_adjudication_ready_v1"
        elif adapter == "gen_raw_identity_source_v1":
            required_adapter = "gen_terminal_ready_v1"
        elif adapter in {"hotspot_program_semantics_v1", "hotspot_program_membership_v1"}:
            required_adapter = "hotspot_semantic_ready_v1"
        elif adapter in {"plan13_program_context_v1", "dataset_status_v1"}:
            if terminal_id not in by_id or by_id[terminal_id].row["adapter_id"] not in {
                "plan13_final_ready_v1", "plan13_terminal_ready_v1"
            }:
                raise contract.PassportContractError(
                    "TERMINAL_LINK_MISMATCH",
                    f"{item.input_id} requires selected Plan13 final READY",
                )
            continue
        elif adapter.startswith("myojin_") and adapter != "myojin_terminal_ready_v1":
            required_adapter = "myojin_terminal_ready_v1"
        if required_adapter:
            if terminal_id not in by_id or by_id[terminal_id].row["adapter_id"] != required_adapter:
                raise contract.PassportContractError(
                    "TERMINAL_LINK_MISMATCH",
                    f"{item.input_id} requires terminal adapter {required_adapter}",
                )


def adapt_selection(
    selected: Sequence[SelectedInput], signature: Mapping[str, Any], allow_fixture: bool
) -> tuple[Assembly, IdentityReference]:
    _validate_gate_artifacts(selected, signature, allow_fixture)
    _validate_terminal_links(selected)
    identity = _load_identity(selected)
    assembly = Assembly(analysis_release_id=str(signature["analysis_release_id"]))
    ordered_adapters = [
        "gen_identity_adjudication_ready_v1",
        "gen_identity_adjudication_quarantine_v1",
        "gen_adjudicated_classes_v1", "gen_frozen_classes_v1",
        "accepted_gene_evidence_v1",
        "synthetic_gene_evidence_v2", "hotspot_program_semantics_v1",
        "hotspot_program_membership_v1", "plan13_program_context_v1",
        "dataset_status_v1", "myojin_assay_verdict_v1",
        "myojin_class_nonsupport_v1", "myojin_program_nonsupport_v1",
    ]
    for adapter in ordered_adapters:
        for item in sorted(
            [value for value in selected if value.row["adapter_id"] == adapter],
            key=lambda value: value.input_id,
        ):
            if adapter == "gen_identity_adjudication_ready_v1":
                adapt_gen_identity_status(assembly, item)
            elif adapter == "gen_identity_adjudication_quarantine_v1":
                adapt_gen_identity_quarantine(assembly, item)
            elif adapter == "gen_adjudicated_classes_v1":
                adapt_gen_classes(assembly, identity, item, adjudicated=True)
            elif adapter == "gen_frozen_classes_v1":
                adapt_gen_classes(assembly, identity, item)
            elif adapter in {"accepted_gene_evidence_v1", "synthetic_gene_evidence_v2"}:
                adapt_accepted_gene_evidence(assembly, identity, item)
            elif adapter == "hotspot_program_semantics_v1":
                adapt_hotspot_semantics(assembly, item)
            elif adapter == "hotspot_program_membership_v1":
                adapt_hotspot_membership(assembly, identity, item)
            elif adapter == "plan13_program_context_v1":
                adapt_plan13_context(assembly, item)
            elif adapter == "dataset_status_v1":
                adapt_dataset_status(assembly, item)
            elif adapter.startswith("myojin_"):
                adapt_myojin(assembly, item)
    symbol_counts = Counter(row["symbol"] for row in assembly.genes.values())
    for row in assembly.genes.values():
        row["symbol_collision"] = symbol_counts[row["symbol"]] > 1
    fatal_quarantine = [
        row for row in assembly.quarantine
        if row["reason_code"] not in {
            "PROGRAM_MEMBERSHIP_UNMAPPED_SOURCE",
            "GEN_IRREDUCIBLE_MULTI_ENSEMBL",
            "GEN_ENSEMBL_ABSENT_V49",
            "GEN_IDENTITY_CONFLICT",
            "GEN_ENSEMBL_UNRESOLVED",
            "GEN_INVALID_ENSEMBL_TOKEN",
        }
    ]
    if fatal_quarantine:
        codes = Counter(row["reason_code"] for row in fatal_quarantine)
        raise contract.PassportContractError(
            "IDENTITY_QUARANTINE_NONEMPTY",
            f"{len(fatal_quarantine)} selected gene-evidence rows require coordinator resolution: {dict(codes)}",
        )
    if not assembly.genes:
        raise contract.PassportContractError("EMPTY_GENE_UNIVERSE", "no accepted gene-grain evidence")
    return assembly, identity


def _source_graph(
    selected: Sequence[SelectedInput], assembly: Assembly
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    nodes: dict[str, dict[str, str]] = {}
    edges: dict[str, dict[str, str]] = {}

    def add_node(row: dict[str, str]) -> None:
        key = row["source_node_id"]
        if key in nodes and nodes[key] != row:
            raise contract.PassportContractError("SOURCE_NODE_CONFLICT", key)
        nodes[key] = row

    def add_edge(source: str, target: str, edge_type: str, result_id: str, note: str) -> None:
        edge_id = _stable_id("edge", source, target, edge_type, result_id)
        edges[edge_id] = {
            "source_edge_id": edge_id, "from_node_id": source, "to_node_id": target,
            "edge_type": edge_type, "evidence_result_id": result_id, "note": note,
        }

    by_input = {item.input_id: item for item in selected}
    for item in selected:
        release_node = _stable_id("node:release", item.input_id, item.row["source_release_id"])
        artifact_uri = (
            f"artifact://selected-input/{item.input_id}/{item.path.name}"
        )
        release_uri = f"release://{item.row['source_release_id']}"
        add_node(
            {
                "source_node_id": item.source_node_id, "node_type": "derived_object",
                "label": item.row["artifact_role"], "source_release_id": item.row["source_release_id"],
                "uri": artifact_uri, "artifact_sha256": item.sha256,
                "access_status": "coordinator_selected_immutable_artifact",
                "biological_unit": item.row["biological_unit"],
                "source_datasets": item.row["source_datasets"],
            }
        )
        add_node(
            {
                "source_node_id": release_node, "node_type": "frozen_release",
                "label": item.row["source_release_id"], "source_release_id": item.row["source_release_id"],
                "uri": release_uri, "artifact_sha256": item.sha256,
                "access_status": "selected_frozen_release", "biological_unit": item.row["biological_unit"],
                "source_datasets": item.row["source_datasets"],
            }
        )
        add_edge(item.source_node_id, release_node, "derived_from", "", "Selected artifact derives from frozen release.")
        dataset_nodes = []
        for dataset in _split_values(item.row["source_datasets"]):
            node_id = _stable_id("node:dataset", item.input_id, dataset, item.row["source_url"])
            dataset_nodes.append(node_id)
            add_node(
                {
                    "source_node_id": node_id, "node_type": "dataset", "label": dataset,
                    "source_release_id": item.row["source_release_id"], "uri": item.row["source_url"],
                    "artifact_sha256": item.sha256, "access_status": "selected_source",
                    "biological_unit": item.row["biological_unit"], "source_datasets": dataset,
                }
            )
            add_edge(release_node, node_id, "derived_from", "", "Frozen release derives from source dataset.")
            for cohort in _split_values(item.row["source_cohorts"]):
                cohort_id = _stable_id("node:cohort", item.input_id, dataset, cohort)
                add_node(
                    {
                        "source_node_id": cohort_id, "node_type": "cohort", "label": cohort,
                        "source_release_id": item.row["source_release_id"], "uri": item.row["source_url"],
                        "artifact_sha256": item.sha256, "access_status": "selected_source",
                        "biological_unit": item.row["biological_unit"], "source_datasets": dataset,
                    }
                )
                add_edge(node_id, cohort_id, "derived_from", "", "Dataset derives from cohort.")
            if item.row["source_publication"]:
                publication_id = _stable_id(
                    "node:publication", item.input_id, dataset,
                    item.row["source_publication"],
                )
                add_node(
                    {
                        "source_node_id": publication_id, "node_type": "publication",
                        "label": item.row["source_publication"], "source_release_id": item.row["source_release_id"],
                        "uri": item.row["source_url"], "artifact_sha256": item.sha256,
                        "access_status": "selected_source", "biological_unit": "not applicable",
                        "source_datasets": dataset,
                    }
                )
                add_edge(node_id, publication_id, "derived_from", "", "Dataset is described by publication.")
    for item in selected:
        for parent_id in _split_values(item.row["parent_input_ids"]):
            parent = by_input[parent_id]
            if parent.input_id == item.input_id:
                raise contract.PassportContractError(
                    "SOURCE_GRAPH_SELF_PARENT", item.input_id
                )
            add_edge(
                item.source_node_id,
                parent.source_node_id,
                "derived_from",
                "",
                "Coordinator-attested selected-artifact dependency.",
            )
    for result in assembly.evidence:
        item = by_input[result["source_input_row_id"].split(":", 1)[0]]
        for dataset in _split_values(item.row["source_datasets"]):
            add_edge(
                result["source_node_id"], _stable_id("node:dataset", item.input_id, dataset, item.row["source_url"]),
                "tested_by", result["evidence_result_id"], "Gene result was tested by selected dataset.",
            )
    for result in assembly.program_context:
        item = by_input[result["source_input_row_id"].split(":", 1)[0]]
        for dataset in _split_values(item.row["source_datasets"]):
            add_edge(
                result["source_node_id"],
                _stable_id("node:dataset", item.input_id, dataset, item.row["source_url"]),
                "tested_by", result["program_context_id"],
                "Program-grain context result was tested by selected dataset; no member-gene expansion is authorized.",
            )
    return sorted(nodes.values(), key=lambda row: row["source_node_id"]), sorted(edges.values(), key=lambda row: row["source_edge_id"])


def _experiments(assembly: Assembly) -> list[dict[str, str]]:
    rows = []
    for gene in sorted(assembly.genes.values(), key=lambda row: row["passport_id"]):
        fields = contract.EXPERIMENT_RULES[gene["next_experiment_rule_id"]]
        rows.append(
            {
                "analysis_release_id": assembly.analysis_release_id,
                "passport_id": gene["passport_id"], "ensembl_id": gene["ensembl_id"],
                "symbol": gene["symbol"], "experiment_rule_id": gene["next_experiment_rule_id"],
                "role_rule_id": gene["role_rule_id"], **fields, "disclaimer": contract.DISCLAIMER,
            }
        )
    return rows


def _write_json_records(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    payload = []
    for row in rows:
        payload.append(
            {
                key: None if pd.isna(value) else value.item() if hasattr(value, "item") else value
                for key, value in row.items()
            }
        )
    contract.atomic_write_json(path, {"records": payload})


def _program_context_domain(assay: str) -> str:
    value = assay.lower()
    if any(token in value for token in ("visium", "geomx", "cosmx", "spatial", "lipid")):
        return "spatial"
    if any(token in value for token in ("protein", "proteom", "dia_ms", "olink")):
        return "proteomics"
    if any(token in value for token in ("atac", "chromatin", "accessibility")):
        return "chromatin"
    return "cell_context"


def _domain_coverage(
    selected: Sequence[SelectedInput], assembly: Assembly
) -> list[dict[str, Any]]:
    """Expose native-grain coverage without projecting program calls to genes."""

    by_input = {item.input_id: item.row["adapter_id"] for item in selected}
    gene_counts: Counter[str] = Counter()
    context_counts: Counter[str] = Counter()
    adapter_ids: dict[str, set[str]] = {
        domain: set() for domain in contract.EVIDENCE_DOMAINS
    }
    for row in assembly.evidence:
        domain = str(row["evidence_domain"])
        gene_counts[domain] += 1
        source_id = str(row["source_input_row_id"]).split(":", 1)[0]
        if source_id in by_input:
            adapter_ids[domain].add(by_input[source_id])
    for row in assembly.program_context:
        domain = _program_context_domain(str(row["assay"]))
        context_counts[domain] += 1
        source_id = str(row["source_input_row_id"]).split(":", 1)[0]
        if source_id in by_input:
            adapter_ids[domain].add(by_input[source_id])
    program_counts: Counter[str] = Counter()
    membership_counts: Counter[str] = Counter()
    if assembly.programs:
        program_counts["cell_context"] = len(assembly.programs)
        adapter_ids["cell_context"].update(
            item.row["adapter_id"]
            for item in selected
            if item.row["adapter_id"] == "hotspot_program_semantics_v1"
        )
    if assembly.memberships:
        membership_counts["cell_context"] = len(assembly.memberships)
        adapter_ids["cell_context"].update(
            item.row["adapter_id"]
            for item in selected
            if item.row["adapter_id"] == "hotspot_program_membership_v1"
        )
    assay_counts: Counter[str] = Counter()
    if assembly.assay_status:
        assay_counts["functional"] = len(assembly.assay_status)
        adapter_ids["functional"].update(
            item.row["adapter_id"]
            for item in selected
            if item.row["adapter_id"].startswith("myojin_")
            and item.row["adapter_id"] != "myojin_terminal_ready_v1"
        )

    rows: list[dict[str, Any]] = []
    for domain in contract.EVIDENCE_DOMAINS:
        n_gene = gene_counts[domain]
        n_program = program_counts[domain]
        n_membership = membership_counts[domain]
        n_context = context_counts[domain]
        n_assay = assay_counts[domain]
        if n_gene:
            boundary = "gene_level"
            wording = f"Accepted {domain} evidence is populated at gene grain."
        elif n_context:
            boundary = "program_context_only"
            wording = (
                f"Accepted {domain} evidence is displayed only at program grain; "
                "member-gene evidence calls are not authorized."
            )
        elif n_program or n_membership:
            boundary = "program_definition_only"
            wording = (
                "Cell-program definitions and membership are displayed for context only; "
                "membership is not a member-gene evidence call."
            )
        elif n_assay:
            boundary = "assay_or_program_status_only"
            wording = (
                "Functional non-support is displayed only at its frozen assay, class, "
                "or program grain; no member-gene result is created."
            )
        else:
            boundary = "not_populated"
            wording = f"No accepted {domain} evidence is populated in this candidate."
        rows.append(
            {
                "analysis_release_id": assembly.analysis_release_id,
                "evidence_domain": domain,
                "gene_evidence_status": "populated" if n_gene else "not_populated",
                "program_context_status": (
                    "populated" if n_context or n_program or n_membership else "not_populated"
                ),
                "assay_status": "populated" if n_assay else "not_populated",
                "gene_evidence_authorized": str(bool(n_gene)).lower(),
                "member_gene_expansion_authorized": "false",
                "n_gene_evidence_rows": n_gene,
                "n_programs": n_program,
                "n_program_membership_rows": n_membership,
                "n_program_context_rows": n_context,
                "n_assay_status_rows": n_assay,
                "source_adapter_ids": ";".join(sorted(adapter_ids[domain])),
                "coverage_boundary": boundary,
                "allowed_wording": wording,
                "limitation": (
                    "Protein, chromatin, spatial, functional, cell-context, and "
                    "cross-species program or assay products never expand to gene calls; "
                    "live-atlas columns are outside the signed selection."
                ),
            }
        )
    return rows


def _review_environment(output_root: Path) -> str:
    manifest = output_root / "portal_candidate/review/review_manifest.tsv"
    if not manifest.is_file() or manifest.is_symlink():
        raise contract.PassportContractError(
            "MANIFEST_RENDER_ENVIRONMENT", str(manifest)
        )
    rows = contract.read_tsv_rows(manifest)
    if not rows:
        raise contract.PassportContractError(
            "MANIFEST_RENDER_ENVIRONMENT", "empty review manifest"
        )
    fields = (
        "browser_label", "browser_version", "browser_sha256", "mesa_module"
    )
    observed = {tuple(row[field] for field in fields) for row in rows}
    if len(observed) != 1:
        raise contract.PassportContractError(
            "MANIFEST_RENDER_ENVIRONMENT", "inconsistent renderer provenance"
        )
    browser_label, browser_version, browser_sha, mesa_module = next(iter(observed))
    return (
        f"renderer:{browser_label};version={browser_version};"
        f"browser_sha256={browser_sha};mesa={mesa_module};"
        "fonts=system_fonts_unfrozen"
    )


def manifest_producer_binding(
    relative_path: str, *, finalized: bool
) -> tuple[str, Path]:
    """Map an artifact to a stable logical producer and executing file."""

    if relative_path in {
        "portal_candidate/index.html",
        "portal_candidate/ui_contract.json",
        "portal_candidate/ui_source_manifest.tsv",
    }:
        return UI_GENERATOR_LOGICAL_PRODUCER_ID, UI_GENERATOR_PRODUCER_PATH
    if relative_path in {
        "portal_candidate/review/overview.png",
        "portal_candidate/review/THRB.png",
        "portal_candidate/review/HKDC1.png",
        "portal_candidate/review/GLP1R.png",
        "portal_candidate/review/MTARC1.png",
        "portal_candidate/review/review_manifest.tsv",
    }:
        return UI_RENDERER_LOGICAL_PRODUCER_ID, UI_RENDERER_PRODUCER_PATH
    if relative_path in {
        "preflight/gen_identity_adjudication_v1/GEN_IDENTITY_ADJUDICATION_READY",
        "preflight/gen_identity_adjudication_v1/gen_frozen_classes_ensembl_adjudicated.tsv",
        "preflight/gen_identity_adjudication_v1/gen_identity_adjudication_audit.tsv",
        "preflight/gen_identity_adjudication_v1/gen_identity_quarantine.tsv",
    }:
        return GEN_IDENTITY_LOGICAL_PRODUCER_ID, GEN_IDENTITY_PRODUCER_PATH
    if relative_path == "passport_input_selection.tsv":
        return SELECTION_LOGICAL_PRODUCER_ID, SELECTION_PRODUCER_PATH
    if relative_path in {
        "passport_input_selection.signature.json",
        "passport_manual_acceptance.tsv",
        "passport_manual_acceptance.signature.json",
    }:
        return ATTESTATION_LOGICAL_PRODUCER_ID, ATTESTATION_PRODUCER_PATH
    if finalized and relative_path in {
        "passport_gate_status.tsv", "passport_build_status.tsv"
    }:
        return VALIDATOR_LOGICAL_PRODUCER_ID, VALIDATOR_PRODUCER_PATH
    return BUNDLE_LOGICAL_PRODUCER_ID, SCRIPT_PATH


def manifest_producer_paths() -> dict[str, Path]:
    """Return the exact logical producer allowlist for the executing tool copy."""

    return {
        BUNDLE_LOGICAL_PRODUCER_ID: SCRIPT_PATH,
        SELECTION_LOGICAL_PRODUCER_ID: SELECTION_PRODUCER_PATH,
        ATTESTATION_LOGICAL_PRODUCER_ID: ATTESTATION_PRODUCER_PATH,
        VALIDATOR_LOGICAL_PRODUCER_ID: VALIDATOR_PRODUCER_PATH,
        GEN_IDENTITY_LOGICAL_PRODUCER_ID: GEN_IDENTITY_PRODUCER_PATH,
        UI_GENERATOR_LOGICAL_PRODUCER_ID: UI_GENERATOR_PRODUCER_PATH,
        UI_RENDERER_LOGICAL_PRODUCER_ID: UI_RENDERER_PRODUCER_PATH,
    }


def _write_manifest(
    output_root: Path, selection_sha256: str, release_status: str,
    *, finalized: bool = False,
) -> Path:
    manifest_path = output_root / "passport_release_manifest.tsv"
    temp_residue = sorted(
        path.relative_to(output_root).as_posix()
        for path in output_root.rglob("*")
        if path.is_file()
        and path.name.startswith(".")
        and path.name.endswith(".tmp")
        and "fixtures" not in path.relative_to(output_root).parts
        and "superseded_bundles" not in path.relative_to(output_root).parts
    )
    if temp_residue:
        raise contract.PassportContractError(
            "MANIFEST_TEMP_RESIDUE", ", ".join(temp_residue)
        )
    excluded = {
        manifest_path.resolve(),
        (output_root / "passport_validation_report.tsv").resolve(),
        (output_root / "passport_plan60_handoff.tsv").resolve(),
        (output_root / "passport_terminal_provenance.tsv").resolve(),
        (output_root / "PASS06_VALIDATED").resolve(),
        (output_root / "FIXTURE_PASS06_VALIDATED").resolve(),
        (output_root / "plan60_terminal_artifacts.tsv").resolve(),
        (output_root / "plan60_terminal_artifacts.signature.json").resolve(),
    }

    rows = []
    review_environment = _review_environment(output_root)
    for path in sorted(value for value in output_root.rglob("*") if value.is_file()):
        relative_parts = path.relative_to(output_root).parts
        if (
            path.resolve() in excluded
            or "fixtures" in relative_parts
            or "superseded_bundles" in relative_parts
        ):
            continue
        relative_path = path.relative_to(output_root).as_posix()
        producer_id, producer_path = manifest_producer_binding(
            relative_path, finalized=finalized
        )
        if not producer_path.is_file() or producer_path.is_symlink():
            raise contract.PassportContractError(
                "MANIFEST_PRODUCER_MISSING", str(producer_path)
            )
        rows.append(
            {
                "relative_path": relative_path,
                "sha256": contract.sha256_file(path), "bytes": path.stat().st_size,
                "producer": producer_id,
                "environment": (
                    review_environment
                    if relative_path.startswith("portal_candidate/review/")
                    else f"spatial:pandas-{pd.__version__}:pyarrow-{pa.__version__}"
                ),
                "producer_sha256": contract.sha256_file(producer_path),
                "upstream_artifact_sha256": selection_sha256,
                "release_status": release_status,
            }
        )
    contract.atomic_write_tsv(manifest_path, rows, contract.TSV_SCHEMAS[manifest_path.name])
    return manifest_path


def write_bundle(
    output_root: Path,
    selection_path: Path,
    selected: Sequence[SelectedInput],
    signature: Mapping[str, Any],
    assembly: Assembly,
    fixture_mode: bool,
) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    selection_sha = contract.sha256_file(selection_path)
    selection_destination = output_root / "passport_input_selection.tsv"
    signature_destination = output_root / "passport_input_selection.signature.json"
    if selection_path.resolve() != selection_destination.resolve():
        shutil.copy2(selection_path, selection_destination)
        shutil.copy2(selection_path.with_name("passport_input_selection.signature.json"), signature_destination)
    nodes, edges = _source_graph(selected, assembly)
    genes = sorted(assembly.genes.values(), key=lambda row: row["passport_id"])
    evidence = sorted(assembly.evidence, key=lambda row: (row["passport_id"], row["evidence_result_id"]))
    coverage = sorted(assembly.coverage, key=lambda row: (row["passport_id"], row["coverage_result_id"]))
    programs = sorted(assembly.programs.values(), key=lambda row: row["program_uid"])
    memberships = sorted(assembly.memberships, key=lambda row: (row["program_uid"], row["ensembl_id"], row["source_input_row_id"]))
    contexts = sorted(assembly.program_context, key=lambda row: (row["program_uid"], row["dataset_id"], row["program_context_id"]))
    for filename, rows in (
        ("passport_gene_index.parquet", genes),
        ("passport_evidence_long.parquet", evidence),
        ("passport_coverage_long.parquet", coverage),
        ("passport_program_index.parquet", programs),
        ("passport_program_membership.parquet", memberships),
        ("passport_program_context.parquet", contexts),
    ):
        schema = contract.PRODUCTION_PARQUET_SCHEMAS[filename]
        contract.atomic_write_parquet(output_root / filename, contract.coerce_dataframe(rows, schema), schema)
    contract.atomic_write_tsv(output_root / "passport_source_nodes.tsv", nodes, contract.TSV_SCHEMAS["passport_source_nodes.tsv"])
    contract.atomic_write_tsv(output_root / "passport_source_edges.tsv", edges, contract.TSV_SCHEMAS["passport_source_edges.tsv"])
    contract.atomic_write_tsv(output_root / "passport_next_experiment.tsv", _experiments(assembly), contract.TSV_SCHEMAS["passport_next_experiment.tsv"])
    contract.atomic_write_tsv(output_root / "passport_rulebook.tsv", contract.rulebook_rows(), contract.TSV_SCHEMAS["passport_rulebook.tsv"])
    contract.atomic_write_tsv(output_root / "passport_controlled_vocabularies.tsv", contract.controlled_vocabulary_rows(), contract.TSV_SCHEMAS["passport_controlled_vocabularies.tsv"])
    contract.atomic_write_tsv(output_root / "passport_data_dictionary.tsv", contract.data_dictionary_rows(include_production=True), contract.TSV_SCHEMAS["passport_data_dictionary.tsv"])
    accepted_audit = []
    for item in selected:
        accepted_audit.append(
            {
                **{column: item.row[column] for column in contract.PRODUCTION_ONLY_TSV_SCHEMAS["accepted_input_audit.tsv"] if column in item.row},
                "selection_sha256": selection_sha, "validation_status": "pass",
                "validation_reason": "signed selection, artifact checksum, producer checksum, and terminal links validated",
            }
        )
    contract.atomic_write_tsv(output_root / "accepted_input_audit.tsv", accepted_audit, contract.PRODUCTION_ONLY_TSV_SCHEMAS["accepted_input_audit.tsv"])
    contract.atomic_write_tsv(output_root / "passport_adapter_audit.tsv", assembly.adapter_audit, contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_adapter_audit.tsv"])
    contract.atomic_write_tsv(output_root / "passport_identity_quarantine.tsv", assembly.quarantine, contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_identity_quarantine.tsv"])
    contract.atomic_write_tsv(output_root / "passport_gen_identity_status.tsv", assembly.gen_identity_status, contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_gen_identity_status.tsv"])
    domain_coverage = _domain_coverage(selected, assembly)
    contract.atomic_write_tsv(output_root / "passport_domain_coverage.tsv", domain_coverage, contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_domain_coverage.tsv"])
    contract.atomic_write_tsv(output_root / "passport_dataset_status.tsv", assembly.dataset_status, contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_dataset_status.tsv"])
    contract.atomic_write_tsv(output_root / "passport_assay_status.tsv", assembly.assay_status, contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_assay_status.tsv"])
    gate_rows = [
        {"gate_id": "PASS00_SIGNED_SELECTION", "status": "pass", "promotion_allowed": "false", "reason": "Coordinator selection and immutable inputs validated."},
        {"gate_id": "PASS02_SOURCE_GRAPH", "status": "pass", "promotion_allowed": "false", "reason": "Source nodes, edges, and provenance emitted."},
        {"gate_id": "PASS03_EVIDENCE_TABLES", "status": "pass", "promotion_allowed": "false", "reason": "Gene, evidence, coverage, and program-grain tables assembled without call recomputation."},
        {"gate_id": "PASS04_EXPERIMENT_ROUTING", "status": "pass", "promotion_allowed": "false", "reason": "Deterministic rulebook generated one discriminating experiment per gene."},
        {"gate_id": "PASS05_PORTAL_EXPORT", "status": "pass", "promotion_allowed": "false", "reason": "Candidate exports and local review UI generated; PASS06 validation remains required."},
        {"gate_id": "PASS06_TERMINAL", "status": "pending_validation", "promotion_allowed": "false", "reason": "Automated and signed manual acceptance must complete."},
    ]
    contract.atomic_write_tsv(output_root / "passport_gate_status.tsv", gate_rows, contract.TSV_SCHEMAS["passport_gate_status.tsv"])
    build_rows = [
        {"stage_id": f"PASS0{stage}", "status": "pass", "promotion_allowed": "false", "selection_sha256": selection_sha, "reason": "candidate assembly complete"}
        for stage in range(2, 6)
    ] + [
        {"stage_id": "PASS06", "status": "pending_validation", "promotion_allowed": "false", "selection_sha256": selection_sha, "reason": "terminal validation not yet run"}
    ]
    contract.atomic_write_tsv(output_root / "passport_build_status.tsv", build_rows, contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_build_status.tsv"])
    portal_root = output_root / "portal_export"
    portal_root.mkdir(parents=True, exist_ok=True)
    summary_schema = contract.PORTAL_PARQUET_SCHEMAS["portal_export/passport_summary.parquet"]
    summary_rows = [{**row, "alphabetical_key": f"{row['symbol'].upper()}|{row['ensembl_id']}"} for row in genes]
    contract.atomic_write_parquet(portal_root / "passport_summary.parquet", contract.coerce_dataframe(summary_rows, summary_schema), summary_schema)
    for source_name, destination_name in (
        ("passport_evidence_long.parquet", "passport_evidence.parquet"),
        ("passport_program_index.parquet", "passport_program_index.parquet"),
        ("passport_program_membership.parquet", "passport_program_membership.parquet"),
        ("passport_program_context.parquet", "passport_program_context.parquet"),
    ):
        shutil.copy2(output_root / source_name, portal_root / destination_name)
    dictionary_rows = contract.data_dictionary_rows(include_production=True)
    _write_json_records(portal_root / "passport_dictionary.json", dictionary_rows)
    contract.atomic_write_json(portal_root / "passport_source_graph.json", {"nodes": nodes, "edges": edges})
    _write_json_records(portal_root / "passport_domain_coverage.json", domain_coverage)
    contract.atomic_write_json(
        portal_root / "passport_bundle_status.json",
        {
            "analysis_release_id": assembly.analysis_release_id,
            "selection_sha256": selection_sha,
            "fixture_only": fixture_mode,
            "default_order": "alphabetical",
            "combined_score_constructed": False,
            "scientific_call_recomputed": False,
            "gene_level_domain_scope": ["genetics", "transcriptomics"],
            "program_or_assay_only_domains": [
                "cell_context", "spatial", "proteomics", "chromatin", "functional"
            ],
            "cross_species_populated": False,
            "gene_rows": len(genes), "evidence_rows": len(evidence),
            "program_rows": len(programs), "program_context_rows": len(contexts),
        },
    )
    try:
        passport_ui.build_ui(output_root, fixture_mode=fixture_mode)
        passport_ui_renderer.render_ui(
            output_root,
            mesa_module=os.environ.get(
                "PASSPORT_MESA_MODULE", passport_ui_renderer.EXPECTED_MESA_MODULE
            ),
        )
    except (passport_ui.PassportUIError, passport_ui_renderer.PassportUIRenderError) as exc:
        raise contract.PassportContractError("PASS05_UI_BUILD", str(exc)) from exc
    _write_manifest(output_root, selection_sha, "synthetic_fixture_not_for_promotion" if fixture_mode else "assembled_candidate_prevalidation")
    return output_root


def build_from_selection(
    selection_path: Path, output_root: Path, fixture_mode: bool = False
) -> Path:
    selected, signature = _resolve_selected_inputs(selection_path, allow_fixture=fixture_mode)
    assembly, _ = adapt_selection(selected, signature, allow_fixture=fixture_mode)
    return write_bundle(output_root, selection_path, selected, signature, assembly, fixture_mode)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, default=contract.DEFAULT_CANDIDATE_ROOT / "passport_input_selection.tsv")
    parser.add_argument("--output-root", type=Path, default=contract.DEFAULT_CANDIDATE_ROOT)
    parser.add_argument("--fixture-mode", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output_root.resolve()
    candidate = contract.DEFAULT_CANDIDATE_ROOT.resolve()
    if args.fixture_mode:
        contract.assert_fixture_path(output, candidate)
    elif output != candidate:
        raise contract.PassportContractError("OUTPUT_ROOT", f"production output must be {candidate}")
    result = build_from_selection(args.selection.resolve(), output, fixture_mode=args.fixture_mode)
    print(f"PASS: PASS-02--PASS-05 candidate assembled at {result}")


if __name__ == "__main__":
    try:
        main()
    except contract.PassportContractError as exc:
        raise SystemExit(str(exc)) from exc
