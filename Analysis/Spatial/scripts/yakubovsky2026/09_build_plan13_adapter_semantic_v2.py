#!/usr/bin/env python3
"""Build the semantic-v2 Plan 13 view of immutable Yakubovsky outputs.

This script is deliberately outside the frozen Plan 11 analysis contract.  It
does not refit a model, rescore a spot, or rewrite the historical Plan 13
adapter.  Instead, it verifies the byte-identical historical bundle, carries
its assay-native estimates forward, and corrects evidence-state semantics in a
new write-once candidate directory.
"""

from __future__ import annotations

import argparse
import math
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT_DEFAULT = SCRIPT_PATH.parents[4]
CONTRACT_ROOT_RELATIVE = Path(
    "Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
)
sys.path.insert(0, str(PROJECT_ROOT_DEFAULT / CONTRACT_ROOT_RELATIVE))

from contract_lib import (  # noqa: E402
    ADAPTER_REGISTRY_COLUMNS,
    ContractError,
    RELEASE_ID,
    SCHEMAS,
    SEMANTIC_CONTRACT_ID,
    parse_bool,
    read_tsv,
    sha256_file,
    validate_candidate_contract,
    write_tsv,
)


DATASET = "Yakubovsky_2026"
ADAPTER_ID = "yakubovsky2026_binary_lipid"
FROZEN_ADAPTER_PRODUCER_RELATIVE = Path(
    "Analysis/Spatial/scripts/yakubovsky2026/05_build_plan13_adapter.py"
)
FROZEN_ADAPTER_PRODUCER_SHA256 = (
    "f66c305dd2afbf25dce4975ac3a17a09b711d092f035ddf6814ca2f677728c1d"
)

# These hashes identify the reviewed, immutable pre-remediation adapter.  A
# semantic rebuild fails closed if any historical byte changes.
LEGACY_ADAPTER_SHA256 = {
    "PLAN13_ADAPTER_READY": "81622102e392641ddca9eae04fd02ec04e95d65567fc46daccff3e7eabd3c036",
    "adapter_registry.tsv": "491c9265dd4799b4ea5c4055c5e0695fed234121cdcd1d0607f1cf9554ea0b05",
    "design_audit.tsv": "69564d8ce031c113b394dc8ca1fe86331467b1c55d1a2d5ad9782505039b2b2a",
    "execution_manifest.tsv": "e9be729c5bcd1fe34353a3cf20708d843b87f6aab9e0b9e9fa511f2791278410",
    "gate_status.tsv": "b26efc0f86ae5f52ea1033a658be987fb94015498ec96714b468f44e50f2e0df",
    "gene_mapping_audit.tsv": "beb956ca82bc8d060f73e3218d7258a4c1136908ab6e1ec3b8b4eb01c78d950a",
    "per_sample_program_scores.tsv": "4020ddcea0910c7a1d2275dadf7a62e3f37a9594bdc66b7ee92907d096890e1c",
    "program_effects.tsv": "5ca3bcdd6529fa1ffaee26adeedc3bd16ed576cce3c10cc4f00e9b9770ab24da",
    "program_testability.tsv": "828fc3769ce2b34571c800f821ac402b38354bd3ce9ac4ddc7cda356d1d4a02d",
    "sample_manifest.tsv": "52d024c6133b89b1fbf90d5508588af904636e0be8a5eff02131824acaa074dc",
    "sensitivity.tsv": "e61a5e75109c8bf9145f706d82d382b716f3d73752310d27e8fde9be4889bbd7",
    "source_manifest.tsv": "3da0437520dcb1208f2eb79d5c2dcbbb084e4f05b83351fd96d8f5ceb40e05c9",
    "validation_report.tsv": "5dc4704da4477744bdb6cebba3e3a55bfb95fa7f17954c5719f99ceb8ca9d6c5",
}

NATIVE_INPUT_SHA256 = {
    "READY": "3dbf18fb10db55717b816e9cab30480ec1007e71928dec9d6d70a08e1acb1716",
    "analysis_freeze_manifest.tsv": "1ddee881942a4e14ba98ecd24639d22289cefefbce9cfcf54448eef984d283b5",
    "donor_program_effects.tsv": "f9b6c825d8c6e1416467954e43ff9c011b7cf51da76a70c1d5477125295d4d79",
    "gate_status.tsv": "39e9539b975797b6e333c7a71c1d273674cf46d135b13238f9ca501ba4070424",
    "program_effects.tsv": "20ff894781dd78609dadd3743569fd6f3a0ec999d1d7fcd429f3cb35e4fb1231",
    "program_testability.tsv": "3341489e7c671446e02828c026903ce644fa39e52fbf4b231ce2d0451ef70324",
    "sensitivity.tsv": "a7ed790d72a9dc6ef34fb063fef76a2462bc52515b327d9fc9fcad310ecd215e",
}

LEGACY_READY_CONTRACT_LIB_SHA256 = (
    "93cbab090e3b53906b034d60cccb4608757e42d183e1a850dd7fc0392d6186c1"
)
LEGACY_READY_VALIDATOR_SHA256 = (
    "aab3d5cd0d4f88ef7e377acdd75094048335cc54dc8b7d2d53b7cf430d2490b0"
)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def one_row(path: Path, required: tuple[str, ...] = ()) -> dict[str, str]:
    _, rows = read_tsv(path, required or None)
    if len(rows) != 1:
        raise ContractError(f"expected exactly one row in {path}, found {len(rows)}")
    return rows[0]


def direction_from_number(value: str | float) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"direction input is not numeric: {value!r}") from exc
    if not (-float("inf") < numeric < float("inf")):
        raise ContractError(f"direction input is not finite: {value!r}")
    return "positive" if numeric > 0 else "negative" if numeric < 0 else "zero"


def semantic_evidence_state(
    native_state: str,
    robustness_pass: bool,
    negative_call_rule_id: str = "",
) -> str:
    """Map a historical label to semantic-v2 without inventing adequacy.

    A nonsignificant result is indeterminate unless an explicit, prespecified
    adequate-negative rule exists.  This Yakubovsky bridge supplies no such
    rule, so its historical ``tested_negative`` rows become ``indeterminate``.
    """

    if native_state == "robust":
        if not robustness_pass:
            raise ContractError("historical robust state lacks robustness_pass=TRUE")
        if negative_call_rule_id:
            raise ContractError("robust state cannot carry a negative-call rule")
        return "robust"
    if native_state == "tested_negative":
        if robustness_pass:
            raise ContractError("historical tested_negative has robustness_pass=TRUE")
        return "tested_negative" if negative_call_rule_id.strip() else "indeterminate"
    if native_state == "indeterminate":
        if robustness_pass or negative_call_rule_id:
            raise ContractError("indeterminate state has incompatible robustness/rule")
        return "indeterminate"
    if native_state in {"untestable", "not_applicable", "skipped"}:
        if robustness_pass or negative_call_rule_id:
            raise ContractError("terminal state has incompatible robustness/rule")
        return native_state
    raise ContractError(f"unsupported historical evidence state: {native_state!r}")


def semantic_effect_row(
    legacy: Mapping[str, str],
    native: Mapping[str, str],
    *,
    producer: str,
    producer_sha256: str,
    source_manifest_sha256: str,
) -> dict[str, object]:
    """Add semantic-v2 direction/heterogeneity fields without changing effects."""

    uid = legacy.get("program_uid", "")
    if not uid or native.get("program_uid") != uid:
        raise ContractError(f"legacy/native program mismatch: {uid!r}")
    for field in ("estimate", "pvalue", "qvalue", "membership_sha256"):
        native_field = "qvalue" if field == "qvalue" else field
        legacy_field = "padj" if field == "qvalue" else field
        if field == "membership_sha256":
            if legacy.get(field) != native.get(field):
                raise ContractError(f"legacy/native membership mismatch for {uid}")
        elif not math.isclose(
            float(legacy[legacy_field]),
            float(native[native_field]),
            rel_tol=1e-12,
            abs_tol=1e-15,
        ):
            raise ContractError(f"legacy/native {field} mismatch for {uid}")

    robustness_pass = parse_bool(
        legacy["robustness_pass"], f"legacy.robustness_pass[{uid}]"
    )
    negative_rule = ""
    state = semantic_evidence_state(
        legacy["evidence_state"], robustness_pass, negative_rule
    )
    testable = parse_bool(legacy["testable"], f"legacy.testable[{uid}]")

    row: dict[str, object] = {column: legacy.get(column, "") for column in SCHEMAS["program_effects.tsv"]}
    row.update(
        {
            "matched_null_sd": "",
            "negative_call_rule_id": negative_rule,
            "evidence_state": state,
            "producer": producer,
            "producer_sha256": producer_sha256,
            "source_manifest_sha256": source_manifest_sha256,
        }
    )
    if testable:
        descriptive = direction_from_number(legacy["estimate"])
        inferential = direction_from_number(native["combined_z"])
        native_descriptive = native.get("median_donor_slope_direction", "")
        native_inferential = native.get("combined_direction", "")
        if native_descriptive != descriptive or native_inferential != inferential:
            raise ContractError(f"native direction audit drift for {uid}")
        agreement = descriptive == inferential
        if parse_bool(
            native["combined_and_median_direction_agree"],
            f"native.combined_and_median_direction_agree[{uid}]",
        ) != agreement:
            raise ContractError(f"native direction-agreement flag drift for {uid}")
        if legacy["direction_observed"] != descriptive:
            raise ContractError(f"legacy descriptive direction drift for {uid}")
        row.update(
            {
                "descriptive_effect_direction": descriptive,
                "inferential_test_direction": inferential,
                "direction_agreement": agreement,
                "heterogeneity_statistic": native["cochran_q_descriptive"],
                "heterogeneity_df": native["cochran_q_df"],
                "heterogeneity_pvalue": native["cochran_q_p_descriptive"],
            }
        )
    else:
        row.update(
            {
                "descriptive_effect_direction": "",
                "inferential_test_direction": "",
                "direction_agreement": "",
                "heterogeneity_statistic": "",
                "heterogeneity_df": "",
                "heterogeneity_pvalue": "",
            }
        )
    return row


def _verify_manifest(root: Path, project_root: Path, filename: str) -> None:
    _, rows = read_tsv(root / filename, SCHEMAS[filename])
    for row in rows:
        relative = Path(row["relative_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ContractError(f"unsafe path in historical {filename}: {relative}")
        if filename == "source_manifest.tsv":
            scope = row["path_scope"]
            if scope == "project_relative":
                path = project_root / relative
            elif scope == "adapter_relative":
                path = root / relative
            else:
                raise ContractError(f"invalid historical source scope: {scope!r}")
        else:
            path = root / relative
        if (
            not path.is_file()
            or path.stat().st_size != int(row["bytes"])
            or sha256_file(path) != row["sha256"]
        ):
            raise ContractError(f"historical manifest linkage drift: {path}")


def verify_immutable_inputs(
    project_root: Path, native_root: Path, legacy_root: Path
) -> dict[str, str]:
    frozen_producer = project_root / FROZEN_ADAPTER_PRODUCER_RELATIVE
    if sha256_file(frozen_producer) != FROZEN_ADAPTER_PRODUCER_SHA256:
        raise ContractError("frozen 05_build_plan13_adapter.py hash drift")
    for name, expected in NATIVE_INPUT_SHA256.items():
        path = native_root / name
        if not path.is_file() or sha256_file(path) != expected:
            raise ContractError(f"immutable native Yakubovsky input drift: {path}")
    for name, expected in LEGACY_ADAPTER_SHA256.items():
        path = legacy_root / name
        if not path.is_file() or sha256_file(path) != expected:
            raise ContractError(f"immutable historical Yakubovsky adapter drift: {path}")

    legacy_ready = one_row(legacy_root / "PLAN13_ADAPTER_READY")
    expected_ready = {
        "status": "validated_plan13_adapter_candidate",
        "native_ready_sha256": NATIVE_INPUT_SHA256["READY"],
        "adapter_registry_sha256": LEGACY_ADAPTER_SHA256["adapter_registry.tsv"],
        "gate_status_sha256": LEGACY_ADAPTER_SHA256["gate_status.tsv"],
        "validation_report_sha256": LEGACY_ADAPTER_SHA256["validation_report.tsv"],
        "contract_lib_sha256": LEGACY_READY_CONTRACT_LIB_SHA256,
        "contract_validator_sha256": LEGACY_READY_VALIDATOR_SHA256,
        "canonical_promotion_authorized": "False",
    }
    for field, expected in expected_ready.items():
        if legacy_ready.get(field) != expected:
            raise ContractError(f"historical adapter READY drift for {field}")

    legacy_registry = one_row(legacy_root / "adapter_registry.tsv")
    if (
        legacy_registry.get("adapter_id") != ADAPTER_ID
        or legacy_registry.get("gate_expectation") != "valid_null"
    ):
        raise ContractError("historical adapter registry semantic baseline drift")
    legacy_gate = one_row(legacy_root / "gate_status.tsv")
    if (
        legacy_gate.get("gate_status") != "valid_null"
        or legacy_gate.get("source_manifest_sha256")
        != LEGACY_ADAPTER_SHA256["source_manifest.tsv"]
        or legacy_gate.get("execution_manifest_sha256")
        != LEGACY_ADAPTER_SHA256["execution_manifest.tsv"]
    ):
        raise ContractError("historical adapter gate linkage drift")
    # The historical source manifest and validation report are themselves
    # byte-pinned above.  Do not re-read multi-gigabyte raw sources here: this
    # semantic layer consumes the immutable adapter/native result tables, not
    # the raw expression object.
    _verify_manifest(legacy_root, project_root, "execution_manifest.tsv")

    _, legacy_effects = read_tsv(legacy_root / "program_effects.tsv")
    if len(legacy_effects) != 2:
        raise ContractError("historical Yakubovsky adapter must contain two effects")
    for row in legacy_effects:
        if (
            row.get("producer") != FROZEN_ADAPTER_PRODUCER_RELATIVE.as_posix()
            or row.get("producer_sha256") != FROZEN_ADAPTER_PRODUCER_SHA256
            or row.get("source_manifest_sha256")
            != LEGACY_ADAPTER_SHA256["source_manifest.tsv"]
        ):
            raise ContractError("historical effect provenance drift")
    return legacy_ready


def _source_manifest_rows(
    project_root: Path,
    native_root: Path,
    legacy_root: Path,
    created: str,
) -> list[dict[str, object]]:
    sources: dict[Path, str] = {}
    for name in sorted(LEGACY_ADAPTER_SHA256):
        sources[legacy_root / name] = f"immutable_legacy_plan13_adapter_{name}"
    for name in sorted(NATIVE_INPUT_SHA256):
        sources[native_root / name] = f"immutable_native_plan11_{name}"
    sources[project_root / FROZEN_ADAPTER_PRODUCER_RELATIVE] = (
        "immutable_historical_adapter_producer"
    )
    sources[SCRIPT_PATH] = "semantic_v2_adapter_producer"
    sources[project_root / CONTRACT_ROOT_RELATIVE / "contract_lib.py"] = (
        "semantic_v2_shared_contract"
    )
    sources[project_root / CONTRACT_ROOT_RELATIVE / "03_validate_contract.py"] = (
        "semantic_v2_contract_validator"
    )
    rows = []
    for path, role in sorted(sources.items(), key=lambda item: item[0].as_posix()):
        if not path.is_file() or not path.resolve().is_relative_to(project_root):
            raise ContractError(f"invalid semantic-v2 source: {path}")
        rows.append(
            {
                "path_scope": "project_relative",
                "relative_path": path.resolve().relative_to(project_root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "source_role": role,
                "public_access": False,
                "retrieved_utc": created,
            }
        )
    return rows


def _gate_from_states(states: list[str], outcomes_tested: bool) -> tuple[str, str]:
    if not outcomes_tested:
        return "untestable", "immutable_native_plan11_spatial_model_unestimable"
    if "robust" in states:
        return "pass", "semantic_v2_at_least_one_prespecified_robust_program"
    if states and all(state == "tested_negative" for state in states):
        return "valid_null", "semantic_v2_all_tested_programs_pass_adequate_negative_rules"
    return (
        "indeterminate",
        "semantic_v2_nonrobust_family_without_prespecified_adequate_negative_rule",
    )


def _semantic_rows(
    project_root: Path,
    native_root: Path,
    legacy_root: Path,
    source_manifest_sha256: str,
) -> tuple[dict[str, list[dict[str, object]]], str, str, bool]:
    tables: dict[str, list[dict[str, object]]] = {}
    for filename in (
        "sample_manifest.tsv",
        "gene_mapping_audit.tsv",
        "design_audit.tsv",
        "per_sample_program_scores.tsv",
        "sensitivity.tsv",
    ):
        _, rows = read_tsv(legacy_root / filename)
        tables[filename] = [dict(row) for row in rows]

    _, native_effect_rows = read_tsv(native_root / "program_effects.tsv")
    native_by_uid = {row["program_uid"]: row for row in native_effect_rows}
    if len(native_by_uid) != len(native_effect_rows):
        raise ContractError("duplicate native Yakubovsky program effect")

    producer = SCRIPT_PATH.relative_to(project_root).as_posix()
    producer_sha = sha256_file(SCRIPT_PATH)
    _, legacy_effect_rows = read_tsv(legacy_root / "program_effects.tsv")
    effect_rows = [
        semantic_effect_row(
            row,
            native_by_uid.get(row["program_uid"], {}),
            producer=producer,
            producer_sha256=producer_sha,
            source_manifest_sha256=source_manifest_sha256,
        )
        for row in legacy_effect_rows
    ]
    tables["program_effects.tsv"] = effect_rows

    state_by_uid = {str(row["program_uid"]): str(row["evidence_state"]) for row in effect_rows}
    _, legacy_testability = read_tsv(legacy_root / "program_testability.tsv")
    testability_rows: list[dict[str, object]] = []
    for row in legacy_testability:
        uid = row["program_uid"]
        if uid not in state_by_uid:
            raise ContractError(f"historical testability lacks matched effect: {uid}")
        transformed = dict(row)
        transformed["evidence_state"] = state_by_uid[uid]
        testability_rows.append(transformed)
    tables["program_testability.tsv"] = testability_rows

    legacy_gate = one_row(legacy_root / "gate_status.tsv")
    outcomes_tested = parse_bool(legacy_gate["outcomes_tested"], "legacy.outcomes_tested")
    gate_status, gate_reason = _gate_from_states(list(state_by_uid.values()), outcomes_tested)
    return tables, gate_status, gate_reason, outcomes_tested


def build_semantic_adapter(
    project_root: Path, native_root: Path, output_root: Path
) -> None:
    legacy_root = native_root / "plan13_adapter"
    verify_immutable_inputs(project_root, native_root, legacy_root)
    if output_root.exists() and any(output_root.iterdir()):
        raise ContractError(f"refusing to overwrite semantic-v2 adapter: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)

    contract_root = project_root / CONTRACT_ROOT_RELATIVE
    hotspot_root = (
        project_root
        / "Analysis/Multimodal_Program_Projection/candidates"
        / RELEASE_ID
        / "hotspot"
    )
    created = utc_now()
    source_rows = _source_manifest_rows(
        project_root, native_root, legacy_root, created
    )
    write_tsv(
        output_root / "source_manifest.tsv",
        SCHEMAS["source_manifest.tsv"],
        source_rows,
    )
    source_hash = sha256_file(output_root / "source_manifest.tsv")
    tables, gate_status, gate_reason, outcomes_tested = _semantic_rows(
        project_root, native_root, legacy_root, source_hash
    )
    for filename, rows in tables.items():
        write_tsv(output_root / filename, SCHEMAS[filename], rows)

    execution_files = (
        "sample_manifest.tsv",
        "gene_mapping_audit.tsv",
        "design_audit.tsv",
        "program_testability.tsv",
        "per_sample_program_scores.tsv",
        "program_effects.tsv",
        "sensitivity.tsv",
        "source_manifest.tsv",
    )
    write_tsv(
        output_root / "execution_manifest.tsv",
        SCHEMAS["execution_manifest.tsv"],
        [
            {
                "role": "immutable_plan11_to_semantic_v2_adapter_artifact",
                "relative_path": filename,
                "bytes": (output_root / filename).stat().st_size,
                "sha256": sha256_file(output_root / filename),
                "created_utc": created,
            }
            for filename in execution_files
        ],
    )

    legacy_registry = one_row(legacy_root / "adapter_registry.tsv")
    adapter_row = {column: legacy_registry.get(column, "") for column in ADAPTER_REGISTRY_COLUMNS}
    adapter_row.update(
        {
            "gate_expectation": gate_status,
            "status": "validated_semantic_v2_candidate_adapter",
        }
    )
    write_tsv(output_root / "adapter_registry.tsv", ADAPTER_REGISTRY_COLUMNS, [adapter_row])
    write_tsv(
        output_root / "gate_status.tsv",
        SCHEMAS["gate_status.tsv"],
        [
            {
                "release_id": RELEASE_ID,
                "dataset": DATASET,
                "analysis_set_id": legacy_registry["analysis_set_id"],
                "gate_status": gate_status,
                "gate_reason": gate_reason,
                "outcomes_tested": outcomes_tested,
                "registry_sha256": legacy_registry["registry_sha256"],
                "source_manifest_sha256": source_hash,
                "execution_manifest_sha256": sha256_file(
                    output_root / "execution_manifest.tsv"
                ),
            }
        ],
    )

    checks = validate_candidate_contract(
        output_root,
        output_root / "adapter_registry.tsv",
        hotspot_root,
        project_root,
    )
    contract_sha = sha256_file(contract_root / "contract_lib.py")
    validator_sha = sha256_file(contract_root / "03_validate_contract.py")
    write_tsv(
        output_root / "validation_report.tsv",
        (
            "release_id",
            "semantic_contract_id",
            "check_id",
            "status",
            "detail",
            "contract_lib_sha256",
            "validated_utc",
        ),
        [
            {
                "release_id": RELEASE_ID,
                "semantic_contract_id": SEMANTIC_CONTRACT_ID,
                "check_id": check.check_id,
                "status": check.status,
                "detail": check.detail,
                "contract_lib_sha256": contract_sha,
                "validated_utc": utc_now(),
            }
            for check in checks
        ],
    )
    write_tsv(
        output_root / "PLAN13_ADAPTER_READY",
        (
            "release_id",
            "semantic_contract_id",
            "status",
            "native_ready_sha256",
            "legacy_adapter_ready_sha256",
            "adapter_registry_sha256",
            "gate_status_sha256",
            "source_manifest_sha256",
            "execution_manifest_sha256",
            "validation_report_sha256",
            "contract_lib_sha256",
            "contract_validator_sha256",
            "producer_sha256",
            "canonical_promotion_authorized",
            "completed_utc",
            "python",
            "platform",
            "conda_prefix",
        ),
        [
            {
                "release_id": RELEASE_ID,
                "semantic_contract_id": SEMANTIC_CONTRACT_ID,
                "status": "validated_plan13_adapter_candidate",
                "native_ready_sha256": NATIVE_INPUT_SHA256["READY"],
                "legacy_adapter_ready_sha256": LEGACY_ADAPTER_SHA256[
                    "PLAN13_ADAPTER_READY"
                ],
                "adapter_registry_sha256": sha256_file(
                    output_root / "adapter_registry.tsv"
                ),
                "gate_status_sha256": sha256_file(output_root / "gate_status.tsv"),
                "source_manifest_sha256": source_hash,
                "execution_manifest_sha256": sha256_file(
                    output_root / "execution_manifest.tsv"
                ),
                "validation_report_sha256": sha256_file(
                    output_root / "validation_report.tsv"
                ),
                "contract_lib_sha256": contract_sha,
                "contract_validator_sha256": validator_sha,
                "producer_sha256": sha256_file(SCRIPT_PATH),
                "canonical_promotion_authorized": False,
                "completed_utc": utc_now(),
                "python": sys.version.replace("\n", " "),
                "platform": platform.platform(),
                "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
            }
        ],
    )


def check_semantic_adapter(
    project_root: Path, native_root: Path, output_root: Path
) -> dict[str, object]:
    legacy_root = native_root / "plan13_adapter"
    verify_immutable_inputs(project_root, native_root, legacy_root)
    hotspot_root = (
        project_root
        / "Analysis/Multimodal_Program_Projection/candidates"
        / RELEASE_ID
        / "hotspot"
    )
    checks = validate_candidate_contract(
        output_root,
        output_root / "adapter_registry.tsv",
        hotspot_root,
        project_root,
    )
    ready = one_row(output_root / "PLAN13_ADAPTER_READY")
    contract_root = project_root / CONTRACT_ROOT_RELATIVE
    expected = {
        "release_id": RELEASE_ID,
        "semantic_contract_id": SEMANTIC_CONTRACT_ID,
        "status": "validated_plan13_adapter_candidate",
        "native_ready_sha256": NATIVE_INPUT_SHA256["READY"],
        "legacy_adapter_ready_sha256": LEGACY_ADAPTER_SHA256[
            "PLAN13_ADAPTER_READY"
        ],
        "adapter_registry_sha256": sha256_file(output_root / "adapter_registry.tsv"),
        "gate_status_sha256": sha256_file(output_root / "gate_status.tsv"),
        "source_manifest_sha256": sha256_file(output_root / "source_manifest.tsv"),
        "execution_manifest_sha256": sha256_file(
            output_root / "execution_manifest.tsv"
        ),
        "validation_report_sha256": sha256_file(
            output_root / "validation_report.tsv"
        ),
        "contract_lib_sha256": sha256_file(contract_root / "contract_lib.py"),
        "contract_validator_sha256": sha256_file(
            contract_root / "03_validate_contract.py"
        ),
        "producer_sha256": sha256_file(SCRIPT_PATH),
        "canonical_promotion_authorized": "FALSE",
    }
    for field, value in expected.items():
        if ready.get(field) != value:
            raise ContractError(f"semantic-v2 READY drift for {field}")
    _, effects = read_tsv(
        output_root / "program_effects.tsv", SCHEMAS["program_effects.tsv"]
    )
    if len(effects) != 2 or any(row["evidence_state"] != "indeterminate" for row in effects):
        raise ContractError("semantic-v2 Yakubovsky effects are not both indeterminate")
    if not any(not parse_bool(row["direction_agreement"], "direction_agreement") for row in effects):
        raise ContractError("semantic-v2 Yakubovsky direction disagreement was hidden")
    return {
        "release_id": RELEASE_ID,
        "semantic_contract_id": SEMANTIC_CONTRACT_ID,
        "status": "validated_plan13_adapter_candidate",
        "n_contract_checks": len(checks),
        "n_programs": len(effects),
        "n_indeterminate": sum(row["evidence_state"] == "indeterminate" for row in effects),
        "native_outputs_modified": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--native-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.check_only == args.write:
        parser.error("choose exactly one of --check-only or --write")
    project_root = args.base.resolve()
    native_root = (
        args.native_dir
        or (
            project_root
            / "Analysis/Spatial/candidates"
            / RELEASE_ID
            / "yakubovsky2026"
        )
    ).resolve()
    output_root = (
        args.output_dir
        or (
            project_root
            / "Analysis/Spatial/candidates"
            / RELEASE_ID
            / "semantic_v2_2026-08-08/yakubovsky_plan13_adapter"
        )
    ).resolve()
    try:
        if args.write:
            build_semantic_adapter(project_root, native_root, output_root)
        result = check_semantic_adapter(project_root, native_root, output_root)
        import json

        print(json.dumps(result, indent=2, sort_keys=True))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
