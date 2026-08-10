#!/usr/bin/env python3
"""Independent SP-INT-06/07 validation and candidate READY sealing."""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

from contract_lib import (
    ContractError,
    RELEASE_ID,
    SCHEMAS,
    parse_bool,
    parse_float,
    parse_int,
    read_tsv,
    sha256_file,
    write_tsv,
)
from final_integration_lib import (
    FIGURE4_MATRIX_COLUMNS,
    HOTFIX_EXPECTED_SHA256,
    INPUT_SEAL_COLUMNS,
    INTEGRATED_EFFECT_COLUMNS,
    SOURCE_TABLE_MANIFEST_COLUMNS,
    VERDICT_COLUMNS,
    canonical_row_sha256,
    expected_programs,
    figure4_matrix,
    independent_bh,
    make_verdicts,
    require_complete_groups,
    require_no_cross_assay_construct,
    verify_hotfix,
)


REPORT_COLUMNS = ("scope", "check_id", "status", "detail")
BUNDLE_COLUMNS = ("relative_path", "bytes", "sha256", "role")
READY_COLUMNS = (
    "release_id",
    "status",
    "mode",
    "input_seal_manifest_sha256",
    "integrated_effects_sha256",
    "figure4_matrix_sha256",
    "figure4_verdict_sha256",
    "source_table_manifest_sha256",
    "validation_report_sha256",
    "bundle_manifest_sha256",
    "code_freeze_sha256",
    "validator_hotfix_manifest_sha256",
    "n_programs",
    "n_dataset_assay_groups",
    "n_integrated_program_rows",
    "n_robust_rows",
    "n_indeterminate_rows",
    "n_tested_negative_rows",
    "n_terminal_rows",
    "v1_preservation_unchanged",
    "cross_assay_score_constructed",
    "pdf_rendered",
    "canonical_promotion_authorized",
    "plan13_complete",
    "validated_utc",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def check(condition: bool, check_id: str, detail: str, checks: list[dict[str, str]]) -> None:
    if not condition:
        raise ContractError(f"{check_id}: {detail}")
    checks.append({"scope": "SP-INT-06/07", "check_id": check_id, "status": "pass", "detail": detail})


def numeric_empty(row: Mapping[str, str]) -> bool:
    return not any(
        row[column].strip()
        for column in (
            "estimate", "std_error", "matched_null_sd", "interval_low", "interval_high",
            "pvalue", "padj", "heterogeneity_statistic", "heterogeneity_df",
            "heterogeneity_pvalue",
        )
    )


def validate_input_seals(project_root: Path, root: Path, mode: str, checks: list[dict[str, str]]) -> list[dict[str, str]]:
    _, rows = read_tsv(root / "input_seal_manifest.tsv", INPUT_SEAL_COLUMNS)
    check(bool(rows), "input_seal_nonempty", f"{len(rows)} pinned inputs", checks)
    seen = set()
    roles = set()
    for row in rows:
        path = project_root / row["relative_path"]
        check(row["relative_path"] not in seen, "unique_input_path", row["relative_path"], checks)
        seen.add(row["relative_path"])
        roles.add(row["source_role"])
        if not path.is_file() or path.stat().st_size != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
            raise ContractError(f"input seal drift: {path}")
    if mode == "real":
        required_roles = {
            "real_ready", "protein_ready", "plan12_gate", "yak_ready", "yak_adapter_ready",
            "yak_hotfix", "native_v1_ready", "native_v2_ready",
            "native_spatial_v2_program_effects", "native_spatial_v2_output_manifest",
            "outcome_blind_final_assembly_code_freeze",
        }
        missing = required_roles - roles
        check(not missing, "required_real_seals", f"all required roles present; missing={sorted(missing)}", checks)
    else:
        check(
            roles == {"fixture_registry_ready_only_no_external_outcome", "fixture_registry_only_no_external_outcome"},
            "fixture_inputs_outcome_free",
            "only Plan 20 registry/READY are read",
            checks,
        )
    return rows


def validate_effect_semantics(rows: Sequence[Mapping[str, str]], expected_uids: set[str], mode: str, checks: list[dict[str, str]]) -> None:
    require_complete_groups(rows, expected_uids)
    require_no_cross_assay_construct(rows, INTEGRATED_EFFECT_COLUMNS)
    unique_keys = {
        (row["dataset"], row["assay"], row["analysis_set_id"], row["program_uid"])
        for row in rows
    }
    check(len(unique_keys) == len(rows), "unique_integrated_key", f"{len(rows)} unique rows", checks)
    for row in rows:
        if row["program_uid"] not in expected_uids:
            raise ContractError(f"out-of-family program in integrated table: {row['program_uid']}")
        if not parse_bool(row["figure4_include"], f"include[{row['dataset']}/{row['program_uid']}]"):
            raise ContractError("a preselected program row was removed from Figure 4 source data")
        state = row["evidence_state"]
        testable = parse_bool(row["testable"], f"testable[{row['dataset']}/{row['program_uid']}]")
        robust = parse_bool(row["robustness_pass"], f"robust[{row['dataset']}/{row['program_uid']}]")
        comparable = parse_bool(row["cross_assay_comparable"], f"comparable[{row['dataset']}/{row['program_uid']}]")
        if comparable:
            raise ContractError("cross-assay comparability is prohibited")
        if state in {"untestable", "not_applicable", "skipped"}:
            if testable or not numeric_empty(row) or row["interval_type"] != "not_applicable" or row["pvalue_method"] != "not_applicable":
                raise ContractError(f"terminal integrated state carries inference: {row['dataset']}/{row['program_uid']}")
        elif state in {"robust", "indeterminate", "tested_negative"}:
            if not testable:
                raise ContractError("tested integrated state is marked untestable")
            estimate = parse_float(row["estimate"], "estimate")
            pvalue = parse_float(row["pvalue"], "pvalue")
            padj = parse_float(row["padj"], "padj")
            if None in (estimate, pvalue, padj) or not 0 <= float(pvalue) <= 1 or not 0 <= float(padj) <= 1:
                raise ContractError("tested integrated state lacks finite p/q/effect")
            if (
                not row["std_error"].strip()
                and not row["matched_null_sd"].strip()
                and not (row["interval_low"].strip() and row["interval_high"].strip())
            ):
                raise ContractError("tested integrated state lacks assay-native uncertainty")
            if row["std_error"].strip() and row["matched_null_sd"].strip():
                raise ContractError("one dispersion is labeled as both sampling SE and matched-null SD")
            descriptive = "positive" if float(estimate) > 0 else "negative" if float(estimate) < 0 else "zero"
            if row["direction_observed"] != descriptive or row["descriptive_effect_direction"] != descriptive:
                raise ContractError("descriptive effect direction does not match estimate")
            inferential = row["inferential_test_direction"]
            if inferential not in {"positive", "negative", "zero"}:
                raise ContractError("tested integrated state lacks inferential direction")
            if parse_bool(row["direction_agreement"], "direction_agreement") != (descriptive == inferential):
                raise ContractError("descriptive/inferential direction-agreement drift")
            heterogeneity = tuple(
                row[column].strip()
                for column in ("heterogeneity_statistic", "heterogeneity_df", "heterogeneity_pvalue")
            )
            if any(heterogeneity) and not all(heterogeneity):
                raise ContractError("heterogeneity summary is only partially populated")
            if state == "robust" and not robust:
                raise ContractError("robust state lacks robustness flag")
            if state == "tested_negative" and robust:
                raise ContractError("tested-negative state carries robustness flag")
            if state == "tested_negative" and not row["negative_call_rule_id"].strip():
                raise ContractError("tested-negative state lacks an explicit adequate-negative rule")
            if state != "tested_negative" and row["negative_call_rule_id"].strip():
                raise ContractError("non-negative state carries a negative-call rule")
            if state == "indeterminate" and robust:
                raise ContractError("indeterminate state carries robustness flag")
        elif state == "source_dependent":
            if row["source_dependence"] != "source_dependent":
                raise ContractError("source-dependent state/label mismatch")
        else:
            raise ContractError(f"unsupported integrated state: {state}")
        if mode == "real" and row["source_dependence"] not in {"independent", "source_dependent", "descriptive"}:
            raise ContractError("invalid source-dependence class")
    counts = Counter(row["evidence_state"] for row in rows)
    check(True, "state_semantics", f"states={dict(sorted(counts.items()))}; positivity not required", checks)
    check(all(parse_bool(row["figure4_include"], "include") for row in rows), "all_preselected_rows_retained", f"{len(rows)}/{len(rows)}", checks)
    check(True, "no_cross_assay_score", "assay-native units retained; no score/rank/count composite", checks)


def validate_matrix_and_verdict(root: Path, rows: list[dict[str, str]], checks: list[dict[str, str]]) -> None:
    _, matrix = read_tsv(root / "figure4_program_matrix.tsv", FIGURE4_MATRIX_COLUMNS)
    expected_matrix = figure4_matrix(rows)
    matrix_key = lambda row: (row["dataset"], row["assay"], row["analysis_set_id"], row["program_uid"])
    observed_by_key = {matrix_key(row): row for row in matrix}
    expected_by_key = {matrix_key(row): {column: str(row[column]) for column in FIGURE4_MATRIX_COLUMNS} for row in expected_matrix}
    if set(observed_by_key) != set(expected_by_key) or len(observed_by_key) != len(matrix):
        raise ContractError("Figure 4 matrix key grid drift")
    for key in expected_by_key:
        if observed_by_key[key] != expected_by_key[key]:
            raise ContractError(f"Figure 4 matrix does not exactly project integrated effects: {key}")
    check(True, "figure4_matrix_exact", f"{len(matrix)} exact rows", checks)

    _, verdicts = read_tsv(root / "figure4_dataset_verdict.tsv", VERDICT_COLUMNS)
    expected = make_verdicts(rows)
    verdict_key = lambda row: (row["dataset"], row["assay"], row["analysis_set_id"])
    observed_v = {verdict_key(row): row for row in verdicts}
    expected_v = {verdict_key(row): {column: str(row[column]) for column in VERDICT_COLUMNS} for row in expected}
    if observed_v != expected_v or len(observed_v) != len(verdicts):
        raise ContractError("Figure 4 dataset verdict does not independently rederive")
    check(True, "figure4_verdict_rederived", f"{len(verdicts)} dataset/assay verdicts", checks)


def validate_source_tables(project_root: Path, root: Path, mode: str, checks: list[dict[str, str]]) -> None:
    _, rows = read_tsv(root / "figure4_source_table_manifest.tsv", SOURCE_TABLE_MANIFEST_COLUMNS)
    expected_count = 8 if mode == "real" else 3
    check(len(rows) == expected_count, "source_table_count", f"{len(rows)} source tables", checks)
    ids = set()
    for row in rows:
        if row["table_id"] in ids:
            raise ContractError("duplicate Figure 4 source-table ID")
        ids.add(row["table_id"])
        path = root / row["relative_path"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
            raise ContractError(f"Figure 4 source-table manifest drift: {path}")
        with path.open("r", encoding="utf-8") as handle:
            n_rows = max(sum(1 for _ in handle) - 1, 0)
        if n_rows != int(row["data_rows"]):
            raise ContractError(f"Figure 4 source-table row count drift: {path}")
    if mode == "real":
        candidate_root = root.parent
        expected_copies = {
            "source_native/govaere_cosmx_il32_context.tsv": candidate_root / "real_adapters/govaere2026_cosmx_il32_context/source_native_context.tsv",
            "source_native/fixed_25_protein_context.tsv": candidate_root / "protein_atac/adapters/pxd051911_diams_program/fixed_25_protein_context.tsv",
            "source_native/static_atac_context.tsv": candidate_root / "protein_atac/source_native/atac_static_context/static_atac_context.tsv",
            "source_native/static_atac_hits.tsv": candidate_root / "protein_atac/source_native/atac_static_context/static_atac_hits.tsv",
            "source_native/static_atac_study_coverage.tsv": candidate_root / "protein_atac/source_native/atac_static_context/static_study_coverage.tsv",
        }
        for relative, source in expected_copies.items():
            target = root / relative
            if target.stat().st_size != source.stat().st_size or sha256_file(target) != sha256_file(source):
                raise ContractError(f"source-native Figure 4 table is not byte-identical: {relative}")
        check(True, "source_native_byte_identity", "CosMx, fixed-25 protein, and static ATAC tables preserved", checks)


def validate_adapter_imports(project_root: Path, rows: list[dict[str, str]], checks: list[dict[str, str]]) -> None:
    imported = [
        row for row in rows
        if row["import_method"] in {
            "exact_adapter_effect_row",
            "adjudicated_unresolved_array_not_biological_replicate",
        }
    ]
    by_file: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in imported:
        by_file[row["source_effect_path"]].append(row)
        source = project_root / row["source_effect_path"]
        ready = project_root / row["source_ready_path"]
        if sha256_file(source) != row["source_effect_sha256"] or sha256_file(ready) != row["source_ready_sha256"]:
            raise ContractError("adapter source/READY hash drift in integrated row")
    for relative, integrated in by_file.items():
        source = project_root / relative
        _, source_rows = read_tsv(source, SCHEMAS["program_effects.tsv"])
        source_by_uid = {row["program_uid"]: row for row in source_rows}
        if len(source_by_uid) != len(source_rows):
            raise ContractError("duplicate source adapter program row")
        for row in integrated:
            uid = row["program_uid"]
            if uid not in source_by_uid:
                raise ContractError("integrated adapter row absent from source")
            source_row = source_by_uid[uid]
            adjudicated = row["import_method"] == "adjudicated_unresolved_array_not_biological_replicate"
            exact_columns = tuple(
                column for column in SCHEMAS["program_effects.tsv"]
                if not adjudicated or column not in {"biological_unit", "n_biological"}
            )
            if any(row[column] != source_row[column] for column in exact_columns):
                raise ContractError(f"integrated adapter values drift from source: {relative}/{uid}")
            if adjudicated and not (
                row["dataset"] == "Govaere2026_CosMx"
                and row["biological_unit_resolution"] == "unresolved"
                and row["biological_unit"] == "unknown_public_biological_unit"
                and row["n_biological"] == ""
                and row["technical_unit"] == "physical_array"
                and row["n_technical"] == source_row["n_technical"]
            ):
                raise ContractError("CosMx unresolved-array biological-count adjudication drift")
            if row["source_effect_row_sha256"] != canonical_row_sha256(source_row, SCHEMAS["program_effects.tsv"]):
                raise ContractError("integrated adapter row hash drift")
    check(
        len(imported) == 14,
        "adapter_imports",
        "12 exact plus 2 CosMx unresolved-array count-adjudicated rows",
        checks,
    )


def validate_native_spatial(project_root: Path, root: Path, rows: list[dict[str, str]], checks: list[dict[str, str]]) -> None:
    native = [row for row in rows if row["import_method"] == "deterministic_native_spatial_adapter"]
    if len(native) != 4:
        raise ContractError("final integration must contain four native spatial rows")
    source_path = root.parent / "native_spatial/v2_candidate/spatial_program_results.tsv"
    _, source_rows = read_tsv(
        source_path,
        (
            "program_id", "dataset", "testable", "residual_moran_i", "residual_null_mean",
            "residual_null_sd", "residual_pvalue", "residual_padj", "sensitivity_sign_agree", "robust",
        ),
    )
    source = {(row["dataset"], row["program_id"]): row for row in source_rows}
    for row in native:
        key = (row["dataset"], row["program_uid"])
        src = source[key]
        producer = project_root / row["producer"]
        if not producer.is_file() or sha256_file(producer) != row["producer_sha256"]:
            raise ContractError("native spatial candidate-producer provenance drift")
        if parse_bool(src["testable"], "native.testable"):
            expected_estimate = float(src["residual_moran_i"]) - float(src["residual_null_mean"])
            for column, expected in (
                ("estimate", expected_estimate),
                ("matched_null_sd", float(src["residual_null_sd"])),
                ("pvalue", float(src["residual_pvalue"])),
                ("padj", float(src["residual_padj"])),
            ):
                observed = parse_float(row[column], f"native.{column}")
                if observed is None or not math.isclose(observed, expected, rel_tol=1e-12, abs_tol=1e-14):
                    raise ContractError(f"native spatial {column} does not rederive: {key}")
            if row["uncertainty_semantics"] != "matched_gene_null_standard_deviation_not_sampling_standard_error":
                raise ContractError("native spatial null SD is mislabeled as sampling uncertainty")
            if row["std_error"].strip() or row["interval_type"] != "none":
                raise ContractError("matched-null SD was mislabeled as sampling SE/CI")
            expected_robust = parse_bool(src["robust"], "native.robust")
            if parse_bool(row["robustness_pass"], "integrated.robust") != expected_robust:
                raise ContractError("native spatial robust flag drift")
        if row["source_dependence"] != ("independent" if row["dataset"] == "GSE192741" else "source_dependent"):
            raise ContractError("native spatial source-dependence label drift")
        if row["dataset"] == "Vu_et_al_2025" and not (
            row["biological_unit_resolution"] == "unresolved"
            and row["biological_unit"] == "unknown_public_biological_unit"
            and row["n_biological"] == ""
            and int(row["n_technical"]) > 0
        ):
            raise ContractError("Vu arrays were reported as biological replicates")
    for dataset in ("GSE192741", "Vu_et_al_2025"):
        src_part = {row["program_id"]: float(row["residual_pvalue"]) for row in source_rows if row["dataset"] == dataset}
        expected_q = independent_bh(src_part)
        for row in native:
            if row["dataset"] == dataset:
                observed = float(row["padj"])
                if not math.isclose(observed, expected_q[row["program_uid"]], rel_tol=1e-12, abs_tol=1e-14):
                    raise ContractError(f"native spatial complete-family BH drift: {dataset}")
    check(True, "native_spatial_rederived", "4 effects, matched-null uncertainty, two complete BH families", checks)


def validate_real_specific(project_root: Path, root: Path, rows: list[dict[str, str]], checks: list[dict[str, str]]) -> None:
    hotfix = project_root / "Analysis/Spatial/candidates" / RELEASE_ID / "yakubovsky2026/validator_hotfix_manifest.tsv"
    verify_hotfix(hotfix)
    check(sha256_file(hotfix) == HOTFIX_EXPECTED_SHA256, "yak_hotfix_provenance", HOTFIX_EXPECTED_SHA256, checks)
    validate_adapter_imports(project_root, rows, checks)
    validate_native_spatial(project_root, root, rows, checks)
    groups = {(row["dataset"], row["assay"], row["analysis_set_id"]) for row in rows}
    check(len(groups) == 9 and len(rows) == 18, "real_complete_grid", "nine assay groups x two frozen programs", checks)
    sources = Counter(row["source_dependence"] for row in rows)
    check(sources["independent"] > 0 and sources["source_dependent"] > 0, "source_dependence_explicit", str(dict(sources)), checks)
    summary = root.parent / "v1_preservation_summary.tsv"
    _, summary_rows = read_tsv(summary, ("n_unchanged", "n_changed", "n_missing", "n_added", "aggregate_sha256"))
    if len(summary_rows) != 1 or summary_rows[0]["n_unchanged"] != "591" or any(
        summary_rows[0][column] != "0" for column in ("n_changed", "n_missing", "n_added")
    ):
        raise ContractError("v1 preservation summary is not 591/591 unchanged")
    check(True, "v1_preservation", f"591/591; {summary_rows[0]['aggregate_sha256']}", checks)


def validate_assembly_status(root: Path, mode: str, checks: list[dict[str, str]]) -> None:
    _, rows = read_tsv(root / "assembly_status.tsv")
    if len(rows) != 1:
        raise ContractError("assembly status must contain one row")
    row = rows[0]
    expected_status = "assembled_pending_independent_validation" if mode == "real" else "synthetic_fixture_assembled_pending_validation"
    if row.get("status") != expected_status or row.get("mode") != mode:
        raise ContractError("assembly status/mode drift")
    links = {
        "input_seal_manifest_sha256": root / "input_seal_manifest.tsv",
        "integrated_effects_sha256": root / "integrated_program_effects.tsv",
        "figure4_matrix_sha256": root / "figure4_program_matrix.tsv",
        "figure4_verdict_sha256": root / "figure4_dataset_verdict.tsv",
        "source_table_manifest_sha256": root / "figure4_source_table_manifest.tsv",
    }
    for field, path in links.items():
        if row.get(field) != sha256_file(path):
            raise ContractError(f"assembly status hash link drift: {field}")
    for field in ("cross_assay_score_constructed", "pdf_rendered", "canonical_write_authorized", "plan13_complete"):
        if parse_bool(row[field], f"status.{field}"):
            raise ContractError(f"pre-validation assembly status overstates {field}")
    check(True, "assembly_status_links", "all source-table hashes linked; no PDF/canonical write", checks)


def run_checks(project_root: Path, root: Path, mode: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    checks: list[dict[str, str]] = []
    validate_input_seals(project_root, root, mode, checks)
    hotspot_root = root.parent.parent / "hotspot"
    selected = expected_programs(hotspot_root)
    _, rows = read_tsv(root / "integrated_program_effects.tsv", INTEGRATED_EFFECT_COLUMNS)
    validate_effect_semantics(rows, set(selected), mode, checks)
    expected_rows = 18 if mode == "real" else 8
    expected_groups = 9 if mode == "real" else 4
    check(len(rows) == expected_rows, "integrated_row_count", str(expected_rows), checks)
    check(
        len({(row["dataset"], row["assay"], row["analysis_set_id"]) for row in rows}) == expected_groups,
        "integrated_group_count",
        str(expected_groups),
        checks,
    )
    validate_matrix_and_verdict(root, rows, checks)
    validate_source_tables(project_root, root, mode, checks)
    validate_assembly_status(root, mode, checks)
    if mode == "real":
        validate_real_specific(project_root, root, rows, checks)
    else:
        states = Counter(row["evidence_state"] for row in rows)
        check(
            states == Counter({"indeterminate": 2, "skipped": 2, "untestable": 2, "robust": 1, "tested_negative": 1}),
            "fixture_state_coverage",
            str(dict(states)),
            checks,
        )
    return checks, rows


def bundle_rows(root: Path) -> list[dict[str, object]]:
    excluded = {"validation_report.tsv", "bundle_manifest.tsv", "READY", "FIXTURE_READY"}
    rows = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name in excluded:
            continue
        rows.append(
            {
                "relative_path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "role": "sp_int_06_07_candidate_artifact",
            }
        )
    return rows


def write_ready(project_root: Path, root: Path, mode: str, checks: list[dict[str, str]], rows: list[dict[str, str]]) -> None:
    report = root / "validation_report.tsv"
    bundle = root / "bundle_manifest.tsv"
    ready = root / ("READY" if mode == "real" else "FIXTURE_READY")
    for path in (report, bundle, ready):
        if path.exists():
            raise ContractError(f"refusing to overwrite final validation artifact: {path}")
    write_tsv(report, REPORT_COLUMNS, checks)
    write_tsv(bundle, BUNDLE_COLUMNS, bundle_rows(root))
    states = Counter(row["evidence_state"] for row in rows)
    code_freeze = root.parent / "final_integration_code_freeze.tsv"
    ready_row = {
        "release_id": RELEASE_ID,
        "status": (
            "ready_plan13_candidate_source_tables_no_pdf_no_canonical_write"
            if mode == "real"
            else "validated_synthetic_null_fixture_no_real_outcomes"
        ),
        "mode": mode,
        "input_seal_manifest_sha256": sha256_file(root / "input_seal_manifest.tsv"),
        "integrated_effects_sha256": sha256_file(root / "integrated_program_effects.tsv"),
        "figure4_matrix_sha256": sha256_file(root / "figure4_program_matrix.tsv"),
        "figure4_verdict_sha256": sha256_file(root / "figure4_dataset_verdict.tsv"),
        "source_table_manifest_sha256": sha256_file(root / "figure4_source_table_manifest.tsv"),
        "validation_report_sha256": sha256_file(report),
        "bundle_manifest_sha256": sha256_file(bundle),
        "code_freeze_sha256": sha256_file(code_freeze) if code_freeze.is_file() else "not_applicable_fixture_pre_freeze",
        "validator_hotfix_manifest_sha256": HOTFIX_EXPECTED_SHA256 if mode == "real" else "not_applicable_fixture",
        "n_programs": 2,
        "n_dataset_assay_groups": 9 if mode == "real" else 4,
        "n_integrated_program_rows": len(rows),
        "n_robust_rows": states["robust"],
        "n_indeterminate_rows": states["indeterminate"],
        "n_tested_negative_rows": states["tested_negative"],
        "n_terminal_rows": states["untestable"] + states["not_applicable"] + states["skipped"],
        "v1_preservation_unchanged": mode == "real",
        "cross_assay_score_constructed": False,
        "pdf_rendered": False,
        "canonical_promotion_authorized": False,
        "plan13_complete": mode == "real",
        "validated_utc": utc_now(),
    }
    write_tsv(ready, READY_COLUMNS, [ready_row])


def check_ready(root: Path, mode: str) -> None:
    path = root / ("READY" if mode == "real" else "FIXTURE_READY")
    _, rows = read_tsv(path, READY_COLUMNS)
    if len(rows) != 1:
        raise ContractError("final READY must contain one row")
    row = rows[0]
    expected_status = (
        "ready_plan13_candidate_source_tables_no_pdf_no_canonical_write"
        if mode == "real"
        else "validated_synthetic_null_fixture_no_real_outcomes"
    )
    if row["status"] != expected_status or row["mode"] != mode:
        raise ContractError("final READY status/mode drift")
    for field, file_name in (
        ("input_seal_manifest_sha256", "input_seal_manifest.tsv"),
        ("integrated_effects_sha256", "integrated_program_effects.tsv"),
        ("figure4_matrix_sha256", "figure4_program_matrix.tsv"),
        ("figure4_verdict_sha256", "figure4_dataset_verdict.tsv"),
        ("source_table_manifest_sha256", "figure4_source_table_manifest.tsv"),
        ("validation_report_sha256", "validation_report.tsv"),
        ("bundle_manifest_sha256", "bundle_manifest.tsv"),
    ):
        if row[field] != sha256_file(root / file_name):
            raise ContractError(f"final READY hash drift: {field}")
    if parse_bool(row["cross_assay_score_constructed"], "READY.score") or parse_bool(row["pdf_rendered"], "READY.pdf"):
        raise ContractError("final READY claims a prohibited score or PDF")
    if parse_bool(row["canonical_promotion_authorized"], "READY.promotion"):
        raise ContractError("final READY authorizes canonical promotion")
    if parse_bool(row["plan13_complete"], "READY.plan13") != (mode == "real"):
        raise ContractError("final READY Plan 13 state drift")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("fixture", "real"), required=True)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--candidate-root", type=Path, default=None)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    candidate_root = args.candidate_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "spatial_context"
    )
    root = candidate_root / ("final_integration_fixture" if args.mode == "fixture" else "final_integration")
    try:
        checks, rows = run_checks(project_root, root, args.mode)
        if args.check_only:
            check_ready(root, args.mode)
            print(f"PASS: rederived {len(checks)} {args.mode} final-integration checks and READY hashes")
        else:
            write_ready(project_root, root, args.mode, checks, rows)
            check_ready(root, args.mode)
            print(f"PASS: {len(checks)} {args.mode} final-integration checks; candidate READY")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
