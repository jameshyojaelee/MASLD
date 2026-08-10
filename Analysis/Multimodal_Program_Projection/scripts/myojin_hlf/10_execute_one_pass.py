#!/usr/bin/env python3
"""Execute the sealed Myojin HLF analysis exactly once after all gates pass."""

from __future__ import annotations

import argparse
import datetime as dt
import math
import platform
import re
import traceback
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import numpy as np
import pandas as pd
import scipy
import statsmodels

from myojin_firewall_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    RELEASE_ID,
    atomic_write_text,
    md5_file,
    parse_bool,
    read_tsv,
    require_within,
    sha256_file,
    stable_bundle_sha256,
    write_tsv,
)
from phase_c_stats import (
    add_frozen_strata,
    apply_program_bh,
    fit_class_binary,
    fit_class_continuous,
    matched_program_test,
    mechanical_figure_verdict,
    rank_first_bins,
)


SPECIFICATION_SHA256 = (
    "1f99b7a7a95d9aebbbb911bc2317860f4a3fada0d1da2a4932d5c1ab190c4cdd"
)
RAW_SHA256 = "2663a811676a1abe2778375e6d28f6397eaded1f895d274f3c241733669c1c66"
CLASS_SEED = 2026080740
BINARY_SEED = 2026080741
PROGRAM_SEED = 2026080742
CLASS_DRAWS = 100_000
BINARY_DRAWS = 100_000
PROGRAM_DRAWS = 10_000
EXPECTED_VERSIONS = {
    "python": "3.10",
    "numpy": "2.2.6",
    "pandas": "2.3.3",
    "scipy": "1.15.2",
    "statsmodels": "0.14.5",
}
EXPECTED_DEPMAP = {
    "Model.csv": (645696, "675210d17675f3517b0ce39a3c274f16"),
    "CRISPRGeneEffect.csv": (428678699, "6edf7ade09b9b34199210b559d4745d3"),
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": (
        506628654,
        "71794802b750ce77c422dad0720a40af",
    ),
}
EXPECTED_GATES = {
    "masked_source_ready",
    "external_outcomes_unread",
    "frozen_evidence_class_hash",
    "frozen_hotspot_registry_hashes",
    "gencode_v49_hash",
    "version_matched_depmap_triplet",
    "HLF_identity_unique_ACH_000393",
    "HLF_expression_available",
    "depmap_entrez_conflicts_explicitly_excluded",
    "mapping_unambiguous_for_primary_universe",
    "HLF_TPM_ge_1_eligibility_attrition_reported",
    "nonreference_class_complete_covariates_at_least_80pct",
}
SPEC_COMPONENTS = [
    "blind_analysis_spec.yaml",
    "prediction_manifest.tsv",
    "known_hit_exclusion.tsv",
    "screen_gene_universe.tsv",
    "gene_mapping_audit.tsv",
    "covariate_coverage_audit.tsv",
    "program_testability.tsv",
    "gate_status.tsv",
    "firewall_input_manifest.tsv",
]
RESULT_NAMES = [
    "screen_results.tsv",
    "class_effects.tsv",
    "class_protective_hits.tsv",
    "program_effects.tsv",
    "permutation_audit.tsv",
    "matched_null_audit.tsv",
    "sensitivity.tsv",
    "fig5_verdict.tsv",
    "phase_c_release_manifest.tsv",
    "PHASE_C_COMPLETE",
    "PHASE_C_VALIDATED",
]
NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS = {"m": NS_MAIN}
EXPECTED_HEADERS = {
    "A3": "id",
    "B3": "num",
    "C3": "pos|score",
    "D3": "pos|p-value",
    "E3": "pos|fdr",
    "F3": "pos|rank",
    "G3": "pos|goodsgrna",
    "H3": "pos|lfc",
    "I3": "-Log10pvalue_Pos",
    "J3": "num",
    "K3": "pos|score",
    "L3": "pos|p-value",
    "M3": "pos|fdr",
    "N3": "pos|rank",
    "O3": "pos|goodsgrna",
    "P3": "pos|lfc",
    "Q3": "-Log10pvalue_Pos",
}
SCREEN_FIELDS = [
    "source_row",
    "source_gene_id",
    "gene_symbol",
    "mapping_state",
    "gene_biotype",
    "primary_evidence_class",
    "expected_direction",
    "depmap_entrez_conflict",
    "neutral_source_base",
    "primary_covariates_numeric_complete",
    "HLF_expression_floor_pass",
    "direct_primary_fields_complete",
    "primary_eligible",
    "guide5_sensitivity_eligible",
    "exclusion_reasons",
    "guide_count_d21_d0",
    "guide_count_pa_vehicle",
    "guide_count_min",
    "HLF_TPM",
    "log1p_HLF_TPM",
    "HLF_Chronos",
    "D21_vs_D0_positive_selection_score",
    "D21_vs_D0_positive_selection_p",
    "D21_vs_D0_positive_selection_FDR",
    "D21_vs_D0_positive_selection_LFC",
    "PA_vs_vehicle_positive_selection_score",
    "PA_vs_vehicle_positive_selection_p",
    "PA_vs_vehicle_positive_selection_FDR",
    "PA_vs_vehicle_positive_selection_LFC",
    "protective_hit_FDR_lt_0.05_LFC_ge_0.25",
    "positive_selection_FDR_hit_without_LFC_floor",
    "source_significance_interpretation",
    "analysis_guide_bin",
    "analysis_expression_quintile",
    "analysis_chronos_quintile",
    "actual_vehicle_fitness_observed",
    "outcome_values_read",
    "specification_sha256",
    "release_id",
]


class PreflightError(RuntimeError):
    pass


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def one_row(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    if len(rows) != 1:
        raise PreflightError(f"expected_one_row::{path.name}::{len(rows)}")
    return rows[0]


def preflight_check(
    audit: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    audit.append(
        {
            "check": name,
            "passed": str(passed).upper(),
            "detail": detail,
            "external_outcomes_read": "FALSE",
            "release_id": RELEASE_ID,
        }
    )
    if not passed:
        raise PreflightError(f"{name}::{detail}")


def validate_environment() -> dict[str, str]:
    return {
        "python": ".".join(platform.python_version().split(".")[:2]),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "statsmodels": statsmodels.__version__,
    }


def verify_code_manifest(outdir: Path, audit: list[dict[str, object]]) -> str:
    marker = one_row(outdir / "PHASE_C_CODE_READY")
    manifest_path = outdir / "phase_c_code_manifest.tsv"
    preflight_check(
        audit,
        "phase_c_code_manifest_hash",
        marker.get("manifest_sha256") == sha256_file(manifest_path),
        marker.get("manifest_sha256", "missing"),
    )
    rows = read_tsv(manifest_path)
    code_paths: list[Path] = []
    input_paths: list[Path] = []
    all_match = bool(rows)
    environment_rows: dict[str, str] = {}
    for row in rows:
        role = row["role"]
        if role == "environment_version":
            environment_rows[row["relative_path"].split("::", 1)[1]] = row["sha256"]
            continue
        path = PROJECT_ROOT / row["relative_path"]
        all_match &= (
            path.is_file()
            and str(path.stat().st_size) == row["bytes"]
            and sha256_file(path) == row["sha256"]
            and row["specification_sha256"] == SPECIFICATION_SHA256
        )
        if role == "phase_c_code":
            code_paths.append(path)
        elif role == "frozen_input":
            input_paths.append(path)
        else:
            all_match = False
    preflight_check(
        audit, "all_frozen_code_and_inputs_unchanged", all_match, f"rows={len(rows)}"
    )
    code_hash, _ = stable_bundle_sha256(code_paths, PROJECT_ROOT)
    input_hash, _ = stable_bundle_sha256(input_paths, PROJECT_ROOT)
    preflight_check(
        audit,
        "phase_c_code_bundle_rederived",
        code_hash == marker.get("code_bundle_sha256"),
        code_hash,
    )
    preflight_check(
        audit,
        "phase_c_input_bundle_rederived",
        input_hash == marker.get("input_bundle_sha256"),
        input_hash,
    )
    observed_environment = validate_environment()
    preflight_check(
        audit,
        "phase_c_environment_exact",
        observed_environment == EXPECTED_VERSIONS == environment_rows,
        str(observed_environment),
    )
    return code_hash


def preflight(outdir: Path) -> tuple[list[dict[str, object]], Path, str]:
    audit: list[dict[str, object]] = []
    required = [
        "SEALED",
        "FIREWALL_READY",
        "DEPMAP_READY",
        "PHASE_C_CODE_READY",
        "phase_c_code_manifest.tsv",
        "unseal_record.tsv",
        "masked_screen_schema.tsv",
        *SPEC_COMPONENTS,
        "source_manifest.tsv",
        "depmap_source_manifest.tsv",
        "hlf_identity_audit.tsv",
    ]
    missing = [name for name in required if not (outdir / name).is_file()]
    preflight_check(
        audit, "required_artifacts_present", not missing, f"missing={missing}"
    )
    forbidden = [
        name
        for name in [
            "BLOCKED",
            "UNSEALED",
            "phase_c_terminal_blocker.tsv",
            "phase_c_terminal_failure.tsv",
            *RESULT_NAMES,
        ]
        if (outdir / name).exists()
    ]
    preflight_check(
        audit, "phase_c_never_started", not forbidden, f"present={forbidden}"
    )
    preflight_check(
        audit,
        "unseal_record_empty",
        len(read_tsv(outdir / "unseal_record.tsv")) == 0,
        "append-only record has zero data rows",
    )

    sealed = one_row(outdir / "SEALED")
    ready = one_row(outdir / "FIREWALL_READY")
    authorized = (
        sealed.get("status") == "sealed_ready_for_unseal"
        and ready.get("status") == "ready_for_one_pass_unseal"
        and sealed.get("specification_sha256") == SPECIFICATION_SHA256
        and ready.get("specification_sha256") == SPECIFICATION_SHA256
        and sealed.get("external_outcomes_read") == "FALSE"
        and ready.get("external_outcomes_read") == "FALSE"
    )
    preflight_check(
        audit, "sealed_firewall_authorization", authorized, SPECIFICATION_SHA256
    )
    spec_hash, _ = stable_bundle_sha256(
        [outdir / name for name in SPEC_COMPONENTS], outdir
    )
    preflight_check(
        audit,
        "specification_bundle_rederived",
        spec_hash == SPECIFICATION_SHA256,
        spec_hash,
    )
    code_hash = verify_code_manifest(outdir, audit)

    gates = read_tsv(outdir / "gate_status.tsv")
    gate_names = {row["gate"] for row in gates}
    gates_pass = (
        len(gates) == 12
        and gate_names == EXPECTED_GATES
        and all(
            parse_bool(row["passed"]) and row["external_outcomes_read"] == "FALSE"
            for row in gates
        )
    )
    preflight_check(
        audit, "all_12_source_gates_pass", gates_pass, f"gates={sorted(gate_names)}"
    )
    deviations = read_tsv(outdir / "deviations.tsv")
    preflight_check(
        audit,
        "no_prespecification_deviations",
        not deviations,
        f"rows={len(deviations)}",
    )
    design_detail = preflight_design_gates(outdir)
    preflight_check(
        audit,
        "all_prespecified_designs_estimable_before_outcome_access",
        design_detail["passed"],
        design_detail["detail"],
    )

    source = one_row(outdir / "source_manifest.tsv")
    workbook = outdir / source["relative_path"]
    raw_ok = (
        workbook.is_file()
        and source.get("sha256") == RAW_SHA256
        and sha256_file(workbook) == RAW_SHA256
        and workbook.stat().st_size == int(source["bytes"])
        and workbook.stat().st_mode & 0o222 == 0
    )
    preflight_check(audit, "raw_workbook_authenticated_readonly", raw_ok, RAW_SHA256)

    depmap_manifest = read_tsv(outdir / "depmap_source_manifest.tsv")
    depmap_by_name = {row["filename"]: row for row in depmap_manifest}
    depmap_ok = set(depmap_by_name) == set(EXPECTED_DEPMAP)
    for filename, (expected_bytes, expected_md5) in EXPECTED_DEPMAP.items():
        path = outdir / "source/depmap_24q4" / filename
        row = depmap_by_name.get(filename, {})
        depmap_ok &= (
            path.is_file()
            and path.stat().st_size == expected_bytes
            and md5_file(path) == expected_md5
            and sha256_file(path) == row.get("observed_sha256")
            and row.get("published_md5") == expected_md5
            and row.get("observed_md5") == expected_md5
            and row.get("doi") == "10.25452/figshare.plus.27993248.v1"
            and row.get("depmap_release") == "DepMap Public 24Q4"
        )
    preflight_check(
        audit, "official_depmap_24q4_triplet_reauthenticated", depmap_ok, "3 files"
    )
    identity = one_row(outdir / "hlf_identity_audit.tsv")
    identity_ok = (
        identity.get("model_id") == "ACH-000393"
        and identity.get("cell_line_name") == "HLF"
        and identity.get("oncotree_lineage") == "Liver"
        and identity.get("hlf_unique") == "TRUE"
        and identity.get("hlfa_conflated") == "FALSE"
        and identity.get("actual_vehicle_fitness_observed") == "FALSE"
    )
    preflight_check(
        audit, "HLF_identity_reauthenticated", identity_ok, "ACH-000393 HLF liver"
    )
    return audit, workbook, code_hash


def frozen_covariate_frame(universe: list[dict[str, str]], predicate) -> pd.DataFrame:
    rows = []
    for row in universe:
        if not predicate(row):
            continue
        chronos = row["HLF_Chronos"]
        rows.append(
            {
                "source_row": int(row["source_row"]),
                "gene_symbol": row["gene_symbol"],
                "primary_evidence_class": row["primary_evidence_class"],
                "guide_count_min": int(row["guide_count_min"]),
                "log1p_HLF_TPM": float(row["log1p_HLF_TPM"]),
                "HLF_Chronos": np.nan if chronos == "" else float(chronos),
                "HLF_Chronos_missing": int(chronos == ""),
            }
        )
    frame = pd.DataFrame.from_records(rows)
    if frame.empty or frame["gene_symbol"].duplicated().any():
        raise PreflightError("empty_or_duplicate_frozen_design_universe")
    return frame.sort_values("gene_symbol", kind="stable").reset_index(drop=True)


def design_ranks(
    frame: pd.DataFrame, include_missing: bool
) -> tuple[int, int, int, int]:
    classes = frame["primary_evidence_class"].to_numpy()
    columns = [
        np.ones(len(frame)),
        (classes == "genetic_only").astype(float),
        (classes == "disease_state_only").astype(float),
        (classes == "convergent").astype(float),
        frame["log1p_HLF_TPM"].to_numpy(dtype=float),
        frame["guide_count_min"].to_numpy(dtype=float),
        frame["HLF_Chronos"].to_numpy(dtype=float),
    ]
    if include_missing:
        columns.append(frame["HLF_Chronos_missing"].to_numpy(dtype=float))
    full = np.column_stack(columns)
    reduced = np.column_stack([columns[0], *columns[4:]])
    return (
        int(np.linalg.matrix_rank(full)),
        full.shape[1],
        int(np.linalg.matrix_rank(reduced)),
        reduced.shape[1],
    )


def preflight_design_gates(outdir: Path) -> dict[str, object]:
    universe = read_tsv(outdir / "screen_gene_universe.tsv")
    primary = frozen_covariate_frame(
        universe, lambda row: parse_bool(row["primary_eligible"])
    )
    expanded = frozen_covariate_frame(
        universe,
        lambda row: (
            parse_bool(row["neutral_source_base"])
            and not parse_bool(row["depmap_entrez_conflict"])
            and parse_bool(row["HLF_expression_floor_pass"])
            and row["log1p_HLF_TPM"] != ""
        ),
    )
    primary, expanded = freeze_analysis_strata(primary, expanded)
    known_hits = {
        row["canonical_symbol"] for row in read_tsv(outdir / "known_hit_exclusion.tsv")
    }
    d21_rows = {
        int(row["source_row"])
        for row in read_tsv(outdir / "masked_screen_schema.tsv")
        if parse_bool(row["d21_d0_lfc_present"])
    }
    variants = {
        "primary": (primary, False),
        "guide_floor_5": (primary[primary["guide_count_min"] >= 5], False),
        "remove_all_known_hits": (
            primary[~primary["gene_symbol"].isin(known_hits)],
            False,
        ),
        "remove_ACSL3": (primary[primary["gene_symbol"] != "ACSL3"], False),
        "missing_Chronos_indicator": (expanded, True),
        "PA_D21_vs_D0": (primary[primary["source_row"].isin(d21_rows)], False),
        "absolute_PA_vs_vehicle_LFC": (primary, False),
    }
    failures = []
    details = []
    if len(primary) != 7_762:
        failures.append(f"primary_n={len(primary)}")
    if int(expanded["HLF_Chronos_missing"].sum()) == 0:
        failures.append("missing_Chronos_sensitivity_has_no_missing_rows")
    for name, (frame, include_missing) in variants.items():
        rank_full, columns_full, rank_reduced, columns_reduced = design_ranks(
            frame, include_missing
        )
        details.append(
            f"{name}:{len(frame)}:{rank_full}/{columns_full}:{rank_reduced}/{columns_reduced}"
        )
        if rank_full != columns_full or rank_reduced != columns_reduced:
            failures.append(f"{name}_rank_failure")
        strata = add_frozen_strata(frame)["permutation_stratum"]
        class_frame = frame.assign(permutation_stratum=strata.to_numpy())
        for class_name in ("genetic_only", "disease_state_only", "convergent"):
            shared = 0
            for _, group in class_frame[
                class_frame["primary_evidence_class"].isin(["neither", class_name])
            ].groupby("permutation_stratum"):
                labels = set(group["primary_evidence_class"])
                shared += labels == {"neither", class_name}
            if shared == 0:
                failures.append(f"{name}_{class_name}_no_shared_binary_stratum")
    return {
        "passed": not failures,
        "detail": "failures="
        + (",".join(failures) if failures else "none")
        + "|"
        + ";".join(details),
    }


def column_letters(reference: str) -> str:
    match = re.match(r"^([A-Z]+)", reference)
    if not match:
        raise ValueError(f"Invalid cell reference {reference}")
    return match.group(1)


def shared_strings(archive: ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return [
        "".join(node.text or "" for node in item.iterfind(".//m:t", NS))
        for item in root.findall("m:si", NS)
    ]


def cell_value(cell: ET.Element, strings: list[str]) -> str:
    cell_type = cell.attrib.get("t", "")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.iterfind(".//m:t", NS))
    value = cell.find("m:v", NS)
    raw = "" if value is None else (value.text or "")
    if cell_type == "s" and raw:
        return strings[int(raw)]
    return raw


def numeric(value: str, label: str) -> float | None:
    if value == "":
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite numeric result at {label}")
    return result


def parse_outcomes_once(workbook: Path) -> list[dict[str, object]]:
    with ZipFile(workbook) as archive:
        strings = shared_strings(archive)
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    dimension = sheet.find("m:dimension", NS)
    if dimension is None or dimension.attrib.get("ref") != "A1:Q18346":
        raise ValueError("Held-out workbook dimension changed")
    rows: dict[int, dict[str, str]] = {}
    for node in sheet.findall(".//m:row", NS):
        row_number = int(node.attrib["r"])
        values: dict[str, str] = {}
        for cell in node.findall("m:c", NS):
            values[column_letters(cell.attrib["r"])] = cell_value(cell, strings)
        rows[row_number] = values
    for reference, expected in EXPECTED_HEADERS.items():
        row_number = int(re.sub(r"^[A-Z]+", "", reference))
        if rows.get(row_number, {}).get(column_letters(reference), "") != expected:
            raise ValueError(f"Held-out workbook header changed at {reference}")
    parsed: list[dict[str, object]] = []
    for row_number in range(4, 18347):
        row = rows.get(row_number, {})
        source_gene_id = row.get("A", "").strip()
        if not source_gene_id:
            raise ValueError(f"Blank source gene at row {row_number}")
        values = {
            "source_row": row_number,
            "source_gene_id": source_gene_id,
            "D21_vs_D0_positive_selection_score": numeric(
                row.get("C", ""), f"C{row_number}"
            ),
            "D21_vs_D0_positive_selection_p": numeric(
                row.get("D", ""), f"D{row_number}"
            ),
            "D21_vs_D0_positive_selection_FDR": numeric(
                row.get("E", ""), f"E{row_number}"
            ),
            "D21_vs_D0_positive_selection_LFC": numeric(
                row.get("H", ""), f"H{row_number}"
            ),
            "PA_vs_vehicle_positive_selection_score": numeric(
                row.get("K", ""), f"K{row_number}"
            ),
            "PA_vs_vehicle_positive_selection_p": numeric(
                row.get("L", ""), f"L{row_number}"
            ),
            "PA_vs_vehicle_positive_selection_FDR": numeric(
                row.get("M", ""), f"M{row_number}"
            ),
            "PA_vs_vehicle_positive_selection_LFC": numeric(
                row.get("P", ""), f"P{row_number}"
            ),
        }
        for field in (
            "D21_vs_D0_positive_selection_p",
            "D21_vs_D0_positive_selection_FDR",
            "PA_vs_vehicle_positive_selection_p",
            "PA_vs_vehicle_positive_selection_FDR",
        ):
            value = values[field]
            if value is not None and not 0 <= float(value) <= 1:
                raise ValueError(
                    f"Probability outside [0,1] at row {row_number} field {field}"
                )
        parsed.append(values)
    if len(parsed) != 18_343:
        raise ValueError(f"Expected 18,343 outcome rows; observed {len(parsed)}")
    return parsed


def build_screen_results(
    universe: list[dict[str, str]],
    prediction: list[dict[str, str]],
    outcomes: list[dict[str, object]],
) -> list[dict[str, object]]:
    if len(universe) != len(outcomes) or len(universe) != 18_343:
        raise ValueError("Universe/outcome row-count mismatch")
    directions = {
        row["gene_symbol"]: row["expected_direction"]
        for row in prediction
        if row["object_type"] == "evidence_class_gene"
    }
    result: list[dict[str, object]] = []
    for frozen, outcome in zip(universe, outcomes):
        if (
            int(frozen["source_row"]) != int(outcome["source_row"])
            or frozen["source_gene_id"] != outcome["source_gene_id"]
        ):
            raise ValueError(
                f"Held-out source join mismatch at row {outcome['source_row']}"
            )
        lfc = outcome["PA_vs_vehicle_positive_selection_LFC"]
        fdr = outcome["PA_vs_vehicle_positive_selection_FDR"]
        protective = lfc is not None and fdr is not None and fdr < 0.05 and lfc >= 0.25
        fdr_only = fdr is not None and fdr < 0.05
        interpretation = (
            "no_direct_result"
            if lfc is None
            else "continuous_negative_only_not_source_significant_sensitizing"
            if lfc < 0
            else "positive_selection_direction_continuous"
        )
        # Keep the full frozen universe row in memory for mechanical sensitivity
        # selection.  write_tsv later emits only the explicit SCREEN_FIELDS
        # outcome-release contract.
        row: dict[str, object] = {
            **frozen,
            **outcome,
            "expected_direction": directions.get(
                frozen["gene_symbol"], "not_applicable"
            ),
            "protective_hit_FDR_lt_0.05_LFC_ge_0.25": str(protective).upper(),
            "positive_selection_FDR_hit_without_LFC_floor": str(fdr_only).upper(),
            "source_significance_interpretation": interpretation,
            "analysis_guide_bin": "",
            "analysis_expression_quintile": "",
            "analysis_chronos_quintile": "",
            "actual_vehicle_fitness_observed": "FALSE",
            "outcome_values_read": "TRUE",
            "specification_sha256": SPECIFICATION_SHA256,
            "release_id": RELEASE_ID,
        }
        result.append(row)
    return result


def analysis_frame(rows: list[dict[str, object]], predicate) -> pd.DataFrame:
    selected = [row for row in rows if predicate(row)]
    records = []
    for row in selected:
        chronos_raw = row.get("HLF_Chronos", "")
        records.append(
            {
                "source_row": int(row["source_row"]),
                "gene_symbol": str(row["gene_symbol"]),
                "primary_evidence_class": str(row["primary_evidence_class"]),
                "guide_count_min": int(row["guide_count_min"]),
                "log1p_HLF_TPM": float(row["log1p_HLF_TPM"]),
                "HLF_Chronos": np.nan if chronos_raw == "" else float(chronos_raw),
                "HLF_Chronos_missing": int(chronos_raw == ""),
                "PA_vs_vehicle_LFC": float(row["PA_vs_vehicle_positive_selection_LFC"]),
                "PA_vs_vehicle_absolute_LFC": abs(
                    float(row["PA_vs_vehicle_positive_selection_LFC"])
                ),
                "D21_vs_D0_LFC": (
                    np.nan
                    if row["D21_vs_D0_positive_selection_LFC"] is None
                    else float(row["D21_vs_D0_positive_selection_LFC"])
                ),
                "protective_hit": parse_bool(
                    row["protective_hit_FDR_lt_0.05_LFC_ge_0.25"]
                ),
                "positive_selection_FDR_hit": parse_bool(
                    row["positive_selection_FDR_hit_without_LFC_floor"]
                ),
            }
        )
    frame = pd.DataFrame.from_records(records)
    if frame.empty or frame["gene_symbol"].duplicated().any():
        raise ValueError(
            "Analysis universe is empty or contains duplicate canonical symbols"
        )
    return frame.sort_values("gene_symbol", kind="stable").reset_index(drop=True)


def freeze_analysis_strata(
    primary: pd.DataFrame, expanded: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = add_frozen_strata(primary)
    base_bins = primary.set_index("gene_symbol")[
        [
            "guide_bin",
            "expression_quintile",
            "chronos_missing",
            "chronos_quintile",
            "permutation_stratum",
        ]
    ]
    expanded = expanded.copy()
    complete = expanded["HLF_Chronos_missing"].eq(0)
    if set(expanded.loc[complete, "gene_symbol"]) != set(primary["gene_symbol"]):
        raise ValueError(
            "Complete rows in missing-Chronos sensitivity differ from primary universe"
        )
    expanded["guide_bin"] = np.minimum(
        expanded["guide_count_min"].to_numpy(dtype=int), 8
    )
    expanded["expression_quintile"] = 0
    expanded["chronos_missing"] = expanded["HLF_Chronos_missing"].astype(int)
    expanded["chronos_quintile"] = 0
    for column in ("expression_quintile", "chronos_quintile"):
        expanded.loc[complete, column] = expanded.loc[complete, "gene_symbol"].map(
            base_bins[column]
        )
    missing = ~complete
    if missing.any():
        expanded.loc[missing, "expression_quintile"] = rank_first_bins(
            expanded.loc[missing, "log1p_HLF_TPM"].to_numpy(),
            expanded.loc[missing, "gene_symbol"].astype(str).to_numpy(),
        )
    median_chronos = float(primary["HLF_Chronos"].median())
    expanded.loc[missing, "HLF_Chronos"] = median_chronos
    expanded["permutation_stratum"] = [
        f"g{g}_e{e}_m{m}"
        for g, e, m in zip(
            expanded["guide_bin"],
            expanded["expression_quintile"],
            expanded["chronos_missing"],
        )
    ]
    expanded = add_frozen_strata(expanded)
    return primary, expanded


def run_class_analyses(
    primary: pd.DataFrame, expanded: pd.DataFrame, known_hits: set[str]
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    continuous_variants = [
        ("primary", primary, "PA_vs_vehicle_LFC", False),
        (
            "guide_floor_5",
            primary[primary["guide_count_min"] >= 5],
            "PA_vs_vehicle_LFC",
            False,
        ),
        (
            "remove_all_known_hits",
            primary[~primary["gene_symbol"].isin(known_hits)],
            "PA_vs_vehicle_LFC",
            False,
        ),
        (
            "remove_ACSL3",
            primary[primary["gene_symbol"] != "ACSL3"],
            "PA_vs_vehicle_LFC",
            False,
        ),
        ("missing_Chronos_indicator", expanded, "PA_vs_vehicle_LFC", True),
        (
            "PA_D21_vs_D0",
            primary[primary["D21_vs_D0_LFC"].notna()],
            "D21_vs_D0_LFC",
            False,
        ),
        ("absolute_PA_vs_vehicle_LFC", primary, "PA_vs_vehicle_absolute_LFC", False),
    ]
    class_rows: list[dict[str, object]] = []
    permutation_rows: list[dict[str, object]] = []
    sensitivity_rows: list[dict[str, object]] = []
    for variant, frame, outcome, include_missing in continuous_variants:
        rows, audit = fit_class_continuous(
            frame.copy(),
            outcome,
            variant,
            CLASS_SEED,
            accepted_draws=CLASS_DRAWS,
            include_missing_indicator=include_missing,
        )
        class_rows.extend(rows)
        permutation_rows.append(dict(audit, analysis_family="class_continuous"))
        sensitivity_rows.append(
            {
                "analysis_family": "class_continuous",
                "analysis_variant": variant,
                "outcome": outcome,
                "n_tested": len(frame),
                "accepted_draws": CLASS_DRAWS,
                "status": "completed",
            }
        )

    binary_variants = [
        ("primary_protective_hit", primary, "protective_hit", False),
        (
            "guide_floor_5",
            primary[primary["guide_count_min"] >= 5],
            "protective_hit",
            False,
        ),
        (
            "remove_all_known_hits",
            primary[~primary["gene_symbol"].isin(known_hits)],
            "protective_hit",
            False,
        ),
        (
            "remove_ACSL3",
            primary[primary["gene_symbol"] != "ACSL3"],
            "protective_hit",
            False,
        ),
        ("missing_Chronos_indicator", expanded, "protective_hit", True),
        (
            "positive_selection_FDR_without_LFC_floor",
            primary,
            "positive_selection_FDR_hit",
            False,
        ),
    ]
    binary_rows: list[dict[str, object]] = []
    for variant, frame, event, _ in binary_variants:
        rows, audits = fit_class_binary(
            frame.copy(), event, variant, BINARY_SEED, accepted_draws=BINARY_DRAWS
        )
        binary_rows.extend(rows)
        permutation_rows.extend(
            dict(row, analysis_family="class_binary") for row in audits
        )
        sensitivity_rows.append(
            {
                "analysis_family": "class_binary",
                "analysis_variant": variant,
                "outcome": event,
                "n_tested": len(frame),
                "accepted_draws": BINARY_DRAWS,
                "status": "completed",
            }
        )
    return class_rows, binary_rows, permutation_rows, sensitivity_rows


def run_program_analyses(
    primary: pd.DataFrame,
    prediction: list[dict[str, str]],
    known_hits: set[str],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    membership: dict[str, list[dict[str, object]]] = defaultdict(list)
    source_effects: dict[str, set[float]] = defaultdict(set)
    for row in prediction:
        if row["object_type"] != "hepatocyte_program_gene":
            continue
        membership[row["object_id"]].append(
            {
                "gene_symbol": row["gene_symbol"],
                "original_l1_weight": float(row["frozen_weight"]),
            }
        )
        source_effects[row["object_id"]].add(float(row["frozen_source_effect"]))
    if len(membership) != 2 or any(
        len(values) != 1 for values in source_effects.values()
    ):
        raise ValueError("Frozen external hepatocyte program family changed")

    variant_frames = {
        "primary": primary,
        "guide_floor_5": primary[primary["guide_count_min"] >= 5],
        "remove_all_known_hits": primary[~primary["gene_symbol"].isin(known_hits)],
        "remove_ACSL3": primary[primary["gene_symbol"] != "ACSL3"],
    }
    task_keys: list[tuple[str, str, str]] = []
    for uid in sorted(membership):
        for mode in ("original_weight", "equal_weight", "leave_highest_weight_out"):
            task_keys.append(("primary", mode, uid))
        for variant in ("guide_floor_5", "remove_all_known_hits", "remove_ACSL3"):
            task_keys.append((variant, "original_weight", uid))
    task_keys.sort()
    substreams = np.random.SeedSequence(PROGRAM_SEED).spawn(len(task_keys))
    program_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    sensitivity_rows: list[dict[str, object]] = []
    standardized: dict[str, pd.DataFrame] = {}
    for variant, frame in variant_frames.items():
        work = frame.copy()
        sd = float(work["PA_vs_vehicle_LFC"].std(ddof=1))
        if not math.isfinite(sd) or sd <= 0:
            raise ValueError(f"Nonpositive program outcome SD in {variant}")
        work["outcome_z"] = (
            work["PA_vs_vehicle_LFC"] - work["PA_vs_vehicle_LFC"].mean()
        ) / sd
        standardized[variant] = add_frozen_strata(work)
    for task_index, ((variant, mode, uid), seed_sequence) in enumerate(
        zip(task_keys, substreams)
    ):
        effect = next(iter(source_effects[uid]))
        if effect == 0:
            raise ValueError(f"Zero frozen direction for {uid}")
        members = membership[uid]
        result, audit = matched_program_test(
            standardized[variant],
            members,
            {str(row["gene_symbol"]) for row in members},
            1 if effect > 0 else -1,
            uid,
            variant,
            mode,
            seed_sequence,
            draws=PROGRAM_DRAWS,
        )
        result["rng_substream_index"] = task_index
        audit["rng_substream_index"] = task_index
        program_rows.append(result)
        audit_rows.append(audit)
    apply_program_bh(program_rows)
    seen_sensitivity: set[tuple[str, str]] = set()
    for variant, mode, _ in task_keys:
        key = (variant, mode)
        if key in seen_sensitivity:
            continue
        seen_sensitivity.add(key)
        rows = [
            row
            for row in program_rows
            if row["analysis_variant"] == variant and row["weight_mode"] == mode
        ]
        sensitivity_rows.append(
            {
                "analysis_family": "program_matched_null",
                "analysis_variant": variant,
                "outcome": mode,
                "n_tested": sum(row["testable"] == "TRUE" for row in rows),
                "accepted_draws": PROGRAM_DRAWS,
                "status": "completed_including_valid_untestable_rows",
            }
        )
    return program_rows, audit_rows, sensitivity_rows


def write_phase_c_outputs(
    outdir: Path,
    screen_rows: list[dict[str, object]],
    class_rows: list[dict[str, object]],
    binary_rows: list[dict[str, object]],
    program_rows: list[dict[str, object]],
    permutation_rows: list[dict[str, object]],
    matched_rows: list[dict[str, object]],
    sensitivity_rows: list[dict[str, object]],
    source_gate_pass: bool,
    code_hash: str,
) -> None:
    write_tsv(outdir / "screen_results.tsv", screen_rows, SCREEN_FIELDS)
    write_tsv(outdir / "class_effects.tsv", class_rows, list(class_rows[0]))
    write_tsv(outdir / "class_protective_hits.tsv", binary_rows, list(binary_rows[0]))
    permutation_fields: list[str] = []
    for row in permutation_rows:
        for field in row:
            if field not in permutation_fields:
                permutation_fields.append(field)
    write_tsv(outdir / "permutation_audit.tsv", permutation_rows, permutation_fields)
    write_tsv(outdir / "program_effects.tsv", program_rows, list(program_rows[0]))
    write_tsv(outdir / "matched_null_audit.tsv", matched_rows, list(matched_rows[0]))
    sensitivity_fields = list(sensitivity_rows[0]) + [
        "specification_sha256",
        "release_id",
    ]
    write_tsv(
        outdir / "sensitivity.tsv",
        [
            dict(row, specification_sha256=SPECIFICATION_SHA256, release_id=RELEASE_ID)
            for row in sensitivity_rows
        ],
        sensitivity_fields,
    )
    verdict = mechanical_figure_verdict(source_gate_pass, class_rows, program_rows)
    verdict.update(
        {
            "specification_sha256": SPECIFICATION_SHA256,
            "phase_c_code_bundle_sha256": code_hash,
            "biological_positivity_required_for_software_acceptance": "FALSE",
            "release_id": RELEASE_ID,
        }
    )
    write_tsv(outdir / "fig5_verdict.tsv", [verdict], list(verdict))

    release_paths = [
        outdir / "unseal_record.tsv",
        outdir / "UNSEALED",
        outdir / "phase_c_preflight_audit.tsv",
        outdir / "phase_c_code_manifest.tsv",
        outdir / "PHASE_C_CODE_READY",
        outdir / "screen_results.tsv",
        outdir / "class_effects.tsv",
        outdir / "class_protective_hits.tsv",
        outdir / "program_effects.tsv",
        outdir / "permutation_audit.tsv",
        outdir / "matched_null_audit.tsv",
        outdir / "sensitivity.tsv",
        outdir / "fig5_verdict.tsv",
    ]
    manifest_rows = []
    for path in release_paths:
        outcome_bearing = path.name not in {
            "phase_c_preflight_audit.tsv",
            "phase_c_code_manifest.tsv",
            "PHASE_C_CODE_READY",
        }
        manifest_rows.append(
            {
                "role": "phase_c_result"
                if path.name.endswith(".tsv")
                else "phase_c_state",
                "relative_path": str(path.relative_to(outdir)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "line_count": sum(1 for _ in path.open(encoding="utf-8")),
                "external_outcomes_read": str(outcome_bearing).upper(),
                "canonical": "FALSE",
                "release_id": RELEASE_ID,
            }
        )
    write_tsv(
        outdir / "phase_c_release_manifest.tsv",
        manifest_rows,
        [
            "role",
            "relative_path",
            "bytes",
            "sha256",
            "line_count",
            "external_outcomes_read",
            "canonical",
            "release_id",
        ],
    )
    atomic_write_text(
        outdir / "PHASE_C_COMPLETE",
        "release_id\tstatus\tspecification_sha256\tphase_c_code_bundle_sha256\t"
        "release_manifest_sha256\texternal_outcomes_read\tbiological_positivity_required\n"
        f"{RELEASE_ID}\tone_pass_execution_complete\t{SPECIFICATION_SHA256}\t{code_hash}\t"
        f"{sha256_file(outdir / 'phase_c_release_manifest.tsv')}\tTRUE\tFALSE\n",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=CANDIDATE_ROOT)
    args = parser.parse_args()
    outdir = require_within(args.outdir, CANDIDATE_ROOT)
    unsealed = False
    try:
        audit, workbook, code_hash = preflight(outdir)
        write_tsv(
            outdir / "phase_c_preflight_audit.tsv",
            audit,
            ["check", "passed", "detail", "external_outcomes_read", "release_id"],
        )
    except Exception as error:
        write_tsv(
            outdir / "phase_c_terminal_blocker.tsv",
            [
                {
                    "blocked_at_utc": utc_now(),
                    "stage": "pre_unseal",
                    "reason": str(error),
                    "specification_sha256": SPECIFICATION_SHA256,
                    "external_outcomes_read": "FALSE",
                    "terminal": "TRUE",
                    "release_id": RELEASE_ID,
                }
            ],
            [
                "blocked_at_utc",
                "stage",
                "reason",
                "specification_sha256",
                "external_outcomes_read",
                "terminal",
                "release_id",
            ],
        )
        raise

    unsealed_at = utc_now()
    executor = str(Path(__file__).resolve().relative_to(PROJECT_ROOT))
    write_tsv(
        outdir / "unseal_record.tsv",
        [
            {
                "unsealed_at_utc": unsealed_at,
                "specification_sha256": SPECIFICATION_SHA256,
                "executor": executor,
                "raw_outcome_source_sha256": RAW_SHA256,
                "one_pass_release_id": RELEASE_ID,
            }
        ],
        [
            "unsealed_at_utc",
            "specification_sha256",
            "executor",
            "raw_outcome_source_sha256",
            "one_pass_release_id",
        ],
    )
    atomic_write_text(
        outdir / "UNSEALED",
        "release_id\tstatus\tunsealed_at_utc\tspecification_sha256\t"
        "phase_c_code_bundle_sha256\traw_outcome_source_sha256\texternal_outcomes_read\n"
        f"{RELEASE_ID}\tone_pass_outcome_access_started\t{unsealed_at}\t"
        f"{SPECIFICATION_SHA256}\t{code_hash}\t{RAW_SHA256}\tTRUE\n",
    )
    unsealed = True
    try:
        outcomes = parse_outcomes_once(workbook)
        universe = read_tsv(outdir / "screen_gene_universe.tsv")
        prediction = read_tsv(outdir / "prediction_manifest.tsv")
        screen_rows = build_screen_results(universe, prediction, outcomes)
        primary = analysis_frame(
            screen_rows, lambda row: parse_bool(row["primary_eligible"])
        )
        if len(primary) != 7_762:
            raise ValueError(f"Frozen primary universe changed: {len(primary)} != 7762")
        expanded = analysis_frame(
            screen_rows,
            lambda row: (
                parse_bool(row["neutral_source_base"])
                and not parse_bool(row["depmap_entrez_conflict"])
                and parse_bool(row["HLF_expression_floor_pass"])
                and row["log1p_HLF_TPM"] != ""
            ),
        )
        primary, expanded = freeze_analysis_strata(primary, expanded)
        bins = primary.set_index("source_row")
        for row in screen_rows:
            source_row = int(row["source_row"])
            if source_row in bins.index:
                row["analysis_guide_bin"] = int(bins.loc[source_row, "guide_bin"])
                row["analysis_expression_quintile"] = int(
                    bins.loc[source_row, "expression_quintile"]
                )
                row["analysis_chronos_quintile"] = int(
                    bins.loc[source_row, "chronos_quintile"]
                )
        known_hits = {
            row["canonical_symbol"]
            for row in read_tsv(outdir / "known_hit_exclusion.tsv")
        }
        class_rows, binary_rows, permutation_rows, class_sensitivity = (
            run_class_analyses(primary, expanded, known_hits)
        )
        program_rows, matched_rows, program_sensitivity = run_program_analyses(
            primary, prediction, known_hits
        )
        write_phase_c_outputs(
            outdir,
            screen_rows,
            class_rows,
            binary_rows,
            program_rows,
            permutation_rows,
            matched_rows,
            class_sensitivity + program_sensitivity,
            source_gate_pass=True,
            code_hash=code_hash,
        )
        print(
            "PHASE_C_ONE_PASS_COMPLETE "
            f"specification_sha256={SPECIFICATION_SHA256} tested_universe={len(primary)} "
            "biological_positivity_not_required=TRUE"
        )
    except Exception as error:
        if unsealed:
            failure_type = type(error).__name__
            write_tsv(
                outdir / "phase_c_terminal_failure.tsv",
                [
                    {
                        "failed_at_utc": utc_now(),
                        "stage": "post_unseal_no_rerun",
                        "error_type": failure_type,
                        "reason": str(error),
                        "traceback_sha256": __import__("hashlib")
                        .sha256(traceback.format_exc().encode("utf-8"))
                        .hexdigest(),
                        "specification_sha256": SPECIFICATION_SHA256,
                        "external_outcomes_read": "TRUE",
                        "terminal": "TRUE",
                        "release_id": RELEASE_ID,
                    }
                ],
                [
                    "failed_at_utc",
                    "stage",
                    "error_type",
                    "reason",
                    "traceback_sha256",
                    "specification_sha256",
                    "external_outcomes_read",
                    "terminal",
                    "release_id",
                ],
            )
        raise


if __name__ == "__main__":
    main()
