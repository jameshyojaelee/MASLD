#!/usr/bin/env python3
"""Pure helpers for the SP-INT-06/07 candidate-only final integration.

The integration is deliberately a union of assay-native rows.  It never
standardizes, sums, ranks, or otherwise makes effects from different assays
commensurate.
"""

from __future__ import annotations

import csv
import hashlib
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

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


HOTFIX_EXPECTED_SHA256 = "22e34649838bfea2ce36cf324559b088dfc8a1fed5ca7475a91889b188a5e09a"

INTEGRATION_EXTRA_COLUMNS = (
    "source_adapter_id",
    "evidence_role",
    "claim_scope",
    "uncertainty_semantics",
    "import_method",
    "source_ready_path",
    "source_ready_sha256",
    "source_effect_path",
    "source_effect_sha256",
    "source_effect_row_sha256",
    "figure4_include",
    "figure4_interpretation",
)
INTEGRATED_EFFECT_COLUMNS = SCHEMAS["program_effects.tsv"] + INTEGRATION_EXTRA_COLUMNS

FIGURE4_MATRIX_COLUMNS = (
    "release_id",
    "dataset",
    "assay",
    "analysis_set_id",
    "program_uid",
    "legacy_program_id",
    "program_label",
    "source_dependence",
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
    "n_genes_measured",
    "retained_l1_weight",
    "testable",
    "direction_observed",
    "descriptive_effect_direction",
    "inferential_test_direction",
    "direction_agreement",
    "heterogeneity_statistic",
    "heterogeneity_df",
    "heterogeneity_pvalue",
    "sensitivity_sign_agree",
    "robustness_pass",
    "evidence_state",
    "negative_call_rule_id",
    "evidence_role",
    "claim_scope",
    "uncertainty_semantics",
    "figure4_interpretation",
)

VERDICT_COLUMNS = (
    "release_id",
    "dataset",
    "assay",
    "analysis_set_id",
    "source_dependence",
    "biological_unit",
    "n_biological",
    "technical_unit",
    "n_technical",
    "effect_unit",
    "n_programs",
    "n_robust",
    "n_indeterminate",
    "n_tested_negative",
    "n_untestable",
    "n_not_applicable",
    "n_skipped",
    "dataset_gate",
    "figure4_verdict",
    "figure4_role",
    "interpretation",
)

INPUT_SEAL_COLUMNS = (
    "release_id",
    "relative_path",
    "bytes",
    "sha256",
    "source_role",
    "seal_state",
)

SOURCE_TABLE_MANIFEST_COLUMNS = (
    "release_id",
    "table_id",
    "relative_path",
    "data_rows",
    "bytes",
    "sha256",
    "figure4_role",
    "interpretation",
)

FORBIDDEN_INTEGRATION_TOKENS = (
    "universal_score",
    "combined_score",
    "cross_assay_score",
    "cross_assay_rank",
    "overall_score",
    "overall_rank",
    "modality_count",
    "evidence_count",
)


def project_relative(project_root: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(project_root.resolve()).as_posix()
    except ValueError as exc:
        raise ContractError(f"path is outside project root: {path}") from exc


def canonical_row_sha256(row: Mapping[str, object], columns: Sequence[str]) -> str:
    payload = "\t".join(str(row.get(column, "")) for column in columns) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def bool_text(value: object) -> str:
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return "TRUE" if parse_bool(value, "boolean") else "FALSE"


def direction(value: float) -> str:
    return "positive" if value > 0 else "negative" if value < 0 else "zero"


def independent_bh(values: Mapping[str, float]) -> dict[str, float]:
    ordered = sorted(values.items(), key=lambda item: (item[1], item[0]))
    adjusted: dict[str, float] = {}
    running = 1.0
    total = len(ordered)
    for index in range(total - 1, -1, -1):
        key, value = ordered[index]
        running = min(running, value * total / (index + 1), 1.0)
        adjusted[key] = running
    return adjusted


def expected_programs(hotspot_root: Path) -> dict[str, dict[str, str]]:
    _, rows = read_tsv(
        hotspot_root / "program_registry_v2.tsv",
        (
            "cell_type",
            "module",
            "program_uid",
            "membership_sha256",
            "module_name",
            "primary_direction",
            "external_test_eligible",
        ),
    )
    selected = {
        row["program_uid"]: row
        for row in rows
        if parse_bool(row["external_test_eligible"], f"external[{row['program_uid']}]")
    }
    if len(selected) != 2 or any(row["cell_type"] != "hepatocytes" for row in selected.values()):
        raise ContractError(f"final integration expects two frozen hepatocyte programs, found {len(selected)}")
    return selected


def require_complete_groups(
    rows: Sequence[Mapping[str, str]],
    expected_uids: Iterable[str],
    group_columns: Sequence[str] = ("dataset", "assay", "analysis_set_id"),
) -> None:
    expected = set(expected_uids)
    groups: dict[tuple[str, ...], list[str]] = defaultdict(list)
    for row in rows:
        groups[tuple(row[column] for column in group_columns)].append(row["program_uid"])
    if not groups:
        raise ContractError("integrated effect table has no dataset/assay groups")
    for key, observed in groups.items():
        if len(observed) != len(set(observed)) or set(observed) != expected:
            raise ContractError(
                f"incomplete/duplicate frozen program family for {key}: observed={sorted(observed)}, expected={sorted(expected)}"
            )


def require_no_cross_assay_construct(rows: Sequence[Mapping[str, str]], columns: Sequence[str]) -> None:
    lower_columns = {column.lower() for column in columns}
    for token in FORBIDDEN_INTEGRATION_TOKENS:
        if token in lower_columns:
            raise ContractError(f"prohibited cross-assay integration column: {token}")
    for row in rows:
        comparable_value = row["cross_assay_comparable"]
        comparable = (
            comparable_value
            if isinstance(comparable_value, bool)
            else parse_bool(str(comparable_value), f"cross_assay[{row['dataset']}]")
        )
        if comparable:
            raise ContractError("integrated row claims cross-assay comparability")
        if row["effect_unit"] in {"universal_score", "combined_score", "standardized_cross_assay_score"}:
            raise ContractError(f"prohibited integrated effect unit: {row['effect_unit']}")


def adapter_row(
    row: Mapping[str, str],
    *,
    adapter_id: str,
    evidence_role: str,
    claim_scope: str,
    uncertainty_semantics: str,
    ready_path: str,
    ready_sha256: str,
    effect_path: str,
    effect_sha256: str,
    interpretation: str,
) -> dict[str, object]:
    output: dict[str, object] = {column: row[column] for column in SCHEMAS["program_effects.tsv"]}
    unresolved_array = (
        row["dataset"] == "Govaere2026_CosMx"
        and row["biological_unit_resolution"] == "unresolved"
        and row["technical_unit"] == "physical_array"
    )
    if unresolved_array:
        output["biological_unit"] = "unknown_public_biological_unit"
        output["n_biological"] = ""
    output.update(
        {
            "source_adapter_id": adapter_id,
            "evidence_role": evidence_role,
            "claim_scope": claim_scope,
            "uncertainty_semantics": uncertainty_semantics,
            "import_method": (
                "adjudicated_unresolved_array_not_biological_replicate"
                if unresolved_array
                else "exact_adapter_effect_row"
            ),
            "source_ready_path": ready_path,
            "source_ready_sha256": ready_sha256,
            "source_effect_path": effect_path,
            "source_effect_sha256": effect_sha256,
            "source_effect_row_sha256": canonical_row_sha256(row, SCHEMAS["program_effects.tsv"]),
            "figure4_include": True,
            "figure4_interpretation": interpretation,
        }
    )
    return output


def native_spatial_rows(
    project_root: Path,
    native_root: Path,
    hotspot_root: Path,
) -> list[dict[str, object]]:
    selected = expected_programs(hotspot_root)
    result_path = native_root / "spatial_program_results.tsv"
    universe_path = native_root / "tested_universe.tsv"
    design_path = native_root / "native_design_audit.tsv"
    source_manifest_path = native_root / "source_manifest.tsv"
    execution_manifest_path = native_root / "execution_manifest.tsv"
    _, results = read_tsv(
        result_path,
        (
            "program_id",
            "dataset",
            "program_name",
            "n_measured",
            "retained_l1_weight",
            "testable",
            "residual_moran_i",
            "residual_null_mean",
            "residual_null_sd",
            "residual_moran_z",
            "residual_pvalue",
            "residual_padj",
            "sensitivity_sign_agree",
            "robust",
        ),
    )
    _, universe_rows = read_tsv(
        universe_path,
        (
            "program_id",
            "legacy_program_id",
            "cell_type",
            "module",
            "program_name",
            "membership_sha256",
            "primary_direction",
        ),
    )
    universe = {row["program_id"]: row for row in universe_rows}
    if set(universe) != set(selected) or len(universe_rows) != 2:
        raise ContractError("native v2 tested universe is not the sealed two-program family")
    _, designs = read_tsv(
        design_path,
        (
            "dataset",
            "biological_unit",
            "n_biological",
            "technical_unit",
            "n_technical",
            "design_status",
        ),
    )
    design = {row["dataset"]: row for row in designs}
    expected_datasets = {"GSE192741", "Vu_et_al_2025"}
    if set(design) != expected_datasets or len(designs) != 2:
        raise ContractError("native v2 design audit must cover GSE192741 and Vu_et_al_2025")
    if len(results) != 4 or {(row["dataset"], row["program_id"]) for row in results} != {
        (dataset, uid) for dataset in expected_datasets for uid in selected
    }:
        raise ContractError("native v2 spatial program table is not a complete 2 x 2 grid")

    _, native_sources = read_tsv(
        source_manifest_path,
        ("relative_path", "source_role", "bytes", "sha256"),
    )
    candidate_producers = [
        row for row in native_sources if row["source_role"] == "candidate_producer"
    ]
    if len(candidate_producers) != 1:
        raise ContractError("native spatial source manifest must identify one candidate producer")
    producer_row = candidate_producers[0]
    producer = project_root / producer_row["relative_path"]
    if not producer.is_file() or sha256_file(producer) != producer_row["sha256"]:
        raise ContractError("native spatial candidate-producer provenance drift")
    _, execution_rows = read_tsv(execution_manifest_path, ("parameter", "value"))
    execution = {row["parameter"]: row["value"] for row in execution_rows}
    if execution.get("candidate_producer_sha256") != producer_row["sha256"]:
        raise ContractError("native spatial execution/source producer hashes disagree")

    source_ready = native_root / "READY"
    ready_sha = sha256_file(source_ready)
    source_effect_sha = sha256_file(result_path)
    source_manifest_sha = sha256_file(source_manifest_path)
    producer_sha = sha256_file(producer)
    output = []
    for row in results:
        uid = row["program_id"]
        reg = selected[uid]
        univ = universe[uid]
        ds = row["dataset"]
        dsg = design[ds]
        testable = parse_bool(row["testable"], f"native.testable[{ds}/{uid}]")
        robust = parse_bool(row["robust"], f"native.robust[{ds}/{uid}]")
        sign_agree = parse_bool(row["sensitivity_sign_agree"], f"native.sign[{ds}/{uid}]")
        n_measured = parse_int(row["n_measured"], f"native.n[{ds}/{uid}]")
        retained = parse_float(row["retained_l1_weight"], f"native.weight[{ds}/{uid}]")
        if n_measured is None or retained is None:
            raise ContractError("native spatial testability fields are missing")
        common = {
            "release_id": RELEASE_ID,
            "registry_sha256": sha256_file(hotspot_root / "program_registry_v2.tsv"),
            "membership_sha256": reg["membership_sha256"],
            "program_uid": uid,
            "legacy_program_id": f"{reg['cell_type']}:{reg['module']}",
            "program_label": reg["module_name"],
            "cell_type": reg["cell_type"],
            "dataset": ds,
            "assay": "Visium_spatial_transcriptomics",
            "source_publication": (
                "GSE192741_NCBI_GEO_public_record"
                if ds == "GSE192741"
                else "Vu_et_al_2025_public_spatial_liver_dataset"
            ),
            "source_dependence": "independent" if ds == "GSE192741" else "source_dependent",
            "analysis_set_id": "matched_null_residual_spatial_autocorrelation",
            "biological_unit": (
                dsg["biological_unit"] if ds == "GSE192741" else "unknown_public_biological_unit"
            ),
            "biological_unit_resolution": "resolved" if ds == "GSE192741" else "unresolved",
            "n_biological": dsg["n_biological"] if ds == "GSE192741" else "",
            "technical_unit": dsg["technical_unit"],
            "n_technical": dsg["n_technical"],
            "contrast_or_exposure": "residual_spatial_autocorrelation_against_matched_gene_null",
            "effect_unit": "centered_residual_moran_i_against_matched_gene_null",
            "n_genes_measured": n_measured,
            "retained_l1_weight": retained,
            "testable": testable,
            "testability_reason": (
                "passed_8_gene_20pct_L1_and_matched_null_spatial_gate"
                if testable
                else "failed_native_spatial_testability_gate"
            ),
            "direction_expected": reg["primary_direction"],
            "cross_assay_comparable": False,
            "producer": project_relative(project_root, producer),
            "producer_sha256": producer_sha,
            "source_manifest_sha256": source_manifest_sha,
        }
        if testable:
            residual = parse_float(row["residual_moran_i"], f"native.residual[{ds}/{uid}]")
            null_mean = parse_float(row["residual_null_mean"], f"native.null_mean[{ds}/{uid}]")
            null_sd = parse_float(row["residual_null_sd"], f"native.null_sd[{ds}/{uid}]")
            pvalue = parse_float(row["residual_pvalue"], f"native.p[{ds}/{uid}]")
            padj = parse_float(row["residual_padj"], f"native.q[{ds}/{uid}]")
            if None in (residual, null_mean, null_sd, pvalue, padj) or (null_sd or 0.0) <= 0:
                raise ContractError("testable native spatial row lacks a valid matched-null effect")
            estimate = float(residual) - float(null_mean)
            common.update(
                {
                    "estimate": estimate,
                    "std_error": "",
                    "matched_null_sd": null_sd,
                    "interval_low": "",
                    "interval_high": "",
                    "interval_type": "none",
                    "pvalue": pvalue,
                    "padj": padj,
                    "pvalue_method": "empirical_matched_gene_null_two_sided_9999_draws",
                    "multiplicity_family": "complete_two_program_family_within_dataset",
                    "direction_observed": direction(estimate),
                    "descriptive_effect_direction": direction(estimate),
                    "inferential_test_direction": direction(estimate),
                    "direction_agreement": True,
                    "heterogeneity_statistic": "",
                    "heterogeneity_df": "",
                    "heterogeneity_pvalue": "",
                    "sensitivity_sign_agree": sign_agree,
                    "robustness_pass": robust,
                    "negative_call_rule_id": "",
                    "evidence_state": "robust" if robust else "indeterminate",
                }
            )
        else:
            common.update(
                {
                    "estimate": "",
                    "std_error": "",
                    "matched_null_sd": "",
                    "interval_low": "",
                    "interval_high": "",
                    "interval_type": "not_applicable",
                    "pvalue": "",
                    "padj": "",
                    "pvalue_method": "not_applicable",
                    "multiplicity_family": "complete_two_program_family_within_dataset",
                    "direction_observed": "",
                    "descriptive_effect_direction": "",
                    "inferential_test_direction": "",
                    "direction_agreement": False,
                    "heterogeneity_statistic": "",
                    "heterogeneity_df": "",
                    "heterogeneity_pvalue": "",
                    "sensitivity_sign_agree": False,
                    "robustness_pass": False,
                    "negative_call_rule_id": "",
                    "evidence_state": "untestable",
                }
            )
        native_row = {column: common[column] for column in SCHEMAS["program_effects.tsv"]}
        native_row.update(
            {
                "source_adapter_id": f"{ds.lower()}_visium_native_v2",
                "evidence_role": "spatial_organization",
                "claim_scope": (
                    "donor_supported_spatial_organization"
                    if ds == "GSE192741"
                    else "physical_array_spatial_pattern_source_dependent"
                ),
                "uncertainty_semantics": "matched_gene_null_standard_deviation_not_sampling_standard_error",
                "import_method": "deterministic_native_spatial_adapter",
                "source_ready_path": project_relative(project_root, source_ready),
                "source_ready_sha256": ready_sha,
                "source_effect_path": project_relative(project_root, result_path),
                "source_effect_sha256": source_effect_sha,
                "source_effect_row_sha256": canonical_row_sha256(row, tuple(row)),
                "figure4_include": True,
                "figure4_interpretation": (
                    "retain_complete_donor_supported_spatial_panel"
                    if ds == "GSE192741"
                    else "retain_complete_array_level_source_dependent_spatial_panel"
                ),
            }
        )
        output.append(native_row)
    return output


def make_verdicts(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, str], list[Mapping[str, object]]] = defaultdict(list)
    for row in rows:
        groups[(str(row["dataset"]), str(row["assay"]), str(row["analysis_set_id"]))].append(row)
    verdicts = []
    for key in sorted(groups):
        part = groups[key]
        states = Counter(str(row["evidence_state"]) for row in part)
        if states["robust"]:
            gate = "pass"
        elif states["tested_negative"]:
            gate = "valid_null"
        elif states["indeterminate"]:
            gate = "indeterminate"
        elif states["untestable"] == len(part):
            gate = "untestable"
        elif states["skipped"] == len(part):
            gate = "skipped"
        elif states["not_applicable"] == len(part):
            gate = "not_applicable"
        else:
            gate = "mixed_terminal"
        dataset = key[0]
        if dataset == "Govaere2026_CosMx":
            verdict = "retain_source_native_context_only"
            role = "cell_proximity_context"
            interpretation = "whole-program inference not applicable; retain IL32 array-level source-native panel"
        elif dataset == "Govaere2026_GeoMx":
            verdict = "retain_exclusion_boundary"
            role = "transparent_assay_exclusion"
            interpretation = "local GeoMx reanalysis remains dropped; no numeric program effect"
        elif dataset == "PXD051911":
            verdict = "retain_fixed25_descriptive_and_program_untestable_rows"
            role = "protein_observability_boundary"
            interpretation = "two frozen programs fail DIA observability; fixed 25 remains selection-conditioned descriptive context"
        elif gate == "skipped":
            verdict = "retain_transparent_skip_rows"
            role = "source_gate_boundary"
            interpretation = "preserve both preselected program rows without inventing inference"
        else:
            verdict = "retain_complete_program_panel"
            role = str(part[0]["evidence_role"])
            interpretation = "show both preselected programs regardless of significance"
        first = part[0]
        verdicts.append(
            {
                "release_id": RELEASE_ID,
                "dataset": dataset,
                "assay": key[1],
                "analysis_set_id": key[2],
                "source_dependence": first["source_dependence"],
                "biological_unit": first["biological_unit"],
                "n_biological": first["n_biological"],
                "technical_unit": first["technical_unit"],
                "n_technical": first["n_technical"],
                "effect_unit": first["effect_unit"],
                "n_programs": len(part),
                "n_robust": states["robust"],
                "n_indeterminate": states["indeterminate"],
                "n_tested_negative": states["tested_negative"],
                "n_untestable": states["untestable"],
                "n_not_applicable": states["not_applicable"],
                "n_skipped": states["skipped"],
                "dataset_gate": gate,
                "figure4_verdict": verdict,
                "figure4_role": role,
                "interpretation": interpretation,
            }
        )
    return verdicts


def figure4_matrix(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    return [{column: row[column] for column in FIGURE4_MATRIX_COLUMNS} for row in rows]


def count_data_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def seal_rows(project_root: Path, paths: Sequence[tuple[Path, str, str]]) -> list[dict[str, object]]:
    result = []
    seen: set[Path] = set()
    for path, role, state in paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if not path.is_file():
            raise ContractError(f"required final-integration input is missing: {path}")
        result.append(
            {
                "release_id": RELEASE_ID,
                "relative_path": project_relative(project_root, path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "source_role": role,
                "seal_state": state,
            }
        )
    return result


def verify_hotfix(path: Path) -> dict[str, str]:
    if sha256_file(path) != HOTFIX_EXPECTED_SHA256:
        raise ContractError("Yakubovsky validator-hotfix provenance hash drift")
    _, rows = read_tsv(
        path,
        (
            "status",
            "repair_scope",
            "model_or_result_recomputed",
            "acceptance_rule_changed",
            "post_outcome_validation_only",
        ),
    )
    if len(rows) != 1:
        raise ContractError("Yakubovsky validator hotfix manifest must contain one row")
    row = rows[0]
    if row["status"] != "validation_only_iterator_hotfix_applied":
        raise ContractError("Yakubovsky hotfix status drift")
    if row["repair_scope"] != "source_abundance_fixture_first_two_references_against_first_two_certified_axis_rows":
        raise ContractError("Yakubovsky hotfix repair scope drift")
    if parse_bool(row["model_or_result_recomputed"], "hotfix.model"):
        raise ContractError("Yakubovsky hotfix unexpectedly recomputed model/results")
    if parse_bool(row["acceptance_rule_changed"], "hotfix.acceptance"):
        raise ContractError("Yakubovsky hotfix unexpectedly changed an acceptance rule")
    if not parse_bool(row["post_outcome_validation_only"], "hotfix.post_outcome"):
        raise ContractError("Yakubovsky hotfix lost its post-outcome validation-only disclosure")
    return row


def write_effect_outputs(root: Path, rows: Sequence[Mapping[str, object]]) -> None:
    ordered = sorted(
        rows,
        key=lambda row: (str(row["dataset"]), str(row["assay"]), str(row["analysis_set_id"]), str(row["program_uid"])),
    )
    write_tsv(root / "integrated_program_effects.tsv", INTEGRATED_EFFECT_COLUMNS, ordered)
    write_tsv(root / "figure4_program_matrix.tsv", FIGURE4_MATRIX_COLUMNS, figure4_matrix(ordered))
    write_tsv(root / "figure4_dataset_verdict.tsv", VERDICT_COLUMNS, make_verdicts(ordered))


def finite_or_empty(value: str) -> bool:
    if not value.strip():
        return True
    try:
        return math.isfinite(float(value))
    except ValueError:
        return False
