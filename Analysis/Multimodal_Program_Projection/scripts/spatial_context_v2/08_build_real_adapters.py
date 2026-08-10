#!/usr/bin/env python3
"""Build the bounded real-adapter seam for Plan 13.

The producer performs no new biological inference.  It closes GSE287826 as a
donor-key skip, propagates the repository's explicit Govaere GeoMx exclusion,
and preserves the accepted Govaere CosMx IL32 proximity result at its native
physical-array unit.  The CosMx source-native rows are not converted into
whole-program effects; the required v2 program rows are therefore explicitly
``not_applicable`` while assay coverage is audited.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from contract_lib import (
    ADAPTER_REGISTRY_COLUMNS,
    ContractError,
    RELEASE_ID,
    SCHEMAS,
    canonical_id_hash,
    parse_bool,
    program_universe,
    read_tsv,
    sha256_file,
    validate_hotspot_seal,
    write_tsv,
)


SOURCE_NATIVE_COLUMNS = (
    "release_id",
    "dataset",
    "assay",
    "source_publication",
    "source_dependence",
    "reporting_unit",
    "reporting_id",
    "is_mash",
    "include_direction_concordance",
    "n_hepatocytes",
    "n_macrophages",
    "n_metmac",
    "il32_partition",
    "contrast_or_exposure",
    "effect_unit",
    "estimate",
    "secondary_effect_unit",
    "secondary_estimate",
    "direction_concordant",
    "inferential_status",
    "interpretation",
)


GOVAERE_CITATION = (
    "Boesch_et_al_Nature_Genetics_2026_"
    "doi_10.1038/s41588-026-02600-3"
)


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _project_relative(project_root: Path, path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(project_root.resolve()).as_posix()
    except ValueError as exc:
        raise ContractError(f"source path is outside the project root: {path}") from exc


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ContractError(f"missing source CSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ContractError(f"missing CSV header: {path}")
        return [{key: (value if value is not None else "") for key, value in row.items()} for row in reader]


def _source_manifest(
    project_root: Path,
    adapter_root: Path,
    source_rows: Sequence[tuple[Path, str, bool, str]],
) -> str:
    rows = []
    for path, role, public_access, retrieved_utc in source_rows:
        if not path.is_file():
            raise ContractError(f"missing source-manifest input: {path}")
        rows.append(
            {
                "path_scope": "project_relative",
                "relative_path": _project_relative(project_root, path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "source_role": role,
                "public_access": public_access,
                "retrieved_utc": retrieved_utc,
            }
        )
    write_tsv(adapter_root / "source_manifest.tsv", SCHEMAS["source_manifest.tsv"], rows)
    return sha256_file(adapter_root / "source_manifest.tsv")


def _execution_manifest(
    adapter_root: Path,
    filenames: Iterable[str],
    created_utc: str,
) -> str:
    rows = []
    for filename in filenames:
        path = adapter_root / filename
        if not path.is_file():
            raise ContractError(f"cannot seal missing adapter artifact: {path}")
        rows.append(
            {
                "role": "real_adapter_artifact",
                "relative_path": filename,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "created_utc": created_utc,
            }
        )
    write_tsv(
        adapter_root / "execution_manifest.tsv",
        SCHEMAS["execution_manifest.tsv"],
        rows,
    )
    return sha256_file(adapter_root / "execution_manifest.tsv")


def _base_program_rows(
    spec: Mapping[str, object],
    seal,
    source_manifest_sha256: str,
    producer_relative: str,
    producer_sha256: str,
    coverage: Mapping[str, tuple[int, float, str]],
    n_biological: int | None,
    n_technical: int,
    state: str,
    reason: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    universe = program_universe(seal, str(spec["program_universe"]))
    mapping_rows = []
    testability_rows = []
    effect_rows = []
    sensitivity_rows = []
    for uid in sorted(universe):
        reg = universe[uid]
        n_genes, retained_weight, mapping_status = coverage.get(
            uid,
            (0, 0.0, "not_measured_or_not_opened"),
        )
        mapping_rows.append(
            {
                "release_id": RELEASE_ID,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "program_uid": uid,
                "membership_sha256": reg["membership_sha256"],
                "n_source_genes": reg["n_source_genes"],
                "n_genes_measured": n_genes,
                "retained_l1_weight": retained_weight,
                "mapping_status": mapping_status,
            }
        )
        testability_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_sha256": seal.registry_sha256,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "program_uid": uid,
                "membership_sha256": reg["membership_sha256"],
                "testable": False,
                "n_genes_measured": n_genes,
                "retained_l1_weight": retained_weight,
                "testability_reason": reason,
                "evidence_state": state,
            }
        )
        effect_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_sha256": seal.registry_sha256,
                "membership_sha256": reg["membership_sha256"],
                "program_uid": uid,
                "legacy_program_id": f"{reg['cell_type']}:{reg['module']}",
                "program_label": reg["module_name"],
                "cell_type": reg["cell_type"],
                "dataset": spec["dataset"],
                "assay": spec["assay"],
                "source_publication": spec["source_publication"],
                "source_dependence": spec["source_dependence"],
                "analysis_set_id": spec["analysis_set_id"],
                "biological_unit": spec["biological_unit"],
                "biological_unit_resolution": spec["biological_unit_resolution"],
                "n_biological": n_biological,
                "technical_unit": spec["technical_unit"],
                "n_technical": n_technical,
                "contrast_or_exposure": spec["contrast_or_exposure"],
                "effect_unit": spec["effect_unit"],
                "estimate": "",
                "std_error": "",
                "matched_null_sd": "",
                "interval_low": "",
                "interval_high": "",
                "interval_type": "not_applicable",
                "pvalue": "",
                "padj": "",
                "pvalue_method": "not_applicable",
                "multiplicity_family": "not_tested",
                "n_genes_measured": n_genes,
                "retained_l1_weight": retained_weight,
                "testable": False,
                "testability_reason": reason,
                "direction_expected": reg.get("primary_direction", "not_prespecified"),
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
                "cross_assay_comparable": False,
                "evidence_state": state,
                "producer": producer_relative,
                "producer_sha256": producer_sha256,
                "source_manifest_sha256": source_manifest_sha256,
            }
        )
        sensitivity_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_sha256": seal.registry_sha256,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "program_uid": uid,
                "membership_sha256": reg["membership_sha256"],
                "sensitivity_id": "not_run_no_source_native_program_estimand",
                "estimate": "",
                "effect_unit": spec["effect_unit"],
                "direction_observed": "",
                "sign_agree": False,
                "status": state,
            }
        )
    return mapping_rows, testability_rows, effect_rows, sensitivity_rows


def _write_adapter_tables(
    adapter_root: Path,
    sample_rows: Sequence[Mapping[str, object]],
    design_row: Mapping[str, object],
    mapping_rows: Sequence[Mapping[str, object]],
    testability_rows: Sequence[Mapping[str, object]],
    effect_rows: Sequence[Mapping[str, object]],
    sensitivity_rows: Sequence[Mapping[str, object]],
) -> None:
    write_tsv(adapter_root / "sample_manifest.tsv", SCHEMAS["sample_manifest.tsv"], sample_rows)
    write_tsv(adapter_root / "design_audit.tsv", SCHEMAS["design_audit.tsv"], [design_row])
    write_tsv(adapter_root / "gene_mapping_audit.tsv", SCHEMAS["gene_mapping_audit.tsv"], mapping_rows)
    write_tsv(adapter_root / "program_testability.tsv", SCHEMAS["program_testability.tsv"], testability_rows)
    write_tsv(adapter_root / "per_sample_program_scores.tsv", SCHEMAS["per_sample_program_scores.tsv"], [])
    write_tsv(adapter_root / "program_effects.tsv", SCHEMAS["program_effects.tsv"], effect_rows)
    write_tsv(adapter_root / "sensitivity.tsv", SCHEMAS["sensitivity.tsv"], sensitivity_rows)


def _empty_design(spec: Mapping[str, object], status: str) -> dict[str, object]:
    unresolved = spec["biological_unit_resolution"] == "unresolved"
    return {
        "release_id": RELEASE_ID,
        "dataset": spec["dataset"],
        "analysis_set_id": spec["analysis_set_id"],
        "contrast_or_exposure": spec["contrast_or_exposure"],
        "biological_unit": spec["biological_unit"],
        "biological_unit_resolution": spec["biological_unit_resolution"],
        "n_biological": "" if unresolved else 0,
        "technical_unit": spec["technical_unit"],
        "n_technical": 0,
        "design_status": status,
        "biological_ids_sha256": canonical_id_hash([]),
        "technical_ids_sha256": canonical_id_hash([]),
    }


def _write_gate(
    adapter_root: Path,
    spec: Mapping[str, object],
    seal,
    source_hash: str,
    execution_hash: str,
    reason: str,
    outcomes_tested: bool,
) -> None:
    write_tsv(
        adapter_root / "gate_status.tsv",
        SCHEMAS["gate_status.tsv"],
        [
            {
                "release_id": RELEASE_ID,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "gate_status": spec["gate_expectation"],
                "gate_reason": reason,
                "outcomes_tested": outcomes_tested,
                "registry_sha256": seal.registry_sha256,
                "source_manifest_sha256": source_hash,
                "execution_manifest_sha256": execution_hash,
            }
        ],
    )


def _gse_spec(seal) -> dict[str, object]:
    return {
        "adapter_id": "gse287826_skip",
        "adapter_root": "real_adapters/gse287826_skip",
        "dataset": "GSE287826",
        "assay": "GeoMx_WTA_CK8_18_positive",
        "source_publication": "GSE287826_NCBI_GEO_public_record",
        "source_dependence": "independent",
        "biological_unit": "unknown_public_biological_unit",
        "biological_unit_resolution": "unresolved",
        "technical_unit": "AOI",
        "analysis_set_id": "CK8_18_positive_MASH_vs_healthy",
        "contrast_or_exposure": "MASH_minus_healthy",
        "effect_unit": "donor_score_difference",
        "cross_assay_comparable": False,
        "program_universe": "external_test_eligible",
        "min_genes_testable": 8,
        "min_retained_l1_weight": 0.2,
        "gate_expectation": "skipped",
        "registry_sha256": seal.registry_sha256,
        "ready_sha256": seal.ready_sha256,
        "status": "real_skip_closed_no_donor_key",
    }


def _geomx_spec(seal) -> dict[str, object]:
    return {
        "adapter_id": "govaere2026_geomx_dropped",
        "adapter_root": "real_adapters/govaere2026_geomx_dropped",
        "dataset": "Govaere2026_GeoMx",
        "assay": "GeoMx_WTA",
        "source_publication": GOVAERE_CITATION,
        "source_dependence": "source_dependent",
        "biological_unit": "unknown_public_biological_unit",
        "biological_unit_resolution": "unresolved",
        "technical_unit": "segment",
        "analysis_set_id": "repository_exclusion_decision",
        "contrast_or_exposure": "all_local_GeoMx_reanalysis_contrasts",
        "effect_unit": "dropped_geomx_program_effect_not_applicable",
        "cross_assay_comparable": False,
        "program_universe": "external_test_eligible",
        "min_genes_testable": 8,
        "min_retained_l1_weight": 0.2,
        "gate_expectation": "skipped",
        "registry_sha256": seal.registry_sha256,
        "ready_sha256": seal.ready_sha256,
        "status": "real_skip_existing_geomx_exclusion",
    }


def _cosmx_spec(seal) -> dict[str, object]:
    return {
        "adapter_id": "govaere2026_cosmx_il32_context",
        "adapter_root": "real_adapters/govaere2026_cosmx_il32_context",
        "dataset": "Govaere2026_CosMx",
        "assay": "CosMx_SMI",
        "source_publication": GOVAERE_CITATION,
        "source_dependence": "source_dependent",
        "biological_unit": "unknown_public_biological_unit",
        "biological_unit_resolution": "unresolved",
        "technical_unit": "physical_array",
        "analysis_set_id": "IL32_hepatocyte_macrophage_proximity_MASH_arrays",
        "contrast_or_exposure": "IL32_high_vs_low_hepatocyte_nearest_macrophage_distance",
        "effect_unit": "per_array_delta_nearest_macrophage_distance_um",
        "cross_assay_comparable": False,
        "program_universe": "external_test_eligible",
        "min_genes_testable": 8,
        "min_retained_l1_weight": 0.2,
        "gate_expectation": "pass",
        "registry_sha256": seal.registry_sha256,
        "ready_sha256": seal.ready_sha256,
        "status": "real_source_native_context_no_program_refit",
    }


def _build_gse(
    project_root: Path,
    adapter_root: Path,
    spec: Mapping[str, object],
    seal,
    producer_relative: str,
    producer_sha256: str,
    created: str,
) -> None:
    acquisition = project_root / "Analysis/Spatial/candidates" / RELEASE_ID / "acquisition/gse287826"
    closed = project_root / "Analysis/Spatial/candidates" / RELEASE_ID / "gse287826"
    _, source_gate_rows = read_tsv(
        acquisition / "gate_status.tsv",
        ("status", "inference_authorized", "pairing_inference_prohibited"),
    )
    if len(source_gate_rows) != 1:
        raise ContractError("GSE287826 source gate must contain one row")
    source_gate = source_gate_rows[0]
    if (
        source_gate["status"] != "skipped_no_donor_key"
        or source_gate["inference_authorized"].strip().lower() != "false"
        or source_gate["pairing_inference_prohibited"].strip().lower() != "true"
    ):
        raise ContractError("GSE287826 no-donor-key terminal gate no longer holds")
    _, closed_rows = read_tsv(closed / "gate_status.tsv", ("status", "outcomes_read"))
    if len(closed_rows) != 1 or closed_rows[0]["status"] != "skipped_no_donor_key":
        raise ContractError("Plan 12 GSE287826 skip closure is absent or drifted")
    if closed_rows[0]["outcomes_read"].strip().lower() != "false":
        raise ContractError("Plan 12 GSE287826 closure says outcomes were read")

    adapter_root.mkdir(parents=True)
    source_hash = _source_manifest(
        project_root,
        adapter_root,
        (
            (acquisition / "gate_status.tsv", "terminal_public_donor_gate", True, "2026-08-07T22:36:23Z"),
            (acquisition / "source_files.tsv", "public_source_file_manifest", True, "2026-08-07T22:32:04Z"),
            (acquisition / "gse287826_record_reconciliation.tsv", "public_record_reconciliation", True, "2026-08-07T22:36:23Z"),
            (closed / "gate_status.tsv", "plan12_outcome_blind_skip_closure", True, "2026-08-07T22:36:23Z"),
            (closed / "execution_manifest.tsv", "plan12_skip_execution_manifest", True, "2026-08-07T22:36:23Z"),
        ),
    )
    rows = _base_program_rows(
        spec,
        seal,
        source_hash,
        producer_relative,
        producer_sha256,
        {},
        None,
        0,
        "skipped",
        "skipped_no_authoritative_30_AOI_to_15_donor_key;outcomes_not_opened",
    )
    _write_adapter_tables(adapter_root, [], _empty_design(spec, "skipped_no_donor_key"), *rows)
    execution_hash = _execution_manifest(
        adapter_root,
        (
            "sample_manifest.tsv",
            "gene_mapping_audit.tsv",
            "design_audit.tsv",
            "program_testability.tsv",
            "per_sample_program_scores.tsv",
            "program_effects.tsv",
            "sensitivity.tsv",
            "source_manifest.tsv",
        ),
        created,
    )
    _write_gate(
        adapter_root,
        spec,
        seal,
        source_hash,
        execution_hash,
        "skipped_no_donor_key;30_AOIs_never_treated_as_30_biological_replicates",
        False,
    )


def _build_geomx_drop(
    project_root: Path,
    adapter_root: Path,
    spec: Mapping[str, object],
    seal,
    producer_relative: str,
    producer_sha256: str,
    created: str,
) -> None:
    dropped = project_root / "Analysis/Spatial/results/govaere2026/GEOMX_DROPPED.md"
    text = dropped.read_text(encoding="utf-8")
    required_phrases = (
        "NOT used as an evidence arm anywhere in the paper",
        "0 significant genes",
        "result is not trustworthy",
    )
    if any(phrase not in text for phrase in required_phrases):
        raise ContractError("Govaere GeoMx exclusion decision drifted")
    adapter_root.mkdir(parents=True)
    source_hash = _source_manifest(
        project_root,
        adapter_root,
        (
            (dropped, "authoritative_local_geomx_exclusion_decision", True, "2026-06-16T17:10:20Z"),
            (project_root / "data/external/govaere2026_natgenetics/README.md", "public_source_and_license_inventory", True, "2026-05-21"),
        ),
    )
    rows = _base_program_rows(
        spec,
        seal,
        source_hash,
        producer_relative,
        producer_sha256,
        {},
        None,
        0,
        "skipped",
        "repository_excludes_local_geomx_reanalysis_due_unresolved_LOQ_or_orientation_problem",
    )
    _write_adapter_tables(adapter_root, [], _empty_design(spec, "skipped_repository_exclusion"), *rows)
    execution_hash = _execution_manifest(
        adapter_root,
        (
            "sample_manifest.tsv",
            "gene_mapping_audit.tsv",
            "design_audit.tsv",
            "program_testability.tsv",
            "per_sample_program_scores.tsv",
            "program_effects.tsv",
            "sensitivity.tsv",
            "source_manifest.tsv",
        ),
        created,
    )
    _write_gate(
        adapter_root,
        spec,
        seal,
        source_hash,
        execution_hash,
        "historical_local_GeoMx_outputs_exist_but_are_explicitly_dropped_and_not_imported",
        False,
    )


def measured_program_coverage(
    membership_rows: Sequence[Mapping[str, str]],
    measured_symbols: set[str],
    uids: Iterable[str],
) -> dict[str, tuple[int, float, str]]:
    result: dict[str, tuple[int, float, str]] = {}
    for uid in uids:
        rows = [row for row in membership_rows if row["program_uid"] == uid]
        seen: set[str] = set()
        retained = 0.0
        for row in rows:
            symbol = row["mapped_symbol"].strip()
            if not symbol or symbol.upper() == "NA" or symbol in seen:
                continue
            seen.add(symbol)
            if symbol in measured_symbols:
                retained += float(row["original_l1_weight"])
        measured = seen & measured_symbols
        result[uid] = (
            len(measured),
            retained,
            "panel_coverage_audited_but_source_native_result_is_IL32_specific_no_program_refit",
        )
    return result


def source_native_cosmx_rows(source_rows: Sequence[Mapping[str, str]]) -> list[dict[str, object]]:
    required_slides = {"Leuven_1", "Leuven_2", "Leuven_3", "Leuven_4"}
    observed = {row["slide"] for row in source_rows}
    if observed != required_slides or len(source_rows) != 4:
        raise ContractError(f"CosMx IL32 source must contain the four frozen arrays, found {sorted(observed)}")
    result = []
    for row in sorted(source_rows, key=lambda item: item["slide"]):
        is_mash = parse_bool(row["is_mash"], f"CosMx.is_mash[{row['slide']}]")
        concordant = parse_bool(
            row["il32high_closer"],
            f"CosMx.il32high_closer[{row['slide']}]",
        )
        result.append(
            {
                "release_id": RELEASE_ID,
                "dataset": "Govaere2026_CosMx",
                "assay": "CosMx_SMI",
                "source_publication": GOVAERE_CITATION,
                "source_dependence": "source_dependent",
                "reporting_unit": "physical_array",
                "reporting_id": row["slide"],
                "is_mash": is_mash,
                "include_direction_concordance": is_mash,
                "n_hepatocytes": row["n_hep"],
                "n_macrophages": row["n_mac"],
                "n_metmac": row["n_metmac"],
                "il32_partition": row["il32_mode"],
                "contrast_or_exposure": "IL32_high_vs_low_hepatocyte_nearest_macrophage_distance",
                "effect_unit": "per_array_delta_nearest_macrophage_distance_um",
                "estimate": row["delta_high_minus_low_um"],
                "secondary_effect_unit": "per_array_spearman_rho_IL32_vs_mean_5NN_macrophage_CD74",
                "secondary_estimate": row["rho_il32_vs_knn_cd74"],
                "direction_concordant": concordant,
                "inferential_status": (
                    "MASH_array_direction_concordance" if is_mash
                    else "mixed_normal_reference_reported_not_in_concordance"
                ),
                "interpretation": (
                    "source_dependent_proximity_context_not_whole_program_validation;"
                    "cell_level_p_values_excluded_as_pseudoreplicated"
                ),
            }
        )
    mash = [row for row in result if row["include_direction_concordance"]]
    if len(mash) != 3 or not all(bool(row["direction_concordant"]) for row in mash):
        raise ContractError("accepted CosMx source must rederive 3/3 MASH-array direction concordance")
    return result


def _build_cosmx(
    project_root: Path,
    adapter_root: Path,
    spec: Mapping[str, object],
    seal,
    producer_relative: str,
    producer_sha256: str,
    created: str,
) -> None:
    source_dir = project_root / "Analysis/Spatial/results/govaere2026/il32_colocalization"
    per_slide_path = source_dir / "il32_coloc_per_slide.csv"
    summary_path = source_dir / "il32_coloc_summary.txt"
    original_producer = project_root / "Analysis/Spatial/scripts/42f_govaere2026_cosmx_il32_coloc.py"
    panel_path = project_root / "Analysis/Spatial/results/govaere2026/cosmx_de_Hepatocyte_MASH_vs_noMASH.csv"
    membership_path = (
        project_root
        / "Analysis/Multimodal_Program_Projection/candidates"
        / RELEASE_ID
        / "hotspot/program_membership_v2.tsv"
    )
    source_rows = _read_csv(per_slide_path)
    native_rows = source_native_cosmx_rows(source_rows)
    panel_rows = _read_csv(panel_path)
    measured_symbols = {row["gene"].strip() for row in panel_rows if row["gene"].strip()}
    if len(measured_symbols) != 968:
        raise ContractError(f"CosMx measured panel must contain 968 unique genes, found {len(measured_symbols)}")
    _, membership_rows = read_tsv(
        membership_path,
        ("program_uid", "mapped_symbol", "original_l1_weight", "membership_sha256"),
    )
    universe = program_universe(seal, str(spec["program_universe"]))
    coverage = measured_program_coverage(membership_rows, measured_symbols, universe)

    adapter_root.mkdir(parents=True)
    source_hash = _source_manifest(
        project_root,
        adapter_root,
        (
            (per_slide_path, "accepted_source_native_per_array_context", True, "2026-06-16T17:09:07Z"),
            (summary_path, "accepted_source_native_direction_summary", True, "2026-06-16T17:09:08Z"),
            (original_producer, "source_native_producer", True, "2026-07-13T17:19:04Z"),
            (panel_path, "frozen_968_gene_panel_universe_only", True, "2026-06-16"),
            (project_root / "data/external/govaere2026_natgenetics/README.md", "public_source_and_license_inventory", True, "2026-05-21"),
            (project_root / "data/external/govaere2026_natgenetics/cosmx/manifest.tsv", "public_source_checksum_manifest", True, "2026-05-21T15:53:18Z"),
        ),
    )
    sample_rows = []
    for row in native_rows:
        include = bool(row["include_direction_concordance"])
        sample_rows.append(
            {
                "release_id": RELEASE_ID,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "biological_id": "",
                "technical_id": row["reporting_id"],
                "biological_unit": spec["biological_unit"],
                "technical_unit": spec["technical_unit"],
                "include_primary": include,
                "gate_state": (
                    "included_MASH_direction_concordance" if include
                    else "mixed_normal_reference_reported_not_in_concordance"
                ),
            }
        )
    included = [row for row in sample_rows if bool(row["include_primary"])]
    biological_ids: list[str] = []
    technical_ids = [str(row["technical_id"]) for row in included]
    design = {
        "release_id": RELEASE_ID,
        "dataset": spec["dataset"],
        "analysis_set_id": spec["analysis_set_id"],
        "contrast_or_exposure": spec["contrast_or_exposure"],
        "biological_unit": spec["biological_unit"],
        "biological_unit_resolution": spec["biological_unit_resolution"],
        "n_biological": "",
        "technical_unit": spec["technical_unit"],
        "n_technical": len(technical_ids),
        "design_status": "source_dependent_three_MASH_array_direction_concordance",
        "biological_ids_sha256": canonical_id_hash(biological_ids),
        "technical_ids_sha256": canonical_id_hash(technical_ids),
    }
    rows = _base_program_rows(
        spec,
        seal,
        source_hash,
        producer_relative,
        producer_sha256,
        coverage,
        None,
        len(technical_ids),
        "not_applicable",
        (
            "lineage_relevant_but_source_native_result_is_IL32_specific;"
            "no_frozen_whole_program_estimand;program_refit_prohibited"
        ),
    )
    _write_adapter_tables(adapter_root, sample_rows, design, *rows)
    write_tsv(adapter_root / "source_native_context.tsv", SOURCE_NATIVE_COLUMNS, native_rows)
    execution_hash = _execution_manifest(
        adapter_root,
        (
            "sample_manifest.tsv",
            "gene_mapping_audit.tsv",
            "design_audit.tsv",
            "program_testability.tsv",
            "per_sample_program_scores.tsv",
            "program_effects.tsv",
            "sensitivity.tsv",
            "source_native_context.tsv",
            "source_manifest.tsv",
        ),
        created,
    )
    _write_gate(
        adapter_root,
        spec,
        seal,
        source_hash,
        execution_hash,
        (
            "accepted_existing_IL32_hepatocyte_macrophage_array_context;"
            "no_whole_program_refit_or_program_effect_claim"
        ),
        True,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--candidate-root", type=Path, default=None)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    candidate_root = args.candidate_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "spatial_context"
    )
    hotspot_root = (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "hotspot"
    )
    try:
        if not (candidate_root / "map_selection_manifest.tsv").is_file():
            raise ContractError("map_selection_manifest.tsv must be frozen before real adapters")
        target = candidate_root / "real_adapters"
        registry_target = candidate_root / "real_adapter_registry.tsv"
        if target.exists() or registry_target.exists():
            raise ContractError("refusing to overwrite an existing real-adapter candidate")
        seal = validate_hotspot_seal(hotspot_root)
        producer = Path(__file__).resolve()
        producer_relative = producer.relative_to(project_root).as_posix()
        producer_sha256 = sha256_file(producer)
        created = utc_now()
        staging = candidate_root / f".real_adapters.incomplete.{os.getpid()}"
        registry_staging = candidate_root / f".real_adapter_registry.incomplete.{os.getpid()}.tsv"
        if staging.exists() or registry_staging.exists():
            raise ContractError("staging path already exists")
        staging.mkdir(parents=True)

        specs = (_gse_spec(seal), _geomx_spec(seal), _cosmx_spec(seal))
        _build_gse(project_root, staging / "gse287826_skip", specs[0], seal, producer_relative, producer_sha256, created)
        _build_geomx_drop(project_root, staging / "govaere2026_geomx_dropped", specs[1], seal, producer_relative, producer_sha256, created)
        _build_cosmx(project_root, staging / "govaere2026_cosmx_il32_context", specs[2], seal, producer_relative, producer_sha256, created)
        write_tsv(registry_staging, ADAPTER_REGISTRY_COLUMNS, specs)

        os.replace(staging, target)
        os.replace(registry_staging, registry_target)
        print("PASS: built three isolated real adapters without program refitting")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
