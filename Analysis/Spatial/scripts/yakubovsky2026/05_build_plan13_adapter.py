#!/usr/bin/env python3
"""Adapt sealed Plan 11 results to the Plan 13 assay-native long contract.

This bridge never refits a lipid model.  It keeps native donor slopes and
signed-Stouffer P values in their native unit, supplies the descriptive donor
effect range (not a confidence interval), and exports donor mean frozen program
scores only as sample-level context.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from yakubovsky_common import (
    DATASET,
    RELEASE_ID,
    STAGE_ARTIFACTS,
    ContractError,
    atomic_write_frame,
    bool_value,
    default_paths,
    sha256_file,
    validate_registry_contract,
    validate_stage_chain,
)


CONTRACT_LIB_SHA256 = "93cbab090e3b53906b034d60cccb4608757e42d183e1a850dd7fc0392d6186c1"
CONTRACT_VALIDATOR_SHA256 = "aab3d5cd0d4f88ef7e377acdd75094048335cc54dc8b7d2d53b7cf430d2490b0"
ANALYSIS_SET_ID = "yakubovsky_binary_lipid_context"
ADAPTER_ID = "yakubovsky2026_binary_lipid"
ASSAY = "Visium_lipid_context"
SOURCE_PUBLICATION = "Yakubovsky_et_al_Nature_2026_s41586-026-10377-y"
SOURCE_DEPENDENCE = "source_dependent"
CONTRAST = "source_defined_lipid_zone_vs_non_lipid_zone"
EFFECT_UNIT = "source_lipid_zone_minus_non_lipid_zone_score"
ROBUST_VARIANTS = ("equal_weight", "leave_top_gene", "zonation_df3", "zonation_df5")


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_contract_module(base: Path):
    script_dir = base / "Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
    contract_path = script_dir / "contract_lib.py"
    validator_path = script_dir / "03_validate_contract.py"
    if sha256_file(contract_path) != CONTRACT_LIB_SHA256:
        raise ContractError("Plan 13 contract_lib.py hash drift")
    if sha256_file(validator_path) != CONTRACT_VALIDATOR_SHA256:
        raise ContractError("Plan 13 validator hash drift")
    spec = importlib.util.spec_from_file_location("plan13_contract_lib", contract_path)
    if spec is None or spec.loader is None:
        raise ContractError("Cannot import the sealed Plan 13 contract")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_native_ready(native: Path, registry) -> pd.Series:
    validate_stage_chain(native, "terminal")
    ready = pd.read_csv(native / "READY", sep="\t", dtype=str, keep_default_na=False)
    if len(ready) != 1:
        raise ContractError("Native Plan 11 READY must contain exactly one row")
    row = ready.iloc[0]
    expected = {
        "release_id": RELEASE_ID,
        "status": "ready_for_plan13_candidate_integration",
        "registry_sha256": registry.registry_sha256,
        "membership_sha256": registry.membership_sha256,
        "canonical_promotion_authorized": "False",
    }
    for field, value in expected.items():
        if row.get(field) != value:
            raise ContractError(
                f"Native Plan 11 READY mismatch for {field}: {row.get(field)!r} != {value!r}"
            )
    if sha256_file(native / "gate_status.tsv") != row["gate_status_sha256"]:
        raise ContractError("Native Plan 11 READY/gate hash mismatch")
    return row


def native_evidence_state(value: str) -> str:
    if value == "robust":
        return "robust"
    if value == "tested_negative":
        return "tested_negative"
    if value == "untestable":
        return "untestable"
    raise ContractError(f"Unexpected native Plan 11 evidence state: {value!r}")


def contract_direction(value: float) -> str:
    if not np.isfinite(value):
        return "not_applicable"
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "zero"


def source_rows(
    base: Path,
    native: Path,
    paths: dict[str, Path],
    created: str,
) -> list[dict[str, Any]]:
    native_roles = {
        "READY": "sealed_native_plan11_ready",
        "gate_status.tsv": "native_plan11_gate",
        "validation_report.tsv": "native_plan11_validation",
        "program_effects.tsv": "native_plan11_program_effects",
        "donor_program_effects.tsv": "native_plan11_donor_program_effects",
        "sensitivity.tsv": "native_plan11_sensitivities",
        "multiplicity_manifest.tsv": "native_plan11_multiplicity_family",
        "program_testability.tsv": "native_plan11_program_testability",
        "per_sample_program_scores.tsv": "native_plan11_spot_program_scores",
        "model_sample_manifest.tsv": "native_plan11_spot_inclusion",
        "source_manifest.tsv": "native_plan11_source_hash_chain",
    }
    for stage, artifacts in STAGE_ARTIFACTS.items():
        for name in artifacts:
            native_roles.setdefault(name, f"sealed_native_plan11_{stage}_artifact")
        native_roles[f"{stage}_stage_seal.tsv"] = (
            f"sealed_native_plan11_{stage}_seal"
        )
    plan20_sources = {
        "plan20_registry": paths["registry"],
        "plan20_membership": paths["membership"],
        "plan20_external_family": paths["external"],
        "plan20_tested_universe": paths["tested_universe"],
        "plan20_validation": paths["hotspot_validation"],
        "plan20_release_manifest": paths["hotspot_release_manifest"],
        "plan20_ready": paths["hotspot_ready"],
    }
    rows = []
    for name, role in sorted(native_roles.items()):
        path = native / name
        if not path.is_file():
            raise ContractError(f"Missing consumed native Plan 11 artifact: {path}")
        rows.append(
            {
                "path_scope": "project_relative",
                "relative_path": str(path.relative_to(base)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "source_role": role,
                "public_access": False,
                "retrieved_utc": created,
            }
        )
    for role, path in sorted(plan20_sources.items()):
        if not path.is_file():
            raise ContractError(f"Missing consumed Plan 20 artifact: {path}")
        rows.append(
            {
                "path_scope": "project_relative",
                "relative_path": str(path.relative_to(base)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "source_role": role,
                "public_access": False,
                "retrieved_utc": created,
            }
        )
    return rows


def check_existing(output: Path, base: Path, contract, native: Path, paths) -> dict[str, Any]:
    registry = validate_registry_contract(paths)
    native_ready = read_native_ready(native, registry)
    checks = contract.validate_candidate_contract(
        output,
        output / "adapter_registry.tsv",
        paths["hotspot"],
        base,
    )
    ready_path = output / "PLAN13_ADAPTER_READY"
    ready = pd.read_csv(ready_path, sep="\t", dtype=str, keep_default_na=False)
    if len(ready) != 1 or ready.iloc[0]["status"] != "validated_plan13_adapter_candidate":
        raise ContractError("Plan 13 adapter READY is absent or invalid")
    row = ready.iloc[0]
    expected_hashes = {
        "native_ready_sha256": sha256_file(native / "READY"),
        "adapter_registry_sha256": sha256_file(output / "adapter_registry.tsv"),
        "gate_status_sha256": sha256_file(output / "gate_status.tsv"),
        "validation_report_sha256": sha256_file(output / "validation_report.tsv"),
        "contract_lib_sha256": CONTRACT_LIB_SHA256,
        "contract_validator_sha256": CONTRACT_VALIDATOR_SHA256,
    }
    for field, value in expected_hashes.items():
        if row.get(field) != value:
            raise ContractError(f"Plan 13 adapter READY hash drift for {field}")
    return {
        "release_id": RELEASE_ID,
        "status": "validated_plan13_adapter_candidate",
        "n_contract_checks": len(checks),
        "native_ready_sha256": sha256_file(native / "READY"),
        "native_terminal_status": native_ready["status"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path)
    parser.add_argument("--native-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.check_only == args.write:
        raise SystemExit("Choose exactly one of --check-only or --write")
    paths = default_paths(args.base)
    base = paths["base"]
    native = (args.native_dir or paths["candidate"]).resolve()
    output = (args.output_dir or (native / "plan13_adapter")).resolve()
    contract = load_contract_module(base)
    if args.check_only:
        print(json.dumps(check_existing(output, base, contract, native, paths), indent=2, sort_keys=True))
        return

    protected = [
        *contract.SCHEMAS,
        "adapter_registry.tsv",
        "validation_report.tsv",
        "PLAN13_ADAPTER_READY",
    ]
    existing = [str(output / name) for name in protected if (output / name).exists()]
    if existing:
        raise ContractError("Refusing to overwrite Plan 13 adapter artifact(s): " + ", ".join(existing))
    output.mkdir(parents=True, exist_ok=True)

    registry = validate_registry_contract(paths)
    read_native_ready(native, registry)
    gate = pd.read_csv(native / "gate_status.tsv", sep="\t", dtype=str, keep_default_na=False)
    if len(gate) != 1:
        raise ContractError("Native Plan 11 gate must contain one row")
    native_status = gate.iloc[0]["status"]
    programs = pd.read_csv(
        native / "program_effects.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    testability = pd.read_csv(
        native / "program_testability.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    donor_effects = pd.read_csv(
        native / "donor_program_effects.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    sensitivities = pd.read_csv(
        native / "sensitivity.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    multiplicity = pd.read_csv(
        native / "multiplicity_manifest.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    model_samples = pd.read_csv(
        native / "model_sample_manifest.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    spot_scores = pd.read_csv(
        native / "per_sample_program_scores.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    family = registry.family.set_index("program_uid", drop=False)
    if set(programs["program_uid"]) != set(family.index):
        raise ContractError("Native Plan 11 effect family drift")

    testability["testable"] = testability["testable"].str.lower().isin(["true", "1"])
    donor_effects["estimable"] = donor_effects["estimable"].str.lower().isin(["true", "1"])
    programs["robust"] = programs["robust"].str.lower().isin(["true", "1"])
    model_samples["primary_model_eligible"] = model_samples[
        "primary_model_eligible"
    ].str.lower().isin(["true", "1"])
    for frame, columns in (
        (
            programs,
            ("estimate", "ci_lower", "ci_upper", "pvalue", "qvalue", "n_donors"),
        ),
        (donor_effects, ("observed_beta", "bootstrap_pvalue")),
        (sensitivities, ("combined_z", "median_donor_slope")),
        (spot_scores, ("primary_score",)),
    ):
        for column in columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")

    if multiplicity["variant"].duplicated().any():
        raise ContractError("Native multiplicity manifest has duplicate variants")
    family_complete_by_variant = {
        str(row.variant): bool_value(row.family_complete)
        for row in multiplicity.itertuples(index=False)
    }
    missing_multiplicity = sorted(
        set(sensitivities["variant"]).difference(family_complete_by_variant)
    )
    if missing_multiplicity:
        raise ContractError(
            "Native multiplicity manifest lacks sensitivity variant(s): "
            + ", ".join(missing_multiplicity)
        )

    native_testable_uids = set(
        testability.loc[testability["testable"], "program_uid"]
    )
    common_donors: set[str] = set()
    complete_native = native_status == "complete_binary_lipid_analysis"
    if complete_native:
        donor_sets = []
        for uid in sorted(native_testable_uids):
            donor_sets.append(
                set(
                    donor_effects.loc[
                        (donor_effects["program_uid"] == uid)
                        & (donor_effects["variant"] == "primary")
                        & donor_effects["estimable"],
                        "donor",
                    ]
                )
            )
        if not donor_sets or len({tuple(sorted(values)) for values in donor_sets}) != 1:
            raise ContractError(
                "Plan 13 one-analysis-set contract requires identical valid primary donors "
                "across every assay-testable frozen program"
            )
        common_donors = donor_sets[0]
        if len(common_donors) < 3:
            raise ContractError("Native complete branch has fewer than three common donors")
        gate_expectation = "pass" if programs["robust"].any() else "valid_null"
        outcomes_tested = True
    elif native_status == "complete_zonation_reference_spatial_model_unestimable":
        gate_expectation = "untestable"
        outcomes_tested = False
    else:
        raise ContractError(f"Unsupported native Plan 11 terminal status: {native_status}")

    included_spots = model_samples[
        model_samples["primary_model_eligible"]
        & model_samples["donor"].isin(common_donors)
    ].copy()
    included_spot_ids = set(included_spots["spot_id"])
    created = utc_now()
    source_manifest_rows = source_rows(base, native, paths, created)
    contract.write_tsv(
        output / "source_manifest.tsv",
        contract.SCHEMAS["source_manifest.tsv"],
        source_manifest_rows,
    )
    source_manifest_hash = sha256_file(output / "source_manifest.tsv")

    sample_rows = [
        {
            "release_id": RELEASE_ID,
            "dataset": DATASET,
            "analysis_set_id": ANALYSIS_SET_ID,
            "biological_id": row.donor,
            "technical_id": row.spot_id,
            "biological_unit": "donor",
            "technical_unit": "spot",
            "include_primary": True,
            "gate_state": "included_primary_model",
        }
        for row in included_spots.itertuples(index=False)
    ]
    contract.write_tsv(
        output / "sample_manifest.tsv",
        contract.SCHEMAS["sample_manifest.tsv"],
        sample_rows,
    )
    biological_ids = sorted(common_donors)
    technical_ids = sorted(included_spot_ids)
    design_row = {
        "release_id": RELEASE_ID,
        "dataset": DATASET,
        "analysis_set_id": ANALYSIS_SET_ID,
        "contrast_or_exposure": CONTRAST,
        "biological_unit": "donor",
        "biological_unit_resolution": "resolved",
        "n_biological": len(biological_ids),
        "technical_unit": "spot",
        "n_technical": len(technical_ids),
        "design_status": (
            "passed_native_plan11_binary_spatial_model"
            if complete_native
            else "native_plan11_spatial_model_unestimable"
        ),
        "biological_ids_sha256": contract.canonical_id_hash(biological_ids),
        "technical_ids_sha256": contract.canonical_id_hash(technical_ids),
    }
    contract.write_tsv(
        output / "design_audit.tsv",
        contract.SCHEMAS["design_audit.tsv"],
        [design_row],
    )

    producer = str(Path(__file__).resolve().relative_to(base))
    producer_hash = sha256_file(Path(__file__).resolve())
    mapping_rows: list[dict[str, Any]] = []
    testability_rows: list[dict[str, Any]] = []
    effect_rows: list[dict[str, Any]] = []
    score_rows: list[dict[str, Any]] = []
    sensitivity_rows: list[dict[str, Any]] = []
    technical_counts = included_spots.groupby("donor").size().to_dict()
    for uid, family_row in family.iterrows():
        test_row = testability[testability["program_uid"] == uid]
        native_row = programs[programs["program_uid"] == uid]
        if len(test_row) != 1 or len(native_row) != 1:
            raise ContractError(f"Missing or duplicate native program row: {uid}")
        test_row = test_row.iloc[0]
        native_row = native_row.iloc[0]
        native_gene_testable = bool(test_row["testable"])
        generic_testable = complete_native and native_gene_testable
        if not complete_native:
            evidence_state = "untestable"
            reason = "native_plan11_spatial_model_unestimable"
        elif not native_gene_testable:
            evidence_state = "untestable"
            reason = test_row["testability_reason"]
        else:
            evidence_state = native_evidence_state(native_row["evidence_state"])
            reason = test_row["testability_reason"]

        n_source_genes = int(test_row["n_original_membership_rows"])
        n_measured = int(test_row["n_mapped_genes"])
        retained = float(test_row["retained_original_l1_weight"])
        mapping_rows.append(
            {
                "release_id": RELEASE_ID,
                "dataset": DATASET,
                "analysis_set_id": ANALYSIS_SET_ID,
                "program_uid": uid,
                "membership_sha256": family_row["membership_sha256"],
                "n_source_genes": n_source_genes,
                "n_genes_measured": n_measured,
                "retained_l1_weight": retained,
                "mapping_status": reason,
            }
        )
        testability_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_sha256": registry.registry_sha256,
                "dataset": DATASET,
                "analysis_set_id": ANALYSIS_SET_ID,
                "program_uid": uid,
                "membership_sha256": family_row["membership_sha256"],
                "testable": generic_testable,
                "n_genes_measured": n_measured,
                "retained_l1_weight": retained,
                "testability_reason": reason,
                "evidence_state": evidence_state,
            }
        )

        sensitivity_sign_agree = False
        if generic_testable:
            estimate = float(native_row["estimate"])
            observed_direction = contract_direction(estimate)
            native_sensitivity = sensitivities[
                (sensitivities["program_uid"] == uid)
                & sensitivities["variant"].isin(ROBUST_VARIANTS)
            ]
            sensitivity_sign_agree = len(native_sensitivity) == len(ROBUST_VARIANTS) and all(
                family_complete_by_variant.get(str(row.variant), False)
                and contract_direction(float(row.combined_z)) == observed_direction
                for row in native_sensitivity.itertuples(index=False)
            )
            interval_low = float(native_row["ci_lower"])
            interval_high = float(native_row["ci_upper"])
            if not interval_low <= estimate <= interval_high:
                raise ContractError(f"Native donor slope range does not contain estimate: {uid}")
            numeric = {
                "estimate": estimate,
                "std_error": "",
                "interval_low": interval_low,
                "interval_high": interval_high,
                "interval_type": "donor_effect_range",
                "pvalue": float(native_row["pvalue"]),
                "padj": float(native_row["qvalue"]),
                "pvalue_method": "equal_donor_signed_stouffer_of_spatial_block_bootstrap_sign_tail_pvalues",
            }
            donor_score_source = spot_scores[
                (spot_scores["program_uid"] == uid)
                & spot_scores["spot_id"].isin(included_spot_ids)
            ]
            for donor in biological_ids:
                donor_values = donor_score_source.loc[
                    donor_score_source["donor"] == donor, "primary_score"
                ].to_numpy(float)
                if len(donor_values) != technical_counts[donor]:
                    raise ContractError(f"Donor score/spot count mismatch for {uid}/{donor}")
                score_rows.append(
                    {
                        "release_id": RELEASE_ID,
                        "registry_sha256": registry.registry_sha256,
                        "dataset": DATASET,
                        "analysis_set_id": ANALYSIS_SET_ID,
                        "program_uid": uid,
                        "membership_sha256": family_row["membership_sha256"],
                        "biological_id": donor,
                        "technical_unit_count": technical_counts[donor],
                        "program_score": float(np.mean(donor_values)),
                        "score_unit": "donor_mean_frozen_weighted_within_donor_gene_z_score",
                    }
                )
        else:
            observed_direction = "not_applicable"
            numeric = {
                "estimate": "",
                "std_error": "",
                "interval_low": "",
                "interval_high": "",
                "interval_type": "not_applicable",
                "pvalue": "",
                "padj": "",
                "pvalue_method": "not_applicable",
            }

        effect_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_sha256": registry.registry_sha256,
                "membership_sha256": family_row["membership_sha256"],
                "program_uid": uid,
                "legacy_program_id": f"{family_row['cell_type']}:{family_row['module']}",
                "program_label": family_row["module_name"],
                "cell_type": family_row["cell_type"],
                "dataset": DATASET,
                "assay": ASSAY,
                "source_publication": SOURCE_PUBLICATION,
                "source_dependence": SOURCE_DEPENDENCE,
                "analysis_set_id": ANALYSIS_SET_ID,
                "biological_unit": "donor",
                "biological_unit_resolution": "resolved",
                "n_biological": len(biological_ids),
                "technical_unit": "spot",
                "n_technical": len(technical_ids),
                "contrast_or_exposure": CONTRAST,
                "effect_unit": EFFECT_UNIT,
                **numeric,
                "multiplicity_family": (
                    f"BH_complete_assay_testable_external_hepatocyte_family_n{len(native_testable_uids)}"
                ),
                "n_genes_measured": n_measured,
                "retained_l1_weight": retained,
                "testable": generic_testable,
                "testability_reason": reason,
                "direction_expected": family_row["expected_direction"],
                "direction_observed": observed_direction,
                "sensitivity_sign_agree": sensitivity_sign_agree,
                "robustness_pass": bool(native_row["robust"]) if generic_testable else False,
                "cross_assay_comparable": False,
                "evidence_state": evidence_state,
                "producer": producer,
                "producer_sha256": producer_hash,
                "source_manifest_sha256": source_manifest_hash,
            }
        )
        if generic_testable:
            primary_direction = contract_direction(float(native_row["estimate"]))
            for sensitivity_row in sensitivities[sensitivities["program_uid"] == uid].itertuples(
                index=False
            ):
                estimate = float(sensitivity_row.median_donor_slope)
                observed = contract_direction(estimate)
                sensitivity_rows.append(
                    {
                        "release_id": RELEASE_ID,
                        "registry_sha256": registry.registry_sha256,
                        "dataset": DATASET,
                        "analysis_set_id": ANALYSIS_SET_ID,
                        "program_uid": uid,
                        "membership_sha256": family_row["membership_sha256"],
                        "sensitivity_id": sensitivity_row.variant,
                        "estimate": estimate,
                        "effect_unit": EFFECT_UNIT,
                        "direction_observed": observed,
                        "sign_agree": observed == primary_direction,
                        "status": (
                            "complete_prespecified_family"
                            if family_complete_by_variant.get(
                                str(sensitivity_row.variant), False
                            )
                            and np.isfinite(sensitivity_row.combined_z)
                            else "unestimable_complete_family"
                            if family_complete_by_variant.get(
                                str(sensitivity_row.variant), False
                            )
                            else "incomplete_prespecified_family"
                        ),
                    }
                )

    for filename, rows in (
        ("gene_mapping_audit.tsv", mapping_rows),
        ("program_testability.tsv", testability_rows),
        ("per_sample_program_scores.tsv", score_rows),
        ("program_effects.tsv", effect_rows),
        ("sensitivity.tsv", sensitivity_rows),
    ):
        contract.write_tsv(output / filename, contract.SCHEMAS[filename], rows)

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
    execution_rows = [
        {
            "role": "plan11_to_plan13_adapter_artifact",
            "relative_path": filename,
            "bytes": (output / filename).stat().st_size,
            "sha256": sha256_file(output / filename),
            "created_utc": created,
        }
        for filename in execution_files
    ]
    contract.write_tsv(
        output / "execution_manifest.tsv",
        contract.SCHEMAS["execution_manifest.tsv"],
        execution_rows,
    )
    gate_reason = (
        "native_plan11_complete_with_at_least_one_robust_program"
        if gate_expectation == "pass"
        else "native_plan11_complete_valid_all_null_or_nonrobust_family"
        if gate_expectation == "valid_null"
        else "native_plan11_spatial_model_unestimable"
    )
    gate_row = {
        "release_id": RELEASE_ID,
        "dataset": DATASET,
        "analysis_set_id": ANALYSIS_SET_ID,
        "gate_status": gate_expectation,
        "gate_reason": gate_reason,
        "outcomes_tested": outcomes_tested,
        "registry_sha256": registry.registry_sha256,
        "source_manifest_sha256": source_manifest_hash,
        "execution_manifest_sha256": sha256_file(output / "execution_manifest.tsv"),
    }
    contract.write_tsv(
        output / "gate_status.tsv", contract.SCHEMAS["gate_status.tsv"], [gate_row]
    )
    adapter_row = {
        "adapter_id": ADAPTER_ID,
        "adapter_root": ".",
        "dataset": DATASET,
        "assay": ASSAY,
        "source_publication": SOURCE_PUBLICATION,
        "source_dependence": SOURCE_DEPENDENCE,
        "biological_unit": "donor",
        "biological_unit_resolution": "resolved",
        "technical_unit": "spot",
        "analysis_set_id": ANALYSIS_SET_ID,
        "contrast_or_exposure": CONTRAST,
        "effect_unit": EFFECT_UNIT,
        "cross_assay_comparable": False,
        "program_universe": "external_test_eligible",
        "min_genes_testable": 8,
        "min_retained_l1_weight": 0.2,
        "gate_expectation": gate_expectation,
        "registry_sha256": registry.registry_sha256,
        "ready_sha256": registry.ready_sha256,
        "status": "validated_candidate_adapter",
    }
    contract.write_tsv(
        output / "adapter_registry.tsv", contract.ADAPTER_REGISTRY_COLUMNS, [adapter_row]
    )
    checks = contract.validate_candidate_contract(
        output,
        output / "adapter_registry.tsv",
        paths["hotspot"],
        base,
    )
    validation = pd.DataFrame(
        [
            {
                "release_id": RELEASE_ID,
                "check_id": check.check_id,
                "status": check.status,
                "detail": check.detail,
                "contract_lib_sha256": CONTRACT_LIB_SHA256,
                "validated_utc": utc_now(),
            }
            for check in checks
        ]
    )
    atomic_write_frame(output / "validation_report.tsv", validation)
    adapter_ready = pd.DataFrame(
        [
            {
                "release_id": RELEASE_ID,
                "status": "validated_plan13_adapter_candidate",
                "native_ready_sha256": sha256_file(native / "READY"),
                "adapter_registry_sha256": sha256_file(output / "adapter_registry.tsv"),
                "gate_status_sha256": sha256_file(output / "gate_status.tsv"),
                "validation_report_sha256": sha256_file(output / "validation_report.tsv"),
                "contract_lib_sha256": CONTRACT_LIB_SHA256,
                "contract_validator_sha256": CONTRACT_VALIDATOR_SHA256,
                "canonical_promotion_authorized": False,
                "completed_utc": utc_now(),
                "python": sys.version.replace("\n", " "),
                "platform": platform.platform(),
                "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
            }
        ]
    )
    atomic_write_frame(output / "PLAN13_ADAPTER_READY", adapter_ready)
    print(
        json.dumps(
            {
                "release_id": RELEASE_ID,
                "status": "validated_plan13_adapter_candidate",
                "gate_status": gate_expectation,
                "n_biological": len(biological_ids),
                "n_technical_spots": len(technical_ids),
                "n_contract_checks": len(checks),
                "source_dependence": SOURCE_DEPENDENCE,
                "cross_assay_comparable": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
