#!/usr/bin/env python3
"""Freeze the model-specific GSE244832 observed-ATAC execution authority."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
from typing import Any, Mapping

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from scripts.materialize_gse244832_label_free_query_atac import tn5_positions


SCHEMA = "masld-bench-gse244832-observed-atac-execution-registration-v1"
ARTIFACT_CLASS = "gse244832_observed_atac_execution_registration"
DATASET_ID = "gse244832"
SOURCE_DATASET_ID = "gse296875"
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage")
ALLOWED_MODELS = (
    "lsi",
    "observed_atac_glm",
    "observed_atac_only",
    "peakvi",
)
ALLOWED_EXECUTION = {
    "lsi": "eligible_training_fitted_idf_svd_fixed_window_projection_only",
    "observed_atac_glm": "eligible_mandatory_fixed_window_baseline",
    "observed_atac_only": "eligible_mandatory_fixed_window_baseline",
    "peakvi": "eligible_source_fitted_fixed_window_decoder_without_target_adaptation",
}
ALLOWED_OUTPUT_FAMILIES = {
    model_id: ["masked_accessibility_count"] for model_id in ALLOWED_MODELS
}
QUERY_ROLES = {
    "valid": "test",
    "test": "valid",
}


class GSE244832ExecutionRegistrationError(RuntimeError):
    """Raised when the execution registration is not outcome-blind and exact."""


def _digest(path: Path) -> str:
    from hashlib import sha256

    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise GSE244832ExecutionRegistrationError(f"JSON object required: {path}")
    return value


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, fields: tuple[str, ...], rows: list[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _verify_authority(
    root: Path,
    spec: Mapping[str, Any],
    expected_class: str,
) -> Path:
    path = (root / str(spec.get("path", ""))).resolve(strict=True)
    manifest = verify_frozen_tree(path)
    if (
        path.is_symlink()
        or _digest(path / "ARTIFACTS.json") != spec.get("artifacts_sha256")
        or manifest.get("metadata", {}).get("artifact_class") != expected_class
        or manifest.get("metadata", {}).get("condition_values_read") is not False
        or manifest.get("metadata", {}).get("rna_assay_opened") is not False
        or manifest.get("metadata", {}).get("sequence_execution_authorized") is not False
    ):
        raise GSE244832ExecutionRegistrationError(f"authority differs: {path}")
    return path


def _read_axis_states(axis_root: Path) -> tuple[list[dict[str, Any]], Any]:
    import numpy as np

    rows = _read_tsv(axis_root / "donor_lineage_axis.tsv")
    if len(rows) != 90:
        raise GSE244832ExecutionRegistrationError("donor-lineage axis must be 18 x 5")
    result: list[dict[str, Any]] = []
    states = np.full((18, 4), 255, dtype=np.uint8)
    seen: set[tuple[int, int]] = set()
    for row in rows:
        if row["dataset_id"] != DATASET_ID:
            raise GSE244832ExecutionRegistrationError("axis dataset differs")
        donor_index = int(row["donor_index"])
        lineage_index = int(row["lineage_index"])
        if lineage_index == 4:
            if row["lineage_id"] != "t_cell" or row["evidence_state"] != "not_applicable":
                raise GSE244832ExecutionRegistrationError("T-cell missingness differs")
            continue
        if not (0 <= donor_index < 18 and 0 <= lineage_index < 4):
            raise GSE244832ExecutionRegistrationError("primary axis index differs")
        key = donor_index, lineage_index
        if key in seen or row["lineage_id"] != LINEAGES[lineage_index]:
            raise GSE244832ExecutionRegistrationError("primary axis duplicated or reordered")
        seen.add(key)
        cells = int(row["cells"])
        eligible = row["eligible_min_50_cells"] == "true"
        evidence = row["evidence_state"]
        if eligible and cells >= 50 and evidence == "observed":
            state = 0
        elif not eligible and 0 < cells < 50 and evidence == "below_qc":
            state = 3
        elif not eligible and cells == 0 and evidence == "structurally_missing":
            state = 1
        else:
            raise GSE244832ExecutionRegistrationError("axis missingness semantics differ")
        states[key] = state
        result.append(
            {
                "donor_index": donor_index,
                "donor_id": row["donor_id"],
                "lineage_index": lineage_index,
                "lineage_id": row["lineage_id"],
                "outer_fold": int(row["outer_fold"]),
                "cells": cells,
                "evidence_state": evidence,
                "missing_state_code": state,
                "eligible_min_50_cells": str(eligible).lower(),
            }
        )
    if len(seen) != 72 or np.any(states == 255) or int((states == 0).sum()) != 48:
        raise GSE244832ExecutionRegistrationError("18-donor primary missingness rectangle differs")
    result.sort(key=lambda row: (row["donor_index"], row["lineage_index"]))
    return result, states


def _read_axis_windows(axis_root: Path) -> dict[str, list[tuple[str, str, int, int]]]:
    rows = _read_tsv(axis_root / "windows.tsv")
    if len(rows) != 32000:
        raise GSE244832ExecutionRegistrationError("window axis must contain 32,000 rows")
    by_role = {role: [] for role in QUERY_ROLES}
    for row in rows:
        role = row["role"]
        if role not in by_role:
            raise GSE244832ExecutionRegistrationError("window role differs")
        start, end = int(row["start"]), int(row["end"])
        if end - start != 1000:
            raise GSE244832ExecutionRegistrationError("fixed-window width differs")
        by_role[role].append((row["window_id"], row["contig"], start, end))
    if any(len(rows_for_role) != 16000 for rows_for_role in by_role.values()):
        raise GSE244832ExecutionRegistrationError("window role cardinality differs")
    valid_contigs = {row[1] for row in by_role["valid"]}
    test_contigs = {row[1] for row in by_role["test"]}
    if not valid_contigs or not test_contigs or valid_contigs & test_contigs:
        raise GSE244832ExecutionRegistrationError("whole-contig target firewall differs")
    return by_role


def _read_query_windows(path: Path) -> list[tuple[str, str, int, int]]:
    rows = _read_tsv(path)
    result: list[tuple[str, str, int, int]] = []
    for expected, row in enumerate(rows):
        if int(row["context_index"]) != expected:
            raise GSE244832ExecutionRegistrationError("query window order differs")
        result.append((row["window_id"], row["contig"], int(row["start"]), int(row["end"])))
    return result


def _verify_query(
    query_root: Path,
    role: str,
    axis_rows: list[dict[str, Any]],
    unit_states: Any,
    axis_windows: Mapping[str, list[tuple[str, str, int, int]]],
) -> dict[str, Any]:
    import numpy as np

    contract = _load_json(query_root / "query_contract.json")
    authorities = _load_json(query_root / "input_authorities.json")
    if (
        contract.get("context_role") != role
        or contract.get("scored_target_role") != QUERY_ROLES[role]
        or contract.get("donors") != 18
        or contract.get("lineages") != list(LINEAGES)
        or contract.get("windows") != 16000
        or contract.get("profile_bins") != 20
        or contract.get("fragment_cut_sites") != "start_plus_4_and_end_minus_5"
        or contract.get("scored_target_contigs_exported_in_query_values") is not False
        or contract.get("rna_state") != "structurally_missing_for_atac_exchange"
        or authorities.get("condition_or_phenotype_authority_included") is not False
        or authorities.get("rna_or_cross_assay_authority_included") is not False
        or authorities.get("evaluator_outcome_authority_included") is not False
        or len(authorities.get("fragment_sha256_by_donor", {})) != 18
    ):
        raise GSE244832ExecutionRegistrationError("query contract or authority differs")
    query_axis = _read_tsv(query_root / "query_axis.tsv")
    expected_axis = [
        {
            "donor_index": str(row["donor_index"]),
            "donor_id": row["donor_id"],
            "lineage_index": str(row["lineage_index"]),
            "lineage_id": row["lineage_id"],
        }
        for row in axis_rows
    ]
    if query_axis != expected_axis:
        raise GSE244832ExecutionRegistrationError("query donor-lineage axis differs")
    query_windows = _read_query_windows(query_root / "query_windows.tsv")
    if query_windows != axis_windows[role]:
        raise GSE244832ExecutionRegistrationError("query windows differ from context role")
    target_contigs = {row[1] for row in axis_windows[QUERY_ROLES[role]]}
    if target_contigs & {row[1] for row in query_windows}:
        raise GSE244832ExecutionRegistrationError("scored target contig entered query values")
    profile = np.load(
        query_root / "query_unit_fragment_profile.uint32.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    count = np.load(
        query_root / "query_unit_fragment_total.uint32.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    if profile.shape != (18, 4, 16000, 20) or count.shape != (18, 4, 16000):
        raise GSE244832ExecutionRegistrationError("query array geometry differs")
    if profile.dtype != np.uint32 or count.dtype != np.uint32:
        raise GSE244832ExecutionRegistrationError("query array dtype differs")
    for donor in range(18):
        collapsed = profile[donor].sum(axis=2, dtype=np.uint64)
        if not np.array_equal(collapsed, count[donor]):
            raise GSE244832ExecutionRegistrationError("query profile/count conservation differs")
    structurally_missing = unit_states == 1
    if np.any(np.asarray(count)[structurally_missing] != 0):
        raise GSE244832ExecutionRegistrationError("structurally missing query unit is nonzero")
    return {
        "context_role": role,
        "scored_target_role": QUERY_ROLES[role],
        "profile_shape": list(profile.shape),
        "count_shape": list(count.shape),
        "query_windows": len(query_windows),
        "query_contigs": sorted({row[1] for row in query_windows}),
        "target_contigs_absent": True,
        "profile_count_conservation": True,
    }


def _execution_family(axis_root: Path) -> tuple[list[dict[str, Any]], list[str]]:
    rows = _read_tsv(axis_root / "family_eligibility.tsv")
    if len(rows) != 25 or len({row["model_id"] for row in rows}) != 25:
        raise GSE244832ExecutionRegistrationError("family roster differs")
    if not set(ALLOWED_MODELS) <= {row["model_id"] for row in rows}:
        raise GSE244832ExecutionRegistrationError("allowed family is absent")
    result: list[dict[str, Any]] = []
    blocked: list[str] = []
    for row in rows:
        model = row["model_id"]
        allowed = model in ALLOWED_MODELS
        if allowed and row["input_regime"] != "observed_atac":
            raise GSE244832ExecutionRegistrationError("allowed family is not observed ATAC")
        if not allowed:
            blocked.append(model)
        result.append(
            {
                "model_id": model,
                "input_regime": row["input_regime"],
                "output_family": row["output_family"],
                "axis_eligibility": row["eligibility"],
                "execution_eligibility": (
                    ALLOWED_EXECUTION[model]
                    if allowed
                    else "ineligible_model_specific_reference_motif_sequence_or_decoder_contract"
                ),
                "execution_allowed": str(allowed).lower(),
                "query_lineage_mask_required": "true",
            }
        )
    if tuple(sorted(model for model in ALLOWED_MODELS)) != tuple(sorted(ALLOWED_EXECUTION)):
        raise GSE244832ExecutionRegistrationError("allowed execution map differs")
    return result, sorted(blocked)


def validate_config(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("dataset_id") != DATASET_ID
        or config.get("source_dataset_id") != SOURCE_DATASET_ID
        or tuple(config.get("allowed_model_ids", ())) != ALLOWED_MODELS
        or config.get("direct_axis_family_table_execution_authorized") is not False
        or config.get("direct_query_artifact_execution_authorized") is not False
        or config.get("condition_values_read") is not False
        or config.get("phenotype_values_read") is not False
        or config.get("rna_assay_opened") is not False
        or config.get("cross_assay_join_inferred") is not False
        or config.get("biological_outcome_arrays_read") is not False
        or config.get("sequence_or_reference_features_authorized") is not False
        or config.get("metrics_calculated") is not False
        or config.get("champion_claim_allowed") is not False
    ):
        raise GSE244832ExecutionRegistrationError("execution registration config differs")
    paths: dict[str, Path] = {}
    implementation = config.get("implementation", {})
    for name in ("registration_builder", "registration_test", "commit", "commit_test", "scorer", "scorer_test"):
        path = (root / str(implementation.get(f"{name}_path", ""))).resolve(strict=True)
        if _digest(path) != implementation.get(f"{name}_sha256"):
            raise GSE244832ExecutionRegistrationError(f"implementation differs: {name}")
        paths[name] = path
    return paths


def register(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    import numpy as np

    if output.exists():
        raise GSE244832ExecutionRegistrationError("refusing to overwrite registration")
    root = root.resolve(strict=True)
    config_path = config_path.resolve(strict=True)
    config = _load_json(config_path)
    implementation = validate_config(root, config)
    axis_spec = config["axis_authority"]
    axis_root = _verify_authority(
        root,
        axis_spec,
        "gse244832_reference_guarded_atac_exchange_axis",
    )
    axis_rows, unit_states = _read_axis_states(axis_root)
    axis_windows = _read_axis_windows(axis_root)
    families, blocked = _execution_family(axis_root)
    queries: dict[str, dict[str, Any]] = {}
    for spec in config["query_authorities"]:
        role = str(spec["context_role"])
        if role not in QUERY_ROLES or role in queries:
            raise GSE244832ExecutionRegistrationError("query role roster differs")
        query_root = _verify_authority(
            root,
            spec,
            "gse244832_label_free_query_atac",
        )
        receipt = _verify_query(query_root, role, axis_rows, unit_states, axis_windows)
        receipt.update(
            {
                "path": str(query_root),
                "artifacts_sha256": spec["artifacts_sha256"],
            }
        )
        queries[role] = receipt
    if set(queries) != set(QUERY_ROLES) or tn5_positions(100, 120) != (104, 115):
        raise GSE244832ExecutionRegistrationError("query roles or Tn5 geometry differ")

    output.mkdir(parents=True, mode=0o750)
    np.save(output / "query_lineage_missing_state.uint8.npy", unit_states)
    _write_tsv(
        output / "query_axis_with_missingness.tsv",
        (
            "donor_index",
            "donor_id",
            "lineage_index",
            "lineage_id",
            "outer_fold",
            "cells",
            "evidence_state",
            "missing_state_code",
            "eligible_min_50_cells",
        ),
        axis_rows,
    )
    _write_tsv(
        output / "family_execution_eligibility.tsv",
        (
            "model_id",
            "input_regime",
            "output_family",
            "axis_eligibility",
            "execution_eligibility",
            "execution_allowed",
            "query_lineage_mask_required",
        ),
        families,
    )
    implementation_hashes = {
        name: {"path": str(path.relative_to(root)), "sha256": _digest(path)}
        for name, path in sorted(implementation.items())
    }
    contract = {
        "schema_version": "masld-bench-gse244832-observed-atac-execution-contract-v1",
        "dataset_id": DATASET_ID,
        "source_dataset_id": SOURCE_DATASET_ID,
        "axis_path": str(axis_root),
        "axis_artifacts_sha256": axis_spec["artifacts_sha256"],
        "query_authorities": queries,
        "allowed_model_ids": list(ALLOWED_MODELS),
        "allowed_output_families_by_model": ALLOWED_OUTPUT_FAMILIES,
        "blocked_model_ids": blocked,
        "model_manifest_requirements": {
            "model_id_matches_prediction_bundle": True,
            "registration_artifacts_sha256_required": True,
            "path_relative_to_frozen_prediction_root": True,
            "fit_scope": "gse296875_source_training_only",
            "sequence_features_consumed": False,
            "rna_assay_consumed": False,
            "target_dataset_fit_or_adaptation": False,
            "lsi_glm_atac_only_query_transform_scope": "source_training_only",
            "per_donor_context_only_if": (
                "deterministic_normalization_without_fitted_target_state"
            ),
        },
        "query_lineage_missing_state_path": "query_lineage_missing_state.uint8.npy",
        "query_lineage_missing_state_sha256": _digest(
            output / "query_lineage_missing_state.uint8.npy"
        ),
        "missingness": {
            "shape": [18, 4],
            "observed_eligible_code": 0,
            "structurally_missing_code": 1,
            "below_qc_code": 3,
            "eligible_units": int((unit_states == 0).sum()),
            "structurally_missing_units": int((unit_states == 1).sum()),
            "below_qc_units": int((unit_states == 3).sum()),
            "mask_before_any_query_transform": True,
            "missing_as_zero": False,
        },
        "target_contig_firewall": {
            "whole_contig_roles_disjoint": True,
            "scored_role_absent_from_query_values": True,
            "mask_before_query_transform": True,
        },
        "fragment_cut_sites": "start_plus_4_and_end_minus_5",
        "direct_axis_family_table_execution_authorized": False,
        "direct_query_artifact_execution_authorized": False,
        "query_values_reusable_only_through_this_registration": True,
        "prediction_commit_required_before_outcome_open": True,
        "condition_values_read": False,
        "phenotype_values_read": False,
        "rna_assay_opened": False,
        "cross_assay_join_inferred": False,
        "biological_outcome_arrays_read": False,
        "metrics_calculated": False,
        "sequence_execution_authorized": False,
        "external_or_sealed_evaluation": False,
        "champion_claim_allowed": False,
        "implementation": implementation_hashes,
    }
    write_json_exclusive(output / "execution_contract.json", contract)
    receipt = {
        "schema_version": "masld-bench-gse244832-observed-atac-registration-receipt-v1",
        "status": "passed",
        "donors": 18,
        "primary_lineages": 4,
        "eligible_units": int((unit_states == 0).sum()),
        "structurally_missing_units": int((unit_states == 1).sum()),
        "below_qc_units": int((unit_states == 3).sum()),
        "family_rows": len(families),
        "allowed_family_rows": len(ALLOWED_MODELS),
        "blocked_family_rows": len(blocked),
        "query_roles": sorted(queries),
        "windows_per_query_role": 16000,
        "profile_bins": 20,
        "tn5_fixture": {"start": 100, "end": 120, "cuts": [104, 115]},
        "condition_values_read": False,
        "biological_outcome_arrays_read": False,
        "prediction_files_read": False,
        "metrics_calculated": False,
        "champion_claim_allowed": False,
    }
    write_json_exclusive(output / "validation_receipt.json", receipt)
    (output / "python_version.txt").write_text(sys.version + "\n", encoding="utf-8")
    (output / "pip_freeze.txt").write_text(
        subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout,
        encoding="utf-8",
    )
    (output / "source.sha256").write_text(
        "".join(
            f"{_digest(path)}  {path}\n"
            for path in [config_path, *implementation.values()]
        ),
        encoding="utf-8",
    )
    freeze_tree(
        output,
        {
            "artifact_class": ARTIFACT_CLASS,
            "dataset_id": DATASET_ID,
            "source_dataset_id": SOURCE_DATASET_ID,
            "axis_artifacts_sha256": axis_spec["artifacts_sha256"],
            "allowed_model_ids": list(ALLOWED_MODELS),
            "query_lineage_mask_required": True,
            "condition_values_read": False,
            "rna_assay_opened": False,
            "cross_assay_join_inferred": False,
            "biological_outcome_arrays_read": False,
            "prediction_files_read": False,
            "metrics_calculated": False,
            "sequence_execution_authorized": False,
            "champion_claim_allowed": False,
            "status": "passed",
        },
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    print(json.dumps(register(arguments.root, arguments.config, arguments.output), sort_keys=True))


if __name__ == "__main__":
    main()
