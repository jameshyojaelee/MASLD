#!/usr/bin/env python3
"""Validate the native-task boundaries for local accessibility and variant models."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import tomllib
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import (
    ArtifactError,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
    write_text_exclusive,
)


SCHEMA = "masld-bench-local-variant-family-reconciliation-v1"
STATUS = "pass_reconciled_native_tasks_no_model_promoted"
MODELS = (
    "chrombpnet",
    "bpnet",
    "sequence_cnn_control",
    "sequence_transformer_control",
    "mpralegnet",
    "gkmsvm",
    "deltasvm",
    "abc",
    "re2g",
    "epinformer",
)
TASKS = (
    "local_profile_accessibility",
    "mpra_reporter_activity",
    "local_accessibility_variant_effect",
    "enhancer_gene_link",
    "observed_atac_native_link",
)
AUDIT_GROUPS = (
    "chrombpnet",
    "bpnet",
    "mpralegnet",
    "deltasvm",
    "gkmsvm",
    "abc",
    "re2g",
    "epinformer",
    "sequence_controls",
)
BASE_FILES = {
    "regulatory_local_registry",
    "variant_capability_registry",
    "variant_task",
}
EXPECTED_FILES = BASE_FILES | {
    f"{group}_{suffix}"
    for group in AUDIT_GROUPS
    for suffix in ("checkpoints", "crosswalk", "exposure")
}
EXPECTED_TREES = {
    "sequence_training_probe",
    "local_profile_one_split",
    "local_profile_incomplete_matrix",
    "bpnet_terminal_bundle",
    "sequence_cnn_terminal_bundle",
    "mpralegnet_admission",
    "mpralegnet_native_probe",
    "mpralegnet_development_evaluation",
}
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class LocalVariantReconciliationError(RuntimeError):
    """Raised when a model, endpoint, or read-only authority differs."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise LocalVariantReconciliationError(f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise LocalVariantReconciliationError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise LocalVariantReconciliationError(f"{label} must be a JSON object")
    return value


def safe_project_path(root: Path, relative: str) -> Path:
    value = Path(relative)
    if value.is_absolute() or not value.parts or ".." in value.parts:
        raise LocalVariantReconciliationError(f"unsafe project path: {relative}")
    try:
        path = reject_symlink_components(root / value, label="local variant authority")
        path.resolve(strict=True).relative_to(root)
    except (ArtifactError, OSError, ValueError) as error:
        raise LocalVariantReconciliationError(
            f"local variant authority is missing or escapes root: {relative}"
        ) from error
    return path


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != STATUS
        or config.get("family_id") != "local_accessibility_and_variant_effect"
        or tuple(config.get("model_order", ())) != MODELS
        or tuple(config.get("native_task_order", ())) != TASKS
        or set(config.get("models", {})) != set(MODELS)
        or set(config.get("native_tasks", {})) != set(TASKS)
    ):
        raise LocalVariantReconciliationError("family identity or roster differs")

    tasks = config["native_tasks"]
    if (
        tasks["local_profile_accessibility"].get("models")
        != list(MODELS[:4])
        or tasks["mpra_reporter_activity"].get("models") != ["mpralegnet"]
        or tasks["mpra_reporter_activity"].get("biological_unit")
        != "none_zero_donors"
        or tasks["mpra_reporter_activity"].get(
            "replicates_used_as_independent_donors"
        )
        is not False
        or tasks["enhancer_gene_link"].get("models")
        != ["abc", "re2g", "epinformer"]
        or tasks["observed_atac_native_link"].get("models") != ["abc", "re2g"]
        or tasks["local_profile_accessibility"].get("current_matrix_rankable")
        is not False
        or any(
            row.get("additional_development_screen_justified_now") is not False
            for row in tasks.values()
        )
    ):
        raise LocalVariantReconciliationError("native-task boundary differs")

    models = config["models"]
    if (
        models["epinformer"].get("code_license")
        != "NO_OPERATIVE_LICENSE_DETECTED_AT_SELECTED_COMMIT"
        or not models["epinformer"].get("runtime", "").startswith("blocked")
        or models["gkmsvm"].get("code_license") != "GPL-3.0"
        or models["deltasvm"].get("code_license") != "GPL-3.0"
        or models["abc"].get("exposure") != "reference_only"
        or models["re2g"].get("exposure") != "clean_declared"
        or not models["abc"].get("query_time_observed_atac")
        or not models["re2g"].get("query_time_observed_atac")
        or any(model.get("signed_gene_effect") for model in models.values())
        or any(model.get("open_champion_eligible_now") for model in models.values())
    ):
        raise LocalVariantReconciliationError("license, exposure, or endpoint differs")

    boundaries = config.get("family_boundaries", {})
    if not boundaries or any(value is not False for value in boundaries.values()):
        raise LocalVariantReconciliationError("claim or execution boundary differs")

    ask = config.get("bundled_cpu_validation_ask", {})
    if (
        ask.get("job_name") != "model-check-218"
        or ask.get("partition") != "cpu"
        or ask.get("qos") != "nslab"
        or ask.get("cpus_per_task") != 1
        or ask.get("memory_gb") != 8
        or ask.get("walltime") != "02:00:00"
        or ask.get("array") is not False
        or ask.get("job_count") != 1
        or ask.get("total_cpu_hours") != 2
        or ask.get("total_gpu_hours") != 0
        or ask.get("submitted") is not False
    ):
        raise LocalVariantReconciliationError("CPU validation ask differs")


def validate_authorities(
    root: Path, config: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    files = config.get("frozen_file_authorities", {})
    trees = config.get("frozen_tree_authorities", {})
    if set(files) != EXPECTED_FILES or set(trees) != EXPECTED_TREES:
        raise LocalVariantReconciliationError("frozen authority census differs")
    result: dict[str, dict[str, Any]] = {}
    for name, binding in files.items():
        path = safe_project_path(root, str(binding.get("path", "")))
        expected = str(binding.get("sha256", ""))
        if not SHA256.fullmatch(expected) or sha256_file(path) != expected:
            raise LocalVariantReconciliationError(f"frozen file changed: {name}")
        result[name] = {
            "kind": "file",
            "path": path.relative_to(root).as_posix(),
            "sha256": expected,
        }
    for name, binding in trees.items():
        path = safe_project_path(root, str(binding.get("path", "")))
        expected = str(binding.get("artifacts_sha256", ""))
        if not SHA256.fullmatch(expected):
            raise LocalVariantReconciliationError(f"tree hash differs: {name}")
        try:
            verify_frozen_tree(path)
        except ArtifactError as error:
            raise LocalVariantReconciliationError(
                f"frozen tree verification failed: {name}: {error}"
            ) from error
        if sha256_file(path / "ARTIFACTS.json") != expected:
            raise LocalVariantReconciliationError(f"frozen tree changed: {name}")
        result[name] = {
            "kind": "tree",
            "path": path.relative_to(root).as_posix(),
            "artifacts_sha256": expected,
        }
    return result


def validate_registry(root: Path, config: Mapping[str, Any]) -> list[dict[str, Any]]:
    binding = config["frozen_file_authorities"]["regulatory_local_registry"]
    path = root / binding["path"]
    try:
        registry = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise LocalVariantReconciliationError(
            f"invalid regulatory registry: {error}"
        ) from error
    wanted = {
        "chrombpnet": ("MIT", "target_label_unexposed", "candidate"),
        "bpnet": ("MIT", "target_label_unexposed", "candidate"),
        "mpralegnet": ("code_MIT_weights_MIT", "target_label_unexposed", "candidate"),
        "deltasvm": ("GPL-3.0", "target_label_unexposed", "candidate"),
        "gkmsvm": ("GPL-3.0", "target_label_unexposed", "candidate"),
        "abc": ("MIT", "reference_only", "candidate"),
        "re2g": ("MIT", "clean_declared", "candidate"),
        "epinformer": (
            "code_NO_LICENSE_DETECTED_weights_MIT",
            "target_label_unexposed",
            "blocked_terms",
        ),
    }
    rows = {
        row.get("model_id"): row
        for row in registry.get("models", ())
        if row.get("model_id") in wanted
    }
    if set(rows) != set(wanted):
        raise LocalVariantReconciliationError("regulatory registry roster differs")
    output: list[dict[str, Any]] = []
    for model_id, expected in wanted.items():
        row = rows[model_id]
        observed = (
            row.get("license_status"),
            row.get("exposure_status"),
            row.get("status"),
        )
        if observed != expected:
            raise LocalVariantReconciliationError(
                f"regulatory registry state differs: {model_id}"
            )
        output.append(
            {
                "model_id": model_id,
                "license_status": observed[0],
                "exposure_status": observed[1],
                "registry_status": observed[2],
                "admission_blocking": row.get("admission_blocking"),
            }
        )
    return output


def validate_capability_and_task(root: Path, config: Mapping[str, Any]) -> None:
    files = config["frozen_file_authorities"]
    capability_path = root / files["variant_capability_registry"]["path"]
    task_path = root / files["variant_task"]["path"]
    try:
        capability = tomllib.loads(capability_path.read_text(encoding="utf-8"))
        task = tomllib.loads(task_path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise LocalVariantReconciliationError(
            f"invalid variant capability or task contract: {error}"
        ) from error
    rows = {
        row.get("model_id"): row
        for row in capability.get("models", ())
        if row.get("model_id") in MODELS
    }
    if set(rows) != set(MODELS):
        raise LocalVariantReconciliationError("variant capability roster differs")
    expected_roles = {
        "chrombpnet": "baseline",
        "bpnet": "baseline",
        "sequence_cnn_control": "baseline",
        "sequence_transformer_control": "baseline",
        "mpralegnet": "secondary_only",
        "gkmsvm": "secondary_only",
        "deltasvm": "secondary_only",
        "abc": "link_only",
        "re2g": "link_only",
        "epinformer": "link_only",
    }
    if any(
        rows[model_id].get("role") != expected_roles[model_id]
        or rows[model_id].get("primary_eligible") is not False
        for model_id in MODELS
    ):
        raise LocalVariantReconciliationError("variant native role differs")
    mandatory = {
        "chrombpnet",
        "bpnet",
        "sequence_cnn_control",
        "sequence_transformer_control",
        "gkmsvm",
        "deltasvm",
        "abc",
        "re2g",
    }
    if {
        model_id
        for model_id, row in rows.items()
        if row.get("is_mandatory_baseline") is True
    } != mandatory:
        raise LocalVariantReconciliationError("mandatory baseline role differs")
    if (
        task.get("task_id") != "variant_to_regulation"
        or task.get("unit_of_inference") != "independent_locus"
        or task.get("split_id") != "locus_outer"
        or task.get("missingness_policy") != "explicit_state_and_mask"
        or task.get("uncertainty_method") != "paired_ld_block_bootstrap"
        or task.get("bootstrap_replicates") != 10000
        or "deltaSVM and direct gkm-SVM consume the same fold/state/seed LS-GKM weights"
        not in " ".join(task.get("admission_gates", ()))
        or any(
            model_id not in task.get("baseline_model_ids", ())
            for model_id in mandatory
        )
    ):
        raise LocalVariantReconciliationError("variant TaskSpec differs")


def validate_frozen_evidence(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    trees = config["frozen_tree_authorities"]
    matrix = load_json(
        root / trees["local_profile_incomplete_matrix"]["path"] / "matrix_summary.json",
        label="local-profile matrix summary",
    )
    one_split = load_json(
        root / trees["local_profile_one_split"]["path"] / "evaluation/evaluation.json",
        label="local-profile one-split evaluation",
    )
    mpr = load_json(
        root
        / trees["mpralegnet_development_evaluation"]["path"]
        / "evaluation/receipt.json",
        label="MPRALegNet development evaluation",
    )
    probe = load_json(
        root
        / trees["mpralegnet_native_probe"]["path"]
        / "predictions/probe_receipt.json",
        label="MPRALegNet runtime probe",
    )
    bpnet_manifest = load_json(
        root / trees["bpnet_terminal_bundle"]["path"] / "ARTIFACTS.json",
        label="BPNet terminal bundle",
    )
    cnn_manifest = load_json(
        root / trees["sequence_cnn_terminal_bundle"]["path"] / "ARTIFACTS.json",
        label="sequence CNN terminal bundle",
    )
    bpnet_bundle = bpnet_manifest.get("metadata", bpnet_manifest)
    cnn_bundle = cnn_manifest.get("metadata", cnn_manifest)
    if not isinstance(bpnet_bundle, Mapping) or not isinstance(cnn_bundle, Mapping):
        raise LocalVariantReconciliationError("terminal bundle metadata differs")
    if (
        matrix.get("dataset_id") != "gse296875"
        or matrix.get("evaluation_role") != "valid"
        or matrix.get("evaluation_groups") != 7
        or matrix.get("cross_fold_candidate_ranking_calculated") is not False
        or matrix.get("incomplete_coverage_used_for_comparison") is not False
        or matrix.get("test_outcomes_read") is not False
        or one_split.get("one_seed_one_split_screen") is not True
        or one_split.get("promotion_gate_evaluated") is not False
        or one_split.get("test_outcomes_read") is not False
    ):
        raise LocalVariantReconciliationError("local-profile evidence boundary differs")
    expected_profile = config["development_evidence_summary"]["local_profile_one_split"]
    for model, expected in expected_profile["relative_block_deviance_reduction"].items():
        observed = one_split["relative_block_deviance_reduction"][model]["valid"]
        if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-15):
            raise LocalVariantReconciliationError(f"profile metric differs: {model}")
    if (
        mpr.get("donor_count") != 0
        or mpr.get("experimental_replicates_per_context") != 4
        or mpr.get("replicates_used_as_independent_donors") is not False
        or mpr.get("models_ranked") is not False
        or mpr.get("champion_eligible") is not False
        or mpr.get("sealed_outcomes_read") is not False
        or probe.get("checkpoint_sha256")
        != "470dc7bfd3f0912c91f307a5f4019598072b8786294c46fc501b826c030cbe39"
        or probe.get("strict_checkpoint_restore") is not True
        or probe.get("deterministic_repeat_max_abs_diff") != 0.0
        or probe.get("sealed_outcomes_loaded") is not False
        or bpnet_bundle.get("logical_tasks") != 8
        or bpnet_bundle.get("terminal_status") != "passed"
        or cnn_bundle.get("logical_tasks") != 8
        or cnn_bundle.get("terminal_status") != "passed"
    ):
        raise LocalVariantReconciliationError("runtime or MPRA evidence differs")
    expected_mpr = config["development_evidence_summary"]["mpralegnet_descriptive_transfer"]
    summaries = mpr["metric_summary"]
    for context, pearson_key, spearman_key in (
        ("HepG2_control", "HepG2_control_pearson", "HepG2_control_spearman"),
        ("HepG2_PAOA", "HepG2_PAOA_pearson", "HepG2_PAOA_spearman"),
    ):
        row = summaries[f"mpralegnet_hepg2_test1_val2::{context}"]
        if not math.isclose(
            row["pearson"], expected_mpr[pearson_key], abs_tol=1e-15
        ):
            raise LocalVariantReconciliationError(f"MPRA Pearson differs: {context}")
        if not math.isclose(
            row["spearman"], expected_mpr[spearman_key], abs_tol=1e-15
        ):
            raise LocalVariantReconciliationError(f"MPRA Spearman differs: {context}")
    return {
        "local_profile_evaluation_groups": matrix["evaluation_groups"],
        "local_profile_matrix_rankable": False,
        "bpnet_frozen_logical_tasks": bpnet_bundle["logical_tasks"],
        "sequence_cnn_frozen_logical_tasks": cnn_bundle["logical_tasks"],
        "mpralegnet_elements": mpr["mpralegnet_elements"],
        "mpralegnet_locus_groups": mpr["locus_sequence_groups"],
        "mpralegnet_biological_donors": mpr["donor_count"],
    }


def tsv(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in fields})
    return buffer.getvalue()


def preflight(
    *, root: Path, config_path: Path, output: Path | None = None
) -> dict[str, Any]:
    root = root.resolve(strict=True)
    config = load_json(config_path, label="local variant reconciliation config")
    validate_config(config)
    authorities = validate_authorities(root, config)
    registry_rows = validate_registry(root, config)
    validate_capability_and_task(root, config)
    evidence = validate_frozen_evidence(root, config)
    task_rows = [
        {
            "task_id": task_id,
            "models": "|".join(config["native_tasks"][task_id]["models"]),
            "current_evidence": config["native_tasks"][task_id]["current_evidence"],
            "additional_screen_justified_now": config["native_tasks"][task_id][
                "additional_development_screen_justified_now"
            ],
            "next_action": config["native_tasks"][task_id]["next_action"],
        }
        for task_id in TASKS
    ]
    receipt = {
        "schema_version": "masld-bench-local-variant-family-reconciliation-receipt-v1",
        "status": STATUS,
        "family_id": config["family_id"],
        "model_rows": len(MODELS),
        "native_task_rows": len(TASKS),
        "authority_rows": len(authorities),
        "registry_rows_checked": len(registry_rows),
        "frozen_evidence": evidence,
        "additional_development_screens_allowed_now": [],
        "conditional_next_screens": [
            "gkmsvm_plus_deltasvm_after_ATAC_runtime_and_variant_fixture_activation",
            "ABC_plus_rE2G_after_observed_ATAC_link_task_and_denominator_activation",
            "additional_fixed_seeds_only_for_complete_local_profile_matrix_finalists",
        ],
        "active_campaign_outputs_bound_or_read": False,
        "sealed_assets_read": False,
        "test_outcomes_read": False,
        "new_training_scoring_or_model_import_performed": False,
        "primary_signed_eQTL_or_ieQTL_candidate_present": False,
        "champion_claim_allowed": False,
        "universal_model_claim_allowed": False,
        "clinical_claim_allowed": False,
        "bundled_cpu_validation_ask": config["bundled_cpu_validation_ask"],
    }
    if output is not None:
        output.mkdir(parents=True, exist_ok=False)
        write_json_exclusive(output / "receipt.json", receipt)
        write_text_exclusive(
            output / "native_task_dispositions.tsv",
            tsv(
                task_rows,
                (
                    "task_id",
                    "models",
                    "current_evidence",
                    "additional_screen_justified_now",
                    "next_action",
                ),
            ),
        )
        write_text_exclusive(
            output / "registry_census.tsv",
            tsv(
                registry_rows,
                (
                    "model_id",
                    "license_status",
                    "exposure_status",
                    "registry_status",
                    "admission_blocking",
                ),
            ),
        )
    return receipt


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    receipt = preflight(root=args.root, config_path=args.config, output=args.output)
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
