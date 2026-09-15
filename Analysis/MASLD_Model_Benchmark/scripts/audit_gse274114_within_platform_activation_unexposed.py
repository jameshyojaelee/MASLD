#!/usr/bin/env python3
"""Audit GSE274114 within-platform TaskSpecs without opening evaluator labels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.hashing import sha256_file
from masld_bench.registry import load_dataset_manifest, load_task_spec


FIXTURE_ARTIFACTS_SHA256 = "b334ff4fe77903be82bb228fc479a28d92d74f4cec0b06ea98c3d6a6100c9300"
PROMOTION_GATE_SHA256 = "2e45f65809ce927c5ef63e892b9434f41dd77bd2c28913815783c3d5a4577ab1"
SELECTIVELY_HASHED_FIXTURE_PATHS = frozenset(
    {
        "activation_audit.json",
        "authorities/exposure_audit.json",
        "authorities/instrument_confounding.json",
        "authorities/metadata.json",
        "authorities/qc.json",
        "authorities/reference.json",
        "authorities/rights.json",
        "authorities/topology.json",
        "model_inputs/gene_crosswalk.tsv",
        "model_inputs/quant_manifest.tsv",
        "source.sha256",
    }
)


class GSE274114WithinPlatformAuditError(ValueError):
    """Raised when a source requirement or within-platform TaskSpec differs."""


def require_fields(
    value: Mapping[str, Any], expected: Mapping[str, Any], label: str
) -> None:
    if any(value.get(key) != required for key, required in expected.items()):
        raise GSE274114WithinPlatformAuditError(f"{label} contract differs")


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise GSE274114WithinPlatformAuditError(f"JSON authority is not an object: {path}")
    return value


def audit_fixture_controls_without_protected_bytes(fixture_root: Path) -> dict[str, Any]:
    """Verify the pinned manifest and allowed members, never protected contents."""

    manifest_path = fixture_root / "ARTIFACTS.json"
    completion_path = fixture_root / "COMPLETE"
    if sha256_file(manifest_path) != FIXTURE_ARTIFACTS_SHA256:
        raise GSE274114WithinPlatformAuditError("fixture ARTIFACTS identity differs")
    manifest = load_json(manifest_path)
    completion = load_json(completion_path)
    if (
        set(manifest) != {"schema_version", "metadata", "artifacts"}
        or manifest.get("schema_version") != "masld-bench-artifacts-v1"
        or completion
        != {
            "schema_version": "masld-bench-complete-v1",
            "manifest_sha256": FIXTURE_ARTIFACTS_SHA256,
            "artifact_count": len(manifest.get("artifacts", [])),
        }
    ):
        raise GSE274114WithinPlatformAuditError("fixture freeze controls differ")
    require_fields(
        manifest.get("metadata", {}),
        {
            "artifact_class": "gse274114_no_fit_etiology_activation",
            "participants": 39,
            "sra_technical_runs": 71,
            "status": "passed_label_separated_no_fit_fixture",
            "four_class_performance_allowed": False,
            "MASH_vs_non_MASH_performance_allowed": False,
            "normalization_run": False,
            "model_fit_or_scoring_run": False,
        },
        "fixture metadata",
    )
    records: dict[str, dict[str, Any]] = {}
    for raw in manifest["artifacts"]:
        if (
            not isinstance(raw, dict)
            or set(raw) != {"path", "sha256", "size_bytes"}
            or raw["path"] in records
        ):
            raise GSE274114WithinPlatformAuditError("fixture artifact inventory differs")
        records[raw["path"]] = raw
    if set(SELECTIVELY_HASHED_FIXTURE_PATHS) - set(records):
        raise GSE274114WithinPlatformAuditError("fixture lacks an allowed audit member")
    protected = {
        path for path in records if path.startswith("evaluator_only/")
    }
    molecular = {
        path for path in records if path.startswith("model_inputs/quant_v49/")
    }
    if protected != {
        "evaluator_only/participant_folds.tsv",
        "evaluator_only/source_labels.tsv",
    } or len(molecular) != 39:
        raise GSE274114WithinPlatformAuditError("protected fixture census differs")
    for relative, record in records.items():
        path = fixture_root / relative
        if not path.is_file() or path.stat().st_size != record["size_bytes"]:
            raise GSE274114WithinPlatformAuditError("fixture member size or presence differs")
    for relative in SELECTIVELY_HASHED_FIXTURE_PATHS:
        if sha256_file(fixture_root / relative) != records[relative]["sha256"]:
            raise GSE274114WithinPlatformAuditError("allowed fixture member identity differs")
    return {
        "artifact_count": len(records),
        "evaluator_only_files_declared_not_opened": len(protected),
        "molecular_quant_files_declared_not_opened": len(molecular),
        "selectively_hashed_allowed_files": len(SELECTIVELY_HASHED_FIXTURE_PATHS),
    }


def validate_authorities(fixture_root: Path) -> dict[str, Any]:
    activation = load_json(fixture_root / "activation_audit.json")
    rights = load_json(fixture_root / "authorities/rights.json")
    reference = load_json(fixture_root / "authorities/reference.json")
    topology = load_json(fixture_root / "authorities/topology.json")
    instrument = load_json(fixture_root / "authorities/instrument_confounding.json")
    metadata = load_json(fixture_root / "authorities/metadata.json")
    qc = load_json(fixture_root / "authorities/qc.json")
    exposure = load_json(fixture_root / "authorities/exposure_audit.json")
    require_fields(
        activation,
        {
            "status": "pass_no_fit_fixture_two_within_instrument_contrasts_only",
            "participants": 39,
            "sra_runs": 71,
            "source_group_counts": {"CTRL": 9, "ENEG": 11, "ENEG_NASH": 9, "NASH": 10},
            "model_input_labels_separated": True,
            "normalization_run": False,
            "model_fit_or_scoring_run": False,
            "four_class_performance_allowed": False,
            "MASH_vs_non_MASH_performance_allowed": False,
            "HBV_only_is_MASLD_negative_control": False,
        },
        "activation",
    )
    require_fields(
        rights,
        {
            "access_tier": "public",
            "internal_nonclinical_research_use": True,
            "new_DUA_or_controlled_access_required": False,
            "redistribution": "prohibited_pending_source_specific_review",
            "derivative_weight_release": "requires_model_specific_terms_and_review",
        },
        "rights",
    )
    require_fields(
        reference,
        {
            "source_annotation": "GENCODE_v40",
            "target_annotation": "GENCODE_v49",
            "source_features": 61598,
            "mapped_source_rows": 60368,
            "unique_v49_target_genes": 60324,
            "par_y_x_two_row_target_genes": 44,
            "unmapped_source_rows": 1230,
            "rounding_allowed": False,
        },
        "reference",
    )
    require_fields(
        topology,
        {
            "participants": 39,
            "biopsies": 39,
            "geo_samples": 39,
            "biosamples": 39,
            "sra_experiments": 39,
            "sra_runs": 71,
            "pairing_topology": "same_study_unpaired",
            "technical_runs_are_biological_replicates": False,
        },
        "topology",
    )
    require_fields(
        instrument,
        {
            "perfect_MASH_superclass_instrument_confounding": True,
            "four_class_performance_allowed": False,
            "MASH_vs_non_MASH_performance_allowed": False,
            "OOD_accuracy_using_MASH_vs_healthy_allowed": False,
            "group_instrument": {
                "CTRL": "Illumina HiSeq 4000",
                "ENEG": "Illumina HiSeq 4000",
                "ENEG_NASH": "Illumina NovaSeq 6000",
                "NASH": "Illumina NovaSeq 6000",
            },
        },
        "instrument confounding",
    )
    require_fields(
        metadata,
        {
            "level": "source_group_aggregate_only",
            "participant_level_metadata_available": False,
            "participant_input_allowed": False,
        },
        "metadata",
    )
    require_fields(
        qc,
        {
            "quant_files": 39,
            "output_v49_quant_files": 39,
            "features_per_file": 61598,
            "output_v49_genes_per_file": 60324,
            "normalization_run": False,
            "model_fit_or_scoring_run": False,
        },
        "QC",
    )
    require_fields(
        exposure,
        {"project_exposure": "downstream_demo", "sealed_eligible": False},
        "exposure",
    )
    return {
        "participants": 39,
        "technical_runs": 71,
        "model_input_genes": 60324,
        "participant_metadata_inputs_available": False,
        "project_sealed": False,
    }


def validate_task_specs(
    dataset_path: Path,
    hiseq_task_path: Path,
    novaseq_task_path: Path,
    promotion_gate_path: Path,
) -> dict[str, Any]:
    dataset = load_dataset_manifest(dataset_path)
    hiseq = load_task_spec(hiseq_task_path)
    novaseq = load_task_spec(novaseq_task_path)
    gate = load_json(promotion_gate_path)
    if (
        dataset.dataset_id != "gse274114_mash_hbv"
        or dataset.biological_unit != "participant"
        or dataset.expected_biological_units != 39
        or dataset.exposure_status.value != "downstream_demo"
        or sha256_file(promotion_gate_path) != PROMOTION_GATE_SHA256
        or gate.get("champion_eligible") is not False
        or gate.get("four_class_performance_allowed") is not False
        or gate.get("mash_vs_non_mash_performance_allowed") is not False
        or gate.get("external_holdout_dataset_ids") != []
    ):
        raise GSE274114WithinPlatformAuditError("dataset or promotion gate differs")
    expected = {
        "gse274114_hiseq_healthy_vs_hbv": {
            "participants": 20,
            "counts": {"CTRL": 9, "ENEG": 11},
            "instrument": "Illumina HiSeq 4000",
        },
        "gse274114_novaseq_mash_vs_mash_hbv": {
            "participants": 19,
            "counts": {"NASH": 10, "ENEG_NASH": 9},
            "instrument": "Illumina NovaSeq 6000",
        },
    }
    tasks = {hiseq.task_id: hiseq, novaseq.task_id: novaseq}
    if set(tasks) != set(expected):
        raise GSE274114WithinPlatformAuditError("within-platform TaskSpec roster differs")
    for task_id, task in tasks.items():
        parameters = task.evaluator_parameters
        task_expected = expected[task_id]
        gates = " ".join(task.admission_gates).lower()
        claim = task.claim_gate.lower()
        if (
            task.status.value != "candidate"
            or task.unit_of_inference != "participant"
            or task.datasets_train != ("gse274114_mash_hbv",)
            or task.datasets_development
            or task.datasets_sealed
            or task.required_pairing_levels[0].value != "same_study_unpaired"
            or task.resampling_units != ("participant",)
            or task.bootstrap_replicates != 10000
            or parameters.get("participants") != task_expected["participants"]
            or parameters.get("source_group_counts") != task_expected["counts"]
            or parameters.get("instrument") != task_expected["instrument"]
            or parameters.get("outer_folds") != 5
            or parameters.get("technical_runs_are_replicates") is not False
            or parameters.get("evaluator_labels_available_to_model") is not False
            or parameters.get("four_class_performance_allowed") is not False
            or parameters.get("global_ood_accuracy_allowed") is not False
            or parameters.get("mash_vs_non_mash_performance_allowed") is not False
            or parameters.get("model_fit_or_scoring_run") is not False
            or "global four-class" not in gates
            or "cross-platform ood" not in claim
            or "champion" not in claim
        ):
            raise GSE274114WithinPlatformAuditError(f"TaskSpec differs: {task_id}")
    if sum(item["participants"] for item in expected.values()) != 39:
        raise GSE274114WithinPlatformAuditError("within-platform participant census differs")
    return {
        "task_ids": sorted(tasks),
        "hiseq_participants": 20,
        "novaseq_participants": 19,
        "outer_folds_per_task": 5,
        "global_four_class_performance_allowed": False,
        "global_ood_accuracy_allowed": False,
        "champion_eligible": False,
    }


def run(
    *,
    fixture_root: Path,
    dataset: Path,
    hiseq_task: Path,
    novaseq_task: Path,
    promotion_gate: Path,
    output: Path,
) -> None:
    if output.exists():
        raise GSE274114WithinPlatformAuditError("refusing to overwrite activation audit")
    fixture = audit_fixture_controls_without_protected_bytes(fixture_root)
    authorities = validate_authorities(fixture_root)
    tasks = validate_task_specs(dataset, hiseq_task, novaseq_task, promotion_gate)
    output.mkdir(mode=0o750)
    receipt = {
        "schema_version": "masld-bench-gse274114-within-platform-activation-audit-v1",
        "status": "passed_outcome_blind_within_platform_task_activation",
        "fixture": fixture,
        "authorities": authorities,
        "tasks": tasks,
        "fixture_artifacts_sha256": FIXTURE_ARTIFACTS_SHA256,
        "dataset_sha256": sha256_file(dataset),
        "hiseq_task_sha256": sha256_file(hiseq_task),
        "novaseq_task_sha256": sha256_file(novaseq_task),
        "promotion_gate_sha256": sha256_file(promotion_gate),
        "legacy_combined_task_activated_by_this_receipt": False,
        "global_four_class_task_activated": False,
        "global_ood_task_activated": False,
        "mash_vs_non_mash_task_activated": False,
        "evaluator_label_files_opened": False,
        "evaluator_fold_file_opened": False,
        "molecular_quant_values_opened": False,
        "participant_label_map_read": False,
        "model_fit_or_scoring_run": False,
        "sealed_or_champion_claim_eligible": False,
    }
    write_json_exclusive(output / "receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse274114_within_platform_activation_audit",
            "status": receipt["status"],
            "evaluator_labels_opened": False,
            "molecular_values_opened": False,
            "model_fit_or_scoring_run": False,
            "global_ood_task_activated": False,
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-root", required=True, type=Path)
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--hiseq-task", required=True, type=Path)
    parser.add_argument("--novaseq-task", required=True, type=Path)
    parser.add_argument("--promotion-gate", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(**vars(args))
    print(json.dumps({"output": str(args.output), "status": "passed"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
