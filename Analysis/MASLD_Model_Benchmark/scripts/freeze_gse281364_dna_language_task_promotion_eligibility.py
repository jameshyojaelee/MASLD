#!/usr/bin/env python3
"""Freeze performance-blind task eligibility for two GSE281364 DNA LMs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping

from masld_bench.artifacts import (
    ArtifactError,
    freeze_tree,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
    write_text_exclusive,
)


SCHEMA = "masld-bench-gse281364-dna-language-task-promotion-eligibility-v1"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
MODELS = {"hyenadna", "nucleotide_transformer"}
REQUIRED_COMPARATORS = [
    "deltaSVM",
    "gkm-SVM",
    "sequence_CNN_control",
    "sequence_transformer_control",
    "rekeyed_signed_MPRALegNet",
]
FALSE_SCOPE_FIELDS = (
    "external_or_champion_role",
    "performance_read",
    "prediction_values_read",
    "reporter_outcomes_read",
    "sealed_assets_read",
    "metrics_calculated",
)
FORBIDDEN_AUTHORITY_FRAGMENTS = (
    "21099008",
    "/evaluation/",
    "/selection/",
    "gse289173",
    "sealed",
)


class EligibilityAuditError(RuntimeError):
    """Raised when the task-specific eligibility authority differs."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise EligibilityAuditError(f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise EligibilityAuditError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise EligibilityAuditError(f"{label} must be a JSON object")
    return value


def project_path(root: Path, relative: str) -> Path:
    value = Path(relative)
    if value.is_absolute() or not value.parts or ".." in value.parts:
        raise EligibilityAuditError(f"unsafe authority path: {relative}")
    lowered = f"/{value.as_posix().casefold()}/"
    if any(fragment in lowered for fragment in FORBIDDEN_AUTHORITY_FRAGMENTS):
        raise EligibilityAuditError(f"performance or sealed authority forbidden: {relative}")
    path = reject_symlink_components(root / value, label="eligibility authority")
    try:
        path.resolve(strict=True).relative_to(root)
    except (OSError, ValueError) as error:
        raise EligibilityAuditError(
            f"authority is missing or escapes project root: {relative}"
        ) from error
    return path


def validate_config(config: Mapping[str, Any]) -> None:
    if config.get("schema_version") != SCHEMA:
        raise EligibilityAuditError("eligibility schema differs")
    if config.get("status") != "performance_blind_task_specific_eligibility_frozen":
        raise EligibilityAuditError("eligibility status differs")

    scope = config.get("scope")
    if not isinstance(scope, dict):
        raise EligibilityAuditError("scope is missing")
    if scope.get("dataset_id") != "gse281364" or scope.get("task_id") != "mpra_allelic_direction":
        raise EligibilityAuditError("task scope differs")
    if set(scope.get("models", [])) != MODELS:
        raise EligibilityAuditError("model census differs")
    if scope.get("biological_donors") != 0 or scope.get("experimental_replicates_per_assay_context") != 4:
        raise EligibilityAuditError("replication-unit contract differs")
    if any(scope.get(field) is not False for field in FALSE_SCOPE_FIELDS):
        raise EligibilityAuditError("performance, outcome, sealed, or claim firewall differs")

    gate = config.get("mandatory_comparator_gate")
    if not isinstance(gate, dict) or gate.get("promotion_allowed") is not False:
        raise EligibilityAuditError("promotion must remain blocked")
    if gate.get("required_comparators_not_yet_complete_in_one_authoritative_same_universe_comparison") != REQUIRED_COMPARATORS:
        raise EligibilityAuditError("mandatory comparator census differs")

    models = config.get("models")
    if not isinstance(models, dict) or set(models) != MODELS:
        raise EligibilityAuditError("model eligibility census differs")

    hyena = models["hyenadna"]
    if hyena["identity"]["checkpoint_sha256"] != "deafb53209bfafb314d14bd81108546a34d473c07ed7ab7a355674376b56ad12":
        raise EligibilityAuditError("HyenaDNA checkpoint differs")
    if hyena["terms"]["code_license"] != "Apache-2.0" or hyena["terms"]["checkpoint_repository_declared_license"] != "BSD-3-Clause":
        raise EligibilityAuditError("HyenaDNA code/weight license split differs")
    if not hyena["terms"]["upstream_checkpoint_redistribution"].startswith("blocked_"):
        raise EligibilityAuditError("HyenaDNA redistribution must fail closed")
    if hyena["contamination"]["base_checkpoint_GSE281364_target_labels"] != "target_label_unexposed":
        raise EligibilityAuditError("HyenaDNA target-label exposure differs")
    if hyena["contamination"]["adapted_GSE281364_head"] != "continual_seen_training_folds_only":
        raise EligibilityAuditError("HyenaDNA adapted-head exposure differs")
    if hyena["action_eligibility"]["task_specific_tournament_promotion"] != "blocked_required_comparators_missing":
        raise EligibilityAuditError("HyenaDNA promotion gate differs")
    if not hyena["action_eligibility"]["open_champion"].startswith("conditionally_eligible_only_after_"):
        raise EligibilityAuditError("HyenaDNA open-champion status differs")

    nt = models["nucleotide_transformer"]
    if nt["identity"]["source_checkpoint_sha256"] != "971cb721bd8cc8134d3665b5d94b195fffb18a1480ad5925d6dcc9290bf1c17f":
        raise EligibilityAuditError("Nucleotide Transformer source checkpoint differs")
    if nt["identity"]["converted_checkpoint_sha256"] != "05f00d75d230a4e31bb0fc51565180355efd46c9f749aa9ae470c04290028136":
        raise EligibilityAuditError("Nucleotide Transformer conversion differs")
    if nt["terms"]["checkpoint_repository_declared_license"] != "CC-BY-NC-SA-4.0":
        raise EligibilityAuditError("Nucleotide Transformer weight license differs")
    if nt["terms"]["commercial_restriction"] is not True or nt["terms"]["share_alike"] is not True:
        raise EligibilityAuditError("Nucleotide Transformer restriction differs")
    if nt["declared_pretraining"]["exact_pinned_model_card_pretraining_tokens"] != 300_000_000_000:
        raise EligibilityAuditError("pinned Nucleotide Transformer 300B token declaration differs")
    if nt["declared_pretraining"]["prior_generic_manifest_pretraining_tokens"] != 50_000_000_000:
        raise EligibilityAuditError("Nucleotide Transformer conflicting prior token field differs")
    if nt["contamination"]["base_checkpoint_GSE281364_target_labels"] != "target_label_unexposed":
        raise EligibilityAuditError("Nucleotide Transformer target-label exposure differs")
    if nt["contamination"]["adapted_GSE281364_head"] != "continual_seen_training_folds_only":
        raise EligibilityAuditError("Nucleotide Transformer adapted-head exposure differs")
    if nt["action_eligibility"]["open_champion"] != "ineligible_under_current_terms":
        raise EligibilityAuditError("restricted comparator cannot be an open champion")
    if nt["action_eligibility"]["conditional_model_trigger_input"] != "ineligible_restricted_terms":
        raise EligibilityAuditError("restricted comparator cannot trigger the conditional model")

    firewall = config.get("firewall")
    if not isinstance(firewall, dict):
        raise EligibilityAuditError("firewall is missing")
    if firewall.get("eligibility_must_not_depend_on_performance") is not True:
        raise EligibilityAuditError("eligibility/performance separation differs")
    if firewall.get("no_ranking_or_shortlist_created") is not True:
        raise EligibilityAuditError("ranking firewall differs")
    forbidden = set(firewall.get("forbidden_input_classes", []))
    if not {"completed_head_predictions", "candidate_metrics", "GSE289173_features_or_labels"}.issubset(forbidden):
        raise EligibilityAuditError("forbidden input classes differ")

    authorities = config.get("frozen_authorities")
    if not isinstance(authorities, dict) or len(authorities) != 14:
        raise EligibilityAuditError("frozen authority census differs")
    for name, binding in authorities.items():
        if not isinstance(binding, dict):
            raise EligibilityAuditError(f"authority binding differs: {name}")
        if binding.get("kind") not in {"file", "tree"}:
            raise EligibilityAuditError(f"authority kind differs: {name}")
        if not SHA256.fullmatch(str(binding.get("sha256", ""))):
            raise EligibilityAuditError(f"authority SHA-256 differs: {name}")
        relative = str(binding.get("path", ""))
        lowered = f"/{Path(relative).as_posix().casefold()}/"
        if any(fragment in lowered for fragment in FORBIDDEN_AUTHORITY_FRAGMENTS):
            raise EligibilityAuditError(f"performance or sealed binding forbidden: {name}")


def verify_authorities(root: Path, config: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for name, binding in config["frozen_authorities"].items():
        path = project_path(root, str(binding["path"]))
        if binding["kind"] == "file":
            if path.is_symlink() or not path.is_file():
                raise EligibilityAuditError(f"authority is not a regular file: {name}")
            observed = sha256_file(path)
        else:
            try:
                verify_frozen_tree(path)
            except ArtifactError as error:
                raise EligibilityAuditError(f"frozen tree differs: {name}: {error}") from error
            observed = sha256_file(path / "ARTIFACTS.json")
        if observed != binding["sha256"]:
            raise EligibilityAuditError(f"frozen authority changed: {name}")
        records[name] = {
            "kind": binding["kind"],
            "path": str(path.relative_to(root)),
            "sha256": observed,
        }
    return records


def audit(*, root: Path, config_path: Path) -> dict[str, Any]:
    root = reject_symlink_components(root, label="project root").resolve(strict=True)
    config_path = reject_symlink_components(config_path, label="eligibility config")
    config = load_json(config_path, label="eligibility config")
    validate_config(config)
    authorities = verify_authorities(root, config)
    return {
        "schema_version": SCHEMA,
        "status": "pass_performance_blind_task_specific_eligibility",
        "dataset_id": "gse281364",
        "task_id": "mpra_allelic_direction",
        "biological_donors": 0,
        "experimental_replicates_per_assay_context": 4,
        "models": {
            model_id: {
                "action_eligibility": config["models"][model_id]["action_eligibility"],
                "terms": config["models"][model_id]["terms"],
                "contamination": config["models"][model_id]["contamination"],
                "declared_pretraining": config["models"][model_id]["declared_pretraining"],
            }
            for model_id in sorted(MODELS)
        },
        "mandatory_comparator_gate": config["mandatory_comparator_gate"],
        "frozen_authorities": authorities,
        "performance_read": False,
        "prediction_values_read": False,
        "reporter_outcomes_read": False,
        "sealed_assets_read": False,
        "metrics_calculated": False,
        "ranking_created": False,
        "shortlist_created": False,
        "promotion_executed": False,
    }


def eligibility_tsv(receipt: Mapping[str, Any]) -> str:
    fields = [
        "model_id",
        "benchmark_comparator",
        "development_result_reporting",
        "task_specific_tournament_promotion",
        "conditional_model_trigger_input",
        "project_sealed_evaluation",
        "open_champion",
        "universal_model_claim",
    ]
    rows = ["\t".join(fields)]
    for model_id, record in sorted(receipt["models"].items()):
        action = record["action_eligibility"]
        rows.append("\t".join([model_id] + [str(action[field]) for field in fields[1:]]))
    return "\n".join(rows) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    receipt = audit(root=args.root, config_path=args.config)
    output = reject_symlink_components(args.output, label="eligibility output")
    if output.exists():
        raise EligibilityAuditError(f"refusing to overwrite output: {output}")
    output.mkdir(parents=True, mode=0o750)
    write_json_exclusive(output / "receipt.json", receipt, mode=0o640)
    write_text_exclusive(output / "eligibility.tsv", eligibility_tsv(receipt), mode=0o640)
    manifest_sha256 = freeze_tree(
        output,
        {
            "artifact_class": "gse281364_dna_language_task_promotion_eligibility",
            "dataset_id": "gse281364",
            "task_id": "mpra_allelic_direction",
            "models": sorted(MODELS),
            "performance_read": False,
            "reporter_outcomes_read": False,
            "sealed_assets_read": False,
            "metrics_calculated": False,
            "promotion_executed": False,
            "status": "pass_performance_blind_task_specific_eligibility",
        },
    )
    verify_frozen_tree(output)
    print(json.dumps({"output": str(output), "artifacts_sha256": manifest_sha256, "status": receipt["status"]}, sort_keys=True))


if __name__ == "__main__":
    main()
