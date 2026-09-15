#!/usr/bin/env python3
"""Validate and freeze the restricted Enformer/Sei MPRA head TaskSpec."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-gse281364-enformer-sei-head-campaign-v1"
SEEDS = (1103, 2909, 4721, 6673, 8111)
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
CANDIDATES = (
    ("available_simple_controls", "zero"),
    ("available_simple_controls", "outer_training_mean"),
    ("available_simple_controls", "allele_identity_ridge"),
    ("enformer_crested_restricted_port", "rc_ensemble_sad_ridge"),
    ("enformer_crested_restricted_port", "rc_ensemble_sar_ridge"),
    ("sei", "native_40class_ridge"),
)
MISSING_BASELINES = (
    "deltaSVM",
    "gkm-SVM",
    "sequence_CNN",
    "sequence_transformer",
    "MPRALegNet_signed_rekeyed_to_1033_elements_239_blocks",
)
PREDICTION_FIELDS = (
    "seed",
    "row_hash",
    "unit_hash",
    "block_hash",
    "stratum",
    "outer_fold",
    "study_id",
    "assay_context_id",
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
    "model_id",
    "head_id",
    "prediction",
    "experimental_replicates",
    "biological_donors",
    "outcome_role",
)


class EnformerSeiTaskSpecError(RuntimeError):
    """Raised when the campaign requirement differs from the frozen design."""


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EnformerSeiTaskSpecError(f"invalid campaign config: {error}") from error
    if not isinstance(value, dict):
        raise EnformerSeiTaskSpecError("campaign config must be a JSON object")
    return value


def _binding_path(root: Path, binding: Mapping[str, Any], *, label: str) -> Path:
    value = binding.get("tree_path")
    if not isinstance(value, str) or not value:
        raise EnformerSeiTaskSpecError(f"{label} path missing")
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise EnformerSeiTaskSpecError(f"unsafe {label} path")
    try:
        tree = (root / relative).resolve(strict=True)
        tree.relative_to(root)
    except (OSError, ValueError) as error:
        raise EnformerSeiTaskSpecError(f"{label} path missing or escapes root") from error
    expected = binding.get("artifacts_sha256")
    if not isinstance(expected, str) or file_sha256(tree / "ARTIFACTS.json") != expected:
        raise EnformerSeiTaskSpecError(f"{label} ARTIFACTS identity differs")
    verify_frozen_tree(tree)
    return tree


def validate_config(config: Mapping[str, Any]) -> None:
    authorization = config.get("authorization", {})
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != "prespecified_restricted_development_head_campaign"
        or config.get("dataset_id") != "gse281364"
        or config.get("task_id") != "variant_to_regulation"
        or authorization.get("outcome_access") is not True
        or authorization.get("static_head_fit") is not True
        or any(
            authorization.get(field) is not False
            for field in (
                "base_checkpoint_tuning",
                "stack_fit",
                "residual_complementarity",
                "sealed_asset_access",
                "global_census_mutation",
            )
        )
    ):
        raise EnformerSeiTaskSpecError("campaign authorization boundary differs")
    task = config.get("task", {})
    if (
        task.get("role") != "exposed_development_mpra_proxy_only"
        or task.get("biological_donors") != 0
        or task.get("experimental_replicates_per_context") != 4
        or tuple(task.get("contexts", ())) != CONTEXTS
        or task.get("elements") != 1033
        or task.get("long_range_blocks") != 239
        or task.get("outer_folds") != 5
        or tuple(task.get("seeds", ())) != SEEDS
        or task.get("largest_receptive_field_buffer_bp") != 524288
        or task.get("prediction_rows_per_candidate") != 10330
    ):
        raise EnformerSeiTaskSpecError("task topology differs")
    fit = config.get("fit", {})
    if (
        fit.get("five_seed_strategy") != "outer_training_long_range_block_bootstrap"
        or fit.get("seeds_are_biological_replicates") is not False
        or fit.get("standardization") != "fit_on_seed_specific_training_bootstrap_only"
        or fit.get("feature_selection") != "absolute_training_Pearson_top_k_only"
        or fit.get("calibration") != "ridge_intercept_only_no_posthoc_calibration"
        or tuple(fit.get("alpha_grid", ()))
        != (0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)
        or tuple(fit.get("enformer_top_k_grid", ())) != (32, 128, 512, 2048, 5313)
        or tuple(fit.get("sei_top_k_grid", ())) != (40,)
        or tuple(fit.get("allele_top_k_grid", ())) != (16,)
        or fit.get("held_fold_outcomes_used_for_fit_or_tuning") is not False
        or fit.get("held_fold_features_used_for_preprocessing_selection_or_tuning") is not False
    ):
        raise EnformerSeiTaskSpecError("head-fit recipe differs")
    roster = tuple((row.get("model_id"), row.get("head_id")) for row in config.get("candidates", ()))
    if roster != CANDIDATES:
        raise EnformerSeiTaskSpecError("candidate roster differs")
    firewall = config.get("completion_firewall", {})
    if (
        firewall.get("mandatory_task_native_baselines_complete") is not False
        or tuple(firewall.get("missing_baselines", ())) != MISSING_BASELINES
        or firewall.get("family_native_results_reported_separately") is not True
        or any(
            firewall.get(field) is not False
            for field in (
                "models_ranked_as_interchangeable",
                "shortlist_authorized",
                "champion_claim_authorized",
                "native_equivalence_claim_authorized",
                "external_claim_authorized",
                "clinical_claim_authorized",
            )
        )
    ):
        raise EnformerSeiTaskSpecError("completion firewall differs")


def freeze_task_spec(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise EnformerSeiTaskSpecError("TaskSpec output already exists")
    config = load_config(config_path)
    validate_config(config)
    authorities = {}
    for name in ("static_authority", "outcome_authority", "allele_authority"):
        binding = config[name]
        tree = _binding_path(root, binding, label=name)
        authorities[name] = {
            "tree_path": binding["tree_path"],
            "artifacts_sha256": binding["artifacts_sha256"],
        }
        for member_key, hash_key in (
            ("base_member", "base_sha256"),
            ("row_member", "row_sha256"),
            ("contract_member", "contract_sha256"),
            ("member", "member_sha256"),
        ):
            if member_key not in binding:
                continue
            member = tree / binding[member_key]
            if file_sha256(member) != binding[hash_key]:
                raise EnformerSeiTaskSpecError(f"{name} member identity differs")
    output.mkdir(parents=True, mode=0o750)
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_prespecified_outcome_blind_taskspec",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "candidate_count": len(CANDIDATES),
        "prediction_rows_per_candidate": 10330,
        "fixed_seeds": list(SEEDS),
        "contexts": list(CONTEXTS),
        "elements": 1033,
        "long_range_blocks": 239,
        "config_sha256": file_sha256(config_path),
        "authorities": authorities,
        "outcomes_read": False,
        "prediction_values_read": False,
        "model_fit": False,
        "sealed_assets_read": False,
        "models_ranked": False,
        "champion_claim": False,
    }
    write_json_exclusive(output / "task_spec.json", receipt)
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    freeze_task_spec(args.root.resolve(strict=True), args.config.resolve(strict=True), args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
