#!/usr/bin/env python3
"""Freeze outcome-blind authorization for the 24 LS-GKM fits after the scale check."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


SCHEMA = "masld-bench-lsgkm-gse281364-dinucleotide-full-campaign-input-v1"
SEEDS = [1103, 2909, 4721, 6673, 8111]
SPLIT_COUNTS = {
    "donor0_genomic0": 175,
    "donor1_genomic1": 254,
    "donor2_genomic2": 205,
    "donor3_genomic3": 154,
    "donor4_genomic4": 245,
}
REUSED_FIT = ("donor0_genomic0", 1103)


class FullCampaignAuthorizationError(RuntimeError):
    """Raised when an authority or campaign invariant differs."""


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise FullCampaignAuthorizationError(f"JSON object differs: {path}")
    return value


def resolve_file(root: Path, relative_text: str, expected_sha256: str) -> Path:
    relative = Path(relative_text)
    if relative.is_absolute() or ".." in relative.parts:
        raise FullCampaignAuthorizationError("unsafe project-relative file")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    if path.is_symlink() or not path.is_file() or file_sha256(path) != expected_sha256:
        raise FullCampaignAuthorizationError(f"file authority differs: {relative_text}")
    return path


def manifest_members(manifest: Mapping[str, Any]) -> dict[str, str]:
    members = manifest.get("artifacts")
    if not isinstance(members, list):
        raise FullCampaignAuthorizationError("artifact manifest differs")
    roster: dict[str, str] = {}
    for member in members:
        if not isinstance(member, dict) or not isinstance(member.get("path"), str):
            raise FullCampaignAuthorizationError("artifact member differs")
        if member["path"] in roster:
            raise FullCampaignAuthorizationError("duplicate artifact member")
        roster[member["path"]] = str(member.get("sha256"))
    return roster


def bind_member(
    root: Path,
    authority_path: Path,
    manifest: Mapping[str, Any],
    member_path: Path,
    expected_sha256: str,
) -> None:
    relative = member_path.relative_to(authority_path.parent).as_posix()
    if manifest_members(manifest).get(relative) != expected_sha256:
        raise FullCampaignAuthorizationError(f"member absent from authority: {relative}")


def expected_run_fits(config: Mapping[str, Any]) -> list[tuple[str, int]]:
    result: list[tuple[str, int]] = []
    for split_id, split in config["fit_grid"]["splits"].items():
        result.extend((split_id, int(seed)) for seed in split["run_seeds"])
    return result


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status")
        != "freeze_remaining_24_outcome_blind_fits_after_scale_pass"
        or config.get("design_id")
        != "lsgkm_hepatocyte_atac_peak_exact_dinucleotide_null_v1"
    ):
        raise FullCampaignAuthorizationError("campaign identity differs")
    runtime = config["runtime"]
    if (
        runtime.get("gkmtrain_argv")
        != ["-t", "2", "-l", "11", "-k", "7", "-d", "3", "-c", "1", "-e", "0.001", "-w", "1", "-m", "4096", "-T", "1", "-z"]
        or runtime.get("gkmpredict_argv") != ["-T", "1"]
        or runtime.get("native_fit_seed_available") is not False
        or runtime.get("resume_supported") is not False
    ):
        raise FullCampaignAuthorizationError("runtime differs")
    grid = config["fit_grid"]
    if (
        grid.get("fixed_seeds") != SEEDS
        or set(grid.get("splits", {})) != set(SPLIT_COUNTS)
        or grid.get("total_fit_count") != 25
        or grid.get("reused_fit_count") != 1
        or grid.get("remaining_fit_count") != 24
        or grid.get("positive_count_per_fit") != 10000
        or grid.get("negative_count_per_fit") != 10000
        or grid.get("sequence_length_bp") != 300
        or grid.get("model_state_id") != "hepatocyte"
    ):
        raise FullCampaignAuthorizationError("fit grid differs")
    all_fits = {(split_id, seed) for split_id in SPLIT_COUNTS for seed in SEEDS}
    run_fits = expected_run_fits(config)
    if (
        len(run_fits) != len(set(run_fits))
        or len(run_fits) != 24
        or set(run_fits) != all_fits - {REUSED_FIT}
        or any(
            config["fit_grid"]["splits"][split_id]["expected_held_elements"]
            != count
            for split_id, count in SPLIT_COUNTS.items()
        )
    ):
        raise FullCampaignAuthorizationError("remaining fit rectangle differs")
    scale = config["scale_probe"]
    if (
        (scale["reused_fit"]["split_id"], scale["reused_fit"]["model_seed"])
        != REUSED_FIT
        or scale.get("required_status")
        != "pass_one_fit_production_command_scale_probe"
    ):
        raise FullCampaignAuthorizationError("probe reuse contract differs")
    execution = config["execution"]
    if (
        execution.get("bundle_by_split") is not True
        or execution.get("bundle_count") != 5
        or execution.get("one_native_thread_per_fit") is not True
        or execution.get("requested_memory_gib_per_bundle") != 32
        or execution.get("requested_wall_seconds_per_bundle") != 14400
        or execution.get("requested_cpu_count_by_split")
        != {split_id: len(grid["splits"][split_id]["run_seeds"]) for split_id in SPLIT_COUNTS}
        or execution.get("maximum_requested_cpu_hours") != 96
        or execution.get("maximum_requested_memory_gib_hours") != 640
    ):
        raise FullCampaignAuthorizationError("execution contract differs")
    firewall = config["action_firewall"]
    if (
        firewall.get("remaining_24_fits_authorized") is not True
        or firewall.get("held_sequence_scoring_authorized") is not True
        or any(
            firewall.get(field) is not False
            for field in (
                "probe_fit_retraining_authorized",
                "outcome_access_authorized",
                "reporter_count_access_authorized",
                "sealed_asset_access_authorized",
                "control_prediction_value_access_authorized",
                "benchmark_metric_calculation_authorized",
                "aggregation_or_evaluation_authorized",
            )
        )
    ):
        raise FullCampaignAuthorizationError("action firewall differs")


def freeze_authorization(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    config = load_json(config_path)
    validate_config(config)
    if output.exists():
        raise FullCampaignAuthorizationError("authorization output exists")
    output.mkdir(parents=True, mode=0o750)

    inputs = config["inputs"]
    input_authority = resolve_file(root, inputs["artifacts_path"], inputs["artifacts_sha256"])
    input_manifest = load_json(input_authority)
    for path_key, sha_key in (
        ("independent_audit_path", "independent_audit_sha256"),
        ("materialization_contract_path", "materialization_contract_sha256"),
        ("scoring_manifest_path", "scoring_manifest_sha256"),
        ("scoring_fasta_path", "scoring_fasta_sha256"),
    ):
        member = resolve_file(root, inputs[path_key], inputs[sha_key])
        bind_member(root, input_authority, input_manifest, member, inputs[sha_key])
    audit = load_json(root / inputs["independent_audit_path"])
    materialization = load_json(root / inputs["materialization_contract_path"])
    if (
        audit.get("status") != "pass_independent_emitted_fasta_invariant_audit"
        or audit.get("fit_count") != 25
        or audit.get("audited_pair_rows") != 250000
        or audit.get("outcomes_read") is not False
        or materialization.get("shared_fit_count") != 25
        or materialization.get("production_fits_executed") != 0
        or materialization.get("outcomes_read") is not False
    ):
        raise FullCampaignAuthorizationError("input readiness differs")

    scale = config["scale_probe"]
    scale_authority = resolve_file(root, scale["artifacts_path"], scale["artifacts_sha256"])
    scale_manifest = load_json(scale_authority)
    for path_key, sha_key in (
        ("receipt_path", "receipt_sha256"),
        ("model_path", "model_sha256"),
        ("scores_path", "scores_sha256"),
    ):
        member = resolve_file(root, scale[path_key], scale[sha_key])
        bind_member(root, scale_authority, scale_manifest, member, scale[sha_key])
    receipt = load_json(root / scale["receipt_path"])
    if (
        receipt.get("status") != scale["required_status"]
        or receipt.get("production_fits_executed") != 1
        or receipt.get("active_fit") != {
            "split_id": "donor0_genomic0",
            "model_seed": 1103,
            "positive_count": 10000,
            "negative_count": 10000,
            "sequence_length_bp": 300,
            "positive_class_weight": 1,
            "model_state_id": "hepatocyte",
        }
        or receipt.get("model_receipt", {}).get("sha256") != scale["model_sha256"]
        or receipt.get("held_element_scores_sha256") != scale["scores_sha256"]
        or receipt.get("benchmark_metrics_calculated") is not False
        or any(
            receipt.get(field) is not False
            for field in (
                "outcomes_read",
                "reporter_counts_read",
                "control_prediction_values_read",
                "sealed_assets_read",
            )
        )
    ):
        raise FullCampaignAuthorizationError("scale gate differs")
    execution = config["execution"]
    scale_wall = float(receipt["total_native_wall_seconds"])
    scale_maxrss_gib = int(receipt["maximum_resident_set_kbytes"]) / 1024 / 1024
    if (
        execution["requested_wall_seconds_per_bundle"]
        < scale_wall * execution["scale_wall_buffer_multiplier_minimum"]
        or execution["requested_memory_gib_per_bundle"]
        < scale_maxrss_gib * execution["scale_maxrss_buffer_multiplier_minimum"]
    ):
        raise FullCampaignAuthorizationError("scale-derived resource buffer differs")

    result = {
        "schema_version": "masld-bench-lsgkm-gse281364-full-campaign-authorization-v1",
        "status": "authorized_remaining_24_outcome_blind_fits",
        "design_id": config["design_id"],
        "config_sha256": file_sha256(config_path),
        "runtime_artifacts_sha256": config["runtime"]["artifacts_sha256"],
        "input_artifacts_sha256": inputs["artifacts_sha256"],
        "scale_artifacts_sha256": scale["artifacts_sha256"],
        "scale_receipt_sha256": scale["receipt_sha256"],
        "reused_fit": scale["reused_fit"],
        "reused_model_sha256": scale["model_sha256"],
        "reused_scores_sha256": scale["scores_sha256"],
        "authorized_run_fits": [
            {"split_id": split_id, "model_seed": seed}
            for split_id, seed in expected_run_fits(config)
        ],
        "remaining_fit_count": 24,
        "bundle_count": 5,
        "scale_total_native_wall_seconds": scale_wall,
        "scale_maximum_resident_set_kbytes": receipt["maximum_resident_set_kbytes"],
        "execution": execution,
        "action_firewall": config["action_firewall"],
        "benchmark_metrics_calculated": False,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "control_prediction_values_read": False,
        "sealed_assets_read": False,
    }
    (output / "full_campaign_authorization.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    freeze_authorization(
        arguments.project_root.resolve(strict=True),
        arguments.config.resolve(strict=True),
        arguments.output,
    )


if __name__ == "__main__":
    main()
