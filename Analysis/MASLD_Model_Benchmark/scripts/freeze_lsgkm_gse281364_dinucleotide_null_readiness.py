#!/usr/bin/env python3
"""Freeze an outcome-blind LS-GKM exact-dinucleotide-null redesign."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import random
from typing import Any, Iterable


DNA = frozenset("ACGT")
SEEDS = [1103, 2909, 4721, 6673, 8111]


class DinucleotideNullReadinessError(RuntimeError):
    """Raised when an authority or redesign invariant differs."""


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise DinucleotideNullReadinessError(f"JSON object differs: {path}")
    return value


def resolve_project_path(root: Path, relative_text: str) -> Path:
    relative = Path(relative_text)
    if relative.is_absolute() or ".." in relative.parts:
        raise DinucleotideNullReadinessError("unsafe project-relative path")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    return path


def manifest_members(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    members = manifest.get("artifacts")
    if not isinstance(members, list):
        raise DinucleotideNullReadinessError("artifact manifest differs")
    result: dict[str, dict[str, Any]] = {}
    for member in members:
        if not isinstance(member, dict) or not isinstance(member.get("path"), str):
            raise DinucleotideNullReadinessError("artifact member differs")
        if member["path"] in result:
            raise DinucleotideNullReadinessError("duplicate artifact member")
        result[member["path"]] = member
    return result


def bind_manifest(root: Path, binding: dict[str, Any]) -> tuple[Path, dict[str, Any]]:
    path = resolve_project_path(root, binding["path"])
    if path.name != "ARTIFACTS.json" or file_sha256(path) != binding["sha256"]:
        raise DinucleotideNullReadinessError(f"top artifact authority differs: {path}")
    manifest = load_json(path)
    if manifest.get("schema_version") != "masld-bench-artifacts-v1":
        raise DinucleotideNullReadinessError("artifact schema differs")
    return path, manifest


def require_member(
    manifest: dict[str, Any], member_path: str, expected_sha256: str
) -> dict[str, Any]:
    member = manifest_members(manifest).get(member_path)
    if member is None or member.get("sha256") != expected_sha256:
        raise DinucleotideNullReadinessError(f"artifact member differs: {member_path}")
    return member


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical_sequence(sequence: str) -> str:
    reverse = reverse_complement(sequence)
    return min(sequence, reverse)


def monomer_counts(sequence: str) -> Counter[str]:
    return Counter(sequence)


def dinucleotide_counts(sequence: str) -> Counter[str]:
    return Counter(sequence[index : index + 2] for index in range(len(sequence) - 1))


def _seed_from_parts(parts: Iterable[str | int]) -> int:
    material = "|".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(material).digest(), "big")


def exact_dinucleotide_shuffle(
    sequence: str,
    *,
    design_id: str,
    split_id: str,
    positive_id: str,
    model_seed: int,
    attempt: int,
) -> str:
    """Return one deterministic randomized Euler trail over the dinucleotide graph."""
    sequence = sequence.upper()
    if len(sequence) < 2 or set(sequence) - DNA:
        raise DinucleotideNullReadinessError("shuffle input must be unambiguous DNA")
    if model_seed not in SEEDS or attempt < 0:
        raise DinucleotideNullReadinessError("shuffle seed material differs")
    sequence_sha256 = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    rng = random.Random(
        _seed_from_parts(
            (design_id, split_id, positive_id, sequence_sha256, model_seed, attempt)
        )
    )
    outgoing: dict[str, list[str]] = defaultdict(list)
    for left, right in zip(sequence, sequence[1:]):
        outgoing[left].append(right)
    for edges in outgoing.values():
        rng.shuffle(edges)

    stack = [sequence[0]]
    trail: list[str] = []
    while stack:
        node = stack[-1]
        edges = outgoing[node]
        if edges:
            stack.append(edges.pop())
        else:
            trail.append(stack.pop())
    shuffled = "".join(reversed(trail))
    if (
        len(shuffled) != len(sequence)
        or shuffled[0] != sequence[0]
        or shuffled[-1] != sequence[-1]
        or monomer_counts(shuffled) != monomer_counts(sequence)
        or dinucleotide_counts(shuffled) != dinucleotide_counts(sequence)
    ):
        raise DinucleotideNullReadinessError("Euler-trail shuffle invariant failed")
    return shuffled


def five_seed_shuffles(
    sequence: str,
    *,
    design_id: str,
    split_id: str,
    positive_id: str,
    max_attempts: int,
) -> dict[int, str]:
    """Find one collision-free canonical shuffle for every frozen seed."""
    forbidden = {canonical_sequence(sequence)}
    result: dict[int, str] = {}
    for model_seed in SEEDS:
        for attempt in range(max_attempts):
            candidate = exact_dinucleotide_shuffle(
                sequence,
                design_id=design_id,
                split_id=split_id,
                positive_id=positive_id,
                model_seed=model_seed,
                attempt=attempt,
            )
            canonical = canonical_sequence(candidate)
            if canonical not in forbidden:
                result[model_seed] = candidate
                forbidden.add(canonical)
                break
        else:
            raise DinucleotideNullReadinessError(
                f"five-seed distinct shuffle gate failed: {positive_id} seed={model_seed}"
            )
    return result


def validate_config_boundary(config: dict[str, Any]) -> None:
    if (
        config.get("schema_version")
        != "masld-bench-lsgkm-gse281364-dinucleotide-null-readiness-input-v1"
        or config.get("status")
        != "freeze_outcome_blind_independent_redesign_readiness"
        or config.get("design_id")
        != "lsgkm_hepatocyte_atac_peak_exact_dinucleotide_null_v1"
    ):
        raise DinucleotideNullReadinessError("readiness config identity differs")
    relation = config["campaign_relation"]
    if (
        relation.get("independent_redesign") is not True
        or relation.get("failed_campaign_evidence_reused") is not False
        or relation.get("failed_campaign_thresholds_reinterpreted") is not False
        or relation.get("failed_campaign_matching_tolerances_reused") is not False
    ):
        raise DinucleotideNullReadinessError("independent redesign boundary differs")
    if any(value is not False for value in config["action_firewall"].values()):
        raise DinucleotideNullReadinessError("action firewall differs")

    runtime = config["runtime"]
    if (
        runtime.get("gkmtrain_argv")
        != ["-t", "2", "-l", "11", "-k", "7", "-d", "3", "-c", "1", "-e", "0.001", "-w", "1", "-m", "4096", "-T", "1", "-z"]
        or runtime.get("native_fit_seed_available") is not False
        or runtime.get("internal_cross_validation_enabled") is not False
        or runtime.get("resume_supported") is not False
    ):
        raise DinucleotideNullReadinessError("runtime contract differs")

    scoring = config["scoring_universe"]
    if (
        scoring.get("element_count") != 1033
        or scoring.get("long_range_block_count") != 239
        or scoring.get("outer_fold_count") != 5
        or scoring.get("fixed_seeds") != SEEDS
        or scoring.get("rows_per_readout") != 10330
        or scoring.get("context_specific_prediction") is not False
        or scoring.get("allele_score_sign") != "ALT_minus_REF"
    ):
        raise DinucleotideNullReadinessError("scoring universe differs")

    negative = config["negative_design"]
    exact_true = (
        "positive_roster_shared_across_seeds",
        "preserve_sequence_length_exactly",
        "preserve_monomer_counts_exactly",
        "preserve_dinucleotide_counts_exactly",
        "preserve_first_and_last_base",
        "reject_identity_or_reverse_complement_of_source",
        "reject_canonical_collision_with_any_fit_positive_or_negative",
        "reject_canonical_collision_with_held_scoring_alleles",
        "reject_ambiguous_or_incomplete_windows",
        "reject_blacklist_overlap",
    )
    if (
        negative.get("positive_count_per_fit") != 10000
        or negative.get("negative_count_per_fit") != 10000
        or negative.get("minimum_positive_count_per_fit") != 10000
        or negative.get("positive_class_weight") != 1
        or negative.get("max_shuffle_attempts_per_positive_seed") != 256
        or negative.get("required_distinct_canonical_shuffles_across_five_seeds") != 5
        or any(negative.get(field) is not True for field in exact_true)
        or negative.get("old_genomic_matching_thresholds_present") is not False
        or negative.get("old_minimum_positive_coverage_present") is not False
    ):
        raise DinucleotideNullReadinessError("negative design differs")
    if any(
        token in json.dumps(negative, sort_keys=True)
        for token in ("gc_tolerance", "mappability_tolerance", "repeat_tolerance", "atac_decile", "tss_distance_bin")
    ):
        raise DinucleotideNullReadinessError("old genomic-match tolerance leaked into redesign")

    gate = config["production_gate"]
    if (
        gate.get("current_state")
        != "ready_for_outcome_blind_input_materialization_and_scale_probe_only"
        or gate.get("require_exactly_10000_pairs_in_every_split_seed") is not True
        or gate.get("production_training_currently_authorized") is not False
    ):
        raise DinucleotideNullReadinessError("production gate differs")
    claims = config["claim_limits"]
    if claims.get("assay_native_sequence_transfer_baseline") is not True or any(
        claims.get(field) is not False
        for field in (
            "champion_eligible",
            "external_claim_eligible",
            "accessibility_classifier_claim_supported",
            "condition_specific_effect_supported",
            "donor_conditioning_supported",
            "target_gene_supported",
            "signed_eQTL_or_ieQTL_supported",
            "causal_claim_supported",
        )
    ):
        raise DinucleotideNullReadinessError("claim boundary differs")


def build_readiness(root: Path, config_path: Path) -> dict[str, Any]:
    config = load_json(config_path)
    validate_config_boundary(config)

    historical_path, historical = bind_manifest(root, config["campaign_relation"]["historical_boundary"])
    if historical.get("metadata", {}).get("status") != "terminally_blocked_300bp_matched_negative_design":
        raise DinucleotideNullReadinessError("historical terminal boundary differs")

    runtime_path, runtime_manifest = bind_manifest(root, config["runtime"])
    require_member(
        runtime_manifest,
        config["runtime"]["receipt_member"],
        config["runtime"]["receipt_sha256"],
    )
    runtime_receipt = load_json(runtime_path.parent / config["runtime"]["receipt_member"])
    if (
        runtime_receipt.get("status") != "pass_cpu_runtime_and_synthetic_fixture"
        or runtime_receipt.get("binary_sha256") != config["runtime"]["binary_sha256"]
        or runtime_receipt.get("reverse_complement_fixture_pass") is not True
        or runtime_receipt.get("direct_gkmsvm_direction_fixture_pass") is not True
        or runtime_receipt.get("deltasvm_sign_conversion_once") is not True
        or runtime_receipt.get("outcomes_read") is not False
    ):
        raise DinucleotideNullReadinessError("pinned runtime receipt differs")

    fasta_path, fasta = bind_manifest(root, config["reference_and_split"]["sequence_fasta_artifacts"])
    split_path, split = bind_manifest(root, config["reference_and_split"]["sequence_split_artifacts"])
    fasta_metadata = fasta.get("metadata", {})
    split_metadata = split.get("metadata", {})
    if (
        fasta_metadata.get("artifact_class") != "sequence_fasta"
        or fasta_metadata.get("build") != "GRCh38.p14"
        or fasta_metadata.get("indexed") is not True
        or fasta_metadata.get("pyfaidx_random_access_passed") is not True
        or split_metadata.get("artifact_class") != "sequence_split_contract"
        or split_metadata.get("dataset_id") != "gse296875"
        or split_metadata.get("genomic_folds") != 5
        or split_metadata.get("held_atac_used_for_window_selection") is not False
    ):
        raise DinucleotideNullReadinessError("reference or split authority differs")

    training_receipts = []
    for source in config["training_sources"]:
        manifest_path, manifest = bind_manifest(root, source)
        metadata = manifest.get("metadata", {})
        require_member(manifest, source["peak_member"], source["peak_member_sha256"])
        require_member(manifest, source["summary_member"], source["summary_member_sha256"])
        summary = load_json(manifest_path.parent / source["summary_member"])
        if (
            metadata.get("artifact_class") != "chrombpnet_training_inputs"
            or metadata.get("dataset_id") != "gse296875"
            or metadata.get("lineage_id") != "hepatocyte"
            or metadata.get("split_id") != source["split_id"]
            or metadata.get("held_atac_used_for_peak_discovery") is not False
            or metadata.get("test_outcomes_available_to_training") is not False
            or summary.get("retained") != source["expected_retained_peaks"]
            or summary.get("held_atac_test_used_for_peak_selection") is not False
            or summary.get("held_atac_valid_used_for_peak_selection") is not False
        ):
            raise DinucleotideNullReadinessError(f"training source differs: {source['split_id']}")
        training_receipts.append(
            {
                "split_id": source["split_id"],
                "artifacts_sha256": source["sha256"],
                "peak_member_sha256": source["peak_member_sha256"],
                "retained_training_peaks": summary["retained"],
                "held_atac_used_for_peak_selection": False,
            }
        )

    row_path, row_universe = bind_manifest(root, config["scoring_universe"]["row_universe_artifacts"])
    row_metadata = row_universe.get("metadata", {})
    if (
        row_metadata.get("dataset_id") != "gse281364"
        or row_metadata.get("source_locus_group_count") != 1033
        or row_metadata.get("long_range_block_count") != 239
        or row_metadata.get("fixed_seed_count") != 5
        or row_metadata.get("seeded_row_count") != 10330
        or row_metadata.get("outcomes_read") is not False
        or row_metadata.get("prediction_values_read") is not False
    ):
        raise DinucleotideNullReadinessError("outcome-blind row universe differs")
    allele_path, _allele_fixture = bind_manifest(root, config["scoring_universe"]["allele_fixture_artifacts"])

    control_path, control = bind_manifest(root, config["comparison"]["shared_control_artifacts"])
    require_member(
        control,
        config["comparison"]["control_prediction_member"],
        config["comparison"]["control_prediction_member_sha256"],
    )
    # Intentionally do not open the prediction member or any reporter-owned source.

    synthetic = ("ACGTTGCAAGTCGATCGGATCCGATGCTAGCTAGGCTA" * 8)[:300]
    shuffled = five_seed_shuffles(
        synthetic,
        design_id=config["design_id"],
        split_id="synthetic_split",
        positive_id="synthetic_positive",
        max_attempts=config["negative_design"]["max_shuffle_attempts_per_positive_seed"],
    )
    if len({canonical_sequence(value) for value in shuffled.values()}) != 5:
        raise DinucleotideNullReadinessError("synthetic five-seed diversity fixture differs")

    return {
        "schema_version": "masld-bench-lsgkm-gse281364-dinucleotide-null-readiness-v1",
        "status": "ready_for_outcome_blind_input_materialization_and_scale_probe_only",
        "design_id": config["design_id"],
        "scientific_estimand": config["negative_design"]["estimand"],
        "not_an_estimand": config["negative_design"]["not_an_estimand"],
        "campaign_relation": config["campaign_relation"],
        "runtime_artifacts_sha256": config["runtime"]["sha256"],
        "runtime_capabilities": {
            "cpu_runtime_passed": True,
            "balanced_separate_positive_negative_FASTA_supported": True,
            "fixed_length_300bp_supported": True,
            "native_fit_seed_available": False,
            "seed_role": "deterministic_negative_generation_and_record_order_only",
            "resume_supported": False,
            "shared_fit_readouts": config["readouts"]["shared_fit_readouts"],
            "independent_model_families": 1,
        },
        "reference_artifacts_sha256": file_sha256(fasta_path),
        "split_artifacts_sha256": file_sha256(split_path),
        "training_sources": training_receipts,
        "scoring_universe": {
            "row_universe_artifacts_sha256": file_sha256(row_path),
            "allele_fixture_artifacts_sha256": file_sha256(allele_path),
            "element_count": 1033,
            "long_range_block_count": 239,
            "fixed_seeds": SEEDS,
            "rows_per_readout": 10330,
            "context_specific_prediction": False,
            "allele_score_sign": "ALT_minus_REF",
        },
        "negative_design": config["negative_design"],
        "comparison": {
            "shared_control_artifacts_sha256": file_sha256(control_path),
            "control_prediction_member": config["comparison"]["control_prediction_member"],
            "control_prediction_member_sha256": config["comparison"]["control_prediction_member_sha256"],
            "control_model_id": config["comparison"]["control_model_id"],
            "control_head_id": config["comparison"]["control_head_id"],
            "canonical_control_rows_sha256": config["comparison"]["canonical_control_rows_sha256"],
            "prediction_member_opened": False,
        },
        "synthetic_shuffle_fixture": {
            "sequence_length": len(synthetic),
            "five_seed_canonical_shuffle_count": 5,
            "length_preserved": all(len(value) == len(synthetic) for value in shuffled.values()),
            "monomer_counts_preserved": all(monomer_counts(value) == monomer_counts(synthetic) for value in shuffled.values()),
            "dinucleotide_counts_preserved": all(dinucleotide_counts(value) == dinucleotide_counts(synthetic) for value in shuffled.values()),
        },
        "production_gate": config["production_gate"],
        "claim_limits": config["claim_limits"],
        "input_fastas_written": 0,
        "production_fits_executed": 0,
        "production_predictions_generated": 0,
        "benchmark_metrics_calculated": False,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "prediction_values_read": False,
        "sealed_assets_read": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise DinucleotideNullReadinessError("refusing to overwrite readiness receipt")
    readiness = build_readiness(
        arguments.project_root.resolve(strict=True),
        arguments.config.resolve(strict=True),
    )
    arguments.output.write_text(
        json.dumps(readiness, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(readiness, sort_keys=True))


if __name__ == "__main__":
    main()
