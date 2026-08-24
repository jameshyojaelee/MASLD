#!/usr/bin/env python3
"""Build outcome-free identities for the five conditional-model ablations."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import tomllib

from masld_bench.artifacts import verify_frozen_tree
from masld_bench.context_preconditions import require_precondition_only_campaign
from masld_bench.hashing import canonical_sha256, sha256_file


ABLATION_NAMES = (
    "sequence_only",
    "trans_only",
    "permuted_context",
    "reverse_complement",
    "matched_shuffle",
)


def build(campaign_path: Path, context_fixture: Path) -> dict[str, object]:
    with campaign_path.open("rb") as handle:
        campaign = tomllib.load(handle)
    require_precondition_only_campaign(campaign)
    fixture = context_fixture.resolve(strict=True)
    manifest = verify_frozen_tree(fixture)
    expected_metadata = {
        "artifact_class": "context_conditioning_precondition_fixture",
        "status": "pass",
        "conditional_model_authorized": False,
        "architecture_built": False,
        "model_training_started": False,
        "sealed_outcomes_exposed": False,
    }
    if manifest.get("metadata") != expected_metadata:
        raise ValueError("context precondition artifact metadata differs")
    shared = {
        "status": "precondition_only_parent_unresolved",
        "parent_binding_required": True,
        "selected_parent_identity": None,
        "same_parent_rows_output_axis_splits_seeds_and_evaluator_required": True,
        "outcome_blind": True,
        "sealed_outcomes_allowed": False,
        "observed_query_atac_allowed": False,
        "independent_champion_identity": False,
        "executable": False,
        "model_training_allowed": False,
        "future_realization_requires_new_frozen_parent_bound_spec": True,
    }
    ablations = {
        "sequence_only": {
            **shared,
            "scientific_role": "parent_bound_cis_pathway_ablation",
            "preserve": [
                "parent_sequence_window_and_reference",
                "parent_sequence_backbone_and_head",
                "parent_output_crosswalk",
                "parent_row_split_and_seed_ids",
            ],
            "intervention": (
                "Disable every query-context-conditioned FiLM gate and LoRA delta. "
                "Do not read query RNA. Preserve the sequence path and compatible head; "
                "if the selected architecture cannot emit the registered outputs under "
                "that intervention, mark the ablation not_applicable rather than fitting "
                "a different hidden parent."
            ),
            "variant_sign": "ALT_minus_REF_after_forward_orientation_restoration",
            "missing_context": "not_read_not_zero_filled",
        },
        "trans_only": {
            **shared,
            "scientific_role": "parent_bound_trans_pathway_ablation",
            "preserve": [
                "parent_context_mapper_order_masks_and_context_encoder",
                "parent_output_crosswalk",
                "parent_row_split_and_seed_ids",
            ],
            "intervention": (
                "Remove cis sequence information at the named parent attachment point "
                "using a constant neutral latent fixed before outcomes. No genomic "
                "coordinate, REF, ALT, sequence-derived token, or sequence prototype may "
                "enter the trans path. If a full registered output cannot be emitted, "
                "mark not_applicable rather than add a new sequence-bearing head."
            ),
            "variant_sign": "exact_zero_ALT_minus_REF_by_construction",
            "neutral_latent_fit": "outer_training_partition_only_if_parent_requires_it",
        },
        "permuted_context": {
            **shared,
            "scientific_role": "parent_bound_context_assignment_negative_control",
            "preserve": [
                "parent_sequence_and_weights",
                "parent_context_values_masks_and_missingness_states",
                "parent_output_crosswalk",
                "parent_row_split_and_seed_ids",
            ],
            "intervention": (
                "Derange context row assignments with a precommitted seed inside the "
                "same dataset, outer partition, assay, lineage, pairing topology, and "
                "complete missingness-pattern stratum. Never move a context across donors' "
                "outer folds or use outcomes to construct or select a permutation."
            ),
            "self_assignment_allowed": False,
            "minimum_stratum_size": 2,
            "ineligible_small_stratum_state": "not_applicable",
            "permutation_receipt_required": [
                "ordered_row_ids",
                "stratum_ids",
                "seed",
                "source_to_destination_row_ids",
                "assignment_sha256",
            ],
        },
        "reverse_complement": {
            **shared,
            "scientific_role": "parent_bound_strand_equivariance_diagnostic",
            "preserve": [
                "parent_context_vector_and_mask",
                "parent_weights_and_heads",
                "parent_row_split_and_seed_ids",
            ],
            "intervention": (
                "Reverse-complement the complete DNA window, reverse output-bin order, "
                "swap every declared strand-paired output channel, restore genomic "
                "orientation, and compare or ensemble only after restoration. Context "
                "gene order is unchanged."
            ),
            "variant_sign": (
                "ALT_minus_REF_is_preserved_after_orientation_restoration; REF_ALT "
                "swap_must_negate"
            ),
            "required_future_fixture": [
                "sequence_round_trip",
                "output_bin_round_trip",
                "strand_channel_round_trip",
                "variant_sign_preservation",
                "ambiguous_base_mask_preservation",
            ],
        },
        "matched_shuffle": {
            **shared,
            "scientific_role": "parent_bound_cis_sequence_negative_control",
            "preserve": [
                "parent_context_vector_and_mask",
                "parent_window_length_and_output_bin_alignment",
                "parent_row_split_and_seed_ids",
                "ambiguous_base_positions",
            ],
            "intervention": (
                "Apply an outcome-blind deterministic dinucleotide-preserving shuffle "
                "independently within each 64-bp sequence bin and each contiguous ACGT "
                "run. Preserve N-mask positions and bin boundaries. For variant rows, "
                "keep the normalized REF or ALT span fixed and shuffle the remaining "
                "left and right run segments separately with the same seed contract."
            ),
            "matching_constraints": [
                "length",
                "per_bin_dinucleotide_counts_where_mathematically_feasible",
                "ambiguous_base_mask",
                "variant_anchor_and_allele",
                "output_bin_alignment",
            ],
            "failure_state": (
                "not_applicable_if_a_segment_has_no_nonidentity_valid_shuffle; never "
                "relax matching after outcomes"
            ),
            "shuffle_receipt_required": [
                "row_id",
                "seed",
                "segment_boundaries",
                "input_sequence_sha256",
                "shuffled_sequence_sha256",
                "matching_check_sha256",
            ],
        },
    }
    if tuple(ablations) != ABLATION_NAMES:
        raise ValueError("conditional ablation roster differs")
    identity = {
        "schema_version": "masld-bench-conditional-ablation-preconditions-v1",
        "status": "precondition_only_parent_unresolved",
        "ablation_names": list(ABLATION_NAMES),
        "conditional_campaign": {
            "campaign_id": campaign["campaign_id"],
            "status": campaign["status"],
            "submit_enabled": campaign["submit_enabled"],
            "conditional_model_spec_path": campaign["conditional_model_spec_path"],
            "campaign_sha256": sha256_file(campaign_path),
        },
        "context_precondition_binding": {
            "path": fixture.as_posix(),
            "manifest_sha256": sha256_file(fixture / "ARTIFACTS.json"),
            "fixture_sha256": sha256_file(fixture / "fixture.json"),
        },
        "architecture_count_cap": 1,
        "selected_parent_identity": None,
        "architecture_built": False,
        "model_training_started": False,
        "conditional_model_authorized": False,
        "sealed_outcomes_exposed": False,
        "ablations": ablations,
    }
    return {"precondition_spec_id": canonical_sha256(identity), **identity}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--context-fixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = build(args.campaign, args.context_fixture)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(args.output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o440)
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({
        "precondition_spec_id": payload["precondition_spec_id"],
        "ablation_names": payload["ablation_names"],
        "executable": False,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
