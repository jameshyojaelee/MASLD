"""Replay-plus-EWC expansion of strict7 with external healthy donors."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .ewc import EWCRegularizer, build_anchor_and_masks, compute_empirical_fisher
from .external_reference_common import (
    _enable_inference_without_optimization,
    load_external26_policy,
    verify_external26_common,
)
from .firewall import validate_program_firewall
from .sampling import concatenate_query_and_replay


class ExternalHealthyContinualError(ContractError):
    """Raised when the V22 external-healthy continual pilot changes identity."""


POLICY_SCHEMA = "masld-cl-external-healthy-continual-policy-v22"


def load_external_healthy_continual_policy(
    config: dict[str, Any], value: str | Path,
) -> tuple[Path, dict[str, Any], dict[str, Any]]:
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "external_healthy_continual_policy_v22.json"
    )
    if path != expected:
        raise ExternalHealthyContinualError("V22 requires its source-controlled policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("status") != "paper_setting_stop_go_pilot"
        or policy.get("case_stage_program_hero_gene_umap_cas13_used") is not False
    ):
        raise ExternalHealthyContinualError("V22 policy identity differs")
    decision_spec = policy["source_reference_decision"]
    decision_path = (path.parent / decision_spec["path"]).resolve()
    with decision_path.open() as handle:
        decision = json.load(handle)
    if (
        sha256_path(decision_path) != decision_spec["file_sha256"]
        or decision.get("lock_sha256") != decision_spec["lock_sha256"]
        or decision.get("decision") != decision_spec["required_decision"]
        or decision.get("selected_roster") != "common_strict7"
        or decision.get("expanded_reference_allowed_downstream") is not False
    ):
        raise ExternalHealthyContinualError("V22 source reference decision differs")
    return path, policy, decision


def _verify_strict7_source_model(
    config: dict[str, Any], decision: dict[str, Any], model_value: str | Path,
) -> tuple[Path, dict[str, Any]]:
    from .training import _model_bundle_identity

    model = Path(model_value).resolve()
    manifest_path = model.parent / "reference_manifest.json"
    expected = decision["selected_manifest"]
    if (
        str(manifest_path) != expected["path"]
        or sha256_path(manifest_path) != expected["sha256"]
    ):
        raise ExternalHealthyContinualError("V22 source reference manifest differs")
    with manifest_path.open() as handle:
        manifest = json.load(handle)
    observed = _model_bundle_identity(model)
    if (
        manifest.get("schema_version")
        != "masld-cl-external26-common-reference-v21"
        or manifest.get("config_sha256") != config["_config_sha256"]
        or manifest.get("roster_name") != "common_strict7"
        or manifest.get("model_kind") != "all_lineage"
        or manifest.get("model_bundle") != observed
    ):
        raise ExternalHealthyContinualError("V22 strict7 source model differs")
    return manifest_path, manifest


def train_external_healthy_continual_pilot(
    config: dict[str, Any], prepared_value: str | Path, lock_value: str | Path,
    common_policy_value: str | Path, policy_value: str | Path,
    reference_model_value: str | Path, output_value: str | Path,
    *, seed: int, replay_fraction: float, ewc_lambda: float,
) -> dict[str, Any]:
    """Run the V22 paper-setting stop/go pilot with no outcome access."""
    import anndata as ad
    import torch
    from scvi.model import SCANVI

    from .scvi_adapter import (
        assert_runtime_versions, attach_training_plan, scvi_reconstruction_loss,
    )
    from .training import (
        _capped_indices, _drop_unused_categories, _ensure_output,
        _export_embedding, _external_donor_split, _filter_fisher_from_eligibility_pool,
        _model_bundle_identity, _registered_loader,
        _sample_identity, _sample_replay, _set_model_labels, _train_kwargs,
        _trained_epochs, _validate_update_trainability, _write_environment,
        set_all_seeds,
    )

    assert_runtime_versions(config)
    set_all_seeds(seed)
    validate_program_firewall(config)
    common_lock = verify_external26_common(
        config, prepared_value, lock_value, common_policy_value
    )
    load_external26_policy(config, common_policy_value)
    policy_path, policy, decision = load_external_healthy_continual_policy(
        config, policy_value
    )
    pilot = policy["pilot"]
    if (
        seed != pilot["seed"]
        or replay_fraction != pilot["replay_fraction"]
        or ewc_lambda != pilot["ewc_lambda"]
        or pilot["maximum_epochs"] != config["architecture"]["all_lineage_max_epochs"]
        or pilot["early_stopping_patience"]
        != config["architecture"]["early_stopping_patience"]
        or pilot["control_fisher_fraction"]
        != config["sampling"]["control_fisher_fraction"]
    ):
        raise ExternalHealthyContinualError("V22 pilot hyperparameters differ")
    source_manifest_path, source_manifest = _verify_strict7_source_model(
        config, decision, reference_model_value
    )

    full = ad.read_h5ad(prepared_value)
    contract = full.uns.get("masld_external26_common", {})
    if (
        full.n_obs != 1374354 or full.n_vars != 4000
        or contract.get("common_gene_order_sha256") != common_lock["common_gene_order_sha256"]
        or contract.get("combined_obs_sha256") != common_lock["combined_obs_sha256"]
        or contract.get("combined_counts_sha256") != common_lock["combined_counts_sha256"]
    ):
        raise ExternalHealthyContinualError("V22 common-universe input differs")
    reference_indices = np.flatnonzero(
        full.obs["strict_reference"].astype(bool).to_numpy()
    )
    query_indices = np.flatnonzero(
        full.obs["external_reference_candidate"].astype(bool).to_numpy()
    )
    if (
        len(reference_indices) != 216957
        or full.obs.iloc[reference_indices]["donor_id"].nunique() != 7
        or len(query_indices) != 142036
        or full.obs.iloc[query_indices]["donor_id"].nunique() != 19
        or full.obs.iloc[query_indices]["analysis_eligible"].astype(bool).any()
    ):
        raise ExternalHealthyContinualError("V22 reference/query roster differs")
    reference_pool = _capped_indices(
        full, reference_indices, "all_lineage", config, seed + 1
    )
    query_pool = _capped_indices(
        full, query_indices, "all_lineage", config, seed + 2
    )
    replay = _sample_replay(
        full, reference_pool, replay_fraction, seed + 3
    )
    expected_replay = int(np.floor(len(reference_pool) * replay_fraction + 0.5))
    if len(replay) != expected_replay:
        raise ExternalHealthyContinualError("V22 replay size differs")
    training_indices = concatenate_query_and_replay(query_pool, replay)
    training = full[training_indices].copy()
    batch_key = config["features"]["selection_batch_key"]
    _drop_unused_categories(training, (batch_key,))
    _set_model_labels(
        training, query_unknown=True,
        unlabeled=config["features"]["unlabeled_category"],
    )

    reference = SCANVI.load(reference_model_value)
    reference_state = {
        name: value.detach().clone()
        for name, value in reference.module.named_parameters()
    }
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    model = SCANVI.load_query_data(
        training, reference, unfrozen=True,
        accelerator=accelerator, device=1,
    )
    trainable = {
        name: bool(parameter.requires_grad)
        for name, parameter in model.module.named_parameters()
    }
    _validate_update_trainability(
        trainable, model_kind="all_lineage", method="continual_learning"
    )
    fisher_reference = _filter_fisher_from_eligibility_pool(
        full, replay, reference_pool,
        config["sampling"]["lineage_fisher_min_cells"],
    )
    # External cells are healthy controls but intentionally retain primary_query=False.
    # Sample them with the same donor/lineage/preparation implementation directly.
    from .training import _sample_fisher
    controls = _sample_fisher(
        full, query_pool, pilot["control_fisher_fraction"], seed + 5, config
    )
    control_eligible = _filter_fisher_from_eligibility_pool(
        full, query_pool, query_pool,
        config["sampling"]["lineage_fisher_min_cells"],
    )
    expected_controls = int(np.floor(
        len(control_eligible) * pilot["control_fisher_fraction"] + 0.5
    ))
    if len(controls) != expected_controls or not len(fisher_reference):
        raise ExternalHealthyContinualError("V22 Fisher roster differs")
    fisher_reference_adata = full[fisher_reference].copy()
    controls_adata = full[controls].copy()
    _drop_unused_categories(fisher_reference_adata, (batch_key,))
    _drop_unused_categories(controls_adata, (batch_key,))
    _set_model_labels(
        fisher_reference_adata, query_unknown=False,
        unlabeled=config["features"]["unlabeled_category"],
    )
    _set_model_labels(
        controls_adata, query_unknown=True,
        unlabeled=config["features"]["unlabeled_category"],
    )
    anchors, masks = build_anchor_and_masks(model.module, reference_state)
    batch_size = config["sampling"]["fisher_batch_size"]
    reference_fisher, reference_summary = compute_empirical_fisher(
        model.module,
        _registered_loader(model, fisher_reference_adata, batch_size, drop_last=False),
        scvi_reconstruction_loss, expected_batch_size=batch_size,
        include_partial=True,
    )
    control_fisher, control_summary = compute_empirical_fisher(
        model.module,
        _registered_loader(model, controls_adata, batch_size, drop_last=False),
        scvi_reconstruction_loss, expected_batch_size=batch_size,
        include_partial=True,
    )
    if (
        reference_summary.n_observations != len(fisher_reference)
        or control_summary.n_observations != len(controls)
    ):
        raise ExternalHealthyContinualError("V22 Fisher did not consume every cell")
    regularizer = EWCRegularizer(
        anchors, masks, reference_fisher, control_fisher
    )
    regularizer.assert_zero_at_anchor(model.module)
    novel_batches = sorted(
        set(full.obs.iloc[query_pool][batch_key].astype(str))
        - set(full.obs.iloc[reference_pool][batch_key].astype(str))
    )
    masked_elements = int(sum((value == 0).sum().item() for value in masks.values()))
    if not novel_batches or masked_elements == 0:
        raise ExternalHealthyContinualError("V22 novel batch slices are not unpenalized")
    attach_training_plan(model, semi_supervised=True)
    external_split = _external_donor_split(training, seed)
    kwargs = _train_kwargs(config, "all_lineage", external_split)
    kwargs["plan_kwargs"] = {
        "ewc_regularizer": regularizer, "ewc_lambda": ewc_lambda,
    }
    model.train(**kwargs)

    output = _ensure_output(output_value)
    model.save(output / "model", save_anndata=True)
    model_bundle = _model_bundle_identity(output / "model")
    np.save(output / "training_latent.npy", model.get_latent_representation(training))
    ewc_path = output / "ewc_anchor_mask_fishers.pt"
    torch.save({
        "anchors": {name: value.detach().cpu() for name, value in anchors.items()},
        "masks": {name: value.detach().cpu() for name, value in masks.items()},
        "reference_fisher": {
            name: value.detach().cpu() for name, value in reference_fisher.items()
        },
        "query_control_fisher": {
            name: value.detach().cpu() for name, value in control_fisher.items()
        },
    }, ewc_path)

    _set_model_labels(
        full, query_unknown=True,
        unlabeled=config["features"]["unlabeled_category"],
    )
    mapping_model = SCANVI.load_query_data(
        full, model, unfrozen=False, accelerator=accelerator, device=1
    )
    mapping_inference = _enable_inference_without_optimization(mapping_model)
    mapping_model.save(output / "mapping_model", save_anndata=False)
    embedding = _export_embedding(
        mapping_model, full, np.ones(full.n_obs, dtype=bool),
        "all_lineage", config, output,
    )
    result = {
        "schema_version": "masld-cl-external-healthy-continual-pilot-v22",
        "model_kind": "all_lineage",
        "method": "continual_learning",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "common_prepared_lock": {
            "path": str(Path(lock_value).resolve()), "sha256": sha256_path(lock_value),
            "lock_sha256": common_lock["lock_sha256"],
        },
        "source_reference": {
            "manifest": str(source_manifest_path),
            "manifest_sha256": sha256_path(source_manifest_path),
            "model": source_manifest["model_bundle"],
        },
        "seed": seed,
        "replay_fraction": replay_fraction,
        "ewc_lambda": ewc_lambda,
        "n_reference_pool_cells": int(len(reference_pool)),
        "n_external_query_pool_cells": int(len(query_pool)),
        "n_replay_cells": int(len(replay)),
        "n_training_cells": int(len(training_indices)),
        "samples": {
            "reference_pool": _sample_identity(full, reference_pool),
            "external_query_pool": _sample_identity(full, query_pool),
            "replay": _sample_identity(full, replay),
            "reference_fisher": _sample_identity(full, fisher_reference),
            "external_control_fisher": _sample_identity(full, controls),
        },
        "fisher": {
            "estimator": "empirical_diagonal_minibatch_mean_gradient",
            "evaluated_at_same_post_load_anchor": True,
            "reference": asdict(reference_summary),
            "external_healthy_control": asdict(control_summary),
            "reference_is_exact_replay_subset": bool(
                set(map(int, fisher_reference)).issubset(set(map(int, replay)))
            ),
            "novel_batch_categories": novel_batches,
            "masked_parameter_elements": masked_elements,
            "tensor_bundle": {"file": ewc_path.name, "sha256": sha256_path(ewc_path)},
        },
        "stopping_epoch": _trained_epochs(model),
        "model_bundle": model_bundle,
        "mapping_model_bundle": _model_bundle_identity(output / "mapping_model"),
        "mapping_inference": mapping_inference,
        "embedding": embedding,
        "query_labels_hidden": True,
        "case_stage_program_hero_gene_umap_cas13_used": False,
        "sensitivity_only": True,
    }
    write_json_exclusive(output / "continual_reference_manifest.json", result)
    _write_environment(output)
    return result


def _external_healthy_centroids(bundle):
    from .metrics import donor_centroids

    _, latent, cells = bundle
    reference = cells["strict_reference"].to_numpy(dtype=bool)
    external = cells["cell_id"].astype(str).str.startswith("HLiCA|").to_numpy()
    selected_mask = reference | external
    selected = cells.loc[selected_mask].reset_index(drop=True)
    selected["_external_healthy"] = selected["cell_id"].astype(str).str.startswith(
        "HLiCA|"
    )
    donors, values = donor_centroids(
        np.asarray(latent[selected_mask], dtype=float),
        selected["donor_id"].astype(str),
    )
    metadata = selected.drop_duplicates("donor_id").set_index("donor_id").loc[donors]
    is_reference = metadata["strict_reference"].to_numpy(dtype=bool)
    is_external = metadata["_external_healthy"].to_numpy(dtype=bool)
    if is_reference.sum() != 7 or is_external.sum() != 19:
        raise ExternalHealthyContinualError(
            "V22 external-alignment donor roster differs"
        )
    return (
        donors.astype(str), values, metadata["dataset"].astype(str).to_numpy(),
        is_reference, is_external,
    )


def _v22_gate_decision(
    policy: dict[str, Any], reference_f1: dict[str, float],
    jaccard_loss: float, query_macro_change: float,
    query_lineage_changes: dict[str, float],
    external_macro_change: float, external_lineage_changes: dict[str, float],
    primary_alignment: dict[str, float], external_alignment: dict[str, float],
) -> tuple[list[dict[str, Any]], str]:
    thresholds = policy["stop_go_gates"]
    observations = [
        (
            "strict7_macro_f1_change_ci_low", reference_f1["ci_low"],
            thresholds["strict7_macro_f1_change_ci_low_greater_than"],
            reference_f1["ci_low"]
            > thresholds["strict7_macro_f1_change_ci_low_greater_than"],
        ),
        (
            "strict7_neighborhood_jaccard_loss", jaccard_loss,
            thresholds["strict7_neighborhood_jaccard_loss_at_most"],
            jaccard_loss <= thresholds["strict7_neighborhood_jaccard_loss_at_most"],
        ),
        (
            "canonical_query_macro_f1_change", query_macro_change,
            thresholds["canonical_query_macro_f1_change_at_least"],
            query_macro_change >= thresholds["canonical_query_macro_f1_change_at_least"],
        ),
        (
            "canonical_query_worst_major_lineage_f1_change",
            min(query_lineage_changes.values()),
            thresholds["worst_major_lineage_f1_change_at_least"],
            min(query_lineage_changes.values())
            >= thresholds["worst_major_lineage_f1_change_at_least"],
        ),
        (
            "external_healthy_macro_f1_change", external_macro_change,
            thresholds["external_healthy_macro_f1_change_at_least"],
            external_macro_change
            >= thresholds["external_healthy_macro_f1_change_at_least"],
        ),
        (
            "external_healthy_worst_major_lineage_f1_change",
            min(external_lineage_changes.values()),
            thresholds["external_healthy_worst_major_lineage_f1_change_at_least"],
            min(external_lineage_changes.values())
            >= thresholds["external_healthy_worst_major_lineage_f1_change_at_least"],
        ),
        (
            "primary_query_control_alignment_improvement",
            primary_alignment["improvement"],
            thresholds["primary_query_control_alignment_improvement_greater_than"],
            primary_alignment["improvement"]
            > thresholds["primary_query_control_alignment_improvement_greater_than"]
            and primary_alignment["ci_low"]
            > thresholds["primary_query_control_alignment_ci_low_greater_than"],
        ),
        (
            "external_healthy_alignment_improvement",
            external_alignment["improvement"],
            thresholds["external_healthy_alignment_improvement_at_least"],
            external_alignment["improvement"]
            >= thresholds["external_healthy_alignment_improvement_at_least"]
            and external_alignment["ci_low"]
            > thresholds["external_healthy_alignment_ci_low_greater_than"],
        ),
    ]
    gates = [
        {"gate": name, "observed": observed, "threshold": threshold, "pass": passed}
        for name, observed, threshold, passed in observations
    ]
    decision = "pass_paper_setting_pilot" if all(x["pass"] for x in gates) else "reject_paper_setting_pilot"
    return gates, decision


def evaluate_external_healthy_continual_pilot(
    config: dict[str, Any], policy_value: str | Path,
    common_value: str | Path, candidate_value: str | Path,
    output_value: str | Path, anchor_policy_value: str | Path | None = None,
) -> dict[str, Any]:
    """Evaluate the V22 pilot without reading disease-state outcomes."""
    from .benchmark import _bundle_centroids, _label_scores, _paired_improvement
    from .control_evaluation import _neighbor_jaccard_loss, _reference_f1_change_ci
    from .embedding import load_embedding, matched_rows
    from .external_reference_evaluation import _held_external_scores, _load_owned_reference

    validate_program_firewall(config)
    policy_path, policy, decision_lock = load_external_healthy_continual_policy(
        config, policy_value
    )
    common, common_manifest_path, common_manifest = _load_owned_reference(
        config, common_value, "common_strict7"
    )
    candidate_path = Path(candidate_value).resolve()
    candidate_manifest_path = candidate_path.parent / (
        "continual_reference_manifest.json"
        if anchor_policy_value is None else "anchor_export_manifest.json"
    )
    with candidate_manifest_path.open() as handle:
        candidate_manifest = json.load(handle)
    candidate = load_embedding(candidate_path)
    if anchor_policy_value is None:
        valid_candidate = (
            candidate_manifest.get("schema_version")
            == "masld-cl-external-healthy-continual-pilot-v22"
            and candidate_manifest.get("config_sha256") == config["_config_sha256"]
            and candidate_manifest.get("policy", {}).get("sha256") == sha256_path(policy_path)
            and candidate_manifest.get("embedding") == candidate[0]
            and candidate_manifest.get("query_labels_hidden") is True
            and candidate_manifest.get("case_stage_program_hero_gene_umap_cas13_used") is False
            and candidate_manifest.get("common_prepared_lock", {}).get("lock_sha256")
            == common_manifest["common_prepared_lock"]["lock_sha256"]
        )
    else:
        with Path(anchor_policy_value).open() as handle:
            transform_schema = json.load(handle).get("schema_version")
        if transform_schema == "masld-cl-external-healthy-anchor-policy-v23":
            from .external_healthy_anchor import load_external_healthy_anchor_policy

            anchor_policy_path, _, _ = load_external_healthy_anchor_policy(
                config, anchor_policy_value
            )
            expected_schema = "masld-cl-external-healthy-anchor-export-v23"
            extra_valid = True
        elif transform_schema == "masld-cl-external-healthy-direct-freeze-policy-v24":
            from .external_healthy_direct_freeze import (
                load_external_healthy_direct_freeze_policy,
            )

            anchor_policy_path, _, _ = load_external_healthy_direct_freeze_policy(
                config, anchor_policy_value
            )
            expected_schema = "masld-cl-external-healthy-direct-freeze-export-v24"
            extra_valid = (
                candidate_manifest.get("query_coordinates_byte_identical_to_v22") is True
                and candidate_manifest.get("fitted_transform") is False
            )
        else:
            raise ExternalHealthyContinualError("unknown external-healthy transform policy")
        valid_candidate = (
            candidate_manifest.get("schema_version") == expected_schema
            and candidate_manifest.get("config_sha256") == config["_config_sha256"]
            and candidate_manifest.get("policy", {}).get("sha256")
            == sha256_path(anchor_policy_path)
            and candidate_manifest.get("embedding") == candidate[0]
            and candidate_manifest.get("strict7_coordinates_and_predictions_exact") is True
            and candidate_manifest.get("query_predictions_changed") is False
            and candidate_manifest.get("case_stage_program_hero_gene_umap_cas13_read") is False
            and Path(candidate_manifest.get("strict7_embedding", {}).get("path", "")).resolve()
            == Path(common_value).resolve()
            and extra_valid
        )
    if not valid_candidate:
        raise ExternalHealthyContinualError("external-healthy candidate provenance differs")
    left, right = matched_rows(common[2], candidate[2])
    if len(left) != len(common[2]) or len(right) != len(candidate[2]):
        raise ExternalHealthyContinualError("V22 candidate cell roster differs")
    common_latent = np.asarray(common[1])[left]
    common_cells = common[2].iloc[left].reset_index(drop=True)
    candidate_latent = np.asarray(candidate[1])[right]
    candidate_cells = candidate[2].iloc[right].reset_index(drop=True)
    if not np.array_equal(
        common_cells["audit_cell_type"].astype(str).to_numpy(),
        candidate_cells["audit_cell_type"].astype(str).to_numpy(),
    ):
        raise ExternalHealthyContinualError("V22 frozen audit labels changed")

    current7 = common_cells["strict_reference"].to_numpy(dtype=bool)
    common_current7 = common_cells.loc[current7].reset_index(drop=True)
    candidate_current7 = candidate_cells.loc[current7].reset_index(drop=True)
    reference_f1 = _reference_f1_change_ci(
        common_current7, candidate_current7,
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"],
    )
    jaccard_loss = _neighbor_jaccard_loss(
        common_latent[current7], candidate_latent[current7], common_current7,
        k=config["evaluation"]["reference_neighborhood_k"],
        maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
        seed=config["screen"]["seed"],
    )
    common_bundle = (common[0], common_latent, common_cells)
    candidate_bundle = (candidate[0], candidate_latent, candidate_cells)
    common_query_macro, common_query_lineages = _label_scores(
        common_bundle, config["lineages"]
    )
    candidate_query_macro, candidate_query_lineages = _label_scores(
        candidate_bundle, config["lineages"]
    )
    query_macro_change = candidate_query_macro - common_query_macro
    query_lineage_changes = {
        lineage: candidate_query_lineages[lineage] - common_query_lineages[lineage]
        for lineage in config["lineages"]
    }
    common_external = _held_external_scores(common_cells, config["lineages"])
    candidate_external = _held_external_scores(candidate_cells, config["lineages"])
    external_macro_change = (
        candidate_external["mean_donor_macro_f1"]
        - common_external["mean_donor_macro_f1"]
    )
    external_lineage_changes = {
        lineage: (
            candidate_external["major_lineages"][lineage]["mean_donor_positive_f1"]
            - common_external["major_lineages"][lineage]["mean_donor_positive_f1"]
        )
        for lineage in config["lineages"]
    }
    primary_alignment = _paired_improvement(
        _bundle_centroids(candidate_bundle), _bundle_centroids(common_bundle),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 2200,
    )
    external_alignment = _paired_improvement(
        _external_healthy_centroids(candidate_bundle),
        _external_healthy_centroids(common_bundle),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 2201,
    )
    gates, decision = _v22_gate_decision(
        policy, reference_f1, jaccard_loss, query_macro_change,
        query_lineage_changes, external_macro_change, external_lineage_changes,
        primary_alignment, external_alignment,
    )
    result = {
        "schema_version": "masld-cl-external-healthy-continual-decision-v22",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "source_reference_decision_lock": {
            "path": str((
                policy_path.parent / policy["source_reference_decision"]["path"]
            ).resolve()),
            "file_sha256": policy["source_reference_decision"]["file_sha256"],
            "lock_sha256": decision_lock["lock_sha256"],
        },
        "sources": {
            "common_embedding": {"path": str(Path(common_value).resolve()), "sha256": sha256_path(common_value)},
            "common_manifest": {"path": str(common_manifest_path), "sha256": sha256_path(common_manifest_path)},
            "candidate_embedding": {"path": str(candidate_path), "sha256": sha256_path(candidate_path)},
            "candidate_manifest": {"path": str(candidate_manifest_path), "sha256": sha256_path(candidate_manifest_path)},
        },
        "control_only": True,
        "case_stage_program_hero_gene_umap_cas13_read": False,
        "anchored_export": anchor_policy_value is not None,
        "reference_retention": {
            "macro_f1_change": reference_f1,
            "neighborhood_jaccard_loss": jaccard_loss,
        },
        "canonical_query_label_concordance": {
            "macro_f1_change": query_macro_change,
            "major_lineage_changes": query_lineage_changes,
        },
        "external_healthy_label_transfer": {
            "common": common_external,
            "candidate": candidate_external,
            "macro_f1_change": external_macro_change,
            "major_lineage_changes": external_lineage_changes,
        },
        "primary_query_control_alignment": primary_alignment,
        "external_healthy_alignment": external_alignment,
        "gates": gates,
        "decision": decision,
        "pilot_pass": decision == "pass_paper_setting_pilot",
        "sensitivity_only": True,
    }
    write_json_exclusive(output_value, result)
    return result
