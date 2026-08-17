#!/usr/bin/env python3
"""Tiny CPU smoke test of the scvi-tools 1.3.3 EWC adapter."""

from __future__ import annotations

import json
import anndata as ad
import numpy as np
import pandas as pd
import torch
from scvi import REGISTRY_KEYS
from scvi.data.fields import NumericalObsField
from scvi.model import SCANVI, SCVI

from masld_cl.bi_replay import compute_bi_scores
from masld_cl.ewc import EWCRegularizer, build_anchor_and_masks, compute_empirical_fisher
from masld_cl.distillation import (
    LatentDistillationRegularizer, posterior_means_at_anchor,
)
from masld_cl.scvi_adapter import (
    assert_scvi_version,
    attach_distillation_training_plan,
    attach_ewc_distillation_training_plan,
    attach_training_plan,
    scvi_reconstruction_loss,
)
from masld_cl.training import (
    _set_model_labels, _validate_scanvi_label_codes,
    _validate_update_trainability, set_all_seeds,
)


def make_adata(seed: int = 17) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    reference = rng.poisson(2.0, size=(32, 12)).astype(np.int32)
    query_control = rng.poisson(2.0, size=(16, 12)).astype(np.int32)
    query_case = rng.poisson(2.0, size=(16, 12)).astype(np.int32)
    query_case[:, :2] += rng.poisson(3.0, size=(16, 2)).astype(np.int32)
    values = np.vstack([reference, query_control, query_case])
    obs = pd.DataFrame({
        "batch": ["reference"] * 32 + ["query"] * 32,
        "audit_cell_type": (["A"] * 16 + ["pDCs"] * 16) * 2,
        "strict_reference": [True] * 32 + [False] * 32,
        "role": ["reference"] * 32 + ["control"] * 16 + ["case"] * 16,
    })
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    result = ad.AnnData(values, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(12)]))
    result.layers["counts"] = result.X.copy()
    _set_model_labels(result, query_unknown=True, unlabeled="Unknown")
    return result


def loader(model, adata, indices):
    registered = model._validate_anndata(adata)
    return model._make_data_loader(
        registered, indices=indices, batch_size=8, shuffle=False, drop_last=True
    )


def main() -> None:
    assert_scvi_version()
    set_all_seeds(17)
    full = make_adata()
    reference = full[full.obs["role"] == "reference"].copy()
    _set_model_labels(reference, query_unknown=False, unlabeled="Unknown")
    _validate_scanvi_label_codes(reference, unlabeled="Unknown")
    SCVI.setup_anndata(reference, layer="counts", batch_key="batch", labels_key="model_label")
    vae = SCVI(
        reference, n_hidden=16, n_latent=4, n_layers=1, dropout_rate=0,
        gene_likelihood="nb", dispersion="gene",
    )
    vae.train(max_epochs=2, train_size=1.0, batch_size=8, accelerator="cpu", devices=1)
    scanvi = SCANVI.from_scvi_model(vae, unlabeled_category="Unknown")
    scanvi.train(max_epochs=2, train_size=1.0, batch_size=8, accelerator="cpu", devices=1)
    reference_state = {
        name: value.detach().clone() for name, value in scanvi.module.named_parameters()
    }
    combined = full.copy()

    def load_query(seed):
        set_all_seeds(seed)
        model = SCANVI.load_query_data(
            combined.copy(), scanvi, unfrozen=True, accelerator="cpu", device=1
        )
        _validate_update_trainability(
            {
                name: bool(parameter.requires_grad)
                for name, parameter in model.module.named_parameters()
            },
            model_kind="all_lineage", method="continual_learning",
        )
        return model

    # Exercise lambda-zero equivalence on the actual scvi-tools adapter.
    ordinary = load_query(41)
    ordinary_repeat = load_query(41)
    zero_ewc = load_query(41)
    zero_state = {name: value.detach().clone() for name, value in zero_ewc.module.named_parameters()}
    zero_anchors, zero_masks = build_anchor_and_masks(zero_ewc.module, zero_state)
    ones = {name: torch.ones_like(value) for name, value in zero_anchors.items()}
    zero_regularizer = EWCRegularizer(zero_anchors, zero_masks, ones, ones)
    attach_training_plan(zero_ewc, semi_supervised=True)
    common_train = {
        "max_epochs": 1, "train_size": 1.0, "batch_size": 8,
        "accelerator": "cpu", "devices": 1,
    }
    set_all_seeds(73)
    ordinary.train(**common_train)
    set_all_seeds(73)
    ordinary_repeat.train(**common_train)
    set_all_seeds(73)
    zero_ewc.train(
        **common_train,
        plan_kwargs={"ewc_regularizer": zero_regularizer, "ewc_lambda": 0.0},
    )
    maximum_lambda_zero_difference = max(
        float(torch.max(torch.abs(left - right)).detach().cpu())
        for left, right in zip(ordinary.module.parameters(), zero_ewc.module.parameters())
    )
    if maximum_lambda_zero_difference > 1e-6:
        raise RuntimeError(
            f"lambda-zero adapter differs from ordinary fine tuning: {maximum_lambda_zero_difference}"
        )
    maximum_same_seed_parameter_difference = max(
        float(torch.max(torch.abs(left - right)).detach().cpu())
        for left, right in zip(
            ordinary.module.parameters(), ordinary_repeat.module.parameters()
        )
    )
    ordinary_latent = ordinary.get_latent_representation(combined)
    repeated_latent = ordinary_repeat.get_latent_representation(combined)
    maximum_same_seed_latent_difference = float(
        np.max(np.abs(ordinary_latent - repeated_latent), initial=0)
    )
    if max(
        maximum_same_seed_parameter_difference,
        maximum_same_seed_latent_difference,
    ) > 1e-6:
        raise RuntimeError(
            "same-seed CPU runs exceed 1e-6: "
            f"parameters={maximum_same_seed_parameter_difference}, "
            f"latent={maximum_same_seed_latent_difference}"
        )

    query = load_query(101)
    bi_scores = compute_bi_scores(
        query.module,
        loader(
            query, combined,
            np.flatnonzero(combined.obs["role"].to_numpy() == "reference"),
        ),
        seed=109, n_augmentations=200, mask_fraction=0.50,
    )
    if bi_scores.shape != (32,) or not np.isfinite(bi_scores).all() or np.any(bi_scores < 0):
        raise RuntimeError(f"actual scANVI BI adapter produced invalid scores: {bi_scores.shape}")
    anchors, masks = build_anchor_and_masks(query.module, reference_state)
    ref_indices = np.flatnonzero(combined.obs["role"].to_numpy() == "reference")
    control_indices = np.flatnonzero(combined.obs["role"].to_numpy() == "control")
    reference_fisher, reference_summary = compute_empirical_fisher(
        query.module, loader(query, combined, ref_indices), scvi_reconstruction_loss,
        expected_batch_size=8,
    )
    control_fisher, control_summary = compute_empirical_fisher(
        query.module, loader(query, combined, control_indices), scvi_reconstruction_loss,
        expected_batch_size=8,
    )
    regularizer = EWCRegularizer(anchors, masks, reference_fisher, control_fisher)
    regularizer.assert_zero_at_anchor(query.module)
    attach_training_plan(query, semi_supervised=True)
    query.train(
        max_epochs=1,
        train_size=1.0,
        batch_size=8,
        accelerator="cpu",
        devices=1,
        plan_kwargs={"ewc_regularizer": regularizer, "ewc_lambda": 0.5},
    )
    penalty = float(regularizer.penalty(query.module).detach().cpu())
    if not np.isfinite(penalty) or penalty <= 0:
        raise RuntimeError(f"EWC smoke penalty did not become positive: {penalty}")
    latent = query.get_latent_representation(combined)
    if latent.shape != (64, 4) or not np.isfinite(latent).all():
        raise RuntimeError(f"unexpected smoke latent: {latent.shape}")

    # Exercise exact replay-index transfer and the scANVI distillation plan.
    distillation = load_query(131)
    distillation.adata.obs["_masld_distill_index"] = np.arange(
        distillation.adata.n_obs, dtype=np.int64
    )
    distillation.get_anndata_manager(
        distillation.adata, required=True
    ).register_new_fields([
        NumericalObsField(REGISTRY_KEYS.INDICES_KEY, "_masld_distill_index")
    ])
    targets = posterior_means_at_anchor(
        distillation.module,
        distillation._make_data_loader(
            adata=distillation.adata, batch_size=8, shuffle=False, drop_last=False
        ),
    ).numpy().astype(np.float32, copy=False)
    replay_mask = np.asarray(
        distillation.adata.obs["role"].astype(str) == "reference", dtype=bool
    )
    distillation_regularizer = LatentDistillationRegularizer(
        torch.from_numpy(targets), torch.from_numpy(replay_mask)
    ).to(next(distillation.module.parameters()).device)
    distillation_loader = distillation._make_data_loader(
        adata=distillation.adata, batch_size=8, shuffle=False, drop_last=False
    )
    anchor_maximum = distillation_regularizer.assert_zero_at_anchor(
        distillation.module, distillation_loader
    )
    attach_distillation_training_plan(distillation, semi_supervised=True)
    set_all_seeds(137)
    distillation.train(
        **common_train,
        plan_kwargs={
            "distillation_regularizer": distillation_regularizer,
            "distillation_weight": 100.0,
        },
    )
    post_distillation_penalties = [
        float(distillation_regularizer.penalty(distillation.module, value).detach().cpu())
        for value in distillation._make_data_loader(
            adata=distillation.adata, batch_size=8, shuffle=False, drop_last=False
        )
    ]
    post_distillation_penalty = max(post_distillation_penalties)
    if not np.isfinite(post_distillation_penalty) or post_distillation_penalty <= 0:
        raise RuntimeError(
            f"distillation smoke penalty did not become positive: {post_distillation_penalty}"
        )

    # Exercise the combined V13 block-normalized EWC plus replay-distillation plan.
    combined_model = load_query(149)
    combined_model.adata.obs["_masld_distill_index"] = np.arange(
        combined_model.adata.n_obs, dtype=np.int64
    )
    combined_model.get_anndata_manager(
        combined_model.adata, required=True
    ).register_new_fields([
        NumericalObsField(REGISTRY_KEYS.INDICES_KEY, "_masld_distill_index")
    ])
    combined_targets = posterior_means_at_anchor(
        combined_model.module,
        combined_model._make_data_loader(
            adata=combined_model.adata, batch_size=8, shuffle=False, drop_last=False
        ),
    ).numpy().astype(np.float32, copy=False)
    combined_distillation = LatentDistillationRegularizer(
        torch.from_numpy(combined_targets), torch.from_numpy(replay_mask)
    ).to(next(combined_model.module.parameters()).device)
    combined_anchors, combined_masks = build_anchor_and_masks(
        combined_model.module, reference_state
    )
    combined_reference_fisher, _ = compute_empirical_fisher(
        combined_model.module,
        loader(combined_model, combined_model.adata, ref_indices),
        scvi_reconstruction_loss, expected_batch_size=8,
    )
    combined_control_fisher, _ = compute_empirical_fisher(
        combined_model.module,
        loader(combined_model, combined_model.adata, control_indices),
        scvi_reconstruction_loss, expected_batch_size=8,
    )
    for name, value in combined_masks.items():
        if name.startswith("classifier."):
            value.zero_()
    combined_ewc = EWCRegularizer(
        combined_anchors, combined_masks,
        combined_reference_fisher, combined_control_fisher,
        normalization="block_mean_product",
    )
    combined_ewc.assert_zero_at_anchor(combined_model.module)
    combined_distillation.assert_zero_at_anchor(
        combined_model.module,
        combined_model._make_data_loader(
            adata=combined_model.adata, batch_size=8, shuffle=False, drop_last=False
        ),
    )
    attach_ewc_distillation_training_plan(combined_model, semi_supervised=True)
    set_all_seeds(151)
    combined_model.train(
        **common_train,
        plan_kwargs={
            "ewc_regularizer": combined_ewc,
            "ewc_lambda": 1.0,
            "distillation_regularizer": combined_distillation,
            "distillation_weight": 100.0,
        },
    )
    combined_ewc_penalty = float(combined_ewc.penalty(combined_model.module).detach().cpu())
    combined_distillation_penalty = max(
        float(combined_distillation.penalty(combined_model.module, value).detach().cpu())
        for value in combined_model._make_data_loader(
            adata=combined_model.adata, batch_size=8, shuffle=False, drop_last=False
        )
    )
    if min(combined_ewc_penalty, combined_distillation_penalty) <= 0:
        raise RuntimeError("combined V13 penalties did not both become positive")
    print(json.dumps({
        "status": "pass",
        "post_step_ewc_penalty": penalty,
        "latent_shape": list(latent.shape),
        "reference_fisher_batches": reference_summary.n_batches,
        "control_fisher_batches": control_summary.n_batches,
        "maximum_lambda_zero_difference": maximum_lambda_zero_difference,
        "maximum_same_seed_parameter_difference": maximum_same_seed_parameter_difference,
        "maximum_same_seed_latent_difference": maximum_same_seed_latent_difference,
        "bi_score_count": int(len(bi_scores)),
        "bi_augmentations": 200,
        "bi_mask_fraction": 0.50,
        "cpu_determinism_tolerance": 1e-6,
        "distillation_zero_at_anchor_maximum": anchor_maximum,
        "post_step_distillation_penalty_maximum": post_distillation_penalty,
        "combined_v13_ewc_penalty": combined_ewc_penalty,
        "combined_v13_distillation_penalty_maximum": combined_distillation_penalty,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
