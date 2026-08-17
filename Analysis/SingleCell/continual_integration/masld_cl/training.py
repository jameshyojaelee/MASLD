"""Reference, baseline, and corrected continual-learning training workflows."""

from __future__ import annotations

import hashlib
import json
import os
import random
import subprocess
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import torch

from .config import write_json_exclusive
from .contracts import (
    ContractError, DeterministicGzipTextWriter, sha256_path, sha256_tree,
    verify_contract_lock,
)
from .data import verify_prepared, verify_prepared_file
from .ewc import EWCRegularizer, build_anchor_and_masks, compute_empirical_fisher
from .distillation import LatentDistillationRegularizer, posterior_means_at_anchor
from .firewall import validate_program_firewall
from .bi_replay import BI_MODES, compute_bi_scores, load_bi_replay_unlock
from .sampling import (
    balanced_exact_fraction_indices, capped_group_indices,
    concatenate_query_and_replay,
    donor_train_validation_split,
    stratified_exact_bi_replay,
)


class TrainingError(RuntimeError):
    pass


METHODS = {
    "architecture_surgery", "fine_tune", "replay_only", "ewc_only",
    "continual_learning", "latent_distillation", "continual_distillation",
}

DISTILLATION_INDEX_KEY = "_masld_distill_index"


def set_all_seeds(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True, warn_only=True)
    torch.backends.cudnn.benchmark = False
    try:
        import scvi
        scvi.settings.seed = seed
    except ImportError:
        pass


def _ensure_output(path: str | Path) -> Path:
    output = Path(path)
    output.mkdir(parents=True, exist_ok=False)
    return output


def _model_bundle_identity(model_dir: str | Path) -> dict[str, Any]:
    path = Path(model_dir).resolve()
    tree_hash, files = sha256_tree(path)
    return {
        "realpath": str(path),
        "tree_sha256": tree_hash,
        "files": files,
    }


def _verify_reference_model_bundle(
    model_dir: str | Path, config: dict[str, Any], model_kind: str,
) -> dict[str, Any]:
    path = Path(model_dir).resolve()
    manifests = [
        candidate for candidate in (
            path.parent / "reference_manifest.json",
            path.parent / "update_manifest.json",
        ) if candidate.is_file()
    ]
    if len(manifests) != 1:
        raise TrainingError(
            f"reference model requires exactly one sibling training manifest: {path}"
        )
    with manifests[0].open() as handle:
        manifest = json.load(handle)
    if manifest.get("config_sha256") != config["_config_sha256"]:
        raise TrainingError("reference model manifest has a different config")
    if manifest.get("model_kind") != model_kind:
        raise TrainingError("reference model kind differs from requested update")
    observed = _model_bundle_identity(path)
    expected = manifest.get("model_bundle", {})
    if (
        expected.get("realpath") != observed["realpath"]
        or expected.get("tree_sha256") != observed["tree_sha256"]
        or expected.get("files") != observed["files"]
    ):
        raise TrainingError("reference model bundle changed after training")
    return {
        "manifest_realpath": str(manifests[0].resolve()),
        "manifest_sha256": sha256_path(manifests[0]),
        "model_bundle": observed,
    }


def _model_mask(adata, model_kind: str, config: dict[str, Any]) -> np.ndarray:
    if model_kind == "all_lineage":
        return np.ones(adata.n_obs, dtype=bool)
    if model_kind not in config["lineages"]:
        raise TrainingError(f"model_kind must be all_lineage or one of {config['lineages']}")
    return adata.obs["audit_cell_type"].astype(str).to_numpy() == model_kind


def _capped_indices(adata, indices: np.ndarray, model_kind: str, config: dict[str, Any], seed: int) -> np.ndarray:
    obs = adata.obs.iloc[indices]
    if model_kind == "all_lineage":
        groups = list(zip(obs["donor_id"].astype(str), obs["audit_cell_type"].astype(str)))
        cap = config["sampling"]["all_lineage_cap_per_donor_label"]
    else:
        groups = [(x,) for x in obs["donor_id"].astype(str)]
        cap = config["sampling"]["lineage_cap_per_donor"]
    local = capped_group_indices(groups, cap, seed)
    return indices[local]


def _set_model_labels(adata, *, query_unknown: bool, unlabeled: str) -> None:
    labels = adata.obs["audit_cell_type"].astype(str).copy()
    if query_unknown:
        labels.loc[~adata.obs["strict_reference"].astype(bool)] = unlabeled
    # scANVI removes the unlabeled class from its classifier output, so its
    # registry code must be the final categorical code.
    import pandas as pd
    categories = sorted(set(labels) - {unlabeled}) + [unlabeled]
    adata.obs["model_label"] = pd.Categorical(labels, categories=categories)


def _validate_scanvi_label_codes(adata, *, unlabeled: str) -> None:
    labels = adata.obs.get("model_label")
    if labels is None or not hasattr(labels.dtype, "categories"):
        raise TrainingError("scANVI model_label must be categorical")
    categories = list(map(str, labels.cat.categories))
    if not categories or categories[-1] != unlabeled:
        raise TrainingError("scANVI unlabeled category must have the final registry code")
    codes = labels.cat.codes.to_numpy()
    labeled = labels.astype(str).to_numpy() != unlabeled
    if np.any(codes[labeled] < 0) or np.any(codes[labeled] >= len(categories) - 1):
        raise TrainingError("scANVI labeled target lies outside the classifier range")


def _validate_update_trainability(
    trainable_by_name: dict[str, bool], *, model_kind: str, method: str,
) -> None:
    """Distinguish scANVI's fixed class prior from unexpectedly frozen weights."""
    inherently_fixed = {"y_prior"} if model_kind == "all_lineage" else set()
    missing = sorted(inherently_fixed - set(trainable_by_name))
    incorrectly_trainable = sorted(
        name for name in inherently_fixed if trainable_by_name.get(name) is True
    )
    if missing or incorrectly_trainable:
        raise TrainingError(
            "scANVI fixed-parameter contract failed: "
            f"missing={missing}, trainable={incorrectly_trainable}"
        )
    optimizable = {
        name: value for name, value in trainable_by_name.items()
        if name not in inherently_fixed
    }
    if not optimizable:
        raise TrainingError("update model contains no optimizable parameters")
    if method == "architecture_surgery":
        if not any(optimizable.values()) or all(optimizable.values()):
            raise TrainingError(
                "architecture surgery did not preserve a frozen reference core"
            )
        return
    frozen = sorted(name for name, value in optimizable.items() if not value)
    if frozen:
        raise TrainingError(f"fully unfrozen update contains frozen parameters: {frozen}")


def _drop_unused_categories(adata, fields: tuple[str, ...]) -> None:
    """Prevent unobserved query categories from entering a reference registry."""
    for field in fields:
        if field in adata.obs and hasattr(adata.obs[field].dtype, "categories"):
            adata.obs[field] = adata.obs[field].cat.remove_unused_categories()


def _sample_identity(adata, indices: np.ndarray) -> dict[str, Any]:
    digest = hashlib.sha256()
    for value in adata.obs_names[indices].astype(str):
        digest.update(value.encode("utf-8"))
        digest.update(b"\0")
    obs = adata.obs.iloc[indices]
    return {
        "n_cells": int(len(indices)),
        "n_donors": int(obs["donor_id"].nunique()) if len(indices) else 0,
        "n_lineages": int(obs["audit_cell_type"].nunique()) if len(indices) else 0,
        "n_preparations": int(obs["preparation_method"].nunique()) if len(indices) else 0,
        "cell_order_sha256": digest.hexdigest(),
    }


def _validate_held_study_rosters(
    config: dict[str, Any], query_datasets: list[str] | None,
    control_fisher_datasets: list[str] | None,
    held_out_datasets: list[str] | None,
) -> set[str]:
    held_out = set(held_out_datasets or [])
    if not held_out:
        return held_out
    powered = set(config["evaluation"]["powered_query_studies"])
    if not held_out.issubset(powered):
        raise TrainingError("held-out datasets must be powered query studies")
    if query_datasets is None or held_out & set(query_datasets):
        raise TrainingError("held-study fitting requires an explicit disjoint query dataset roster")
    if control_fisher_datasets is None or held_out & set(control_fisher_datasets):
        raise TrainingError("held-study control Fisher must explicitly exclude held-out studies")
    return held_out


def _external_donor_split(adata, seed: int) -> list[np.ndarray]:
    train, validation = donor_train_validation_split(
        adata.obs["donor_id"].astype(str).to_numpy(), 0.20, seed
    )
    return [train, validation, np.asarray([], dtype=np.int64)]


def _train_kwargs(
    config: dict[str, Any], model_kind: str, external: list[np.ndarray],
    fixed_epochs: int | None = None,
) -> dict[str, Any]:
    common = {
        "max_epochs": fixed_epochs if fixed_epochs is not None else (
            config["architecture"]["all_lineage_max_epochs"]
            if model_kind == "all_lineage" else config["architecture"]["lineage_max_epochs"]
        ),
        "batch_size": config["architecture"]["batch_size"],
        "accelerator": "gpu" if torch.cuda.is_available() else "cpu",
        "devices": 1,
    }
    if fixed_epochs is not None:
        if fixed_epochs <= 0:
            raise TrainingError("fixed refit epochs must be positive")
        common.update({"early_stopping": False, "train_size": 1.0})
    else:
        common.update({
            "early_stopping": True,
            "early_stopping_patience": config["architecture"]["early_stopping_patience"],
            "early_stopping_monitor": "validation_loss",
            "datasplitter_kwargs": {"external_indexing": external, "drop_last": True},
        })
    return common


def _trained_epochs(model) -> int:
    history = getattr(model, "history_", None)
    if history is None or "elbo_train" not in history:
        raise TrainingError("scvi model did not expose an epoch history")
    epochs = len(history["elbo_train"])
    if epochs <= 0:
        raise TrainingError("scvi model completed zero training epochs")
    return int(epochs)


def _architecture(config: dict[str, Any]) -> dict[str, Any]:
    keys = ("n_hidden", "n_latent", "n_layers", "dropout_rate", "dispersion", "gene_likelihood")
    return {key: config["architecture"][key] for key in keys}


def train_reference(
    config: dict[str, Any], prepared_path: str | Path, contract_lock: str | Path,
    prepared_lock: str | Path, output: str | Path, model_kind: str, seed: int,
) -> dict[str, Any]:
    import anndata as ad
    import scvi
    from scvi.model import SCANVI, SCVI

    from .scvi_adapter import assert_runtime_versions

    assert_runtime_versions(config)
    set_all_seeds(seed)
    validate_program_firewall(config)
    lock = verify_contract_lock(config, contract_lock, full_hash=False)
    prepared_identity = verify_prepared_file(
        prepared_path, prepared_lock, config, lock
    )
    adata = ad.read_h5ad(prepared_path)
    verify_prepared(adata, config, lock, prepared_identity)
    mask = adata.obs["strict_reference"].astype(bool).to_numpy() & _model_mask(adata, model_kind, config)
    indices = _capped_indices(adata, np.flatnonzero(mask), model_kind, config, seed)
    reference = adata[indices].copy()
    _drop_unused_categories(reference, (config["features"]["selection_batch_key"],))
    _set_model_labels(reference, query_unknown=False, unlabeled=config["features"]["unlabeled_category"])
    if model_kind == "all_lineage":
        _validate_scanvi_label_codes(
            reference, unlabeled=config["features"]["unlabeled_category"]
        )
    batch_key = config["features"]["selection_batch_key"]
    labels_key = config["features"]["labels_key"]
    external = _external_donor_split(reference, seed)
    SCVI.setup_anndata(reference, layer="counts", batch_key=batch_key, labels_key=labels_key)
    vae = SCVI(reference, **_architecture(config))
    vae.train(**_train_kwargs(config, model_kind, external))
    if model_kind == "all_lineage":
        model = SCANVI.from_scvi_model(
            vae, unlabeled_category=config["features"]["unlabeled_category"]
        )
        model.train(**_train_kwargs(config, model_kind, external))
        model_class = "SCANVI"
    else:
        model = vae
        model_class = "SCVI"
    output = _ensure_output(output)
    model.save(output / "model", save_anndata=True)
    model_bundle = _model_bundle_identity(output / "model")
    np.save(output / "reference_latent.npy", model.get_latent_representation(reference))
    evaluation_mask = (
        adata.obs["strict_reference"].astype(bool).to_numpy()
        & _model_mask(adata, model_kind, config)
    )
    embedding = _export_embedding(
        model, adata, evaluation_mask, model_kind, config, output
    )
    result = {
        "schema_version": "masld-cl-reference-v1",
        "model_kind": model_kind,
        "model_class": model_class,
        "seed": seed,
        "n_training_cells": reference.n_obs,
        "n_donors": int(reference.obs["donor_id"].nunique()),
        "genes": list(map(str, reference.var_names)),
        "contract_raw_counts_sha256": lock["raw_counts_sha256"],
        "config_sha256": config["_config_sha256"],
        "embedding": embedding,
        "stopping_epoch": _trained_epochs(model),
        "model_bundle": model_bundle,
    }
    write_json_exclusive(output / "reference_manifest.json", result)
    _write_environment(output)
    return result


def train_reference_sensitivity(
    config: dict[str, Any], prepared_path: str | Path, contract_lock: str | Path,
    prepared_lock: str | Path, policy_value: str | Path, output: str | Path,
    model_kind: str, seed: int,
) -> dict[str, Any]:
    """Fit the frozen-feature local nine-donor healthy-reference sensitivity."""
    import anndata as ad
    from scvi.model import SCANVI, SCVI

    from .scvi_adapter import assert_runtime_versions

    assert_runtime_versions(config)
    set_all_seeds(seed)
    validate_program_firewall(config)
    pipeline_root = Path(config["_config_path"]).resolve().parent
    policy_path = Path(policy_value).resolve()
    if policy_path != pipeline_root / "reference" / "local9_sensitivity_policy_v11.json":
        raise TrainingError("reference sensitivity requires the source-controlled local9 policy")
    with policy_path.open() as handle:
        policy = json.load(handle)
    assessment_path = pipeline_root / policy.get("assessment_relative_path", "")
    if (
        policy.get("schema_version") != "masld-cl-local9-reference-sensitivity-v11"
        or policy.get("config_sha256") != config["_config_sha256"]
        or sha256_path(assessment_path) != policy.get("assessment_sha256")
        or policy.get("sensitivity_only") is not True
    ):
        raise TrainingError("local9 reference policy or assessment differs")
    lock = verify_contract_lock(config, contract_lock, full_hash=False)
    prepared_identity = verify_prepared_file(prepared_path, prepared_lock, config, lock)
    adata = ad.read_h5ad(prepared_path)
    verify_prepared(adata, config, lock, prepared_identity)
    roster = set(map(str, policy["donors"]))
    donor_values = adata.obs["donor_id"].astype(str)
    sensitivity_reference = (
        adata.obs["analysis_eligible"].astype(bool).to_numpy()
        & donor_values.isin(roster).to_numpy()
    )
    observed_donors = set(donor_values[sensitivity_reference])
    if (
        observed_donors != roster
        or len(observed_donors) != int(policy["expected_donors"])
        or int(sensitivity_reference.sum()) != int(policy["expected_cells"])
    ):
        raise TrainingError("local9 reference roster or census differs from policy")
    original_strict = adata.obs["strict_reference"].astype(bool).to_numpy(copy=True)
    if not np.all(sensitivity_reference[original_strict]):
        raise TrainingError("local9 sensitivity dropped a strict seven-donor reference cell")
    adata.obs["strict_reference"] = sensitivity_reference
    mask = sensitivity_reference & _model_mask(adata, model_kind, config)
    indices = _capped_indices(adata, np.flatnonzero(mask), model_kind, config, seed)
    reference = adata[indices].copy()
    _drop_unused_categories(reference, (config["features"]["selection_batch_key"],))
    _set_model_labels(reference, query_unknown=False, unlabeled=config["features"]["unlabeled_category"])
    if model_kind == "all_lineage":
        _validate_scanvi_label_codes(reference, unlabeled=config["features"]["unlabeled_category"])
    batch_key = config["features"]["selection_batch_key"]
    labels_key = config["features"]["labels_key"]
    external = _external_donor_split(reference, seed)
    SCVI.setup_anndata(reference, layer="counts", batch_key=batch_key, labels_key=labels_key)
    vae = SCVI(reference, **_architecture(config))
    vae.train(**_train_kwargs(config, model_kind, external))
    if model_kind == "all_lineage":
        model = SCANVI.from_scvi_model(
            vae, unlabeled_category=config["features"]["unlabeled_category"]
        )
        model.train(**_train_kwargs(config, model_kind, external))
        model_class = "SCANVI"
    else:
        model, model_class = vae, "SCVI"
    output = _ensure_output(output)
    model.save(output / "model", save_anndata=True)
    model_bundle = _model_bundle_identity(output / "model")
    np.save(output / "reference_latent.npy", model.get_latent_representation(reference))
    evaluation_mask = sensitivity_reference & _model_mask(adata, model_kind, config)
    embedding = _export_embedding(model, adata, evaluation_mask, model_kind, config, output)
    result = {
        "schema_version": "masld-cl-reference-v1", "model_kind": model_kind,
        "model_class": model_class, "seed": seed, "sensitivity_only": True,
        "reference_roster": policy["roster_name"], "reference_donors": sorted(roster),
        "reference_policy": str(policy_path), "reference_policy_sha256": sha256_path(policy_path),
        "n_reference_cells_uncapped": int(sensitivity_reference.sum()),
        "n_training_cells": int(reference.n_obs), "n_donors": int(reference.obs["donor_id"].nunique()),
        "genes": list(map(str, reference.var_names)),
        "contract_raw_counts_sha256": lock["raw_counts_sha256"],
        "config_sha256": config["_config_sha256"], "embedding": embedding,
        "stopping_epoch": _trained_epochs(model), "model_bundle": model_bundle,
    }
    write_json_exclusive(output / "reference_manifest.json", result)
    _write_environment(output)
    return result


def train_de_novo(
    config: dict[str, Any], prepared_path: str | Path, contract_lock: str | Path,
    prepared_lock: str | Path, output: str | Path, model_kind: str, seed: int, production: bool,
    selection_lock: str | Path | None = None,
    query_datasets: list[str] | None = None,
) -> dict[str, Any]:
    """Fit the prespecified de novo joint scVI/scANVI baseline."""
    import anndata as ad
    from scvi.model import SCANVI, SCVI

    from .scvi_adapter import assert_runtime_versions
    from .firewall import load_outcome_selection_lock

    assert_runtime_versions(config)
    set_all_seeds(seed)
    validate_program_firewall(config)
    lock = verify_contract_lock(config, contract_lock, full_hash=False)
    prepared_identity = verify_prepared_file(
        prepared_path, prepared_lock, config, lock
    )
    frozen_selection = None
    if selection_lock is not None:
        frozen_selection = load_outcome_selection_lock(selection_lock, config)
    if production and frozen_selection is None:
        raise TrainingError("production baseline fitting requires the immutable selection lock")
    full = ad.read_h5ad(prepared_path)
    verify_prepared(full, config, lock, prepared_identity)
    reference_indices, query_indices = _reference_and_query_indices(
        full, model_kind, config, production, query_datasets
    )
    reference_pool = _capped_indices(full, reference_indices, model_kind, config, seed + 1)
    query_pool = _capped_indices(full, query_indices, model_kind, config, seed + 2)
    indices = concatenate_query_and_replay(query_pool, reference_pool)
    training = full[indices].copy()
    _drop_unused_categories(training, (config["features"]["selection_batch_key"],))
    _set_model_labels(
        training,
        query_unknown=model_kind == "all_lineage",
        unlabeled=config["features"]["unlabeled_category"],
    )
    if model_kind == "all_lineage":
        _validate_scanvi_label_codes(
            training, unlabeled=config["features"]["unlabeled_category"]
        )
    external = _external_donor_split(training, seed)
    batch_key = config["features"]["selection_batch_key"]
    labels_key = config["features"]["labels_key"]
    SCVI.setup_anndata(training, layer="counts", batch_key=batch_key, labels_key=labels_key)
    vae = SCVI(training, **_architecture(config))
    vae.train(**_train_kwargs(config, model_kind, external))
    if model_kind == "all_lineage":
        model = SCANVI.from_scvi_model(
            vae, unlabeled_category=config["features"]["unlabeled_category"]
        )
        model.train(**_train_kwargs(config, model_kind, external))
        model_class = "SCANVI"
    else:
        model = vae
        model_class = "SCVI"
    output = _ensure_output(output)
    model.save(output / "model", save_anndata=True)
    model_bundle = _model_bundle_identity(output / "model")
    np.save(output / "training_latent.npy", model.get_latent_representation(training))
    evaluation_mask = _evaluation_mask(
        full, model_kind, config, production, query_datasets
    )
    embedding = _export_embedding(
        model, full, evaluation_mask, model_kind, config, output
    )
    result = {
        "schema_version": "masld-cl-denovo-v1",
        "model_kind": model_kind,
        "model_class": model_class,
        "production": production,
        "seed": seed,
        "n_reference_cells": int(len(reference_pool)),
        "n_query_cells": int(len(query_pool)),
        "n_training_cells": int(len(indices)),
        "query_datasets": sorted(set(full.obs.iloc[query_pool]["dataset"].astype(str))),
        "contract_raw_counts_sha256": lock["raw_counts_sha256"],
        "config_sha256": config["_config_sha256"],
        "selection_lock_sha256": None if frozen_selection is None else frozen_selection["lock_sha256"],
        "embedding": embedding,
        "stopping_epoch": _trained_epochs(model),
        "model_bundle": model_bundle,
    }
    write_json_exclusive(output / "denovo_manifest.json", result)
    _write_environment(output)
    return result


def _reference_and_query_indices(
    adata, model_kind: str, config: dict[str, Any], production: bool,
    query_datasets: list[str] | None = None,
):
    lineage = _model_mask(adata, model_kind, config)
    analyzed = adata.obs["analysis_eligible"].astype(bool).to_numpy()
    reference = np.flatnonzero(lineage & analyzed & adata.obs["strict_reference"].astype(bool).to_numpy())
    if production:
        query = np.flatnonzero(lineage & analyzed & ~adata.obs["strict_reference"].astype(bool).to_numpy())
    else:
        query = np.flatnonzero(lineage & analyzed & adata.obs["primary_query"].astype(bool).to_numpy())
    if query_datasets is not None:
        requested = set(query_datasets)
        observed = set(adata.obs.iloc[query]["dataset"].astype(str))
        unknown = sorted(requested - observed)
        if unknown:
            raise TrainingError(f"requested query datasets are absent from eligible query: {unknown}")
        keep = adata.obs.iloc[query]["dataset"].astype(str).isin(requested).to_numpy()
        query = query[keep]
    if not len(reference) or not len(query):
        raise TrainingError("reference and query training pools must both be nonempty")
    return reference, query


def _sample_replay(adata, reference: np.ndarray, fraction: float, seed: int) -> np.ndarray:
    obs = adata.obs.iloc[reference]
    strata = list(zip(
        obs["donor_id"].astype(str), obs["audit_cell_type"].astype(str),
        obs["preparation_method"].astype(str),
    ))
    return balanced_exact_fraction_indices(reference, strata, fraction, seed)


def _filter_fisher_minimum(
    adata, indices: np.ndarray, minimum: int,
) -> np.ndarray:
    obs = adata.obs.iloc[indices]
    groups = list(zip(obs["donor_id"].astype(str), obs["audit_cell_type"].astype(str)))
    counts: dict[tuple[str, str], int] = {}
    for group in groups:
        counts[group] = counts.get(group, 0) + 1
    keep = np.fromiter((counts[group] >= minimum for group in groups), dtype=bool)
    return indices[keep]


def _filter_fisher_from_eligibility_pool(
    adata, sampled: np.ndarray, eligibility_pool: np.ndarray, minimum: int,
) -> np.ndarray:
    """Filter sampled cells using donor-lineage counts in the source pool."""
    pool_obs = adata.obs.iloc[eligibility_pool]
    counts: dict[tuple[str, str], int] = {}
    for group in zip(
        pool_obs["donor_id"].astype(str), pool_obs["audit_cell_type"].astype(str)
    ):
        counts[group] = counts.get(group, 0) + 1
    sample_obs = adata.obs.iloc[sampled]
    groups = list(zip(
        sample_obs["donor_id"].astype(str), sample_obs["audit_cell_type"].astype(str)
    ))
    keep = np.fromiter((counts.get(group, 0) >= minimum for group in groups), dtype=bool)
    return sampled[keep]


def _sample_fisher(
    adata, candidates: np.ndarray, fraction: float, seed: int,
    config: dict[str, Any],
) -> np.ndarray:
    eligible = _filter_fisher_minimum(
        adata, candidates, config["sampling"]["lineage_fisher_min_cells"]
    )
    obs = adata.obs.iloc[eligible]
    strata = list(zip(
        obs["donor_id"].astype(str), obs["audit_cell_type"].astype(str),
        obs["preparation_method"].astype(str),
    ))
    return balanced_exact_fraction_indices(eligible, strata, fraction, seed)


def _sample_controls(
    adata, query: np.ndarray, fraction: float, seed: int, config: dict[str, Any],
) -> np.ndarray:
    eligible = query[
        adata.obs.iloc[query]["primary_query"].astype(bool).to_numpy()
        & adata.obs.iloc[query]["query_control"].astype(bool).to_numpy()
    ]
    return _sample_fisher(adata, eligible, fraction, seed, config)


def _reference_fisher_indices(
    adata, method: str, replay: np.ndarray, reference_pool: np.ndarray,
    config: dict[str, Any], seed: int,
) -> np.ndarray:
    if method in {"continual_learning", "continual_distillation"}:
        result = _filter_fisher_from_eligibility_pool(
            adata, replay, reference_pool,
            config["sampling"]["lineage_fisher_min_cells"],
        )
        if not set(map(int, result)).issubset(set(map(int, replay))):
            raise TrainingError("reference Fisher contains cells outside replay")
        return result
    if method == "ewc_only":
        return _sample_fisher(
            adata, reference_pool, config["sampling"]["primary_replay_fraction"],
            seed, config,
        )
    raise TrainingError(f"reference Fisher is undefined for method: {method}")


def _registered_loader(model, adata, batch_size: int, *, drop_last: bool):
    registered = model._validate_anndata(adata)
    return model._make_data_loader(
        adata=registered, batch_size=batch_size, shuffle=False, drop_last=drop_last
    )


def train_update(
    config: dict[str, Any], prepared_path: str | Path, contract_lock: str | Path,
    prepared_lock: str | Path, reference_model: str | Path, output: str | Path,
    model_kind: str, method: str,
    ewc_lambda: float, replay_fraction: float, seed: int, production: bool,
    selection_lock: str | Path | None = None,
    query_datasets: list[str] | None = None,
    control_fisher_datasets: list[str] | None = None,
    held_out_datasets: list[str] | None = None,
    export_embedding: bool = True,
    refit_lock: str | Path | None = None,
    replay_mode: str = "random",
    bi_replay_unlock: str | Path | None = None,
    distillation_weight: float = 0.0,
    reference_embedding_path: str | Path | None = None,
    policy_value: str | Path | None = None,
) -> dict[str, Any]:
    import anndata as ad
    from scvi.model import SCANVI, SCVI

    from .scvi_adapter import (
        assert_runtime_versions, attach_distillation_training_plan,
        attach_ewc_distillation_training_plan, attach_training_plan,
        scvi_reconstruction_loss,
    )
    from .firewall import load_outcome_selection_lock

    if method not in METHODS:
        raise TrainingError(f"unknown update method: {method}")
    if not np.isfinite(distillation_weight) or distillation_weight < 0:
        raise TrainingError("distillation weight must be finite and non-negative")
    if replay_mode != "random" and replay_mode not in BI_MODES:
        raise TrainingError(f"unknown replay mode: {replay_mode}")
    if replay_mode != "random" and method != "continual_learning":
        raise TrainingError("BI replay is defined only for continual learning")
    if replay_mode == "random" and bi_replay_unlock is not None:
        raise TrainingError("random replay cannot consume a BI replay unlock")
    if replay_mode in BI_MODES and bi_replay_unlock is None:
        raise TrainingError("BI replay requires a promotion-derived unlock")
    if method == "architecture_surgery" and (ewc_lambda != 0 or replay_fraction != 0):
        raise TrainingError("architecture surgery must use zero EWC and zero replay")
    if method == "fine_tune" and (ewc_lambda != 0 or replay_fraction != 0):
        raise TrainingError("fine tuning must use zero EWC and zero replay")
    if method == "replay_only" and ewc_lambda != 0:
        raise TrainingError("replay-only must use zero EWC")
    if method == "ewc_only" and replay_fraction != 0:
        raise TrainingError("EWC-only cannot put replay cells in the ELBO")
    if method == "continual_learning" and (ewc_lambda <= 0 or replay_fraction <= 0):
        raise TrainingError("continual learning requires positive EWC and replay")
    if method == "continual_distillation" and (
        ewc_lambda <= 0 or replay_fraction <= 0 or distillation_weight <= 0
    ):
        raise TrainingError(
            "continual distillation requires positive EWC, replay, and distillation"
        )
    if method == "latent_distillation":
        if ewc_lambda != 0 or replay_fraction <= 0 or distillation_weight <= 0:
            raise TrainingError(
                "latent distillation requires zero EWC, positive replay, and positive weight"
            )
    elif method != "continual_distillation" and distillation_weight != 0:
        raise TrainingError("distillation weight is allowed only for latent distillation")
    if method == "continual_distillation" and reference_embedding_path is None:
        raise TrainingError("continual distillation requires the frozen reference embedding")
    if method != "continual_distillation" and reference_embedding_path is not None:
        raise TrainingError("reference-embedding export freezing is V13-only")
    method_policy = None
    method_policy_path = None
    if method == "continual_distillation":
        pipeline_root = Path(config["_config_path"]).resolve().parent
        method_policy_path = Path(policy_value or "").resolve()
        expected_policy = pipeline_root / "reference" / "continual_distillation_policy_v13.json"
        if method_policy_path != expected_policy:
            raise TrainingError("continual distillation requires the source-controlled V13 policy")
        with method_policy_path.open() as handle:
            method_policy = json.load(handle)
        if (
            method_policy.get("schema_version") != "masld-cl-continual-distillation-policy-v13"
            or method_policy.get("config_sha256") != config["_config_sha256"]
            or float(method_policy.get("replay_fraction", -1)) != replay_fraction
            or float(method_policy.get("distillation_weight", -1)) != distillation_weight
            or ewc_lambda not in set(map(float, method_policy.get("ewc_lambda_grid", [])))
            or method_policy.get("case_stage_program_hero_gene_umap_and_cas13_not_used") is not True
        ):
            raise TrainingError("continual-distillation run differs from the prospective V13 policy")
    elif policy_value is not None:
        raise TrainingError("method policy is accepted only for continual distillation")
    frozen_selection = None
    fixed_epochs = None
    if selection_lock is not None:
        frozen_selection = load_outcome_selection_lock(selection_lock, config)
    if production:
        if frozen_selection is None:
            raise TrainingError("production fitting requires the immutable selection lock")
        if method == "continual_learning":
            if frozen_selection.get("schema_version") != "masld-cl-selection-v1":
                raise TrainingError("replay-plus-EWC fitting requires its own frozen hyperparameter selection")
            selected = frozen_selection["selected"]
            if (ewc_lambda, replay_fraction) != (
                float(selected["ewc_lambda"]), float(selected["replay_fraction"])
            ):
                raise TrainingError("production hyperparameters differ from the frozen selection")
            if refit_lock is None:
                raise TrainingError("production continual learning requires the median-epoch refit lock")
            from .refit import load_refit_lock
            frozen_refit = load_refit_lock(refit_lock, config, frozen_selection)
            fixed_epochs = int(
                frozen_refit["median_stopping_epochs"][f"continual_learning|{model_kind}"]
            )
    elif frozen_selection is not None and method == "continual_learning":
        if frozen_selection.get("schema_version") != "masld-cl-selection-v1":
            raise TrainingError("replay-plus-EWC confirmation requires its own frozen selection")
        allowed_pairs = {
            (
                float(frozen_selection["selected"]["ewc_lambda"]),
                float(frozen_selection["selected"]["replay_fraction"]),
            ),
            (100.0, 0.20),
        }
        runner = frozen_selection.get("pareto_runner_up")
        if runner is not None:
            allowed_pairs.add((
                float(runner["ewc_lambda"]), float(runner["replay_fraction"])
            ))
        if method == "continual_learning" and (ewc_lambda, replay_fraction) not in allowed_pairs:
            raise TrainingError("post-selection confirmation setting is not selected, runner-up, or paper")
    if replay_mode in BI_MODES and (not production or frozen_selection is None):
        raise TrainingError("BI replay is sensitivity-only after random-replay promotion")
    bi_unlock = None
    if replay_mode in BI_MODES:
        bi_unlock = load_bi_replay_unlock(
            bi_replay_unlock, config, frozen_selection
        )
    assert_runtime_versions(config)
    set_all_seeds(seed)
    validate_program_firewall(config)
    lock = verify_contract_lock(config, contract_lock, full_hash=False)
    prepared_identity = verify_prepared_file(
        prepared_path, prepared_lock, config, lock
    )
    reference_identity = _verify_reference_model_bundle(
        reference_model, config, model_kind
    )
    full = ad.read_h5ad(prepared_path)
    verify_prepared(full, config, lock, prepared_identity)
    held_out = _validate_held_study_rosters(
        config, query_datasets, control_fisher_datasets, held_out_datasets
    )
    reference_indices, query_indices = _reference_and_query_indices(
        full, model_kind, config, production, query_datasets
    )
    reference_pool = _capped_indices(full, reference_indices, model_kind, config, seed + 1)
    query_pool = _capped_indices(full, query_indices, model_kind, config, seed + 2)
    random_replay = (
        _sample_replay(full, reference_pool, replay_fraction, seed + 3)
        if replay_fraction > 0 else np.asarray([], dtype=np.int64)
    )
    expected_replay = int(
        np.floor(len(reference_pool) * replay_fraction + 0.5)
    )
    if len(random_replay) != expected_replay:
        raise TrainingError("random replay did not contain the exact requested fraction")
    batch_key = config["features"]["selection_batch_key"]
    model_cls = SCANVI if model_kind == "all_lineage" else SCVI
    reference = model_cls.load(reference_model)
    replay = random_replay
    bi_scores = None
    bi_scores_identity = None
    if replay_mode in BI_MODES:
        reference_for_bi = full[reference_pool].copy()
        _drop_unused_categories(reference_for_bi, (batch_key,))
        _set_model_labels(
            reference_for_bi, query_unknown=False,
            unlabeled=config["features"]["unlabeled_category"],
        )
        bi_scores = compute_bi_scores(
            reference.module,
            _registered_loader(reference, reference_for_bi, 8, drop_last=False),
            seed=seed + 6, n_augmentations=200, mask_fraction=0.50,
        )
        if len(bi_scores) != len(reference_pool):
            raise TrainingError("BI scoring did not return one value per reference cell")
        replay_obs = full.obs.iloc[reference_pool]
        strata = list(zip(
            replay_obs["donor_id"].astype(str),
            replay_obs["audit_cell_type"].astype(str),
            replay_obs["preparation_method"].astype(str),
        ))
        random_obs = full.obs.iloc[random_replay]
        random_strata = list(zip(
            random_obs["donor_id"].astype(str),
            random_obs["audit_cell_type"].astype(str),
            random_obs["preparation_method"].astype(str),
        ))
        replay = stratified_exact_bi_replay(
            reference_pool, bi_scores, strata, len(random_replay),
            replay_mode.removeprefix("bi_"),
            target_counts=Counter(random_strata),
        )
        bi_scores_identity = {
            "n_cells": int(len(bi_scores)),
            "sha256_float64_le": hashlib.sha256(
                np.asarray(bi_scores, dtype="<f8").tobytes()
            ).hexdigest(),
            "minimum": float(bi_scores.min()),
            "maximum": float(bi_scores.max()),
            "n_augmentations": 200,
            "mask_fraction": 0.50,
            "batch_size": 8,
            "latent_statistic": "posterior_sample_z_as_upstream",
        }
    training_indices = concatenate_query_and_replay(query_pool, replay)
    training = full[training_indices].copy()
    _drop_unused_categories(training, (batch_key,))
    _set_model_labels(training, query_unknown=True, unlabeled=config["features"]["unlabeled_category"])
    if method in {"latent_distillation", "continual_distillation"}:
        training.obs[DISTILLATION_INDEX_KEY] = np.arange(
            training.n_obs, dtype=np.int64
        )
    reference_state = {name: value.detach().clone() for name, value in reference.module.named_parameters()}
    unfrozen = method != "architecture_surgery"
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    model = model_cls.load_query_data(
        training, reference, unfrozen=unfrozen, accelerator=accelerator, device=1
    )
    trainable_by_name = {
        name: bool(parameter.requires_grad) for name, parameter in model.module.named_parameters()
    }
    _validate_update_trainability(
        trainable_by_name, model_kind=model_kind, method=method
    )
    external = _external_donor_split(training, seed)
    plan_kwargs: dict[str, Any] = {}
    fisher_summaries = None
    ewc_tensor_bundle = None
    distillation_summary = None
    distillation_tensor_bundle = None
    if method in {"ewc_only", "continual_learning", "continual_distillation"}:
        fisher_reference = _reference_fisher_indices(
            full, method, replay, reference_pool, config, seed + 4
        )
        control_candidates = np.flatnonzero(
            _model_mask(full, model_kind, config)
            & full.obs["analysis_eligible"].astype(bool).to_numpy()
            & full.obs["primary_query"].astype(bool).to_numpy()
        )
        powered = set(config["evaluation"]["powered_query_studies"])
        if control_fisher_datasets is not None:
            requested_controls = set(control_fisher_datasets)
            if not requested_controls or not requested_controls.issubset(powered):
                raise TrainingError(
                    "control Fisher datasets must be a nonempty subset of powered query studies"
                )
            control_candidates = control_candidates[
                full.obs.iloc[control_candidates]["dataset"].astype(str)
                .isin(requested_controls).to_numpy()
            ]
        controls = _sample_controls(
            full, control_candidates, config["sampling"]["control_fisher_fraction"], seed + 5,
            config,
        )
        control_eligible = control_candidates[
            full.obs.iloc[control_candidates]["primary_query"].astype(bool).to_numpy()
            & full.obs.iloc[control_candidates]["query_control"].astype(bool).to_numpy()
        ]
        control_eligible = _filter_fisher_minimum(
            full, control_eligible,
            config["sampling"]["lineage_fisher_min_cells"],
        )
        expected_controls = int(np.floor(
            len(control_eligible)
            * config["sampling"]["control_fisher_fraction"] + 0.5
        ))
        if len(controls) != expected_controls:
            raise TrainingError("control Fisher did not contain the exact requested fraction")
        if not len(fisher_reference) or not len(controls):
            raise TrainingError("EWC requires nonempty reference and query-control Fisher samples")
        fisher_reference_adata = full[fisher_reference].copy()
        controls_adata = full[controls].copy()
        _drop_unused_categories(fisher_reference_adata, (batch_key,))
        _drop_unused_categories(controls_adata, (batch_key,))
        _set_model_labels(fisher_reference_adata, query_unknown=False, unlabeled=config["features"]["unlabeled_category"])
        _set_model_labels(controls_adata, query_unknown=True, unlabeled=config["features"]["unlabeled_category"])
        anchors, masks = build_anchor_and_masks(model.module, reference_state)
        batch_size = config["sampling"]["fisher_batch_size"]
        reference_fisher, reference_summary = compute_empirical_fisher(
            model.module,
            _registered_loader(
                model, fisher_reference_adata, batch_size, drop_last=False
            ),
            scvi_reconstruction_loss,
            expected_batch_size=batch_size,
            include_partial=True,
        )
        # The model has not stepped: both Fishers are evaluated at the same anchor.
        control_fisher, control_summary = compute_empirical_fisher(
            model.module,
            _registered_loader(
                model, controls_adata, batch_size, drop_last=False
            ),
            scvi_reconstruction_loss,
            expected_batch_size=batch_size,
            include_partial=True,
        )
        if reference_summary.n_observations != len(fisher_reference):
            raise TrainingError("reference Fisher did not consume every sampled cell")
        if control_summary.n_observations != len(controls):
            raise TrainingError("control Fisher did not consume every sampled cell")
        ewc_normalization = "none"
        ewc_excluded_parameter_prefixes: list[str] = []
        if method == "continual_distillation":
            ewc_normalization = "block_mean_product"
            ewc_excluded_parameter_prefixes = ["classifier."]
            for name, mask in masks.items():
                if any(name.startswith(prefix) for prefix in ewc_excluded_parameter_prefixes):
                    mask.zero_()
        regularizer = EWCRegularizer(
            anchors, masks, reference_fisher, control_fisher,
            normalization=ewc_normalization,
        )
        regularizer.assert_zero_at_anchor(model.module)
        novel_batches = sorted(
            set(full.obs.iloc[query_pool][batch_key].astype(str))
            - set(full.obs.iloc[reference_pool][batch_key].astype(str))
        )
        masked_elements = int(sum((value == 0).sum().item() for value in masks.values()))
        if novel_batches and masked_elements == 0:
            raise TrainingError("query batch categories did not create unpenalized parameter slices")
        attach_training_plan(model, semi_supervised=model_kind == "all_lineage")
        plan_kwargs = {"ewc_regularizer": regularizer, "ewc_lambda": ewc_lambda}
        fisher_summaries = {
            "estimator": "empirical_diagonal_minibatch_mean_gradient",
            "evaluated_at_same_post_load_anchor": True,
            "reference": asdict(reference_summary),
            "query_control": asdict(control_summary),
            "reference_sample": _sample_identity(full, fisher_reference),
            "query_control_sample": _sample_identity(full, controls),
            "query_control_eligibility": _sample_identity(
                full, control_eligible
            ),
            "query_control_requested_fraction": config["sampling"][
                "control_fisher_fraction"
            ],
            "reference_is_exact_replay_subset": bool(
                method not in {"continual_learning", "continual_distillation"}
                or set(map(int, fisher_reference)).issubset(set(map(int, replay)))
            ),
            "control_fisher_datasets": sorted(
                powered if control_fisher_datasets is None
                else set(control_fisher_datasets)
            ),
            "novel_batch_categories": novel_batches,
            "masked_parameter_elements": masked_elements,
            "normalization": ewc_normalization,
            "excluded_parameter_prefixes": ewc_excluded_parameter_prefixes,
        }
        ewc_tensor_bundle = {
            "anchors": {name: value.detach().cpu() for name, value in anchors.items()},
            "masks": {name: value.detach().cpu() for name, value in masks.items()},
            "reference_fisher": {
                name: value.detach().cpu() for name, value in reference_fisher.items()
            },
            "query_control_fisher": {
                name: value.detach().cpu() for name, value in control_fisher.items()
            },
        }
    if method in {"latent_distillation", "continual_distillation"}:
        from scvi import REGISTRY_KEYS
        from scvi.data.fields import NumericalObsField

        manager = model.get_anndata_manager(training, required=True)
        if REGISTRY_KEYS.INDICES_KEY in manager.data_registry:
            raise TrainingError("distillation cell-index registry already exists")
        manager.register_new_fields([
            NumericalObsField(REGISTRY_KEYS.INDICES_KEY, DISTILLATION_INDEX_KEY)
        ])
        anchor_latent = posterior_means_at_anchor(
            model.module,
            _registered_loader(
                model, training, config["architecture"]["batch_size"],
                drop_last=False,
            ),
        ).numpy().astype(np.float32, copy=False)
        replay_mask = np.isin(training_indices, replay, assume_unique=True)
        if int(replay_mask.sum()) != len(replay):
            raise TrainingError("distillation targets do not match the exact replay roster")
        if np.any(full.obs.iloc[training_indices[replay_mask]]["strict_reference"] != True):
            raise TrainingError("a distillation target is not a strict-reference replay cell")
        targets = np.zeros_like(anchor_latent, dtype=np.float32)
        targets[replay_mask] = anchor_latent[replay_mask]
        regularizer = LatentDistillationRegularizer(
            torch.from_numpy(targets), torch.from_numpy(replay_mask)
        ).to(next(model.module.parameters()).device)
        zero_at_anchor = regularizer.assert_zero_at_anchor(
            model.module,
            _registered_loader(
                model, training, config["architecture"]["batch_size"],
                drop_last=False,
            ),
        )
        if method == "continual_distillation":
            attach_ewc_distillation_training_plan(
                model, semi_supervised=model_kind == "all_lineage"
            )
            plan_kwargs.update({
                "distillation_regularizer": regularizer,
                "distillation_weight": distillation_weight,
            })
        else:
            attach_distillation_training_plan(
                model, semi_supervised=model_kind == "all_lineage"
            )
            plan_kwargs = {
                "distillation_regularizer": regularizer,
                "distillation_weight": distillation_weight,
            }
        replay_targets = np.asarray(targets[replay_mask], dtype="<f4")
        distillation_summary = {
            "weight": float(distillation_weight),
            "target_statistic": "qz.loc",
            "encoder_mode": "evaluation",
            "anchor_after_reference_load_and_batch_expansion": True,
            "target_cells": "exact_replay_only",
            "n_targets": int(replay_mask.sum()),
            "n_query_targets": 0,
            "n_latent": int(targets.shape[1]),
            "zero_at_anchor_maximum": float(zero_at_anchor),
            "target_sha256_float32_le": hashlib.sha256(
                replay_targets.tobytes()
            ).hexdigest(),
            "replay_mask_sha256": hashlib.sha256(
                np.asarray(replay_mask, dtype=np.uint8).tobytes()
            ).hexdigest(),
            "index_obs_key": DISTILLATION_INDEX_KEY,
        }
        distillation_tensor_bundle = {
            "targets": torch.from_numpy(targets),
            "replay_mask": torch.from_numpy(replay_mask),
            "training_cell_ids": list(map(str, training.obs_names)),
        }
    kwargs = _train_kwargs(config, model_kind, external, fixed_epochs)
    if plan_kwargs:
        kwargs["plan_kwargs"] = plan_kwargs
    model.train(**kwargs)
    output = _ensure_output(output)
    if bi_scores is not None:
        bi_scores_path = output / "bi_reference_scores.npy"
        np.save(bi_scores_path, np.asarray(bi_scores, dtype="<f8"))
        bi_scores_identity.update({
            "file": bi_scores_path.name,
            "file_sha256": sha256_path(bi_scores_path),
        })
    model.save(output / "model", save_anndata=True)
    model_bundle = _model_bundle_identity(output / "model")
    if ewc_tensor_bundle is not None:
        ewc_path = output / "ewc_anchor_mask_fishers.pt"
        torch.save(ewc_tensor_bundle, ewc_path)
        fisher_summaries["tensor_bundle"] = {
            "file": ewc_path.name,
            "sha256": sha256_path(ewc_path),
            "parameter_names": sorted(ewc_tensor_bundle["anchors"]),
        }
    if distillation_tensor_bundle is not None:
        distillation_path = output / "latent_distillation_targets.pt"
        torch.save(distillation_tensor_bundle, distillation_path)
        distillation_summary["tensor_bundle"] = {
            "file": distillation_path.name,
            "sha256": sha256_path(distillation_path),
        }
    np.save(output / "training_latent.npy", model.get_latent_representation(training))
    embedding = None
    if export_embedding:
        evaluation_mask = _evaluation_mask(
            full, model_kind, config, production, query_datasets
        )
        embedding = _export_embedding(
            model, full, evaluation_mask, model_kind, config, output,
            frozen_reference_embedding=(
                reference_embedding_path if method == "continual_distillation" else None
            ),
        )
    result = {
        "schema_version": "masld-cl-update-v1",
        "model_kind": model_kind,
        "method": method,
        "production": production,
        "seed": seed,
        "ewc_lambda": ewc_lambda,
        "distillation_weight": float(distillation_weight),
        "replay_fraction": replay_fraction,
        "replay_mode": replay_mode,
        "sensitivity_only": replay_mode in BI_MODES,
        "bi_replay_unlock_sha256": (
            None if bi_unlock is None else bi_unlock["lock_sha256"]
        ),
        "bi_scores": bi_scores_identity,
        "n_query_cells": int(len(query_pool)),
        "n_random_replay_cells": int(len(random_replay)),
        "n_replay_cells": int(len(replay)),
        "n_training_cells": int(len(training_indices)),
        "n_training_unique_cells": int(len(np.unique(training_indices))),
        "n_descriptive_only_training_cells": int(
            (~full.obs.iloc[training_indices]["analysis_eligible"].astype(bool)).sum()
        ),
        "samples": {
            "reference_pool": _sample_identity(full, reference_pool),
            "query_pool": _sample_identity(full, query_pool),
            "random_replay": _sample_identity(full, random_replay),
            "replay": _sample_identity(full, replay),
        },
        "query_datasets": sorted(set(full.obs.iloc[query_pool]["dataset"].astype(str))),
        "held_out_datasets": sorted(held_out),
        "fisher": fisher_summaries,
        "distillation": distillation_summary,
        "trainable_parameters": int(sum(trainable_by_name.values())),
        "total_parameters": int(len(trainable_by_name)),
        "contract_raw_counts_sha256": lock["raw_counts_sha256"],
        "config_sha256": config["_config_sha256"],
        "selection_lock_sha256": None if frozen_selection is None else frozen_selection["lock_sha256"],
        "reference_identity": reference_identity,
        "frozen_reference_embedding": (
            None if reference_embedding_path is None else {
                "path": str(Path(reference_embedding_path).resolve()),
                "sha256": sha256_path(reference_embedding_path),
            }
        ),
        "method_policy": (
            None if method_policy_path is None else {
                "path": str(method_policy_path),
                "sha256": sha256_path(method_policy_path),
            }
        ),
        "model_bundle": model_bundle,
        "embedding": embedding,
        "stopping_epoch": _trained_epochs(model),
        "fixed_refit_epoch": fixed_epochs,
    }
    if result["n_training_cells"] != result["n_training_unique_cells"]:
        raise TrainingError("replay entered the training object more than once")
    if result["n_descriptive_only_training_cells"] != 0:
        raise TrainingError("descriptive-only cells entered optimization")
    write_json_exclusive(output / "update_manifest.json", result)
    _write_environment(output)
    return result


def train_sequence(
    config: dict[str, Any], prepared_path: str | Path, contract_lock: str | Path,
    prepared_lock: str | Path, reference_model: str | Path, output: str | Path,
    model_kind: str, order_name: str, ewc_lambda: float, replay_fraction: float,
    seed: int, selection_lock: str | Path,
    refit_lock: str | Path,
) -> dict[str, Any]:
    if order_name not in config["acquisition_orders"]:
        raise TrainingError(f"unknown acquisition order: {order_name}")
    order = list(config["acquisition_orders"][order_name])
    if len(order) != len(set(order)):
        raise TrainingError(f"acquisition order contains duplicate datasets: {order_name}")
    expected = set(config["acquisition_orders"]["accession"])
    if set(order) != expected:
        raise TrainingError(f"acquisition order does not contain the frozen query studies: {order_name}")
    output = _ensure_output(output)
    current_model = Path(reference_model)
    stages = []
    for index, dataset in enumerate(order):
        stage = output / f"stage_{index + 1:02d}_{dataset}"
        result = train_update(
            config, prepared_path, contract_lock, prepared_lock, current_model,
            stage, model_kind, "continual_learning", ewc_lambda, replay_fraction,
            seed + index, True,
            selection_lock=selection_lock,
            query_datasets=[dataset],
            control_fisher_datasets=None,
            held_out_datasets=None,
            export_embedding=index == len(order) - 1,
            refit_lock=refit_lock,
        )
        stages.append({
            "stage": index + 1, "dataset": dataset,
            "output": str(stage.resolve()), "manifest": result,
        })
        current_model = stage / "model"
    result = {
        "schema_version": "masld-cl-sequence-v1",
        "order_name": order_name,
        "order": order,
        "model_kind": model_kind,
        "seed": seed,
        "ewc_lambda": ewc_lambda,
        "replay_fraction": replay_fraction,
        "stages": stages,
        "final_model": str(current_model.resolve()),
        "config_sha256": config["_config_sha256"],
    }
    write_json_exclusive(output / "sequence_manifest.json", result)
    return result


def _evaluation_mask(
    full, model_kind: str, config: dict[str, Any], production: bool,
    query_datasets: list[str] | None = None,
) -> np.ndarray:
    lineage = _model_mask(full, model_kind, config)
    if production:
        return lineage
    requested = np.zeros(full.n_obs, dtype=bool)
    if query_datasets is not None:
        requested = full.obs["dataset"].astype(str).isin(query_datasets).to_numpy()
    return lineage & (
        full.obs["strict_reference"].astype(bool).to_numpy()
        | full.obs["primary_query"].astype(bool).to_numpy()
        | requested
    )


def _export_embedding(
    model, full, mask: np.ndarray, model_kind: str,
    config: dict[str, Any], output: Path,
    frozen_reference_embedding: str | Path | None = None,
) -> dict[str, Any]:
    """Infer the complete evaluation roster; optimization caps never alter it."""
    import csv

    mapped = full[mask].copy()
    from scvi import REGISTRY_KEYS
    if REGISTRY_KEYS.INDICES_KEY in model.adata_manager.data_registry:
        mapped.obs[DISTILLATION_INDEX_KEY] = np.arange(mapped.n_obs, dtype=np.int64)
    _set_model_labels(
        mapped, query_unknown=True, unlabeled=config["features"]["unlabeled_category"]
    )
    latent = model.get_latent_representation(mapped, batch_size=config["architecture"]["batch_size"])
    frozen_predictions: dict[str, str] = {}
    if frozen_reference_embedding is not None:
        from .embedding import load_embedding, matched_rows
        frozen = load_embedding(frozen_reference_embedding)
        reference_positions = np.flatnonzero(mapped.obs["strict_reference"].astype(bool).to_numpy())
        mapped_reference = mapped.obs.iloc[reference_positions].copy()
        mapped_reference["cell_id"] = mapped_reference.index.astype(str)
        left, right = matched_rows(frozen[2], mapped_reference.reset_index(drop=True))
        if len(left) != len(frozen[2]) or len(right) != len(reference_positions):
            raise TrainingError("frozen reference export has a different cell roster")
        latent[reference_positions[right]] = np.asarray(frozen[1])[left]
        frozen_predictions = dict(zip(
            frozen[2].iloc[left]["cell_id"].astype(str),
            frozen[2].iloc[left]["predicted_cell_type"].astype(str),
        ))
    latent_path = output / "embedding_latent.npy"
    cells_path = output / "embedding_cells.tsv.gz"
    np.save(latent_path, latent)
    with DeterministicGzipTextWriter(cells_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "row_index", "cell_id", "library_id", "donor_id", "dataset",
            "audit_cell_type", "predicted_cell_type", "preparation", "technical_batch",
            "strict_reference", "primary_query", "query_control", "analysis_eligible",
        ])
        if model_kind == "all_lineage" and hasattr(model, "predict"):
            predictions = list(map(str, model.predict(mapped)))
        else:
            predictions = mapped.obs["audit_cell_type"].astype(str).tolist()
        for index, ((cell_id, row), prediction) in enumerate(zip(mapped.obs.iterrows(), predictions)):
            prediction = frozen_predictions.get(str(cell_id), prediction)
            writer.writerow([
                index, cell_id, row["library_id"], row["donor_id"], row["dataset"],
                row["audit_cell_type"], prediction, row["preparation_method"], row["technical_batch"],
                bool(row["strict_reference"]), bool(row["primary_query"]),
                bool(row["query_control"]), bool(row["analysis_eligible"]),
            ])
    bundle = {
        "schema_version": "masld-cl-embedding-v1",
        "model_kind": model_kind,
        "n_cells": int(mapped.n_obs),
        "n_latent": int(latent.shape[1]),
        "n_analysis_ineligible": int((~mapped.obs["analysis_eligible"].astype(bool)).sum()),
        "latent_file": latent_path.name,
        "cells_file": cells_path.name,
        "latent_sha256": sha256_path(latent_path),
        "cells_sha256": sha256_path(cells_path),
        "config_sha256": config["_config_sha256"],
        "reference_coordinates_and_predictions_frozen": frozen_reference_embedding is not None,
    }
    write_json_exclusive(output / "embedding_manifest.json", bundle)
    return bundle


def _write_environment(output: Path) -> None:
    import platform
    record = {
        "python": platform.python_version(),
        "torch": torch.__version__,
        "cuda_available": torch.cuda.is_available(),
        "cuda": torch.version.cuda,
        "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu",
    }
    try:
        import anndata, scanpy, scvi
        record.update({"anndata": anndata.__version__, "scanpy": scanpy.__version__, "scvi_tools": scvi.__version__})
    except ImportError:
        pass
    write_json_exclusive(output / "environment.json", record)
    with (output / "pip_freeze.txt").open("x") as handle:
        subprocess.run([os.sys.executable, "-m", "pip", "freeze"], check=True, stdout=handle, text=True)
