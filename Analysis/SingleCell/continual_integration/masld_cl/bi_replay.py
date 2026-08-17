"""Promotion-gated Bregman-information replay sensitivity."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import torch

from .config import canonical_json_bytes, write_json_exclusive
from .contracts import ContractError, sha256_path
from .firewall import load_selection_lock
from .metric_bundle import verify_metric_bundle
from .sampling import bregman_information_lse, masked_augmentations


BI_MODES = {"bi_bottom", "bi_top", "bi_step"}


def compute_bi_scores(
    module: torch.nn.Module, batches: Iterable[dict[str, torch.Tensor]],
    *, seed: int, n_augmentations: int = 200, mask_fraction: float = 0.50,
) -> np.ndarray:
    """Compute BI from exactly 200 half-gene-masked encoder passes per cell."""
    from scvi import REGISTRY_KEYS
    from scvi.module._constants import MODULE_KEYS

    device = next(module.parameters()).device
    generator = torch.Generator(device=device)
    generator.manual_seed(seed)
    prior_mode = module.training
    module.eval()
    values = []
    cuda_devices = (
        [device.index if device.index is not None else torch.cuda.current_device()]
        if device.type == "cuda" else []
    )
    try:
        with torch.random.fork_rng(devices=cuda_devices), torch.inference_mode():
            torch.manual_seed(seed)
            if device.type == "cuda":
                torch.cuda.manual_seed_all(seed)
            for batch in batches:
                tensors = {
                    key: value.to(device) if isinstance(value, torch.Tensor) else value
                    for key, value in batch.items()
                }
                expression = tensors[REGISTRY_KEYS.X_KEY]
                augmented_expression = masked_augmentations(
                    expression, n_augmentations=n_augmentations,
                    mask_fraction=mask_fraction, generator=generator,
                )
                augmented = {}
                for key, value in tensors.items():
                    if not isinstance(value, torch.Tensor):
                        augmented[key] = value
                    elif key == REGISTRY_KEYS.X_KEY:
                        augmented[key] = augmented_expression.reshape(
                            n_augmentations * expression.shape[0], expression.shape[1]
                        )
                    elif value.ndim == 0:
                        augmented[key] = value
                    else:
                        augmented[key] = value.repeat(
                            n_augmentations, *([1] * (value.ndim - 1))
                        )
                outputs = module.inference(**module._get_inference_input(augmented))
                encoded = outputs[MODULE_KEYS.Z_KEY]
                if encoded is None:
                    raise ContractError("BI replay requires an encoded latent sample")
                latent = encoded.reshape(
                    n_augmentations, expression.shape[0], -1
                )
                score = bregman_information_lse(latent).clamp_min(0)
                values.append(score.detach().cpu().numpy())
    finally:
        module.train(prior_mode)
    if not values:
        raise ContractError("BI replay received no reference batches")
    result = np.concatenate(values).astype(np.float64, copy=False)
    if not np.isfinite(result).all() or np.any(result < 0):
        raise ContractError("BI replay scores are non-finite or negative")
    return result


def write_bi_replay_unlock(
    config: dict[str, Any], selection_lock: str | Path,
    promotion_decision: str | Path, metrics: str | Path,
    metrics_lock: str | Path, output: str | Path,
) -> dict[str, Any]:
    selection = load_selection_lock(selection_lock, config)
    verified = verify_metric_bundle(config, metrics, metrics_lock, selection)
    with Path(promotion_decision).open() as handle:
        decision = json.load(handle)
    if (
        decision.get("schema_version") != "masld-cl-promotion-v1"
        or decision.get("config_sha256") != config["_config_sha256"]
        or decision.get("selection_lock_sha256") != selection["lock_sha256"]
        or decision.get("metric_bundle_lock_sha256") != verified["lock_sha256"]
        or decision.get("decision") != "promote_main"
        or decision.get("all_gates_pass") is not True
        or decision.get("claim_allowed") is not True
    ):
        raise ContractError("random-replay model did not pass the locked promotion decision")
    unlock = {
        "schema_version": "masld-cl-bi-replay-unlock-v1",
        "config_sha256": config["_config_sha256"],
        "selection_lock_sha256": selection["lock_sha256"],
        "selected_setting": selection["selected"],
        "promotion_decision_realpath": str(Path(promotion_decision).resolve()),
        "promotion_decision_sha256": sha256_path(promotion_decision),
        "metrics_realpath": str(Path(metrics).resolve()),
        "metrics_sha256": sha256_path(metrics),
        "metrics_lock_realpath": str(Path(metrics_lock).resolve()),
        "metrics_lock_sha256": sha256_path(metrics_lock),
        "metric_bundle_lock_sha256": verified["lock_sha256"],
        "sensitivity_only": True,
    }
    unlock["lock_sha256"] = hashlib.sha256(canonical_json_bytes(unlock)).hexdigest()
    write_json_exclusive(output, unlock)
    return unlock


def load_bi_replay_unlock(
    path: str | Path, config: dict[str, Any], selection: dict[str, Any],
) -> dict[str, Any]:
    with Path(path).open() as handle:
        unlock = json.load(handle)
    payload = {key: value for key, value in unlock.items() if key != "lock_sha256"}
    if (
        unlock.get("schema_version") != "masld-cl-bi-replay-unlock-v1"
        or unlock.get("config_sha256") != config["_config_sha256"]
        or unlock.get("selection_lock_sha256") != selection["lock_sha256"]
        or unlock.get("selected_setting") != selection["selected"]
        or unlock.get("sensitivity_only") is not True
        or hashlib.sha256(canonical_json_bytes(payload)).hexdigest()
        != unlock.get("lock_sha256")
    ):
        raise ContractError("BI replay unlock is invalid")
    for path_key, hash_key in (
        ("promotion_decision_realpath", "promotion_decision_sha256"),
        ("metrics_realpath", "metrics_sha256"),
        ("metrics_lock_realpath", "metrics_lock_sha256"),
    ):
        if sha256_path(unlock[path_key]) != unlock[hash_key]:
            raise ContractError("BI replay unlock source changed")
    return unlock
