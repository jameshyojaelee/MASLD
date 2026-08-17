"""Replay-cell latent distillation for continual scVI/scANVI updates."""

from __future__ import annotations

from collections.abc import Mapping

import torch
from scvi import REGISTRY_KEYS
from scvi.module._constants import MODULE_KEYS


class LatentDistillationError(RuntimeError):
    """Raised when a latent-distillation invariant is violated."""


@torch.inference_mode()
def posterior_means_at_anchor(module: torch.nn.Module, batches) -> torch.Tensor:
    """Read qz.loc directly before scvi marks a query-loaded model as trained."""
    device = next(module.parameters()).device
    values = []
    module_was_training = module.training
    try:
        module.eval()
        for batch in batches:
            tensors = _full_batch(batch)
            inference_inputs = {
                key: value.to(device) if isinstance(value, torch.Tensor) else value
                for key, value in module._get_inference_input(tensors).items()
            }
            values.append(
                module.inference(**inference_inputs)[MODULE_KEYS.QZ_KEY].loc.detach().cpu()
            )
    finally:
        module.train(module_was_training)
    if not values:
        raise LatentDistillationError("anchor loader yielded no cells")
    result = torch.cat(values, dim=0)
    if result.ndim != 2 or not torch.isfinite(result).all():
        raise LatentDistillationError("anchor posterior means are invalid")
    return result


def _full_batch(batch):
    """Return the full-data tensors from ordinary or semi-supervised batches."""
    if isinstance(batch, Mapping):
        if REGISTRY_KEYS.X_KEY in batch:
            return batch
        if 0 in batch:
            return batch[0]
    if isinstance(batch, (tuple, list)) and batch:
        return batch[0]
    raise LatentDistillationError("cannot identify the full dataset in training batch")


class LatentDistillationRegularizer(torch.nn.Module):
    """Penalize replay-cell posterior-mean drift from the loaded reference."""

    def __init__(self, targets: torch.Tensor, replay_mask: torch.Tensor):
        super().__init__()
        targets = torch.as_tensor(targets).detach().clone()
        replay_mask = torch.as_tensor(replay_mask, dtype=torch.bool).detach().clone()
        if targets.ndim != 2 or targets.shape[0] == 0 or targets.shape[1] == 0:
            raise LatentDistillationError("targets must be a nonempty cell-by-latent matrix")
        if replay_mask.ndim != 1 or replay_mask.shape[0] != targets.shape[0]:
            raise LatentDistillationError("replay mask and target rows must match exactly")
        if not torch.isfinite(targets).all():
            raise LatentDistillationError("latent targets contain non-finite values")
        if not replay_mask.any():
            raise LatentDistillationError("latent distillation requires replay targets")
        self.register_buffer("targets", targets)
        self.register_buffer("replay_mask", replay_mask)

    def penalty(self, module: torch.nn.Module, batch) -> torch.Tensor:
        tensors = _full_batch(batch)
        if REGISTRY_KEYS.INDICES_KEY not in tensors:
            raise LatentDistillationError("training batch lacks exact cell indices")
        device = next(module.parameters()).device
        indices = tensors[REGISTRY_KEYS.INDICES_KEY].to(device=device).long().ravel()
        if indices.numel() == 0:
            raise LatentDistillationError("empty training batch")
        if torch.unique(indices).numel() != indices.numel():
            raise LatentDistillationError("training batch contains duplicate cell indices")
        if int(indices.min()) < 0 or int(indices.max()) >= self.targets.shape[0]:
            raise LatentDistillationError("training cell index lies outside target roster")
        selected = self.replay_mask[indices]
        if not selected.any():
            return next(module.parameters()).sum() * 0.0

        module_was_training = module.training
        try:
            module.eval()
            inference_inputs = {
                key: value.to(device) if isinstance(value, torch.Tensor) else value
                for key, value in module._get_inference_input(tensors).items()
            }
            outputs = module.inference(**inference_inputs)
            current = outputs[MODULE_KEYS.QZ_KEY].loc
        finally:
            module.train(module_was_training)
        if current.ndim != 2 or current.shape[0] != indices.shape[0]:
            raise LatentDistillationError("posterior mean does not match batch rows")
        target = self.targets[indices]
        if current.shape[1] != target.shape[1]:
            raise LatentDistillationError("posterior and target latent dimensions differ")
        value = (current[selected] - target[selected]).square().sum(dim=1).mean()
        if not torch.isfinite(value):
            raise LatentDistillationError("latent-distillation penalty is non-finite")
        return value

    def assert_zero_at_anchor(
        self, module: torch.nn.Module, batches, *, atol: float = 1e-7,
    ) -> float:
        maximum = 0.0
        seen = torch.zeros(self.replay_mask.shape, dtype=torch.bool, device="cpu")
        for batch in batches:
            tensors = _full_batch(batch)
            indices = tensors[REGISTRY_KEYS.INDICES_KEY].long().ravel().cpu()
            selected = self.replay_mask.cpu()[indices]
            seen[indices[selected]] = True
            value = float(self.penalty(module, batch).detach().cpu())
            maximum = max(maximum, value)
        missing = self.replay_mask.cpu() & ~seen
        if missing.any():
            raise LatentDistillationError("anchor check did not visit every replay target")
        if maximum > atol:
            raise LatentDistillationError(
                f"latent-distillation penalty is not zero at anchor: {maximum}"
            )
        return maximum
