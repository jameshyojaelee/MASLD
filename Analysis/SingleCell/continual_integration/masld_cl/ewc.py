"""Corrected empirical-Fisher EWC implementation.

The implementation is deliberately independent of scvi-tools. The narrow scvi
adapter lives in :mod:`masld_cl.scvi_adapter`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Iterable, Mapping

import torch
from torch import nn


class EWCError(RuntimeError):
    """Raised when an EWC invariant is violated."""


TensorMap = Mapping[str, torch.Tensor]


def _shared_slices(reference: torch.Tensor, current: torch.Tensor) -> tuple[slice, ...]:
    if reference.ndim != current.ndim:
        raise EWCError(f"rank changed from {reference.ndim} to {current.ndim}")
    if any(old > new for old, new in zip(reference.shape, current.shape)):
        raise EWCError(f"current shape {tuple(current.shape)} cannot contain {tuple(reference.shape)}")
    return tuple(slice(0, int(size)) for size in reference.shape)


@dataclass(frozen=True)
class FisherSummary:
    n_batches: int
    n_observations: int
    expected_batch_size: int | None
    skipped_partial_batches: int
    weighting: str
    squared_norm_by_parameter: dict[str, float]


def build_anchor_and_masks(
    model: nn.Module, reference_state: TensorMap, *, require_bitwise_loaded: bool = True
) -> tuple[dict[str, torch.Tensor], dict[str, torch.Tensor]]:
    """Snapshot the post-load model and mark only reference-owned slices.

    `reference_state` is the state dict saved before query category expansion.
    The model must already contain the loaded reference values. New category
    columns are present in the anchor for serialization but receive a zero mask.
    """
    current = dict(model.named_parameters())
    anchors: dict[str, torch.Tensor] = {}
    masks: dict[str, torch.Tensor] = {}
    missing = sorted(set(reference_state) - set(current))
    if missing:
        raise EWCError(f"reference parameters missing from query model: {missing}")
    for name, parameter in current.items():
        anchor = parameter.detach().clone()
        mask = torch.zeros_like(anchor)
        if name in reference_state:
            old = reference_state[name].detach().to(device=parameter.device, dtype=parameter.dtype)
            shared = _shared_slices(old, parameter)
            if require_bitwise_loaded and not torch.equal(parameter.detach()[shared], old):
                raise EWCError(f"shared slice is not the loaded reference state: {name}")
            mask[shared] = 1
        anchors[name] = anchor
        masks[name] = mask
    return anchors, masks


def compute_empirical_fisher(
    model: nn.Module,
    batches: Iterable[Any],
    loss_fn: Callable[[nn.Module, Any], torch.Tensor],
    *,
    expected_batch_size: int | None = None,
    include_partial: bool = False,
    max_batches: int | None = None,
) -> tuple[dict[str, torch.Tensor], FisherSummary]:
    """Estimate observation-weighted squared minibatch-mean gradients.

    Partial batches are excluded when `expected_batch_size` is supplied unless
    `include_partial` is true. Weighting by batch size preserves the ordinary
    full-minibatch estimate while allowing every presampled cell to contribute.
    This remains an empirical diagonal minibatch-gradient approximation, not an
    exact per-cell Fisher.
    """
    parameters = dict(model.named_parameters())
    trainable = {name: p for name, p in parameters.items() if p.requires_grad}
    if not trainable:
        raise EWCError("model has no trainable parameters")
    fisher = {name: torch.zeros_like(p, memory_format=torch.preserve_format) for name, p in parameters.items()}
    n_batches = 0
    n_observations = 0
    total_weight = 0
    skipped = 0
    prior_mode = model.training
    model.eval()
    try:
        for batch in batches:
            batch_size = infer_batch_size(batch)
            if (
                expected_batch_size is not None
                and batch_size != expected_batch_size
                and not include_partial
            ):
                skipped += 1
                continue
            if max_batches is not None and n_batches >= max_batches:
                break
            model.zero_grad(set_to_none=True)
            loss = loss_fn(model, batch)
            if loss.ndim != 0 or not torch.isfinite(loss):
                raise EWCError("Fisher loss must be a finite scalar")
            loss.backward()
            for name, parameter in trainable.items():
                if parameter.grad is not None:
                    fisher[name].add_(
                        parameter.grad.detach().square(), alpha=batch_size
                    )
            n_batches += 1
            n_observations += batch_size
            total_weight += batch_size
    finally:
        model.train(prior_mode)
        model.zero_grad(set_to_none=True)
    if n_batches == 0:
        raise EWCError("no eligible batches available for Fisher estimation")
    for value in fisher.values():
        value.div_(total_weight)
        if not torch.isfinite(value).all() or torch.any(value < 0):
            raise EWCError("Fisher tensor is non-finite or negative")
    summary = FisherSummary(
        n_batches=n_batches,
        n_observations=n_observations,
        expected_batch_size=expected_batch_size,
        skipped_partial_batches=skipped,
        weighting="observations_over_squared_minibatch_mean_gradients",
        squared_norm_by_parameter={name: float(value.square().sum().cpu()) for name, value in fisher.items()},
    )
    return fisher, summary


def infer_batch_size(batch: Any) -> int:
    if isinstance(batch, torch.Tensor):
        return int(batch.shape[0])
    if isinstance(batch, Mapping):
        sizes = {int(value.shape[0]) for value in batch.values() if isinstance(value, torch.Tensor) and value.ndim}
        if len(sizes) != 1:
            raise EWCError(f"cannot infer one batch size from mapping: {sorted(sizes)}")
        return sizes.pop()
    if isinstance(batch, (tuple, list)) and batch:
        return infer_batch_size(batch[0])
    raise EWCError(f"cannot infer batch size from {type(batch)!r}")


class EWCRegularizer(nn.Module):
    """Named, masked product-Fisher EWC penalty registered as module buffers."""

    def __init__(
        self,
        anchors: TensorMap,
        masks: TensorMap,
        reference_fisher: TensorMap,
        control_fisher: TensorMap,
        *,
        normalization: str = "none",
    ) -> None:
        super().__init__()
        if normalization not in {"none", "block_mean_product"}:
            raise EWCError(f"unknown EWC normalization: {normalization}")
        self.normalization = normalization
        names = list(anchors)
        if not names or any(set(x) != set(names) for x in (masks, reference_fisher, control_fisher)):
            raise EWCError("anchor, mask, and Fisher parameter names must match exactly")
        self._parameter_names = names
        self._buffer_keys: dict[str, str] = {}
        for index, name in enumerate(names):
            tensors = [anchors[name], masks[name], reference_fisher[name], control_fisher[name]]
            if len({tuple(x.shape) for x in tensors}) != 1:
                raise EWCError(f"EWC tensor shape mismatch for {name}")
            if not all(torch.isfinite(x).all() for x in tensors):
                raise EWCError(f"non-finite EWC tensor for {name}")
            if torch.any(reference_fisher[name] < 0) or torch.any(control_fisher[name] < 0):
                raise EWCError(f"negative Fisher tensor for {name}")
            key = f"p{index}"
            self._buffer_keys[name] = key
            self.register_buffer(f"{key}_anchor", anchors[name].detach().clone())
            self.register_buffer(f"{key}_mask", masks[name].detach().clone())
            self.register_buffer(f"{key}_reference", reference_fisher[name].detach().clone())
            self.register_buffer(f"{key}_control", control_fisher[name].detach().clone())
            importance = masks[name] * reference_fisher[name] * control_fisher[name]
            positive = importance > 0
            scale = importance[positive].mean() if positive.any() else importance.new_tensor(1.0)
            active = positive.sum().clamp_min(1).to(dtype=importance.dtype)
            self.register_buffer(f"{key}_importance_scale", scale.detach().clone())
            self.register_buffer(f"{key}_active", active.detach().clone())

    def penalty(self, model: nn.Module) -> torch.Tensor:
        current = dict(model.named_parameters())
        if set(current) != set(self._parameter_names):
            missing = sorted(set(self._parameter_names) - set(current))
            extra = sorted(set(current) - set(self._parameter_names))
            raise EWCError(f"model parameters changed after anchoring; missing={missing}, extra={extra}")
        first = current[self._parameter_names[0]]
        total = torch.zeros((), device=first.device, dtype=first.dtype)
        contributing_blocks = 0
        for name in self._parameter_names:
            parameter = current[name]
            key = self._buffer_keys[name]
            anchor = getattr(self, f"{key}_anchor")
            mask = getattr(self, f"{key}_mask")
            reference = getattr(self, f"{key}_reference")
            control = getattr(self, f"{key}_control")
            importance = mask * reference * control
            contribution = importance * (parameter - anchor).square()
            if self.normalization == "block_mean_product" and torch.any(importance > 0):
                scale = getattr(self, f"{key}_importance_scale")
                active = getattr(self, f"{key}_active")
                total = total + (contribution / scale).sum() / active
                contributing_blocks += 1
            elif self.normalization == "none":
                total = total + contribution.sum()
        if self.normalization == "block_mean_product" and contributing_blocks:
            total = total / contributing_blocks
        return total

    def assert_zero_at_anchor(self, model: nn.Module, atol: float = 1e-12) -> None:
        value = float(self.penalty(model).detach().cpu())
        if abs(value) > atol:
            raise EWCError(f"EWC penalty is not zero at anchor: {value}")


def align_fisher_to_model(
    model: nn.Module, fisher: TensorMap, *, fill: float = 0.0
) -> dict[str, torch.Tensor]:
    """Require named tensors and expand old-shaped Fisher into current slices."""
    output: dict[str, torch.Tensor] = {}
    for name, parameter in model.named_parameters():
        value = torch.full_like(parameter, fill)
        if name in fisher:
            old = fisher[name].to(device=parameter.device, dtype=parameter.dtype)
            shared = _shared_slices(old, parameter)
            value[shared] = old
        output[name] = value
    extra = sorted(set(fisher) - set(output))
    if extra:
        raise EWCError(f"Fisher has parameters absent from model: {extra}")
    return output
