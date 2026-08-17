"""Narrow scvi-tools 1.3.3 adapters for audited continual-learning losses."""

from __future__ import annotations

import platform
from typing import Any

import torch

from .ewc import EWCRegularizer
from .distillation import LatentDistillationRegularizer

try:
    import scvi
    from scvi.train import SemiSupervisedTrainingPlan, TrainingPlan
except ImportError as exc:  # pragma: no cover - exercised only in incomplete environments
    raise ImportError("masld_cl.scvi_adapter requires pinned scvi-tools==1.3.3") from exc


def assert_scvi_version() -> None:
    if scvi.__version__ != "1.3.3":
        raise RuntimeError(f"production port requires scvi-tools==1.3.3, found {scvi.__version__}")


def assert_runtime_versions(config: dict[str, Any]) -> None:
    """Fail closed when the production runtime differs from the lock."""
    assert_scvi_version()
    expected = config["software"]
    python = ".".join(platform.python_version().split(".")[:2])
    torch_version = torch.__version__.split("+", 1)[0]
    observed = {
        "python": python,
        "scvi_tools": scvi.__version__,
        "torch": torch_version,
    }
    mismatches = {
        key: {"expected": str(expected[key]), "observed": value}
        for key, value in observed.items() if value != str(expected[key])
    }
    if mismatches:
        raise RuntimeError(f"production runtime differs from config: {mismatches}")


class _EWCPlanMixin:
    def _initialize_ewc(self, regularizer: EWCRegularizer, ewc_lambda: float) -> None:
        if ewc_lambda < 0:
            raise ValueError("ewc_lambda must be non-negative")
        self.ewc_regularizer = regularizer
        self.ewc_lambda = float(ewc_lambda)
        self.ewc_regularizer.assert_zero_at_anchor(self.module)

    def _add_ewc(self, base_loss: torch.Tensor) -> torch.Tensor:
        penalty = self.ewc_regularizer.penalty(self.module)
        total = base_loss + self.ewc_lambda * penalty
        self.log("ewc_penalty", penalty, on_step=False, on_epoch=True, batch_size=1)
        self.log("objective_with_ewc", total, on_step=False, on_epoch=True, batch_size=1)
        return total


class EWCTrainingPlan(_EWCPlanMixin, TrainingPlan):
    """scVI training plan with the audited EWC term."""

    def __init__(self, module, *, ewc_regularizer: EWCRegularizer, ewc_lambda: float, **kwargs):
        super().__init__(module, **kwargs)
        self._initialize_ewc(ewc_regularizer, ewc_lambda)

    def training_step(self, batch, batch_idx):
        return self._add_ewc(super().training_step(batch, batch_idx))


class EWCSemiSupervisedTrainingPlan(_EWCPlanMixin, SemiSupervisedTrainingPlan):
    """scANVI training plan with the audited EWC term."""

    def __init__(
        self, module, n_classes: int, *, ewc_regularizer: EWCRegularizer,
        ewc_lambda: float, **kwargs,
    ):
        # SCANVI.train passes n_classes positionally in scvi-tools 1.3.3.
        super().__init__(module, n_classes, **kwargs)
        self._initialize_ewc(ewc_regularizer, ewc_lambda)

    def training_step(self, batch, batch_idx):
        return self._add_ewc(super().training_step(batch, batch_idx))


def scvi_reconstruction_loss(module: torch.nn.Module, batch: Any) -> torch.Tensor:
    """Return the ordinary unlabeled scVI/scANVI objective for Fisher estimation."""
    device = next(module.parameters()).device
    batch = {
        key: value.to(device) if isinstance(value, torch.Tensor) else value
        for key, value in batch.items()
    }
    _, _, output = module(batch)
    return output.loss


def attach_training_plan(model: Any, *, semi_supervised: bool) -> None:
    """Set the instance's training-plan class without modifying scvi internals."""
    assert_scvi_version()
    model._training_plan_cls = (
        EWCSemiSupervisedTrainingPlan if semi_supervised else EWCTrainingPlan
    )


class _LatentDistillationPlanMixin:
    def _initialize_distillation(
        self, regularizer: LatentDistillationRegularizer,
        distillation_weight: float,
    ) -> None:
        if distillation_weight <= 0:
            raise ValueError("distillation_weight must be positive")
        self.distillation_regularizer = regularizer
        self.distillation_weight = float(distillation_weight)

    def _add_distillation(self, base_loss: torch.Tensor, batch) -> torch.Tensor:
        penalty = self.distillation_regularizer.penalty(self.module, batch)
        total = base_loss + self.distillation_weight * penalty
        self.log(
            "latent_distillation_penalty", penalty, on_step=False, on_epoch=True,
            batch_size=1,
        )
        self.log(
            "objective_with_latent_distillation", total, on_step=False, on_epoch=True,
            batch_size=1,
        )
        return total


class LatentDistillationTrainingPlan(_LatentDistillationPlanMixin, TrainingPlan):
    """scVI training plan with replay-cell posterior-mean anchoring."""

    def __init__(
        self, module, *, distillation_regularizer: LatentDistillationRegularizer,
        distillation_weight: float, **kwargs,
    ):
        super().__init__(module, **kwargs)
        self._initialize_distillation(distillation_regularizer, distillation_weight)

    def training_step(self, batch, batch_idx):
        return self._add_distillation(
            super().training_step(batch, batch_idx), batch
        )


class LatentDistillationSemiSupervisedTrainingPlan(
    _LatentDistillationPlanMixin, SemiSupervisedTrainingPlan
):
    """scANVI training plan with replay-cell posterior-mean anchoring."""

    def __init__(
        self, module, n_classes: int, *,
        distillation_regularizer: LatentDistillationRegularizer,
        distillation_weight: float, **kwargs,
    ):
        super().__init__(module, n_classes, **kwargs)
        self._initialize_distillation(distillation_regularizer, distillation_weight)

    def training_step(self, batch, batch_idx):
        return self._add_distillation(
            super().training_step(batch, batch_idx), batch
        )


def attach_distillation_training_plan(model: Any, *, semi_supervised: bool) -> None:
    """Attach the audited latent-distillation plan to one model instance."""
    assert_scvi_version()
    model._training_plan_cls = (
        LatentDistillationSemiSupervisedTrainingPlan
        if semi_supervised else LatentDistillationTrainingPlan
    )


class EWCAndLatentDistillationTrainingPlan(
    _EWCPlanMixin, _LatentDistillationPlanMixin, TrainingPlan
):
    """scVI plan with block-normalized EWC and replay-latent anchoring."""

    def __init__(
        self, module, *, ewc_regularizer: EWCRegularizer, ewc_lambda: float,
        distillation_regularizer: LatentDistillationRegularizer,
        distillation_weight: float, **kwargs,
    ):
        super().__init__(module, **kwargs)
        self._initialize_ewc(ewc_regularizer, ewc_lambda)
        self._initialize_distillation(distillation_regularizer, distillation_weight)

    def training_step(self, batch, batch_idx):
        base = super().training_step(batch, batch_idx)
        return self._add_ewc(self._add_distillation(base, batch))


class EWCAndLatentDistillationSemiSupervisedTrainingPlan(
    _EWCPlanMixin, _LatentDistillationPlanMixin, SemiSupervisedTrainingPlan
):
    """scANVI plan with block-normalized EWC and replay-latent anchoring."""

    def __init__(
        self, module, n_classes: int, *, ewc_regularizer: EWCRegularizer,
        ewc_lambda: float, distillation_regularizer: LatentDistillationRegularizer,
        distillation_weight: float, **kwargs,
    ):
        super().__init__(module, n_classes, **kwargs)
        self._initialize_ewc(ewc_regularizer, ewc_lambda)
        self._initialize_distillation(distillation_regularizer, distillation_weight)

    def training_step(self, batch, batch_idx):
        base = super().training_step(batch, batch_idx)
        return self._add_ewc(self._add_distillation(base, batch))


def attach_ewc_distillation_training_plan(model: Any, *, semi_supervised: bool) -> None:
    """Attach the audited combined V13 plan to one model instance."""
    assert_scvi_version()
    model._training_plan_cls = (
        EWCAndLatentDistillationSemiSupervisedTrainingPlan
        if semi_supervised else EWCAndLatentDistillationTrainingPlan
    )
