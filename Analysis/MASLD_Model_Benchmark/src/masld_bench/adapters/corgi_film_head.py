"""Restricted Regular Corgi FiLM-plus-ATAC-head adaptation surface.

This module does not load checkpoints, datasets, or outcomes.  Callers must
strictly restore the released 22-track model first, then pass that model to
``CorgiFilmAtacModel``.  Only the two released FiLM MLPs and the replacement
one-channel ATAC head are trainable.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import torch
from torch import nn
from torch.nn import functional as F


CHECKPOINT_SHA256 = "cd51539f01de1e66a3a4f9b89e772bdeae07df2f0178d5b692f2368db2f4d2a4"
ATAC_CHANNEL_INDEX = 1
EXPECTED_TRAINABLE_NUMEL = 5_869_229
TRAINABLE_PREFIXES = (
    "film_mlp_conv.",
    "film_mlp_transformer.",
    "masld_atac_head.",
)
TRANSFERRED_MODULES = (
    "film_mlp_conv",
    "film_mlp_transformer",
    "film_layers",
    "max_pool",
    "conv_0",
    "conv_1",
    "conv_2",
    "conv_3",
    "conv_4",
    "conv_5",
    "transformer_layers",
    "crop",
    "final_conv",
)


class CorgiFilmHeadError(ValueError):
    """Raised when the restricted adaptation surface differs."""


class CorgiFilmAtacModel(nn.Module):
    """Corgi backbone with top-level FiLM MLPs and one softplus ATAC head."""

    def __init__(self, released_model: nn.Module, *, require_exact: bool = True) -> None:
        super().__init__()
        config = getattr(released_model, "config", None)
        if not isinstance(config, dict):
            raise CorgiFilmHeadError("released model config differs")
        self.config = dict(config)
        if require_exact and (
            self.config.get("input_trans_regulators") != 2_891
            or self.config.get("output_channels") != 22
            or self.config.get("dim") != 1_536
            or self.config.get("output_central_bins") != 6_144
        ):
            raise CorgiFilmHeadError("released Corgi dimensions differ")

        for name in TRANSFERRED_MODULES:
            module = getattr(released_model, name, None)
            if not isinstance(module, nn.Module):
                raise CorgiFilmHeadError(f"released module is absent: {name}")
            delattr(released_model, name)
            self.add_module(name, module)

        pretrained_head = getattr(released_model, "output_head", None)
        if not isinstance(pretrained_head, nn.Conv1d):
            raise CorgiFilmHeadError("released regular Corgi head is not direct Conv1d")
        if tuple(pretrained_head.weight.shape) != (22, 1_920, 1) or tuple(
            pretrained_head.bias.shape
        ) != (22,):
            raise CorgiFilmHeadError("released regular Corgi head shape differs")
        self.masld_atac_head = nn.Conv1d(1_920, 1, kernel_size=1, bias=True)
        with torch.no_grad():
            self.masld_atac_head.weight.copy_(
                pretrained_head.weight[ATAC_CHANNEL_INDEX : ATAC_CHANNEL_INDEX + 1]
            )
            self.masld_atac_head.bias.copy_(
                pretrained_head.bias[ATAC_CHANNEL_INDEX : ATAC_CHANNEL_INDEX + 1]
            )
        delattr(released_model, "output_head")
        self._configure_trainable_surface(require_exact=require_exact)

    @staticmethod
    def _split_film_params(
        film_params: torch.Tensor, film_dimensions: list[int]
    ) -> tuple[dict[int, torch.Tensor], dict[int, torch.Tensor]]:
        offset = 0
        scales: dict[int, torch.Tensor] = {}
        shifts: dict[int, torch.Tensor] = {}
        for index, channel_dim in enumerate(film_dimensions):
            scales[index] = film_params[:, offset : offset + channel_dim]
            shifts[index] = film_params[
                :, offset + channel_dim : offset + 2 * channel_dim
            ]
            offset += 2 * channel_dim
        if offset != film_params.shape[1]:
            raise CorgiFilmHeadError("FiLM parameter width differs")
        return scales, shifts

    def _convolutions(
        self,
        sequence: torch.Tensor,
        scales: Mapping[int, torch.Tensor],
        shifts: Mapping[int, torch.Tensor],
    ) -> torch.Tensor:
        value = self.conv_0(sequence)
        value = self.max_pool(value)
        value = self.film_layers[0](value, scales[0], shifts[0])
        value = self.conv_1(value)
        value = self.max_pool(value)
        value = self.conv_2(value)
        value = self.max_pool(value)
        value = self.film_layers[1](value, scales[1], shifts[1])
        value = self.conv_3(value)
        value = self.max_pool(value)
        value = self.conv_4(value)
        value = self.max_pool(value)
        value = self.conv_5(value)
        return self.max_pool(value)

    @staticmethod
    def _transformer_film(
        value: torch.Tensor,
        transformer_layer: nn.Module,
        film_layer: nn.Module,
        scale: torch.Tensor,
        shift: torch.Tensor,
    ) -> torch.Tensor:
        attended = transformer_layer.block[0](value)
        feed_forward = transformer_layer.block[1].fn
        value = feed_forward[0](attended)
        value = film_layer(value.permute(0, 2, 1), scale, shift).permute(0, 2, 1)
        for layer in feed_forward[1:]:
            value = layer(value)
        return value + attended

    def _transformers(
        self,
        value: torch.Tensor,
        scales: Mapping[int, torch.Tensor],
        shifts: Mapping[int, torch.Tensor],
    ) -> torch.Tensor:
        value = self.transformer_layers[0](value)
        value = self.transformer_layers[1](value)
        value = self._transformer_film(
            value,
            self.transformer_layers[2],
            self.film_layers[2],
            scales[0],
            shifts[0],
        )
        value = self.transformer_layers[3](value)
        value = self.transformer_layers[4](value)
        value = self._transformer_film(
            value,
            self.transformer_layers[5],
            self.film_layers[3],
            scales[1],
            shifts[1],
        )
        value = self.transformer_layers[6](value)
        value = self.transformer_layers[7](value)
        return self.transformer_layers[8](value)

    def forward_raw(
        self, sequence: torch.Tensor, trans_regulator_context: torch.Tensor
    ) -> torch.Tensor:
        conv_parameters = self.film_mlp_conv(trans_regulator_context)
        conv_scales, conv_shifts = self._split_film_params(
            conv_parameters, list(self.config["film_dimensions_conv"])
        )
        transformer_parameters = self.film_mlp_transformer(
            trans_regulator_context
        )
        transformer_scales, transformer_shifts = self._split_film_params(
            transformer_parameters,
            list(self.config["film_dimensions_transformer"]),
        )
        value = self._convolutions(
            sequence.permute(0, 2, 1), conv_scales, conv_shifts
        )
        value = self._transformers(
            value.permute(0, 2, 1), transformer_scales, transformer_shifts
        )
        value = self.crop(value.permute(0, 2, 1))
        value = self.final_conv(value)
        return self.masld_atac_head(value)

    def forward(
        self, sequence: torch.Tensor, trans_regulator_context: torch.Tensor
    ) -> torch.Tensor:
        return F.softplus(self.forward_raw(sequence, trans_regulator_context))

    def _configure_trainable_surface(self, *, require_exact: bool) -> None:
        self.eval()
        for parameter in self.parameters():
            parameter.requires_grad = False
        for name, parameter in self.named_parameters():
            if name.startswith(TRAINABLE_PREFIXES):
                parameter.requires_grad = True
        manifest = trainable_parameter_manifest(self)
        if require_exact and manifest["trainable_parameter_numel"] != EXPECTED_TRAINABLE_NUMEL:
            raise CorgiFilmHeadError("exact trainable parameter count differs")
        if any(module.training for module in self.modules()):
            raise CorgiFilmHeadError("adapted model must remain in eval mode")


def trainable_parameter_manifest(model: nn.Module) -> dict[str, Any]:
    trainable: list[dict[str, Any]] = []
    frozen_numel = 0
    for name, parameter in model.named_parameters():
        if parameter.requires_grad:
            if not name.startswith(TRAINABLE_PREFIXES):
                raise CorgiFilmHeadError(f"unexpected trainable parameter: {name}")
            trainable.append(
                {
                    "name": name,
                    "shape": list(parameter.shape),
                    "numel": int(parameter.numel()),
                    "dtype": str(parameter.dtype),
                }
            )
        else:
            frozen_numel += int(parameter.numel())
    observed_prefixes = {
        prefix for prefix in TRAINABLE_PREFIXES if any(row["name"].startswith(prefix) for row in trainable)
    }
    if observed_prefixes != set(TRAINABLE_PREFIXES):
        raise CorgiFilmHeadError("one or more trainable parameter families are absent")
    batchnorm_buffers = [
        name
        for name, _buffer in model.named_buffers()
        if name.endswith(("running_mean", "running_var", "num_batches_tracked"))
    ]
    return {
        "schema_version": "masld-bench-corgi-film-head-parameter-manifest-v1",
        "trainable_prefixes": list(TRAINABLE_PREFIXES),
        "trainable_parameters": trainable,
        "trainable_parameter_tensors": len(trainable),
        "trainable_parameter_numel": sum(row["numel"] for row in trainable),
        "frozen_parameter_numel": frozen_numel,
        "batchnorm_buffer_names": batchnorm_buffers,
        "batchnorm_running_statistics_frozen_by_eval_mode": not any(
            module.training for module in model.modules()
        ),
    }


def build_optimizer(
    model: nn.Module,
    *,
    film_learning_rate: float = 1.0e-5,
    head_learning_rate: float = 3.0e-5,
    film_weight_decay: float = 1.0e-3,
    head_weight_decay: float = 1.0e-4,
) -> torch.optim.AdamW:
    if min(film_learning_rate, head_learning_rate) <= 0 or min(
        film_weight_decay, head_weight_decay
    ) < 0:
        raise CorgiFilmHeadError("optimizer hyperparameters differ")
    film: list[nn.Parameter] = []
    head: list[nn.Parameter] = []
    for name, parameter in model.named_parameters():
        if not parameter.requires_grad:
            continue
        if name.startswith(("film_mlp_conv.", "film_mlp_transformer.")):
            film.append(parameter)
        elif name.startswith("masld_atac_head."):
            head.append(parameter)
        else:
            raise CorgiFilmHeadError(f"optimizer received unexpected parameter: {name}")
    if not film or not head:
        raise CorgiFilmHeadError("optimizer parameter groups are incomplete")
    return torch.optim.AdamW(
        [
            {
                "params": film,
                "lr": film_learning_rate,
                "weight_decay": film_weight_decay,
                "group_name": "film_mlp",
            },
            {
                "params": head,
                "lr": head_learning_rate,
                "weight_decay": head_weight_decay,
                "group_name": "atac_head",
            },
        ]
    )


def adaptation_delta(model: nn.Module) -> dict[str, Any]:
    state = model.state_dict()
    names = [
        name
        for name, parameter in model.named_parameters()
        if parameter.requires_grad
    ]
    return {
        "schema_version": "masld-bench-corgi-film-head-delta-v1",
        "source_checkpoint_sha256": CHECKPOINT_SHA256,
        "trainable_parameter_names": names,
        "trainable_state_dict": {
            name: state[name].detach().cpu().clone() for name in names
        },
    }


def load_adaptation_delta(model: nn.Module, delta: Mapping[str, Any]) -> None:
    if (
        delta.get("schema_version") != "masld-bench-corgi-film-head-delta-v1"
        or delta.get("source_checkpoint_sha256") != CHECKPOINT_SHA256
    ):
        raise CorgiFilmHeadError("adaptation delta identity differs")
    expected = [
        name
        for name, parameter in model.named_parameters()
        if parameter.requires_grad
    ]
    if delta.get("trainable_parameter_names") != expected:
        raise CorgiFilmHeadError("adaptation delta parameter roster differs")
    values = delta.get("trainable_state_dict")
    if not isinstance(values, Mapping) or set(values) != set(expected):
        raise CorgiFilmHeadError("adaptation delta tensor roster differs")
    parameters = dict(model.named_parameters())
    with torch.no_grad():
        for name in expected:
            value = values[name]
            if not isinstance(value, torch.Tensor) or value.shape != parameters[name].shape:
                raise CorgiFilmHeadError(f"adaptation delta shape differs: {name}")
            parameters[name].copy_(value.to(parameters[name]))


__all__ = [
    "ATAC_CHANNEL_INDEX",
    "CHECKPOINT_SHA256",
    "CorgiFilmAtacModel",
    "CorgiFilmHeadError",
    "EXPECTED_TRAINABLE_NUMEL",
    "TRAINABLE_PREFIXES",
    "adaptation_delta",
    "build_optimizer",
    "load_adaptation_delta",
    "trainable_parameter_manifest",
]
