#!/usr/bin/env python3
"""Small task-native sequence controls for GSE281364 MPRA transfer."""

from __future__ import annotations

from typing import Mapping

import torch
from torch import nn
from torch.nn import functional as F


INPUT_LENGTH = 230
CHANNEL_ORDER = ("A", "G", "C", "T")
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
MODEL_IDS = ("sequence_cnn_control", "sequence_transformer_control")


class SequenceControlError(ValueError):
    """Raised when input geometry or model identity differs."""


def reverse_complement_one_hot(values: torch.Tensor) -> torch.Tensor:
    if values.ndim != 3 or values.shape[1:] != (4, INPUT_LENGTH):
        raise SequenceControlError("one-hot tensor geometry differs")
    return values[:, (3, 2, 1, 0), :].flip(-1)


def _validate_input(values: torch.Tensor) -> None:
    if (
        values.ndim != 3
        or values.shape[1:] != (4, INPUT_LENGTH)
        or not torch.isfinite(values).all()
    ):
        raise SequenceControlError("sequence input differs")


class ResidualDilatedBlock(nn.Module):
    def __init__(self, channels: int, kernel: int, dilation: int, dropout: float):
        super().__init__()
        padding = dilation * (kernel // 2)
        self.norm1 = nn.GroupNorm(8, channels)
        self.conv1 = nn.Conv1d(
            channels, channels, kernel, padding=padding, dilation=dilation
        )
        self.norm2 = nn.GroupNorm(8, channels)
        self.conv2 = nn.Conv1d(
            channels, channels, kernel, padding=padding, dilation=dilation
        )
        self.dropout = nn.Dropout(dropout)

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        residual = values
        values = self.conv1(F.gelu(self.norm1(values)))
        values = self.dropout(values)
        values = self.conv2(F.gelu(self.norm2(values)))
        return residual + self.dropout(values)


class SmallSequenceCNN(nn.Module):
    def __init__(self, dropout: float = 0.1):
        super().__init__()
        self.stem = nn.Conv1d(4, 64, 15, padding=7)
        self.blocks = nn.ModuleList(
            ResidualDilatedBlock(64, 7, dilation, dropout)
            for dilation in (1, 2, 4, 8)
        )
        self.final_norm = nn.GroupNorm(8, 64)
        self.head = nn.Sequential(
            nn.Linear(128, 64), nn.GELU(), nn.Dropout(dropout), nn.Linear(64, 2)
        )

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        _validate_input(values)
        values = self.stem(values)
        for block in self.blocks:
            values = block(values)
        values = F.gelu(self.final_norm(values))
        pooled = torch.cat((values.mean(-1), values.amax(-1)), dim=1)
        return self.head(pooled)


class SmallSequenceTransformer(nn.Module):
    def __init__(self, dropout: float = 0.1):
        super().__init__()
        width = 64
        self.token_projection = nn.Linear(4, width)
        self.cls_token = nn.Parameter(torch.zeros(1, 1, width))
        self.position = nn.Parameter(torch.zeros(1, INPUT_LENGTH + 1, width))
        layer = nn.TransformerEncoderLayer(
            d_model=width,
            nhead=4,
            dim_feedforward=192,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(layer, num_layers=6)
        self.final_norm = nn.LayerNorm(width)
        self.head = nn.Linear(width, 2)
        nn.init.normal_(self.cls_token, std=0.02)
        nn.init.normal_(self.position, std=0.02)

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        _validate_input(values)
        tokens = self.token_projection(values.transpose(1, 2))
        cls = self.cls_token.expand(tokens.shape[0], -1, -1)
        tokens = torch.cat((cls, tokens), dim=1) + self.position
        encoded = self.encoder(tokens)
        return self.head(self.final_norm(encoded[:, 0]))


def build_model(model_id: str, *, seed: int, dropout: float = 0.1) -> nn.Module:
    if model_id not in MODEL_IDS or seed < 0 or not 0.0 <= dropout < 1.0:
        raise SequenceControlError("model identity, seed, or dropout differs")
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    if model_id == "sequence_cnn_control":
        return SmallSequenceCNN(dropout=dropout)
    return SmallSequenceTransformer(dropout=dropout)


def allele_pair_activity(
    model: nn.Module,
    reference: torch.Tensor,
    alternative: torch.Tensor,
    *,
    reverse_complement_average: bool,
) -> Mapping[str, torch.Tensor]:
    _validate_input(reference)
    _validate_input(alternative)
    if reference.shape != alternative.shape:
        raise SequenceControlError("REF and ALT batch geometry differs")
    reference_forward = model(reference)
    alternative_forward = model(alternative)
    if reverse_complement_average:
        reference_reverse = model(reverse_complement_one_hot(reference))
        alternative_reverse = model(reverse_complement_one_hot(alternative))
        reference_mean = 0.5 * (reference_forward + reference_reverse)
        alternative_mean = 0.5 * (alternative_forward + alternative_reverse)
    else:
        reference_reverse = torch.full_like(reference_forward, torch.nan)
        alternative_reverse = torch.full_like(alternative_forward, torch.nan)
        reference_mean = reference_forward
        alternative_mean = alternative_forward
    return {
        "reference_forward": reference_forward,
        "alternative_forward": alternative_forward,
        "reference_reverse_complement": reference_reverse,
        "alternative_reverse_complement": alternative_reverse,
        "reference_orientation_mean": reference_mean,
        "alternative_orientation_mean": alternative_mean,
        "alt_minus_ref": alternative_mean - reference_mean,
    }


def joint_orientation_augmentation(
    reference: torch.Tensor,
    alternative: torch.Tensor,
    reverse_mask: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    _validate_input(reference)
    _validate_input(alternative)
    if (
        reference.shape != alternative.shape
        or reverse_mask.shape != (reference.shape[0],)
        or reverse_mask.dtype != torch.bool
    ):
        raise SequenceControlError("joint orientation mask differs")
    ref_rc = reverse_complement_one_hot(reference)
    alt_rc = reverse_complement_one_hot(alternative)
    mask = reverse_mask[:, None, None]
    return torch.where(mask, ref_rc, reference), torch.where(mask, alt_rc, alternative)
