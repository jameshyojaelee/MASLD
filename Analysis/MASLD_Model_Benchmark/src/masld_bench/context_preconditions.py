"""Outcome-blind preconditions for future sequence-plus-RNA conditioning.

This module does not build or authorize the conditional model.  It provides a
small executable requirements for fold-fitted context mapping, missingness, and
neutral FiLM/LoRA attachment behavior while the complementarity check is shut.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Mapping, Sequence

from .topology import MISSING_STATES, TopologyError, assert_pairing


class ContextPreconditionError(ValueError):
    """Raised when a precondition fixture could leak held or target evidence."""


MAPPER_ARMS = (
    "released_rank_quantile_from_counts",
    "length_adjusted_tpm_then_released_rank_quantile",
    "learned_count_to_rank_context_encoder",
)
REDUNDANT_RANK_ARMS = (
    "cpm_before_released_rank_quantile",
    "log_cpm_before_released_rank_quantile",
    "custom_monotone_quantile_before_released_rank_quantile",
)
_SPLIT_ROLES = frozenset({"train", "valid", "test", "sealed_query"})


@dataclass(frozen=True, slots=True)
class ContextVector:
    """RNA context with per-feature state and mask; missing values stay null."""

    values: tuple[float | None, ...]
    observed_mask: tuple[bool, ...]
    states: tuple[str, ...]

    def validate(self, width: int) -> None:
        if not (
            len(self.values) == len(self.observed_mask) == len(self.states) == width
        ):
            raise ContextPreconditionError("RNA context width differs")
        for value, mask, state in zip(
            self.values, self.observed_mask, self.states, strict=True
        ):
            if state not in MISSING_STATES:
                raise ContextPreconditionError(f"unknown RNA missingness state: {state}")
            if state == "observed":
                if mask is not True or value is None or isinstance(value, bool):
                    raise ContextPreconditionError(
                        "observed RNA features require a true mask and numeric value"
                    )
                numeric = float(value)
                if not isfinite(numeric) or numeric < 0:
                    raise ContextPreconditionError(
                        "observed RNA counts must be finite and nonnegative"
                    )
            elif mask is not False or value is not None:
                raise ContextPreconditionError(
                    "missing RNA features require false masks and null values, never zero"
                )

    def complete_values(self) -> tuple[float, ...]:
        if not all(self.observed_mask):
            raise ContextPreconditionError(
                "the released rank mapper requires complete observed RNA; route "
                "missing context to the explicit sequence-only path"
            )
        return tuple(float(value) for value in self.values if value is not None)


@dataclass(frozen=True, slots=True)
class ContextRow:
    row_id: str
    dataset_id: str
    donor_id: str
    study_id: str
    split_role: str
    genomic_block: str
    pairing_level: str
    rna: ContextVector
    observed_atac_input: bool = False
    sealed_outcome_input: bool = False

    def validate(self, width: int) -> None:
        if not all(
            isinstance(value, str) and value
            for value in (
                self.row_id,
                self.dataset_id,
                self.donor_id,
                self.study_id,
                self.genomic_block,
            )
        ):
            raise ContextPreconditionError("context row identifiers must be non-empty")
        if self.split_role not in _SPLIT_ROLES:
            raise ContextPreconditionError("context row split role differs")
        try:
            assert_pairing(self.dataset_id, self.pairing_level)
        except TopologyError as error:
            raise ContextPreconditionError(str(error)) from error
        if self.observed_atac_input:
            raise ContextPreconditionError(
                "observed ATAC may be a training target but never a context-model input"
            )
        if self.sealed_outcome_input:
            raise ContextPreconditionError("sealed outcomes may never enter context mapping")
        self.rna.validate(width)


@dataclass(frozen=True, slots=True)
class ContextMapperPrecondition:
    """Fold identity for a mapper; it contains no learned conditional model."""

    ordered_gene_ids: tuple[str, ...]
    reference_values: tuple[float, ...]
    mapper_arm: str
    training_donor_ids: frozenset[str]
    training_genomic_blocks: frozenset[str]
    receptive_field_bp: int
    boundary_buffer_bp: int

    @classmethod
    def fit(
        cls,
        rows: Sequence[ContextRow],
        *,
        ordered_gene_ids: Sequence[str],
        reference_values: Sequence[float],
        mapper_arm: str,
        held_donor_ids: Sequence[str],
        held_genomic_blocks: Sequence[str],
        receptive_field_bp: int,
        boundary_buffer_bp: int,
    ) -> "ContextMapperPrecondition":
        genes = tuple(ordered_gene_ids)
        reference = tuple(float(value) for value in reference_values)
        if (
            not rows
            or not genes
            or len(set(genes)) != len(genes)
            or len(reference) != len(genes)
            or any(not isfinite(value) for value in reference)
        ):
            raise ContextPreconditionError("context mapper axes differ")
        if mapper_arm not in MAPPER_ARMS:
            raise ContextPreconditionError("context mapper arm is not prespecified")
        if mapper_arm == "learned_count_to_rank_context_encoder":
            raise ContextPreconditionError(
                "learned context fitting remains blocked until conditional authorization"
            )
        if (
            not isinstance(receptive_field_bp, int)
            or receptive_field_bp <= 0
            or not isinstance(boundary_buffer_bp, int)
            or boundary_buffer_bp < receptive_field_bp
        ):
            raise ContextPreconditionError(
                "boundary buffer must cover the complete selected receptive field"
            )
        held_donors = frozenset(held_donor_ids)
        held_blocks = frozenset(held_genomic_blocks)
        if not held_donors or not held_blocks:
            raise ContextPreconditionError("crossed held donor and genomic blocks are required")
        for row in rows:
            row.validate(len(genes))
            if row.split_role != "train":
                raise ContextPreconditionError("context mapper fit received a nontraining row")
            if row.donor_id in held_donors or row.genomic_block in held_blocks:
                raise ContextPreconditionError("held donor or genomic block entered mapper fit")
        donors = frozenset(row.donor_id for row in rows)
        blocks = frozenset(row.genomic_block for row in rows)
        if donors & held_donors or blocks & held_blocks:
            raise ContextPreconditionError("context mapper split overlap differs")
        return cls(
            ordered_gene_ids=genes,
            reference_values=reference,
            mapper_arm=mapper_arm,
            training_donor_ids=donors,
            training_genomic_blocks=blocks,
            receptive_field_bp=receptive_field_bp,
            boundary_buffer_bp=boundary_buffer_bp,
        )

    def transform_query(
        self,
        row: ContextRow,
        *,
        gene_lengths: Sequence[float] | None = None,
    ) -> tuple[float, ...] | None:
        row.validate(len(self.ordered_gene_ids))
        if row.split_role == "train":
            raise ContextPreconditionError("query transform received a training row")
        if (
            row.donor_id in self.training_donor_ids
            or row.genomic_block in self.training_genomic_blocks
        ):
            raise ContextPreconditionError(
                "query must cross both held donor and held genomic block"
            )
        if not all(row.rna.observed_mask):
            return None
        values = row.rna.complete_values()
        if self.mapper_arm == "length_adjusted_tpm_then_released_rank_quantile":
            if gene_lengths is None or len(gene_lengths) != len(values):
                raise ContextPreconditionError("length-adjusted mapper lacks gene lengths")
            lengths = tuple(float(value) for value in gene_lengths)
            if any(not isfinite(value) or value <= 0 for value in lengths):
                raise ContextPreconditionError("gene lengths must be finite and positive")
            values = tuple(value / length for value, length in zip(values, lengths, strict=True))
        return stable_rank_quantile(values, self.reference_values)


def stable_rank_quantile(
    values: Sequence[float], reference_values: Sequence[float]
) -> tuple[float, ...]:
    """Match the released stable argsort/sorted-reference transform."""

    numeric = tuple(float(value) for value in values)
    reference = tuple(float(value) for value in reference_values)
    if len(numeric) != len(reference) or not numeric:
        raise ContextPreconditionError("rank-quantile axes differ")
    if any(not isfinite(value) for value in (*numeric, *reference)):
        raise ContextPreconditionError("rank-quantile inputs must be finite")
    order = sorted(range(len(numeric)), key=lambda index: numeric[index])
    normalized = [0.0] * len(numeric)
    for index, reference_value in zip(order, sorted(reference), strict=True):
        normalized[index] = reference_value
    return tuple(normalized)


def apply_film_lora_precondition(
    features: Sequence[Sequence[float]],
    *,
    gamma_delta: Sequence[float],
    beta: Sequence[float],
    lora_a: Sequence[Sequence[float]],
    lora_b: Sequence[Sequence[float]],
    alpha: float,
    context_observed: bool,
    gate: float,
) -> tuple[tuple[float, ...], ...]:
    """Apply a minimal shape fixture for residual FiLM plus a LoRA delta."""

    matrix = tuple(tuple(float(value) for value in row) for row in features)
    if not matrix or not matrix[0] or any(len(row) != len(matrix[0]) for row in matrix):
        raise ContextPreconditionError("sequence feature matrix differs")
    width = len(matrix[0])
    gamma = tuple(float(value) for value in gamma_delta)
    offsets = tuple(float(value) for value in beta)
    down = tuple(tuple(float(value) for value in row) for row in lora_a)
    up = tuple(tuple(float(value) for value in row) for row in lora_b)
    if (
        len(gamma) != width
        or len(offsets) != width
        or not down
        or any(len(row) != width for row in down)
        or len(up) != width
        or any(len(row) != len(down) for row in up)
    ):
        raise ContextPreconditionError("FiLM/LoRA attachment shapes differ")
    gate_value = float(gate)
    scale = float(alpha) / len(down)
    if not all(
        isfinite(value)
        for value in (
            gate_value,
            scale,
            *gamma,
            *offsets,
            *(value for row in matrix for value in row),
            *(value for row in down for value in row),
            *(value for row in up for value in row),
        )
    ):
        raise ContextPreconditionError("FiLM/LoRA fixture values must be finite")
    if not 0.0 <= gate_value <= 1.0:
        raise ContextPreconditionError("context gate must be in [0, 1]")
    if context_observed is False and gate_value != 0.0:
        raise ContextPreconditionError("missing context must force the context gate to zero")
    result = []
    for row in matrix:
        rank_state = [sum(a * value for a, value in zip(a_row, row, strict=True)) for a_row in down]
        delta = [
            scale * sum(weight * value for weight, value in zip(b_row, rank_state, strict=True))
            for b_row in up
        ]
        result.append(tuple(
            value * (1.0 + gate_value * gain)
            + gate_value * offset
            + gate_value * residual
            for value, gain, offset, residual in zip(
                row, gamma, offsets, delta, strict=True
            )
        ))
    return tuple(result)


def require_precondition_only_campaign(campaign: Mapping[str, object]) -> None:
    """Prove a fixture cannot turn the blocked conditional campaign executable."""

    if (
        campaign.get("wave") != "conditional_model"
        or campaign.get("status") != "blocked_trigger"
        or campaign.get("submit_enabled") is not False
        or campaign.get("conditional_model_spec_path") != "UNRESOLVED"
    ):
        raise ContextPreconditionError(
            "precondition fixtures require the unresolved, blocked conditional campaign"
        )


__all__ = [
    "ContextMapperPrecondition",
    "ContextPreconditionError",
    "ContextRow",
    "ContextVector",
    "MAPPER_ARMS",
    "REDUNDANT_RANK_ARMS",
    "apply_film_lora_precondition",
    "require_precondition_only_campaign",
    "stable_rank_quantile",
]
