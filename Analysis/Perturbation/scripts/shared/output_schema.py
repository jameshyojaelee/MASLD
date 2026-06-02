"""Standardized output schema for perturbation predictions.

Every Model Runner subagent emits predictions in this schema. Peer Reviewer
validates schema compliance before passing to Consensus Aggregator.

Using pydantic for runtime validation — catches NaN, missing fields, malformed
predictions before they reach downstream consensus.

P0-X5 fix (2026-05-21): predictions now validated per-arm via dispatch in
ModelRunOutput.check_predictions_match_arm — was previously `list` (untyped),
so malformed dicts would silently pass through to consensus.
"""
from __future__ import annotations

from datetime import datetime
from typing import Any, Literal

from pydantic import BaseModel, Field, model_validator


Arm = Literal["D1", "D2", "D3", "D4", "D5"]
Direction = Literal["up", "down", "either", "—"]
Modality = Literal["zero_shot", "fine_tuned"]


class GenePrediction(BaseModel):
    """Per-gene downstream prediction (D1 mechanism, D5 mouse)."""

    target_gene: str = Field(description="HGNC (or mouse) symbol perturbed")
    downstream_gene: str = Field(description="Predicted downstream-affected gene")
    logFC_predicted: float
    abs_rank: int = Field(description="Rank within target's downstream list")
    direction: Direction
    confidence: float = Field(ge=0, le=1)


class ReversalPrediction(BaseModel):
    """Per-gene Diseased→Healthy reversal score (D2)."""

    gene: str
    reversal_score: float
    reversal_rank: int
    stage_specific: Literal["F0", "F1", "F2", "F3", "F4", "—"]
    cell_type: str  # "hepatocyte_progressor" etc.
    reference_signature: Literal["ref_a", "ref_b", "ref_c", "consensus"]


class SynergyPrediction(BaseModel):
    """Combinatorial KO synergy (D3)."""

    gene1: str
    gene2: str
    gene3: str | None = None
    gene4: str | None = None
    additive_baseline: float
    observed_double_or_higher: float
    synergy_magnitude: float  # = ||observed - additive||
    synergy_class: Literal["synergistic", "antagonistic", "epistatic", "additive"]
    sigma_above_additive: float
    k562_bias_confidence: float = Field(ge=0, le=1, description="0=high transfer risk")


class CircuitPrediction(BaseModel):
    """Hep→receiver paracrine perturbation (D4)."""

    sender_gene: str  # hepatocyte ligand
    receiver_cell_type: Literal["Macrophage", "Stellate", "LSEC", "Cholangiocyte"]
    receptor_gene: str
    receiver_response_genes: list[str]
    response_magnitude: float
    propagation_method: Literal["state_2stage", "niches", "commot"]


# Map arm → expected Prediction model. Used by ModelRunOutput validator.
ARM_TO_PRED_MODEL: dict[str, type[BaseModel]] = {}  # populated below


class ModelRunOutput(BaseModel):
    """Standardized envelope around any per-prediction list.

    Every Model Runner emits ONE ModelRunOutput per (arm × model × modality × context).

    Predictions are validated against the arm's expected schema (D1+D5 → GenePrediction,
    D2 → ReversalPrediction, D3 → SynergyPrediction, D4 → CircuitPrediction).
    Each prediction dict is round-tripped through pydantic to catch schema errors at
    write-time rather than at consensus-time.
    """

    arm: Arm
    model: str
    modality: Modality
    context: str = Field(description="Stratification context, e.g. 'F__Steatohepatitis__Progressor'")
    timestamp: datetime
    checkpoint_hash: str  # SHA256 of model weights (reproducibility)
    n_predictions: int
    predictions: list[dict[str, Any]]  # validated per-arm in check_predictions_match_arm
    runtime_seconds: float
    slurm_job_id: str | None = None
    notes: str = ""
    model_ejected: bool = Field(
        default=False,
        description="True if Backtester ejected this model from arm consensus (per F20 SOP). "
                    "Predictions remain in JSON for audit but downstream consensus skips them.",
    )

    @model_validator(mode="after")
    def check_predictions_match_arm(self):
        # Count check
        if self.n_predictions != len(self.predictions):
            raise ValueError(
                f"n_predictions={self.n_predictions} doesn't match len(predictions)={len(self.predictions)}"
            )
        # Per-arm schema validation
        pred_model = ARM_TO_PRED_MODEL.get(self.arm)
        if pred_model is None:
            raise ValueError(f"No prediction schema registered for arm={self.arm}")
        # Validate each item (raises ValidationError with line-level detail)
        for i, p in enumerate(self.predictions):
            try:
                pred_model.model_validate(p)
            except Exception as e:
                raise ValueError(
                    f"Prediction[{i}] failed {pred_model.__name__} validation: "
                    f"{type(e).__name__}: {str(e)[:200]}"
                ) from e
        return self


class ConsensusResult(BaseModel):
    """Per-gene consensus result from Consensus Aggregator."""

    arm: Arm
    gene: str
    n_models_agree: int
    n_models_total: int
    consensus_direction: Direction
    consensus_magnitude: float
    pass_consensus: bool
    contributing_models: list[str]
    disagreement_note: str = ""


# Populate the arm → prediction schema dispatch (must be at end of module so all
# Prediction classes are defined before registration).
ARM_TO_PRED_MODEL.update({
    "D1": GenePrediction,
    "D2": ReversalPrediction,
    "D3": SynergyPrediction,
    "D4": CircuitPrediction,
    "D5": GenePrediction,  # mouse target/downstream pair; same schema as D1
})
