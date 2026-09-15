"""Biology-free RNA-context FiLM/LoRA fixture for the blocked conditional slot.

This is not a model implementation. It has no selected parent, checkpoint,
trainable state, biological input, or GPU path. The fixture only proves shape,
mask, ablation, permutation, and JSON-export behavior while the conditional
complementarity trigger remains closed.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass
import hashlib
import json
from math import isfinite
from pathlib import Path
from typing import Mapping, Sequence

from masld_bench.context_preconditions import (
    ContextPreconditionError,
    apply_film_lora_precondition,
)
from masld_bench.topology import MISSING_STATES


FIXTURE_SCHEMA_VERSION = "masld-bench-context-borzoi-rna-adapter-fixture-v1"
OUTPUT_SCHEMA_VERSION = "masld-bench-context-borzoi-rna-adapter-output-v1"
REQUIRED_ARMS = (
    "conditioned",
    "sequence_only",
    "trans_only",
    "permuted_context",
)
PROHIBITED_EXPORT_SUFFIXES = {
    ".ckpt",
    ".h5",
    ".hdf5",
    ".joblib",
    ".pkl",
    ".pickle",
    ".pt",
    ".pth",
}


class ContextBorzoiFixtureError(ValueError):
    """Raised when the synthetic adapter fixture crosses a closed boundary."""


def _canonical_json(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True, slots=True)
class RNAContextMask:
    """Synthetic RNA values with typed feature-level missingness."""

    values: tuple[float | None, ...]
    observed_mask: tuple[bool, ...]
    states: tuple[str, ...]

    def validate(self, width: int) -> None:
        if width <= 0 or not (
            len(self.values) == len(self.observed_mask) == len(self.states) == width
        ):
            raise ContextBorzoiFixtureError("RNA context axes differ")
        for value, observed, state in zip(
            self.values, self.observed_mask, self.states, strict=True
        ):
            if state not in MISSING_STATES:
                raise ContextBorzoiFixtureError(f"unknown missingness state: {state}")
            if state == "observed":
                if observed is not True or value is None or isinstance(value, bool):
                    raise ContextBorzoiFixtureError(
                        "observed RNA requires a true mask and numeric value"
                    )
                numeric = float(value)
                if not isfinite(numeric) or numeric < 0:
                    raise ContextBorzoiFixtureError(
                        "observed RNA fixture values must be finite and nonnegative"
                    )
            elif observed is not False or value is not None:
                raise ContextBorzoiFixtureError(
                    "missing RNA requires a false mask and null value, never zero"
                )

    @property
    def complete(self) -> bool:
        return bool(self.observed_mask) and all(self.observed_mask)

    def summary(self) -> float:
        if not self.complete:
            raise ContextBorzoiFixtureError(
                "incomplete RNA must route to the sequence-only identity"
            )
        numeric = tuple(float(value) for value in self.values if value is not None)
        weights = tuple(range(1, len(numeric) + 1))
        bounded = tuple(value / (1.0 + value) for value in numeric)
        return sum(weight * value for weight, value in zip(weights, bounded, strict=True)) / sum(
            weights
        )


@dataclass(frozen=True, slots=True)
class SyntheticContextRow:
    """A nonbiological row carrying only RNA context and permutation strata."""

    row_id: str
    outer_partition: str
    dataset_id: str
    assay: str
    lineage: str
    topology: str
    rna: RNAContextMask
    source_kind: str = "synthetic_fixture_only"
    observed_atac_input: bool = False
    biological_data_input: bool = False
    sealed_input: bool = False

    def validate(self, width: int) -> None:
        identifiers = (
            self.row_id,
            self.outer_partition,
            self.dataset_id,
            self.assay,
            self.lineage,
            self.topology,
        )
        if not all(isinstance(value, str) and value for value in identifiers):
            raise ContextBorzoiFixtureError("synthetic context identifiers must be nonempty")
        if self.source_kind != "synthetic_fixture_only":
            raise ContextBorzoiFixtureError("fixture cannot consume a biological source")
        if self.observed_atac_input:
            raise ContextBorzoiFixtureError(
                "observed ATAC belongs to the separate observed-multiome task"
            )
        if self.biological_data_input or self.sealed_input:
            raise ContextBorzoiFixtureError("biological or sealed fixture input is forbidden")
        self.rna.validate(width)

    def permutation_stratum(self) -> tuple[object, ...]:
        return (
            self.outer_partition,
            self.dataset_id,
            self.assay,
            self.lineage,
            self.topology,
            self.rna.states,
            self.rna.observed_mask,
        )


@dataclass(frozen=True, slots=True)
class AdapterFixtureSpec:
    """Small synthetic geometry, never a parent-bound model specification."""

    feature_width: int = 4
    lora_rank: int = 2
    lora_alpha: float = 2.0
    permutation_seed: int = 20260825
    complementarity_trigger_state: str = "NOT_EVALUATED"
    parent_component: str = "UNSELECTED"
    checkpoint_identity: str = "NONE"
    model_weights_present: bool = False
    architecture_built: bool = False
    training_authorized: bool = False
    prediction_authorized: bool = False
    gpu_execution_authorized: bool = False

    def validate(self) -> None:
        if (
            not isinstance(self.feature_width, int)
            or self.feature_width <= 0
            or not isinstance(self.lora_rank, int)
            or not 0 < self.lora_rank <= self.feature_width
            or not isfinite(float(self.lora_alpha))
            or self.lora_alpha <= 0
            or not isinstance(self.permutation_seed, int)
        ):
            raise ContextBorzoiFixtureError("synthetic adapter geometry differs")
        if (
            self.complementarity_trigger_state != "NOT_EVALUATED"
            or self.parent_component != "UNSELECTED"
            or self.checkpoint_identity != "NONE"
            or self.model_weights_present
            or self.architecture_built
            or self.training_authorized
            or self.prediction_authorized
            or self.gpu_execution_authorized
        ):
            raise ContextBorzoiFixtureError(
                "synthetic fixture requires an unselected parent and closed trigger"
            )


@dataclass(frozen=True, slots=True)
class AdapterFixtureOutput:
    schema_version: str
    arm: str
    effective_arm: str
    target_row_id: str
    context_source_row_id: str | None
    context_observed_mask: tuple[bool, ...] | None
    context_states: tuple[str, ...] | None
    values: tuple[tuple[float, ...], ...]
    sequence_values_consumed: bool
    rna_context_values_consumed: bool
    observed_atac_consumed: bool
    biological_data_consumed: bool
    model_weights_consumed: bool
    gpu_executed: bool

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def derange_context_assignments(
    rows: Sequence[SyntheticContextRow],
    *,
    width: int,
    seed: int,
) -> dict[str, SyntheticContextRow]:
    """Derange contexts within every complete, identical permutation stratum."""

    if not rows or len({row.row_id for row in rows}) != len(rows):
        raise ContextBorzoiFixtureError("permutation rows are empty or duplicated")
    groups: dict[tuple[object, ...], list[SyntheticContextRow]] = defaultdict(list)
    for row in rows:
        row.validate(width)
        if not row.rna.complete:
            raise ContextBorzoiFixtureError(
                "permuted-context fixture requires a complete shared missingness pattern"
            )
        groups[row.permutation_stratum()].append(row)
    result: dict[str, SyntheticContextRow] = {}
    for group in groups.values():
        if len(group) < 2:
            raise ContextBorzoiFixtureError(
                "each permutation stratum needs at least two synthetic rows"
            )
        ordered = sorted(
            group,
            key=lambda row: hashlib.sha256(
                f"{seed}\0{row.row_id}".encode("utf-8")
            ).hexdigest(),
        )
        for index, target in enumerate(ordered):
            source = ordered[(index + 1) % len(ordered)]
            if source.row_id == target.row_id:
                raise ContextBorzoiFixtureError("permutation contains a fixed point")
            result[target.row_id] = source
    return result


def _validate_features(
    features: Sequence[Sequence[float]], width: int
) -> tuple[tuple[float, ...], ...]:
    matrix = tuple(tuple(float(value) for value in row) for row in features)
    if (
        not matrix
        or any(len(row) != width for row in matrix)
        or any(not isfinite(value) for row in matrix for value in row)
    ):
        raise ContextBorzoiFixtureError("synthetic sequence feature geometry differs")
    return matrix


def _synthetic_adapter_terms(
    context: RNAContextMask,
    spec: AdapterFixtureSpec,
) -> tuple[
    tuple[float, ...],
    tuple[float, ...],
    tuple[tuple[float, ...], ...],
    tuple[tuple[float, ...], ...],
]:
    """Create deterministic coefficients, not learned or exported model weights."""

    summary = context.summary()
    gamma = tuple(summary * 0.02 * (index + 1) for index in range(spec.feature_width))
    beta = tuple(summary * 0.01 * (spec.feature_width - index) for index in range(spec.feature_width))
    down = tuple(
        tuple((((rank + 1) * (channel + 2)) % 5 - 2) / 10.0 for channel in range(spec.feature_width))
        for rank in range(spec.lora_rank)
    )
    up = tuple(
        tuple(summary * (channel + 1) * (rank + 1) / 100.0 for rank in range(spec.lora_rank))
        for channel in range(spec.feature_width)
    )
    return gamma, beta, down, up


def run_adapter_arm(
    spec: AdapterFixtureSpec,
    sequence_features: Sequence[Sequence[float]],
    target: SyntheticContextRow,
    *,
    arm: str,
    permuted_source: SyntheticContextRow | None = None,
) -> AdapterFixtureOutput:
    """Run a same-parent synthetic arm without weights, outcomes, or ATAC input."""

    spec.validate()
    if arm not in REQUIRED_ARMS:
        raise ContextBorzoiFixtureError(f"unregistered fixture arm: {arm}")
    target.validate(spec.feature_width)
    matrix = _validate_features(sequence_features, spec.feature_width)
    if arm == "sequence_only":
        output = apply_film_lora_precondition(
            matrix,
            gamma_delta=(0.0,) * spec.feature_width,
            beta=(0.0,) * spec.feature_width,
            lora_a=((0.0,) * spec.feature_width,) * spec.lora_rank,
            lora_b=((0.0,) * spec.lora_rank,) * spec.feature_width,
            alpha=spec.lora_alpha,
            context_observed=False,
            gate=0.0,
        )
        return AdapterFixtureOutput(
            OUTPUT_SCHEMA_VERSION,
            arm,
            arm,
            target.row_id,
            None,
            None,
            None,
            output,
            True,
            False,
            False,
            False,
            False,
            False,
        )

    context_source = target
    if arm == "permuted_context":
        if (
            permuted_source is None
            or permuted_source.row_id == target.row_id
            or permuted_source.permutation_stratum() != target.permutation_stratum()
        ):
            raise ContextBorzoiFixtureError(
                "permuted context must be a deranged source from the identical stratum"
            )
        permuted_source.validate(spec.feature_width)
        context_source = permuted_source
    elif permuted_source is not None:
        raise ContextBorzoiFixtureError("permuted source supplied to a nonpermuted arm")

    if not context_source.rna.complete:
        if arm != "conditioned":
            raise ContextBorzoiFixtureError(
                "incomplete context is valid only through conditioned-to-sequence routing"
            )
        sequence_only = run_adapter_arm(spec, matrix, target, arm="sequence_only")
        return AdapterFixtureOutput(
            OUTPUT_SCHEMA_VERSION,
            arm,
            "sequence_only_missing_context",
            target.row_id,
            target.row_id,
            target.rna.observed_mask,
            target.rna.states,
            sequence_only.values,
            True,
            False,
            False,
            False,
            False,
            False,
        )

    gamma, beta, down, up = _synthetic_adapter_terms(context_source.rna, spec)
    input_matrix = matrix
    sequence_consumed = True
    if arm == "trans_only":
        input_matrix = tuple((1.0,) * spec.feature_width for _row in matrix)
        sequence_consumed = False
    conditioned = apply_film_lora_precondition(
        input_matrix,
        gamma_delta=gamma,
        beta=beta,
        lora_a=down,
        lora_b=up,
        alpha=spec.lora_alpha,
        context_observed=True,
        gate=1.0,
    )
    if arm == "trans_only":
        conditioned = tuple(
            tuple(value - 1.0 for value in row) for row in conditioned
        )
    return AdapterFixtureOutput(
        OUTPUT_SCHEMA_VERSION,
        arm,
        arm,
        target.row_id,
        context_source.row_id,
        context_source.rna.observed_mask,
        context_source.rna.states,
        conditioned,
        sequence_consumed,
        True,
        False,
        False,
        False,
        False,
    )


def export_fixture_bundle(
    output: Path,
    spec: AdapterFixtureSpec,
    outputs: Mapping[str, AdapterFixtureOutput],
    assignments: Mapping[str, SyntheticContextRow],
) -> dict[str, object]:
    """Export canonical JSON only; no state dict, checkpoint, or pickle is written."""

    spec.validate()
    if set(outputs) != set(REQUIRED_ARMS) | {"missing_context"}:
        raise ContextBorzoiFixtureError("fixture output arm census differs")
    output.mkdir(parents=True, exist_ok=False)
    contract = {
        "schema_version": FIXTURE_SCHEMA_VERSION,
        "fixture_spec": asdict(spec),
        "required_arms": list(REQUIRED_ARMS),
        "remaining_parent_bound_ablations": ["reverse_complement", "matched_shuffle"],
        "conditional_model_authorized": False,
        "parent_component_selected": False,
        "architecture_built": False,
        "model_weights_present": False,
        "checkpoint_opened": False,
        "biological_data_read": False,
        "observed_atac_read": False,
        "observed_multiome_task_modified": False,
        "gpu_executed": False,
        "export_format": "canonical_json_utf8",
        "pickle_allowed": False,
    }
    payload = {
        "schema_version": FIXTURE_SCHEMA_VERSION,
        "outputs": {name: value.to_dict() for name, value in sorted(outputs.items())},
        "permuted_context_assignment": {
            key: value.row_id for key, value in sorted(assignments.items())
        },
    }
    (output / "fixture_contract.json").write_bytes(_canonical_json(contract))
    (output / "fixture_outputs.json").write_bytes(_canonical_json(payload))
    files = []
    for path in sorted(output.iterdir()):
        files.append(
            {
                "path": path.name,
                "sha256": _sha256(path),
                "size_bytes": path.stat().st_size,
                "format": "canonical_json_utf8",
            }
        )
    manifest = {
        "schema_version": "masld-bench-nonpickle-export-manifest-v1",
        "files": files,
        "pickle_present": False,
        "model_weights_present": False,
        "tensor_state_present": False,
    }
    (output / "export_manifest.json").write_bytes(_canonical_json(manifest))
    (output / "COMPLETE").write_text("complete\n", encoding="utf-8")
    verify_nonpickle_export(output)
    return manifest


def verify_nonpickle_export(output: Path) -> None:
    """Reject pickle-like extensions, magic, or an incomplete JSON manifest."""

    expected_files = {
        "COMPLETE",
        "export_manifest.json",
        "fixture_contract.json",
        "fixture_outputs.json",
    }
    observed_files = {path.name for path in output.iterdir() if path.is_file()}
    if observed_files != expected_files:
        raise ContextBorzoiFixtureError("export file census differs")
    for path in output.iterdir():
        if path.is_symlink() or not path.is_file():
            raise ContextBorzoiFixtureError(f"invalid export member: {path.name}")
        if path.suffix.lower() in PROHIBITED_EXPORT_SUFFIXES:
            raise ContextBorzoiFixtureError(f"prohibited export suffix: {path.name}")
        if path.read_bytes()[:1] == b"\x80":
            raise ContextBorzoiFixtureError(f"pickle protocol magic detected: {path.name}")
    manifest_path = output / "export_manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContextBorzoiFixtureError("export manifest is absent or invalid") from error
    if (
        manifest.get("pickle_present") is not False
        or manifest.get("model_weights_present") is not False
        or manifest.get("tensor_state_present") is not False
    ):
        raise ContextBorzoiFixtureError("export contains pickle or model state")
    for row in manifest.get("files", []):
        path = output / str(row["path"])
        if _sha256(path) != row.get("sha256"):
            raise ContextBorzoiFixtureError(f"export hash differs: {path.name}")
        if row.get("format") != "canonical_json_utf8":
            raise ContextBorzoiFixtureError(f"export format differs: {path.name}")
        json.loads(path.read_text(encoding="utf-8"))
    if (output / "COMPLETE").read_text(encoding="utf-8") != "complete\n":
        raise ContextBorzoiFixtureError("export completion marker differs")


__all__ = [
    "AdapterFixtureOutput",
    "AdapterFixtureSpec",
    "ContextBorzoiFixtureError",
    "FIXTURE_SCHEMA_VERSION",
    "OUTPUT_SCHEMA_VERSION",
    "REQUIRED_ARMS",
    "RNAContextMask",
    "SyntheticContextRow",
    "derange_context_assignments",
    "export_fixture_bundle",
    "run_adapter_arm",
    "verify_nonpickle_export",
]
