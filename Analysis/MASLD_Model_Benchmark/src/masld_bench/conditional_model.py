"""Fail-closed authorization and specification for the one conditional model."""

from __future__ import annotations

from dataclasses import dataclass
import json
from math import isfinite
from pathlib import Path
from typing import Any, Mapping, Sequence

from .artifacts import (
    ArtifactError,
    canonical_hash,
    reject_symlink_components,
    verify_frozen_tree,
)
from .hashing import sha256_file


class ConditionalModelError(RuntimeError):
    """Raised when novelty is attempted without the prespecified evidence trigger."""


STACKING_EVIDENCE_SCHEMA_VERSION = "masld-bench-stacking-evidence-v1"
STACKING_COMPUTE_BASIS_POLICY = (
    "tie_break_requires_one_hash_bound_measured_compute_basis_v1"
)
_STACKING_EVIDENCE_FILENAME = "stacking_evidence.json"
_STACKING_EVIDENCE_SOURCE_CONTRACTS = {
    "development_gain_bundle": {
        "count": 2,
        "artifact_class": "development_stack_gain_bundle",
        "document_filename": "development_stack_gain_bundle.json",
    },
    "development_residual_bundle": {
        "count": 2,
        "artifact_class": "development_residual_bundle",
        "document_filename": "development_residual_bundle.json",
    },
    "development_shortlist": {
        "count": 2,
        "artifact_class": "development_shortlist",
        "document_filename": "development_shortlist.json",
    },
}
_VERIFIED_STACKING_EVIDENCE = object()


def _verified_source_binding(raw: object, *, index: int) -> dict[str, str]:
    if not isinstance(raw, Mapping) or set(raw) != {
        "artifact_class",
        "document_filename",
        "document_sha256",
        "manifest_sha256",
        "path",
        "role",
    }:
        raise ConditionalModelError(
            f"stacking evidence source binding {index} has the wrong schema"
        )
    role = str(raw["role"])
    contract = _STACKING_EVIDENCE_SOURCE_CONTRACTS.get(role)
    if contract is None:
        raise ConditionalModelError(
            f"stacking evidence source binding {index} has unsupported role {role!r}"
        )
    if (
        raw["artifact_class"] != contract["artifact_class"]
        or raw["document_filename"] != contract["document_filename"]
    ):
        raise ConditionalModelError(
            f"stacking evidence source binding {index} differs from its role contract"
        )
    try:
        source = reject_symlink_components(
            Path(str(raw["path"])), label=f"stacking evidence source {index}"
        )
        if not source.is_absolute() or source.is_symlink() or not source.is_dir():
            raise ConditionalModelError(
                "stacking evidence sources must be absolute frozen directories"
            )
        source = source.resolve(strict=True)
        manifest = verify_frozen_tree(source)
    except ConditionalModelError:
        raise
    except (ArtifactError, OSError, ValueError) as exc:
        raise ConditionalModelError(
            f"stacking evidence source {index} is not a valid frozen tree: {exc}"
        ) from exc
    artifact_class = str(raw["artifact_class"])
    metadata = manifest.get("metadata")
    if (
        not isinstance(metadata, Mapping)
        or metadata.get("artifact_class") != artifact_class
    ):
        raise ConditionalModelError(
            f"stacking evidence source {index} changed artifact class"
        )
    document_filename = str(raw["document_filename"])
    document = source / document_filename
    if (
        Path(document_filename).name != document_filename
        or document.is_symlink()
        or not document.is_file()
    ):
        raise ConditionalModelError(
            f"stacking evidence source {index} has an invalid document"
        )
    normalized = {
        "role": role,
        "artifact_class": artifact_class,
        "path": source.as_posix(),
        "manifest_sha256": sha256_file(source / "ARTIFACTS.json"),
        "document_filename": document_filename,
        "document_sha256": sha256_file(document),
    }
    if dict(raw) != normalized:
        raise ConditionalModelError(
            f"stacking evidence source binding {index} changed"
        )
    return normalized


def _load_stacking_evidence(path: str | Path) -> dict[str, Any]:
    if not isinstance(path, (str, Path)):
        raise ConditionalModelError(
            "stacking evidence must be supplied as a frozen-tree path"
        )
    try:
        root = reject_symlink_components(Path(path), label="stacking evidence")
        if not root.is_absolute() or root.is_symlink() or not root.is_dir():
            raise ConditionalModelError(
                "stacking evidence must be an absolute frozen directory"
            )
        root = root.resolve(strict=True)
        manifest = verify_frozen_tree(root)
        document = root / _STACKING_EVIDENCE_FILENAME
        if document.is_symlink() or not document.is_file():
            raise ConditionalModelError(
                f"stacking evidence lacks {_STACKING_EVIDENCE_FILENAME}"
            )
        payload = json.loads(document.read_text(encoding="utf-8"))
    except ConditionalModelError:
        raise
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ConditionalModelError(f"invalid frozen stacking evidence: {exc}") from exc
    expected_fields = {
        "created_before_outcome_unblind",
        "derived_evidence",
        "evidence_id",
        "schema_version",
        "sealed_results_used",
        "source_bindings",
    }
    if not isinstance(payload, Mapping) or set(payload) != expected_fields:
        raise ConditionalModelError("stacking evidence has the wrong exact schema")
    if payload.get("schema_version") != STACKING_EVIDENCE_SCHEMA_VERSION:
        raise ConditionalModelError(
            f"stacking evidence requires {STACKING_EVIDENCE_SCHEMA_VERSION}"
        )
    if payload.get("created_before_outcome_unblind") is not True:
        raise ConditionalModelError(
            "stacking evidence was not frozen before outcome unblinding"
        )
    if payload.get("sealed_results_used") is not False:
        raise ConditionalModelError("stacking evidence may not use sealed results")
    raw_sources = payload.get("source_bindings")
    if not isinstance(raw_sources, list):
        raise ConditionalModelError("stacking evidence source_bindings must be an array")
    sources = [
        _verified_source_binding(raw, index=index)
        for index, raw in enumerate(raw_sources)
    ]
    observed_roles = {
        role: sum(source["role"] == role for source in sources)
        for role in _STACKING_EVIDENCE_SOURCE_CONTRACTS
    }
    expected_roles = {
        role: int(contract["count"])
        for role, contract in _STACKING_EVIDENCE_SOURCE_CONTRACTS.items()
    }
    if observed_roles != expected_roles:
        raise ConditionalModelError(
            "stacking evidence must bind exactly two shortlists, two residual "
            "bundles, and two development gain bundles"
        )
    if sources != sorted(
        sources, key=lambda item: (item["role"], item["path"])
    ):
        raise ConditionalModelError("stacking evidence source bindings must be sorted")
    derived = payload.get("derived_evidence")
    if not isinstance(derived, Mapping):
        raise ConditionalModelError("stacking evidence lacks derived evidence")
    identity = dict(payload)
    claimed = identity.pop("evidence_id", None)
    if not isinstance(claimed, str) or canonical_hash(identity) != claimed:
        raise ConditionalModelError("stacking evidence identity does not rederive")
    metadata = manifest.get("metadata")
    if metadata != {
        "artifact_class": "stacking_evidence",
        "evidence_id": claimed,
    }:
        raise ConditionalModelError("stacking evidence manifest metadata mismatch")
    return {
        "root": root.as_posix(),
        "manifest_sha256": sha256_file(root / "ARTIFACTS.json"),
        "document_sha256": sha256_file(document),
        "payload": dict(payload),
    }


@dataclass(frozen=True, slots=True, init=False)
class ComplementarityEvidence:
    """Verified trigger inputs loaded only from one recursively frozen tree."""

    residual_correlation: float
    relative_deviance_gains: Mapping[str, float]
    absolute_correlation_or_f1_gains: Mapping[str, float]
    qualifying_seed_count: int
    evaluated_seed_count: int
    study_count: int
    cross_fitted: bool
    nonnegative_stack: bool
    best_open_models_compared: bool
    _artifact_root: str
    _manifest_sha256: str
    _document_sha256: str
    _verification_token: object

    @classmethod
    def from_stacking_evidence(
        cls, path: str | Path
    ) -> "ComplementarityEvidence":
        """Build evidence only after recursive verification of a frozen tree."""

        verified = _load_stacking_evidence(path)
        derived = verified["payload"]["derived_evidence"]
        required = {
            "absolute_correlation_or_f1_gains",
            "best_open_models_compared",
            "cross_fitted",
            "evaluated_seed_count",
            "nonnegative_stack",
            "qualifying_seed_count",
            "relative_deviance_gains",
            "residual_correlation",
            "study_count",
        }
        if set(derived) != required:
            raise ConditionalModelError(
                "stacking derived evidence has the wrong exact schema"
            )
        instance = object.__new__(cls)
        values = {
            "residual_correlation": derived["residual_correlation"],
            "relative_deviance_gains": dict(derived["relative_deviance_gains"]),
            "absolute_correlation_or_f1_gains": dict(
                derived["absolute_correlation_or_f1_gains"]
            ),
            "qualifying_seed_count": derived["qualifying_seed_count"],
            "evaluated_seed_count": derived["evaluated_seed_count"],
            "study_count": derived["study_count"],
            "cross_fitted": derived["cross_fitted"],
            "nonnegative_stack": derived["nonnegative_stack"],
            "best_open_models_compared": derived["best_open_models_compared"],
            "_artifact_root": verified["root"],
            "_manifest_sha256": verified["manifest_sha256"],
            "_document_sha256": verified["document_sha256"],
            "_verification_token": _VERIFIED_STACKING_EVIDENCE,
        }
        for field_name, value in values.items():
            object.__setattr__(instance, field_name, value)
        return instance

    @property
    def stacking_evidence_binding(self) -> dict[str, str]:
        return {
            "path": self._artifact_root,
            "manifest_sha256": self._manifest_sha256,
            "document_sha256": self._document_sha256,
        }

    @property
    def derivation_source(self) -> str:
        return STACKING_EVIDENCE_SCHEMA_VERSION

    @property
    def qualifying_endpoints(self) -> tuple[str, ...]:
        relative = {
            endpoint
            for endpoint, gain in self.relative_deviance_gains.items()
            if float(gain) >= 0.05
        }
        absolute = {
            endpoint
            for endpoint, gain in self.absolute_correlation_or_f1_gains.items()
            if float(gain) >= 0.02
        }
        return tuple(sorted(relative | absolute))

    def require_verified_artifact(self) -> None:
        """Re-verify the evidence tree without requiring the trigger to pass."""

        failures = []
        if self._verification_token is not _VERIFIED_STACKING_EVIDENCE:
            failures.append(
                "complementarity evidence was not derived from a frozen "
                "stacking-evidence artifact"
            )
        else:
            try:
                verified = _load_stacking_evidence(self._artifact_root)
                rebound = {
                    "path": verified["root"],
                    "manifest_sha256": verified["manifest_sha256"],
                    "document_sha256": verified["document_sha256"],
                }
                if rebound != self.stacking_evidence_binding:
                    failures.append("stacking-evidence artifact binding changed")
                derived = verified["payload"]["derived_evidence"]
                observed = {
                    "residual_correlation": self.residual_correlation,
                    "relative_deviance_gains": dict(
                        self.relative_deviance_gains
                    ),
                    "absolute_correlation_or_f1_gains": dict(
                        self.absolute_correlation_or_f1_gains
                    ),
                    "qualifying_seed_count": self.qualifying_seed_count,
                    "evaluated_seed_count": self.evaluated_seed_count,
                    "study_count": self.study_count,
                    "cross_fitted": self.cross_fitted,
                    "nonnegative_stack": self.nonnegative_stack,
                    "best_open_models_compared": self.best_open_models_compared,
                }
                if dict(derived) != observed:
                    failures.append("stacking derived evidence changed")
            except ConditionalModelError as exc:
                failures.append(str(exc))
        if failures:
            raise ConditionalModelError(
                "invalid stacking-evidence authorization: " + "; ".join(failures)
            )

    def require_trigger(self) -> None:
        self.require_verified_artifact()
        failures = []
        try:
            residual_correlation = float(self.residual_correlation)
        except (TypeError, ValueError):
            residual_correlation = float("nan")
        if (
            isinstance(self.residual_correlation, bool)
            or not isfinite(residual_correlation)
            or not -1.0 <= residual_correlation <= 1.0
        ):
            failures.append("held residual correlation is not finite in [-1, 1]")
        elif residual_correlation >= 0.80:
            failures.append("held residual correlation is not below 0.80")
        relative_gains = self.relative_deviance_gains
        absolute_gains = self.absolute_correlation_or_f1_gains
        gains_are_mappings = isinstance(relative_gains, Mapping) and isinstance(
            absolute_gains, Mapping
        )
        overlap = set(relative_gains) & set(absolute_gains) if gains_are_mappings else set()
        if overlap:
            failures.append(
                "endpoint gain scale is ambiguous for " + ", ".join(sorted(overlap))
            )
        gains = {**relative_gains, **absolute_gains} if gains_are_mappings else {}
        valid_gains = bool(gains)
        if valid_gains:
            for endpoint, gain in gains.items():
                try:
                    numeric_gain = float(gain)
                except (TypeError, ValueError):
                    valid_gains = False
                    break
                if (
                    not str(endpoint).strip()
                    or isinstance(gain, bool)
                    or not isfinite(numeric_gain)
                ):
                    valid_gains = False
                    break
        if not valid_gains:
            failures.append("endpoint gains must be non-empty finite numbers")
        qualifying_endpoints = self.qualifying_endpoints if valid_gains else ()
        if len(qualifying_endpoints) < 2:
            failures.append("fewer than two mechanistic endpoints meet the gain threshold")
        valid_seed_counts = all(
            isinstance(value, int) and not isinstance(value, bool)
            for value in (self.qualifying_seed_count, self.evaluated_seed_count)
        )
        if not valid_seed_counts or self.evaluated_seed_count != 5:
            failures.append("the trigger must evaluate exactly five frozen seeds")
        if (
            not valid_seed_counts
            or self.qualifying_seed_count < 4
            or self.qualifying_seed_count > 5
        ):
            failures.append("gain does not persist in at least four of five seeds")
        if (
            not isinstance(self.study_count, int)
            or isinstance(self.study_count, bool)
            or self.study_count < 2
        ):
            failures.append("gain does not span more than one study")
        if self.cross_fitted is not True:
            failures.append("stacking evidence is not cross-fitted")
        if self.nonnegative_stack is not True:
            failures.append("stacking weights are not constrained nonnegative")
        if self.best_open_models_compared is not True:
            failures.append("residuals are not from the best open sequence and context models")
        if failures:
            raise ConditionalModelError("conditional model is not authorized: " + "; ".join(failures))


def _license_bound_candidate(
    candidate: Mapping[str, Any], *, require_clean_declared: bool
) -> bool:
    """Require an admitted candidate to carry its immutable open-license decision."""

    decision = candidate.get("open_champion_license_decision")
    decision_sha256 = candidate.get("open_champion_license_decision_sha256")
    if not isinstance(decision, Mapping) or not isinstance(decision_sha256, str):
        return False
    if canonical_hash(decision) != decision_sha256:
        return False
    if decision.get("eligible") is not True:
        return False
    if candidate.get("open_champion_eligible") is not True:
        return False
    if candidate.get("admitted") is not True:
        return False
    if require_clean_declared and decision.get("exposure_status") != "clean_declared":
        return False
    return True


def _rank_one(
    candidates: Sequence[Mapping[str, Any]], label: str
) -> Mapping[str, Any]:
    """Rank on the development endpoint; consult compute only on an exact tie.

    Compute was previously the second sort key with a ``float("inf")`` default,
    so an absent or incomparable cost silently ordered candidates.  Every
    producer in this tree records
    ``compute_cost_status = "unavailable_no_comparable_measured_basis"``, so a
    genuine tie now fails closed instead of being broken by an unmeasured
    number.
    """

    ranked = sorted(
        candidates,
        key=lambda item: (
            -float(item.get("development_score", float("-inf"))),
            str(item["model_id"]),
        ),
    )
    best = ranked[0]
    tied = [
        item
        for item in ranked
        if float(item.get("development_score", float("-inf")))
        == float(best.get("development_score", float("-inf")))
    ]
    if len(tied) == 1:
        return best
    bases = {
        (
            item.get("compute_cost_status"),
            item.get("compute_basis_sha256"),
        )
        for item in tied
    }
    if len(bases) != 1 or next(iter(bases))[0] != "measured":
        raise ConditionalModelError(
            f"{label} selection is tied without one hash-bound comparable "
            f"compute basis ({STACKING_COMPUTE_BASIS_POLICY})"
        )
    return min(
        tied,
        key=lambda item: (
            float(item["compute_cost"]),
            str(item["model_id"]),
        ),
    )


def _selected_candidate_binding(
    candidate: Mapping[str, Any], label: str
) -> dict[str, Any]:
    required_strings = (
        "candidate_id",
        "adaptation_regime",
        "candidate_configuration_sha256",
    )
    values: dict[str, str] = {}
    for field in required_strings:
        value = candidate.get(field)
        if not isinstance(value, str) or not value:
            raise ConditionalModelError(f"{label} lacks frozen {field}")
        if field.endswith("sha256") or field == "candidate_id":
            if len(value) != 64 or any(
                character not in "0123456789abcdef" for character in value
            ):
                raise ConditionalModelError(
                    f"{label}.{field} must be a lowercase SHA-256 identifier"
                )
        values[field] = value
    fit_state = candidate.get("selection_artifact_ensemble_sha256s")
    expected_fit_fields = {
        "checkpoint_sha256",
        "preprocessing_sha256",
        "task_head_sha256",
        "calibration_sha256",
    }
    if not isinstance(fit_state, Mapping) or set(fit_state) != expected_fit_fields:
        raise ConditionalModelError(f"{label} lacks exact fitted-state bindings")
    normalized_fit = {}
    for field in sorted(expected_fit_fields):
        value = fit_state[field]
        if not isinstance(value, str) or len(value) != 64 or any(
            character not in "0123456789abcdef" for character in value
        ):
            raise ConditionalModelError(
                f"{label}.{field} must be a lowercase SHA-256 identifier"
            )
        normalized_fit[field] = value
    return {**values, "selection_artifact_ensemble_sha256s": normalized_fit}


def build_context_model_spec(
    *,
    evidence: ComplementarityEvidence,
    sequence_candidates: Sequence[Mapping[str, Any]],
    cell_candidates: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Select lower-compute eligible backbones within one SE, then freeze design."""

    evidence.require_trigger()
    eligible_sequence = [
        candidate
        for candidate in sequence_candidates
        if candidate.get("within_one_se") is True
        and _license_bound_candidate(candidate, require_clean_declared=False)
    ]
    eligible_cell = [
        candidate
        for candidate in cell_candidates
        if _license_bound_candidate(candidate, require_clean_declared=True)
    ]
    if not eligible_sequence:
        raise ConditionalModelError("no open long-range sequence backbone is within one SE")
    if not eligible_cell:
        raise ConditionalModelError("no clean open cell-state encoder is admitted")
    sequence = _rank_one(eligible_sequence, "open long-range sequence backbone")
    cell = _rank_one(eligible_cell, "open cell-state encoder")
    sequence_binding = _selected_candidate_binding(
        sequence, "open long-range sequence backbone"
    )
    cell_binding = _selected_candidate_binding(cell, "open cell-state encoder")
    identity = {
        "schema_version": "masld-context-model-spec-v1",
        "authorization": {
            "residual_correlation": evidence.residual_correlation,
            "relative_deviance_gains": dict(
                sorted(evidence.relative_deviance_gains.items())
            ),
            "absolute_correlation_or_f1_gains": dict(
                sorted(evidence.absolute_correlation_or_f1_gains.items())
            ),
            "qualifying_endpoints": list(evidence.qualifying_endpoints),
            "qualifying_seed_count": evidence.qualifying_seed_count,
            "evaluated_seed_count": evidence.evaluated_seed_count,
            "study_count": evidence.study_count,
            "cross_fitted": evidence.cross_fitted,
            "nonnegative_stack": evidence.nonnegative_stack,
            "best_open_models_compared": evidence.best_open_models_compared,
            "derivation_source": evidence.derivation_source,
            "stacking_evidence_binding": evidence.stacking_evidence_binding,
            "compute_tie_break_policy": STACKING_COMPUTE_BASIS_POLICY,
        },
        "sequence_backbone": sequence["model_id"],
        "sequence_candidate_binding": sequence_binding,
        "sequence_license_decision_sha256": sequence[
            "open_champion_license_decision_sha256"
        ],
        "cell_encoder": cell["model_id"],
        "cell_candidate_binding": cell_binding,
        "cell_license_decision_sha256": cell[
            "open_champion_license_decision_sha256"
        ],
        "competing_context_path": "training_fold_hvg_pseudobulk_encoder",
        "conditioning": ["gated_film", "lora_adapters"],
        "heads": ["atac_profile", "atac_count", "rna_coverage"],
        "forbidden_heads": ["masld_histone"],
        "missing_modality_masks": True,
        "outputs": ["gene_score", "peak_score", "variant_score", "five_seed_uncertainty"],
        "required_ablations": [
            "sequence_only",
            "trans_only",
            "permuted_context",
            "reverse_complement",
            "matched_shuffle",
        ],
        "hotspot_117_usage": "read_only_post_hoc_projection_only",
        "retrospective_diagnostics": ["mpra", "crop_seq"],
        "architecture_count_cap": 1,
    }
    identity["spec_sha256"] = canonical_hash(identity)
    return identity
