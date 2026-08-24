"""Frozen campaign authority for the single conditional MASLD architecture."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import json
from pathlib import Path
import shutil
import tempfile
from typing import Any

from .artifacts import (
    ArtifactError,
    freeze_tree,
    publish_directory_noreplace,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)
from .conditional_model import (
    ComplementarityEvidence,
    ConditionalModelError,
    build_context_model_spec,
)
from .hashing import canonical_sha256, sha256_file
from .tournament import (
    CELL_TASK,
    VARIANT_TASK,
    TournamentError,
    verify_conditional_model_decision,
    verify_development_shortlist,
    verify_selection_candidate_ledger,
)


class ContextSpecError(RuntimeError):
    """Raised when a novel-model campaign lacks frozen authorization."""


CONTEXT_MODEL_SPEC_ARTIFACT_SCHEMA_VERSION = (
    "masld-bench-context-model-spec-artifact-v1"
)
_SEQUENCE_FAMILY = "regulatory_sequence"
_CELL_FAMILY = "cell_foundation"


def _safe_root(path: str | Path, label: str) -> Path:
    try:
        root = reject_symlink_components(Path(path), label=label)
        if not root.is_absolute() or root.is_symlink() or not root.is_dir():
            raise ContextSpecError(f"{label} must be an absolute frozen directory")
        return root.resolve(strict=True)
    except ContextSpecError:
        raise
    except (ArtifactError, OSError) as exc:
        raise ContextSpecError(f"invalid {label}: {exc}") from exc


def _binding(
    root: Path, *, artifact_class: str, filename: str, identity: tuple[str, str]
) -> dict[str, str]:
    try:
        manifest = verify_frozen_tree(root)
    except (ArtifactError, OSError, ValueError, KeyError, TypeError) as exc:
        raise ContextSpecError(f"invalid frozen {artifact_class}: {exc}") from exc
    metadata = manifest.get("metadata")
    identity_field, identity_value = identity
    if not isinstance(metadata, Mapping) or metadata.get("artifact_class") != artifact_class:
        raise ContextSpecError(f"source is not classified as {artifact_class}")
    document = root / filename
    if document.is_symlink() or not document.is_file():
        raise ContextSpecError(f"frozen {artifact_class} lacks {filename}")
    return {
        "path": root.as_posix(),
        "manifest_sha256": sha256_file(root / "ARTIFACTS.json"),
        "document_sha256": sha256_file(document),
        identity_field: identity_value,
    }


def _within_one_se_ids(shortlist: Mapping[str, Any]) -> set[str]:
    candidates = [
        item
        for item in shortlist["candidates"]
        if isinstance(item, Mapping) and item.get("eligible") is True
    ]
    if not candidates:
        return set()
    metric = str(shortlist["primary_metric"])
    higher = bool(shortlist["higher_is_better"])

    def score(candidate: Mapping[str, Any]) -> float:
        return float(candidate["metrics"][metric])

    ordered = sorted(
        candidates,
        key=lambda item: (
            -score(item) if higher else score(item), str(item["candidate_id"])
        ),
    )
    best = ordered[0]
    error = float(best["standard_errors"][metric])
    threshold = score(best) - error if higher else score(best) + error
    return {
        str(item["candidate_id"])
        for item in candidates
        if (score(item) >= threshold if higher else score(item) <= threshold)
    }


def _candidate_pool(
    shortlist: Mapping[str, Any], *, family_id: str, require_within_one_se: bool
) -> list[dict[str, Any]]:
    ledger_binding = shortlist.get("selection_candidate_ledger_binding")
    if not isinstance(ledger_binding, Mapping):
        raise ContextSpecError("development shortlist lacks its scientific ledger")
    try:
        ledger = verify_selection_candidate_ledger(ledger_binding["path"])
    except (TournamentError, KeyError, TypeError) as exc:
        raise ContextSpecError(f"invalid shortlist scientific ledger: {exc}") from exc
    dispositions = {
        str(item["model_id"]): item
        for item in ledger["model_dispositions"]
        if isinstance(item, Mapping)
    }
    ledger_candidates = {
        str(item["candidate_id"]): item
        for item in ledger["model_candidates"]
        if isinstance(item, Mapping)
    }
    within_one_se = _within_one_se_ids(shortlist)
    selected = set(str(item) for item in shortlist["selected_candidate_ids"])
    pool = []
    for raw in shortlist["candidates"]:
        if not isinstance(raw, Mapping):
            continue
        candidate_id = str(raw["candidate_id"])
        if candidate_id not in selected or raw.get("family") != family_id:
            continue
        model_id = str(raw["model_id"])
        disposition = dispositions.get(model_id)
        if not isinstance(disposition, Mapping):
            raise ContextSpecError(f"shortlisted model {model_id} lacks a disposition")
        ledger_candidate = ledger_candidates.get(candidate_id)
        if not isinstance(ledger_candidate, Mapping):
            raise ContextSpecError(
                f"shortlisted candidate {candidate_id} lacks a scientific ledger record"
            )
        decision = disposition.get("open_champion_license_decision")
        decision_sha256 = disposition.get(
            "open_champion_license_decision_sha256"
        )
        if not isinstance(decision, Mapping) or not isinstance(
            decision_sha256, str
        ):
            raise ContextSpecError(f"shortlisted model {model_id} lacks a license decision")
        metric = str(shortlist["primary_metric"])
        pool.append(
            {
                "candidate_id": candidate_id,
                "model_id": model_id,
                "adaptation_regime": str(raw["adaptation_regime"]),
                "candidate_configuration_sha256": str(
                    raw["candidate_configuration_sha256"]
                ),
                "selection_artifact_ensemble_sha256s": dict(
                    ledger_candidate["selection_artifact_ensemble_sha256s"]
                ),
                "admitted": raw.get("eligible") is True
                and disposition.get("champion_eligible") is True,
                "within_one_se": candidate_id in within_one_se,
                "development_score": float(raw["metrics"][metric]),
                "compute_cost_status": raw.get("complexity_status"),
                "open_champion_eligible": disposition.get(
                    "open_champion_eligible"
                ),
                "open_champion_license_decision": dict(decision),
                "open_champion_license_decision_sha256": decision_sha256,
            }
        )
    if require_within_one_se:
        pool = [item for item in pool if item["within_one_se"] is True]
    return sorted(pool, key=lambda item: (item["model_id"], item["candidate_id"]))


def _source_state(
    *,
    conditional_decision_dir: str | Path,
    sequence_shortlist_dir: str | Path,
    cell_shortlist_dir: str | Path,
) -> dict[str, Any]:
    decision_root = _safe_root(conditional_decision_dir, "conditional decision")
    sequence_root = _safe_root(sequence_shortlist_dir, "sequence shortlist")
    cell_root = _safe_root(cell_shortlist_dir, "cell shortlist")
    try:
        decision = verify_conditional_model_decision(decision_root)
        sequence_shortlist = verify_development_shortlist(sequence_root)
        cell_shortlist = verify_development_shortlist(cell_root)
    except TournamentError as exc:
        raise ContextSpecError(f"invalid context-model source authority: {exc}") from exc
    if decision.get("triggered") is not True:
        raise ContextSpecError("conditional decision did not authorize a novel model")
    if sequence_shortlist.get("task_id") != VARIANT_TASK:
        raise ContextSpecError("sequence shortlist must be variant_to_regulation")
    if cell_shortlist.get("task_id") != CELL_TASK:
        raise ContextSpecError("cell shortlist must be cell_state_mapping")
    for label, shortlist in (
        ("sequence", sequence_shortlist),
        ("cell", cell_shortlist),
    ):
        if shortlist.get("source_wave") != "full_specialist_screen" or shortlist.get(
            "source_task_seeds"
        ) != [1103, 2909, 4721, 6673, 8111]:
            raise ContextSpecError(
                f"{label} shortlist must be the five-seed full specialist screen"
            )
    decision_sequence_bindings = [
        item
        for item in decision["development_shortlist_bindings"]
        if isinstance(item, Mapping) and item.get("shortlist_id") == sequence_shortlist["shortlist_id"]
    ]
    if len(decision_sequence_bindings) != 1:
        raise ContextSpecError(
            "sequence shortlist is not bound by the triggered conditional decision"
        )
    sequence_binding = _binding(
        sequence_root,
        artifact_class="development_shortlist",
        filename="development_shortlist.json",
        identity=("shortlist_id", str(sequence_shortlist["shortlist_id"])),
    )
    if dict(decision_sequence_bindings[0]) != sequence_binding:
        raise ContextSpecError("conditional decision sequence binding changed")
    evidence_binding = decision.get("complementarity_evidence", {}).get(
        "stacking_evidence_binding"
    )
    if not isinstance(evidence_binding, Mapping):
        raise ContextSpecError("conditional decision lacks stacking evidence")
    try:
        evidence = ComplementarityEvidence.from_stacking_evidence(
            evidence_binding["path"]
        )
    except (ConditionalModelError, KeyError, TypeError) as exc:
        raise ContextSpecError(f"invalid stacking evidence: {exc}") from exc
    sequence_candidates = _candidate_pool(
        sequence_shortlist,
        family_id=_SEQUENCE_FAMILY,
        require_within_one_se=True,
    )
    cell_candidates = _candidate_pool(
        cell_shortlist,
        family_id=_CELL_FAMILY,
        require_within_one_se=False,
    )
    try:
        model_spec = build_context_model_spec(
            evidence=evidence,
            sequence_candidates=sequence_candidates,
            cell_candidates=cell_candidates,
        )
    except ConditionalModelError as exc:
        raise ContextSpecError(f"context model cannot be specified: {exc}") from exc
    decision_binding = _binding(
        decision_root,
        artifact_class="conditional_decision",
        filename="conditional_decision.json",
        identity=("decision_id", str(decision["decision_id"])),
    )
    cell_binding = _binding(
        cell_root,
        artifact_class="development_shortlist",
        filename="development_shortlist.json",
        identity=("shortlist_id", str(cell_shortlist["shortlist_id"])),
    )
    return {
        "model_spec": model_spec,
        "conditional_decision_binding": decision_binding,
        "sequence_shortlist_binding": sequence_binding,
        "cell_shortlist_binding": cell_binding,
        "candidate_pool_sha256": canonical_sha256(
            {
                "sequence": sequence_candidates,
                "cell": cell_candidates,
            }
        ),
    }


def freeze_context_model_spec(
    *,
    conditional_decision_dir: str | Path,
    sequence_shortlist_dir: str | Path,
    cell_shortlist_dir: str | Path,
    output_root: str | Path,
) -> Path:
    """Freeze the only artifact that may authorize a novel-model campaign."""

    state = _source_state(
        conditional_decision_dir=conditional_decision_dir,
        sequence_shortlist_dir=sequence_shortlist_dir,
        cell_shortlist_dir=cell_shortlist_dir,
    )
    identity = {
        "schema_version": CONTEXT_MODEL_SPEC_ARTIFACT_SCHEMA_VERSION,
        **state,
        "architecture_count_cap": 1,
        "created_before_outcome_unblind": True,
        "sealed_results_used": False,
        "campaign_authority": "conditional_model_wave_only",
    }
    payload = {"context_model_spec_id": canonical_sha256(identity), **identity}
    output = reject_symlink_components(
        Path(output_root), label="context model spec output root"
    ).resolve()
    output.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".context-model-spec.", dir=output))
    try:
        write_json_exclusive(staging / "context_model_spec.json", payload, mode=0o440)
        freeze_tree(
            staging,
            {
                "artifact_class": "context_model_spec",
                "context_model_spec_id": payload["context_model_spec_id"],
            },
        )
        target = output / f"context-model-spec--{payload['context_model_spec_id']}"
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_context_model_spec(target)
    return target


def verify_context_model_spec(path: str | Path) -> dict[str, Any]:
    """Re-derive a frozen context architecture from all source authorities."""

    root = _safe_root(path, "context model spec")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "context_model_spec.json").read_text(encoding="utf-8")
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContextSpecError(f"invalid frozen context model spec: {exc}") from exc
    fields = {
        "context_model_spec_id",
        "schema_version",
        "model_spec",
        "conditional_decision_binding",
        "sequence_shortlist_binding",
        "cell_shortlist_binding",
        "candidate_pool_sha256",
        "architecture_count_cap",
        "created_before_outcome_unblind",
        "sealed_results_used",
        "campaign_authority",
    }
    if not isinstance(payload, Mapping) or set(payload) != fields:
        raise ContextSpecError("context model spec has the wrong exact schema")
    if payload.get("schema_version") != CONTEXT_MODEL_SPEC_ARTIFACT_SCHEMA_VERSION:
        raise ContextSpecError("unsupported context model spec schema")
    claimed = payload.get("context_model_spec_id")
    identity = dict(payload)
    identity.pop("context_model_spec_id", None)
    if not isinstance(claimed, str) or canonical_sha256(identity) != claimed:
        raise ContextSpecError("context model spec identity does not rederive")
    if manifest.get("metadata") != {
        "artifact_class": "context_model_spec",
        "context_model_spec_id": claimed,
    }:
        raise ContextSpecError("context model spec manifest metadata mismatch")
    for field in (
        "conditional_decision_binding",
        "sequence_shortlist_binding",
        "cell_shortlist_binding",
    ):
        if not isinstance(payload.get(field), Mapping) or not isinstance(
            payload[field].get("path"), str
        ):
            raise ContextSpecError(f"context model spec lacks {field}")
    state = _source_state(
        conditional_decision_dir=payload["conditional_decision_binding"]["path"],
        sequence_shortlist_dir=payload["sequence_shortlist_binding"]["path"],
        cell_shortlist_dir=payload["cell_shortlist_binding"]["path"],
    )
    expected = {
        "context_model_spec_id": claimed,
        "schema_version": CONTEXT_MODEL_SPEC_ARTIFACT_SCHEMA_VERSION,
        **state,
        "architecture_count_cap": 1,
        "created_before_outcome_unblind": True,
        "sealed_results_used": False,
        "campaign_authority": "conditional_model_wave_only",
    }
    if dict(payload) != expected:
        raise ContextSpecError("context model spec does not rederive from frozen sources")
    return dict(payload)


__all__ = [
    "CONTEXT_MODEL_SPEC_ARTIFACT_SCHEMA_VERSION",
    "ContextSpecError",
    "freeze_context_model_spec",
    "verify_context_model_spec",
]
