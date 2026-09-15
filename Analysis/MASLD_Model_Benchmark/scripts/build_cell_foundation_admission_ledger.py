#!/usr/bin/env python3
"""Build a read-only, development-only cell-foundation inclusion ledger."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


class AdmissionLedgerError(ValueError):
    """Raised when an included input or ledger request differs."""


CHECKPOINTS = {
    "geneformer_v1_10m": ("executions/geneformer_v1_10m-safe-inventory-21066210/checkpoint_inventory.json", "a5e33a757431643b3697de7ef6127950cdc49e06e58d4266b3a3ab191b683f14", 41183536),
    "geneformer_v2_104m": ("executions/geneformer_v2_104m-safe-inventory-21066212/checkpoint_inventory.json", "fff5cba29ddd8792991fa77b4872246fbe548a178cebda3775cdc72b67780e7f", 417571156),
    "geneformer_v2_316m": ("executions/geneformer_v2_316m-safe-inventory-21066211/checkpoint_inventory.json", "965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3", 1265455076),
    "scgpt_continual": ("executions/scgpt_continual-weights-only-inspection-21066249/weights_only_inspection.json", "ad0252a1971e0cd619b7116dbab3177432236c4537225d54280a2aa7e5fe402a", 207861754),
    "scgpt_whole_human": ("executions/scgpt_whole_human-weights-only-inspection-21066250/weights_only_inspection.json", "6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922", 205385258),
    "uce_4l": ("executions/uce_4l-weights-only-inspection-21066251/weights_only_inspection.json", "acb28f3f0a1d803e4a4ffe891b9bab38bf93c84762dc06b2452f0d515da91560", 3403514339),
    "uce_33l": ("executions/uce_33l-weights-only-inspection-21066252/weights_only_inspection.json", "3f458726196308e171611ed28394b55865708749f82af68cfd7771d5dfee661e", 5686283829),
    "scimilarity_v1_1": ("executions/scimilarity_v1_1_encoder-weights-only-inspection-21066254/weights_only_inspection.json", "a08f023788bdbded191e2d9344e5d4240ffc7d7505eca4e0b7a37526b858d19f", 124616535),
    "transcriptformer_sapiens": ("executions/transcriptformer_sapiens-weights-only-inspection-21066248/weights_only_inspection.json", "eff027e7393f92a32cadcbb2e13a917e24b62b71341571bf40f38ff6ec6d3021", 1719542526),
    "regformer": ("executions/regformer-weights-only-inspection-21066098/weights_only_inspection.json", "8e8e751a0b05caa4e6e4d6b5c0f8ee847ae607757f0431b4bee90cd9ef26a434", 197677201),
    "scprint_v1_5_medium": ("executions/scprint_v1_5_medium-weights-only-inspection-21066255/weights_only_inspection.json", "a4cf0753270d4ff451a5dbddadd4c59a53e10e4c3e96a8af0db407ad893c36c5", 221775592),
    "scprint2_small_v2": ("executions/scprint2_small_v2-weights-only-inspection-21066099/weights_only_inspection.json", "2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822", 889706878),
}


TERMS = {
    "geneformer": ("Apache-2.0", "Apache-2.0"),
    "scgpt": ("MIT", "UNDECLARED"),
    "uce": ("MIT", "CC-BY-4.0"),
    "scimilarity": ("Apache-2.0", "CC-BY-SA-4.0"),
    "transcriptformer": ("MIT", "MIT"),
    "regformer": ("MIT", "CC-BY-4.0"),
    "scprint": ("MIT", "Apache-2.0"),
    "scprint2": ("GPL-3.0-or-later", "Apache-2.0"),
    "cellplm": ("BSD-2-Clause", "UNDECLARED"),
    "scbert": ("GPL-3.0-or-later", "UNDECLARED"),
    "sccello": ("NO_LICENSE_DETECTED", "NO_LICENSE_DETECTED"),
    "langcell": ("MIT", "UNDECLARED"),
}


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _checkpoint(root: Path, model_id: str) -> dict[str, Any]:
    relative, expected_sha256, expected_size = CHECKPOINTS[model_id]
    path = root / relative
    payload = json.loads(path.read_text(encoding="utf-8"))
    if (
        payload.get("status") != "pass"
        or payload.get("sha256") != expected_sha256
        or payload.get("size_bytes") != expected_size
    ):
        raise AdmissionLedgerError(f"checkpoint inspection differs: {model_id}")
    return {
        "status": "exact_bytes_safe_inventory_passed",
        "sha256": expected_sha256,
        "size_bytes": expected_size,
        "tensor_count": payload.get("tensor_count"),
        "inspection_path": relative,
        "inspection_sha256": _sha256_file(path),
    }


def _family(model_id: str) -> str:
    for family in TERMS:
        if model_id == family or model_id.startswith(family + "_"):
            return family
    raise AdmissionLedgerError(f"unknown family: {model_id}")


def _admitted_record(root: Path, model_id: str) -> dict[str, Any]:
    family = _family(model_id)
    code_terms, weight_terms = TERMS[family]
    runtime = "runtime_and_forward_parity_pending"
    exposure = "unknown"
    head = "not_run"
    blockers: list[str] = []
    blockers.append("immutable_runtime_and_reference_forward_parity_pending")
    if family == "geneformer":
        exposure = "target_label_unexposed_sealed; development_source_overlap_is_source_specific"
        head = "existing_development_only_outputs_not_promoted; symmetric_head_review_pending"
        blockers.append("formal_identical_common_head_receipts_pending")
    elif family == "scgpt":
        exposure = "target_label_unexposed_sealed"
        blockers.append("checkpoint_weight_terms_undeclared")
        if model_id == "scgpt_continual":
            head = "existing_development_only_outputs_not_promoted; symmetric_head_review_pending"
            blockers.append("formal_identical_common_head_receipts_pending")
    elif family == "uce":
        exposure = "target_label_unexposed_sealed"
        head = "existing_development_only_outputs_not_promoted; symmetric_head_review_pending"
        blockers.append("formal_identical_common_head_receipts_pending")
    elif family == "scimilarity":
        exposure = "target_label_unexposed_sealed"
        head = "existing_development_only_outputs_not_promoted; symmetric_head_review_pending"
        blockers.append("formal_identical_common_head_receipts_pending")
    elif family in {"transcriptformer", "regformer"}:
        blockers.append("training_corpus_exposure_unknown")
    elif family == "scprint":
        exposure = "target_label_unexposed_sealed"
        blockers.extend(("exact_native_runtime_pending", "scprint_1_6_4_strict_restore_expects_gene_encoder_embeddings_weight_but_checkpoint_carries_gene_encoder_embedding_weight"))
    elif family == "scprint2":
        blockers.extend(("training_corpus_exposure_unknown", "exact_torch_2_8_scprint2_1_0_3_runtime_pending", "knn_metacell_fold_isolation_pending", "derivative_weight_terms_unresolved"))
    return {
        "model_id": model_id,
        "family": family,
        "checkpoint_gate": _checkpoint(root, model_id),
        "terms_gate": {
            "code": code_terms,
            "weights": weight_terms,
            "status": "pass" if "UNDECLARED" not in weight_terms and "NO_LICENSE" not in weight_terms else "blocked",
        },
        "exposure_gate": exposure,
        "runtime_gate": runtime,
        "development_common_head_gate": head,
        "blockers": blockers,
        "promotion_status": "not_promoted",
    }


def _blocked_record(model_id: str, family: str, blockers: list[str]) -> dict[str, Any]:
    code_terms, weight_terms = TERMS[family]
    return {
        "model_id": model_id,
        "family": family,
        "checkpoint_gate": {"status": "not_acquired"},
        "terms_gate": {"code": code_terms, "weights": weight_terms, "status": "blocked"},
        "exposure_gate": "not_sufficient_for_admission",
        "runtime_gate": "not_started_due_to_prior_gate",
        "development_common_head_gate": "not_run",
        "blockers": blockers,
        "promotion_status": "not_promoted",
    }


def build(root: Path, output: Path) -> dict[str, Any]:
    if output.exists() or not root.is_absolute():
        raise AdmissionLedgerError("ledger request differs")
    model_ids = tuple(CHECKPOINTS)
    records = [_admitted_record(root, model_id) for model_id in model_ids]
    records.extend(
        (
            _blocked_record("cellplm_20230926_85m", "cellplm", ["checkpoint_weight_terms_undeclared", "canonical_checkpoint_identity_unresolved"]),
            _blocked_record("scbert_panglao_pretrained", "scbert", ["checkpoint_weight_terms_undeclared", "authenticated_share_has_no_immutable_digest"]),
            _blocked_record("sccello_zeroshot", "sccello", ["no_detected_code_or_weight_license"]),
            _blocked_record("langcell_annotation_zeroshot", "langcell", ["checkpoint_weight_terms_undeclared", "required_multimodule_checkpoint_hashes_unresolved"]),
        )
    )
    config_sources = []
    for family in TERMS:
        for name in ("checkpoints.json", "exposure_audit.json", "development_crosswalk.json"):
            path = root / "config" / "artifacts" / "models" / family / name
            config_sources.append({"path": path.relative_to(root).as_posix(), "sha256": _sha256_file(path)})
    result = {
        "schema_version": "masld-bench-cell-foundation-admission-ledger-v1",
        "scope": "development_only_cell_state_methodology_candidate",
        "sealed_or_test_outcomes_opened": False,
        "champion_claims_allowed": False,
        "common_lane_rule": "identical frozen encoder representation plus identical donor-grouped nested common-head selection; not yet promoted",
        "native_lane_rule": "separate identity; no KNN/metacell information may cross donor, study, or outer-fold boundaries",
        "records": records,
        "source_config_artifacts": config_sources,
    }
    output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = build(arguments.root, arguments.output)
    print(json.dumps({"records": len(result["records"]), "sealed_or_test_outcomes_opened": False}))


if __name__ == "__main__":
    main()
