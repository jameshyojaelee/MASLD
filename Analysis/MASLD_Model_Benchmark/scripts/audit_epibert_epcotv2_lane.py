#!/usr/bin/env python3
"""Freeze fail-closed EpiBERT and EPCOTv2 lane readiness evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


class LaneAuditError(RuntimeError):
    """Raised when a frozen input or model requirement differs."""


EXPECTED = {
    "admission": "26e12915b05c2727c45eca091f85ec89f96c205e817c6ec9b64f181a563e7577",
    "runtime": "83e42afe350081391204703e9178ac30b52b7dfc15b863ba1af2e10611342aa4",
    "shape": "937e455dcf7bd96081e8bb148925c0fe9b3b286a6a97529c6d48fa31287a6c12",
    "forward": "cff3537ef5cf023db472acff523c791b758575db428bae5d8e9ec6d973eb7a5a",
    "tn5": "01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f",
    "motif": "aaa5f79457e991d470c69efd0066b9db25e395d4c4e8bd07bc4bffa3c2d6768c",
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise LaneAuditError(f"expected JSON object: {path}")
    return value


def require_artifact(root: Path, expected: str) -> dict[str, Any]:
    manifest = root / "ARTIFACTS.json"
    try:
        verified = verify_frozen_tree(root)
    except ArtifactError as error:
        raise LaneAuditError(f"frozen artifact differs: {root}: {error}") from error
    if digest(manifest) != expected:
        raise LaneAuditError(f"frozen artifact differs: {root}")
    return dict(verified)


def api_license_absent(metadata: dict[str, Any], expected_revision: str) -> bool:
    if metadata.get("sha") != expected_revision:
        raise LaneAuditError("Hugging Face revision differs")
    card = metadata.get("cardData") or {}
    if not isinstance(card, dict):
        raise LaneAuditError("Hugging Face cardData differs")
    siblings = metadata.get("siblings") or []
    names = {
        str(item.get("rfilename", "")).lower()
        for item in siblings
        if isinstance(item, dict)
    }
    return not card.get("license") and not any(
        name in {"license", "license.md", "license.txt", "copying"}
        for name in names
    )


def audit(args: argparse.Namespace) -> dict[str, Any]:
    admission = require_artifact(args.admission, EXPECTED["admission"])
    runtime = require_artifact(args.runtime, EXPECTED["runtime"])
    shape = require_artifact(args.shape_audit, EXPECTED["shape"])
    forward = require_artifact(args.forward_probe, EXPECTED["forward"])
    tn5 = require_artifact(args.tn5_bigwigs, EXPECTED["tn5"])
    receipt = load(args.admission / "admission_sources_receipt.json")
    extraction = load(args.runtime / "validation/selected-checkpoint-extraction.json")
    tn5_contract = load(args.tn5_bigwigs / "contract.json")
    forward_receipts = [
        load(args.forward_probe / f"models/{name}/receipt.json")
        for name in ("pretrained_model1", "pretrained_model2", "fine_tuned_rampage")
    ]
    space = load(args.epcot_space_api)
    reference = load(args.epcot_reference_api)

    motif = next(args.runtime.glob("source/EpiBERT-*/data_processing/motif_enrichment/consensus_pwms_vierstra.meme"))
    signal_wdl = next(args.runtime.glob("source/EpiBERT-*/data_processing/create_signal_tracks/fragments_to_bed_scores.wdl"))
    sea_wdl = next(args.runtime.glob("source/EpiBERT-*/data_processing/motif_enrichment/meme_run_sea.wdl"))
    signal_text = signal_wdl.read_text(encoding="utf-8")
    sea_text = sea_wdl.read_text(encoding="utf-8")

    epi = receipt["epibert"]
    if (
        epi["model_archive_sha256"] != "67684dc1f1562955a726f56fc7c53b758afe666f1138d833615fa1def46cf651"
        or not epi["safe_member_inventory_passed"]
        or extraction["status"] != "pass"
        or digest(motif) != EXPECTED["motif"]
        or "scale_factor_inv = 20.0 / scale_factor" not in signal_text
        or "$2-5,$3+5" not in signal_text
        or "memesuite/memesuite:5.4.1" not in sea_text
    ):
        raise LaneAuditError("EpiBERT source or checkpoint contract differs")
    fmeta = forward["metadata"]
    if not (
        fmeta["released_checkpoints_restored"]
        and fmeta["repeat_inference_bit_identical"]
        and not fmeta["project_data_read"]
        and not fmeta["outcomes_read"]
        and all(
            item["restore_existing_objects_matched"]
            and item["checkpoint_extra_training_state_expected"]
            and not item["checkpoint_resume_executed"]
            for item in forward_receipts
        )
    ):
        raise LaneAuditError("EpiBERT executable checkpoint evidence differs")
    if (
        shape["metadata"]["status"] != "passed"
        or runtime["metadata"]["status"] != "passed"
        or tn5["metadata"]["biological_unit"] != "donor"
        or tn5_contract["output_signal_contract"] != "two_insertions_per_unique_fragment_record"
        or not tn5_contract["source_coordinates_are_already_tn5_adjusted"]
    ):
        raise LaneAuditError("runtime or GSE296875 ATAC contract differs")

    epcot = receipt["epcotv2"]
    space_terms_absent = api_license_absent(
        space, "f93a2a944a49797a111f0294709416d38bedbcbe"
    )
    reference_terms_absent = api_license_absent(
        reference, "f445058917800811cb82c38b2c806674cc02b47f"
    )
    if not (
        epcot["checkpoint_lfs_sha256"] == "b824a80ed238e64e15c5a6eeae91329d24f04825c62d2a9dd328658822aba45d"
        and not epcot["checkpoint_object_downloaded"]
        and not epcot["checkpoint_deserialized"]
        and space_terms_absent
        and reference_terms_absent
    ):
        raise LaneAuditError("EPCOTv2 terms or checkpoint-access contract differs")

    return {
        "schema_version": "masld-bench-epibert-epcotv2-lane-audit-v2",
        "status": "pass",
        "input_artifacts_sha256": EXPECTED,
        "epibert": {
            "checkpoint_archive_sha256": epi["model_archive_sha256"],
            "released_checkpoint_prefixes": extraction["explicit_restore_prefixes"],
            "existing_model_objects_matched": True,
            "extra_checkpoint_training_state_permitted": True,
            "restore_consumed_assertion_passed": False,
            "repeat_inference_bit_identical": True,
            "native_inputs": ["hg38_sequence", "observed_ATAC", "693D_SEA_motif_context"],
            "gse296875_native_topology_supported": True,
            "gse296875_signal_transform_exact": False,
            "signal_transform_blocker": "Released WDL defines 20/scale_factor and a +/-5-bp expansion, but no invocation freezes scale_factor units for these donor-lineage tracks.",
            "motif_context_exact": False,
            "motif_context_blocker": "The 693-motif file and MEME 5.4.1 container are frozen, but top_peaks, half_peak_width, background_peaks and SEA threshold are not bound to the released checkpoints.",
            "gpu_execution_ready": False,
            "terminal_disposition": "checkpoint_runtime_ready_native_preprocessing_blocked",
        },
        "epcotv2": {
            "checkpoint_lfs_sha256": epcot["checkpoint_lfs_sha256"],
            "checkpoint_object_downloaded": False,
            "checkpoint_deserialized": False,
            "space_revision": space["sha"],
            "reference_revision": reference["sha"],
            "space_terms_absent": bool(space_terms_absent),
            "reference_terms_absent": bool(reference_terms_absent),
            "gpu_execution_ready": False,
            "terminal_disposition": "blocked_terms_no_execution",
        },
        "queue_disposition": {
            "gpu_wrapper_created": False,
            "gpu_queue_item_created": False,
            "reason": "Neither model currently passes all admission prerequisites.",
        },
        "native_task_contract": {
            "epibert": "masked observed-ATAC reconstruction, caQTL/accessibility deltas, and RAMPAGE proxy",
            "epcotv2": "sequence-plus-observed-ATAC multi-assay profile and contact prediction",
            "rna_conditioned_atac_models": [],
        },
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_labels_read": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--admission", type=Path, required=True)
    parser.add_argument("--runtime", type=Path, required=True)
    parser.add_argument("--shape-audit", type=Path, required=True)
    parser.add_argument("--forward-probe", type=Path, required=True)
    parser.add_argument("--tn5-bigwigs", type=Path, required=True)
    parser.add_argument("--epcot-space-api", type=Path, required=True)
    parser.add_argument("--epcot-reference-api", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise LaneAuditError("output already exists")
    result = audit(args)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "lane_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
