#!/usr/bin/env python3
"""Reconcile Enformer and Sei sequence-native lane evidence fail closed."""

from __future__ import annotations

import argparse
import hashlib
import json
import tomllib
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


class LaneAuditError(RuntimeError):
    """Raised when a frozen output file or lane requirement differs."""


AUTHORITY_SHA256 = {
    "config/artifacts/models/enformer/checkpoints.json": "114b724c94908081f5e22f46731b3192a16947d9621497cd155e6183b8caff6e",
    "config/artifacts/models/enformer/development_crosswalk.json": "64fe8bedf4d7e6a9caeb286c2a67ac60c0bfbfc2c5fd597ae5582bee04575bf1",
    "config/artifacts/models/enformer/exposure_audit.json": "a2fefdaa1b076635fabb6bee0cca9528405fcb900ad4e1e1548c3f7d1b9e2772",
    "config/artifacts/models/sei/checkpoints.json": "76efd0bdd9241724f4e882134fee5b1b1228af9b580cc55ab6e3d987ece9ede8",
    "config/artifacts/models/sei/development_crosswalk.json": "0d37a5438747f654c5a386860ef8d029dde9b8404c300029f7309257a58f1c84",
    "config/artifacts/models/sei/exposure_audit.json": "9c7be7b56b41e7f0108b25521e0809810094e9c20f9c7f17e916d3362d7ccb66",
}


ARTIFACT_SHA256 = {
    "executions/borzoi-enformer-admission-sources-21064438": "c6af0a8723e5a4d8186af52c277d935c5ec4a2f4c46037b3f6954f9f77c39487",
    "executions/enformer-crested-restricted-port-admission-21065579": "e69877c43f47ce7d40962129dae7819a66123ac3f5bef7e0283ac654f77006a5",
    "executions/enformer-crested-safe-keras-inspection-21065759": "155537e853b730d440244581fda83ddff9ba8c13c9dabaa1f57bf841bb69196f",
    "executions/enformer-crested-track-crosswalk-21065831": "a658a159b03fb05c9a468b4fdbc9c4d745916ad3579b7752fcc4254074dbea10",
    "executions/environments/enformer-crested-tensorflow-2.17.0-21065799": "3534aee9e54d8fe6dc09c01c63aaead911590b97a08a7a686d56024520885ec6",
    "executions/enformer-crested-restricted-forward-21066179": "2bb24f274fd58f3011eef9e221e395392594391f57a8a17f8e9ecb7a510771c3",
    "executions/enformer-crested-restricted-forward-21066184": "f2b2d04685789e5b297b7150cc8e47747d54a99ebe287d35e4e51034a97c0341",
    "executions/gse281364-enformer-fixture-21070094": "0faf9d48511a76c282bfc412195deebda63503fb5685cd3fb2bfd77c72505231",
    "executions/gse281364-enformer-fixture-validation-21070139": "20d35712aa9ffbba6708c82731c76b9b4c8b90d214e736646f62c30a4868c198",
    "executions/model-work-203-21070998": "5b4940a0afca2276900c71bc08d247bc699180271caf0e653da5789312e580d6",
    "executions/alphagenome-sei-admission-sources-21064512": "84ae3470dca123c20904e315ef9926d5a139fb9fe6d99c51f6aaf694ef6e6785",
    "executions/sei-safe-native-archive-21064521": "beff13e916cf12eb6db2e29cb301b3151cc519048293a2ada2a66f6fa3691d2d",
    "executions/sei-scoring-sources-21065483": "529a2034a4a8168d8e436c7044364c65735bb274a43b9652d65000110562c5f7",
    "executions/sei-weights-only-conversion-21065468": "2b19a580f496db0a30a6174cc4f7fe65b0cf399e5dde94c8087665d374ba892b",
    "executions/sei-native-probe-21066123": "c0607cba297ad2e3c2889ebd4825b287655e6f0064f0d9137aaa60be34026894",
    "executions/gse281364-sei-fixture-21068925": "620464c0ef7cf378eb2b44a1608cb6b9c273bd70b8770a43b2a95d21e725af4d",
    "executions/sei-gse281364-smoke-21069658": "dd22116a2927712d1c7cafc1a3ca654783efce2173c44f010cd9fc22034c1803",
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise LaneAuditError(f"JSON authority is unreadable: {path}: {error}") from error
    if not isinstance(value, dict):
        raise LaneAuditError(f"expected JSON object: {path}")
    return value


def require_artifact(root: Path, expected: str) -> dict[str, Any]:
    try:
        verified = verify_frozen_tree(root)
    except ArtifactError as error:
        raise LaneAuditError(f"frozen artifact differs: {root}: {error}") from error
    if digest(root / "ARTIFACTS.json") != expected:
        raise LaneAuditError(f"artifact manifest identity differs: {root}")
    return dict(verified)


def require_authority(root: Path, relative: str, expected: str) -> dict[str, Any]:
    path = root / relative
    if digest(path) != expected:
        raise LaneAuditError(f"authority identity differs: {relative}")
    return load_json(path)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise LaneAuditError(message)


def model_registry(root: Path) -> dict[str, dict[str, Any]]:
    path = root / "config/models/regulatory_sequence.toml"
    try:
        raw = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise LaneAuditError(f"model registry is unreadable: {error}") from error
    rows = raw.get("models")
    if not isinstance(rows, list):
        raise LaneAuditError("model registry rows differ")
    result = {
        str(row.get("model_id")): row for row in rows if isinstance(row, dict)
    }
    if len(result) != len(rows):
        raise LaneAuditError("model registry has duplicate or malformed rows")
    return result


def audit(root: Path) -> dict[str, Any]:
    authorities = {
        relative: require_authority(root, relative, expected)
        for relative, expected in AUTHORITY_SHA256.items()
    }
    artifacts = {
        relative: require_artifact(root / relative, expected)
        for relative, expected in ARTIFACT_SHA256.items()
    }
    registry = model_registry(root)

    enformer = authorities["config/artifacts/models/enformer/checkpoints.json"]
    enformer_exposure = authorities[
        "config/artifacts/models/enformer/exposure_audit.json"
    ]
    enformer_registry = registry.get("enformer", {})
    enformer_source = load_json(
        root
        / "executions/borzoi-enformer-admission-sources-21064438/review/source_admission_receipt.json"
    )
    enformer_port = load_json(
        root
        / "executions/enformer-crested-restricted-port-admission-21065579/enformer_crested_port_admission_receipt.json"
    )
    enformer_crosswalk = load_json(
        root
        / "executions/enformer-crested-track-crosswalk-21065831/crosswalk/receipt.json"
    )
    enformer_forward = load_json(
        root
        / "executions/enformer-crested-restricted-forward-21066179/probe/receipt.json"
    )
    enformer_forward_repeat = load_json(
        root
        / "executions/enformer-crested-restricted-forward-21066184/probe/receipt.json"
    )
    enformer_fixture = load_json(
        root / "executions/gse281364-enformer-fixture-21070094/fixture/receipt.json"
    )
    enformer_prediction = load_json(
        root / "executions/model-work-203-21070998/predictions/receipt.json"
    )

    require(enformer["architecture_contract"]["input_length_bp"] == 196608,
            "Enformer input length differs")
    require(enformer["architecture_contract"]["output_bins"] == 896,
            "Enformer output-bin count differs")
    require(enformer["architecture_contract"]["human_output_tracks"] == 5313,
            "Enformer track count differs")
    require(enformer["input_contract"]["donor_or_cell_context"] == "none",
            "Enformer context contract differs")
    require(enformer["registered_converted_checkpoint"]["archive_sha256"]
            == "628f67f540304d4d0e143176dc824ed72b3413f78f2fa2efe5d4f0ab51ea1bcc",
            "Enformer converted archive differs")
    require(enformer["native_vs_converted_parity"]["parity_status"] == "UNRESOLVED",
            "Enformer native parity must remain unresolved")
    require("UNDECLARED" in enformer["terms_audit"]["native_checkpoint_terms"],
            "Enformer native-weight terms disposition differs")
    require(enformer_registry.get("checkpoint_sha256")
            == "628f67f540304d4d0e143176dc824ed72b3413f78f2fa2efe5d4f0ab51ea1bcc",
            "Enformer registry checkpoint identity differs")
    require(enformer_registry.get("status") == "restricted_comparator"
            and enformer_registry.get("admission_blocking") is True,
            "Enformer registry must remain fail closed")
    require(enformer_source["weight_terms"]["enformer"] == "UNDECLARED_BLOCKED"
            and not enformer_source["checkpoint_bytes_downloaded"],
            "official native Enformer admission differs")
    require(enformer_port["safe_member_inventory_passed"]
            and enformer_port["archive_sha256"]
            == "628f67f540304d4d0e143176dc824ed72b3413f78f2fa2efe5d4f0ab51ea1bcc"
            and not enformer_port["native_parity_established"]
            and not enformer_port["open_champion_eligible"],
            "restricted Enformer port evidence differs")
    require(enformer_crosswalk["canonical_target_rows"] == 5313
            and enformer_crosswalk["native_prediction_contract"]["output_bins"] == 896
            and enformer_crosswalk["native_prediction_contract"]["allele_delta_sign"]
            == "ALT_minus_REF",
            "Enformer output crosswalk differs")
    require(enformer_forward["checkpoint_member_sha256"]
            == "29dc3835c1d13a6c6bd93315b554075107cc409a9179b11ee0e00180f67cef56"
            and enformer_forward["deterministic_repeat_bit_identical"]
            and enformer_forward["predictions_artifact_sha256"]
            == enformer_forward_repeat["predictions_artifact_sha256"]
            and not enformer_forward["native_sonnet_parity_established"],
            "Enformer fixed-fixture forward evidence differs")
    require(enformer_fixture["elements"] == 1033
            and enformer_fixture["outer_locus_sequence_groups"] == 1033
            and enformer_fixture["outer_folds"] == 5
            and not enformer_fixture["outcomes_read"],
            "Enformer development fixture topology differs")
    require(enformer_prediction["status"]
            == "pass_outcome_blind_restricted_prediction"
            and enformer_prediction["elements"] == 1033
            and enformer_prediction["tracks"] == 5313
            and enformer_prediction["accessibility_tracks"] == 684
            and enformer_prediction["deterministic_repeat_max_abs"] == 0.0
            and not enformer_prediction["outcomes_read"]
            and not enformer_prediction["model_fitted_or_adapted"]
            and not enformer_prediction["native_sonnet_parity_established"]
            and not enformer_prediction["champion_eligible"],
            "Enformer outcome-blind development smoke differs")
    require(enformer_exposure["task_disposition"]["rna_conditioned_atac"]
            == "static_sequence_baseline_only",
            "Enformer RNA-conditioned disposition differs")

    sei = authorities["config/artifacts/models/sei/checkpoints.json"]
    sei_exposure = authorities["config/artifacts/models/sei/exposure_audit.json"]
    sei_registry = registry.get("sei", {})
    sei_source = load_json(
        root
        / "executions/alphagenome-sei-admission-sources-21064512/review/source_admission_receipt.json"
    )
    sei_archive = load_json(
        root / "executions/sei-safe-native-archive-21064521/inventory/admission_summary.json"
    )
    sei_scoring = load_json(
        root / "executions/sei-scoring-sources-21065483/validation/scoring_contract.json"
    )
    sei_conversion = load_json(
        root
        / "executions/sei-weights-only-conversion-21065468/converted/conversion_receipt.json"
    )
    sei_probe = load_json(
        root / "executions/sei-native-probe-21066123/predictions/probe_receipt.json"
    )
    sei_fixture = load_json(
        root / "executions/gse281364-sei-fixture-21068925/fixture/receipt.json"
    )
    sei_prediction = load_json(
        root / "executions/sei-gse281364-smoke-21069658/predictions/receipt.json"
    )

    require(sei["architecture_contract"]["input_length_bp"] == 4096
            and sei["architecture_contract"]["native_output_dim"] == 21907
            and sei["architecture_contract"]["sequence_class_dim"] == 40,
            "Sei architecture contract differs")
    require(sei["input_contract"]["donor_or_cell_context"] == "none",
            "Sei context contract differs")
    require(sei["native_checkpoint"]["archive"]["sha256"]
            == "3e1c41468d4b01d312aa44ebe0d2a708a705344908efb60daf8ff76548255b46",
            "Sei archive identity differs")
    require(sei["conversion_contract"]["converted_checkpoint_registered"]
            and sei["conversion_contract"]["converted_safetensors_sha256"]
            == "d0168411589bc47265e286b1224ab0564d4439c6d220a929f144924d1840a31e",
            "Sei safe conversion contract differs")
    require(sei["terms_audit"]["result"] == "restricted_comparator"
            and sei["weight_license"] == "CC-BY-4.0",
            "Sei code/weight terms distinction differs")
    require(sei_registry.get("checkpoint_sha256")
            == "3e1c41468d4b01d312aa44ebe0d2a708a705344908efb60daf8ff76548255b46"
            and sei_registry.get("status") == "restricted_comparator"
            and sei_registry.get("admission_blocking") is True,
            "Sei registry must remain fail closed")
    require(sei_source["remote_checkpoint_metadata"]["sei"]["record_license"]
            == "cc-by-4.0",
            "Sei Zenodo weight terms differ")
    require(sei_archive["archive_sha256"]
            == "3e1c41468d4b01d312aa44ebe0d2a708a705344908efb60daf8ff76548255b46"
            and sei_archive["expected_native_members"]["model/sei.pth"]["sha256"]
            == "9d4771d5363e4955eec2a36567fedac986191f64c7874c85bfb9e203c7f60783"
            and not sei_archive["pickle_loaded"],
            "Sei archive admission differs")
    require(sei_conversion["loader"] == "torch.load_weights_only_true_mmap_false"
            and sei_conversion["converted_safetensors_sha256"]
            == "d0168411589bc47265e286b1224ab0564d4439c6d220a929f144924d1840a31e"
            and sei_conversion["converted_tensor_key_count"] == 38
            and sei_conversion["converted_tensor_element_count"] == 889979983
            and not sei_conversion["custom_model_code_imported"],
            "Sei weights-only conversion differs")
    require(sei_scoring["native_variant_output_count"] == 40
            and sei_scoring["projection_array_shape"] == [61, 21907]
            and sei_scoring["native_variant_mapping"]
            == "first_40_projection_rows_after_histone_normalization_alt_minus_ref",
            "Sei native scoring contract differs")
    require(sei_probe["strict_checkpoint_restore"]
            and sei_probe["profile_shape"] == [12, 21907]
            and sei_probe["native_variant_score_shape"] == [3, 40]
            and sei_probe["input_gradient_backward_finite"]
            and sei_probe["deterministic_repeat_max_abs_diff"] == 0.0,
            "Sei native execution probe differs")
    require(sei_fixture["elements"] == 1033
            and sei_fixture["sequences"] == 4132
            and sei_fixture["outer_locus_sequence_groups"] == 1033
            and not sei_fixture["outcomes_read"],
            "Sei development fixture topology differs")
    require(sei_prediction["status"] == "pass_outcome_blind_native_prediction"
            and sei_prediction["elements"] == 1033
            and sei_prediction["native_sequence_classes"] == 40
            and sei_prediction["repeat_max_abs"] == 0.0
            and not sei_prediction["outcomes_read"]
            and not sei_prediction["model_fitted_or_adapted"]
            and not sei_prediction["champion_eligible"],
            "Sei outcome-blind development smoke differs")
    require(sei_exposure["task_disposition"]["rna_conditioned_atac"]
            == "static_sequence_center_bin_comparator_only",
            "Sei RNA-conditioned disposition differs")

    return {
        "schema_version": "masld-bench-enformer-sei-native-lane-audit-v1",
        "status": "pass",
        "authority_sha256": AUTHORITY_SHA256,
        "artifact_manifest_sha256": ARTIFACT_SHA256,
        "enformer": {
            "official_native_checkpoint_admitted": False,
            "official_native_weight_terms_resolved": False,
            "restricted_conversion_executable": True,
            "native_sonnet_parity_established": False,
            "input_length_bp": 196608,
            "native_output": "896x5313 fixed assay-track profiles",
            "development_smoke": "GSE281364 complete 1033-group outcome-blind REF/ALT/RC SAD-SAR fixture",
            "development_smoke_passed": True,
            "rna_conditioned": False,
            "open_champion_eligible": False,
            "terminal_disposition": "restricted_sequence_track_comparator_complete; official_native_lane_blocked",
        },
        "sei": {
            "official_archive_admitted": True,
            "safe_weights_only_conversion": True,
            "strict_native_source_forward_passed": True,
            "input_length_bp": 4096,
            "native_output": "21907 center-bin profile probabilities plus 40 projected sequence classes",
            "development_smoke": "GSE281364 complete 1033-group outcome-blind REF/ALT/RC class-score fixture",
            "development_smoke_passed": True,
            "rna_conditioned": False,
            "open_champion_eligible": False,
            "terminal_disposition": "restricted_sequence_center_bin_comparator_complete",
        },
        "task_separation": {
            "shared_development_source": "gse281364",
            "shared_biological_role": "sequence-only retrospective MPRA allele-effect diagnostic",
            "outputs_are_not_interchangeable": True,
            "observed_atac_input_models": [],
            "rna_conditioned_models": [],
            "signed_gene_effect_models": [],
        },
        "queue_disposition": {
            "new_gpu_job_required": False,
            "reason": "Both restricted development smokes are already complete; official Enformer Sonnet remains blocked by terms and admission rather than compute.",
        },
        "project_data_read_by_audit": False,
        "outcomes_read_by_audit": False,
        "sealed_labels_read_by_audit": False,
        "artifact_metadata_verified": sorted(artifacts),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve(strict=True)
    if args.output.exists():
        raise LaneAuditError("output already exists")
    result = audit(root)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "lane_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
