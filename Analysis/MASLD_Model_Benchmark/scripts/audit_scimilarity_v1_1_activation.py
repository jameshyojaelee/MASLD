#!/usr/bin/env python3
"""Audit the exact SCimilarity v1.1 assets used by the frozen extraction lane."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import tomllib


ARCHIVE_SHA256 = "3904fb66969bca0dd9b162955d9092f497b0769b6975caff70a17257901e68a1"
ARCHIVE_MD5 = "546251b7c435f3b1dbe38e2e420ad57f"
ARCHIVE_SIZE = 30_310_810_843
ENCODER_SHA256 = "a08f023788bdbded191e2d9344e5d4240ffc7d7505eca4e0b7a37526b858d19f"
GENE_ORDER_SHA256 = "0a9aeacf7fc19399426273d8cfab8683253ec56b51be556d42edd115f125b8ca"
CODE_LICENSE_SHA256 = "938fb579905acc3f45d10f20cd18da74f2e754300220f9669d0e5c5698f85cac"
NN_MODELS_SHA256 = "5b3edf1a57ae653202d9655e049b516ec31cfa14e81daaab7c96b647fec6c110"
REPOSITORY_REVISION = "3ce3ec81e122a115eb3e667706725120f7ee252d"


class SCimilarityAuditError(RuntimeError):
    """Raised when activation evidence differs from the frozen requirements."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def audit(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.archive.stat().st_size != ARCHIVE_SIZE or sha256_file(arguments.archive) != ARCHIVE_SHA256:
        raise SCimilarityAuditError("SCimilarity v1.1 archive identity differs")
    required_model_hashes = {
        "encoder.ckpt": ENCODER_SHA256,
        "gene_order.tsv": GENE_ORDER_SHA256,
        "hyperparameters.json": "e0a341792671dbce83b8069cd05552d40c5ce90581136867c8e821128b5dfb34",
        "label_ints.csv": "c822aa1201fa8b3a60cf384b2ee72d74db80bf4954fcd091a6e5ee1b321646de",
        "layer_sizes.json": "f0ef0291ff2ec05d45689cca53a052b3ec20dfda8bf7db94bd74649b1a70b46f",
        "metadata.json": "791ebc4f1636f9327a3c33fa2feb312c3bac0a41f24e57326747e8b222ec4cf6",
    }
    for name, expected in required_model_hashes.items():
        if sha256_file(arguments.model_root / name) != expected:
            raise SCimilarityAuditError(f"staged model file {name} differs")
    genes = (arguments.model_root / "gene_order.tsv").read_text(encoding="utf-8").splitlines()
    if len(genes) != 28_231 or len(set(genes)) != len(genes):
        raise SCimilarityAuditError("SCimilarity gene order cardinality differs")
    if sha256_file(arguments.code_license) != CODE_LICENSE_SHA256:
        raise SCimilarityAuditError("SCimilarity Apache-2.0 license file differs")
    if sha256_file(arguments.nn_models) != NN_MODELS_SHA256:
        raise SCimilarityAuditError("SCimilarity architecture source differs from pinned revision")

    checkpoint_config = json.loads(arguments.checkpoint_config.read_text(encoding="utf-8"))
    exposure = json.loads(arguments.exposure_audit.read_text(encoding="utf-8"))
    crosswalk = json.loads(arguments.crosswalk.read_text(encoding="utf-8"))
    with arguments.model_registry.open("rb") as handle:
        model_registry = tomllib.load(handle)
    registered = {
        item.get("model_id"): item for item in model_registry.get("models", [])
    }.get("scimilarity_v1_1")
    if (
        checkpoint_config.get("repository_revision") != REPOSITORY_REVISION
        or checkpoint_config.get("code_license") != "Apache-2.0"
        or checkpoint_config.get("weight_license") != "CC-BY-SA-4.0"
        or checkpoint_config.get("archive", {}).get("source_md5") != ARCHIVE_MD5
        or checkpoint_config.get("archive", {}).get("sha256") != ARCHIVE_SHA256
        or checkpoint_config.get("weight_content_downloaded") is not True
        or checkpoint_config.get("staged_model_contract", {}).get("encoder_sha256")
        != ENCODER_SHA256
    ):
        raise SCimilarityAuditError("SCimilarity repository/license/source contract differs")
    if (
        not isinstance(registered, dict)
        or registered.get("checkpoint_sha256") != ENCODER_SHA256
        or registered.get("license_status") != "code_Apache-2.0_weights_CC-BY-SA-4.0"
        or registered.get("exposure_status") != "encoder_seen"
        or registered.get("admission_blocking") is not False
        or registered.get("blockers") != []
    ):
        raise SCimilarityAuditError("SCimilarity central model registry differs")
    findings = crosswalk.get("findings", {})
    expected_exposure = {
        "GSE136103": "reference_only",
        "GSE174748": "clean_declared",
        "GSE185477": "encoder_seen",
        "GSE189600": "clean_declared",
        "GSE202379": "clean_declared",
        "GSE244832": "clean_declared",
        "Liver_Atlas": "clean_declared",
    }
    if {key: findings.get(key, {}).get("exposure_state") for key in expected_exposure} != expected_exposure:
        raise SCimilarityAuditError("SCimilarity development exposure crosswalk differs")
    if (
        exposure.get("checkpoint_findings", {}).get("scimilarity_v1_1", {}).get("exposure_state")
        != "encoder_seen"
        or exposure.get("checkpoint_findings", {}).get("scimilarity_v1_1", {}).get(
            "sealed_target_exposure_state"
        )
        != "target_label_unexposed"
        or exposure.get("sealed_champion_eligibility")
        != "ineligible_aggregate_development_encoder_seen"
        or exposure.get("sealed_champion_eligible") is not False
    ):
        raise SCimilarityAuditError("SCimilarity sealed exposure contract differs")
    fixture_receipt = json.loads(
        (arguments.fixture_root / "fixture_receipt.json").read_text(encoding="utf-8")
    )
    if (
        fixture_receipt.get("aggregate_exposure_status") != "encoder_seen"
        or fixture_receipt.get("sealed_outcomes_read") is not False
        or fixture_receipt.get("evaluation_label_columns_read")
        or fixture_receipt.get("study_gene_observability_order")
        != list(expected_exposure)
        or fixture_receipt.get("study_gene_observability_shape") != [7, 28_231]
        or fixture_receipt.get("study_gene_observability_is_global") is not True
        or fixture_receipt.get("gene_overlap") != 20_609
    ):
        raise SCimilarityAuditError("SCimilarity frozen fixture observability differs")
    fixture_contract = exposure.get("frozen_fixture_contract", {})
    if (
        fixture_contract.get("study_gene_observability_shape") != [7, 28_231]
        or fixture_contract.get("study_gene_observability_order") != list(expected_exposure)
        or fixture_contract.get("sealed_outcomes_read") is not False
        or fixture_contract.get("evaluation_labels_read") is not False
        or fixture_contract.get("study_gene_observability_sha256")
        != sha256_file(arguments.fixture_root / "study_gene_observability.npy")
        or fixture_contract.get("observed_gene_mask_sha256")
        != sha256_file(arguments.fixture_root / "observed_gene_mask.npy")
    ):
        raise SCimilarityAuditError("SCimilarity central observability registry differs")
    receipt = {
        "schema_version": "masld-bench-scimilarity-v1-1-activation-audit-v2",
        "status": "pass_exact_assets_terms_and_exposure",
        "model_id": "scimilarity_v1_1",
        "repository_revision": REPOSITORY_REVISION,
        "archive_sha256": ARCHIVE_SHA256,
        "archive_source_md5": ARCHIVE_MD5,
        "archive_size_bytes": ARCHIVE_SIZE,
        "checkpoint_sha256": ENCODER_SHA256,
        "gene_order_sha256": GENE_ORDER_SHA256,
        "gene_order_rows": len(genes),
        "code_license": "Apache-2.0",
        "code_license_sha256": CODE_LICENSE_SHA256,
        "weight_license": "CC-BY-SA-4.0",
        "weight_license_authority": "Zenodo record 10685499",
        "architecture_source_sha256": NN_MODELS_SHA256,
        "development_exposure": expected_exposure,
        "sealed_gse289173_exposure": "target_label_unexposed",
        "aggregate_development_exposure": "encoder_seen",
        "released_reference_index_permitted": False,
        "central_registry_reconciled": True,
        "central_registry_records": {
            "checkpoint_config_sha256": sha256_file(arguments.checkpoint_config),
            "development_crosswalk_sha256": sha256_file(arguments.crosswalk),
            "exposure_audit_sha256": sha256_file(arguments.exposure_audit),
            "model_registry_sha256": sha256_file(arguments.model_registry),
        },
        "fixture_rows": fixture_receipt["rows"],
        "fixture_studies": fixture_receipt["studies"],
        "study_gene_observability_order": fixture_receipt[
            "study_gene_observability_order"
        ],
        "study_gene_observability_shape": fixture_receipt[
            "study_gene_observability_shape"
        ],
        "study_gene_observability_is_global": True,
        "observed_genes": fixture_receipt["gene_overlap"],
        "study_gene_observability_sha256": sha256_file(
            arguments.fixture_root / "study_gene_observability.npy"
        ),
        "sealed_champion_eligible": False,
        "raw_unrefined_embeddings_only": True,
        "downstream_head_fit": False,
        "evaluation_labels_read": False,
        "sealed_outcomes_read": False,
    }
    arguments.output.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--model-root", type=Path, required=True)
    parser.add_argument("--code-license", type=Path, required=True)
    parser.add_argument("--nn-models", type=Path, required=True)
    parser.add_argument("--checkpoint-config", type=Path, required=True)
    parser.add_argument("--exposure-audit", type=Path, required=True)
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--model-registry", type=Path, required=True)
    parser.add_argument("--fixture-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    result = audit(parser.parse_args())
    print(json.dumps({key: result[key] for key in ("status", "checkpoint_sha256", "gene_order_rows")}))


if __name__ == "__main__":
    main()
