#!/usr/bin/env python3
"""Audit the frozen scBasset five-fold development campaign without refitting."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


EXPECTED_ARTIFACTS = {
    "matrix_gate": (
        "executions/scbasset-diagonal-matrix-gate-21067170",
        "1b566cabc831f775c569b573bf2c76191fa66c8ceaf83c69ade9078ff1a3ae13",
    ),
    "source_lock_superseded": (
        "executions/scbasset-crossfold-evaluator-source-21070148",
        "9b3eb31e66840ea79acd905c5a13d159de9191ae43339b3977793e2af5f27f7a",
    ),
    "source_lock_active": (
        "executions/scbasset-crossfold-evaluator-source-21070472",
        "be05e2b17d22de121880ca91b3b9750b0800f049731c05d1e243d7ad0441e9ad",
    ),
    "chain_fold1to4": (
        "executions/scbasset-fold1to4-chain-21069969",
        "9b0869623baab7b7062aa31e9ef794fa4529a404fe07f1fd759ab4a7f85c9056",
    ),
    "prediction_lock_fold1to4": (
        "executions/scbasset-fold1to4-prediction-lock-21070492",
        "9686a37c7bbeeac898cc2c4e1162c7033fcf2e0d09088abc9bb6c6c97bef3f07",
    ),
    "evaluation_fold0": (
        "executions/scbasset-donor0-genomic0-seed11-valid-evaluation-21065092",
        "7ab89892e508c145f7b3bae8819af0f4695624294cfcd22f1200984cdcba8341",
    ),
    "evaluation_fold1to4": (
        "executions/scbasset-fold1to4-valid-crossfold-evaluation-21070600",
        "2fb574e16edc290b19ca3420fc5258a766c1bce5a7e98550449420e0ac5f992c",
    ),
    "aggregate_all5": (
        "executions/scbasset-all5-valid-meta-aggregate-21070675",
        "79ecbd2d473c84ebb75e5e5476dcb4059ac27fe1d061263c733e9d35dc1f9a3f",
    ),
}

FOLD_ARTIFACTS = {
    0: {
        "input": (
            "executions/scbasset-donor0-genomic0-inputs-21064455",
            "02f76f4f7af62063a606943e2b1bcf11ca0c8d8a63f29879ed31bfc1fbc30dc0",
        ),
        "model": (
            "executions/scbasset-donor0-genomic0-seed11-21064694",
            "06db5a37eed92ebb6fd95b47f10cb998121cc0a641cab42686c27d47660e9e00",
        ),
        "prediction": (
            "executions/scbasset-donor0-genomic0-seed11-valid-21064859",
            "c27f080243968ec5815e158fc22280c1718ef9b0cb0ea08012442755d35a9db0",
        ),
        "seed": 11,
    },
    1: {
        "input": (
            "executions/scbasset-donor1_genomic1-inputs-21067171",
            "5a608a40763e0d72e0af1a284690c762f54379a8230e06889361176935bb88c3",
        ),
        "model": (
            "executions/scbasset-donor1_genomic1-seed20260824-21069969",
            "c77bbb004bccf1856008c409377bdb68e81fef93d16e3f860454a7167099aee2",
        ),
        "prediction": (
            "executions/scbasset-donor1_genomic1-seed20260824-valid-21069969",
            "95410886fd5f97f28f9a4d4b84924ae65cf93a246ab6be6d4699b67560878459",
        ),
        "seed": 20260824,
    },
    2: {
        "input": (
            "executions/scbasset-donor2_genomic2-inputs-21067172",
            "1b449406fba08190a94bb2750eaf558e7014b187fb239c998659df68561f610f",
        ),
        "model": (
            "executions/scbasset-donor2_genomic2-seed20260824-21069969",
            "4933e2d2ad64be526da93c25e3b3a723a175026d40fd3affdbc03af5d0d85ff8",
        ),
        "prediction": (
            "executions/scbasset-donor2_genomic2-seed20260824-valid-21069969",
            "eed430b0f559005e9857e6a3f84d16658829a0c2509b20a071b4e0a4528c12e3",
        ),
        "seed": 20260824,
    },
    3: {
        "input": (
            "executions/scbasset-donor3_genomic3-inputs-21067173",
            "bb9b1ce8f4e84ddd4509e97ee593ddca5a39e4b66556aa228c057c50a35f9678",
        ),
        "model": (
            "executions/scbasset-donor3_genomic3-seed20260824-21069969",
            "2f2a43cd57742e57188d758447785591643482bc9bc36df97df4a29c702a0c17",
        ),
        "prediction": (
            "executions/scbasset-donor3_genomic3-seed20260824-valid-21069969",
            "aa7df5227a84535b40a6ed49a2d0809f2c0258bc312c0990021c63059775907e",
        ),
        "seed": 20260824,
    },
    4: {
        "input": (
            "executions/scbasset-donor4_genomic4-inputs-21067174",
            "4796aa202ec5f2a171f7422fcb80544d66f9e288372816a3e856cf688e8d4505",
        ),
        "model": (
            "executions/scbasset-donor4_genomic4-seed20260824-21069969",
            "ac29a200e2d23136451e5e4f349978448320c02fe099b064ca9460b28f93cae3",
        ),
        "prediction": (
            "executions/scbasset-donor4_genomic4-seed20260824-valid-21069969",
            "a91f5f7fe771312b122d50e7f800aaa337d7c5fadec14eff3d7bb7ecf91756c5",
        ),
        "seed": 20260824,
    },
}


class ScBassetCampaignAuditError(RuntimeError):
    """Raised when the frozen campaign closure differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ScBassetCampaignAuditError(f"missing or linked JSON: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ScBassetCampaignAuditError(f"JSON object required: {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    if path.is_symlink() or not path.is_file():
        raise ScBassetCampaignAuditError(f"missing or linked TSV: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def verify_root(root: Path, relative: str, expected: str) -> Path:
    path = (root / relative).resolve(strict=True)
    if sha256_file(path / "ARTIFACTS.json") != expected:
        raise ScBassetCampaignAuditError(f"ARTIFACTS identity differs: {relative}")
    return path


def terminal_disposition(aggregate: dict[str, Any]) -> dict[str, Any]:
    gate = aggregate.get("development_gate", {})
    lineage = aggregate.get("lineage_relative_deviance_reduction", {})
    if (
        aggregate.get("outer_folds") != [0, 1, 2, 3, 4]
        or aggregate.get("n_donors") != 39
        or aggregate.get("mixed_seed_meta_aggregate") is not True
        or aggregate.get("fixed_seed_five_fold_cv") is not False
        or gate.get("overall_threshold_passed") is not False
        or gate.get("improved_lineages") != 0
        or set(lineage)
        != {"cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell"}
        or any(float(value) >= 0.0 for value in lineage.values())
    ):
        raise ScBassetCampaignAuditError("all-fold development disposition differs")
    return {
        "development_gate_passed": False,
        "gpu_production_authorized": False,
        "terminal_disposition": "stop_after_all5_validation_loss_to_training_only_baselines",
    }


def audit(root: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise ScBassetCampaignAuditError("refusing to overwrite campaign audit")
    authorities = {
        name: verify_root(root, relative, expected)
        for name, (relative, expected) in EXPECTED_ARTIFACTS.items()
    }
    test_sets: list[set[str]] = []
    valid_sets: list[set[str]] = []
    fold_rows = []
    for fold, contract in FOLD_ARTIFACTS.items():
        input_root = verify_root(root, *contract["input"])
        model_root = verify_root(root, *contract["model"])
        prediction_root = verify_root(root, *contract["prediction"])
        summary = read_json(input_root / "inputs/summary.json")
        model = read_json(model_root / "model/model_receipt.json")
        prediction = read_json(
            prediction_root / "predictions/prediction_receipt.json"
        )
        test = read_tsv(input_root / "inputs/held_test_donors.tsv")
        valid = read_tsv(input_root / "inputs/held_valid_donors.tsv")
        training = read_tsv(input_root / "inputs/training_cells.tsv")
        regions = read_json(input_root / "regions/fold.json")
        test_donors = {row["donor_id"] for row in test}
        valid_donors = {row["donor_id"] for row in valid}
        training_donors = {row["donor_id"] for row in training}
        training_folds = {int(row["outer_fold"]) for row in training}
        expected_training_folds = set(range(5)) - {fold, (fold + 1) % 5}
        split_id = f"donor{fold}_genomic{fold}"
        if (
            summary.get("status") != "pass"
            or summary.get("split_id") != split_id
            or summary.get("donor_test_fold") != fold
            or summary.get("donor_valid_fold") != (fold + 1) % 5
            or set(summary.get("donor_train_folds", [])) != expected_training_folds
            or summary.get("held_donor_atac_read") is not False
            or summary.get("held_donor_atac_exported") is not False
            or summary.get("genomic_test_atac_used_for_matrix") is not False
            or summary.get("model_training_input_contains_genomic_test_atac")
            is not False
            or summary.get("missing_evidence_encoded_as_zero") is not False
            or test_donors & valid_donors
            or test_donors & training_donors
            or valid_donors & training_donors
            or training_folds != expected_training_folds
            or len(training_donors) != summary.get("training_donors")
            or any(int(row["outer_fold"]) != fold for row in test)
            or any(int(row["outer_fold"]) != (fold + 1) % 5 for row in valid)
            or set(regions) != {"train", "valid", "test"}
            or set(regions["train"]) & set(regions["valid"])
            or set(regions["train"]) & set(regions["test"])
            or set(regions["valid"]) & set(regions["test"])
        ):
            raise ScBassetCampaignAuditError(f"fold {fold} donor/genomic firewall differs")
        if (
            model.get("status") != "pass"
            or model.get("split_id") != split_id
            or model.get("seed") != contract["seed"]
            or model.get("input_artifacts_sha256") != contract["input"][1]
            or model.get("held_donor_atac_used") is not False
            or model.get("genomic_test_atac_used") is not False
            or prediction.get("status") != "pass"
            or prediction.get("role") != "valid"
            or prediction.get("split_id") != split_id
            or prediction.get("input_artifacts_sha256") != contract["input"][1]
            or prediction.get("model_artifacts_sha256") != contract["model"][1]
            or prediction.get("held_donor_atac_used") is not False
            or prediction.get("rna_used") is not False
        ):
            raise ScBassetCampaignAuditError(f"fold {fold} fit/prediction binding differs")
        test_sets.append(test_donors)
        valid_sets.append(valid_donors)
        fold_rows.append(
            {
                "outer_fold": fold,
                "validation_fold": (fold + 1) % 5,
                "seed": contract["seed"],
                "test_donors": len(test_donors),
                "validation_donors": len(valid_donors),
                "training_donors": len(training_donors),
                "input_artifacts_sha256": contract["input"][1],
                "model_artifacts_sha256": contract["model"][1],
                "prediction_artifacts_sha256": contract["prediction"][1],
            }
        )
    if (
        len(set().union(*test_sets)) != 39
        or sum(map(len, test_sets)) != 39
        or len(set().union(*valid_sets)) != 39
        or sum(map(len, valid_sets)) != 39
    ):
        raise ScBassetCampaignAuditError("five-fold donor partition differs")

    active_source = read_json(authorities["source_lock_active"] / "source_lock.json")
    superseded_source = read_json(
        authorities["source_lock_superseded"] / "source_lock.json"
    )
    prediction_lock = read_json(
        authorities["prediction_lock_fold1to4"] / "prediction_lock.json"
    )
    fold0_evaluation = read_json(
        authorities["evaluation_fold0"] / "evaluation/evaluation.json"
    )
    fold1to4_evaluation = read_json(
        authorities["evaluation_fold1to4"] / "evaluation/evaluation.json"
    )
    aggregate = read_json(
        authorities["aggregate_all5"] / "evaluation/evaluation.json"
    )
    if (
        active_source.get("folds") != [1, 2, 3, 4]
        or active_source.get("outcomes_read") is not False
        or "lock_unit_test" not in active_source.get("records", {})
        or "lock_unit_test" in superseded_source.get("records", {})
        or prediction_lock.get("folds") != [1, 2, 3, 4]
        or prediction_lock.get("outcomes_read") is not False
        or prediction_lock.get("outcome_paths_recorded") is not False
        or len(prediction_lock.get("records", [])) != 4
        or fold0_evaluation.get("test_atac_read") is not False
        or fold0_evaluation.get("prediction_frozen_before_outcomes") is not True
        or fold1to4_evaluation.get("prediction_lock_verified_before_outcomes")
        is not True
        or fold1to4_evaluation.get("test_atac_read") is not False
        or aggregate.get("outcome_files_opened_or_recomputed") is not False
        or aggregate.get("test_atac_read") is not False
        or aggregate.get("champion_claim_allowed") is not False
    ):
        raise ScBassetCampaignAuditError("evaluator/lock firewall differs")

    result = {
        "schema_version": "masld-bench-scbasset-all5-campaign-audit-v1",
        "status": "pass",
        "model_id": "scbasset",
        "dataset_id": "gse296875",
        "evaluation_role": "development_validation",
        "biological_unit": "donor",
        "outer_folds_present": [0, 1, 2, 3, 4],
        "donors_partitioned_exactly_once_as_test": True,
        "donors_partitioned_exactly_once_as_validation": True,
        "donor_safe_five_fold_artifact_campaign_complete": True,
        "uniform_prediction_lock_five_fold_campaign_complete": False,
        "uniform_prediction_lock_limitation": "fold_0_is_legacy_seed11_without_the_fold1to4_prediction_lock_protocol",
        "fixed_seed_three_seed_screen_complete": False,
        "observed_models": 5,
        "prespecified_screen_models": 15,
        "seed_by_outer_fold": {str(row["outer_fold"]): row["seed"] for row in fold_rows},
        "folds": fold_rows,
        "active_evaluator_source_lock_artifacts_sha256": EXPECTED_ARTIFACTS[
            "source_lock_active"
        ][1],
        "superseded_evaluator_source_lock_artifacts_sha256": EXPECTED_ARTIFACTS[
            "source_lock_superseded"
        ][1],
        "prediction_lock_fold1to4_artifacts_sha256": EXPECTED_ARTIFACTS[
            "prediction_lock_fold1to4"
        ][1],
        "all5_aggregate_artifacts_sha256": EXPECTED_ARTIFACTS["aggregate_all5"][1],
        "frozen_metric_summary_read": True,
        "raw_outcome_authority_opened": False,
        "outcomes_used_for_fitting": False,
        "sealed_data_read": False,
        "champion_claim_allowed": False,
        **terminal_disposition(aggregate),
    }
    output.mkdir(parents=True, mode=0o750)
    (output / "campaign_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    result = audit(arguments.root, arguments.output)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
