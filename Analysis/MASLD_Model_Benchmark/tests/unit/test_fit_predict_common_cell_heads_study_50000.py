#!/usr/bin/env python3
"""Unit checks for sharded, identity-bound cell-foundation common heads."""

from __future__ import annotations

import ast
import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file
from scripts import fit_predict_common_cell_heads_study_50000 as common


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts/fit_predict_common_cell_heads_study_50000.py"
SHARD_WRAPPER = ROOT / "slurm/run_common_cell_head_bundle_50000.sbatch"
AGGREGATE_WRAPPER = ROOT / "slurm/aggregate_common_cell_heads_study_50000.sbatch"
TRANSCRIPTFORMER_QUEUE = ROOT / "config/campaigns/gpu_bundle_queue/052-model-training-505.json"


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_embedding_root(
    parent: Path,
    *,
    metadata_exposure: str = "unknown",
    receipt_exposure: str = "unknown",
) -> Path:
    root = parent / "embedding_root"
    (root / "embeddings").mkdir(parents=True)
    (root / "embeddings/common_embeddings.npz").write_bytes(b"fixture")
    (root / "embeddings/receipt.json").write_text(
        json.dumps(
            {
                "checkpoint_sha256": "a" * 64,
                "embedding_width": 32,
                "evaluation_label_columns_read": [],
                "exposure_status": receipt_exposure,
                "policies": {"common": {}},
                "rows": 50_000,
                "sealed_outcomes_read": False,
            }
        ),
        encoding="utf-8",
    )
    freeze_tree(
        root,
        {
            "artifact_class": "synthetic_embeddings",
            "dataset_view_id": common.DATASET_VIEW_ID,
            "development_rows": 50_000,
            "evaluation_labels_read": False,
            "exposure_status": metadata_exposure,
            "model_id": "fixture_model",
            "sealed_champion_eligible": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )
    return root


class CommonCellHeadStudy50000Tests(unittest.TestCase):
    def test_donor_class_weights_equalize_class_and_donor_mass(self) -> None:
        donors = np.asarray(["a", "a", "b", "c", "d", "e"])
        labels = np.asarray([0, 0, 0, 1, 1, 2])
        values = common.donor_class_weights(
            donors, labels, require_full_roster=False
        )
        self.assertAlmostEqual(float(values.sum()), len(values), places=5)
        class_mass = [float(values[labels == label].sum()) for label in range(3)]
        self.assertAlmostEqual(class_mass[0], class_mass[1], places=5)
        self.assertAlmostEqual(class_mass[1], class_mass[2], places=5)
        self.assertAlmostEqual(
            float(values[(donors == "a") & (labels == 0)].sum()),
            float(values[(donors == "b") & (labels == 0)].sum()),
            places=5,
        )

    def test_identity_is_derived_from_frozen_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = make_embedding_root(Path(directory))
            identity, path = common.derive_embedding_identity(
                root, Path("embeddings/common_embeddings.npz")
            )
            self.assertEqual(identity["model_id"], "fixture_model")
            self.assertEqual(identity["representation_id"], "common")
            self.assertEqual(identity["checkpoint_sha256"], "a" * 64)
            self.assertEqual(identity["exposure_status"], "unknown")
            self.assertEqual(path.name, "common_embeddings.npz")

    def test_identity_rejects_exposure_spoof(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = make_embedding_root(
                Path(directory),
                metadata_exposure="unknown",
                receipt_exposure="clean_declared",
            )
            with self.assertRaisesRegex(common.CommonHeadError, "identity"):
                common.derive_embedding_identity(
                    root, Path("embeddings/common_embeddings.npz")
                )

    def test_embedding_path_escape_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            root = make_embedding_root(parent)
            (parent / "outside_embeddings.npz").write_bytes(b"outside")
            with self.assertRaisesRegex(common.CommonHeadError, "escapes"):
                common.derive_embedding_identity(
                    root, Path("../outside_embeddings.npz")
                )

    def test_study_and_donor_outer_firewalls(self) -> None:
        donors = np.asarray([f"d{index}" for index in range(10)])
        datasets = np.asarray([f"s{index // 2}" for index in range(10)])
        outer = np.repeat(np.arange(5), 2).astype(np.int8)
        common.validate_study_outer_split(donors, datasets, outer)
        crossed_study = datasets.copy()
        crossed_study[2] = crossed_study[0]
        with self.assertRaisesRegex(common.CommonHeadError, "study crosses"):
            common.validate_study_outer_split(donors, crossed_study, outer)
        crossed_donor = donors.copy()
        crossed_donor[2] = crossed_donor[0]
        with self.assertRaisesRegex(common.CommonHeadError, "donor crosses"):
            common.validate_study_outer_split(crossed_donor, datasets, outer)

    def test_inner_assignment_requires_exact_five_fold_coverage(self) -> None:
        donors = np.asarray([f"d{index}" for index in range(10)])
        assignment = {donor: index // 2 for index, donor in enumerate(donors)}
        observed = common.validate_inner_assignment(donors, assignment)
        np.testing.assert_array_equal(observed, np.repeat(np.arange(5), 2))
        invalid = dict(assignment)
        invalid["d9"] = 5
        with self.assertRaisesRegex(common.CommonHeadError, "outside"):
            common.validate_inner_assignment(donors, invalid)
        empty = {donor: min(index // 2, 3) for index, donor in enumerate(donors)}
        with self.assertRaisesRegex(common.CommonHeadError, "empty"):
            common.validate_inner_assignment(donors, empty)

    def test_fixed_epoch_policy_does_not_select_on_validation_labels(self) -> None:
        self.assertEqual(common.FIXED_EPOCHS, 30)
        self.assertEqual(common.SCREEN_SEEDS, (1103, 1201, 1301))
        tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
        function = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "_fit_fixed_head"
        )
        arguments = {argument.arg for argument in function.args.kwonlyargs}
        self.assertNotIn("validation_targets", arguments)
        self.assertNotIn("validation_y", ast.unparse(function))
        self.assertNotIn("validation_loss", ast.unparse(function))
        text = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("fit_indices=training", text)
        self.assertIn("prediction_indices=test", text)
        self.assertIn('"final_fit_rows": int(final_record["fit_rows"])', text)
        self.assertIn('"inner_models_used_for_outer_prediction": False', text)
        self.assertIn('"outer_prediction_uses_final_refit": True', text)
        self.assertEqual(text.count("from torch import nn"), 1)
        self.assertIn("probabilities /= probabilities.sum(axis=1, keepdims=True)", text)
        self.assertIn("rtol=0.0, atol=1.0e-12", text)

    def test_prediction_schema_has_no_outcome_fields(self) -> None:
        self.assertFalse(common.FORBIDDEN_PREDICTION_FIELDS & set(common.PREDICTION_FIELDS))
        self.assertNotIn("broad_label", common.PREDICTION_FIELDS)
        self.assertNotIn("histology", common.PREDICTION_FIELDS)

    def test_prediction_reader_rejects_outcomes_and_invalid_probabilities(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            labeled = root / "labeled.tsv"
            write_tsv(
                labeled,
                (*common.PREDICTION_FIELDS, "broad_label"),
                [],
            )
            with self.assertRaisesRegex(common.CommonHeadError, "fields"):
                common._read_prediction_shard(labeled)
            invalid = root / "invalid.tsv"
            row = {
                "row_id": "r0",
                "donor_id": "d0",
                "dataset": "s0",
                "outer_fold": 0,
                "predicted_class": common.ROSTER[0],
                **{f"probability::{label}": "0.1" for label in common.ROSTER},
            }
            write_tsv(invalid, common.PREDICTION_FIELDS, [row])
            with self.assertRaisesRegex(common.CommonHeadError, "values"):
                common._read_prediction_shard(invalid)

    def test_synthetic_aggregate_requires_and_joins_all_shards(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            source = parent / "source"
            split = parent / "split"
            source.mkdir()
            split.mkdir()
            rows = [
                {
                    "row_id": f"r{fold}",
                    "donor_id": f"d{fold}",
                    "dataset": f"s{fold}",
                    "outer_fold": fold,
                }
                for fold in common.OUTER_FOLDS
            ]
            write_tsv(
                split / "row_outer_folds.tsv",
                ("row_id", "donor_id", "dataset", "outer_fold"),
                rows,
            )
            freeze_tree(
                source,
                {
                    "subset_id": common.DATASET_VIEW_ID,
                    "row_count": 50_000,
                    "status": "passed",
                },
            )
            freeze_tree(
                split,
                {
                    "sealed_outcomes_read": False,
                    "split_id": common.SPLIT_ID,
                    "target_labels_used_for_assignment": False,
                    "status": "passed",
                },
            )
            source_sha = sha256_file(source / "ARTIFACTS.json")
            split_sha = sha256_file(split / "ARTIFACTS.json")
            identity = {
                "model_id": "fixture_model",
                "representation_id": "common",
                "checkpoint_sha256": "a" * 64,
                "exposure_status": "unknown",
                "embedding_width": 32,
                "embedding_relative": "embeddings/common_embeddings.npz",
                "sealed_champion_eligible": False,
            }
            manifest_rows = []
            for head in common.HEADS:
                for seed in common.SCREEN_SEEDS:
                    for fold in common.OUTER_FOLDS:
                        shard = parent / f"shard-{head}-{seed}-{fold}"
                        shard.mkdir()
                        probabilities = ["1", "0", "0", "0", "0"]
                        prediction = {
                            **rows[fold],
                            "predicted_class": common.ROSTER[0],
                            **{
                                f"probability::{label}": probabilities[index]
                                for index, label in enumerate(common.ROSTER)
                            },
                        }
                        write_tsv(
                            shard / "predictions.tsv",
                            common.PREDICTION_FIELDS,
                            [prediction],
                        )
                        (shard / "prediction_receipt.json").write_text(
                            json.dumps(
                                {
                                    "status": "pass_development_prediction_shard",
                                    **identity,
                                    "dataset_view_id": common.DATASET_VIEW_ID,
                                    "split_id": common.SPLIT_ID,
                                    "head_id": head,
                                    "screen_seed": seed,
                                    "outer_fold": fold,
                                    "training_rows": 4,
                                    "training_donors": 4,
                                    "training_studies": [
                                        f"s{value}" for value in common.OUTER_FOLDS if value != fold
                                    ],
                                    "fixed_epochs": common.FIXED_EPOCHS,
                                    "checkpoint_selection_uses_validation_labels": False,
                                    "inner_models_used_for_outer_prediction": False,
                                    "oof_exact_one_time_coverage": True,
                                    "final_fit_rows": 4,
                                    "final_fit_donors": 4,
                                    "final_fit_studies": 4,
                                    "final_fit_fixed_epochs": common.FIXED_EPOCHS,
                                    "outer_prediction_model": "fixed_epoch_final_outer_training_refit",
                                    "outer_prediction_uses_final_refit": True,
                                    "prediction_tables_contain_observed_labels": False,
                                    "sealed_outcomes_read": False,
                                    "input_artifacts_sha256": {
                                        "embeddings": "b" * 64,
                                        "source": source_sha,
                                        "split": split_sha,
                                    },
                                }
                            ),
                            encoding="utf-8",
                        )
                        freeze_tree(shard, {"status": "passed"})
                        manifest_rows.append(
                            {
                                "shard_root": shard.resolve().as_posix(),
                                "artifacts_sha256": sha256_file(shard / "ARTIFACTS.json"),
                            }
                        )
            manifest = parent / "shards.tsv"
            write_tsv(
                manifest,
                ("shard_root", "artifacts_sha256"),
                manifest_rows,
            )
            output = parent / "aggregate"
            common.aggregate_shards(
                shard_manifest=manifest,
                source=source,
                split=split,
                output=output,
                expected_source_artifacts_sha256=source_sha,
                expected_split_artifacts_sha256=split_sha,
            )
            self.assertEqual(len(list((output / "predictions").glob("*.tsv"))), 6)
            for path in (output / "predictions").glob("*.tsv"):
                self.assertEqual(len(common._read_prediction_shard(path)), 5)

    def test_production_wrappers_are_nslab_masked_and_not_arrays(self) -> None:
        for wrapper in (SHARD_WRAPPER, AGGREGATE_WRAPPER):
            text = wrapper.read_text(encoding="utf-8")
            header = "\n".join(
                line for line in text.splitlines() if line.startswith("#SBATCH")
            )
            self.assertIn("--account=nslab", header)
            self.assertIn("--qos=nslab", header)
            self.assertNotIn("--array", header)
            self.assertNotIn("innovation", text.lower())
            job_name = next(
                line for line in text.splitlines() if line.startswith("#SBATCH --job-name=")
            )
            self.assertNotIn("masld", job_name.lower())
        shard_text = SHARD_WRAPPER.read_text(encoding="utf-8")
        self.assertIn("--partition=gpu", shard_text)
        self.assertIn("--gres=gpu:l40s:1", shard_text)
        self.assertIn("--time=72:00:00", shard_text)
        self.assertIn("CUBLAS_WORKSPACE_CONFIG=:4096:8", shard_text)
        self.assertIn("PYTHONHASHSEED=20260824", shard_text)
        self.assertIn("MASLD_GPU_DISPATCHER_AUTHORITY", shard_text)
        self.assertIn("sequence_gpu_dispatcher_v2", shard_text)
        self.assertIn('${HEAD_IDS:=linear two_layer_mlp}', shard_text)
        self.assertIn('${SCREEN_SEEDS:=1103 1201 1301}', shard_text)
        self.assertIn('${OUTER_FOLDS:=0 1 2 3 4}', shard_text)

    def test_transcriptformer_common_head_queue_is_dispatcher_only(self) -> None:
        value = json.loads(TRANSCRIPTFORMER_QUEUE.read_text(encoding="utf-8"))
        self.assertEqual(value["bundle_id"], "model-training-505")
        self.assertEqual(value["priority"], 52)
        self.assertTrue(value["enabled"])
        self.assertEqual(value["logical_tasks"], 30)
        self.assertEqual(value["family"], "cell_foundation_common_head")
        self.assertEqual(
            value["wrapper_sha256"],
            sha256_file(SHARD_WRAPPER),
        )
        self.assertEqual(
            value["exports"]["MASLD_GPU_DISPATCHER_AUTHORITY"],
            "sequence_gpu_dispatcher_v2",
        )
        self.assertEqual(
            value["exports"]["EMBEDDINGS_ARTIFACTS_SHA256"],
            "fc4601988491b086d5e89f0ff35f2f73d9c27cf74846a6ff7895efc130c825a5",
        )
        self.assertIn(
            "executions/transcriptformer-common-head-queue-validation-v1/ARTIFACTS.json",
            value["required_paths"],
        )


if __name__ == "__main__":
    unittest.main()
