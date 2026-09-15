#!/usr/bin/env python3
"""Focused tests for outcome-blind common embedding preparation."""

from __future__ import annotations

import ast
import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.hashing import sha256_file
from scripts import prepare_common_cell_embeddings_50000 as prepare


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts/prepare_common_cell_embeddings_50000.py"
MODEL_ID = "geneformer_v2_316m"
CHECKPOINT_SHA256 = "a" * 64


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["row_id", "donor_id", "dataset", "outer_fold"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def make_fixture(parent: Path) -> dict[str, object]:
    rows = [
        {
            "row_id": f"r{index}",
            "donor_id": f"d{index // 2}",
            "dataset": f"s{index // 2}",
            "outer_fold": index // 2,
        }
        for index in range(10)
    ]
    split = parent / "split"
    split.mkdir()
    write_tsv(split / "row_outer_folds.tsv", rows)
    freeze_tree(
        split,
        {
            "artifact_class": "cell_state_study_outer_split",
            "dataset_view_id": prepare.DATASET_VIEW_ID,
            "split_id": prepare.SPLIT_ID,
            "rows": 10,
            "donors": 5,
            "studies": 5,
            "target_labels_used_for_assignment": False,
            "sealed_outcomes_read": False,
            "histology_read": False,
            "status": "passed",
        },
    )

    raw = parent / "raw"
    (raw / "extraction").mkdir(parents=True)
    embeddings = np.arange(30, dtype=np.float32).reshape(10, 3)
    np.save(raw / "extraction/embeddings.npy", embeddings, allow_pickle=False)
    (raw / "extraction/embedding_row_order.txt").write_text(
        "".join(f"r{index}\n" for index in range(10)), encoding="utf-8"
    )
    embedding_sha = sha256_file(raw / "extraction/embeddings.npy")
    row_sha = sha256_file(raw / "extraction/embedding_row_order.txt")
    (raw / "extraction/extraction_receipt.json").write_text(
        json.dumps(
            {
                "model_id": MODEL_ID,
                "dataset_view_id": prepare.DATASET_VIEW_ID,
                "rows": 10,
                "embedding_shape": [10, 3],
                "embedding_dtype": "float32",
                "embedding_sha256": embedding_sha,
                "row_order_sha256": row_sha,
                "fixture_row_order_sha256": row_sha,
                "embedding_policy": "fixture_mean_pool",
                "raw_unrefined_embeddings_only": True,
                "downstream_head_fit": False,
                "evaluation_labels_read": False,
                "histology_read": False,
                "sealed_outcomes_read": False,
                "status": "pass_raw_outcome_blind_extraction",
            }
        ),
        encoding="utf-8",
    )
    freeze_tree(
        raw,
        {
            "artifact_class": "geneformer_v2_316m_frozen_screen_raw_embeddings",
            "model_id": MODEL_ID,
            "dataset_view_id": prepare.DATASET_VIEW_ID,
            "rows": 10,
            "raw_unrefined_embeddings_only": True,
            "downstream_head_fit": False,
            "evaluation_labels_read": False,
            "histology_read": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )

    registry = parent / "checkpoints.json"
    registry.write_text(
        json.dumps(
            {
                "architecture_contract": {MODEL_ID: {"hidden_size": 3}},
                "artifacts": [
                    {
                        "path": "Geneformer-V2-316M/model.safetensors",
                        "sha256": CHECKPOINT_SHA256,
                        "training_cutoff": "2024-12",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    exposure = parent / "exposure.json"
    exposure.write_text(
        json.dumps(
            {
                "checkpoint_findings": {
                    MODEL_ID: {"exposure_state": "target_label_unexposed"}
                },
                "sealed_champion_eligibility": (
                    "eligible_on_exposure_only_other_admission_gates_still_apply"
                ),
            }
        ),
        encoding="utf-8",
    )
    return {
        "rows": rows,
        "raw": raw,
        "split": split,
        "registry": registry,
        "exposure": exposure,
    }


def run_fixture(fixture: dict[str, object], output: Path) -> None:
    raw = fixture["raw"]
    split = fixture["split"]
    registry = fixture["registry"]
    exposure = fixture["exposure"]
    assert isinstance(raw, Path)
    assert isinstance(split, Path)
    assert isinstance(registry, Path)
    assert isinstance(exposure, Path)
    prepare.run(
        raw_root=raw,
        raw_embeddings_relative=Path("extraction/embeddings.npy"),
        raw_row_order_relative=Path("extraction/embedding_row_order.txt"),
        raw_receipt_relative=Path("extraction/extraction_receipt.json"),
        split_root=split,
        checkpoint_registry=registry,
        exposure_audit=exposure,
        output=output,
        expected_raw_artifacts_sha256=sha256_file(raw / "ARTIFACTS.json"),
        expected_split_artifacts_sha256=sha256_file(split / "ARTIFACTS.json"),
        expected_checkpoint_registry_sha256=sha256_file(registry),
        expected_exposure_audit_sha256=sha256_file(exposure),
        expected_model_id=MODEL_ID,
        expected_checkpoint_artifact_path="Geneformer-V2-316M/model.safetensors",
        expected_checkpoint_sha256=CHECKPOINT_SHA256,
        expected_exposure_status="target_label_unexposed",
        expected_sealed_champion_eligible=True,
        expected_rows=10,
        expected_width=3,
        expected_donors=5,
        expected_studies=5,
    )


class PrepareCommonCellEmbeddingsTests(unittest.TestCase):
    def test_prepares_exact_outcome_blind_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(Path(directory))
            output = Path(directory) / "prepared"
            run_fixture(fixture, output)
            manifest = verify_frozen_tree(output)
            self.assertEqual(manifest["metadata"]["model_id"], MODEL_ID)
            self.assertTrue(manifest["metadata"]["sealed_champion_eligible"])
            self.assertFalse(manifest["metadata"]["evaluation_labels_read"])
            receipt = json.loads(
                (output / "embeddings/receipt.json").read_text(encoding="utf-8")
            )
            self.assertEqual(receipt["checkpoint_sha256"], CHECKPOINT_SHA256)
            self.assertEqual(receipt["evaluation_label_columns_read"], [])
            self.assertFalse(receipt["source_label_table_read"])
            with np.load(
                output / "embeddings/common_embeddings.npz", allow_pickle=False
            ) as bundle:
                self.assertEqual(
                    set(bundle.files), {"embeddings", "outer_folds", "row_ids"}
                )
                self.assertEqual(bundle["embeddings"].shape, (10, 3))
                self.assertEqual(bundle["outer_folds"].tolist(), [0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
                self.assertEqual(bundle["row_ids"].tolist(), [f"r{i}" for i in range(10)])

    def test_rejects_row_order_not_matching_frozen_split(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            fixture = make_fixture(parent)
            raw = fixture["raw"]
            assert isinstance(raw, Path)
            # Build a new immutable raw fixture with a valid internal receipt but
            # an order that is inconsistent with the independently frozen split.
            second = parent / "raw_permuted"
            (second / "extraction").mkdir(parents=True)
            np.save(
                second / "extraction/embeddings.npy",
                np.arange(30, dtype=np.float32).reshape(10, 3),
                allow_pickle=False,
            )
            order = ["r1", "r0", *[f"r{i}" for i in range(2, 10)]]
            (second / "extraction/embedding_row_order.txt").write_text(
                "".join(f"{value}\n" for value in order), encoding="utf-8"
            )
            receipt = json.loads(
                (raw / "extraction/extraction_receipt.json").read_text(encoding="utf-8")
            )
            receipt["embedding_sha256"] = sha256_file(second / "extraction/embeddings.npy")
            receipt["row_order_sha256"] = sha256_file(
                second / "extraction/embedding_row_order.txt"
            )
            receipt["fixture_row_order_sha256"] = receipt["row_order_sha256"]
            (second / "extraction/extraction_receipt.json").write_text(
                json.dumps(receipt), encoding="utf-8"
            )
            freeze_tree(second, verify_frozen_tree(raw)["metadata"])
            fixture["raw"] = second
            with self.assertRaisesRegex(
                prepare.CommonEmbeddingPreparationError, "row order"
            ):
                run_fixture(fixture, parent / "rejected")

    def test_rejects_exposure_spoof(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            fixture = make_fixture(parent)
            exposure = fixture["exposure"]
            assert isinstance(exposure, Path)
            value = json.loads(exposure.read_text(encoding="utf-8"))
            value["checkpoint_findings"][MODEL_ID]["exposure_state"] = "unknown"
            exposure.write_text(json.dumps(value), encoding="utf-8")
            with self.assertRaisesRegex(
                prepare.CommonEmbeddingPreparationError, "exposure"
            ):
                run_fixture(fixture, parent / "rejected")

    def test_preparation_has_no_source_label_input(self) -> None:
        text = SCRIPT.read_text(encoding="utf-8")
        tree = ast.parse(text)
        run_function = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "run"
        )
        arguments = {argument.arg for argument in run_function.args.kwonlyargs}
        self.assertNotIn("source", arguments)
        self.assertNotIn("selection", arguments)
        self.assertNotIn("broad_label", text)
        self.assertIn('"source_label_table_read": False', text)


if __name__ == "__main__":
    unittest.main()
