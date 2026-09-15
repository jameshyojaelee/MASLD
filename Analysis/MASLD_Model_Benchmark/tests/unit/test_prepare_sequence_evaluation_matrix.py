from __future__ import annotations

import csv
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

from masld_bench.artifacts import freeze_tree
from scripts.prepare_sequence_evaluation_matrix import (
    SequenceMatrixError,
    _prediction_receipts,
    _prediction_views,
)


class SequenceEvaluationMatrixTests(unittest.TestCase):
    def test_chrombpnet_artifact_expands_to_two_prespecified_views(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "chrombpnet"
            for subdir in ("predictions/full_model", "predictions/nobias_model"):
                path = root / subdir
                path.mkdir(parents=True)
                (path / "summary.json").write_text("{}\n", encoding="utf-8")
            freeze_tree(
                root,
                {
                    "artifact_class": "chrombpnet_full_depth_training_and_predictions",
                    "dataset_id": "gse296875",
                    "model_id": "chrombpnet",
                    "split_id": "donor2_genomic2",
                    "lineage_id": "fibroblast",
                    "seed": 11,
                    "status": "passed",
                    "test_outcomes_used": False,
                    "benchmark_metrics_calculated": False,
                },
            )
            digest = sha256((root / "ARTIFACTS.json").read_bytes()).hexdigest()
            views = _prediction_views(root, digest)
        self.assertEqual(
            {view["candidate_id"] for view in views},
            {"chrombpnet_full__s11", "chrombpnet_nobias__s11"},
        )

    def test_bundle_parser_ignores_failed_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary) / "bundle"
            receipt_dir = bundle / "task_receipts"
            receipt_dir.mkdir(parents=True)
            fields = (
                "status",
                "expected_artifact",
                "artifacts_sha256",
            )
            with (receipt_dir / "failed.tsv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "status": "failed",
                        "expected_artifact": "/not/read",
                        "artifacts_sha256": "not_available",
                    }
                )
            with self.assertRaises(SequenceMatrixError):
                _prediction_receipts(bundle)

    def test_bundle_parser_reads_generic_json_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            artifact = base / "prediction"
            artifact.mkdir()
            bundle = base / "bundle"
            receipt_dir = bundle / "task_receipts"
            receipt_dir.mkdir(parents=True)
            (receipt_dir / "passed.json").write_text(
                json.dumps(
                    {
                        "status": "recovered_passed",
                        "expected_artifact": str(artifact),
                        "artifacts_sha256": "b" * 64,
                    }
                ),
                encoding="utf-8",
            )
            values = _prediction_receipts(bundle)
        self.assertEqual(values, [(artifact, "b" * 64)])


if __name__ == "__main__":
    unittest.main()
