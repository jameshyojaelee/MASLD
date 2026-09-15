from __future__ import annotations

import csv
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import tempfile
import unittest

from scripts.evaluate_gse281364_sei_mpralegnet import (
    CONTEXTS,
    GSE281364EvaluationError,
    OUTCOME_FIELDS,
    _matched_endpoint_rows,
    _load_manifest,
    _outcome_deltas,
    activity_delta,
    correlations,
    group_bootstrap_interval,
)


class GSE281364DevelopmentEvaluatorTests(unittest.TestCase):
    def test_frozen_tree_verifier_rejects_tampered_member(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            payload = root / "payload.txt"
            payload.write_text("original\n", encoding="utf-8")
            manifest = {
                "artifacts": [
                    {
                        "path": "payload.txt",
                        "sha256": sha256(payload.read_bytes()).hexdigest(),
                        "size_bytes": payload.stat().st_size,
                    }
                ],
                "metadata": {"fixture": True},
                "schema_version": "masld-bench-artifacts-v1",
            }
            manifest_path = root / "ARTIFACTS.json"
            manifest_path.write_text(
                json.dumps(manifest, separators=(",", ":"), sort_keys=True) + "\n",
                encoding="utf-8",
            )
            complete = {
                "artifact_count": 1,
                "manifest_sha256": sha256(manifest_path.read_bytes()).hexdigest(),
                "schema_version": "masld-bench-complete-v1",
            }
            (root / "COMPLETE").write_text(
                json.dumps(complete, separators=(",", ":"), sort_keys=True) + "\n",
                encoding="utf-8",
            )
            expected = sha256(manifest_path.read_bytes()).hexdigest()
            self.assertEqual(_load_manifest(root, expected)["metadata"], {"fixture": True})
            payload.write_text("tampered\n", encoding="utf-8")
            with self.assertRaises(GSE281364EvaluationError):
                _load_manifest(root, expected)

    def test_activity_delta_is_alt_minus_ref_log_ratio(self) -> None:
        observed = activity_delta(9, 19, 9, 39, pseudocount=1.0)
        self.assertAlmostEqual(observed, 1.0)
        self.assertAlmostEqual(activity_delta(10, 10, 10, 10), 0.0)

    def test_correlations_are_signed_and_exact(self) -> None:
        pearson, spearman = correlations([1.0, 2.0, 3.0], [-1.0, -2.0, -3.0])
        self.assertAlmostEqual(pearson, -1.0)
        self.assertAlmostEqual(spearman, -1.0)
        with self.assertRaises(GSE281364EvaluationError):
            correlations([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])

    def test_group_bootstrap_is_deterministic(self) -> None:
        observed = [float(value) for value in range(1, 21)]
        predicted = [value + 0.1 * (-1) ** index for index, value in enumerate(observed)]
        groups = [f"group_{index // 2}" for index in range(20)]
        first = group_bootstrap_interval(
            observed, predicted, groups, replicates=100, seed=17
        )
        second = group_bootstrap_interval(
            observed, predicted, groups, replicates=100, seed=17
        )
        self.assertEqual(first, second)
        self.assertGreaterEqual(first["valid_bootstrap_replicates"], 95)

    def test_outcome_loader_requires_four_complete_replicate_pairs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "outcomes.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=list(OUTCOME_FIELDS),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for context in CONTEXTS:
                    for replicate in range(1, 5):
                        for allele in ("ref", "alt"):
                            writer.writerow(
                                {
                                    "element_id": "element_1",
                                    "allele": allele,
                                    "context_id": context,
                                    "cell_line": context.split("_", 1)[0],
                                    "condition": context.split("_", 1)[1],
                                    "experimental_replicate": replicate,
                                    "sample_id": f"{context}_r{replicate}",
                                    "DNA": 10,
                                    "RNA": 10 if allele == "ref" else 20,
                                    "n_barcodes": 5,
                                    "assay_state": "observed",
                                    "missing_reason": "not_applicable",
                                    "pairing": "same_sample_different_aliquot",
                                    "biological_unit": "experimental_replicate",
                                    "donor_id": "not_applicable",
                                }
                            )
            outcomes, coverage = _outcome_deltas(path, {"element_1", "element_2"})
            self.assertEqual(coverage, {context: 1 for context in CONTEXTS})
            self.assertEqual(len(outcomes), 4)
            self.assertTrue(
                all(value["replicates"] == 4 for value in outcomes.values())
            )

    def test_native_endpoints_keep_sei_magnitude_and_mpralegnet_sign(self) -> None:
        prediction = [
            {
                "element_id": "element_1",
                "outer_locus_sequence_group_id": "group_1",
                "outer_fold": 0,
                "prediction_score": -2.0,
            }
        ]
        outcomes = {
            ("element_1", "HepG2_control"): {
                "mean_signed_delta": -3.0,
                "sd_signed_delta": 0.25,
                "replicates": 4,
            }
        }
        sei = _matched_endpoint_rows(
            prediction,
            outcomes,
            model_id="sei",
            endpoint_id="magnitude",
            contexts=("HepG2_control",),
            magnitude=True,
            comparability="magnitude",
        )
        mpralegnet = _matched_endpoint_rows(
            prediction,
            outcomes,
            model_id="mpralegnet",
            endpoint_id="signed",
            contexts=("HepG2_control",),
            magnitude=False,
            comparability="signed",
        )
        self.assertEqual(sei[0]["observed_score"], 3.0)
        self.assertEqual(mpralegnet[0]["observed_score"], -3.0)
        self.assertFalse(math.isclose(sei[0]["observed_score"], mpralegnet[0]["observed_score"]))


if __name__ == "__main__":
    unittest.main()
