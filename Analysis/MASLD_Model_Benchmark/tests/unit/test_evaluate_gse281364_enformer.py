from __future__ import annotations

import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

from scripts.evaluate_gse281364_enformer import (
    CONTEXTS,
    EXCLUDED_SCORES,
    SCORES,
    EnformerEvaluationError,
    OUTCOME_FIELDS,
    PREDICTION_FIELDS,
    activity_delta,
    correlations,
    group_bootstrap,
    load_outcomes,
    read_predictions,
    verify_frozen_tree,
)


class EnformerDevelopmentEvaluationTests(unittest.TestCase):
    def test_frozen_tree_verifier_rejects_tamper(self) -> None:
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
            digest = sha256(manifest_path.read_bytes()).hexdigest()
            (root / "COMPLETE").write_text(
                json.dumps(
                    {
                        "artifact_count": 1,
                        "manifest_sha256": digest,
                        "schema_version": "masld-bench-complete-v1",
                    },
                    separators=(",", ":"),
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
            self.assertEqual(verify_frozen_tree(root, digest)["metadata"], {"fixture": True})
            payload.write_text("tampered\n", encoding="utf-8")
            with self.assertRaises(EnformerEvaluationError):
                verify_frozen_tree(root, digest)

    def test_signed_activity_delta(self) -> None:
        self.assertAlmostEqual(activity_delta(9, 19, 9, 39, pseudocount=1.0), 1.0)
        self.assertAlmostEqual(activity_delta(10, 10, 10, 10), 0.0)

    def test_only_hepg2_contexts_and_specific_scores_are_evaluated(self) -> None:
        self.assertEqual(CONTEXTS, ("HepG2_control", "HepG2_PAOA"))
        self.assertEqual(
            set(SCORES),
            {
                "hepg2_accessibility_sad",
                "hepg2_accessibility_sar",
                "liver_accessibility_sad",
                "liver_accessibility_sar",
            },
        )
        self.assertEqual(
            set(EXCLUDED_SCORES),
            {"all_accessibility_mean_sad", "all_accessibility_mean_sar"},
        )

    def test_correlations_remain_signed(self) -> None:
        pearson, spearman = correlations([1.0, 2.0, 3.0], [-1.0, -2.0, -3.0])
        self.assertAlmostEqual(pearson, -1.0)
        self.assertAlmostEqual(spearman, -1.0)

    def test_prediction_receipt_rejects_champion_or_native_parity_claim(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prediction_root = root / "predictions"
            prediction_root.mkdir()
            path = prediction_root / "predictions.tsv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=list(PREDICTION_FIELDS), delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                for index in range(1033):
                    writer.writerow(
                        {
                            "fixture_id": f"fixture_{index}",
                            "element_id": f"element_{index}",
                            "outer_locus_sequence_group_id": f"group_{index}",
                            "outer_fold": index % 5,
                            "hepg2_accessibility_sad": index + 0.1,
                            "hepg2_accessibility_sar": index + 0.2,
                            "liver_accessibility_sad": index + 0.3,
                            "liver_accessibility_sar": index + 0.4,
                            "all_accessibility_mean_sad": index + 0.5,
                            "all_accessibility_mean_sar": index + 0.6,
                        }
                    )
            receipt = {
                "status": "pass_outcome_blind_restricted_prediction",
                "model_id": "enformer_crested_restricted_port",
                "registered_scientific_identity": "restricted_conversion_not_native_sonnet",
                "elements": 1033,
                "outer_locus_sequence_groups": 1033,
                "outer_folds": 5,
                "tracks": 5313,
                "accessibility_tracks": 684,
                "allele_delta_sign": "ALT_minus_REF",
                "deterministic_repeat_max_abs": 0.0,
                "outcomes_read": False,
                "reporter_counts_read": False,
                "sealed_outcomes_read": False,
                "model_fitted_or_adapted": False,
                "native_sonnet_parity_established": False,
                "champion_eligible": False,
                "terms": "internal_academic_noncommercial_nontransferable_restricted_comparator",
                "predictions_sha256": sha256(path.read_bytes()).hexdigest(),
            }
            receipt_path = prediction_root / "receipt.json"
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            self.assertEqual(len(read_predictions(root)), 1033)
            receipt["champion_eligible"] = True
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaises(EnformerEvaluationError):
                read_predictions(root)

    def test_outcomes_require_four_complete_repeat_pairs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "outcomes.tsv.gz"
            rows = []
            for context in CONTEXTS:
                for replicate in range(1, 5):
                    for allele in ("ref", "alt"):
                        rows.append(
                            {
                                "element_id": "element_1",
                                "allele": allele,
                                "context_id": context,
                                "cell_line": "HepG2",
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
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=list(OUTCOME_FIELDS), delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)
            self.assertEqual(len(load_outcomes(path, {"element_1"})), 2)
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=list(OUTCOME_FIELDS), delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows[:-1])
            with self.assertRaises(EnformerEvaluationError):
                load_outcomes(path, {"element_1"})

    def test_group_bootstrap_is_deterministic(self) -> None:
        observed = [float(value) for value in range(20)]
        predicted = [value + (-1) ** value * 0.1 for value in observed]
        groups = [f"group_{value}" for value in range(20)]
        first = group_bootstrap(observed, predicted, groups, seed=17, replicates=100)
        second = group_bootstrap(observed, predicted, groups, seed=17, replicates=100)
        self.assertEqual(first, second)
        self.assertGreaterEqual(first["valid_bootstrap_replicates"], 95)


if __name__ == "__main__":
    unittest.main()
