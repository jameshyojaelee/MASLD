from __future__ import annotations

import unittest

from scripts.correct_geneformer_weight_state import (
    DECLARED_LICENCE,
    GeneformerCorrectionError,
    correct_bundle,
    correct_pinner,
    normalise_licence,
)


def bundle_document(**overrides: object) -> dict:
    document = {
        "weight_content_downloaded": False,
        "license": DECLARED_LICENCE,
        "artifacts": [{"path": "Geneformer-V1-10M/model.safetensors", "sha256": "a" * 64}],
        "code_source_artifacts": [{"path": "setup.py", "sha256": "b" * 64}],
    }
    document.update(overrides)
    return document


def pinner_document() -> dict:
    return {
        "bound_evidence": [
            {"path": "config/artifacts/models/geneformer/checkpoints.json", "sha256": "old"},
            {"path": "executions/other/ARTIFACTS.json", "sha256": "keep"},
        ],
        "recorded_inconsistencies": [
            {
                "action": "recorded_only_not_corrected_by_this_overlay",
                "field": "config/artifacts/models/geneformer/checkpoints.json:weight_content_downloaded",
                "recorded_value": False,
            }
        ],
    }


class LicenceNormalisationTests(unittest.TestCase):
    def test_singular_key_splits_into_the_standard_pair(self) -> None:
        doc = bundle_document()
        normalise_licence(doc)
        self.assertEqual(doc["weight_license"], DECLARED_LICENCE)
        self.assertEqual(doc["code_license"], DECLARED_LICENCE)
        self.assertNotIn("license", doc)

    def test_the_declared_value_is_preserved_exactly(self) -> None:
        """The whole point: do not deny terms that were granted."""

        doc = bundle_document()
        correct_bundle(doc)
        self.assertEqual(doc["weight_license"], "Apache-2.0")
        self.assertNotEqual(doc["weight_license"], "UNDECLARED")

    def test_a_correction_that_would_downgrade_the_licence_is_refused(self) -> None:
        doc = bundle_document(license="UNDECLARED")
        with self.assertRaises(GeneformerCorrectionError) as caught:
            correct_bundle(doc)
        self.assertIn("not the audited", str(caught.exception))

    def test_already_normalised_bundle_is_refused(self) -> None:
        doc = bundle_document(weight_license=DECLARED_LICENCE)
        with self.assertRaises(GeneformerCorrectionError):
            normalise_licence(doc)


class CorrectBundleTests(unittest.TestCase):
    def test_stale_flag_is_flipped(self) -> None:
        doc = bundle_document()
        changes = correct_bundle(doc)
        self.assertIs(doc["weight_content_downloaded"], True)
        self.assertEqual(changes[0]["field"], "weight_content_downloaded")

    def test_already_corrected_bundle_is_refused(self) -> None:
        doc = bundle_document(weight_content_downloaded=True)
        with self.assertRaises(GeneformerCorrectionError):
            correct_bundle(doc)

    def test_recorded_digests_are_not_touched(self) -> None:
        """Unlike scGPT, geneformer's digests were already resolved."""

        doc = bundle_document()
        before = [a["sha256"] for a in doc["artifacts"]]
        correct_bundle(doc)
        self.assertEqual([a["sha256"] for a in doc["artifacts"]], before)


class CorrectPinnerTests(unittest.TestCase):
    def test_only_the_bundle_pin_moves(self) -> None:
        doc = pinner_document()
        correct_pinner(doc, "n" * 64)
        by_path = {e["path"]: e["sha256"] for e in doc["bound_evidence"]}
        self.assertEqual(
            by_path["config/artifacts/models/geneformer/checkpoints.json"], "n" * 64
        )
        self.assertEqual(by_path["executions/other/ARTIFACTS.json"], "keep")

    def test_the_deferred_item_is_marked_closed_not_deleted(self) -> None:
        doc = pinner_document()
        correct_pinner(doc, "n" * 64)
        record = doc["recorded_inconsistencies"][0]
        self.assertEqual(record["action"], "corrected_2026_08_25")
        self.assertIn("deferred item is now closed", record["correction"])
        self.assertIs(record["recorded_value"], False)

    def test_a_pinner_not_referencing_the_bundle_aborts(self) -> None:
        with self.assertRaises(GeneformerCorrectionError):
            correct_pinner({"bound_evidence": [{"path": "x", "sha256": "y"}]}, "n" * 64)


if __name__ == "__main__":
    unittest.main()
