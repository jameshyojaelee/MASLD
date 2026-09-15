from __future__ import annotations

import copy
import json
import unittest

from scripts.correct_scgpt_weight_state import (
    CORRECTED_BLOCKER_TEXT,
    STALE_BLOCKER_TEXT,
    STALE_EVIDENCE,
    ScgptCorrectionError,
    correct_bundle,
    correct_pinner,
)


RESOLVED = {
    "scgpt_continual": "a" * 64,
    "scgpt_whole_human": "b" * 64,
}


def bundle_document() -> dict:
    return {
        "weight_content_downloaded": False,
        "weight_license": "UNDECLARED",
        "source_evidence": "Queried on 2026-08-22. " + STALE_EVIDENCE,
        "bundles": {
            name: {
                "artifacts": [
                    {"path": "best_model.pt", "sha256": "UNRESOLVED", "size_bytes": 1},
                    {"path": "args.json", "sha256": "c" * 64, "size_bytes": 2},
                ]
            }
            for name in RESOLVED
        },
    }


def pinner_document() -> dict:
    return {
        "bound_evidence": [
            {"path": "config/artifacts/models/scgpt/checkpoints.json", "sha256": "old"},
            {"path": "executions/other/ARTIFACTS.json", "sha256": "keep"},
        ],
        "family_specific_blockers": {
            "detail": "No runtime_contract block. " + STALE_BLOCKER_TEXT,
            "runtime_contract_block_present": False,
        },
    }


class CorrectBundleTests(unittest.TestCase):
    def test_flag_and_both_digests_are_corrected(self) -> None:
        doc = bundle_document()
        changes = correct_bundle(doc, RESOLVED)
        self.assertIs(doc["weight_content_downloaded"], True)
        for name, expected in RESOLVED.items():
            artifacts = {a["path"]: a for a in doc["bundles"][name]["artifacts"]}
            self.assertEqual(artifacts["best_model.pt"]["sha256"], expected)
            self.assertEqual(artifacts["args.json"]["sha256"], "c" * 64)
        self.assertEqual(len(changes), 4)

    def test_licence_is_never_upgraded_by_the_correction(self) -> None:
        """Recording the truth must not imply terms nobody granted."""

        doc = bundle_document()
        correct_bundle(doc, RESOLVED)
        self.assertEqual(doc["weight_license"], "UNDECLARED")

    def test_stale_prose_is_replaced_so_the_record_does_not_still_deny_it(self) -> None:
        doc = bundle_document()
        correct_bundle(doc, RESOLVED)
        self.assertNotIn(STALE_EVIDENCE, doc["source_evidence"])
        self.assertIn("acquired on 2026-08-23", doc["source_evidence"])
        self.assertIn("redistribution stays prohibited", doc["source_evidence"])

    def test_already_corrected_bundle_is_refused(self) -> None:
        doc = bundle_document()
        doc["weight_content_downloaded"] = True
        with self.assertRaises(ScgptCorrectionError):
            correct_bundle(doc, RESOLVED)

    def test_a_resolved_digest_is_not_overwritten(self) -> None:
        doc = bundle_document()
        doc["bundles"]["scgpt_continual"]["artifacts"][0]["sha256"] = "d" * 64
        with self.assertRaises(ScgptCorrectionError):
            correct_bundle(doc, RESOLVED)

    def test_missing_stale_evidence_aborts(self) -> None:
        doc = bundle_document()
        doc["source_evidence"] = "something else entirely"
        with self.assertRaises(ScgptCorrectionError):
            correct_bundle(doc, RESOLVED)


class CorrectPinnerTests(unittest.TestCase):
    def test_only_the_bundle_pin_moves(self) -> None:
        doc = pinner_document()
        correct_pinner(doc, "n" * 64)
        by_path = {e["path"]: e["sha256"] for e in doc["bound_evidence"]}
        self.assertEqual(by_path["config/artifacts/models/scgpt/checkpoints.json"], "n" * 64)
        self.assertEqual(by_path["executions/other/ARTIFACTS.json"], "keep")

    def test_the_prose_claim_is_corrected_too(self) -> None:
        doc = pinner_document()
        correct_pinner(doc, "n" * 64)
        detail = doc["family_specific_blockers"]["detail"]
        self.assertNotIn(STALE_BLOCKER_TEXT, detail)
        self.assertIn(CORRECTED_BLOCKER_TEXT, detail)
        self.assertIn("UNDECLARED", detail)

    def test_a_pinner_not_referencing_the_bundle_aborts(self) -> None:
        doc = {"bound_evidence": [{"path": "other", "sha256": "x"}],
               "family_specific_blockers": {"detail": "unrelated"}}
        with self.assertRaises(ScgptCorrectionError):
            correct_pinner(doc, "n" * 64)


if __name__ == "__main__":
    unittest.main()
