from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts import audit_sequence_regulatory_full_rectangle_completion as gate


ROOT = Path(__file__).resolve().parents[2]


class SequenceRegulatoryFullRectangleCompletionTests(unittest.TestCase):
    def test_contract_separates_native_tasks_and_forbids_partial_ranking(self) -> None:
        contract = gate.load_contract(
            ROOT / "config/evaluation/sequence_regulatory_full_rectangle_completion_20260825.json",
            ROOT,
        )
        self.assertEqual(
            contract["task_families"]["task_native_profile"]["expected_fits"],
            500,
        )
        self.assertEqual(
            contract["task_families"]["scbasset_sequence_only"]["expected_fits"],
            25,
        )
        self.assertFalse(contract["gate"]["partial_ranking_allowed"])
        self.assertFalse(contract["claims"]["cross_task_family_ranking_allowed"])

    def test_scbasset_expected_rectangle_is_exactly_five_by_five(self) -> None:
        keys = [(fold, seed) for fold in gate.FOLDS for seed in gate.SEEDS]
        self.assertEqual(len(keys), 25)
        self.assertEqual(len(set(keys)), 25)

    def test_combined_gate_refuses_partial_and_terminal_rectangles(self) -> None:
        self.assertEqual(
            gate.combined_status(
                "ready_full_rectangle_locked", "waiting_for_full_rectangle"
            ),
            "waiting_for_full_rectangles",
        )
        self.assertEqual(
            gate.combined_status(
                "ready_full_rectangle_locked",
                "blocked_terminal_failure_requires_new_immutable_retry",
            ),
            "blocked_rectangle_requires_new_immutable_retry",
        )
        self.assertEqual(
            gate.combined_status(
                "ready_full_rectangle_locked", "ready_full_rectangle_locked"
            ),
            "ready_both_rectangles_locked",
        )

    def test_source_manifest_omits_partial_prediction_rows(self) -> None:
        contract = json.loads(
            (
                ROOT
                / "config/evaluation/sequence_regulatory_full_rectangle_completion_20260825.json"
            ).read_text(encoding="utf-8")
        )
        fit_rows = [
            {
                "lineage_id": "hepatocyte",
                "split_id": "donor0_genomic0",
                "disposition": "not_dispatched_or_pending",
                "canonical_artifact_path": "not_available",
                "canonical_artifacts_sha256": "not_available",
            }
        ]
        with tempfile.TemporaryDirectory(dir=ROOT / "executions") as temporary:
            output = Path(temporary) / "sources.tsv"
            with patch.object(
                gate,
                "_verify_bound_root",
                side_effect=lambda root, value, expected: (root / value).resolve(),
            ):
                gate.write_sequence_source_manifest(
                    root=ROOT,
                    contract=contract,
                    fit_rows=fit_rows,
                    output=output,
                )
            lines = output.read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(lines), 6)
        self.assertTrue(all("prediction_artifact" not in line for line in lines[1:]))


if __name__ == "__main__":
    unittest.main()
