from __future__ import annotations

import json
from pathlib import Path
import unittest

from scripts.register_gse244832_observed_atac_execution import (
    ALLOWED_MODELS,
    ALLOWED_OUTPUT_FAMILIES,
    GSE244832ExecutionRegistrationError,
    _execution_family,
    _read_axis_states,
    _read_axis_windows,
    _read_query_windows,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = (
    ROOT
    / "config/evaluation/gse244832_observed_atac_execution_registration_20260825.json"
)
AXIS = ROOT / "executions/gse244832-reference-guarded-atac-exchange-axis-20260825"


class GSE244832ObservedATACExecutionRegistrationTests(unittest.TestCase):
    def test_real_config_binds_only_four_coordinate_independent_families(self) -> None:
        config = json.loads(CONFIG.read_text(encoding="utf-8"))
        validate_config(ROOT, config)
        self.assertEqual(tuple(config["allowed_model_ids"]), ALLOWED_MODELS)
        self.assertFalse(config["sequence_or_reference_features_authorized"])
        self.assertFalse(config["direct_query_artifact_execution_authorized"])
        self.assertEqual(
            ALLOWED_OUTPUT_FAMILIES,
            {model_id: ["masked_accessibility_count"] for model_id in ALLOWED_MODELS},
        )

    def test_real_axis_has_eighteen_donors_and_explicit_missingness(self) -> None:
        rows, states = _read_axis_states(AXIS)
        self.assertEqual(len(rows), 18 * 4)
        self.assertEqual(states.shape, (18, 4))
        self.assertEqual(int((states == 0).sum()), 48)
        self.assertEqual(int((states == 1).sum()), 4)
        self.assertEqual(int((states == 3).sum()), 20)

    def test_model_specific_family_guard_blocks_reference_dependent_models(self) -> None:
        rows, blocked = _execution_family(AXIS)
        allowed = {
            row["model_id"]
            for row in rows
            if row["execution_allowed"] == "true"
        }
        self.assertEqual(allowed, set(ALLOWED_MODELS))
        for model_id in (
            "epibert",
            "epcotv2",
            "get",
            "epiagent",
            "scbasset_observed_cell_embedding",
        ):
            self.assertIn(model_id, blocked)

    def test_query_window_axes_are_role_exact_and_whole_contig_disjoint(self) -> None:
        axis = _read_axis_windows(AXIS)
        valid = _read_query_windows(
            ROOT
            / "executions/gse244832-label-free-atac-query-valid-context-20260825/query_windows.tsv"
        )
        test = _read_query_windows(
            ROOT
            / "executions/gse244832-label-free-atac-query-test-context-20260825/query_windows.tsv"
        )
        self.assertEqual(valid, axis["valid"])
        self.assertEqual(test, axis["test"])
        self.assertFalse({row[1] for row in valid} & {row[1] for row in test})

    def test_unknown_or_incomplete_axis_fails_closed(self) -> None:
        with self.assertRaises((FileNotFoundError, GSE244832ExecutionRegistrationError)):
            _read_axis_states(ROOT / "executions/not-a-real-axis")


if __name__ == "__main__":
    unittest.main()
