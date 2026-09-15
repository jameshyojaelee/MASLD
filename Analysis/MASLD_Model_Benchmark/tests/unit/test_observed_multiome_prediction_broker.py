from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from masld_bench.adapters.observed_multiome_factorized import FactorizedBaselineError, export, fit, load, predict
from scripts.probe_gse296875_observed_multiome_prediction_broker import (
    PredictionBrokerError,
    _load_npz,
    stage_binding_keys,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_prediction_broker_fixture_20260825.json"


class ObservedMultiomePredictionBrokerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_parent_contracts_pass(self) -> None:
        resolved = validate_config(ROOT, deepcopy(self.config))
        self.assertEqual(len(resolved), 5)

    def test_process_bindings_exclude_forbidden_artifacts(self) -> None:
        self.assertNotIn("evaluator_outcomes", stage_binding_keys("fit"))
        self.assertNotIn("training_labels", stage_binding_keys("predict"))
        self.assertNotIn("fit_state", stage_binding_keys("evaluate"))

    def test_open_firewall_is_rejected_before_parent_reads(self) -> None:
        mutated = deepcopy(self.config)
        mutated["firewall"]["sealed_outcome_read"] = True
        with self.assertRaises(PredictionBrokerError):
            validate_config(ROOT, mutated)

    def test_npz_roster_is_exact(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "values.npz"
            with path.open("xb") as handle:
                np.savez(handle, allowed=np.asarray([1]))
            with self.assertRaises(PredictionBrokerError):
                _load_npz(path, {"different"})

    def test_state_round_trip_keeps_fitted_target_rejection(self) -> None:
        generator = np.random.default_rng(7)
        rows = [f"row-{index}" for index in range(6)]
        targets = [f"target-{index}" for index in range(4)]
        with tempfile.TemporaryDirectory() as directory:
            state = fit("rna_only", row_ids=rows, strata=["a", "b"] * 3, target_ids=targets, rna=generator.normal(size=(6, 2)), observed_atac=generator.poisson(2, size=(6, 2)), sequence_features=generator.normal(size=(4, 2)), target_covariates=generator.normal(size=(4, 2)), target_counts=generator.poisson(2, size=(6, 4)), seed=11)
            path = Path(directory) / "state"
            export(state, path)
            restored = load(path)
            with self.assertRaises(FactorizedBaselineError):
                predict(restored, row_ids=rows, strata=["a", "b"] * 3, target_ids=[targets[0]], rna=generator.normal(size=(6, 2)), observed_atac=generator.poisson(2, size=(6, 2)), sequence_features=generator.normal(size=(1, 2)), target_covariates=generator.normal(size=(1, 2)))


if __name__ == "__main__":
    unittest.main()
