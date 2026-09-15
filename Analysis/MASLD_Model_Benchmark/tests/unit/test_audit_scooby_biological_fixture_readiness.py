#!/usr/bin/env python3
"""Unit tests for the released-Scooby biological-fixture readiness audit."""

from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts/audit_scooby_biological_fixture_readiness.py"
SPEC = importlib.util.spec_from_file_location("scooby_biological_readiness", MODULE)
assert SPEC and SPEC.loader
readiness = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(readiness)


class ScoobyBiologicalFixtureReadinessTest(unittest.TestCase):
    def test_pinned_terms_require_exact_revision(self) -> None:
        metadata = {
            "sha": "abc",
            "cardData": {"license": "mit"},
            "siblings": [{"rfilename": "README.md"}],
        }
        self.assertEqual(readiness.hf_terms(metadata, "abc")["declared_license"], "mit")
        with self.assertRaises(readiness.ScoobyBiologicalReadinessError):
            readiness.hf_terms(metadata, "def")

    def test_adjacent_license_file_does_not_create_declared_terms(self) -> None:
        metadata = {
            "sha": "abc",
            "cardData": {},
            "siblings": [{"rfilename": "README.md"}],
        }
        result = readiness.hf_terms(metadata, "abc")
        self.assertFalse(result["terms_declared"])
        self.assertIsNone(result["declared_license"])

    def test_sbatch_parser_keeps_generic_nslab_l40s_identity(self) -> None:
        source = "\n".join(
            (
                "#SBATCH --job-name=model-probe-284",
                "#SBATCH --partition=gpu",
                "#SBATCH --qos=nslab",
                "#SBATCH --gres=gpu:l40s:1",
            )
        )
        header = readiness.parse_sbatch_header(source)
        self.assertEqual(header["job-name"], "model-probe-284")
        self.assertEqual(header["qos"], "nslab")
        self.assertEqual(header["gres"], "gpu:l40s:1")

    def test_biological_decisions_never_treat_raw_atac_as_context(self) -> None:
        decisions = readiness.biological_decisions()
        self.assertEqual(
            decisions["released_checkpoint_on_gse281367"]["status"],
            "not_applicable",
        )
        self.assertFalse(
            decisions["released_checkpoint_on_gse281367"][
                "biological_execution_authorized"
            ]
        )

    def test_new_liver_latent_is_an_adaptation_not_zero_shot(self) -> None:
        decisions = copy.deepcopy(readiness.biological_decisions())
        adaptation = decisions["project_fitted_scooby_adaptation_on_gse296875"]
        self.assertEqual(
            adaptation["identity"],
            "new_project_adaptation_not_released_checkpoint_inference",
        )
        self.assertFalse(adaptation["biological_execution_authorized"])


if __name__ == "__main__":
    unittest.main()
