#!/usr/bin/env python3
"""Unit checks for outcome-blind elastic-net convergence diagnostics."""

from __future__ import annotations

import unittest

from scripts.diagnose_elastic_net_convergence_study_50000 import (
    ElasticDiagnosticError,
    select_configs,
)


class ElasticNetDiagnosticTests(unittest.TestCase):
    def test_default_preserves_original_five_configurations(self) -> None:
        selected = select_configs(None)
        self.assertEqual(len(selected), 5)
        self.assertNotIn(
            "float64_tol1e-3_c1_max50000",
            {str(config["config_id"]) for config in selected},
        )

    def test_requested_order_is_stable(self) -> None:
        selected = select_configs(
            ["float64_tol1e-3_c1_max50000", "float64_tol1e-2_c1"]
        )
        self.assertEqual(
            [config["max_iter"] for config in selected], [50_000, 10_000]
        )

    def test_unknown_or_duplicate_configuration_is_rejected(self) -> None:
        with self.assertRaises(ElasticDiagnosticError):
            select_configs(["missing"])
        with self.assertRaises(ElasticDiagnosticError):
            select_configs(["float64_tol1e-2_c1", "float64_tol1e-2_c1"])


if __name__ == "__main__":
    unittest.main()
