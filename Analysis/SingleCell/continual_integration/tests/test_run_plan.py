from __future__ import annotations

import unittest
import json
import tempfile
from pathlib import Path
from unittest.mock import patch

from masld_cl.config import ConfigError, load_config
from masld_cl.run_plan import RunPlanError, create_run_plan, screen_settings


class TestRunPlan(unittest.TestCase):
    def setUp(self):
        self.config = {"screen": {
            "lambda_values": [0, 0.1, 0.5, 1, 10, 100],
            "replay_values": [0, 0.1, 0.2, 0.4], "seed": 17,
        }}

    def test_screen_has_six_then_nine_unique_settings(self):
        first = screen_settings(self.config)
        final = screen_settings(self.config, 1.0)
        self.assertEqual(len(first), 6)
        self.assertEqual(len(final), 9)
        self.assertEqual(len({(x["ewc_lambda"], x["replay_fraction"]) for x in final}), 9)
        duplicate = [x for x in final if (x["ewc_lambda"], x["replay_fraction"]) == (1.0, 0.2)]
        self.assertEqual(len(duplicate), 1)

    def test_best_lambda_cannot_leave_grid_or_be_zero(self):
        for value in (0, 3):
            with self.assertRaises(RunPlanError):
                screen_settings(self.config, value)

    def test_condition_or_stage_cannot_be_added_as_covariate(self):
        source = Path(__file__).parents[1] / "config_v1.json"
        config = json.loads(source.read_text())
        config["features"]["categorical_covariate_keys"] = ["stage"]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.json"
            path.write_text(json.dumps(config))
            with self.assertRaisesRegex(ConfigError, "covariates are forbidden"):
                load_config(path)

    def test_secondary_stress_adapts_powered_queries_plus_target(self):
        config = {
            "_config_sha256": "cfg",
            "screen": {
                **self.config["screen"], "confirmation_seeds": [17, 41, 73, 101, 137],
            },
            "lineages": ["H", "M", "F", "C", "T"],
            "evaluation": {
                "powered_query_studies": ["P1", "P2"],
                "secondary_stress_studies": ["S1", "S2"],
                "mapping_only_studies": ["MAP"],
            },
            "acquisition_orders": {
                "accession": ["S1", "S2", "MAP", "P1", "P2"],
                "reverse_accession": ["P2", "P1", "MAP", "S2", "S1"],
            },
        }
        selection = {
            "lock_sha256": "selection",
            "selected": {"ewc_lambda": 1.0, "replay_fraction": 0.2},
            "pareto_runner_up": None,
        }
        with tempfile.TemporaryDirectory() as directory, patch(
            "masld_cl.run_plan.load_selection_lock", return_value=selection
        ):
            plan = create_run_plan(
                config, Path(directory) / "plan.json", selection_lock="selection.json"
            )
        for stress in plan["secondary_stress_tests"]:
            self.assertEqual(
                stress["adaptation_datasets"], ["P1", "P2", stress["stress_dataset"]]
            )
            self.assertEqual(stress["control_fisher_datasets"], ["P1", "P2"])
            self.assertEqual(stress["held_out_datasets"], [])
