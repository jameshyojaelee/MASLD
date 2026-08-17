from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from masld_cl.config import canonical_json_bytes, load_config
from masld_cl.contracts import sha256_path
from masld_cl.firewall import FirewallError, assert_control_only_metric_names, load_selection_lock
from masld_cl.selection import SelectionError, _validate_exact_grid, select_setting


ROOT = Path(__file__).resolve().parents[1]


def row(setting, ewc, replay, metric, value):
    return {
        "setting_id": setting, "ewc_lambda": str(ewc), "replay_fraction": str(replay),
        "metric": metric, "value": str(value),
    }


def setting_rows(setting, ewc, replay, shift, se=0.02):
    return [
        row(setting, ewc, replay, "reference_macro_f1_change_ci_low", -0.01),
        row(setting, ewc, replay, "reference_neighborhood_jaccard_loss", 0.03),
        row(setting, ewc, replay, "shift_control", shift),
        row(setting, ewc, replay, "shift_control_standard_error", se),
    ]


class TestSelection(unittest.TestCase):
    def setUp(self):
        self.config = load_config(ROOT / "config_v1.json")

    def test_outcomes_are_forbidden_before_selection(self):
        for name in ("case_distance", "stage_beta", "program_q", "cas13_hit"):
            with self.assertRaises(FirewallError):
                assert_control_only_metric_names([name])

    def test_one_se_prefers_smaller_lambda_then_replay(self):
        rows = []
        rows += setting_rows("large", 100, 0.4, 0.50, 0.03)
        rows += setting_rows("small", 1, 0.2, 0.52, 0.03)
        rows += setting_rows("zero_ablation", 0, 0.2, 0.1, 0.01)
        selected = select_setting(self.config, rows)["selected"]
        self.assertEqual(selected["setting_id"], "small")

    def test_reference_failure_blocks_setting(self):
        rows = setting_rows("bad", 1, 0.2, 0.4)
        rows[0]["value"] = "-0.5"
        with self.assertRaises(SelectionError):
            select_setting(self.config, rows)

    def test_selection_lock_detects_tampering(self):
        import hashlib
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sources = {}
            for name in ("metrics.tsv", "details.json", "provisional.json"):
                source = root / name
                source.write_text(name)
                sources[name] = source
            records = []
            for index in range(60):
                record = root / f"record_{index}.json"
                record.write_text(str(index))
                records.append({"path": str(record), "sha256": sha256_path(record)})
            lock = {
                "schema_version": "masld-cl-selection-v1",
                "config_sha256": self.config["_config_sha256"],
                "selection_frozen": True,
                "selected": {"ewc_lambda": 1.0, "replay_fraction": 0.2},
                "control_metrics_realpath": str(sources["metrics.tsv"]),
                "control_metrics_sha256": sha256_path(sources["metrics.tsv"]),
                "aggregate_details_realpath": str(sources["details.json"]),
                "aggregate_details_sha256": sha256_path(sources["details.json"]),
                "provisional_lambda_lock_realpath": str(sources["provisional.json"]),
                "provisional_lambda_lock_sha256": sha256_path(sources["provisional.json"]),
                "execution_records": records,
            }
            lock["lock_sha256"] = hashlib.sha256(canonical_json_bytes(lock)).hexdigest()
            path = Path(directory) / "lock.json"
            path.write_text(json.dumps(lock))
            with patch("masld_cl.execution.verify_execution_record", return_value={}):
                load_selection_lock(path, self.config)
            lock["selected"]["ewc_lambda"] = 100
            path.write_text(json.dumps(lock))
            with self.assertRaises(FirewallError):
                load_selection_lock(path, self.config)

    def test_duplicate_model_metrics_must_be_aggregated_first(self):
        rows = setting_rows("a", 1.0, 0.2, shift=0.9, se=0.1)
        rows.append(dict(rows[-1]))
        with self.assertRaisesRegex(SelectionError, "aggregate the six model kinds"):
            select_setting(self.config, rows)

    def test_selection_requires_exact_six_then_nine_setting_grids(self):
        first = []
        for value in self.config["screen"]["lambda_values"]:
            first.extend(setting_rows(f"lambda_{value:g}__replay_0.2", value, 0.2, 1.0))
        _validate_exact_grid(self.config, first, best_lambda=None)
        with self.assertRaisesRegex(SelectionError, "exact prespecified"):
            _validate_exact_grid(self.config, first[:-4], best_lambda=None)
        final = list(first)
        for replay in (0.0, 0.1, 0.4):
            final.extend(setting_rows(f"lambda_1__replay_{replay:g}", 1, replay, 1.0))
        _validate_exact_grid(self.config, final, best_lambda=1.0)
        with self.assertRaisesRegex(SelectionError, "exact prespecified"):
            _validate_exact_grid(self.config, final[:-4], best_lambda=1.0)


if __name__ == "__main__":
    unittest.main()
