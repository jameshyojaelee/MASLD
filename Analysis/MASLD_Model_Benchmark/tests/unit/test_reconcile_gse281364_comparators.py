from __future__ import annotations

import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest

from scripts import reconcile_gse281364_comparators as reconcile


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/gse281364_comparator_reconciliation.toml"


class ComparatorReconciliationTests(unittest.TestCase):
    def test_config_freezes_outcome_and_claim_firewalls(self) -> None:
        config = reconcile.load_config(CONFIG)
        reconcile.validate_config(config)
        self.assertFalse(config["outcome_access_authorized"])
        self.assertFalse(config["metric_calculation_authorized"])
        self.assertFalse(config["shortlist_authorized"])
        self.assertFalse(config["complementarity_authorized"])
        self.assertFalse(
            config["completion_firewall"][
                "deterministic_schema_seed_repeats_are_independent_fits"
            ]
        )
        self.assertTrue(
            config["completion_firewall"][
                "stochastic_stability_claim_blocked_for_deterministic_heads"
            ]
        )

    def test_firewall_drift_is_rejected(self) -> None:
        config = reconcile.load_config(CONFIG)
        changed = copy.deepcopy(config)
        changed["metric_calculation_authorized"] = True
        with self.assertRaises(reconcile.ReconciliationError):
            reconcile.validate_config(changed)
        changed = copy.deepcopy(config)
        changed["completion_firewall"][
            "deterministic_schema_seed_repeats_are_independent_fits"
        ] = True
        with self.assertRaises(reconcile.ReconciliationError):
            reconcile.validate_config(changed)

    def test_candidate_summary_reader_consumes_header_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "summary.tsv"
            path.write_text(
                "model_id\thead_id\tpositive_gain_seeds\tdevelopment_shortlist\tshortlist_blocked\tfinalist_claim\n"
                "do_not_parse\tdo_not_parse\t999\ttrue\tfalse\ttrue\n",
                encoding="utf-8",
            )
            fields = reconcile.read_candidate_summary_header(path)
            self.assertIn("positive_gain_seeds", fields)

    def test_full_reconciliation_has_exact_coverage_and_no_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "reconciliation"
            arguments = type(
                "Arguments",
                (),
                {"root": ROOT, "config": CONFIG, "output": output},
            )()
            receipt = reconcile.reconcile(arguments)
            self.assertEqual(receipt["status"], "pass_outcome_blind_join_coverage_disposition")
            self.assertEqual(receipt["mpralegnet_selected_elements"], 1033)
            self.assertEqual(receipt["row_keys"], 10330)
            self.assertEqual(receipt["mpralegnet_reassigned_from_original_fold"], 830)
            self.assertFalse(receipt["candidate_summary_metric_values_read"])
            self.assertFalse(receipt["metrics_calculated"])
            self.assertFalse(receipt["outcomes_read"])
            self.assertTrue(receipt["shortlist_blocked"])
            with (output / "comparator_disposition.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                dispositions = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(dispositions), 5)
            self.assertEqual(dispositions[0]["joined_row_keys"], "10330")
            self.assertTrue(
                all(row["valid_for_current_shortlist"] == "false" for row in dispositions)
            )
            audit = json.loads(
                (output / "terminal_interpretation_audit.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertFalse(
                audit[
                    "positive_gain_seeds_for_deterministic_heads_is_stability_evidence"
                ]
            )
            self.assertFalse(
                audit["four_of_five_stability_gate_satisfied_by_deterministic_repeats"]
            )


if __name__ == "__main__":
    unittest.main()
