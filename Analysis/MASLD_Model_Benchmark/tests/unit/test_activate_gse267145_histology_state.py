from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from scripts.activate_gse267145_histology_state import (
    HistologyActivationError,
    activation_contract,
    build_endpoint_records,
    stage3,
    validate_task_contract,
)


ROOT = Path(__file__).resolve().parents[2]
JOIN = ROOT / "executions/gse267145-authoritative-join-21064930/participant_join.tsv"
TASK = ROOT / "config/evaluation/gse267145_histology_state_task.toml"
GATE = ROOT / "config/evaluation/gse267145_histology_state_promotion_gate.json"


def join_rows() -> list[dict[str, str]]:
    with JOIN.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class GSE267145HistologyStateTests(unittest.TestCase):
    def test_source_stage_collapse_is_explicit(self) -> None:
        self.assertEqual(stage3("NOR"), "NOR")
        self.assertEqual(stage3("NAFL"), "NAFL")
        for value in ("NASH_F0", "NASH_F1", "NASH_F23"):
            self.assertEqual(stage3(value), "NASH")
        with self.assertRaises(HistologyActivationError):
            stage3("MASH")

    def test_frozen_join_reproduces_histology_census(self) -> None:
        records, summary = build_endpoint_records(join_rows())
        self.assertEqual(len(records), 99)
        self.assertEqual(summary["source_stage5_counts"]["NOR"], 24)
        self.assertEqual(summary["stage3_counts"], {"NAFL": 28, "NASH": 47, "NOR": 24})
        self.assertEqual(summary["nash_crn_component_sum_counts"]["0"], 24)
        self.assertEqual(summary["fibrosis_counts"], {"0": 71, "1": 15, "2": 9, "3": 4})
        self.assertEqual(summary["recorded_sex_counts"], {"F": 85, "M": 14})
        self.assertFalse(summary["lobular_necrosis_in_nas"])
        self.assertTrue(all(record["histology_state"] == "observed" for record in records))

    def test_source_stage_fibrosis_inconsistency_fails_closed(self) -> None:
        rows = join_rows()
        rows[0] = dict(rows[0])
        rows[0]["fibrosis"] = "3"
        with self.assertRaises(HistologyActivationError):
            build_endpoint_records(rows)

    def test_task_is_development_only_and_blocks_published_feature_sets(self) -> None:
        observed = validate_task_contract(TASK, GATE)
        self.assertEqual(observed["task_id"], "paired_bulk_histology_state")
        self.assertEqual(observed["primary_metric"], "participant_macro_f1_stage3")
        self.assertFalse(observed["champion_eligible"])
        text = TASK.read_text(encoding="utf-8").lower()
        self.assertIn("differential regions", text)
        self.assertIn("differential genes", text)
        self.assertIn("fit feature selection inside training folds", text)

    def test_promotable_gate_is_rejected(self) -> None:
        value = json.loads(GATE.read_text(encoding="utf-8"))
        value["champion_eligible"] = True
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "gate.json"
            path.write_text(json.dumps(value), encoding="utf-8")
            with self.assertRaises(HistologyActivationError):
                validate_task_contract(TASK, path)

    def test_activation_keeps_query_modalities_and_claims_separate(self) -> None:
        observed = activation_contract({"task_id": "paired_bulk_histology_state"})
        self.assertTrue(
            observed["allowed_lanes"]["observed_pair_multimodal"]
            ["cannot_support_rna_only_claim"]
        )
        self.assertEqual(
            observed["metadata_contract"]["participant_age"], "structurally_missing"
        )
        self.assertFalse(observed["feature_contract"]["published_outcome_selected_features_allowed"])
        self.assertFalse(observed["claim_contract"]["champion_eligible"])


if __name__ == "__main__":
    unittest.main()
