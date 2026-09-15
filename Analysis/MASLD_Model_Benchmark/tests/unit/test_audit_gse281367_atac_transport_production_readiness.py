from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
import tomllib
import unittest

from scripts.audit_gse281367_atac_transport_production_readiness import (
    ProductionReadinessError,
    audit_axis_header,
    audit_windows,
    run_audit,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = (
    ROOT / "config/evaluation/gse281367_atac_transport_production_readiness.toml"
)


class GSE281367ATACTransportProductionReadinessTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = tomllib.loads(CONTRACT_PATH.read_text(encoding="utf-8"))

    def test_contract_retains_fail_closed_production_boundary(self) -> None:
        validate_contract(self.contract)
        decision = self.contract["decision"]
        self.assertTrue(decision["frozen_outcome_authority_valid"])
        self.assertTrue(decision["donor_grouping_contract_valid"])
        self.assertFalse(decision["label_agnostic_production_ready"])
        self.assertFalse(decision["observed_atac_production_ready"])
        self.assertFalse(decision["existing_outcome_artifact_may_be_given_to_model_jobs"])

    def test_firewall_cannot_be_relaxed(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["production_firewall"]["condition_blind_scorer_registered"] = True
        with self.assertRaises(ProductionReadinessError):
            validate_contract(changed)
        changed = copy.deepcopy(self.contract)
        changed["condition_values_read"] = True
        with self.assertRaises(ProductionReadinessError):
            validate_contract(changed)

    def test_axis_audit_reads_only_header_and_never_emits_label_values(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "axis.tsv"
            fields = list(self.contract["expected_axis"]["fields"])
            path.write_text(
                "\t".join(fields)
                + "\n"
                + "gse281367\t0\tD01\t0\thepatocyte\tSECRET_CASE_LABEL\t0\t51\ttrue\n",
                encoding="utf-8",
            )
            result = audit_axis_header(path, fields)
            self.assertEqual(result["data_rows_read"], 0)
            self.assertFalse(result["condition_values_read"])
            self.assertFalse(result["production_exchange_label_safe"])
            self.assertNotIn("SECRET_CASE_LABEL", json.dumps(result))

    def test_window_contract_rejects_cross_role_contig_leakage(self) -> None:
        expected = {
            "window_n": 2,
            "window_width_bp": 1000,
            "profile_bins": 20,
            "profile_bin_width_bp": 50,
            "valid_window_n": 1,
            "test_window_n": 1,
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            roles = root / "fold.json"
            roles.write_text(
                json.dumps({"train": ["chr3"], "valid": ["chr2"], "test": ["chr1"]}),
                encoding="utf-8",
            )
            windows = root / "windows.tsv"
            windows.write_text(
                "window_index\twindow_id\trole\tcontig\tstart\tend\n"
                "0\tw0\tvalid\tchr2\t0\t1000\n"
                "1\tw1\ttest\tchr1\t100\t1100\n",
                encoding="utf-8",
            )
            result = audit_windows(windows, roles, expected)
            self.assertEqual(result["role_counts"], {"test": 1, "valid": 1})
            self.assertFalse(result["train_windows_present"])
            windows.write_text(
                "window_index\twindow_id\trole\tcontig\tstart\tend\n"
                "0\tw0\tvalid\tchr2\t0\t1000\n"
                "1\tw1\ttest\tchr2\t100\t1100\n",
                encoding="utf-8",
            )
            with self.assertRaises(ProductionReadinessError):
                audit_windows(windows, roles, expected)

    def test_real_authorities_are_hash_bound_without_outcome_or_label_read(self) -> None:
        audit = run_audit(ROOT, self.contract)
        self.assertEqual(audit["artifact_integrity"]["outcome_array_files_opened"], 0)
        self.assertEqual(audit["artifact_integrity"]["outcome_array_files_rehashed"], 0)
        self.assertFalse(audit["read_policy"]["condition_values_read"])
        self.assertFalse(audit["read_policy"]["sealed_outcomes_read"])
        self.assertEqual(audit["genomic_window_contract"]["window_n"], 32000)
        self.assertEqual(audit["donor_safety"]["dataset_donor_n"], 12)
        self.assertFalse(audit["decision"]["family_native_tournament_ready"])


if __name__ == "__main__":
    unittest.main()
