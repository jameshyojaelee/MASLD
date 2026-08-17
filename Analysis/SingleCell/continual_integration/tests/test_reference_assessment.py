from __future__ import annotations

import copy
import csv
import tempfile
import unittest
from pathlib import Path

from masld_cl.reference_assessment import (
    CURRENT_REFERENCE_DONORS,
    ReferenceAssessmentError,
    load_policy,
    validate_local_arms,
    validate_registry,
)


class TestReferenceAssessment(unittest.TestCase):
    def setUp(self):
        self.root = Path(__file__).resolve().parents[1]
        self.policy_path = self.root / "reference_policy_v1.json"
        self.registry_path = self.root / "external_reference_registry_v1.tsv"

    def test_production_policy_and_registry(self):
        policy = load_policy(self.policy_path)
        registry = validate_registry(self.registry_path)
        self.assertTrue(policy["primary_reference_unchanged"])
        self.assertEqual(policy["external_arm"]["expected_reference_donors"], 26)
        self.assertEqual(
            policy["external_arm"]["current_reference_alias_overlap"],
            {
                "GSE185477_D02": "C58",
                "GSE185477_D03": "C70",
                "GSE185477_D04": "C72",
            },
        )
        self.assertFalse(
            policy["external_arm"]["compatibility_audit"]
            ["eligible_for_identical_current_4000_hvgs"]
        )
        self.assertEqual(
            policy["external_arm"]["compatibility_audit"]["mapped_current_hvgs"],
            3802,
        )
        self.assertEqual(
            policy["external_arm"]["compatibility_audit"]["missing_current_hvgs"],
            198,
        )
        self.assertEqual(
            {row["candidate_id"] for row in registry if row["reference_use"] == "include_external_clean"},
            {"andrews_2022", "andrews_2024"},
        )

    def test_whole_hlica_cannot_be_enabled(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "registry.tsv"
            rows = list(csv.DictReader(self.registry_path.open(), delimiter="\t"))
            rows[0]["reference_use"] = "include_external_clean"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ReferenceAssessmentError, "whole-HLiCA"):
                validate_registry(path)

    def test_query_control_leakage_fails(self):
        policy = load_policy(self.policy_path)
        policy = copy.deepcopy(policy)
        policy["local_arms"]["current_strict_7"]["added_donors"] = []
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "donors.tsv"
            fields = [
                "donor_id", "dataset", "harmonized_stage", "strict_reference",
                "query_control", "analysis_eligible", "n_cells_analyzed",
            ]
            rows = []
            for donor in sorted(CURRENT_REFERENCE_DONORS):
                rows.append({
                    "donor_id": donor,
                    "dataset": donor.split("_")[0],
                    "harmonized_stage": "Healthy",
                    "strict_reference": "True",
                    "query_control": "True" if donor == "GSE185477_D02" else "False",
                    "analysis_eligible": "True",
                    "n_cells_analyzed": "1",
                })
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ReferenceAssessmentError, "leaks query controls"):
                validate_local_arms(policy, path)

    def test_external_donor_arithmetic_fails_closed(self):
        policy = load_policy(self.policy_path)
        policy = copy.deepcopy(policy)
        policy["external_arm"]["new_unique_donors"] = 20
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "policy.json"
            import json
            path.write_text(json.dumps(policy))
            with self.assertRaisesRegex(ReferenceAssessmentError, "donor arithmetic"):
                load_policy(path)

    def test_external_feature_compatibility_fails_closed(self):
        policy = copy.deepcopy(load_policy(self.policy_path))
        policy["external_arm"]["compatibility_audit"][
            "eligible_for_identical_current_4000_hvgs"
        ] = True
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "policy.json"
            import json
            path.write_text(json.dumps(policy))
            with self.assertRaisesRegex(ReferenceAssessmentError, "feature-compatibility"):
                load_policy(path)
