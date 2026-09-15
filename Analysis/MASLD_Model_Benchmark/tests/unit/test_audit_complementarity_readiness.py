from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import tomllib
import unittest

from scripts.audit_complementarity_readiness import (
    ComplementarityReadinessError,
    LANES,
    PAIRS,
    SOURCE_FIELDS,
    STATUS,
    preflight,
    validate_authorities,
    validate_config,
    validate_prediction_metadata,
    validate_trigger_authorities,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "config/complementarity_readiness_audit.toml"


class ComplementarityReadinessAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = tomllib.loads(CONFIG_PATH.read_text(encoding="utf-8"))

    def test_trigger_requires_six_unbound_development_sources(self) -> None:
        validate_config(self.config)
        self.assertEqual(
            tuple(row["lane_id"] for row in self.config["prediction_lanes"]), LANES
        )
        self.assertEqual(
            tuple(row["pair_id"] for row in self.config["pair_audits"]), PAIRS
        )
        self.assertEqual(
            set(self.config["required_stacking_sources"]), set(SOURCE_FIELDS)
        )
        self.assertTrue(
            all(
                value == "UNBOUND"
                for value in self.config["required_stacking_sources"].values()
            )
        )
        self.assertFalse(self.config["build_authorized"])
        self.assertFalse(self.config["stack_fit_or_residual_correlation_authorized"])

    def test_no_lane_or_pair_is_prematurely_ready(self) -> None:
        self.assertTrue(
            all(not row["residual_ready"] for row in self.config["prediction_lanes"])
        )
        self.assertTrue(
            all(not row["stack_gain_ready"] for row in self.config["prediction_lanes"])
        )
        self.assertTrue(
            all(not row["trigger_ready"] for row in self.config["pair_audits"])
        )
        gse244832 = next(
            row
            for row in self.config["prediction_lanes"]
            if row["lane_id"] == "gse244832_rna_atac_development"
        )
        self.assertEqual(gse244832["authority_ids"], [])
        self.assertEqual(gse244832["prediction_schema"], "UNBOUND")

    def test_authorities_and_metadata_match_frozen_receipts(self) -> None:
        authorities = validate_authorities(ROOT, self.config)
        self.assertEqual(len(authorities), 19)
        validate_trigger_authorities(ROOT, self.config)
        metadata = validate_prediction_metadata(ROOT, self.config)
        self.assertEqual(metadata["local_profile_evaluation_groups"], 7)
        self.assertEqual(metadata["corgi_lineages"], ["hepatocyte"])
        self.assertEqual(metadata["midas_seed_ids"], [2711])
        self.assertEqual(metadata["scbasset_seed_ids"], [11, 20260824])
        self.assertEqual(metadata["dna_lm_family_promotion"], "none")
        self.assertFalse(metadata["borzoi_prediction_available"])
        self.assertFalse(metadata["gse244832_prediction_authority_bound"])

    def test_claim_or_source_injection_is_rejected(self) -> None:
        changed = deepcopy(self.config)
        changed["build_authorized"] = True
        with self.assertRaises(ComplementarityReadinessError):
            validate_config(changed)

        changed = deepcopy(self.config)
        changed["required_stacking_sources"]["variant_development_shortlist"] = (
            "unverified/path"
        )
        with self.assertRaises(ComplementarityReadinessError):
            validate_config(changed)

        changed = deepcopy(self.config)
        changed["pair_audits"][2]["trigger_ready"] = True
        with self.assertRaises(ComplementarityReadinessError):
            validate_config(changed)

    def test_integrated_audit_reads_only_prediction_row_identifiers(self) -> None:
        receipt = preflight(root=ROOT, config_path=CONFIG_PATH)
        self.assertEqual(receipt["status"], STATUS)
        self.assertEqual(receipt["authority_rows"], 19)
        self.assertEqual(receipt["prediction_lane_rows"], 8)
        self.assertEqual(receipt["pair_audit_rows"], 6)
        self.assertEqual(receipt["required_stacking_source_rows"], 6)
        self.assertEqual(receipt["bound_stacking_source_rows"], 0)
        self.assertEqual(receipt["trigger_ready_pair_rows"], 0)
        self.assertEqual(len(receipt["scbasset_midas_row_overlap"]), 5)
        self.assertTrue(
            all(
                row["same_table_schema_sha256"]
                for row in receipt["scbasset_midas_row_overlap"]
            )
        )
        self.assertFalse(receipt["residual_correlation_calculated"])
        self.assertFalse(receipt["stack_weights_fit"])
        self.assertFalse(receipt["prediction_values_read"])
        self.assertTrue(receipt["row_identifiers_read"])
        self.assertFalse(receipt["raw_or_summary_outcomes_read"])
        self.assertFalse(receipt["sealed_assets_read"])
        self.assertFalse(receipt["novel_model_built_or_fit"])
        self.assertFalse(receipt["conditional_model_build_authorized"])


if __name__ == "__main__":
    unittest.main()
