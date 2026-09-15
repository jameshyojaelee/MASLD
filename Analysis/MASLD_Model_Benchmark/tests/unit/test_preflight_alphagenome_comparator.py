from __future__ import annotations

import json
from pathlib import Path
import unittest

from scripts.preflight_alphagenome_comparator import (
    RELEVANT_CURIES,
    audit_metadata,
    parse_metadata_block,
    preflight,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "config/alphagenome_comparator_preflight.json"


class AlphaGenomeComparatorPreflightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

    def test_current_release_is_restricted_and_fail_closed(self) -> None:
        self.assertEqual(
            self.config["status"],
            "terminal_current_release_restricted_and_execution_blocked",
        )
        terms = self.config["terms_contract"]
        self.assertFalse(terms["open_champion_eligible"])
        self.assertFalse(terms["conditional_open_model_training_allowed"])
        self.assertFalse(terms["weight_or_derived_weight_redistribution_allowed"])
        self.assertFalse(
            self.config["local_checkpoint_route"]["current_execution_allowed"]
        )
        self.assertFalse(self.config["api_route"]["current_execution_allowed"])
        self.assertFalse(self.config["api_route"]["api_key_inspection_allowed"])

    def test_metadata_parser_preserves_assay_and_biosample_fields(self) -> None:
        block = """    metadata {
      name: "CL:0000182 total RNA-seq"
      strand: STRAND_POSITIVE
      ontology_term {
        ontology_type: ONTOLOGY_TYPE_CL
        id: 182
      }
      biosample {
        type: BIOSAMPLE_TYPE_PRIMARY_CELL
        name: "hepatocyte"
        stage: "adult"
      }
      assay: "RNA-seq"
      data_source: "encode"
    }"""
        parsed = parse_metadata_block("OUTPUT_TYPE_RNA_SEQ", block)
        self.assertEqual(parsed["output_type"], "RNA_SEQ")
        self.assertEqual(parsed["ontology_curie"], "CL:0000182")
        self.assertEqual(parsed["track_name"], "CL:0000182 total RNA-seq")
        self.assertEqual(parsed["strand"], "STRAND_POSITIVE")
        self.assertEqual(parsed["biosample_name"], "hepatocyte")
        self.assertEqual(parsed["biosample_type"], "PRIMARY_CELL")
        self.assertEqual(parsed["biosample_stage"], "adult")
        self.assertEqual(parsed["assay"], "RNA-seq")
        self.assertEqual(parsed["data_source"], "encode")

    def test_frozen_human_metadata_roster_is_exact(self) -> None:
        metadata = (
            ROOT
            / self.config["frozen_authorities"]["admission_sources"]["path"]
            / "sources/alphagenome_human_metadata"
        )
        audit = audit_metadata(metadata, self.config)
        self.assertEqual(audit["encoded_output_slots"], 6126)
        self.assertEqual(audit["padding_slots"], 563)
        self.assertEqual(audit["nonpadding_metadata_rows"], 5563)
        self.assertEqual(audit["splice_junction_metadata_rows"], 367)
        self.assertEqual(
            audit["active_tracks_with_second_splice_junction_strands"], 5930
        )
        self.assertEqual(
            audit["relevant_ontology_counts"],
            self.config["metadata_contract"]["relevant_ontology_expected_counts"],
        )
        self.assertEqual(
            {
                curie: audit["relevant_ontology_counts_by_output"][curie]
                for curie in self.config["metadata_contract"][
                    "cell_ontology_expected_counts_by_output"
                ]
            },
            self.config["metadata_contract"][
                "cell_ontology_expected_counts_by_output"
            ],
        )
        self.assertEqual(audit["cholangiocyte_exact_name_tracks"], 0)
        self.assertEqual(audit["kupffer_exact_name_tracks"], 0)
        self.assertEqual(
            len(audit["relevant_rows"]),
            sum(
                self.config["metadata_contract"][
                    "relevant_ontology_expected_counts"
                ].values()
            ),
        )

    def test_api_route_does_not_claim_an_immutable_server_build(self) -> None:
        api = self.config["api_route"]
        self.assertEqual(api["sdk_tag"], "v0.8.0")
        self.assertEqual(
            api["sdk_revision"], "71a6beb8c30832f121309a81c2530efa5af7986a"
        )
        self.assertEqual(api["logical_model_version"], "ALL_FOLDS")
        self.assertFalse(api["default_model_version_allowed"])
        self.assertIn(
            "response_does_not_echo_model_version_or_server_checkpoint_digest",
            api["reproducibility_blockers"],
        )
        self.assertIn(
            "metadata_request_has_no_model_version_field",
            api["reproducibility_blockers"],
        )
        self.assertFalse(api["api_prediction_called"])
        self.assertFalse(api["api_metadata_called"])

    def test_integrated_receipt_is_terminal_without_model_or_api_access(self) -> None:
        receipt = preflight(root=ROOT, config_path=CONFIG_PATH)
        self.assertFalse(receipt["executable"])
        self.assertTrue(receipt["terminal_current_release"])
        self.assertEqual(len(receipt["relevant_track_rows"]), 93)
        self.assertEqual(
            set(receipt["metadata_audit"]["relevant_ontology_counts"]),
            set(RELEVANT_CURIES),
        )
        self.assertEqual(len(receipt["blockers"]), 21)
        self.assertFalse(receipt["local_checkpoint_route"]["checkpoint_downloaded"])
        self.assertFalse(receipt["local_checkpoint_route"]["model_forward_executed"])
        self.assertFalse(receipt["api_route"]["api_connection_attempted"])
        self.assertFalse(receipt["api_route"]["metadata_endpoint_called"])
        self.assertFalse(receipt["api_route"]["prediction_endpoint_called"])
        self.assertFalse(receipt["outcomes_read"])
        self.assertFalse(receipt["sealed_assets_read"])
        self.assertFalse(receipt["metrics_calculated"])
        self.assertFalse(receipt["open_champion_eligible"])


if __name__ == "__main__":
    unittest.main()
