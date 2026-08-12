#!/usr/bin/env python3

from __future__ import annotations

import sys
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from spatial_resource_lib import (  # noqa: E402
    EXPECTED_N_PROGRAMS,
    SpatialResourceError,
    assert_authorized_inference,
    coverage_state,
    load_dataset_registry,
    load_frozen_programs,
    validate_hmsma_metadata_join,
)


class SpatialResourceContractTests(unittest.TestCase):
    def setUp(self):
        self.project = SCRIPT_ROOT.parents[3]
        self.registry = load_dataset_registry(SCRIPT_ROOT / "spatial_resource_datasets.tsv")
        self.by_id = {row["dataset_id"]: row for row in self.registry}

    def test_real_frozen_registry_has_117_programs_and_two_confirmatory(self):
        hotspot = (
            self.project / "Analysis/Multimodal_Program_Projection/candidates"
            / "program-context-v2-candidate-2026-08-07/hotspot"
        )
        programs, weights = load_frozen_programs(hotspot)
        self.assertEqual(len(programs), EXPECTED_N_PROGRAMS)
        self.assertEqual(sum(row["external_test_eligible"] == "TRUE" for row in programs.values()), 2)
        self.assertEqual(set(programs), set(weights))

    def test_coverage_threshold_is_outcome_free_and_exact(self):
        self.assertEqual(coverage_state(8, 0.20)[0], "observable")
        self.assertEqual(coverage_state(7, 0.99)[0], "partial")
        self.assertEqual(coverage_state(100, 0.199999)[0], "partial")
        self.assertEqual(coverage_state(0, 0.0)[0], "untestable")

    def test_unresolved_units_fail_closed(self):
        with self.assertRaises(SpatialResourceError):
            assert_authorized_inference(self.by_id["Vu_et_al_2025"], "donor")
        with self.assertRaises(SpatialResourceError):
            assert_authorized_inference(self.by_id["HRA007511_HMSMA"], "clinical")
        with self.assertRaises(SpatialResourceError):
            assert_authorized_inference(self.by_id["HRA007511_HMSMA"], "histology")
        with self.assertRaises(SpatialResourceError):
            assert_authorized_inference(self.by_id["HRA007511_HMSMA"], "metabolite")
        assert_authorized_inference(self.by_id["GSE192741"], "donor")

    def test_exclusions_remain_explicit(self):
        self.assertEqual(self.by_id["Govaere2026_GeoMx"]["dataset_gate"], "dropped")
        self.assertEqual(self.by_id["GSE287826"]["dataset_gate"], "metadata_pending")
        self.assertEqual(self.by_id["OnTraC"]["dataset_gate"], "skipped")

    @staticmethod
    def authoritative_join_rows():
        base = {
            "array_relationship": "primary",
            "run_aggregation_id": "run_group_1",
            "clinical_label": "MASH",
            "nas_score": "5",
            "fibrosis_score": "2",
            "age": "52",
            "sex": "female",
            "bmi": "31.2",
            "acquisition_batch": "batch_1",
            "histology_registration_status": "registered",
            "histology_image_id": "he_1",
            "histology_orientation": "native_top_left",
            "scrna_overlap_status": "no_overlap",
            "maldi_sample_id": "maldi_1",
            "maldi_registration_status": "registered",
        }
        return [
            {**base, "array_id": "HRA_01", "donor_id": "donor_1"},
            {**base, "array_id": "HRA_02", "donor_id": "donor_2"},
        ]

    def test_authoritative_hmsma_join_accepts_exact_registered_mapping(self):
        rows = self.authoritative_join_rows()
        validate_hmsma_metadata_join(rows, ["HRA_01", "HRA_02"], "clinical")
        validate_hmsma_metadata_join(rows, ["HRA_01", "HRA_02"], "histology")
        validate_hmsma_metadata_join(rows, ["HRA_01", "HRA_02"], "metabolite")

    def test_hmsma_join_rejects_duplicate_or_unmatched_array_ids(self):
        duplicated = self.authoritative_join_rows()
        duplicated[1]["array_id"] = "HRA_01"
        with self.assertRaises(SpatialResourceError):
            validate_hmsma_metadata_join(duplicated, ["HRA_01", "HRA_02"])
        with self.assertRaises(SpatialResourceError):
            validate_hmsma_metadata_join(self.authoritative_join_rows(), ["HRA_01", "HRA_03"])

    def test_hmsma_join_rejects_ambiguous_repeated_donor_relationship(self):
        rows = self.authoritative_join_rows()
        rows[1]["donor_id"] = "donor_1"
        with self.assertRaises(SpatialResourceError):
            validate_hmsma_metadata_join(rows, ["HRA_01", "HRA_02"])
        rows[1]["array_relationship"] = "repeated_section"
        validate_hmsma_metadata_join(rows, ["HRA_01", "HRA_02"])

    def test_hmsma_join_rejects_unregistered_histology_and_maldi(self):
        rows = self.authoritative_join_rows()
        rows[0]["histology_registration_status"] = "unregistered"
        with self.assertRaises(SpatialResourceError):
            validate_hmsma_metadata_join(rows, ["HRA_01", "HRA_02"], "histology")
        rows = self.authoritative_join_rows()
        rows[0]["maldi_registration_status"] = "unregistered"
        with self.assertRaises(SpatialResourceError):
            validate_hmsma_metadata_join(rows, ["HRA_01", "HRA_02"], "metabolite")

    def test_hmsma_join_treats_missing_tsv_cells_as_missing(self):
        rows = self.authoritative_join_rows()
        rows[0]["run_aggregation_id"] = None
        with self.assertRaises(SpatialResourceError):
            validate_hmsma_metadata_join(rows, ["HRA_01", "HRA_02"])


if __name__ == "__main__":
    unittest.main()
