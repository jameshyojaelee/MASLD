from __future__ import annotations

import copy
from pathlib import Path
import tomllib
import unittest

from scripts.audit_gse281574_pooling_topology import (
    PoolingTopologyError,
    build_library_pairs,
    task_decision,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/evaluation/gse281574_pooling_topology.toml"


class GSE281574PoolingTopologyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = tomllib.loads(CONTRACT.read_text(encoding="utf-8"))

    def test_people_libraries_and_geo_records_are_not_conflated(self) -> None:
        pairs, _ = validate_contract(self.contract)
        self.assertEqual(self.contract["reported_person_n"], 15)
        self.assertEqual(self.contract["physical_multiome_library_n"], 6)
        self.assertEqual(self.contract["geo_assay_record_n"], 12)
        self.assertEqual(len(pairs), 6)
        self.assertFalse(self.contract["pooled_libraries_are_donors"])
        self.assertFalse(self.contract["cells_are_biological_replicates"])

    def test_numbered_pool_sizes_and_donor_joins_remain_unresolved(self) -> None:
        validate_contract(self.contract)
        self.assertEqual(sorted(self.contract["pool_size_multiset_per_condition"]), [2, 3])
        self.assertEqual(
            self.contract["numbered_pool_size_assignment"], "UNRESOLVED_PUBLIC_METADATA"
        )
        self.assertEqual(
            self.contract["nucleus_to_individual_join"], "UNAVAILABLE_PUBLIC_METADATA"
        )
        self.assertTrue(all(row["person_n"] == -1 for row in self.contract["assay_record_pair"]))

    def test_only_pool_grouped_descriptive_tasks_are_allowed(self) -> None:
        _, tasks = validate_contract(self.contract)
        decision = task_decision(tasks)
        self.assertFalse(decision["donor_safe_split_available"])
        self.assertTrue(decision["library_grouped_technical_split_available"])
        self.assertEqual(decision["condition_level_independent_unit_n"], 2)
        self.assertFalse(decision["model_selection_allowed"])
        self.assertFalse(decision["external_masld_champion_eligible"])

    def test_normal_pools_cannot_be_masld_negative_controls(self) -> None:
        validate_contract(self.contract)
        self.assertFalse(self.contract["normal_stratum_is_masld_negative_control"])
        self.assertTrue(
            all(not row["masld_negative_control"] for row in self.contract["condition"])
        )
        changed = copy.deepcopy(self.contract)
        changed["condition"][0]["masld_negative_control"] = True
        with self.assertRaises(PoolingTopologyError):
            validate_contract(changed)

    def test_firewalls_cannot_be_relaxed(self) -> None:
        for key in (
            "automatic_activation",
            "sealed_outcomes_read",
            "biological_matrices_downloaded",
            "controlled_data_accessed",
            "pooled_libraries_are_donors",
            "cells_are_biological_replicates",
            "donor_safe_split_available",
            "external_masld_champion_eligible",
        ):
            changed = copy.deepcopy(self.contract)
            changed[key] = True
            with self.assertRaises(PoolingTopologyError):
                validate_contract(changed)

    def test_assay_pair_builder_rejects_a_missing_modality(self) -> None:
        expected, _ = validate_contract(self.contract)
        records = []
        for row in expected:
            for modality, gsm_key in (("ATAC", "atac_gsm"), ("RNA", "rna_gsm")):
                records.append(
                    {
                        "library_id": row["library_id"],
                        "modality": modality,
                        "gsm": row[gsm_key],
                        "condition_id": row["condition_id"],
                        "pool_number": row["pool_number"],
                        "biosample": f"SAMN_{row[gsm_key]}",
                        "srx": f"SRX_{row[gsm_key]}",
                    }
                )
        built = build_library_pairs(records, expected)
        self.assertEqual(len(built), 6)
        with self.assertRaises(PoolingTopologyError):
            build_library_pairs(records[:-1], expected)


if __name__ == "__main__":
    unittest.main()
