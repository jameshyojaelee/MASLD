#!/usr/bin/env python3
"""Unit checks for the complete released-scooby lane audit."""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "audit_scooby_released_lane.py"
SPEC = importlib.util.spec_from_file_location("scooby_released_lane", MODULE)
assert SPEC and SPEC.loader
released = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(released)


class ScoobyReleasedLaneTest(unittest.TestCase):
    def setUp(self) -> None:
        self.path = (
            ROOT
            / "config"
            / "artifacts"
            / "models"
            / "scooby_neurips"
            / "released_variants.json"
        )
        self.authority = json.loads(self.path.read_text(encoding="utf-8"))

    def test_two_unregistered_release_variants_fail_closed(self) -> None:
        variants = released.validate_variants(self.authority)
        self.assertEqual(set(variants), set(released.VARIANTS))
        for value in variants.values():
            self.assertFalse(value["registered"])
            self.assertFalse(value["checkpoint_forward_allowed"])
            self.assertFalse(value["rna_conditioned_atac_eligible"])
            self.assertFalse(value["sealed_masld_champion_eligible"])

    def test_rna_only_name_does_not_imply_rna_conditioned_atac(self) -> None:
        value = self.authority["variants"]["scooby_neurips_rna"]
        self.assertEqual(value["decoder_tracks"], ["rna_plus", "rna_minus"])
        self.assertIn("observed same-nucleus RNA and ATAC", value["held_query_topology"])
        self.assertFalse(value["rna_conditioned_atac_eligible"])

    def test_variant_topology_or_identity_tamper_fails_closed(self) -> None:
        for field, value in (
            ("rna_conditioned_atac_eligible", True),
            ("declared_weight_sha256", "0" * 64),
            ("revision", "main"),
        ):
            with self.subTest(field=field):
                changed = copy.deepcopy(self.authority)
                changed["variants"]["scooby_neurips_flash"][field] = value
                with self.assertRaises(released.ScoobyReleasedLaneError):
                    released.validate_variants(changed)

    def test_checkpoint_roster_is_three_registered_plus_two_variants(self) -> None:
        self.assertEqual(len(released.REGISTERED), 3)
        self.assertEqual(len(released.VARIANTS), 2)
        self.assertFalse(set(released.REGISTERED) & set(released.VARIANTS))


if __name__ == "__main__":
    unittest.main()
