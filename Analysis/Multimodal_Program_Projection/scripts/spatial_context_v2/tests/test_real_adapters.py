#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_ROOT / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


builder = load_module("08_build_real_adapters.py", "real_adapter_builder_test")


class RealAdapterSemantics(unittest.TestCase):
    def test_unresolved_skips_never_encode_zero_as_a_biological_count(self):
        class Seal:
            registry_sha256 = "r" * 64
            ready_sha256 = "s" * 64

        for spec in (builder._gse_spec(Seal()), builder._geomx_spec(Seal())):
            self.assertEqual(spec["biological_unit"], "unknown_public_biological_unit")
            self.assertEqual(spec["biological_unit_resolution"], "unresolved")
            design = builder._empty_design(spec, "skipped_fixture")
            self.assertEqual(design["n_biological"], "")
            self.assertEqual(design["n_technical"], 0)

    def test_program_coverage_uses_unique_measured_symbols_and_original_l1(self):
        rows = [
            {"program_uid": "p1", "mapped_symbol": "A", "original_l1_weight": "0.2"},
            {"program_uid": "p1", "mapped_symbol": "B", "original_l1_weight": "0.3"},
            {"program_uid": "p1", "mapped_symbol": "B", "original_l1_weight": "0.3"},
            {"program_uid": "p1", "mapped_symbol": "NA", "original_l1_weight": "0.2"},
            {"program_uid": "p2", "mapped_symbol": "C", "original_l1_weight": "1.0"},
        ]
        observed = builder.measured_program_coverage(rows, {"A", "B"}, ("p1", "p2"))
        self.assertEqual(observed["p1"][0], 2)
        self.assertAlmostEqual(observed["p1"][1], 0.5)
        self.assertEqual(observed["p2"][0], 0)
        self.assertAlmostEqual(observed["p2"][1], 0.0)

    def test_source_native_context_never_becomes_program_rows(self):
        template = {
            "n_hep": "100",
            "n_mac": "20",
            "n_metmac": "2",
            "il32_mode": "detected_vs_not",
            "delta_high_minus_low_um": "-0.5",
            "rho_il32_vs_knn_cd74": "0.1",
            "il32high_closer": "True",
        }
        rows = []
        for slide, is_mash in (
            ("Leuven_1", "True"),
            ("Leuven_2", "False"),
            ("Leuven_3", "True"),
            ("Leuven_4", "True"),
        ):
            rows.append({**template, "slide": slide, "is_mash": is_mash})
        observed = builder.source_native_cosmx_rows(rows)
        self.assertEqual(len(observed), 4)
        self.assertEqual(sum(bool(row["include_direction_concordance"]) for row in observed), 3)
        self.assertTrue(all("program_uid" not in row for row in observed))
        self.assertTrue(all("pvalue" not in row for row in observed))

    def test_mash_direction_failure_fails_closed(self):
        template = {
            "n_hep": "100",
            "n_mac": "20",
            "n_metmac": "2",
            "il32_mode": "detected_vs_not",
            "delta_high_minus_low_um": "-0.5",
            "rho_il32_vs_knn_cd74": "0.1",
        }
        rows = []
        for slide, is_mash, closer in (
            ("Leuven_1", "True", "True"),
            ("Leuven_2", "False", "True"),
            ("Leuven_3", "True", "False"),
            ("Leuven_4", "True", "True"),
        ):
            rows.append(
                {
                    **template,
                    "slide": slide,
                    "is_mash": is_mash,
                    "il32high_closer": closer,
                }
            )
        with self.assertRaises(Exception):
            builder.source_native_cosmx_rows(rows)


if __name__ == "__main__":
    unittest.main()
