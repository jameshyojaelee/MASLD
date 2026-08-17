from __future__ import annotations

import unittest

from masld_cl.control_adapter_selection import _eligible_summary


class TestControlAdapterSelection(unittest.TestCase):
    def row(self, kind: str, improved: bool = True, protocol: bool = True):
        return {
            "model_kind": kind,
            "gates": {
                "reference_macro_f1": True,
                "reference_neighborhood": True,
                "uniform_offset_invariant": True,
                "unsorted_harmony_no_material_worsening": protocol,
                "unsorted_architecture_no_material_worsening": protocol,
                "harmony_improvement": improved,
                "architecture_improvement": improved,
            },
            "bootstrap": {"harmony": {"candidate_shift": 0.5}},
        }

    def test_four_of_five_lineages_and_all_lineage_are_required(self):
        config = {
            "lineages": ["H", "M", "F", "C", "T"],
            "gates": {"minimum_improved_lineages": 4},
        }
        rows = [self.row("all_lineage")] + [self.row(x) for x in config["lineages"]]
        rows[-1] = self.row("T", improved=False)
        self.assertTrue(_eligible_summary(rows, config)["eligible"])
        rows[-2] = self.row("C", improved=False)
        self.assertFalse(_eligible_summary(rows, config)["eligible"])
        rows[-2] = self.row("C", protocol=False)
        self.assertFalse(_eligible_summary(rows, config)["eligible"])


if __name__ == "__main__":
    unittest.main()
