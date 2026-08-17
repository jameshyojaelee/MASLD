from __future__ import annotations

import unittest

from masld_cl.promotion import GateResult, _decision_from_gates


def gate(gate_id: str, category: str, passed: bool) -> GateResult:
    return GateResult(gate_id, category, passed, 0, "fixture", "fixture")


class TestPromotionDecision(unittest.TestCase):
    def test_nonimproved_but_noninferior_is_supplement_only(self):
        results = [
            gate("integrity", "integrity", True),
            gate("alignment_improvement", "alignment", False),
            gate("maximum_protocol_stratum_worsening", "alignment", True),
        ]
        self.assertEqual(_decision_from_gates(results), "supplement_only")

    def test_material_alignment_worsening_rejects(self):
        results = [
            gate("integrity", "integrity", True),
            gate("alignment_improvement", "alignment", False),
            gate("maximum_protocol_stratum_worsening", "alignment", False),
        ]
        self.assertEqual(_decision_from_gates(results), "reject")


if __name__ == "__main__":
    unittest.main()
