import unittest

from masld_cl.adoption_decision import _minimum_held_protocol_harmony


class TestAdoptionDecision(unittest.TestCase):
    def test_worst_held_protocol_row_is_selected(self):
        artifact = {
            "held_studies": {
                "A": {"held_control_rows": [
                    {"testable": True, "model_kind": "Fibroblasts",
                     "protocol_harmony": {"improvement": -0.114}},
                    {"testable": False, "model_kind": "T cells"},
                ]},
                "B": {"held_control_rows": [
                    {"testable": True, "model_kind": "Hepatocytes",
                     "protocol_harmony": {"improvement": 0.2}},
                ]},
            }
        }
        self.assertEqual(
            _minimum_held_protocol_harmony(artifact),
            {"study": "A", "model_kind": "Fibroblasts", "improvement": -0.114},
        )
