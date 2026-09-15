from __future__ import annotations

from unittest import mock
import unittest

from masld_bench.evaluators import rna_atac_cobolt_5seed_development as evaluator_v1
from masld_bench.evaluators import rna_atac_cobolt_5seed_v6_development as evaluator_v6


class RNAATACCoboltFiveSeedV6EvaluatorTests(unittest.TestCase):
    def test_entry_point_changes_only_campaign_binding_and_restores_base(self) -> None:
        with mock.patch.object(evaluator_v1, "main", return_value=17) as base_main:
            self.assertEqual(evaluator_v6.main(), 17)
        base_main.assert_called_once_with()
        self.assertEqual(evaluator_v1.CAMPAIGN_ID, evaluator_v6.SOURCE_CAMPAIGN_ID)

    def test_entry_point_rejects_changed_base_binding(self) -> None:
        with mock.patch.object(evaluator_v1, "CAMPAIGN_ID", "unexpected"):
            with self.assertRaisesRegex(RuntimeError, "campaign binding changed"):
                evaluator_v6.main()


if __name__ == "__main__":
    unittest.main()
