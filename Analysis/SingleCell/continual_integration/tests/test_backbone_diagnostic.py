import unittest

from masld_cl.backbone_diagnostic import METHODS, _mapping
from masld_cl.contracts import ContractError


class TestBackboneDiagnostic(unittest.TestCase):
    def test_mapping_requires_exact_method_roster(self):
        values = [f"{method}=/{method}.json" for method in sorted(METHODS)]
        self.assertEqual(set(_mapping(values)), METHODS)
        with self.assertRaises(ContractError):
            _mapping(values[:-1])
