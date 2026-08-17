import unittest

import pandas as pd

from masld_cl.local9_reference_pilot import _lineage_counts, _reclassify


class TestLocal9ReferencePilot(unittest.TestCase):
    def test_added_donors_are_reclassified_as_reference_only(self):
        cells = pd.DataFrame({
            "donor_id": ["old", "new"],
            "strict_reference": [True, False],
            "primary_query": [False, True],
            "query_control": [False, True],
        })
        observed = _reclassify(cells, {"new"})
        self.assertTrue(bool(observed.loc[1, "strict_reference"]))
        self.assertFalse(bool(observed.loc[1, "primary_query"]))
        self.assertFalse(bool(observed.loc[1, "query_control"]))

    def test_lineage_counts_use_donors_as_donors(self):
        cells = pd.DataFrame({
            "donor_id": ["a", "a", "b"],
            "audit_cell_type": ["H", "H", "H"],
        })
        self.assertEqual(_lineage_counts(cells, {"a", "b"})["H"], {"cells": 3, "donors": 2})
