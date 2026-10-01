#!/usr/bin/env python3
"""Check signed identity and validation membership at the comparison boundary."""
from io import StringIO
import unittest

import pandas as pd

from model_final import load_adaptation


class AdaptationIdentity(unittest.TestCase):
    def setUp(self):
        self.authority=pd.DataFrame({"key":["14:10:A:C","14:20:G:T"],
            "beta_alt":[0.25,-0.5],"heldout_fold":[1,1],"block_1mb":["14:0","14:0"]})

    def table(self, header="variant_id", first="chr14:20:G:T", effect="-0.5"):
        return StringIO(f"{header}\tobserved_beta\tadapter_example\n"
                        f"{first}\t{effect}\t-0.4\nchr14:10:A:C\t0.25\t0.2\n")

    def test_named_and_historical_unnamed_identity_survive_reordering(self):
        for header in ("variant_id", ""):
            result=load_adaptation(self.table(header),self.authority)
            self.assertEqual(result.key.tolist(),["14:20:G:T","14:10:A:C"])
            self.assertEqual(result.beta_alt.tolist(),[-0.5,0.25])

    def test_missing_reversed_and_duplicate_alleles_rejected(self):
        for first in ("chr14:21:G:T","chr14:20:T:G","chr14:10:A:C"):
            with self.assertRaises(ValueError):
                load_adaptation(self.table(first=first),self.authority)

    def test_changed_effect_and_wrong_fold_rejected(self):
        with self.assertRaises(AssertionError):
            load_adaptation(self.table(effect="0.5"),self.authority)
        self.authority.loc[1,"heldout_fold"]=0
        with self.assertRaises(ValueError):
            load_adaptation(self.table(),self.authority)


if __name__ == "__main__":
    unittest.main()
