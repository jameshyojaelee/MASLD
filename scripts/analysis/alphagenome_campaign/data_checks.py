#!/usr/bin/env python3
"""Focused scientific invariants for endpoint identity and P4 classification."""
import unittest
import numpy as np
import data_p4 as p4


class ScientificInvariants(unittest.TestCase):
    def test_life_stage_requires_adult_metadata(self):
        self.assertEqual(p4.life_stage({"biosample_life_stage":"adult","biosample_type":"tissue"}),"adult_verified")
        self.assertEqual(p4.life_stage({"biosample_life_stage":"unknown","biosample_type":"tissue"}),"life_stage_unresolved")
        self.assertEqual(p4.life_stage({"biosample_life_stage":"embryonic","biosample_type":"in_vitro_differentiated_cells"}),"developmental")
        self.assertEqual(p4.life_stage({"biosample_life_stage":"adult","biosample_type":"cell_line"}),"life_stage_unresolved")

    def test_primary_liver_unknown_not_promoted_to_adult(self):
        var={"ontology_curie":["UBERON:0002107","UBERON:0001114","CL:0000182"],"biosample_name":["liver","right lobe of liver","hepatocyte"],"biosample_type":["tissue","tissue","in_vitro_differentiated_cells"],"biosample_life_stage":["unknown","adult","embryonic"],"histone_mark":["H3K27ac"]*3}
        groups=p4.group_columns("CHIP_HISTONE",var)
        np.testing.assert_array_equal(groups["liver_h3k27ac"],[0,1])
        np.testing.assert_array_equal(groups["liver_h3k27ac_adult_verified"],[1])
        np.testing.assert_array_equal(groups["liver_h3k27ac_developmental"],[2])
        np.testing.assert_array_equal(groups["liver_h3k27ac_life_stage_unresolved"],[0])

    def test_missing_archives_never_complete_a_family(self):
        counts={"U1_h3k27ac":{"present_not_integrity_validated":96459,"archive_incomplete":1}}
        result=p4.family_dispositions([],counts)
        self.assertIn("retrieval_incomplete",result[0]["reason"])
        self.assertTrue(all(r["status"]=="not_run" for r in result))


if __name__=="__main__":
    unittest.main()
