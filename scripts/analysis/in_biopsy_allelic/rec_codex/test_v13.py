import unittest
import numpy as np
from design_v1 import Design,Participant
from lambda_v1 import PairedReliability


class ChecksV13(unittest.TestCase):
    def test_observed_r_squared_is_fixed(self):
        people=[Participant(str(i),'C',i%3,(float(i%2),float(i//3)),.01,20.,2)
                for i in range(12)]
        design=Design(people,[])
        args=({'source':'Y5'},{'source':[.2]}, {str(i):i%5 for i in range(12)},{'C':1.})
        with self.assertRaises(ValueError):PairedReliability(*args)
        fixed=PairedReliability(*args,observed_r_squared=design.stage_r_squared())
        changed=design.with_stage(np.zeros(12,int))
        np.testing.assert_allclose(fixed(changed,0),fixed(design,0))
        sensitive=PairedReliability(*args,recompute_r_squared=True)
        self.assertTrue(np.isnan(sensitive(changed,0)).all())


if __name__=='__main__':unittest.main()
