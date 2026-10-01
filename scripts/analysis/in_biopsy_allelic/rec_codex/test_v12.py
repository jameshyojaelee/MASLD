"""Scientific invariants for Claude's dated v1.2 clarifications."""
import unittest
from unittest.mock import patch
import numpy as np
from design_v1 import Design,Participant
from tests_v1 import p2_tau


class ChecksV12(unittest.TestCase):
    def test_invalid_lambda_quantile(self):
        class Profile:
            calls=[]
            def fit_marginal(self,tau_max,start,fixed_v,**kwargs):
                self.calls.append(fixed_v)
                return {'converged':True,'loglik':-1.}
        model=Profile()
        marginal={'converged':True,'parameters':np.zeros(1),'loglik':0.,'v':0.}
        one_invalid=np.full((50,1),.8);one_invalid[7]=.03
        result=p2_tau(model,marginal,.05,one_invalid)
        self.assertAlmostEqual(result['sources'][0]['tau0'],.04)
        self.assertEqual(len(model.calls),1)
        several_invalid=one_invalid.copy();several_invalid[:3]=np.nan
        result=p2_tau(model,marginal,.05,several_invalid)
        self.assertEqual(result['p_tau'],1.)
        self.assertEqual(len(model.calls),1)

    def test_cross_cohort_cluster(self):
        people=[Participant('A0','A',0,(.1,.2),.01,20.,2),
                Participant('B0','B',1,(.3,.4),.01,20.,2),
                Participant('B1','B',2,(.5,.6),.01,20.,2)]
        design=Design(people,[],{'A0':'family','B0':'family'})
        self.assertEqual(design.cohort_clusters(),{'A':[['A0','B0']],'B':[['B1']]})
        sampled=design.resample({'A':[0],'B':[0]})
        self.assertEqual(len(sampled.participants),3)
        for cohort in ('A','B'):
            self.assertAlmostEqual(sampled.stage_centered[sampled.cohort==cohort].sum(),0.)
        self.assertEqual(sampled.clusters['A0__copy000000'],sampled.clusters['B0__copy000000'])

    def test_adjusted_r_squared(self):
        people=[Participant(str(i),'C',i%3,(float(i%3)+.2*(i%4),.1),.01,20.,2)
                for i in range(12)]
        design=Design(people,[])
        raw=design.stage_r_squared(adjusted=False)['C']
        expected=max(0.,1-(1-raw)*11/10)
        self.assertAlmostEqual(design.stage_r_squared()['C'],expected)

    def test_zero_row_stage_predictors(self):
        people=[Participant(str(i),'C',i%3,(.1,.2),.01,
                            float('nan') if i==0 else 20.,0 if i==0 else 2)
                for i in range(6)]
        design=Design(people,[])
        with patch('design_v1.ordinal_stage_probabilities') as fit:
            design.stage_law()
        x=fit.call_args.args[2]
        self.assertTrue(np.isfinite(x).all())
        np.testing.assert_allclose(x[1]-x[0],[np.log1p(20),np.log1p(2)])


if __name__=='__main__':unittest.main()
