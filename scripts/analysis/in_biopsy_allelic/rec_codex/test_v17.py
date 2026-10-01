"""Checks of inferential coordinates and sampling rules specified in v1.7."""
import unittest
from dataclasses import replace
from unittest.mock import patch
import numpy as np
from scipy.special import logsumexp
from likelihood_v1 import Model,Row,conditional_draw,bb_terms,fixed_bb
from tests_v1 import mc_result


def rows():
    return [Row(g,str(i),'synthetic',30+i%11,7+i%6,.6,.3,.1,.01,.02,.1,(-1)**(g+i),
                float(i%3-1),(),float(i%5-2)/10,float(i%7-3)/10,10,2,10)
            for g in range(3) for i in range(40)]


class ChecksV17(unittest.TestCase):
    def test_zero_covariate_does_not_change_sandwich(self):
        base=Model(rows());extended=Model([replace(r,z=(0.,)) for r in rows()])
        p=base.initial();p[:3]=[.4,.8,-1.]
        q=np.insert(p,2*base.G+1,2.)
        np.testing.assert_allclose(base.evaluate(p)[0],extended.evaluate(q)[0],atol=1e-12)
        np.testing.assert_allclose(base.sandwich(p),extended.sandwich(q),atol=1e-10)
        p[2*base.G]=.9
        self.assertTrue(np.isfinite(base.sandwich(p)))

    def test_missing_gene_coordinate_is_retained(self):
        all_rows=rows();model=Model([r for r in all_rows if r.gene!=1],gene_count=3)
        p=model.initial();p[1]=1.5;p[4]=-2.
        ll,g,h=model.evaluate(p)
        self.assertTrue(np.isfinite(ll));np.testing.assert_array_equal(g[[1,4]],0)
        np.testing.assert_array_equal(h[[1,4]],0)
        self.assertFalse(model.informative()[1]);self.assertFalse(model.informative()[4])

    def test_fixed_offset_in_cold_start(self):
        original=[replace(r,gene=0,orientation=1,omega=0.,stage=0.,z=()) for r in rows()[:8]]
        baseline=Model(original).cold_start()
        offset=Model([replace(r,omega=.4) for r in original]).cold_start()
        np.testing.assert_allclose(offset[0]-baseline[0],.4,atol=1e-5)

    def test_warm_tie_and_higher_cold(self):
        model=Model(rows());p=model.initial()
        with patch.object(model,'cold_start',return_value=p),patch.object(model,'fit_single',
             side_effect=[{'loglik':10.,'converged':True},{'loglik':10.+5e-7,'converged':True}]):
            self.assertEqual(model.fit(p)['kept_start'],'warm')
        with patch.object(model,'cold_start',return_value=p),patch.object(model,'fit_single',
             side_effect=[{'loglik':10.,'converged':True},{'loglik':11.,'converged':True}]):
            self.assertEqual(model.fit(p)['kept_start'],'cold')

    def test_inverse_cdf_includes_equality(self):
        r=rows()[0]
        logp=bb_terms(r.allowed,r.n,.3-r.orientation*r.omega,-3.)[0]
        probability=np.exp(logp-logsumexp(logp))
        u=float(np.cumsum(probability)[2])
        parts=np.stack((logp,fixed_bb(r.allowed,r.n,r.e,r.phi),fixed_bb(r.allowed,r.n,1-r.e,r.phi)))
        weights=np.log(r.prior)+logsumexp(parts,axis=1)
        class_h=np.exp(weights[0]-logsumexp(weights))
        self.assertEqual(conditional_draw(r,.3,-3.,class_h/2,u).a,int(r.allowed[2]))

    def test_mc_standard_error_uses_draw_count(self):
        result=mc_result([True,False,True,False])
        self.assertAlmostEqual(result['mc_se'],np.sqrt(result['p']*(1-result['p'])/4))


if __name__=='__main__':unittest.main()
