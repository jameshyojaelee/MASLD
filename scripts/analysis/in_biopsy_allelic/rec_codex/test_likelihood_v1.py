import unittest
from dataclasses import replace
import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from scipy.stats import binom,norm
from likelihood_v1 import Row,Model,allowed_counts,conditional_terms,bb_terms,mixture_terms,sum_terms,bb_logmass,het_inclusion
from tests_v1 import p2_kappa,holm_three,rng_for
from design_v1 import Design,Participant
from lambda_v1 import marginal_reliability,conditional_reliability,PairedReliability


def fixture():
    return [Row(g,str(i),'synthetic',20+i,8+i//2,.6,.3,.1,.01,.02,.1,(-1)**g,
                float(i-2),(float(i%2)-.4,),.1*(i-2),.1*(i%3-1),10,2,10)
            for g in range(2) for i in range(5)]


class Checks(unittest.TestCase):
    def test_dispersion_floor_precision(self):
        from decimal import Decimal,localcontext
        from math import comb
        with localcontext() as context:
            context.prec=60
            for n,k,eta in ((40,7,.5),(100,20,-1.),(500,250,0.)):
                rho=1e-6+(1-2e-6)*expit(-15.);t=(1-rho)/rho
                aa=expit(eta)*t;bb=expit(-eta)*t
                a=Decimal(aa);b=Decimal(bb);concentration=a+b
                probability=Decimal(comb(n,k))
                for j in range(k):probability*=a+j
                for j in range(n-k):probability*=b+j
                for j in range(n):probability/=concentration+j
                expected=float(probability.ln())
                actual=bb_logmass(np.array([k]),n,aa,bb)[0]
                self.assertAlmostEqual(actual,expected,delta=1e-12)

    def test_relative_cluster_sandwich(self):
        rows=[Row(g,str(i),'synthetic',30+i%11,7+i%6,.6,.3,.1,.01,.02,.1,(-1)**(g+i),
                  float(i%3-1),(),float(i%5-2)/10,float(i%7-3)/10,10,2,10)
              for g in range(3) for i in range(40)]
        model=Model(rows);p=model.initial();p[:3]=[.4,.8,-1.]
        copied=[replace(row,individual=row.individual+suffix) for row in rows for suffix in ('a','b')]
        doubled=Model(copied)
        clusters={row.individual:row.individual[:-1] for row in copied}
        separate=doubled.sandwich(p)
        merged=doubled.sandwich(p,clusters)
        self.assertTrue(np.isfinite(separate) and separate>0)
        np.testing.assert_allclose(merged,np.sqrt(2)*separate,atol=1e-10,rtol=1e-8)
        np.testing.assert_allclose(model.sandwich(p),merged,atol=1e-10,rtol=1e-8)

    def test_label_reliability_invariants(self):
        lower_error=marginal_reliability(2.,1.2,.2)
        higher_error=marginal_reliability(2.,1.2,.4)
        self.assertGreater(lower_error,higher_error)
        self.assertEqual(conditional_reliability(lower_error,0),lower_error)
        self.assertLess(conditional_reliability(lower_error,.2),lower_error)
        self.assertLess(conditional_reliability(.2,.3),0)
        self.assertTrue(np.isnan(conditional_reliability(.7,1.)))
        people=[Participant(str(i),'synthetic',i%3,(float(i%2),float(i//3)),.01,25.,2) for i in range(6)]
        design=Design(people,fixture())
        callback=PairedReliability({'GSE193066':'S3'},{'GSE193066':[.2,.4]},
                                   {str(i):i%5 for i in range(6)},{'synthetic':1.},
                                   observed_r_squared=design.stage_r_squared())
        self.assertGreater(callback(design,0)[0],callback(design,1)[0])

    def test_normalizer_recurrence_and_switch(self):
        from scipy.special import betaln,gammaln,logsumexp
        for n,lo,eta,r in ((10,2,0.,-10),(100,20,3.,-3),(500,50,-1.,-5),(6821,683,.2,-4)):
            rho=1e-6+(1-2e-6)*expit(r);t=(1-rho)/rho
            k=np.arange(n+1);aa=expit(eta)*t;bb=expit(-eta)*t
            direct=gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+betaln(k+aa,n-k+bb)-betaln(aa,bb)
            np.testing.assert_allclose(bb_logmass(k,n,aa,bb),direct,atol=2e-10,rtol=1e-10)
            expected=sum_terms(*bb_terms(np.arange(lo,n-lo+1),n,eta,r))
            actual=het_inclusion(n,lo,eta,r)
            for a,b in zip(actual,expected):
                np.testing.assert_allclose(a,b,atol=2e-9,rtol=1e-8)
        for row in fixture():
            lp,s,h=mixture_terms(row,row.allowed,.4,-3)
            den,ds,dh=sum_terms(lp,s,h)
            j=row.a-row.allowed[0]
            actual=conditional_terms(row,.4,-3)
            for a,b in zip(actual,(lp[j]-den,s[j]-ds,h[j]-dh)):
                np.testing.assert_allclose(a,b,atol=1e-10,rtol=1e-8)

    def test_complete_participant_centering_and_resample(self):
        people=[Participant(str(i),'synthetic',i%3,(float(i),.1, np.nan if i==0 else float(i)),.01,25.,2) for i in range(4)]
        # Participant 3 has no included row but remains in every cohort mean.
        rows=[replace(r,individual=str(int(r.individual)%3)) for r in fixture() if int(r.individual)<3]
        design=Design(people,rows,{'0':'family','1':'family'})
        np.testing.assert_allclose(design.z.mean(axis=0),0,atol=1e-15)
        self.assertEqual(design.schema['indicators'],(2,))
        self.assertTrue(np.all(design.duplication==0))
        model=design.model()
        self.assertAlmostEqual(next(r.stage for r in model.rows if r.individual=='0'),-.75)
        self.assertAlmostEqual(next(r.background for r in model.rows if r.individual=='0'),-1.5)
        sampled=design.resample({'synthetic':np.array([2,2,1])})
        self.assertEqual(len(sampled.participants),5)
        self.assertEqual(len(set(sampled.clusters.values())),3)
        self.assertEqual(sampled.schema,design.schema)
        np.testing.assert_allclose(sampled.z.mean(axis=0),0,atol=1e-15)
        changed=design.with_stage([2,1,0,2])
        self.assertAlmostEqual(changed.stage_centered.sum(),0)

    def test_test_bookkeeping(self):
        result=p2_kappa([0.,.2,-.2,float('nan')],np.ones(4))
        self.assertEqual(result['p_kappa'],3/5)
        np.testing.assert_allclose(holm_three([.01,.04,1.]),[.03,.08,1.])
        np.testing.assert_array_equal(rng_for(1,3).random(10),rng_for(1,3).random(10))

    def test_integer_predicate(self):
        for n in range(2001):
            for d,k,f in ((10,2,10),(20,3,15),(30,5,20)):
                expected=np.array([a for a in range(n+1) if n>=d and min(a,n-a)>=k and 100*min(a,n-a)>=f*n])
                np.testing.assert_array_equal(allowed_counts(n,d,k,f),expected)

    def test_hand_example(self):
        self.assertAlmostEqual(binom.pmf(7,10,expit(.6)),.24973,delta=1e-4)
        row=Row(0,'example','synthetic',10,7,.7,.3,0.,.01,0.,.1,-1,0.,(),0.,0.,10,2,10)
        hetero=bb_terms(np.array([7]),10,.6,-15)[0][0]
        self.assertAlmostEqual(np.exp(hetero),.24973,delta=1e-4)
        prior=np.array([.42,.09,.49])
        np.testing.assert_allclose(row.prior,prior,atol=1e-15)
        k=np.arange(2,9)
        numerator=prior@np.array([binom.pmf(7,10,m) for m in (expit(.6),.01,.99)])
        denominator=prior@np.array([binom.pmf(k,10,m).sum() for m in (expit(.6),.01,.99)])
        self.assertAlmostEqual(conditional_terms(row,.5,-15)[0],np.log(numerator/denominator),delta=1e-4)

    def test_merged_prior_and_common_scale(self):
        row=fixture()[0]
        np.testing.assert_allclose(row.prior,[.36,.48,.15],atol=1e-15)
        opposite=replace(row,orientation=-1)
        np.testing.assert_allclose(opposite.prior,[.36,.15,.48],atol=1e-15)
        scaled=replace(row,p_ref=row.p_ref/2,p_alt=row.p_alt/2,p_x=row.p_x/2)
        self.assertAlmostEqual(conditional_terms(row,.4,-3)[0],conditional_terms(scaled,.4,-3)[0],places=12)

    def test_orientation(self):
        r=fixture()[0]
        rr=replace(r,a=r.n-r.a,orientation=-r.orientation)
        self.assertAlmostEqual(conditional_terms(r,.5,-3)[0],conditional_terms(rr,-.5,-3)[0],places=10)

    def test_gradient_hessian(self):
        m=Model(fixture());p=m.initial();p[:2]=[.4,-.5];p[4:]=[.1,.03,.04,.02,-.03,.01]
        ll,g,h=m.evaluate(p)
        step=1e-5
        for j in range(m.size):
            plus=p.copy();minus=p.copy();plus[j]+=step;minus[j]-=step
            lp,gp,_=m.evaluate(plus);lm,gm,_=m.evaluate(minus)
            np.testing.assert_allclose(g[j],(lp-lm)/(2*step),atol=1e-6,rtol=1e-5)
            np.testing.assert_allclose(h[:,j],(gp-gm)/(2*step),atol=2e-6,rtol=1e-4)

    def test_adaptive_integral_and_gradient(self):
        m=Model(fixture());p=m.initial();p[:2]=[.4,-.5];v=.05**2
        val,g=m.marginal(p,v)
        expected=0.
        for gene in range(m.G):
            rs=[r for r in m.rows if r.gene==gene]
            base=m.evaluate(p,rs)[0]
            mass=quad(lambda u:np.exp(m.evaluate(p,rs,u)[0]-base)*norm.pdf(u,scale=np.sqrt(v)),
                      -12*np.sqrt(v),12*np.sqrt(v),epsabs=1e-12)[0]
            expected+=np.log(mass)+base
        self.assertAlmostEqual(val,expected,places=8)
        self.assertAlmostEqual(val,m.marginal(p,v,40)[0],places=8)
        for j in (0,2*m.G,m.size):
            pars=np.r_[p,v];step=1e-6;plus=pars.copy();minus=pars.copy();plus[j]+=step;minus[j]-=step
            numeric=(m.marginal(plus[:-1],plus[-1])[0]-m.marginal(minus[:-1],minus[-1])[0])/(2*step)
            np.testing.assert_allclose(g[j],numeric,atol=2e-5,rtol=2e-4)

    def test_variance_boundary(self):
        m=Model(fixture());p=m.initial();p[:2]=[.4,-.5]
        ll,g=m.marginal(p,0)
        self.assertEqual(ll,m.evaluate(p)[0])
        numerical=(m.marginal(p,1e-6)[0]-ll)/1e-6
        np.testing.assert_allclose(g[-1],numerical,atol=1e-3,rtol=1e-3)

    def test_marginal_hessian(self):
        m=Model(fixture());p=m.initial();p[:2]=[.4,-.5];v=.01
        _,g,h=m.marginal(p,v,with_hessian=True)
        for j in (0,2*m.G,m.size):
            parameters=np.r_[p,v];step=1e-5 if j<m.size else 1e-6
            plus=parameters.copy();minus=parameters.copy();plus[j]+=step;minus[j]-=step
            gp=m.marginal(plus[:-1],plus[-1])[1];gm=m.marginal(minus[:-1],minus[-1])[1]
            np.testing.assert_allclose(h[:,j],(gp-gm)/(2*step),atol=2e-5,rtol=1e-4)


if __name__=='__main__':unittest.main()
