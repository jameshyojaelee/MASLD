"""Five-gene quadrature accuracy check; entirely synthetic, not Gate B agreement."""
import json
import os
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import betaln,gammaln,logsumexp,expit
from scipy.stats import norm

from likelihood_v1 import Row,Model,conditional_draw


counts=(5,25,75,200,500)
rows=[]
for gene,count in enumerate(counts):
    for j in range(count):
        rows.append(Row(gene,f'{j:04d}','fixture',20+j%17,10,.6,.3,.1,.01,.02,.05,(-1)**(gene+j),
                        float(j%3-1),(float(j%5-2)/10,),float(j%7-3)/10,float(j%5-2)/10,10,2,10))
model=Model(rows);p=model.initial();p[:5]=[-1.2,-.7,.4,.9,1.3];p[10]=.08;p[11]=.04
rng=np.random.default_rng(20260923)
drawn=[]
for row in model.rows:
    eta,r,_,_=model.predictors(row,p)
    drawn.append(conditional_draw(row,eta,r,*rng.random(2)))
model=Model(drawn)
tau_max=.05


def direct_loglik(rs,u):
    """Reference direct betaln enumeration without recurrence or derivatives."""
    total=0.
    for row in rs:
        eta,r,_,_=model.predictors(row,p,u)
        m=expit(eta-row.orientation*row.omega)
        rho=1e-6+(1-2e-6)*expit(r)
        k=row.allowed
        parts=[]
        for mean,disp in ((m,rho),(row.e,row.phi),(1-row.e,row.phi)):
            t=(1-disp)/disp
            parts.append(gammaln(row.n+1)-gammaln(k+1)-gammaln(row.n-k+1)
                         +betaln(k+mean*t,row.n-k+(1-mean)*t)-betaln(mean*t,(1-mean)*t))
        mass=logsumexp(np.asarray(parts)+np.log(row.prior)[:,None],axis=0)
        total+=mass[row.a-k[0]]-logsumexp(mass)
    return float(total)


result={'seed':20260923,'rows_per_gene':counts,'tau_max':tau_max,'cases':[]}
for multiple in (.25,1.,4.):
    tau=multiple*tau_max
    per20={};per40={};per_quad={}
    ll20=model.marginal(p,tau*tau,gene_terms=per20)[0]
    ll40=model.marginal(p,tau*tau,40,gene_terms=per40)[0]
    reference=0.
    for gene in range(5):
        rs=[row for row in model.rows if row.gene==gene]
        base=direct_loglik(rs,0.)
        integral,error=quad(lambda u:np.exp(direct_loglik(rs,u)-base)*norm.pdf(u,scale=tau),
                            -12*tau,12*tau,epsabs=1e-12,epsrel=1e-12,limit=200)
        per_quad[gene]=float(np.log(integral)+base)
        reference+=per_quad[gene]
        assert abs(per20[gene]-per40[gene])<=1e-8,(gene,tau,per20[gene],per40[gene])
        assert abs(per20[gene]-per_quad[gene])<=1e-8*(1+abs(per_quad[gene])),(gene,tau,per20[gene],per_quad[gene])
    case={'tau':tau,'ll20':ll20,'ll40':ll40,'quad':reference,
          'GH20_GH40_abs_difference':abs(ll20-ll40),'quad_abs_difference':abs(ll20-reference),
          'per_gene_GH20':per20,'per_gene_GH40':per40,'per_gene_quad':per_quad}
    print(json.dumps(case),flush=True)
    result['cases'].append(case)
    assert abs(ll20-ll40)<=1e-8,case
    assert abs(ll20-reference)<=1e-8*(1+abs(reference)),case
output=Path(os.environ['CODEX_REC_OUTPUT']) / os.environ.get('CODEX_REC_QUADRATURE_FILE','quadrature_v1.json')
output.write_text(json.dumps(result,indent=2)+'\n')
