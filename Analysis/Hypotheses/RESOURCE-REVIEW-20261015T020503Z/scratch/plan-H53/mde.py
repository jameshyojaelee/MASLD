from math import asin, sqrt, sin
z=1.959964+0.841621
def mde(n1,n2,p0,vif=1.0):
    h=z*sqrt(1/n1+1/n2)*sqrt(vif)
    return sin(asin(sqrt(p0))+h/2)**2-p0, h
for lab,n1,n2 in [("pooled all chr",279,127),("confirmatory excl chr1/19",238,101),("MVP_ALT_EUR vs platelet all",87,127),("MVP_ALT_EUR vs platelet excl chr1/19",75,101)]:
    for p0 in (0.16,0.30,0.46):
        d,h=mde(n1,n2,p0); d2,_=mde(n1,n2,p0,vif=1/0.7)
        print("%-38s n=%d vs %d p0=%.2f  h=%.3f  MDE=%.3f  MDE(VIF 1.43)=%.3f"%(lab,n1,n2,p0,h,d,d2))
def wilson_hw(n,p,z=1.959964):
    return z*sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
for n in (279,238):
    for p in (0.3,0.5):
        print("Wilson half-width n=%d p=%.1f: %.3f"%(n,p,wilson_hw(n,p)))
