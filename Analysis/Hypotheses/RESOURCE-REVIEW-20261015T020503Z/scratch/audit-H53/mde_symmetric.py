# Audit H53: MDE for the symmetric arm definition (liver-only vs platelet-only loci), same method as plan-H53/mde.py
from math import asin, sqrt, sin
z = 1.959964 + 0.841621
def mde(n1, n2, p0, vif=1.0):
    h = z * sqrt(1/n1 + 1/n2) * sqrt(vif)
    return sin(asin(sqrt(p0)) + h/2)**2 - p0, h
for lab, n1, n2 in [("planned 238 vs 101", 238, 101), ("symmetric liver-only 185 vs 101", 185, 101), ("S3 floor 120 vs 60", 120, 60)]:
    for p0 in (0.16, 0.30, 0.46):
        d, h = mde(n1, n2, p0); d2, _ = mde(n1, n2, p0, vif=1/0.7)
        print("%-34s p0=%.2f h=%.3f MDE=%.3f MDE(VIF 1.43)=%.3f" % (lab, p0, h, d, d2))
