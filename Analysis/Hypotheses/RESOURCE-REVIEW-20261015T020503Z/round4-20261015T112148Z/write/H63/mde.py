# Fisher-z MDE for Spearman rho (normal approximation, 1.06 variance inflation for Spearman per Fieller et al. 1957), two-sided alpha 0.05, power 0.8
import numpy as np
from scipy.stats import norm, binom
za, zb = norm.ppf(0.975), norm.ppf(0.8)
for lab,n in [('T1 primary',1337),('anchored',61),('established',27),('unresolved',1208),('untested',40),('protein-zonated nonLM',707),('protein-strong nonLM',142),('LM PC-2',22)]:
    print(lab, n, 'rho_MDE', round(float(np.tanh((za+zb)*np.sqrt(1.06/(n-3)))),3))
print('P(5/5 same sign | p=0.5)', 2*0.5**5, ' one-sided 5/5', 0.5**5)
print('PC-2 binom P(>=18 of 22 | 0.5)', round(float(binom.sf(17,22,0.5)),5), ' P(>=15)', round(float(binom.sf(14,22,0.5)),4))
