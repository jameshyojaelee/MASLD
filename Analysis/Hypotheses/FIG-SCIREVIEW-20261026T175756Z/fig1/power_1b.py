"""Power of the unstratified two-sided exact (Fisher conditional) test, N=3769, G=426, D=220, alpha 0.05."""
import numpy as np
from scipy.stats import hypergeom, nchypergeom_fisher
N, G, D = 3769, 426, 220
x = np.arange(0, min(G, D) + 1)
p0 = hypergeom.pmf(x, N, G, D)
pv = np.array([p0[p0 <= p0[i] * (1 + 1e-7)].sum() for i in range(len(x))])
rej = pv <= 0.05
print("rejection region: O <=", x[rej & (x < 25)].max(), " O >=", x[rej & (x > 25)].min())
for psi in [0.40, 0.44, 0.45, 0.46, 0.47, 0.5, 1.7, 1.75, 1.77, 1.78, 1.8, 1.85]:
    pw = nchypergeom_fisher.pmf(x, N, G, D, psi)[rej].sum()
    print(f"OR {psi:.2f} power {pw:.3f}")
