"""Gate L10-02: how many genome-wide EUR enzyme loci would remain after matching discovery
power to the direct-trait GWAS (lead |z| rescaled by sqrt(N_direct / N_enzyme) must still
reach 5.45, the two-sided z at p = 5e-8). N from GWAS/finemapping/config/gwas_registry.tsv."""
from pathlib import Path
import numpy as np, pandas as pd
S = Path(__file__).parent
N = {'MVP_ALT_EUR': 406803, 'MVP_AST_EUR': 389192, 'UKBB_ALT': 343850, 'UKBB_AST': 343850, 'UKBB_GGT': 343850}
neff = lambda nt, nc: 4 / (1 / nc + 1 / (nt - nc))
ND = {'MVP_NAFLD_EUR N_eff': neff(432709, 29342), 'FinnGen_NAFLD N_eff': neff(438857, 4614),
      'largest PDFF (2022_36402844) N': 44867}
for k, v in ND.items(): print(f'{k} = {v:,.0f}')
def leads(f):
    d = pd.read_csv(S / f'gws_{f}.tsv', sep='\t', header=None, names=['chrom','pos','p','a1','a2','beta','se'])
    d['az'] = (d.beta / d.se).abs(); d = d.sort_values('az', ascending=False)
    taken, out = {}, []
    for r in d.itertuples():
        L = taken.setdefault(r.chrom, [])
        if all(abs(r.pos - x) > 500000 for x in L): L.append(r.pos); out.append(r.az)
    return np.array(out)
print('study\tN\tclumped_loci\t' + '\t'.join(f'kept_vs_{k.split()[0]}' for k in ND))
for f, n in N.items():
    z = leads(f)
    print(f + '\t' + str(n) + '\t' + str(len(z)) + '\t' + '\t'.join(str(int((z * np.sqrt(nd / n) >= 5.45).sum())) for nd in ND.values()))
