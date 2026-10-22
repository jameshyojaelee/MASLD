# Normal-approximation SE from each row's 95% block-bootstrap CI; Bonferroni over the 7 MPRA models in Fig 7B
import csv
from scipy.stats import norm
rows=[r for r in csv.DictReader(open('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/fig7_panels_rebuild-20260901T012245Z/tables/panel_7b_gains.tsv'),delimiter='\t') if r['task'].startswith('MPRA')]
k=len(rows); z=norm.ppf(1-0.05/(2*k))
print('family size',k,'bonferroni z',round(z,3))
for r in rows:
    g,lo,hi=float(r['gain']),float(r['lo']),float(r['hi']); se=(hi-lo)/(2*1.959964)
    p=2*(1-norm.cdf(abs(g)/se))
    print(f"{r['model'][:45]:45s} gain {g:+.4f} 95%CI [{lo:+.4f},{hi:+.4f}] approx p {p:.4f} bonf-lower {g-z*se:+.4f} MDE80 {2.802*se:.4f} gain/MDE {g/(2.802*se):.2f}")
