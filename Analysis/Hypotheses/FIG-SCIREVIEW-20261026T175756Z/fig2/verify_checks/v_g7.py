import pandas as pd, numpy as np
L="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/data/lead_snps/"
studies=["2019_31311600_NAFLD_EUR","2020_32298765_NAFLD_EUR","2021_34128465_PDFF_EUR","2021_34841290_NAFLD_EUR","2021_34957434_PDFF_EUR","2022_36402844_PDFF_EUR","2023_36280732_NAFLD_UKBB_EUR","2023_36280732_NAFLD_deCode_EUR","FinnGen_NAFLD","FinnGen_NASH","UKBB_ALT","UKBB_AST","UKBB_GGT"]
def windows(d):
    n=0
    for c,b in d.sort_values(["CHR","BP"]).groupby("CHR"):
        bp=b.BP.values; n+=1+int((np.diff(bp)>1_000_000).sum())
    return n
tot_all=tot_g=0
for s in studies:
    d=pd.read_csv(L+s+"_leadSNPs_tiered.tsv",sep="\t")
    g=d[d.index_p<5e-8]
    wa,wg=windows(d) if len(d) else 0, windows(g) if len(g) else 0
    print(s,len(d),wa,len(g),wg, d.tier.value_counts().to_dict())
    if not s.startswith("UKBB_A") and not s.startswith("UKBB_G"): tot_all+=wa; tot_g+=wg
print("EUR direct windows all",tot_all,"gws",tot_g)
