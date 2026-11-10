import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/"
reg=pd.read_csv(R+"config/gwas_registry.tsv",sep="\t")
studies=["2019_31311600_NAFLD_EUR","2020_32298765_NAFLD_EUR","2021_34128465_PDFF_EUR","2021_34841290_NAFLD_EUR","2021_34957434_PDFF_EUR","2022_36402844_PDFF_EUR","2023_36280732_NAFLD_UKBB_EUR","2023_36280732_NAFLD_deCode_EUR","2023_36280732_NAFLD_Intermountain_EUR","FinnGen_NAFLD","FinnGen_NASH"]
tot=0
for s in studies:
    p=R+reg.loc[reg.study_name==s,"leadsnps_path"].iloc[0]
    try: d=pd.read_csv(p,sep="\t")
    except Exception as e: print(s,"missing",e); continue
    n=0
    for c,b in d.sort_values(["CHR","BP"]).groupby("CHR"):
        n+=1+int((np.diff(b.BP.values)>1_000_000).sum())
    print(s,len(d),n, list(d.columns)); tot+=n
print("plain lead file EUR direct windows",tot)
