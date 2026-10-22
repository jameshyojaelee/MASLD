# Judge p_informative spot checks (read-only on saved tables). Run with rnaseq env python, PYTHONNOUSERSITE=1.
import pandas as pd
from scipy.stats import spearmanr, chi2_contingency
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/"
# L02-01: morphometric fat fraction vs pathologist-note grade
m=pd.read_csv(R+"jobs/3111/morph_fat.tsv",sep="\t"); g=pd.read_csv(R+"scratch/gate-L02-01/note_grades.tsv",sep="\t")
print("L02-01 slides ok:",(m.status=="ok").sum(),"of",len(m))
d=m[m.status=="ok"].merge(g,on="donor").dropna(subset=["ord"])
print("L02-01 morph vs grade n",len(d),spearmanr(d.ord,d.fat_frac))
# L10-03: GSE213621 accession-order quartile vs fibrosis bin
a=pd.read_csv(R+"scratch/validity-L10-03/accession_order_bins.tsv",sep="\t")
ct=pd.crosstab(a.acc_q,a.bin); print(ct); c=chi2_contingency(ct); print("L10-03 chi2",round(c[0],2),"p",c[1])
# L01-01: GGT1 index
ix=pd.read_csv(R+"scratch/gate-L01-01/index_all_genes.tsv",sep="\t").set_index("gene")
print("L01-01 idx", ix.loc[["GGT1","ALB","APOA1","KRT19","SOX9"],"idx"].round(3).to_dict())
# L10-04: platelet pilot caQTL coloc rate
o=pd.read_csv(R+"scratch/validity-L10-04/pilot_strength.tsv",sep="\t")
for t in (0.5,0.8): print("L10-04 rate PP4>",t,(o.best_pp4>t).groupby(o.study).agg(["mean","size"]).to_dict())
