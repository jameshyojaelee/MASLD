"""audit-H58: which loci were fine-mapped for UKBB ALT/AST/GGT, and are the PNPLA3/TM6SF2 regions among them?
Counts and positions only; no coding fraction."""
import pandas as pd
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
cs=pd.read_csv(ROOT+"GWAS/finemapping/results/credible_sets.csv",usecols=["chromosome","position","study","locus","susie_cs"],dtype={"locus":str,"susie_cs":str})
u=cs[cs.study.isin(["UKBB_ALT","UKBB_AST","UKBB_GGT"])]
print(u.groupby(["study","locus"]).agg(chrom=("chromosome","first"),mn=("position","min"),mx=("position","max"),n=("position","size")).to_string())
for name,(c,p) in {"PNPLA3":(22,44324727),"TM6SF2":(19,19379549)}.items():
    h=cs[(cs.chromosome==c)&((cs.position-p).abs()<500000)]
    print(name,"credible_sets.csv rows within 500kb by study:",h.groupby("study").size().to_dict())
exc=pd.read_csv(ROOT+"figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/excluded_credible_sets.tsv",sep="\t")
print("gate-excluded sets in UKBB studies:",exc[exc.study.astype(str).str.startswith("UKBB_")].groupby("study").size().to_dict())
for t in ["ALT","AST","GGT"]:
    l=pd.read_csv(ROOT+f"GWAS/finemapping/data/lead_snps/UKBB_{t}_leadSNPs_tiered.tsv",sep="\t")
    print(f"UKBB_{t} lead rows",len(l),"tiers",l.tier.value_counts().to_dict())
