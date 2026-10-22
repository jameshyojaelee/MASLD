# Greedy distance clumping (+-500 kb) of genome-wide significant variants; approximate locus count only.
import pandas as pd, sys
O=sys.argv[1]
for f in ["MVP_NAFLD_EUR","MVP_ALT_EUR","MVP_AST_EUR","UKBB_ALT","UKBB_AST","UKBB_GGT"]:
    d=pd.read_csv(f"{O}/gws_{f}.tsv",sep="\t",header=None,names=["chr","pos","p"])
    d=d.sort_values("p"); leads=[]
    taken={}
    for r in d.itertuples():
        L=taken.setdefault(r.chr,[])
        if all(abs(r.pos-x)>500000 for x in L):
            L.append(r.pos); leads.append(r)
    print(f, "gws_variants", len(d), "distance_clumped_loci", len(leads))
