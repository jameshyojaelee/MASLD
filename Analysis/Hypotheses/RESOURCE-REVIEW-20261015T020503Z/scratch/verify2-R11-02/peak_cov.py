import csv
from collections import defaultdict
f="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/consensus_peak_manifest.tsv"
iv=defaultdict(list); n=defaultdict(int)
with open(f) as h:
    for r in csv.DictReader(h,delimiter="\t"):
        if r["supported_gse244832"]=="TRUE" or r["supported_gse281367"]=="TRUE":
            iv[(r["lineage"],r["chrom"])].append((int(r["start0"]),int(r["end"]))); n[r["lineage"]]+=1
bp=defaultdict(int)
for (lin,ch),L in iv.items():
    L.sort(); cs,ce=L[0]
    for s,e in L[1:]:
        if s<=ce: ce=max(ce,e)
        else: bp[lin]+=ce-cs; cs,ce=s,e
    bp[lin]+=ce-cs
G=3.1e9
for lin in sorted(bp): print(lin,"peaks",n[lin],"merged_bp",bp[lin],"frac_of_3.1Gb",round(bp[lin]/G,4))
