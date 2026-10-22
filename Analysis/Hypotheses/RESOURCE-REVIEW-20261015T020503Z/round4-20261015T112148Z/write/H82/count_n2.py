import csv, collections, statistics
S="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/"
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/"
cx={r["sample_id"]:r for r in csv.DictReader(open(R+"scratch/review-R03/complexity_774.tsv"),delimiter="\t")}
man=list(csv.DictReader(open(S+"stage_extension_sample_manifest.tsv"),delimiter="\t"))
c=collections.Counter()
for ax in ["fibrosis","NAS"]:
    u={r["sample_id"]:r["group"] for r in man if r["axis"]==ax}
    c=collections.Counter((u[s],cx[s]["group_binary"]) for s in u)
    print(ax,sorted(c.items()))
t=collections.defaultdict(list)
with open(S+"stage_extension_all_gene_results.tsv") as f:
    for r in csv.DictReader(f,delimiter="\t"):
        if r["axis"] in ("fibrosis_adjacent",) and r["contrast"]=="F0_to_F1": continue
        if r["axis"]=="NAS_adjacent" and r["contrast"]=="NAS1_2_vs_NAS0": continue
        if float(r["FDR"])<0.05: t[r["contrast"]].append(abs(float(r["t"])))
print("contrast n_calls median|t| min|t| frac|t|<4 frac4-5 frac>=5")
for k,v in t.items():
    n=len(v); print(k,n,round(statistics.median(v),3),round(min(v),3),round(sum(x<4 for x in v)/n,3),round(sum(4<=x<5 for x in v)/n,3),round(sum(x>=5 for x in v)/n,3))
