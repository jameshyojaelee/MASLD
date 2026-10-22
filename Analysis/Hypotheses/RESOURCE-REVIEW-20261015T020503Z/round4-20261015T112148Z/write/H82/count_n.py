import csv, collections, gzip
S="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/"
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/"
man=list(csv.DictReader(open(S+"stage_extension_sample_manifest.tsv"),delimiter="\t"))
# roster per axis
for ax in ["fibrosis_adjacent","fibrosis","NAS_adjacent","NAS"]:
    u={}
    for r in man:
        if r["axis"]==ax: u[r["sample_id"]]=(r["dataset"],r["group"])
    c=collections.Counter(u.values())
    print(ax,"units",len(u)); 
    for k in sorted(c): print("  ",k,c[k])
nas_units={r["sample_id"] for r in man if r["axis"].startswith("NAS")}
fib_units={r["sample_id"] for r in man if r["axis"].startswith("fibrosis")}
print("NAS roster",len(nas_units),"fib roster",len(fib_units),"overlap",len(nas_units&fib_units))
cx={r["sample_id"]:r for r in csv.DictReader(open(R+"scratch/review-R03/complexity_774.tsv"),delimiter="\t")}
print("complexity_774 rows",len(cx),"NAS units missing",len(nas_units-set(cx)),"fib missing",len(fib_units-set(cx)))
print("NAS units group_binary",collections.Counter((cx[s]["dataset"],cx[s]["group_binary"]) for s in nas_units))
# released calls per contrast
calls=collections.Counter(); ns={}
seen=collections.defaultdict(dict)
with open(S+"stage_extension_all_gene_results.tsv") as f:
    for r in csv.DictReader(f,delimiter="\t"):
        k=(r["axis"],r["contrast"]); ns[k]=(r["n_reference"],r["n_comparison"])
        if float(r["FDR"])<0.05: calls[k]+=1
        if r["contrast"]=="NAS1_2_vs_NAS0" or r["contrast"] in("F0_to_F1","F1_vs_F0"):
            seen[r["contrast"]+"|"+r["axis"]][r["gene_id_versioned"]]=r["logFC"]
for k in ns: print(k,ns[k],"BH<0.05",calls[k])
a=seen["NAS1_2_vs_NAS0|NAS_adjacent"]; b=seen["NAS1_2_vs_NAS0|NAS"]
print("NAS1_2_vs_NAS0 axes identical logFC:",all(a[g]==b[g] for g in a), len(a))
a=seen["F0_to_F1|fibrosis_adjacent"]; b=seen["F1_vs_F0|fibrosis"]
print("F0_to_F1 vs F1_vs_F0 identical logFC:",all(a[g]==b[g] for g in a), len(a))
