# H82 audit: recount n, released calls, |t| strata, marker values, PC-2a, and list pre-existing outcome files.
import csv, gzip, collections, statistics, math, os
P="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
S=P+"figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/"
R=P+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/"
man=list(csv.DictReader(open(S+"stage_extension_sample_manifest.tsv"),delimiter="\t"))
cx={r["sample_id"]:r for r in csv.DictReader(open(R+"scratch/review-R03/complexity_774.tsv"),delimiter="\t")}
for ax in ["fibrosis","NAS"]:
    u={r["sample_id"]:(r["dataset"],r["group"]) for r in man if r["axis"]==ax}
    print(ax,"units",len(u),sorted(collections.Counter(u.values()).items()))
    print(ax,"group x control",sorted(collections.Counter((g,cx[s]["group_binary"]) for s,(d,g) in u.items()).items()))
t=collections.defaultdict(list); mk={}
G8=["COL1A1","COL1A2","COL3A1","LUM","THBS2","EFEMP1","FBLN5","LOX"]
with open(S+"stage_extension_all_gene_results.tsv") as f:
    for r in csv.DictReader(f,delimiter="\t"):
        if (r["axis"],r["contrast"]) in (("fibrosis_adjacent","F0_to_F1"),("NAS_adjacent","NAS1_2_vs_NAS0")): continue
        if float(r["FDR"])<0.05: t[r["contrast"]].append(abs(float(r["t"])))
        if r["gene_name"] in G8+["GNMT","MAT1A","CYP2C19"]: mk[(r["contrast"],r["gene_name"])]=(float(r["logFC"]),float(r["FDR"]))
for k,v in t.items():
    n=len(v); a=sum(x<4 for x in v); b=sum(4<=x<5 for x in v); c=sum(x>=5 for x in v)
    print(k,n,"median|t| %.4f min %.4f share<4 %.4f strata %d/%d/%d"%(statistics.median(v),min(v),a/n,a,b,c))
e=[mk[("NAS5_8_vs_NAS0",g)] for g in G8]
print("NAS5_8_vs_NAS0 ECM FDR range %.3g to %.3g; logFC %.3f to %.3f"%(min(x[1] for x in e),max(x[1] for x in e),min(x[0] for x in e),max(x[0] for x in e)))
for g in ["GNMT","MAT1A","CYP2C19"]: print("F1_vs_F0",g,"logFC %.4f FDR %.3g"%mk[("F1_vs_F0",g)])
ids={}
with open(S+"stage_extension_all_gene_results.tsv") as f:
    for r in csv.DictReader(f,delimiter="\t"):
        if r["contrast"]=="F4_vs_F0" and r["gene_name"] in G8: ids[r["gene_id_versioned"]]=r["gene_name"]
ok=0
for r in csv.DictReader(gzip.open(R+"jobs/1721/fits/joint__F4_vs_F0.genes.tsv.gz","rt"),delimiter="\t"):
    if r["gene_id_versioned"] in ids: ok+=(float(r["FDR"])<0.05 and float(r["logFC"])>0)
print("PC-2a 1721 joint F4: %d of %d"%(ok,len(ids)))
se=math.sqrt(1/103+1/102); print("alpha_BH %.6f; median d %.4f; min d %.4f"%(0.05*1338/24115,statistics.median(t["F1_vs_F0"])*se,min(t["F1_vs_F0"])*se))
print("pre-existing joint retention (review-G2):")
for r in csv.DictReader(open(R+"scratch/review-G2/refit_vs_released_counts.tsv"),delimiter="\t"): print("  ",r["label"],r["contrast"],r["refit"],r["retained"],r["retained_frac"])
print("pre-existing F3ef check counts:")
for r in csv.DictReader(open(R+"round4-20261015T112148Z/check/F3ef/refit-20261015T121237Z/fit_counts.tsv"),delimiter="\t"):
    if r["fit"] in ("lib","gini","nf") and r["contrast"] in ("F1_to_F2","F2_to_F3","F3_vs_F0"): print("  ",r["fit"],r["contrast"],r["n_bh05"])
