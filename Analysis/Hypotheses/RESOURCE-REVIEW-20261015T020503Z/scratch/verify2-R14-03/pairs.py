import csv, gzip, collections
FM="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
lab={r["study_name"]:(r["tier"],r["tier_label"]) for r in csv.DictReader(open(FM+"/config/gwas_trait_tier.tsv"),delimiter="\t")}
R=FM+"/results/susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz"
rows=list(csv.DictReader(gzip.open(R,"rt"),delimiter="\t"))
print("total pair rows",len(rows))
print("eligible values",collections.Counter(r["eligible"] for r in rows))
c=collections.defaultdict(lambda: collections.Counter())
gs=collections.defaultdict(set); gs_pos=collections.defaultdict(set); g_pos=collections.defaultdict(set)
for r in rows:
    t,l=lab.get(r["gwas"],("NA","NA"))
    k=l
    c[k]["all_pairs"]+=1
    gs[k].add((r["gwas"],r["ensembl"]))
    if r["eligible"]=="TRUE":
        c[k]["elig"]+=1
        try: pp=float(r["PP.H4"])
        except: pp=None
        if pp is None or pp!=pp: c[k]["elig_ppNA"]+=1
        elif pp>0.5:
            c[k]["elig_gt0.5"]+=1; gs_pos[k].add((r["gwas"],r["ensembl"])); g_pos[k].add(r["ensembl"])
for k in c:
    print(k, dict(c[k]), "gene_studies",len(gs[k]), "gene_studies_pos",len(gs_pos[k]), "genes_pos",len(g_pos[k]))
d=c["direct_MASLD"]; e=c["liver_enzyme"]
print("direct rate", d["elig_gt0.5"], "/", d["elig"], round(d["elig_gt0.5"]/d["elig"],4))
print("enzyme rate", e["elig_gt0.5"], "/", e["elig"], round(e["elig_gt0.5"]/e["elig"],4))
