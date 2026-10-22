import csv
S="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/"
G={"COL1A1","COL1A2","COL3A1","LUM","THBS2","EFEMP1","FBLN5","LOX","GNMT","MAT1A","CYP2C19"}
with open(S+"stage_extension_all_gene_results.tsv") as f:
    for r in csv.DictReader(f,delimiter="\t"):
        if r["gene_name"] in G and r["axis"] in ("fibrosis","NAS"):
            print(r["contrast"],r["gene_name"],round(float(r["logFC"]),3),"%.2e"%float(r["FDR"]))
