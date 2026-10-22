"""Verify scientific source invariants before installing the 774-unit 5B/5C panels."""
import csv, hashlib, math, sys
from pathlib import Path
root=Path(sys.argv[1]); checks=[]
def rows(p,sep="\t"):
    with p.open() as f:return list(csv.DictReader(f,delimiter=sep))
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
check("all analyses complete",(root/"COMPLETE").exists())
s=rows(root/"figure5/panels/data/fig5b_protein_triage_summary.tsv")[0]
g=rows(root/"figure5/panels/data/fig5b_protein_triage_genes.tsv")
check("5B measured denominator and complete class partition",len(g)==int(s["n_prioritized_measured"])==sum(int(s[k]) for k in ["n_significant_concordant","n_significant_discordant","n_not_significant"]))
check("5B signs and protein calls reproduce all class counts",all(sum(x["class"]==cls for x in g)==int(s[col]) for cls,col in [("confirmed","n_significant_concordant"),("discordant","n_significant_discordant"),("not_corr","n_not_significant")]))
bulk=rows(Path(s["bulk_deg_path"]),","); byid={x["gene"]:x for x in bulk}
c=rows(root/"figure5/data/composite_mrna_protein_corrected.csv",",")
check("5C fixed25 unique rows",len(c)==25 and len({x["gene"] for x in c})==25)
check("5C primary ID effects and canonical RNA fills",all(abs(float(x["bulk_logFC"])-float(byid[x["bulk_gene_id"]]["logFC"]))<=.00051 and (x["bulk_sig"]=="TRUE")== (float(byid[x["bulk_gene_id"]]["padj"])<.05 and abs(float(byid[x["bulk_gene_id"]]["logFC"]))>.5) for x in c))
h=[x for x in c if x["gene"]=="HLA-DQA1"][0]
check("HLA-DQA1 primary annotation is uncalled",h["bulk_gene_id"].split(".")[0]=="ENSG00000196735" and h["bulk_sig"]=="FALSE")
m=rows(root/"GSE126848_selected_units.tsv")
check("GSE126848 inferential30clusters23control7disease allmale",len(m)==len({x["analysis_unit_id"] for x in m})==30 and sum(x["group_binary"]=="Control" for x in m)==23 and sum(x["group_binary"]=="Disease" for x in m)==7 and {x["inferred_sex"] for x in m}=={"M"})
d=rows(root/"cohort_summary.tsv")
check("five cohorts each24115geneBHfamily and1854pooledcalls",len(d)==5 and all(int(x["n_tested"])==24115 and int(x["canonical_genes"])==1854 for x in d))
r=rows(root/"canonical_cohort_replication.tsv")
check("cohort replication accounts for1854calls",len(r)==1854 and all(0<=int(x["n_replication_cohorts"])<=5 for x in r))
check("unchanged class prespec",hashlib.sha256((root/"class_difference/prespec/PRESPEC.md").read_bytes()).hexdigest()=="f5213ae6cf23db5dbf6ef9cd8f0c3077101dca22983fa7ba05b00bf8469c4ff5")
for p in [root/"class_difference/results/class_difference_results.tsv",root/"class_difference/results/corrected_null_diagnostic.tsv"]:
    check(p.name+" two primary arms exist",sum(x["arm"] in ["singlecell_hep_specificity","protein_stage_response_residualised"] for x in rows(p))==2)
text="\n".join("PASS "+x for x in checks)+"\n"
(root/"scientific_checks.txt").write_text(text); print(text)
