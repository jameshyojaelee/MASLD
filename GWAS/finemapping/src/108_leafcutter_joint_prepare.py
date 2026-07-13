#!/usr/bin/env python3
"""Prepare native joint LeafCutter disease model and LOCO sensitivity inputs."""
import csv, gzip, json, os, re

ROOT=os.environ.get("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PROD=os.path.join(ROOT,"GWAS/finemapping/results/seqfunc/disease_splicing/production")
OUT=os.path.join(PROD,"joint_native"); os.makedirs(OUT,exist_ok=True)
COHORTS=["GSE126848","GSE130970","GSE135251","GSE162694","GSE213621"]
rows=list(csv.DictReader(open(os.path.join(PROD,"production_manifest.tsv")),delimiter="\t"))
rows=[r for r in rows if r["dataset"] in COHORTS and r["group_binary"] in {"Control","Disease"}]
assert len(rows)==860, len(rows)
assert len({r["sample_id"] for r in rows})==len(rows)
assert all(os.path.isfile(r["junc_path"]) and os.path.getsize(r["junc_path"])>0 for r in rows)
for ds in COHORTS:
    z=[r for r in rows if r["dataset"]==ds]
    assert {r["group_binary"] for r in z}=={"Control","Disease"}

with open(os.path.join(OUT,"juncfiles.txt"),"w") as f:
    f.writelines(r["junc_path"]+"\n" for r in rows)

# Headerless by LeafCutter 2.0.3 contract. Dataset is categorical because its
# values are strings; sex/age are intentionally excluded due to >65% missingness.
def write_groups(path, z):
    with open(path,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        for r in z: w.writerow([r["sample_id"],r["group_binary"],r["dataset"]])
write_groups(os.path.join(OUT,"groups_dataset_fixed.tsv"),rows)
for ds in COHORTS:
    write_groups(os.path.join(OUT,f"groups_loco_{ds}.tsv"),[r for r in rows if r["dataset"]!=ds])

# Optional native cluster labels: five-column exon contract documented by
# leafcutter-ds --help (chrom,start,end,strand,gene_name).
gtf="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
pat=re.compile(r'gene_name "([^"]+)"')
with gzip.open(gtf,"rt") as src, gzip.open(os.path.join(OUT,"gencode_v49_exons.tsv.gz"),"wt") as dst:
    dst.write("chr\tstart\tend\tstrand\tgene_name\n")
    for line in src:
        if line.startswith("#"): continue
        x=line.rstrip().split("\t")
        if len(x)==9 and x[2]=="exon":
            m=pat.search(x[8]); gene=m.group(1) if m else "NA"
            dst.write("\t".join((x[0],x[3],x[4],x[6],gene))+"\n")

contract={
 "leafcutter_version":"2.0.3",
 "n_samples":len(rows),"cohorts":COHORTS,
 "primary":"joint LeafCutter clustering; native leafcutter-ds Disease vs Control adjusted for categorical dataset fixed effects",
 "covariates":["dataset"],
 "excluded_covariates":{"sex":"available for only 278/860 production samples","age":"available for only 221/860 production samples"},
 "replication":"existing independently clustered cohort-specific native LeafCutter runs",
 "sensitivity":"five leave-one-cohort-out native models using the joint cluster/count substrate",
 "noncanonical":"Script 107 bespoke Welch-PSI/metafor path is retired",
 "raw_count_pooling_note":"Joint assay-native likelihood is fit to sample-level intron counts with dataset adjustment; counts are not summed across cohorts."
}
with open(os.path.join(OUT,"analysis_contract.json"),"w") as f: json.dump(contract,f,indent=2)
print(f"Prepared {len(rows)} samples for joint native LeafCutter")
