#!/usr/bin/env python3
"""
Build the human kallisto requant driver: the 3 paired-end mega cohorts
(GSE213621, GSE135251, GSE130970) reusing the verified FASTQ paths from the
rna-validation worktree driver. Asserts every R1/R2 exists.

Output: human_driver.tsv  cols: cohort, sample_id, layout, fastq_r1, fastq_r2
"""
import csv, os, sys

WT_DRIVER = ("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/"
             "rna-validation/RNA-seq/scripts/kallisto/sample_driver.tsv")
PE_COHORTS = {"GSE213621", "GSE135251", "GSE130970"}
OUT = ("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
       "RNA-seq/scripts/isoform_diversity/human_driver.tsv")

rows, missing = [], []
with open(WT_DRIVER) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        if row["cohort"] not in PE_COHORTS:
            continue
        r1, r2 = row["fastq_r1"].strip(), row["fastq_r2"].strip()
        for p in (r1, r2):
            if not p or not os.path.exists(p):
                missing.append((row["sample_id"], p))
        rows.append((row["cohort"], row["sample_id"], "PAIRED", r1, r2))

if missing:
    sys.stderr.write(f"[ERROR] {len(missing)} missing FASTQs, e.g. {missing[:3]}\n")
    sys.exit(1)

with open(OUT, "w") as out:
    out.write("cohort\tsample_id\tlayout\tfastq_r1\tfastq_r2\n")
    for r in rows:
        out.write("\t".join(r) + "\n")

from collections import Counter
c = Counter(r[0] for r in rows)
sys.stderr.write(f"[build_human_driver] {len(rows)} samples: {dict(c)}\n")
sys.stderr.write(f"[build_human_driver] wrote {OUT}\n")
