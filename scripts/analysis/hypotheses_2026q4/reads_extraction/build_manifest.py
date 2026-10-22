#!/usr/bin/env python3
"""Join the V2 ten-cohort library list to the STAR BAMs on disk and write the
array manifest for extract_reads_one_library.py.

Libraries with a BAM go into library_manifest.tsv (largest BAM first, so the
array's long tasks start early). Libraries with no BAM go into
libraries_without_bam.tsv; they would need realignment from FASTQ.

Usage: python3 -I build_manifest.py V2_MANIFEST MANIFEST_776 RESULTS_ROOT OUT_DIR
"""
import csv
import os
import sys

v2_path, m776_path, root, out_dir = sys.argv[1:5]
os.makedirs(out_dir, exist_ok=True)
out_manifest = os.path.join(out_dir, "library_manifest.tsv")
out_missing = os.path.join(out_dir, "libraries_without_bam.tsv")
if os.path.exists(out_manifest) or os.path.exists(out_missing):
    sys.exit("refusing to overwrite an existing manifest")

with open(m776_path) as fh:
    m776 = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
with open(v2_path) as fh:
    v2 = list(csv.DictReader(fh, delimiter="\t"))

have, missing = [], []
for r in v2:
    lib, coh = r["library"], r["cohort"]
    d = os.path.join(root, coh, "alignments", "star", lib)
    bam = os.path.join(d, f"{lib}.Aligned.sortedByCoord.out.bam")
    rec = {
        "cohort": coh,
        "library": lib,
        "layout": r["layout"],
        "in_pooled_five": r["in_pooled_five"],
        "in_776_counted": "1" if lib in m776 else "0",
        "strandedness_776": m776[lib]["strandedness"] if lib in m776 else "NA",
        "bam_path": bam,
        "sj_path": os.path.join(d, f"{lib}.SJ.out.tab"),
        "rpg_path": os.path.join(d, f"{lib}.ReadsPerGene.out.tab"),
        "log_final_path": os.path.join(d, f"{lib}.Log.final.out"),
        "fastq1": r["fastq1"],
        "fastq2": r["fastq2"],
    }
    if os.path.exists(bam):
        rec["bam_bytes"] = str(os.path.getsize(bam))
        have.append(rec)
    else:
        rec["bam_bytes"] = "0"
        missing.append(rec)

have.sort(key=lambda x: -int(x["bam_bytes"]))
cols = ["array_index", "cohort", "library", "layout", "in_pooled_five", "in_776_counted", "strandedness_776",
        "bam_bytes", "bam_path", "sj_path", "rpg_path", "log_final_path"]
with open(out_manifest, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for i, rec in enumerate(have, start=1):
        rec["array_index"] = str(i)
        w.writerow([rec[c] for c in cols])
mcols = ["cohort", "library", "layout", "in_pooled_five", "in_776_counted", "fastq1", "fastq2"]
with open(out_missing, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(mcols)
    for rec in missing:
        w.writerow([rec[c] for c in mcols])

by = {}
for rec in have:
    k = rec["cohort"]
    by.setdefault(k, [0, 0, 0])
    by[k][0] += 1
    by[k][1] += int(rec["bam_bytes"])
    by[k][2] += rec["in_776_counted"] == "1"
print("cohort\tn_bam\tbam_GB\tn_in_776")
for k in sorted(by):
    print(f"{k}\t{by[k][0]}\t{by[k][1] / 1e9:.0f}\t{by[k][2]}")
print(f"total\t{len(have)}\t{sum(int(r['bam_bytes']) for r in have) / 1e9:.0f}\t{sum(r['in_776_counted'] == '1' for r in have)}")
print(f"without BAM: {len(missing)}")
