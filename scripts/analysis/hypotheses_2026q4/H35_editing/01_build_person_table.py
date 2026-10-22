#!/usr/bin/env python3
"""H35 step 1: one row per person (primary-selected roster unit) with the Alu editing
counts from the reads extraction and the log2CPM of ADAR, ADARB1 and the fixed ISG list
from STAR ReadsPerGene. Standard library only.

Usage: 01_build_person_table.py OUT_DIR
Reads PRESPEC_H35.md decisions: ALU rows of editing_summary.tsv, library_qc.tsv, run_info.json.
"""
import csv
import gzip
import json
import math
import os
import sys

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
READS = f"{ROOT}/Analysis/Hypotheses/READS-20261009T141017Z/libraries"
ROSTER = (f"{ROOT}/figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/"
          "all_screened_roster.tsv")
GENES = ("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/"
         "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/"
         "gencode_v49_gene_metadata.tsv.gz")
ISG = ["ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "MX1", "MX2",
       "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "XAF1", "USP18", "HERC5", "CMPK2"]
ADAR_ID, ADARB1_ID = "ENSG00000160710", "ENSG00000197381"

out_dir = sys.argv[1]
os.makedirs(out_dir, exist_ok=False)

# symbol -> ensembl base id (GENCODE v49)
sym2id, id2sym = {}, {}
with gzip.open(GENES, "rt") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        sym2id.setdefault(r["gene_name"], []).append(r["ensembl_base"])
        id2sym[r["ensembl_base"]] = r["gene_name"]
isg_ids, isg_problems = {}, []
for g in ISG:
    ids = sym2id.get(g, [])
    if len(ids) != 1:
        isg_problems.append((g, ids))
    else:
        isg_ids[g] = ids[0]
want = {ADAR_ID: "ADAR", ADARB1_ID: "ADARB1"}
want.update({v: k for k, v in isg_ids.items()})
assert id2sym.get(ADAR_ID) == "ADAR" and id2sym.get(ADARB1_ID) == "ADARB1", "ADAR ids do not match GENCODE v49 metadata"


def log2cpm(count, libsize):
    return math.log2((count + 0.5) / (libsize + 1.0) * 1e6)


with open(ROSTER) as fh:
    roster = [r for r in csv.DictReader(fh, delimiter="\t") if r["primary_selected"] == "True"]

cols = (["sample_id", "cohort", "analysis_unit_id", "stage", "stage_resolution", "source_fibrosis_category",
         "group_binary", "age", "inferred_sex", "expr_sex", "reported_sex", "strand_call",
         "alu_ref", "alu_alt", "alu_sites_any", "alu_sites_cov10",
         "alu_ts_ref", "alu_ts_alt", "other_base_alu", "other_base_nonalu", "rpg_libsize", "adar_count", "adarb1_count"]
        + [f"lcpm_{k}" for k in ["ADAR", "ADARB1"] + ISG])
n_missing = 0
rows = []
for r in roster:
    cohort, lib = r["dataset"], r["sample_id"]
    d = f"{READS}/{cohort}/{lib}"
    try:
        es = {}
        with open(f"{d}/editing_summary.tsv") as fh:
            for x in csv.DictReader(fh, delimiter="\t"):
                es[(x["site_class"], x["mode"])] = x
        with open(f"{d}/library_qc.tsv") as fh:
            qc = next(csv.DictReader(fh, delimiter="\t"))
        with open(f"{d}/run_info.json") as fh:
            ri = json.load(fh)
        rpg = ri["args"]["rpg"]
        libsize, counts = 0, {}
        with open(rpg) as fh:
            for line in fh:
                g, c1, _, _ = line.rstrip("\n").split("\t")
                if g.startswith("N_"):
                    continue
                c = int(c1)
                libsize += c
                b = g.split(".")[0]
                if b in want:
                    counts[b] = c
    except (FileNotFoundError, StopIteration, KeyError) as e:
        n_missing += 1
        print("MISSING", cohort, lib, repr(e), file=sys.stderr)
        continue
    a = es[("ALU", "all_fragments")]
    ts = es.get(("ALU", "transcript_strand"))
    row = {
        "sample_id": lib, "cohort": cohort, "analysis_unit_id": r["analysis_unit_id"],
        "stage": r["fibrosis_stage"], "stage_resolution": r["stage_resolution"],
        "source_fibrosis_category": r["source_fibrosis_category"], "group_binary": r["group_binary"],
        "age": r["age"], "inferred_sex": r["inferred_sex"], "expr_sex": r["expr_sex"],
        "reported_sex": r["reported_sex"], "strand_call": qc["strand_call"],
        "alu_ref": a["ref_reads"], "alu_alt": a["alt_reads"], "alu_sites_any": a["n_sites_any"],
        "alu_sites_cov10": a["n_sites_cov10"],
        "alu_ts_ref": ts["ref_reads"] if ts else "NA", "alu_ts_alt": ts["alt_reads"] if ts else "NA",
        "other_base_alu": qc["editing_other_base_alu"], "other_base_nonalu": qc["editing_other_base_nonalu"],
        "rpg_libsize": libsize, "adar_count": counts.get(ADAR_ID, 0), "adarb1_count": counts.get(ADARB1_ID, 0),
    }
    for gid, sym in want.items():
        row[f"lcpm_{sym}"] = f"{log2cpm(counts.get(gid, 0), libsize):.6f}"
    rows.append(row)

with open(f"{out_dir}/person_table.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)
with open(f"{out_dir}/build_notes.json", "w") as fh:
    json.dump({"n_roster_primary": len(roster), "n_rows": len(rows), "n_missing_files": n_missing,
               "isg_genes_used": sorted(isg_ids), "isg_problems": isg_problems,
               "roster": ROSTER, "reads": READS, "python": sys.version}, fh, indent=1)
print("rows", len(rows), "missing", n_missing, "isg problems", isg_problems)
