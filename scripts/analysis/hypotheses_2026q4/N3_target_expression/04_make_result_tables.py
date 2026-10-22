#!/usr/bin/env python
"""N3 step 4: markdown tables for RESULT_N3.md, built from results/*.tsv and summary.json.

Standard library only. Usage: python -I 04_make_result_tables.py <RESULT_DIR>
Writes <RESULT_DIR>/results/result_tables.md. Each table sits between
<!-- BEGIN name --> and <!-- END name --> so RESULT_N3.md can be checked against it.
Also prints the counts quoted in the RESULT_N3.md text.
"""
import csv
import json
import sys
from pathlib import Path

RES = Path(sys.argv[1])
R = RES / "results"
CLASS_LABEL = {
    "stage_change_replicated_in_bulk": "replicated in bulk",
    "stage_change_not_replicated_in_bulk": "not replicated in bulk",
    "no_change_detected": "no change detected",
    "not_detected_in_hepatocytes": "not detected in hepatocytes",
}
STAGES = ["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"]


def read(name):
    with open(R / name, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def signed(x, digits=3):
    v = num(x)
    return "NA" if v is None else f"{v:+.{digits}f}"


def plain(x, digits=3):
    v = num(x)
    return "NA" if v is None else f"{v:.{digits}f}"


def pval(x):
    v = num(x)
    if v is None:
        return "NA"
    return f"{v:.1e}" if v < 0.001 else f"{v:.3f}"


def table(name, headers, rows):
    out = [f"<!-- BEGIN {name} -->", "| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    out.append(f"<!-- END {name} -->")
    return "\n".join(out)


def class_cell(r):
    label = CLASS_LABEL[r["class"]]
    if r["class"] == "stage_change_not_replicated_in_bulk" and r["bulk_same_sign"] == "FALSE":
        label += " (opposite sign)"
    if r["nonhepatocyte_dominant"] == "TRUE":
        label += "; non-hepatocyte-dominant"
    return label


cls = read("target_classes.tsv")
tested = sorted((r for r in cls if r["testable"] == "TRUE"), key=lambda r: num(r["P.Value"]))
undetected = sorted((r for r in cls if r["testable"] != "TRUE"), key=lambda r: r["symbol"])
ordered = tested + undetected
hits = [r for r in tested if num(r["q_bh"]) < 0.05]
summary = json.load(open(R / "summary.json"))
tables = []

tables.append(table(
    "T1a_hits_stage_trend",
    ["Target", "log2 change per step", "SE", "MDE80", "BH q (45 targets)", "Permutation BH q", "maxT adj. p",
     "Cirrhosis vs Healthy log2 (BH q, 135 tests)", "Leave-one-dataset-out same sign", "GSE202379 per F step log2 (BH q)"],
    [[r["symbol"], signed(r["logFC_per_step"]), plain(r["SE"]), plain(r["MDE80"]), pval(r["q_bh"]), pval(r["q_emp"]),
      pval(r["p_maxT"]), f'{signed(r["logFC_Cirrhosis"], 2)} ({pval(r["q_bh_Cirrhosis"])})',
      f'{r["S3_same_sign"]} of {r["S3_runs"]} refits',
      "not testable" if num(r["S2_logFC_per_F"]) is None else f'{signed(r["S2_logFC_per_F"])} ({pval(r["S2_q_bh"])})']
     for r in hits]))

tables.append(table(
    "T1b_hits_measurement_bulk",
    ["Target", "Hepatocytes with at least 1 UMI", "Hepatocyte share of target UMIs", "Bulk log2 (disease vs control)",
     "Bulk target-family BH q", "Bulk genome-wide padj", "Class"],
    [[r["symbol"], plain(r["hep_detect_rate"]), plain(r["hep_share_of_umis"]), signed(r["bulk_logFC"]),
      pval(r["bulk_q_targets"]), pval(r["bulk_genome_padj"]), class_cell(r)] for r in hits]))

class_rows = []
for key, label in CLASS_LABEL.items():
    members = [r for r in ordered if r["class"] == key]
    names = ", ".join(r["symbol"] + ("†" if r["nonhepatocyte_dominant"] == "TRUE" else "") for r in members)
    class_rows.append([label, len(members), names])
tables.append(table("T2_classes", ["Class (pre-specified rule)", "Targets", "Which"], class_rows))

tables.append(table(
    "T3_all_targets",
    ["Target", "log2 change per step", "SE", "MDE80", "BH q", "Bulk log2", "Bulk target-family BH q",
     "Hepatocytes with at least 1 UMI", "Hepatocyte share of UMIs", "Class"],
    [[r["symbol"], signed(r["logFC_per_step"]), plain(r["SE"]), plain(r["MDE80"]), pval(r["q_bh"]),
      signed(r["bulk_logFC"]), pval(r["bulk_q_targets"]), plain(r["hep_detect_rate"]),
      plain(r["hep_share_of_umis"]), class_cell(r)] for r in ordered]))

m2 = read("M2_categorical_targets.tsv")
m2_hits = sorted((r for r in m2 if num(r["q_bh"]) < 0.05), key=lambda r: num(r["P.Value"]))
tables.append(table(
    "T4_M2_contrasts_q_lt_0.05",
    ["Target", "Contrast vs Healthy", "log2 change", "SE", "BH q (135 tests)", "M1 BH q"],
    [[r["symbol"], r["contrast"], signed(r["logFC"]), plain(r["SE"]), pval(r["q_bh"]),
      pval(next(c["q_bh"] for c in cls if c["symbol"] == r["symbol"]))] for r in m2_hits]))

pc = read("positive_control_per_gene.tsv")
pcs = summary["positive_control"]
pc_rows = [[r["symbol"], signed(r["logFC"]), plain(r["SE"]), plain(r["t"], 2), pval(r["P.Value"])] for r in pc]
pc_rows.append(["composite of the 3 genes (OLS, one-sided p)", signed(pcs["estimate"]), "NA", plain(pcs["t"], 2),
                pval(pcs["p_one_sided_negative"])])
tables.append(table("T5_positive_control", ["Gene", "log2 change per step", "SE", "t", "p"], pc_rows))

dist = summary["permutation"]["n_BH_lt_0.05_perm_distribution"]
tables.append(table(
    "T6_permutation_hits",
    ["Targets with BH q < 0.05 in a permuted data set", "Permutations (of 1,000)"],
    [[k, v] for k, v in sorted(dist.items(), key=lambda kv: int(kv[0]))]))

dn = read("donor_counts_by_dataset_stage.tsv")
datasets = sorted({r["dataset"] for r in dn})
cnt = {(r["dataset"], r["stage"]): int(r["N"]) for r in dn}
don_rows = [[d] + [cnt.get((d, s), 0) for s in STAGES] + [sum(cnt.get((d, s), 0) for s in STAGES)] for d in datasets]
don_rows.append(["All"] + [sum(cnt.get((d, s), 0) for d in datasets) for s in STAGES] + [sum(cnt.values())])
tables.append(table("T7_donors", ["Dataset"] + STAGES + ["All"], don_rows))

with open(RES / "target_list.tsv", newline="") as fh:
    tl = {r["symbol"]: r for r in csv.DictReader(fh, delimiter="\t")}
focus = [r["symbol"] for r in hits] + ["KLB", "FGFR4", "GCGR", "GIPR", "GLP1R"]
tables.append(table(
    "T8_database_links",
    ["Target", "Drugs the database links to MASH (highest MASH stage)", "Co-targets listed for the same drug"],
    [[s, tl[s]["drugs"].replace(";", ", "), tl[s]["co_targets_listed_for_same_drug"].replace(";", ", ") or "none"]
     for s in focus]))

(R / "result_tables.md").write_text("\n\n".join(tables) + "\n")

# counts quoted in the text
fl = lambda r, k: num(r[k])
wm = [r for r in tested if fl(r, "hep_share_of_umis") >= 0.5 and fl(r, "hep_detect_rate") >= 0.05]
wm_bulk = [r for r in wm if fl(r, "bulk_logFC") is not None]
s2 = [r for r in tested if fl(r, "S2_q_bh") is not None]
print("testable", len(tested), "undetected", len(undetected), "hits", [r["symbol"] for r in hits])
print("well measured (share>=0.5 and detection>=0.05):", len(wm), "of which hits:", [r["symbol"] for r in wm if r in hits])
print("MDE80 range, well measured:", round(min(fl(r, "MDE80") for r in wm), 3), round(max(fl(r, "MDE80") for r in wm), 3))
print("well measured with bulk:", len(wm_bulk), "same sign:", sum(r["bulk_same_sign"] == "TRUE" for r in wm_bulk))
print("targets with a bulk row:", sum(fl(r, "bulk_logFC") is not None for r in cls), "of", len(cls),
      "; none:", [r["symbol"] for r in cls if fl(r, "bulk_logFC") is None])
print("S2 testable targets:", len(s2), "smallest S2 q:", min((fl(r, "S2_q_bh"), r["symbol"]) for r in s2))
print("S2 hits same sign as M1:", sum((fl(r, "S2_logFC_per_F") > 0) == (fl(r, "logFC_per_step") > 0) for r in hits if fl(r, "S2_logFC_per_F") is not None),
      "of", sum(fl(r, "S2_logFC_per_F") is not None for r in hits))
print("no-change targets:", sum(r["class"] == "no_change_detected" for r in cls),
      "of which hepatocyte share < 0.5:", sum(r["class"] == "no_change_detected" and r["nonhepatocyte_dominant"] == "TRUE" for r in cls))
print("n_BH>=1 permutations:", sum(v for k, v in dist.items() if int(k) >= 1), "; >=4:", sum(v for k, v in dist.items() if int(k) >= 4))
print("fraction nominal p<0.05 permuted: mean", summary["permutation"]["mean_fraction_nominal_p_lt_0.05"],
      "quantiles", summary["permutation"]["quantiles_fraction_nominal"])
