#!/usr/bin/env python3
"""H50 steps 1 and 2: PXD051911 patient and sample counts, then the detection gate.

Usage: python3 -I 01_counts_and_gate.py <new_output_dir>

Reads, as data only: PXD051911 meta_data.txt, plasma and liver quant files, the deposited filtered
plasma matrix, the frozen program membership, and the Olink assay list. Standard library only.
No quantity that relates a protein to fibrosis stage is computed here.
"""
import csv
import json
import math
import os
import sys
from collections import Counter

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
PX = f"{ROOT}/data/PXD051911"
MEMBERSHIP = f"{ROOT}/Analysis/Multimodal_Program_Projection/results/frozen_program_membership.tsv"
OLINK = f"{ROOT}/Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt"
PROGRAMS = {
    "hepatocytes::20": ("Ductular injury (BICC1)", "primary"),
    "fibroblasts::6": ("Biliary-like (PKHD1/ANXA4)", "context_only"),
}
GATE_MIN_PROTEINS = 5
GATE_MAX_MISSING = 0.5
GATE_MIN_PATIENTS = 40
GATE_MIN_PER_CLASS = 10
# Counts stated in Boel 2025 (PMID 40032974): text and Table 1.
PAPER = {
    "liver_proteomes": 58, "scWAT": 58, "oWAT": 27, "plasma_initial": 143,
    "plasma_followup": 41, "followup_no_surgery": 28, "followup_surgery": 13,
}


def na(v):
    return v is None or v.strip() in ("", "NA", "NaN", "nan")


def num(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def genes_of(cell):
    return {g.strip() for g in (cell or "").split(";") if g.strip()}


def main(out_dir):
    if os.path.exists(out_dir):
        sys.exit(f"refusing to overwrite {out_dir}")
    os.makedirs(out_dir)

    # ---- metadata counts -------------------------------------------------
    meta = read_tsv(f"{PX}/meta_data.txt")
    rows = []
    c = Counter()
    for r in meta:
        c["rows"] += 1
        c[f"group:{r['sample_group']}"] += 1
        if not na(r["liver_proteomics_filename"]):
            c["liver_files"] += 1
            c[f"liver_group:{r['sample_group']}"] += 1
        if not na(r["plasma_proteomics_filename"]):
            c["plasma_files"] += 1
            c[f"plasma_group:{r['sample_group']}"] += 1
        if not na(r["plasma_follow_up_proteomics_filename"]):
            c["followup_plasma_files"] += 1
        if not na(r["oWAT_proteomics_filename"]):
            c["oWAT_files"] += 1
        if not na(r["scWAT_proteomics_filename"]):
            c["scWAT_files"] += 1
        if not na(r["kleiner_fibrosis_grade"]):
            c[f"kleiner_group:{r['sample_group']}"] += 1
        rows.append(r)
    patients = {r["patient_name"] for r in rows}
    liver_pat = Counter(r["patient_name"] for r in rows if not na(r["liver_proteomics_filename"]))
    by_pat = {}
    for r in rows:
        by_pat.setdefault(r["patient_name"], []).append(r)
    graded = {p: [r for r in rs if not na(r["kleiner_fibrosis_grade"])] for p, rs in by_pat.items()}
    n_ge2 = sum(1 for g in graded.values() if len(g) >= 2)
    both = only_v2 = only_fu = 0
    v2_age, fu_age = Counter(), Counter()
    for p, g in graded.items():
        if len(g) < 2:
            continue
        grp = {x["sample_group"] for x in g}
        has_v2, has_fu = "at_bariatric_surgery_sample" in grp, "follow_up_sample" in grp
        if has_v2 and has_fu:
            both += 1
        elif has_v2:
            only_v2 += 1
        elif has_fu:
            only_fu += 1
        base = [x for x in g if x["sample_group"] == "initial_sample"]
        if not base:
            continue
        a0 = num(base[0]["alder"])
        for x in g:
            a = num(x["alder"])
            if a0 is None or a is None:
                continue
            if x["sample_group"] == "at_bariatric_surgery_sample":
                v2_age[a - a0] += 1
            elif x["sample_group"] == "follow_up_sample":
                fu_age[a - a0] += 1
    patients_without_initial = sorted(p for p, rs in by_pat.items()
                                      if not any(x["sample_group"] == "initial_sample" for x in rs))

    # S1 sample set: initial visit, plasma file, Kleiner grade, present in the deposited filtered matrix
    with open(f"{PX}/plasma_protein_quant_filtered_and_batch_corrected.txt", encoding="utf-8") as fh:
        filt_header = fh.readline().rstrip("\n").split("\t")
    filt_cols = set(filt_header[1:])
    s1 = [r for r in rows if r["sample_group"] == "initial_sample" and not na(r["plasma_proteomics_filename"])
          and not na(r["kleiner_fibrosis_grade"])]
    s1_in = [r for r in s1 if r["unique_identifier"] in filt_cols]
    grade = Counter(r["kleiner_fibrosis_grade"] for r in s1_in)
    n_lt2 = sum(v for k, v in grade.items() if k in ("F0", "F1"))
    n_ge2_s1 = sum(v for k, v in grade.items() if k in ("F2", "F3", "F4"))
    covars = ["alder", "gender", "bmi", "plasma_batch_effects"]
    n_complete = sum(1 for r in s1_in if all(not na(r[k]) for k in covars))
    g2_pass = len(s1_in) >= GATE_MIN_PATIENTS and n_lt2 >= GATE_MIN_PER_CLASS and n_ge2_s1 >= GATE_MIN_PER_CLASS

    counts = [
        ("rows in meta_data.txt", c["rows"], ""),
        ("distinct patients", len(patients), ""),
        ("rows initial_sample", c["group:initial_sample"], "paper: first visit"),
        ("rows at_bariatric_surgery_sample", c["group:at_bariatric_surgery_sample"], ""),
        ("rows follow_up_sample", c["group:follow_up_sample"], "paper: 2-year follow-up"),
        ("patients without an initial_sample row", len(patients_without_initial), ",".join(patients_without_initial)),
        ("liver proteome files", c["liver_files"], f"paper {PAPER['liver_proteomes']}"),
        ("liver proteome files, initial_sample rows", c["liver_group:initial_sample"], ""),
        ("patients with more than one liver proteome", sum(1 for v in liver_pat.values() if v > 1), ""),
        ("plasma proteome files", c["plasma_files"], f"paper {PAPER['plasma_initial']}"),
        ("plasma files, initial_sample", c["plasma_group:initial_sample"], ""),
        ("plasma files, at_bariatric_surgery_sample", c["plasma_group:at_bariatric_surgery_sample"], ""),
        ("follow-up plasma files", c["followup_plasma_files"], f"paper {PAPER['plasma_followup']} "
         f"({PAPER['followup_no_surgery']} without and {PAPER['followup_surgery']} with bariatric surgery)"),
        ("scWAT files", c["scWAT_files"], f"paper {PAPER['scWAT']}"),
        ("oWAT files", c["oWAT_files"], f"paper {PAPER['oWAT']}"),
        ("rows with a Kleiner grade, initial_sample", c["kleiner_group:initial_sample"], ""),
        ("rows with a Kleiner grade, at_bariatric_surgery_sample", c["kleiner_group:at_bariatric_surgery_sample"], ""),
        ("rows with a Kleiner grade, follow_up_sample", c["kleiner_group:follow_up_sample"], ""),
        ("patients with 2 or more Kleiner grades", n_ge2, "docs/RESULTS.md says 57"),
        ("  of which grades at surgery and follow-up", both, ""),
        ("  of which surgery only", only_v2, ""),
        ("  of which follow-up only", only_fu, ""),
        ("surgery-visit rows, age difference to initial (years): " + json.dumps(dict(sorted(v2_age.items()))), sum(v2_age.values()), "docs/RESULTS.md says 27 at unchanged age"),
        ("follow-up rows, age difference to initial (years): " + json.dumps(dict(sorted(fu_age.items()))), sum(fu_age.values()), ""),
        ("S1 patients (initial visit, plasma, grade, in filtered matrix)", len(s1_in), f"{len(s1)} before the matrix join"),
        ("S1 patients F0-F1", n_lt2, ""),
        ("S1 patients F2-F4", n_ge2_s1, ""),
        ("S1 patients with complete age, sex, BMI, batch", n_complete, ""),
    ]
    with open(f"{out_dir}/pxd051911_counts.tsv", "w", encoding="utf-8") as fh:
        fh.write("metric\tvalue\tnote\n")
        for m, v, n in counts:
            fh.write(f"{m}\t{v}\t{n}\n")
        for k in sorted(grade):
            fh.write(f"S1 patients Kleiner {k}\t{grade[k]}\t\n")

    # ---- detection gate --------------------------------------------------
    prog = {}
    for r in read_tsv(MEMBERSHIP):
        if r["program_id"] in PROGRAMS:
            prog.setdefault(r["program_id"], set()).add(r["gene_symbol"])

    raw_pl = read_tsv(f"{PX}/plasma_protein_quant.txt")
    pl_genes = {}
    for r in raw_pl:
        for g in genes_of(r["Genes"]):
            pl_genes.setdefault(g, set()).add(r["ProteinAccessions"])
    liv_genes = set()
    for r in read_tsv(f"{PX}/liver_protein_quant.txt"):
        liv_genes |= genes_of(r["Genes"])
    with open(OLINK, encoding="utf-8") as fh:
        next(fh)
        olink_assays = {line.split("\t", 1)[0] for line in fh}

    # missingness of each accession across S1 samples in the deposited filtered matrix
    v1_cols = [i for i, h in enumerate(filt_header) if h in {r["unique_identifier"] for r in s1_in}]
    miss = {}
    with open(f"{PX}/plasma_protein_quant_filtered_and_batch_corrected.txt", encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            vals = [num(f[i]) if i < len(f) else None for i in v1_cols]
            miss[f[0]] = sum(v is None for v in vals) / max(len(vals), 1)

    gate_rows = []
    for pid, genes in prog.items():
        name, role = PROGRAMS[pid]
        det_raw = sorted(g for g in genes if g in pl_genes)
        det_filt = sorted(g for g in genes if any(a in miss and miss[a] <= GATE_MAX_MISSING for a in pl_genes.get(g, ())))
        det_liv = sorted(g for g in genes if g in liv_genes)
        det_olk = sorted(g for g in genes if g in olink_assays)
        for src, det in (("S1 plasma, raw matrix (578 protein groups)", det_raw),
                         ("S1 plasma, deposited filtered matrix, at most 50% missing", det_filt),
                         ("S2 liver proteome (tissue)", det_liv),
                         ("S3 Olink assay list (no stage labels)", det_olk)):
            gate_rows.append([pid, name, role, src, len(genes), len(det), ",".join(det),
                              GATE_MIN_PROTEINS, "pass" if len(det) >= GATE_MIN_PROTEINS else "fail"])
    with open(f"{out_dir}/detection_gate.tsv", "w", encoding="utf-8") as fh:
        fh.write("program_id\tprogram_name\trole\tsource\tn_program_genes\tn_detected\tdetected_genes\tmin_required\tgate\n")
        for g in gate_rows:
            fh.write("\t".join(str(x) for x in g) + "\n")

    primary_s1 = next(g for g in gate_rows if g[0] == "hepatocytes::20" and g[3].startswith("S1 plasma, deposited"))
    g1_pass = primary_s1[8] == "pass"
    summary = {
        "G1_ductular_injury_detected_in_S1_filtered": primary_s1[5], "G1_pass": g1_pass,
        "G2_S1_patients": len(s1_in), "G2_F0_F1": n_lt2, "G2_F2_F4": n_ge2_s1, "G2_pass": g2_pass,
        "stage_test_run": bool(g1_pass and g2_pass),
        "decision": "test eligible" if (g1_pass and g2_pass) else "untestable under the prespecified stop rule",
    }
    with open(f"{out_dir}/gate_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
    for g in gate_rows:
        print("\t".join(str(x) for x in g[:6]) + "\t" + g[8])


if __name__ == "__main__":
    main(sys.argv[1])
