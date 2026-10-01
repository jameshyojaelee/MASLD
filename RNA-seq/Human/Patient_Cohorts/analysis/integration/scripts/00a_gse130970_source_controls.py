#!/usr/bin/env python3
"""Per-biopsy source control status for GSE130970 (Hoang et al. 2019, Sci Rep 9:12541).

GEO and SRA carry histology grades but no control flag. The paper defines
controls as "normal liver histology" and its Supplementary Table 1 summarises
them: N=6, steatosis grade 0 in all six, 2 males, age 49.8 +/- 21 years
(MOESM1_ESM.pdf). NAFLD required >5% steatosis, but two NAFLD biopsies are
also listed at steatosis grade 0, so steatosis alone does not identify the
controls.

This script finds every 6-biopsy subset of the steatosis-0 biopsies whose sex
count, age mean and age SD match Table 1 at its printed precision, and fails
unless exactly one subset matches. That subset is the primary control set.
The strict sensitivity set keeps only controls with steatosis 0 and NAS 0 in
GEO; two reconstructed controls (lobular inflammation 1 in GEO) fail that rule.

Replaces the fibrosis_stage == 0 rule in 00_harmonize_metadata.R, which made
16 steatotic F0 biopsies controls.
"""
import csv
import itertools
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[6]
SRC = ROOT / "RNA-seq/Human/Patient_Cohorts/metadata/source/GSE130970"
RUN_TABLE = SRC / "GSE130970_SraRunTable.csv"
OUT = SRC / "GSE130970_source_diagnosis.tsv"

# Supplementary Table 1, histologically normal controls.
TABLE1_N = 6
TABLE1_MALES = 2
TABLE1_AGE_MEAN = 49.8   # printed to 0.1
TABLE1_AGE_SD = 21       # printed to integer
N_BIOPSIES = 78


def matches_table1(rows):
    ages = [int(r["age_at_biopsy"]) for r in rows]
    males = sum(r["sex"] == "male" for r in rows)
    return (males == TABLE1_MALES
            and abs(statistics.mean(ages) - TABLE1_AGE_MEAN) <= 0.05
            and abs(statistics.stdev(ages) - TABLE1_AGE_SD) <= 0.5)


def main():
    rows = list(csv.DictReader(open(RUN_TABLE)))
    if len(rows) != N_BIOPSIES:
        sys.exit(f"expected {N_BIOPSIES} biopsies, found {len(rows)}")
    steat0 = [r for r in rows if r["steatosis_grade"] == "0"]
    hits = [c for c in itertools.combinations(steat0, TABLE1_N) if matches_table1(c)]
    if len(hits) != 1:
        sys.exit(f"Table 1 reconstruction is not unique: {len(hits)} matching subsets")
    controls = {r["Run"] for r in hits[0]}

    fields = ["run", "gsm", "source_control_status", "strict_control_nas0",
              "steatosis_grade", "lobular_inflammation_grade",
              "cytological_ballooning_grade", "nafld_activity_score",
              "fibrosis_stage", "sex_geo", "age_at_biopsy", "control_rule"]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            is_ctrl = r["Run"] in controls
            strict = is_ctrl and r["nafld_activity_score"] == "0"
            if is_ctrl and r["steatosis_grade"] != "0":
                sys.exit(f"steatotic biopsy assigned control: {r['Run']}")
            w.writerow({
                "run": r["Run"], "gsm": r["Sample Name"],
                "source_control_status": "Control" if is_ctrl else "NAFLD",
                "strict_control_nas0": "Control" if strict else "NAFLD",
                "steatosis_grade": r["steatosis_grade"],
                "lobular_inflammation_grade": r["lobular_inflammation_grade"],
                "cytological_ballooning_grade": r["cytological_ballooning_grade"],
                "nafld_activity_score": r["nafld_activity_score"],
                "fibrosis_stage": r["fibrosis_stage"],
                "sex_geo": r["sex"], "age_at_biopsy": r["age_at_biopsy"],
                "control_rule": "hoang2019_table1_unique_reconstruction",
            })
    print(f"wrote {OUT}: {len(controls)} controls "
          f"({', '.join(sorted(controls))}); strict NAS0 subset "
          f"{sum(r['nafld_activity_score'] == '0' for r in hits[0])}")


if __name__ == "__main__":
    main()
