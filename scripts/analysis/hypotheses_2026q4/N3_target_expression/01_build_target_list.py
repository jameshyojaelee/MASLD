#!/usr/bin/env python
"""N3 step 1: derive the MASH drug-target list from the saved Open Targets result.

Implements the rule in PRESPEC_N3.md. Standard library only.
Usage: python -I 01_build_target_list.py <RESULT_DIR>
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

RES = Path(sys.argv[1])
SRC = RES / "inputs" / "ot_mash_drugAndClinicalCandidates_MONDO_0007027.json"
KEEP_STAGES = {"PHASE_2": 2.0, "PHASE_2_3": 2.5, "PHASE_3": 3.0, "APPROVAL": 4.0}
NOT_TARGET_ACTIONS = {"SUBSTRATE", "EXOGENOUS GENE"}
MAX_TARGETS_PER_DRUG = 4
# Target list stated in PRESPEC_N3.md before any expression result was computed.
PRESPEC_LIST = sorted(
    "ACACA ACACB ADORA3 ADRA2A ADRA2B ADRA2C AGTR1 AR CASP1 CASP8 CASP9 CCR2 CCR5 CD3D CD3E CD3G CNR1 "
    "DPP4 ESR1 FGF1 FGFR4 GCGR GIPR GLP1R HMGCR INSR KLB LEPR LOXL2 MAP3K5 NPC1L1 NR1H4 NR3C2 OPRD1 OPRK1 "
    "OPRM1 PDE4A PDE4B PDE4C PDE4D PPARA PPARD PPARG PTGER3 SLC10A2 SLC5A1 SLC5A2 THRA THRB VDR XDH".split())

d = json.load(open(SRC))["data"]["disease"]
rows = d["drugAndClinicalCandidates"]["rows"]
assert d["id"] == "MONDO_0007027" and len(rows) == d["drugAndClinicalCandidates"]["count"]

drug_log, by_target = [], defaultdict(lambda: {"ens": None, "drugs": {}, "co": defaultdict(set)})
for r in rows:
    dr = r["drug"]
    moa = (dr.get("mechanismsOfAction") or {}).get("rows") or []
    genes = {}
    for m in moa:
        if (m["actionType"] or "") in NOT_TARGET_ACTIONS:
            continue
        for t in m["targets"]:
            genes[t["approvedSymbol"]] = t["id"]
    stage = r["maxClinicalStage"]
    if stage not in KEEP_STAGES:
        disp = "excluded_stage"
    elif not genes:
        disp = "no_mechanism_target"
    elif len(genes) > MAX_TARGETS_PER_DRUG:
        disp = f"excluded_more_than_{MAX_TARGETS_PER_DRUG}_targets"
    else:
        disp = "kept"
        for g, ens in genes.items():
            t = by_target[g]
            t["ens"] = ens
            prev = t["drugs"].get(dr["name"], 0.0)
            t["drugs"][dr["name"]] = max(prev, KEEP_STAGES[stage])
            for other in genes:
                if other != g:
                    t["co"][g].add(f"{dr['name']}:{other}")
    drug_log.append((dr["name"], dr["id"], stage, dr["maximumClinicalStage"], len(genes), ";".join(sorted(genes)), disp))

with open(RES / "drug_rows.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["drug", "chembl_id", "mash_row_stage", "drug_max_stage", "n_targets", "targets", "disposition"])
    w.writerows(drug_log)

symbols = sorted(by_target)
if symbols != PRESPEC_LIST:
    sys.exit(f"target list differs from PRESPEC: only_derived={sorted(set(symbols)-set(PRESPEC_LIST))} "
             f"only_prespec={sorted(set(PRESPEC_LIST)-set(symbols))}")

stage_name = {2.0: "PHASE_2", 2.5: "PHASE_2_3", 3.0: "PHASE_3", 4.0: "APPROVAL"}
with open(RES / "target_list.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["symbol", "ensembl", "n_drugs", "best_mash_row_stage", "drugs", "co_targets_listed_for_same_drug"])
    for g in symbols:
        t = by_target[g]
        best = max(t["drugs"].values())
        w.writerow([g, t["ens"], len(t["drugs"]), stage_name[best],
                    ";".join(f"{k}({stage_name[v]})" for k, v in sorted(t["drugs"].items())),
                    ";".join(sorted(t["co"][g]))])
print(f"targets={len(symbols)} drug_rows={len(rows)} kept_rows={sum(1 for x in drug_log if x[-1]=='kept')}")
