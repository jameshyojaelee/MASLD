#!/usr/bin/env python3
"""N4 step one: pair type for every same-person component, from GEO/SRA records only.

Reads the qualified identity roster and the GEO series matrices / SRA run info saved under
<out>/source_metadata. Reads no expression value. Rules are in <out>/PRESPEC.md.
"""
import sys
from pathlib import Path

import pandas as pd

out = Path(sys.argv[1])
repo = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
roster = pd.read_csv(repo / "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/resource_identity_roster.tsv", sep="\t")
src = out / "source_metadata"

SOURCE = {
    "GSE213621": "GEO GSE213621 series matrix and SRA SRP398015 run info (retrieved 2026-10-07); source paper PMID 36535507 as summarised in docs/STATUS.md (368 subjects, one sample per patient, no repeat biopsy described; not re-read in this run)",
    "GSE135251": "GEO GSE135251 series matrix and SRA SRP217231 run info (retrieved 2026-10-07)",
    "GSE130970": "GEO GSE130970 series matrix and SRA SRP197353 run info (retrieved 2026-10-07)",
    "GSE193066": "GEO GSE193066 series matrix (local copy, RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066/metadata/)",
}

rows = []
multi = roster.groupby("screen_component").filter(lambda x: len(x) > 1)
for gse in ("GSE130970", "GSE135251", "GSE213621"):
    meta = pd.read_csv(src / f"{gse}_samples_joined.tsv", sep="\t", dtype=str).set_index("Run")
    for comp, grp in multi[multi["dataset"] == gse].groupby("screen_component"):
        runs = sorted(grp["sample_id"])
        m = meta.loc[runs]
        distinct = (m["geo_accession"].nunique() == len(runs) and m["BioSample"].nunique() == len(runs)
                    and m["Experiment"].nunique() == len(runs))
        base = f"{len(runs)} libraries; distinct GSM, BioSample and SRX: {distinct}; titles: {' | '.join(m['title'])}"
        if gse == "GSE130970":
            ages = [int(x) for x in m["age at biopsy"]]
            if len(set(ages)) > 1:
                ptype = "different_visits"
                ev = f"{base}; GEO 'age at biopsy' {ages} (differs by {max(ages) - min(ages)} y); same recorded sex {sorted(set(m['Sex']))}; GEO design: 78 distinct liver biopsies"
            else:
                ptype, ev = "not_established", f"{base}; identical age at biopsy {ages}"
        elif gse == "GSE135251":
            ptype = "not_established"
            ev = (f"{base}; the source gives each library a different patient number; no age, date or visit field; "
                  f"library names {' | '.join(m['description'])}; GEO design: 216 snap-frozen biopsies")
        else:
            gsm = sorted(int(x[3:]) for x in m["geo_accession"])
            ptype = "not_established"
            ev = (f"{base}; no patient, age, date or visit field in GEO or SRA; GSM numbers "
                  f"{'consecutive' if gsm[-1] - gsm[0] == len(gsm) - 1 else 'not consecutive'}; "
                  f"source describes one sample per patient")
        rows.append({"cohort": gse, "person_id": f"{gse}::{comp}", "definition": "genotype_component",
                     "libraries": ",".join(runs), "n_libraries": len(runs),
                     "stage_codes": ",".join(str(x) for x in grp.set_index("sample_id").loc[runs, "fibrosis_stage"]),
                     "pair_type": ptype, "technical_resequencing_excluded_by_records": distinct,
                     "evidence": ev, "evidence_source": SOURCE[gse]})

# GSE193066: source-declared first/second biopsy pairs
g = roster[roster["dataset"] == "GSE193066"]
gmeta = pd.read_csv(repo / "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066/metadata/metadata.tsv", sep="\t", dtype=str).set_index("Run")
for person, grp in g.groupby("source_person_id"):
    if len(grp) < 2:
        continue
    runs = sorted(grp["sample_id"])
    comps = grp.set_index("sample_id").loc[runs, "a1_component"]
    same = comps.notna().all() and comps.nunique() == 1
    rows.append({"cohort": "GSE193066", "person_id": f"GSE193066::{person}", "definition": "source_declared_person",
                 "libraries": ",".join(runs), "n_libraries": len(runs),
                 "stage_codes": ",".join(str(x) for x in grp.set_index("sample_id").loc[runs, "fibrosis_stage"]),
                 "pair_type": "different_visits", "technical_resequencing_excluded_by_records": True,
                 "evidence": f"GEO 'biopsy' field: {' | '.join(gmeta.loc[runs, 'biopsy'])}; GEO design: 2nd follow-up biopsy for 58 high-risk patients; "
                             f"existing KING table places both libraries in one genetic individual: {bool(same)} "
                             f"(a1_component {' | '.join(comps.fillna('no_call'))})",
                 "evidence_source": SOURCE["GSE193066"], "genotype_same_individual": bool(same)})
for comp, grp in multi[multi["dataset"] == "GSE193066"].groupby("screen_component"):
    if grp["source_person_id"].nunique() > 1:
        runs = sorted(grp["sample_id"])
        rows.append({"cohort": "GSE193066", "person_id": f"GSE193066::genotype::{comp}", "definition": "genotype_component_crossing_source_persons",
                     "libraries": ",".join(runs), "n_libraries": len(runs),
                     "stage_codes": ",".join(str(x) for x in grp.set_index("sample_id").loc[runs, "fibrosis_stage"]),
                     "pair_type": "not_established", "technical_resequencing_excluded_by_records": True,
                     "evidence": f"genotype component joins source persons {' | '.join(grp.set_index('sample_id').loc[runs, 'source_person_id'])}; excluded from analysis",
                     "evidence_source": SOURCE["GSE193066"]})

t = pd.DataFrame(rows)
t.to_csv(out / "pair_type_table.tsv", sep="\t", index=False)
print(t.groupby(["cohort", "definition", "pair_type"]).size())
if "genotype_same_individual" in t:
    print(t[t["cohort"] == "GSE193066"].groupby(["definition", "genotype_same_individual"], dropna=False).size())
