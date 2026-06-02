#!/usr/bin/env python3
"""B5: Extract library_prep_batch from SRA metadata per cohort.

Strategy per cohort (after auditing SraRunTables):
  - PRJNA512027 (Gerhard): `Library Name` prefix (L0/S0) - DOCUMENTED L0/S0 batch confound
  - GSE193066 (Hoshida): `LibraryName` is per-sample (unique), use LoadDate as batch proxy
  - All other cohorts: no usable batch column in SRA → set NA
    (B2 will fall back to dataset-only random effect)

Output:
  - RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata_with_batch.csv
  - agent_status/B5_metadata_with_batch.csv (copy for B2)
"""
import csv
import os
import re
import shutil
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation")
MAIN = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
UNIFIED = MAIN / "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"
COHORT_DIR = MAIN / "RNA-seq/Human/Patient_Cohorts/pipelines/custom"
OUT = ROOT / "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata_with_batch.csv"
COPY = ROOT / "agent_status/B5_metadata_with_batch.csv"

# Per-cohort batch extractors. Each returns dict[Run] -> batch_label or "" if not assignable.
def load_sra(ds):
    p = COHORT_DIR / ds / "metadata/SraRunTable.csv"
    if not p.exists():
        # legacy path
        p2 = MAIN / f"RNA-seq/Human/Patient_Cohorts/archive/old_data_metadata/metadata/{ds}_SraRunTable.csv"
        if p2.exists():
            p = p2
        else:
            return None
    with open(p) as f:
        return list(csv.DictReader(f))

def batch_prjna512027(rows):
    # Library Name like L03346 / S00153 → prefix L0 or S0 (L=liver, S=stomach? no, here it's library-prep batch)
    out = {}
    for r in rows:
        ln = r.get("Library Name", "") or r.get("LibraryName", "")
        if ln:
            m = re.match(r"^([LS])\d", ln)
            out[r["Run"]] = f"PRJNA512027_{m.group(1)}0" if m else f"PRJNA512027_other"
        else:
            out[r["Run"]] = ""
    return out

def batch_gse193066(rows):
    # LibraryName is per-sample (164 unique); use LoadDate floored to day as batch proxy
    out = {}
    for r in rows:
        ld = r.get("LoadDate", "")
        # take date portion only
        d = ld.split(" ")[0] if ld else ""
        out[r["Run"]] = f"GSE193066_{d}" if d else ""
    return out

EXTRACTORS = {
    "PRJNA512027": batch_prjna512027,
    "GSE193066": batch_gse193066,
}

# Build run->batch map
run2batch = {}
for ds, fn in EXTRACTORS.items():
    rows = load_sra(ds)
    if rows is None:
        print(f"WARN: no SRA table for {ds}")
        continue
    m = fn(rows)
    run2batch.update(m)
    uniq = sorted(set(v for v in m.values() if v))
    print(f"{ds}: assigned batch labels for {sum(1 for v in m.values() if v)}/{len(m)} samples; unique batches = {len(uniq)} {uniq[:5]}")

# Also audit all other cohorts and log "no batch proxy available"
audited = set(EXTRACTORS.keys())
for ds_dir in sorted(COHORT_DIR.glob("GSE*")):
    ds = ds_dir.name
    if ds in audited:
        continue
    sra = ds_dir / "metadata/SraRunTable.csv"
    if not sra.exists():
        # try legacy
        sra = MAIN / f"RNA-seq/Human/Patient_Cohorts/archive/old_data_metadata/metadata/{ds}_SraRunTable.csv"
    if sra.exists():
        print(f"{ds}: no batch proxy in SRA metadata → library_prep_batch=NA (B2 falls back to (1|dataset))")

# Merge into unified_metadata
with open(UNIFIED) as f:
    reader = csv.DictReader(f)
    cols = list(reader.fieldnames) + ["library_prep_batch"]
    rows_out = []
    n_assigned = 0
    for r in reader:
        b = run2batch.get(r["sample_id"], "")
        r["library_prep_batch"] = b if b else "NA"
        if b: n_assigned += 1
        rows_out.append(r)

OUT.parent.mkdir(parents=True, exist_ok=True)
COPY.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows_out)
shutil.copy(OUT, COPY)

print(f"\nWROTE {OUT}")
print(f"WROTE {COPY}")
print(f"Total samples: {len(rows_out)}; with library_prep_batch != NA: {n_assigned}")

# Per-dataset summary
from collections import Counter, defaultdict
per_ds = defaultdict(Counter)
for r in rows_out:
    per_ds[r["dataset"]][r["library_prep_batch"]] += 1
print("\nPer-dataset batch summary:")
for ds in sorted(per_ds):
    cnts = per_ds[ds]
    if list(cnts.keys()) == ["NA"]:
        print(f"  {ds}: NA (no batch proxy)")
    else:
        print(f"  {ds}: {dict(cnts)}")
