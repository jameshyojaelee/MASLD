#!/usr/bin/env python3
"""Gate completed extraction on per-donor barcode coverage and depth."""

import csv
import json
import os
from pathlib import Path

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
OUT = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2"
rows = list(csv.DictReader((OUT / "fragment_qc_by_donor.tsv").open(), delimiter="\t"))
audit = []
for row in rows:
    coverage = int(row["barcodes_with_autosomal_fragments"]) / int(row["allowlist_barcodes"])
    depth = int(row["retained_fragment_rows"])
    audit.append({"donor_id": row["donor_id"], "barcode_fragment_coverage": coverage,
                  "retained_fragment_rows": depth,
                  "pass": coverage >= 0.75 and depth >= 1_000_000})
verdict = {
    "n_donors": len(audit),
    "barcode_fragment_coverage_min": min(x["barcode_fragment_coverage"] for x in audit),
    "retained_fragment_rows_min": min(x["retained_fragment_rows"] for x in audit),
    "thresholds": {"barcode_fragment_coverage_min": 0.75, "retained_fragment_rows_min": 1_000_000},
    "failed_donors": [x["donor_id"] for x in audit if not x["pass"]],
    "pass": len(audit) == 18 and all(x["pass"] for x in audit),
    "donors": audit,
}
with (OUT / "donor_extraction_gate.json").open("w") as handle:
    json.dump(verdict, handle, indent=2, sort_keys=True)
    handle.write("\n")
print(json.dumps(verdict, indent=2, sort_keys=True))
if not verdict["pass"]:
    raise SystemExit("donor extraction QC failed")
