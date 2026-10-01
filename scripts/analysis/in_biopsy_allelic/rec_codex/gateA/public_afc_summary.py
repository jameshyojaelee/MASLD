"""Outcome-blind GTEx liver aFC summary for the v2 public tag eGenes."""

from __future__ import annotations

import csv
import gzip
import json
import math
import statistics
from pathlib import Path


ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
TAGS = ROOT / "Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-tags-v2-20260929T134419Z/tags"
GTEX = ROOT / "data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz"


def summarize():
    leads = {}
    for path in sorted(TAGS.glob("chr*.tags.tsv")):
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                gene = row["gene_id"]
                lead = row["lead_variant_id"]
                if gene in leads and leads[gene] != lead:
                    raise ValueError(f"more than one lead for {gene}")
                leads[gene] = lead
    matched = {}
    with gzip.open(GTEX, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["gene_id"]
            if gene not in leads:
                continue
            if gene in matched:
                raise ValueError(f"duplicate GTEx row for {gene}")
            if row["variant_id"] != leads[gene]:
                raise ValueError(f"lead allele mismatch for {gene}")
            value = float(row["log2_aFC"])
            matched[gene] = abs(value * math.log(2)) if math.isfinite(value) else float("nan")
    if set(leads) != set(matched):
        raise ValueError(f"GTEx/tag gene mismatch: {len(set(leads) ^ set(matched))} IDs")
    if any(not math.isfinite(value) for value in matched.values()):
        raise ValueError("nonfinite GTEx aFC")
    values = list(matched.values())
    return {
        "source": str(GTEX),
        "tag_directory": str(TAGS),
        "gene_count": len(values),
        "effect_unit": "absolute natural-log allelic fold change, GTEx lead ALT vs REF",
        "thresholds": [
            {
                "theta": theta,
                "eligible_genes": sum(value >= theta for value in values),
                "median_abs_ln_afc": statistics.median(value for value in values if value >= theta),
            }
            for theta in (0.0, 0.1, 0.2, 0.3)
        ],
    }


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=2))
