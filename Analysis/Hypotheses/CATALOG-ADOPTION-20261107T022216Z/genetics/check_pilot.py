#!/usr/bin/env python3
"""Check exact target identity, terminal status and r2 H4 reproduction."""
import csv
import gzip
import json
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
root = here.parents[3]
output = Path(sys.argv[1])
with (here / "pilot_targets.tsv").open() as handle:
    targets = list(csv.DictReader(handle, delimiter="\t"))
expected = {}
with gzip.open(root / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz", "rt") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["eligible"] == "TRUE" and row["PP.H4"] not in {"", "NA"}:
            key = row["gwas"], row["ensembl"]
            expected[key] = max(expected.get(key, 0), float(row["PP.H4"]))
checks = []
for target in targets:
    path = output / "coloc" / target["gwas_name"] / f"susie_coloc_chr{target['chr']}.csv"
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1 and rows[0]["ensembl"] == target["ensembl"], path
    row = rows[0]
    assert row["method"] == "susie", (target, row["method"])
    value = float(row["PP.H4.susie"])
    reference = expected[target["gwas_name"], target["ensembl"]]
    assert abs(value - reference) <= 1e-4, (target, value, reference)
    prefix = output / "exports" / f"{target['gwas_name']}__{target['ensembl']}"
    for suffix in ["__inputs_and_fits.rds", "__aligned_variants.tsv", "__gwas_credible_sets.tsv"]:
        assert Path(str(prefix) + suffix).is_file(), suffix
    checks.append(dict(**target, pp4_reconstructed=value, pp4_r2=reference,
                       absolute_difference=abs(value-reference), state="original_reference_reproduced"))
(output / "pilot_checks.json").write_text(json.dumps(dict(status="PASS", fits=checks,
    genetic_membership_changes=0, source_matched_MVP_LD_available=False,
    matched_LD_question="unresolved; reproducing the original fit is not a correction",
    scope="three original-reference reconstruction fits with isolated fine-mapping exports"), indent=2) + "\n")
print(json.dumps(checks, indent=2))
