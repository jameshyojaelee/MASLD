#!/usr/bin/env python3
"""Scientific evidence-state invariants for the new missing SuSiE result."""
import csv
import importlib.util
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[4]
path = root / "GWAS/finemapping/src/perf/rebuild_tier12.py"
spec = importlib.util.spec_from_file_location("aggregation", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
rows = []
for gene, methods in [("ENSG00000000001", [("susie_untestable_ambiguous_variant_identity", "")]),
                      ("ENSG00000000002", [("susie_untestable_ambiguous_variant_identity", ""), ("susie", ".3")]),
                      ("ENSG00000000003", [("susie_untestable_ambiguous_variant_identity", ""), ("susie", ".9")])]:
    for number, (method, value) in enumerate(methods):
        rows.append(dict(ensembl=gene, gwas_name=f"s{number}", method=method,
                         **{"PP.H4.susie": value, "PP.H4.abf": ".8"}))
with tempfile.TemporaryDirectory() as directory:
    source = Path(directory) / "fixture.csv"
    with source.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    best, unmapped = module.build(source, {"s0": 2, "s1": 2}, {"s0": "ALT", "s1": "ALT"}, {"s0": "main", "s1": "main"})
assert not unmapped
assert module.susie_state(best["ENSG00000000001"]) == "untestable"
assert module.susie_state(best["ENSG00000000002"]) == "untestable"
assert module.susie_state(best["ENSG00000000003"]) == "positive"
assert all(record["u"] == 1 and record["a"] == .8 for record in best.values())
print("PASS: ambiguous identity survives aggregation; ABF and independent positive support preserved")
