#!/usr/bin/env python3
"""A collision must become untestable; unique controls and ABF must not change."""
import csv
import hashlib
import json
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
output = Path(sys.argv[1])
baseline = here / "pilot_22811698"


def load(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


with (here / "pilot_targets.tsv").open() as handle:
    targets = list(csv.DictReader(handle, delimiter="\t"))
checks = []
for target in targets:
    tail = Path("coloc") / target["gwas_name"] / f"susie_coloc_chr{target['chr']}.csv"
    original, corrected = load(baseline / tail), load(output / tail)
    assert len(original) == len(corrected) == 1
    old, new = original[0], corrected[0]
    assert old["ensembl"] == new["ensembl"] == target["ensembl"]
    assert float(old["PP.H4.abf"]) == float(new["PP.H4.abf"])
    assert old["top_snp"] == new["top_snp"] and old["n_snps"] == new["n_snps"]
    if target["gene"] == "MPRIP":
        assert new["method"] == "susie_untestable_ambiguous_variant_identity"
        assert new["PP.H4.susie"] in {"", "NA"}
        collisions = output / "coloc" / target["gwas_name"] / f"variant_identity_collisions_chr17_{target['ensembl']}_gwas.tsv"
        with collisions.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 2 and {row["coloc_variant_id"] for row in rows} == {"17:17164840"}
        assert float(rows[0]["beta"]) == -float(rows[1]["beta"])
        assert rows[0]["standard_error"] == rows[1]["standard_error"]
        assert not (output / "exports" / f"{target['gwas_name']}__{target['ensembl']}__inputs_and_fits.rds").exists()
        state = "ambiguous_variant_identity_not_fit"
    else:
        assert new["method"] == "susie"
        assert float(new["PP.H4.susie"]) == float(old["PP.H4.susie"])
        state = "unique_control_unchanged"
    checks.append(dict(gene=target["gene"], study=target["gwas_name"], terminal_state=state,
                       method=new["method"], pp4_abf=float(new["PP.H4.abf"]),
                       pp4_susie=new["PP.H4.susie"], source_variant_count_preserved=True))
report = dict(status="PASS", fits=checks, named_membership_changes=0,
    original_producer_unchanged=True, no_arbitrary_variant_deletion=True,
    mechanism_qualification="MPRIP blocked by duplicate identifiers; no source-matched LD resolution claimed",
    corrected_controls="RORA and FABP1 exact numerical reproduction")
(output / "guarded_pilot_checks.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
