#!/usr/bin/env python3
"""Check biological identity, posterior units, and preserved evidence states."""
import csv
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path


def load(path, delimiter="\t"):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


out = Path(sys.argv[1])
report = json.loads((out / "checks.json").read_text())
for entry in report["inputs"].values():
    assert sha(Path(entry["path"])) == entry["sha256"], entry["path"]
for name, checksum in report["outputs"].items():
    assert sha(out / name) == checksum, name
native = load(Path(report["inputs"]["native_tier12"]["path"]), ",")
rows = load(out / "genetics_qualification.tsv")
assert len({row["ensembl"] for row in rows}) == len(native) == len(rows) == 18976
assert {row["ensembl"] for row in native} == {row["ensembl"] for row in rows}
native_positive = {row["ensembl"] for row in native
                   if row["coloc_best_susie_pp4"] and float(row["coloc_best_susie_pp4"]) > 0.5}
for column in ["native_susie_state", "priority_state_after_saved_eligibility_check",
               "any_eligible_state_after_saved_eligibility_check"]:
    assert {row["ensembl"] for row in rows if row[column] == "positive"} == native_positive
for row in rows:
    values = [row[column] for column in ["native_susie_pp4", "susie_pp4_p12_1e_6", "susie_pp4_p12_1e_5"]]
    if not values[0]:
        assert values == ["", "", ""], "missing SuSiE result gained a posterior"
        continue
    primary, lower, upper = map(float, values)
    assert 0 <= lower <= primary <= upper <= 1
    # Independent odds-ratio identity, rather than the builder's reweighting formula.
    # Odds near machine-precision endpoints amplify harmless rounding.
    if 0 < primary < 1 - 1e-8:
        primary_odds = primary / (1 - primary)
        for posterior, multiplier in [(lower, 0.2), (upper, 2)]:
            if posterior < 1 - 1e-8:
                assert math.isclose(posterior / (1 - posterior), primary_odds * multiplier, rel_tol=1e-5)
    assert row["named_class_changed"] == "False"

proposals = load(out / "catalog_eligibility_proposals.tsv")
for row in proposals:
    assert row["named_catalog_class"] not in {"genetically_anchored", "concordant"}
    for rule in ["priority", "any_eligible"]:
        proposed = row[f"proposed_class_{rule}_rule"]
        transition = row[f"{rule}_rule_eligibility_transition"] == "True"
        assert proposed == ("untested" if transition else row["named_catalog_class"])
assert sum(row["priority_rule_class_would_change"] == "True" for row in proposals) == 180
assert sum(row["any_eligible_rule_class_would_change"] == "True" for row in proposals) == 39

ld = load(out / "ld_sensitive_genes.tsv")
assert len(ld) == 14
states = Counter(row["posthoc_sensitivity_state"] for row in ld)
assert states == {"support_below_threshold": 10, "no_credible_set_result": 1, "retains_support": 3}
assert {row["symbol"] for row in ld if row["retains_any_tested_fit_above_0_5"] == "True"} == {"GALK1", "PKN3", "TMPRSS2"}
assert sum(int(row["n_anchor_abs_z_residual_gt_5_45"]) for row in ld) == 14
result = {"status": "PASS", "native_gene_ids": len(rows), "positive_membership_preserved": len(native_positive),
          "named_catalog_changes": 0, "proposed_priority_class_changes": 180,
          "proposed_any_eligible_class_changes": 39, "ld_sensitivity_states": dict(states),
          "source_and_output_hashes_match": True, "posterior_bounds_monotonicity_odds_and_missingness": True,
          "checker_sha256": sha(Path(__file__)),
          "scope": "Saved-output arithmetic and joins only; no new LD fit, replication study, or claim of biological validity."}
destination = out / "independent_checks.json"
with destination.open("x") as handle:
    json.dump(result, handle, indent=2)
    handle.write("\n")
print(json.dumps(result, indent=2))
