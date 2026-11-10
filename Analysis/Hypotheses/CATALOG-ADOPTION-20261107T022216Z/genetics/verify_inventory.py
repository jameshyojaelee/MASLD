#!/usr/bin/env python3
"""Independent identity, family, coverage and posterior-overlap checks."""
import csv
import gzip
import hashlib
import itertools
import json
import math
import sys
from collections import Counter
from pathlib import Path

here = Path(__file__).resolve().parent
root = here.parents[3]
output = Path(sys.argv[1])


def read(path, delimiter="\t"):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


summary = json.loads((output / "summary.json").read_text())
for source, expected in summary["input_sha256"].items():
    assert sha(Path(source)) == expected
comparisons = read(output / "unordered_gene_study_comparisons.tsv")
signals = read(output / "supported_signal_pair_comparisons.tsv")
all_genes = read(output / "all_436_genes_reporting_universe.tsv")
rv9 = root / "Analysis/Hypotheses/RV9-20261102T151012Z/out-20261102T151110Z"
original = read(rv9 / "rows.tsv")
keys = {(row["ensembl"], tuple(sorted([row["discovery"], row["replication"]]))) for row in original}
assert {(row["ensembl"], (row["study1"], row["study2"])) for row in comparisons} == keys
assert len(keys) == len(comparisons) == 383
assert len(all_genes) == len({row["ensembl"] for row in all_genes}) == 436
assert all(row["study1"] < row["study2"] for row in comparisons)
assert Counter(row["trait"] for row in comparisons) == {"ALT": 178, "AST": 192, "NAFLD": 13}
recurring = [row for row in comparisons if row["gene_level_recurrence"] == "True"]
assert len(recurring) == 94 and len({row["ensembl"] for row in recurring}) == 78
assert all(row["trait_group"] == "liver_enzyme" for row in recurring)
for row in signals:
    assert row["same_recorded_eqtl_signal"] == str(row["eqtl_index1"] == row["eqtl_index2"])
    assert row["gwAS_only_credible_sets_available"] == "False"
    assert row["terminal_mechanism_state"].startswith("unresolved_")
    if row["shared_variant_posterior_overlap"]:
        value = float(row["shared_variant_posterior_overlap"])
        assert math.isfinite(value) and 0 <= value <= 1 + 1e-5
    assert float(row["pp4_1"]) > .5 and float(row["pp4_2"]) > .5

# Independently rebuild one nonempty unsigned posterior comparison from original
# exports, including allele-pair identity; no coordinate-only join is accepted.
selected = next(row for row in signals if row["shared_variant_posterior_overlap"])
replay = root / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/genetics"
plan = {(row["gwas_name"], row["ensembl"]): row for row in read(replay / "prepared/replay_plan.tsv")}
vectors = []
for number in [1, 2]:
    study, ensembl = selected[f"study{number}"], selected["ensembl"]
    p = plan[study, ensembl]
    folder = replay / "replay_execution" / f"batch_{p['batch_id']}" / "exports" / study / f"chr{p['chr']}" / ensembl
    pair = [row for row in read(folder / "signal_pairs.tsv") if
            row["idx1"] == selected[f"gwas_index{number}"] and row["idx2"] == selected[f"eqtl_index{number}"]]
    assert len(pair) == 1
    vector = {}
    with gzip.open(folder / "variant_posteriors.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["signal_pair_index"] != pair[0]["signal_pair_index"]:
                continue
            key = (row["chr"], row["hg19_position"], tuple(sorted([row["allele1"].upper(), row["allele2"].upper()])))
            assert key not in vector
            vector[key] = float(row["SNP.PP.H4"])
    assert abs(sum(vector.values()) - 1) < 1e-5
    vectors.append(vector)
common = set(vectors[0]) & set(vectors[1])
independent_overlap = sum(min(vectors[0][key], vectors[1][key]) for key in common)
assert math.isclose(independent_overlap, float(selected["shared_variant_posterior_overlap"]), abs_tol=1e-12)
result = dict(status="PASS", unique_unordered_comparisons=383,
    recurrence_comparisons=94, recurrence_genes=78, direct_disease_recurrence_comparisons=0,
    all_436_membership_preserved=True, trait_families_preserved=True,
    no_final_mechanism_replication_claim=True, variant_allele_pair_overlap_reproduced=True,
    original_input_hashes_unchanged=True, comparison_example=selected,
    limitation="These checks validate saved-output identity and arithmetic, not adequate LD or the historical source-fit hash.")
with (output / "independent_checks.json").open("x") as handle:
    json.dump(result, handle, indent=2)
    handle.write("\n")
print(json.dumps({key: value for key, value in result.items() if key != "comparison_example"}, indent=2))
