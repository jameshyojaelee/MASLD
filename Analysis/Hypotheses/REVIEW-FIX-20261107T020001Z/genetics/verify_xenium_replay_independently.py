#!/usr/bin/env python3
"""Bounded independent scalar reconstruction from saved donor aggregates."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parents[4]
helper = root / "Analysis/Hypotheses/X1-20261028T103410Z/scripts/x1_selection_aware.py"
spec = importlib.util.spec_from_file_location("selection_review", helper)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
output = root / "Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/xenium/results_22810356"
a = dict(np.load(output / "donor_aggregates_selection.npz", allow_pickle=False))
engine = module.SelectionReplay(a)
rng = np.random.default_rng(20261107)
n = len(a["donors"])
hep = list(a["types"]).index("Hepatocyte")
reference = [list(a["genes"]).index(g) for g in ["HPX", "KNG1", "FBP1"]]


def scalar(gene_index, mask, indices):
    use_interior = []
    for k in range(len(a["types"])):
        good = a["n_int"][indices, k] >= 30
        use_interior.append(k == hep or (sum(good[mask]) >= 7 and sum(good[~mask]) >= 7))
    levels, totals, all_levels = [], [], []
    for donor in indices:
        used = np.array([a["m_int" if use_interior[k] else "m_all"][donor, k]
                         for k in range(len(use_interior))])
        used = np.nan_to_num(used, nan=0)
        hepatocyte = a["m_int"][donor, hep]
        used[hep] = hepatocyte
        leak = np.array([np.median(used[k, reference] / hepatocyte[reference]) if k != hep else 0
                         for k in range(len(use_interior))])
        levels.append(np.maximum(used - leak[:, None] * hepatocyte, 0)[:, gene_index])
        totals.append(used.sum(1))
        all_used = np.nan_to_num(a["m_all"][donor].copy(), nan=0)
        all_used[hep] = hepatocyte
        all_leak = np.array([np.median(all_used[k, reference] / hepatocyte[reference]) if k != hep else 0
                             for k in range(len(use_interior))])
        all_levels.append(np.maximum(all_used - all_leak[:, None] * hepatocyte, 0)[:, gene_index])
    levels, totals, all_levels = map(np.array, [levels, totals, all_levels])
    means = levels[~mask].mean(0)
    ordered = sorted(range(len(use_interior)), key=lambda k: (-means[k], k))
    carrier, second = ordered[:2]
    share = a["n_all"][indices, carrier] / a["n_all"][indices].sum(1)
    values = {"within": levels[:, carrier], "second_within": levels[:, second],
              "within_allcells": all_levels[:, carrier], "share": share,
              "within_panelshare": np.divide(levels[:, carrier], totals[:, carrier],
                    out=np.zeros(n), where=totals[:, carrier] > 0)}
    result = {"carrier": carrier, "second_type": second}
    for key, values_ in values.items():
        case, control = values_[mask].mean(), values_[~mask].mean()
        result[key] = 0 if case + control == 0 else (case - control) / (case + control)
    return result


comparisons = 0
genes = ["SCGN", "CYP3A4"]
for gene in genes:
    gene_index = list(a["genes"]).index(gene)
    for mode in ["permutation", "bootstrap"]:
        for iteration in range(24):
            if mode == "permutation":
                mask = np.zeros(n, bool)
                mask[rng.choice(n, 9, replace=False)] = True
                indices = np.arange(n)
            else:
                mask = a["is_mash"].copy()
                indices = np.arange(n)
                for group in [False, True]:
                    original = np.where(mask == group)[0]
                    indices[original] = rng.choice(original, len(original), replace=True)
            observed = engine.statistics(gene, mask[None, :], indices[None, :])
            expected = scalar(gene_index, mask, indices)
            for key, value in expected.items():
                np.testing.assert_allclose(observed[key][0], value, atol=1e-12)
            comparisons += 1

for filename, family_size, correction in [("results_5d.tsv", 14, "within_holm14"),
        ("results_5e.tsv", 12, "within_holm12"), ("results_calls_supplementary.tsv", 103, "within_bh")]:
    with (output / filename).open() as handle:
        table = list(csv.DictReader(handle, delimiter="\t"))
    assert len(table) == family_size and len({row["gene"] for row in table}) == family_size
    p = np.array([float(row["within_p_perm"]) for row in table])
    assert np.isfinite(p).all() and (p >= 0).all() and (p <= 1).all()
    if family_size in [12, 14]:
        ordered = sorted(range(family_size), key=lambda j: p[j])
        q = {j: min(1, max((family_size - k) * p[ordered[k]] for k in range(rank + 1)))
             for rank, j in enumerate(ordered)}
        np.testing.assert_allclose([float(row[correction]) for row in table],
                                   [q[j] for j in range(family_size)], atol=1e-12)
    assert all(abs(float(row["within_p_perm"]) * 48620 - round(float(row["within_p_perm"]) * 48620)) < 1e-8
               for row in table)
with (output / "selection_resampling_diagnostics.tsv").open() as handle:
    diagnostics = list(csv.DictReader(handle, delimiter="\t"))
assert len(diagnostics) == 615
assert all(int(row["n_permutations"]) == 48620 and int(row["n_bootstraps"]) == 10000
           and int(row["n_permutation_nonfinite"]) == 0 for row in diagnostics)
result = {"status": "PASS", "independent_scalar_comparisons": comparisons, "genes": genes,
          "replayed_choices": ["type-level interior eligibility", "stored per-donor fallback",
                               "spillover correction", "primary carrier", "second carrier", "panel and share contrasts"],
          "complete_families": [14, 12, 103], "all_permutation_counts": 48620,
          "all_bootstrap_counts": 10000, "finite_permutation_statistics": True,
          "holm_checked_independently": True, "helper_sha256": hashlib.sha256(helper.read_bytes()).hexdigest(),
          "limits": "Conditional observational donor-label test; does not remove confounding or establish selective/simultaneous interval coverage."}
destination = Path(__file__).with_name("xenium_independent_check.json")
with destination.open("x") as handle:
    json.dump(result, handle, indent=2)
    handle.write("\n")
print(json.dumps(result, indent=2))
