#!/usr/bin/env python3
"""Post hoc X1 correction: replay label-dependent choices in donor resampling.

Input aggregates precede cell-set selection. No raw cells or new gene selection
are needed. Exact label tests assume donor-label exchangeability; they do not
remove the observational age/tissue-source confounding.
"""
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


def bounded_difference(case, control):
    """Monotone in log-ratio, including one-sided zeros; both zero means no signal."""
    total = case + control
    return np.divide(case - control, total, out=np.zeros_like(total), where=total > 0)


def log_ratio(case, control):
    with np.errstate(divide="ignore", invalid="ignore"):
        value = np.log2(case / control)
    return np.where((case == 0) & (control == 0), 0.0, value)


def holm(p):
    # Keep the complete predeclared family, including any unavailable test.
    p = np.asarray(p, float)
    assert np.isfinite(p).all()
    order = np.argsort(p)
    adjusted = np.maximum.accumulate(p[order] * (len(p) - np.arange(len(p)))).clip(0, 1)
    answer = np.empty_like(p)
    answer[order] = adjusted
    return answer


class SelectionReplay:
    def __init__(self, aggregates, minimum_interior=30, minimum_arm=7):
        self.a = aggregates
        self.genes = list(aggregates["genes"])
        self.types = list(aggregates["types"])
        self.hep = self.types.index("Hepatocyte")
        self.n_donors = len(aggregates["donors"])
        self.minimum_interior = minimum_interior
        self.minimum_arm = minimum_arm
        self.share = aggregates["n_all"] / aggregates["n_all"].sum(1)[:, None]
        ref = [self.genes.index(g) for g in ("HPX", "KNG1", "FBP1")]
        hep_ref = aggregates["m_int"][:, self.hep, :]
        self.levels = []
        self.totals = []
        for name in ("m_all", "m_int"):
            used = np.nan_to_num(aggregates[name].copy(), nan=0.0)
            # The all-cell sensitivity still uses interior hepatocytes.
            used[:, self.hep, :] = hep_ref
            with np.errstate(divide="ignore", invalid="ignore"):
                leak = np.nanmedian(used[:, :, ref] / hep_ref[:, None, ref], axis=2)
            leak[:, self.hep] = 0.0
            expected = leak[:, :, None] * hep_ref[:, None, :]
            corrected = np.maximum(used - expected, 0.0)
            if not np.isfinite(corrected).all():
                raise ValueError("Nonfinite corrected levels; no label test is eligible")
            self.levels.append(corrected)
            self.totals.append(used.sum(2))

    def choices(self, masks, indices=None):
        masks = np.asarray(masks, bool)
        if indices is None:
            indices = np.broadcast_to(np.arange(self.n_donors), masks.shape)
        ok = self.a["n_int"][indices] >= self.minimum_interior
        use = ((ok & masks[:, :, None]).sum(1) >= self.minimum_arm) & (
            (ok & ~masks[:, :, None]).sum(1) >= self.minimum_arm)
        use[:, self.hep] = True
        return indices, use

    def statistics(self, gene, masks, indices=None, use=None):
        masks = np.asarray(masks, bool)
        if use is None:
            indices, use = self.choices(masks, indices)
        j = self.genes.index(gene)
        levels = np.where(use[:, None, :], self.levels[1][indices, :, j],
                          self.levels[0][indices, :, j])
        total = np.where(use[:, None, :], self.totals[1][indices], self.totals[0][indices])
        control_means = np.einsum("sd,sdk->sk", ~masks, levels) / (~masks).sum(1)[:, None]
        # Stable order gives the same first-maximum rule as numpy.argmax.
        order = np.argsort(-control_means, axis=1, kind="stable")
        carrier = order[:, 0]
        second = order[:, 1]
        selector = carrier[:, None, None]
        selected = np.take_along_axis(levels, selector, axis=2)[:, :, 0]
        selected_total = np.take_along_axis(total, selector, axis=2)[:, :, 0]
        all_cells = np.take_along_axis(self.levels[0][indices, :, j], selector, axis=2)[:, :, 0]
        second_cells = np.take_along_axis(levels, second[:, None, None], axis=2)[:, :, 0]
        share = np.take_along_axis(self.share[indices], selector, axis=2)[:, :, 0]
        panel = np.divide(selected, selected_total, out=np.zeros_like(selected), where=selected_total > 0)
        values = {"within": selected, "within_panelshare": panel,
                  "within_allcells": all_cells, "second_within": second_cells, "share": share}
        result = {"carrier": carrier, "second_type": second}
        for name, value in values.items():
            case = (value * masks).sum(1) / masks.sum(1)
            control = (value * ~masks).sum(1) / (~masks).sum(1)
            result[name] = bounded_difference(case, control)
            result[name + "_log2"] = log_ratio(case, control)
            result[name + "_both_zero"] = (case == 0) & (control == 0)
        return result


def correct_results(output, n_boot=10000, seed=20261027):
    output = Path(output)
    a = dict(np.load(output / "donor_aggregates_selection.npz", allow_pickle=False))
    engine = SelectionReplay(a)
    observed_masks = a["is_mash"][None, :]
    masks = np.array([np.isin(np.arange(engine.n_donors), c)
                      for c in itertools.combinations(range(engine.n_donors), int(observed_masks.sum()))])
    p_indices, p_use = engine.choices(masks)
    rng = np.random.default_rng(seed)
    case = np.where(a["is_mash"])[0]
    control = np.where(~a["is_mash"])[0]
    b_indices = np.empty((n_boot, engine.n_donors), int)
    b_indices[:, case] = case[rng.integers(0, len(case), (n_boot, len(case))) ]
    b_indices[:, control] = control[rng.integers(0, len(control), (n_boot, len(control))) ]
    b_masks = np.broadcast_to(a["is_mash"], b_indices.shape)
    b_indices, b_use = engine.choices(b_masks, b_indices)
    files = ["results_5d.tsv", "results_5e.tsv", "results_calls_supplementary.tsv"]
    tables = {f: pd.read_csv(output / f, sep="\t") for f in files}
    fixed_genes = sorted(set().union(*(set(t.gene) for t in tables.values())))
    results, checks, diagnostics = {}, [], []
    for i, gene in enumerate(fixed_genes):
        obs = engine.statistics(gene, observed_masks)
        perm = engine.statistics(gene, masks, p_indices, p_use)
        boot = engine.statistics(gene, b_masks, b_indices, b_use)
        result = {"gene": gene}
        for name in ("within", "within_panelshare", "within_allcells", "second_within", "share"):
            # Bounded contrast has the same ordering as |log2 ratio| but is
            # defined at zeros; no permutation or bootstrap replicate is dropped.
            p = np.mean(np.abs(perm[name]) >= abs(obs[name][0]) - 1e-12)
            lo, hi = np.quantile(boot[name + "_log2"], [.025, .975], method="inverted_cdf")
            result.update({name + "_log2_ratio": float(obs[name + "_log2"][0]),
                           name + "_p_perm": float(p), name + "_ci_lo": float(lo), name + "_ci_hi": float(hi)})
            diagnostics.append(dict(gene=gene, statistic=name,
                n_permutations=len(masks), n_bootstraps=n_boot,
                n_permutation_both_zero=int(perm[name + "_both_zero"].sum()),
                n_bootstrap_both_zero=int(boot[name + "_both_zero"].sum()),
                n_bootstrap_infinite=int(np.isinf(boot[name + "_log2"]).sum()),
                n_permutation_nonfinite=int((~np.isfinite(perm[name])).sum()),
                bootstrap_same_carrier_fraction=float(np.mean(boot["carrier"] == obs["carrier"][0]))))
        for f, t in tables.items():
            if gene in set(t.gene):
                old = t.set_index("gene").loc[gene]
                delta = abs(old.within_log2_ratio - result["within_log2_ratio"])
                checks.append(dict(gene=gene, family=f, observed_effect_difference=delta,
                    carrier_identical=old.carrier == engine.types[obs["carrier"][0]]))
                assert np.isclose(old.within_log2_ratio, result["within_log2_ratio"], equal_nan=True)
                assert old.carrier == engine.types[obs["carrier"][0]]
        results[gene] = result
        print(f"selection replay {i + 1}/{len(fixed_genes)}: {gene}", flush=True)
    replay = pd.DataFrame(results.values()).set_index("gene")
    for f, table in tables.items():
        for col in replay.columns:
            if col in table:
                table["selection_conditioned_old_" + col] = table[col]
            table[col] = replay.loc[table.gene, col].to_numpy()
        table["selection_conditioned_old_second_within_p"] = table.second_within_p
        table["second_within_p"] = table.second_within_p_perm
        correction = "within_holm14" if f == files[0] else "within_holm12" if f == files[1] else "within_bh"
        table["selection_conditioned_old_" + correction] = table[correction]
        if correction == "within_bh":
            p = table.within_p_perm.to_numpy()
            order = np.argsort(p)
            q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1].clip(0, 1)
            answer = np.empty_like(q); answer[order] = q
            table[correction] = answer
        else:
            table[correction] = holm(table.within_p_perm)
        table["correction_family_n"] = len(table)
        table["inference"] = "post_hoc_selection_aware_exact_donor_label_permutation; bootstrap_repeats_selection"
        table.to_csv(output / f, sep="\t", index=False)
    pd.DataFrame(checks).to_csv(output / "selection_observed_replay_checks.tsv", sep="\t", index=False)
    pd.DataFrame(diagnostics).to_csv(output / "selection_resampling_diagnostics.tsv", sep="\t", index=False)
    pd.DataFrame(p_use.mean(0), index=engine.types, columns=["permutation_interior_fraction"]).to_csv(
        output / "selection_cellset_permutation_frequencies.tsv", sep="\t")
    summary = {"status": "completed", "post_hoc_correction": True, "seed": seed,
        "n_permutations": len(masks), "n_bootstraps": n_boot, "n_unique_genes": len(fixed_genes),
        "families": {f: len(t) for f, t in tables.items()},
        "zero_rule": "bounded_difference=(case-control)/(case+control); both_zero=0; no_resamples_dropped",
        "interval": "pointwise_percentile_donor_bootstrap_repeating_selection; not_simultaneous_or_selective_CI",
        "assumption": "donor_labels_exchangeable; age_and_tissue_source_confounding_not_removed"}
    (output / "selection_completion.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    correct_results(args.output)
