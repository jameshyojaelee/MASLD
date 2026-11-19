#!/usr/bin/env python3
"""X1 (20261118 N09 remediation) selection-aware exact inference and bootstrap stability.

Input aggregates precede cell-set selection. Remediates N09 circularity and CIs:
  - Primary hepatocyte bleed reference REF = ("KNG1", "FBP1") (HPX removed); emits
    ref_sensitivity_hpx.tsv comparing ("KNG1", "FBP1") vs ("HPX", "KNG1", "FBP1") on the
    exact same 48,620 permutations and 10,000 seeded donor bootstrap draws.
  - Selects carrier via permutation-invariant 18-donor grand mean with made_all spillover
    and donor sufficiency enforced BEFORE selection (guarantees E_{H0}[T] = 0 exactly).
  - Primary 95% CIs are fixed-carrier conditional donor-bootstrap CIs (avoiding cross-lineage
    bimodal mixing) reported alongside carrier_bootstrap_stability, conditional-on-same-carrier
    CIs, and unconditional selection-replaying bootstrap CIs.
  - Primary Holm/BH multiplicity families exclude underpowered (<7 eligible donors, e.g. SPI1)
    and LABEL_DEFINING markers ({"KRT8", "KRT18", "MZB1", "LYVE1"}) while retaining every row
    and reporting *_all_sens alongside *_primary.

Exact label tests assume donor-label exchangeability; they do not remove observational
age/tissue-source confounding between control biopsies/resections and MASH-cirrhosis explants.
"""
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_REF = ("KNG1", "FBP1")
SENS_REF_HPX = ("HPX", "KNG1", "FBP1")
STATS = ("within", "within_panelshare", "within_allcells", "second_within", "share")


def bounded_difference(case, control):
    """Monotone in log-ratio, including one-sided zeros; both zero means no signal."""
    total = case + control
    return np.divide(case - control, total, out=np.zeros_like(total), where=total > 0)


def log_ratio(case, control):
    with np.errstate(divide="ignore", invalid="ignore"):
        value = np.log2(case / control)
    return np.where((case == 0) & (control == 0), 0.0, value)


def holm(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    if len(q) == 0:
        return out
    order = np.argsort(q)
    adjusted = np.maximum.accumulate(q[order] * (len(q) - np.arange(len(q)))).clip(0.0, 1.0)
    answer = np.empty_like(q)
    answer[order] = adjusted
    out[ok] = answer
    return out


def holm_subset(p, mask):
    p = np.asarray(p, float)
    mask = np.asarray(mask, bool)
    out = np.full_like(p, np.nan)
    if mask.any():
        assert np.isfinite(p[mask]).all()
        out[mask] = holm(p[mask])
    return out


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    n = len(q)
    if n == 0:
        return out
    order = np.argsort(q)
    adjusted = np.minimum.accumulate((q[order] * n / np.arange(1, n + 1))[::-1])[::-1].clip(0.0, 1.0)
    answer = np.empty_like(q)
    answer[order] = adjusted
    out[ok] = answer
    return out


def bh_subset(p, mask):
    p = np.asarray(p, float)
    mask = np.asarray(mask, bool)
    out = np.full_like(p, np.nan)
    if mask.any():
        assert np.isfinite(p[mask]).all()
        out[mask] = bh(p[mask])
    return out


class SelectionReplay:
    def __init__(self, aggregates, ref_genes=None, minimum_interior=30, minimum_arm=7):
        self.a = aggregates
        self.genes = [str(g) for g in aggregates["genes"]]
        self.types = [str(t) for t in aggregates["types"]]
        self.hep = self.types.index("Hepatocyte")
        self.n_donors = len(aggregates["donors"])
        self.minimum_interior = minimum_interior
        self.minimum_arm = minimum_arm
        self.minimum_present_donors = 2 * minimum_arm
        self.share = aggregates["n_all"] / aggregates["n_all"].sum(1)[:, None]
        self.n_all = aggregates["n_all"].astype(np.float64)
        self.n_used_int = np.where(aggregates["n_int"] >= minimum_interior, aggregates["n_int"], aggregates["n_all"]).astype(
            np.float64
        )
        self.n_used_int[:, self.hep] = aggregates["n_int"][:, self.hep].astype(np.float64)

        if ref_genes is None:
            if "ref_genes" in aggregates:
                ref_genes = tuple(str(g) for g in aggregates["ref_genes"])
            else:
                ref_genes = DEFAULT_REF
        else:
            ref_genes = tuple(str(g) for g in ref_genes)
        self.ref_genes = ref_genes
        ref = [self.genes.index(g) for g in self.ref_genes]

        hep_ref = aggregates["m_int"][:, self.hep, :]
        self.used_levels = []
        self.expected_levels = []
        self.leaks = []
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
            expected[:, self.hep, :] = 0.0
            corrected = np.maximum(used - expected, 0.0)
            if not np.isfinite(corrected).all():
                raise ValueError("Nonfinite corrected levels; no label test is eligible")
            self.used_levels.append(used)
            self.expected_levels.append(expected)
            self.leaks.append(leak)
            self.levels.append(corrected)
            self.totals.append(used.sum(2))

    def choices(self, masks, indices=None):
        masks = np.asarray(masks, bool)
        if indices is None:
            indices = np.broadcast_to(np.arange(self.n_donors), masks.shape)
        ok = self.a["n_int"][indices] >= self.minimum_interior
        use = ((ok & masks[:, :, None]).sum(1) >= self.minimum_arm) & (
            (ok & ~masks[:, :, None]).sum(1) >= self.minimum_arm
        )
        use[:, self.hep] = True
        return indices, use

    def select_carriers(self, gene, indices, use):
        """Permutation-invariant 18-donor grand mean carrier selection with made_all enforced first.

        Takes no arm-label masks: depends strictly on the pooled 18 donors and symmetric cell-set
        choices, guaranteeing E_{H0}[T] == 0 identically under donor-label exchangeability.
        """
        j = self.genes.index(gene)
        levels = np.where(use[:, None, :], self.levels[1][indices, :, j], self.levels[0][indices, :, j])
        used_g = np.where(use[:, None, :], self.used_levels[1][indices, :, j], self.used_levels[0][indices, :, j])
        exp_g = np.where(use[:, None, :], self.expected_levels[1][indices, :, j], self.expected_levels[0][indices, :, j])
        n_cells = np.where(use[:, None, :], self.n_used_int[indices, :], self.n_all[indices, :])

        mean_used = used_g.mean(axis=1)
        mean_exp = exp_g.mean(axis=1)
        n_present = (n_cells > 0).sum(axis=1)
        n_made_donors = ((used_g > 0.0) & (used_g >= 3.0 * exp_g)).sum(axis=1)
        made_all = (
            (n_present >= self.minimum_present_donors)
            & (mean_used > 0.0)
            & (mean_used >= 3.0 * mean_exp)
            & (n_made_donors >= self.minimum_arm)
        )
        made_all[:, self.hep] = True

        grand_means = levels.mean(axis=1)
        max_gm = np.maximum(grand_means.max(axis=1, keepdims=True), 0.0) + 1.0
        rank_score = np.where(made_all, grand_means, grand_means - max_gm)
        order = np.argsort(-rank_score, axis=1, kind="stable")
        return order[:, 0], order[:, 1], made_all, levels

    def evaluate(self, gene, masks, indices, use, carrier, second, levels=None):
        masks = np.asarray(masks, bool)
        j = self.genes.index(gene)
        if levels is None:
            levels = np.where(use[:, None, :], self.levels[1][indices, :, j], self.levels[0][indices, :, j])
        total = np.where(use[:, None, :], self.totals[1][indices], self.totals[0][indices])
        selector = carrier[:, None, None]
        selected = np.take_along_axis(levels, selector, axis=2)[:, :, 0]
        selected_total = np.take_along_axis(total, selector, axis=2)[:, :, 0]
        all_cells = np.take_along_axis(self.levels[0][indices, :, j], selector, axis=2)[:, :, 0]
        second_cells = np.take_along_axis(levels, second[:, None, None], axis=2)[:, :, 0]
        share = np.take_along_axis(self.share[indices], selector, axis=2)[:, :, 0]
        panel = np.divide(selected, selected_total, out=np.zeros_like(selected), where=selected_total > 0)
        values = {
            "within": selected,
            "within_panelshare": panel,
            "within_allcells": all_cells,
            "second_within": second_cells,
            "share": share,
        }
        result = {"carrier": carrier, "second_type": second}
        for name, value in values.items():
            case = (value * masks).sum(1) / masks.sum(1)
            control = (value * ~masks).sum(1) / (~masks).sum(1)
            result[name] = bounded_difference(case, control)
            result[name + "_log2"] = log_ratio(case, control)
            result[name + "_both_zero"] = (case == 0) & (control == 0)
        return result

    def statistics(self, gene, masks, indices=None, use=None, fixed_carrier=None, fixed_second=None):
        masks = np.asarray(masks, bool)
        if use is None:
            indices, use = self.choices(masks, indices)
        sel_carrier, sel_second, made_all, levels = self.select_carriers(gene, indices, use)
        if fixed_carrier is not None:
            eval_carrier = np.full(len(masks), int(fixed_carrier), dtype=int)
            eval_second = np.full(len(masks), int(fixed_second), dtype=int)
        else:
            eval_carrier, eval_second = sel_carrier, sel_second
        res = self.evaluate(gene, masks, indices, use, eval_carrier, eval_second, levels=levels)
        res["selected_carrier"] = sel_carrier
        res["selected_second"] = sel_second
        res["made_all"] = made_all
        return res

    def detect_ctrl_donors(self, gene, obs_mask, obs_use, carrier_idx):
        j = self.genes.index(gene)
        ctrl = ~obs_mask[0]
        use_int = bool(obs_use[0, carrier_idx])
        det = np.nan_to_num(self.a["det_int" if use_int else "det_all"], nan=0.0)[ctrl, carrier_idx, j]
        ok = det >= 0.05
        if carrier_idx != self.hep:
            used_c = self.used_levels[1 if use_int else 0][ctrl, carrier_idx, j]
            exp_c = self.expected_levels[1 if use_int else 0][ctrl, carrier_idx, j]
            ok &= used_c >= 3.0 * exp_c
        return int(ok.sum())


def _runner_up(selected_arr, obs_idx, types):
    non_obs = selected_arr[selected_arr != obs_idx]
    if len(non_obs) == 0:
        return "none", 0.0
    vals, cnts = np.unique(non_obs, return_counts=True)
    top_i = int(np.argmax(cnts))
    return types[int(vals[top_i])], float(cnts[top_i] / len(selected_arr))


def correct_results(output, n_boot=10000, seed=20261027):
    output = Path(output)
    a = dict(np.load(output / "donor_aggregates_selection.npz", allow_pickle=False))
    engine = SelectionReplay(a)
    observed_masks = a["is_mash"][None, :]
    obs_indices, obs_use = engine.choices(observed_masks)
    masks = np.array([
        np.isin(np.arange(engine.n_donors), c)
        for c in itertools.combinations(range(engine.n_donors), int(observed_masks.sum()))
    ])
    p_indices, p_use = engine.choices(masks)

    rng = np.random.default_rng(seed)
    case = np.where(a["is_mash"])[0]
    control = np.where(~a["is_mash"])[0]
    b_indices = np.empty((n_boot, engine.n_donors), int)
    b_indices[:, case] = case[rng.integers(0, len(case), (n_boot, len(case)))]
    b_indices[:, control] = control[rng.integers(0, len(control), (n_boot, len(control)))]
    b_masks = np.broadcast_to(a["is_mash"], b_indices.shape)
    b_indices, b_use = engine.choices(b_masks, b_indices)

    files = ["results_5d.tsv", "results_5e.tsv", "results_calls_supplementary.tsv"]
    tables = {f: pd.read_csv(output / f, sep="\t") for f in files}
    fixed_genes = sorted(set().union(*(set(t.gene) for t in tables.values())))

    results, checks, diagnostics = {}, [], []
    max_abs_null_mean = 0.0
    for i, gene in enumerate(fixed_genes):
        obs = engine.statistics(gene, observed_masks, obs_indices, obs_use)
        obs_c = int(obs["carrier"][0])
        obs_s = int(obs["second_type"][0])

        perm = engine.statistics(gene, masks, p_indices, p_use)
        boot_uncond = engine.statistics(gene, b_masks, b_indices, b_use)
        boot_fixed = engine.statistics(
            gene, b_masks, b_indices, b_use, fixed_carrier=obs_c, fixed_second=obs_s
        )

        same_carrier_mask = boot_uncond["selected_carrier"] == obs_c
        same_second_mask = boot_uncond["selected_second"] == obs_s
        carrier_stab = float(np.mean(same_carrier_mask))
        second_stab = float(np.mean(same_second_mask))
        ru_type, ru_frac = _runner_up(boot_uncond["selected_carrier"], obs_c, engine.types)

        result = {
            "gene": gene,
            "carrier_bootstrap_stability": carrier_stab,
            "carrier_bootstrap_runner_up": ru_type,
            "carrier_bootstrap_runner_up_frac": ru_frac,
            "second_type_bootstrap_stability": second_stab,
        }
        for name in STATS:
            null_mean = float(np.mean(perm[name]))
            max_abs_null_mean = max(max_abs_null_mean, abs(null_mean))
            assert abs(null_mean) < 1e-12, f"Null asymmetry for {gene}/{name}: mean={null_mean}"

            p = float(np.mean(np.abs(perm[name]) >= abs(obs[name][0]) - 1e-12))
            # Primary 95% CI: fixed-carrier conditional bootstrap (no cross-lineage bimodal mixing)
            fix_lo, fix_hi = np.quantile(boot_fixed[name + "_log2"], [0.025, 0.975], method="inverted_cdf")
            # Unconditional selection-replaying bootstrap 95% CI
            unc_lo, unc_hi = np.quantile(boot_uncond[name + "_log2"], [0.025, 0.975], method="inverted_cdf")
            # Conditional-on-same-selection-event bootstrap 95% CI
            ev_mask = same_second_mask if name == "second_within" else same_carrier_mask
            if ev_mask.any():
                ev_lo, ev_hi = np.quantile(
                    boot_uncond[name + "_log2"][ev_mask], [0.025, 0.975], method="inverted_cdf"
                )
            else:
                ev_lo, ev_hi = np.nan, np.nan

            result.update({
                name + "_log2_ratio": float(obs[name + "_log2"][0]),
                name + "_p_perm": p,
                name + "_ci_lo": float(fix_lo),
                name + "_ci_hi": float(fix_hi),
                name + "_same_carrier_ci_lo": float(ev_lo),
                name + "_same_carrier_ci_hi": float(ev_hi),
                name + "_uncond_ci_lo": float(unc_lo),
                name + "_uncond_ci_hi": float(unc_hi),
            })
            diagnostics.append(dict(
                gene=gene,
                carrier=engine.types[obs_c],
                second_type=engine.types[obs_s],
                statistic=name,
                n_permutations=len(masks),
                n_bootstraps=n_boot,
                perm_null_mean=null_mean,
                perm_null_frac_neg=float(np.mean(perm[name] < 0)),
                perm_null_frac_pos=float(np.mean(perm[name] > 0)),
                n_permutation_carrier_switches=int(np.sum(perm["selected_carrier"] != obs_c)),
                n_permutation_both_zero=int(perm[name + "_both_zero"].sum()),
                n_permutation_nonfinite=int((~np.isfinite(perm[name])).sum()),
                bootstrap_same_carrier_fraction=carrier_stab,
                bootstrap_same_second_fraction=second_stab,
                carrier_bootstrap_runner_up=ru_type,
                carrier_bootstrap_runner_up_frac=ru_frac,
                bootstrap_fixed_ci_lo=float(fix_lo),
                bootstrap_fixed_ci_hi=float(fix_hi),
                bootstrap_same_carrier_ci_lo=float(ev_lo),
                bootstrap_same_carrier_ci_hi=float(ev_hi),
                bootstrap_uncond_ci_lo=float(unc_lo),
                bootstrap_uncond_ci_hi=float(unc_hi),
                n_bootstrap_fixed_both_zero=int(boot_fixed[name + "_both_zero"].sum()),
                n_bootstrap_fixed_infinite=int(np.isinf(boot_fixed[name + "_log2"]).sum()),
                n_bootstrap_uncond_both_zero=int(boot_uncond[name + "_both_zero"].sum()),
                n_bootstrap_uncond_infinite=int(np.isinf(boot_uncond[name + "_log2"]).sum()),
            ))

        for f, t in tables.items():
            if gene in set(t.gene):
                old = t.set_index("gene").loc[gene]
                delta = abs(old.within_log2_ratio - result["within_log2_ratio"])
                checks.append(dict(
                    gene=gene,
                    family=f,
                    observed_effect_difference=delta,
                    carrier_identical=old.carrier == engine.types[obs_c],
                    second_type_identical=old.second_type == engine.types[obs_s],
                ))
                assert np.isclose(old.within_log2_ratio, result["within_log2_ratio"], equal_nan=True), (gene, f)
                assert old.carrier == engine.types[obs_c], (gene, f, old.carrier, engine.types[obs_c])
                assert old.second_type == engine.types[obs_s], (gene, f, old.second_type, engine.types[obs_s])
        results[gene] = result
        print(f"selection replay {i + 1}/{len(fixed_genes)}: {gene}", flush=True)

    replay = pd.DataFrame(results.values()).set_index("gene")
    family_summary = {}
    for f, table in tables.items():
        for col in replay.columns:
            table[col] = replay.loc[table.gene, col].to_numpy()
        table["second_within_p"] = table["second_within_p_perm"]

        if f in (files[0], files[1]):
            prim_mask = table["in_primary_within_holm"].astype(bool).to_numpy()
            table["within_holm_primary"] = holm_subset(table["within_p_perm"], prim_mask)
            table["within_holm_all_sens"] = holm(table["within_p_perm"])
            table["primary_holm_family_n"] = int(prim_mask.sum())
            table["all_rows_family_n"] = len(table)
            table["correction_family_n"] = int(prim_mask.sum())
            family_summary[f] = {
                "all_rows_n": int(len(table)),
                "primary_family_n": int(prim_mask.sum()),
                "excluded_genes": table.loc[~prim_mask, "gene"].tolist(),
            }
        else:
            prim_mask = table["in_primary_within_bh"].astype(bool).to_numpy()
            table["within_bh_primary"] = bh_subset(table["within_p_perm"], prim_mask)
            table["within_bh_all_sens"] = bh(table["within_p_perm"])
            table["primary_bh_family_n"] = int(prim_mask.sum())
            table["all_rows_family_n"] = len(table)
            table["correction_family_n"] = int(prim_mask.sum())
            family_summary[f] = {
                "all_rows_n": int(len(table)),
                "primary_family_n": int(prim_mask.sum()),
                "excluded_n": int((~prim_mask).sum()),
            }

        table["inference"] = (
            "permutation_invariant_18donor_grand_mean_made_all_carrier;"
            "exact_donor_label_permutation_replaying_cellset;"
            "fixed_carrier_conditional_bootstrap_CI_primary_with_carrier_stability_and_uncond_CI"
        )
        table.to_csv(output / f, sep="\t", index=False)

    # ---------------------------------------------------------------- HPX REF sensitivity (exact apples-to-apples)
    engine_hpx = SelectionReplay(a, ref_genes=SENS_REF_HPX)
    sens_rows = []
    for fam_file, fam_label in ((files[0], "5D"), (files[1], "5E")):
        tab = tables[fam_file]
        hpx_p_list = []
        hpx_in_prim_list = []
        fam_records = []
        for row in tab.itertuples(index=False):
            gene = row.gene
            obs_h = engine_hpx.statistics(gene, observed_masks, obs_indices, obs_use)
            obs_ch = int(obs_h["carrier"][0])
            perm_h = engine_hpx.statistics(gene, masks, p_indices, p_use)
            assert abs(float(np.mean(perm_h["within"]))) < 1e-12
            boot_h = engine_hpx.statistics(gene, b_masks, b_indices, b_use)
            p_h = float(np.mean(np.abs(perm_h["within"]) >= abs(obs_h["within"][0]) - 1e-12))
            det_ctrl_h = engine_hpx.detect_ctrl_donors(gene, observed_masks, obs_use, obs_ch)
            in_prim_h = (not bool(row.label_defining)) and (det_ctrl_h >= engine_hpx.minimum_arm)

            use_int_p = bool(obs_use[0, engine.types.index(row.carrier)])
            use_int_h = bool(obs_use[0, obs_ch])
            leak_med_no = float(np.median(engine.leaks[1 if use_int_p else 0][:, engine.types.index(row.carrier)]))
            leak_med_hpx = float(np.median(engine_hpx.leaks[1 if use_int_h else 0][:, obs_ch]))
            j_g = engine.genes.index(gene)
            ctrl_mask = ~a["is_mash"].astype(bool)
            mash_mask = a["is_mash"].astype(bool)
            lvl_h = engine_hpx.levels[1 if use_int_h else 0][:, obs_ch, j_g]

            hpx_p_list.append(p_h)
            hpx_in_prim_list.append(in_prim_h)
            fam_records.append(dict(
                family=fam_label,
                gene=gene,
                label_defining=bool(row.label_defining),
                in_primary_within_holm_no_hpx=bool(row.in_primary_within_holm),
                in_primary_within_holm_with_hpx=bool(in_prim_h),
                carrier_no_hpx=row.carrier,
                carrier_with_hpx=engine_hpx.types[obs_ch],
                carrier_same=(row.carrier == engine_hpx.types[obs_ch]),
                detect_ctrl_donors_no_hpx=int(row.detect_ctrl_donors_pass),
                detect_ctrl_donors_with_hpx=int(det_ctrl_h),
                ctrl_mean_no_hpx=float(row.ctrl_mean_percell_carrier),
                ctrl_mean_with_hpx=float(lvl_h[ctrl_mask].mean()),
                mash_mean_no_hpx=float(row.mash_mean_percell_carrier),
                mash_mean_with_hpx=float(lvl_h[mash_mask].mean()),
                within_log2_no_hpx=float(row.within_log2_ratio),
                within_log2_with_hpx=float(obs_h["within_log2"][0]),
                delta_within_log2=float(row.within_log2_ratio - obs_h["within_log2"][0]),
                within_p_perm_no_hpx=float(row.within_p_perm),
                within_p_perm_with_hpx=p_h,
                delta_within_p_perm=float(row.within_p_perm - p_h),
                within_holm_primary_no_hpx=float(row.within_holm_primary) if pd.notna(row.within_holm_primary) else np.nan,
                carrier_stability_no_hpx=float(row.carrier_bootstrap_stability),
                carrier_stability_with_hpx=float(np.mean(boot_h["selected_carrier"] == obs_ch)),
                carrier_leak_median_no_hpx=leak_med_no,
                carrier_leak_median_with_hpx=leak_med_hpx,
                delta_leak_median=leak_med_no - leak_med_hpx,
            ))
        holm_hpx_prim = holm_subset(np.asarray(hpx_p_list, float), np.asarray(hpx_in_prim_list, bool))
        holm_hpx_all = holm(np.asarray(hpx_p_list, float))
        for rec, hp, ha in zip(fam_records, holm_hpx_prim, holm_hpx_all):
            rec["within_holm_primary_with_hpx"] = float(hp) if np.isfinite(hp) else np.nan
            rec["within_holm_all_sens_with_hpx"] = float(ha) if np.isfinite(ha) else np.nan
            rec["sign_same"] = bool(np.sign(rec["within_log2_no_hpx"]) == np.sign(rec["within_log2_with_hpx"]))
            if rec["in_primary_within_holm_no_hpx"] and rec["in_primary_within_holm_with_hpx"]:
                rec["sig_holm005_same"] = bool(
                    (rec["within_holm_primary_no_hpx"] < 0.05) == (rec["within_holm_primary_with_hpx"] < 0.05)
                )
            else:
                rec["sig_holm005_same"] = True
            sens_rows.append(rec)

    pd.DataFrame(sens_rows).to_csv(output / "ref_sensitivity_hpx.tsv", sep="\t", index=False)
    pd.DataFrame(checks).to_csv(output / "selection_observed_replay_checks.tsv", sep="\t", index=False)
    pd.DataFrame(diagnostics).to_csv(output / "selection_resampling_diagnostics.tsv", sep="\t", index=False)
    pd.DataFrame(p_use.mean(0), index=engine.types, columns=["permutation_interior_fraction"]).to_csv(
        output / "selection_cellset_permutation_frequencies.tsv", sep="\t"
    )
    summary = {
        "status": "completed",
        "remediation": "N09_20261118",
        "ref_genes_primary": list(engine.ref_genes),
        "ref_genes_sensitivity_hpx": list(SENS_REF_HPX),
        "carrier_rule": "permutation_invariant_18donor_grand_mean_with_made_all_spillover_and_donor_sufficiency_prefilter",
        "seed": seed,
        "n_permutations": int(len(masks)),
        "n_bootstraps": int(n_boot),
        "n_unique_genes": int(len(fixed_genes)),
        "max_abs_permutation_null_mean": float(max_abs_null_mean),
        "families": family_summary,
        "zero_rule": "bounded_difference=(case-control)/(case+control); both_zero=0; no_resamples_dropped",
        "primary_interval": "fixed_carrier_conditional_donor_bootstrap_percentile_95_replaying_cellset",
        "secondary_intervals": [
            "same_carrier_event_conditional_bootstrap_95",
            "unconditional_selection_replaying_bootstrap_95",
        ],
        "assumption": "donor_labels_exchangeable; age_and_tissue_source_confounding_not_removed",
    }
    (output / "selection_completion.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"selection replay complete; max |E_H0[T]| = {max_abs_null_mean:.3e}", flush=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    correct_results(args.output)
