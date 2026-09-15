#!/usr/bin/env python3
"""E3 step 3: compare each padding arm to the measured MPRA allele effect.

Spearman rho within a stratum and an arm, never pooled across strata. Uncertainty from a block
bootstrap over the long-range blocks (the resampling unit), 10,000 draws, seed 20260914, with the
native-minus-npad contrast paired inside each draw so the two arms share the draw's blocks.

Outputs (tables/): e3_predicted_effects.tsv, e3_correlations.tsv, e3_predictions.tsv,
e3_summary.json
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

BOOT_SEED = 20260914
# 10,000 is the prespecified draw count; the override exists only so a smoke test on partial data
# can run cheaply, and the value used is written into e3_summary.json.
N_BOOT = int(__import__("os").environ.get("E3_N_BOOT", "10000"))
PRIMARY_CHANNEL = "dnase"
PRIMARY_CONTEXT = "HepG2_control"
PRIMARY_READOUT = "primary"
HYENADNA_REPORTER_RHO = 0.216
CHANNELS = ["atac", "dnase", "rna", "accessibility_mean"]
ARMS = ["native", "npad"]
CONTEXTS = ["HepG2_control", "HepG2_PAOA", "LX2_control", "LX2_TGFb"]
EPS = 1e-6


def load_scores(path: pathlib.Path) -> pd.DataFrame:
    rows = []
    with path.open() as h:
        for line in h:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return pd.DataFrame(rows)


def predicted_effects(sc: pd.DataFrame, readout: str) -> pd.DataFrame:
    """log2((alt + eps) / (ref + eps)) on the fixed readout, per arm and channel."""
    out = sc[["element_id", "in_random_stratum", "in_topeffect_stratum", "long_range_block",
              "mb_block", "contig"]].copy()
    for arm in ARMS:
        for ch in ("atac", "dnase", "rna"):
            a = sc[f"{arm}_alt_{ch}_{readout}"].to_numpy(dtype=float)
            r = sc[f"{arm}_ref_{ch}_{readout}"].to_numpy(dtype=float)
            out[f"{arm}_{ch}_ref"] = r
            out[f"{arm}_{ch}_alt"] = a
            out[f"{arm}_{ch}"] = np.log2((a + EPS) / (r + EPS))
        out[f"{arm}_accessibility_mean"] = out[[f"{arm}_atac", f"{arm}_dnase"]].mean(axis=1)
    return out


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson on midrank-transformed values: identical to scipy's rho, without its p-value cost."""
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 5 or np.unique(x[ok]).size < 3 or np.unique(y[ok]).size < 3:
        return float("nan")
    rx = stats.rankdata(x[ok])
    ry = stats.rankdata(y[ok])
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / denom) if denom > 0 else float("nan")


def macro_fisher_z(x: np.ndarray, ys: list[np.ndarray]) -> float:
    """tanh of the mean Fisher-z Spearman across contexts, the definition the 0.216 figure uses.

    `Analysis/MASLD_Model_Benchmark/config/gse281364_dna_language_seeded_head_campaign.json`
    binds that number to contexts HepG2_control and HepG2_PAOA under this statistic; matching it
    here is what makes the two numbers the same kind of quantity (see prespec AMENDMENT_01).
    """
    zs = []
    for y in ys:
        r = spearman(x, y)
        if not np.isfinite(r):
            return float("nan")
        zs.append(np.arctanh(np.clip(r, -0.999999, 0.999999)))
    return float(np.tanh(float(np.mean(zs))))


def block_bootstrap(df: pd.DataFrame, cols: list[str], label, n_boot: int = N_BOOT,
                    seed: int = BOOT_SEED) -> dict:
    """Resample long-range blocks with replacement; recompute every column's rho in the same draw.

    `label` is one column name (per-context Spearman) or a list of them (macro Fisher-z).
    """
    blocks = df["long_range_block"].to_numpy()
    uniq = np.unique(blocks)
    idx_by_block = {b: np.flatnonzero(blocks == b) for b in uniq}
    labels = [label] if isinstance(label, str) else list(label)
    ys = [df[c].to_numpy(dtype=float) for c in labels]
    xs = {c: df[c].to_numpy(dtype=float) for c in cols}
    rng = np.random.default_rng(seed)
    draws = {c: np.empty(n_boot) for c in cols}
    for i in range(n_boot):
        pick = rng.integers(0, uniq.size, uniq.size)
        take = np.concatenate([idx_by_block[uniq[j]] for j in pick])
        yy = [y[take] for y in ys]
        for c in cols:
            xx = xs[c][take]
            draws[c][i] = spearman(xx, yy[0]) if len(yy) == 1 else macro_fisher_z(xx, yy)
    return {c: draws[c] for c in cols}


def observed_stat(sub: pd.DataFrame, col: str, label) -> float:
    x = sub[col].to_numpy(dtype=float)
    if isinstance(label, str):
        return spearman(x, sub[label].to_numpy(dtype=float))
    return macro_fisher_z(x, [sub[c].to_numpy(dtype=float) for c in label])


def ci(v: np.ndarray) -> tuple[float, float]:
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan"), float("nan")
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def two_sided_p(v: np.ndarray) -> float:
    """Empirical two-sided p that the bootstrap distribution excludes 0, floored at 1/n_draws."""
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan")
    p = 2.0 * min((v <= 0).mean(), (v >= 0).mean())
    return float(max(p, 1.0 / v.size))


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    pe = p[ok]
    n = pe.size
    if n:
        order = np.argsort(pe)
        ranked = pe[order] * n / (np.arange(n) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        out = np.empty(n)
        out[order] = np.minimum(ranked, 1.0)
        q[ok] = out
    return q.tolist()


def main() -> None:
    out_root = pathlib.Path(sys.argv[1])
    tables = out_root / "tables"
    sc = load_scores(out_root / "raw/api_scores.jsonl")
    eff = pd.read_csv(tables / "mpra_allele_effects.tsv", sep="\t")
    meas = eff.pivot_table(index="element_id", columns="context_id", values="allele_effect")
    meas.columns = [f"measured_{c}" for c in meas.columns]

    frames = {}
    for readout in ("primary", "secondary"):
        pe = predicted_effects(sc, readout)
        pe = pe.merge(meas.reset_index(), on="element_id", how="left")
        pe["readout_window"] = readout
        frames[readout] = pe
    allpe = pd.concat(frames.values(), ignore_index=True)
    allpe.to_csv(tables / "e3_predicted_effects.tsv", sep="\t", index=False)

    rows = []
    boot_store = {}
    for readout in ("primary", "secondary"):
        pe = frames[readout]
        for stratum, mask in (("random", pe["in_random_stratum"]),
                              ("topeffect", pe["in_topeffect_stratum"])):
            sub = pe[mask].reset_index(drop=True)
            targets = [(c, f"measured_{c}", "spearman") for c in CONTEXTS
                       if f"measured_{c}" in sub]
            targets.append(("macro_hepg2", ["measured_HepG2_control", "measured_HepG2_PAOA"],
                            "macro_fisher_z_spearman"))
            for context, label, stat_name in targets:
                cols = [f"{arm}_{ch}" for arm in ARMS for ch in CHANNELS]
                draws = block_bootstrap(sub, cols, label)
                for arm in ARMS:
                    for ch in CHANNELS:
                        c = f"{arm}_{ch}"
                        obs = observed_stat(sub, c, label)
                        lo, hi = ci(draws[c])
                        rows.append({"readout_window": readout, "stratum": stratum,
                                     "context": context, "arm": arm, "channel": ch,
                                     "n_elements": int(sub.shape[0]),
                                     "n_blocks": int(sub["long_range_block"].nunique()),
                                     "spearman": obs, "ci_lo": lo, "ci_hi": hi,
                                     "p_boot_two_sided": two_sided_p(draws[c]),
                                     "statistic": stat_name})
                for ch in CHANNELS:
                    d = draws[f"native_{ch}"] - draws[f"npad_{ch}"]
                    obs = (observed_stat(sub, f"native_{ch}", label)
                           - observed_stat(sub, f"npad_{ch}", label))
                    lo, hi = ci(d)
                    rows.append({"readout_window": readout, "stratum": stratum,
                                 "context": context, "arm": "native_minus_npad", "channel": ch,
                                 "n_elements": int(sub.shape[0]),
                                 "n_blocks": int(sub["long_range_block"].nunique()),
                                 "spearman": obs, "ci_lo": lo, "ci_hi": hi,
                                 "p_boot_two_sided": two_sided_p(d),
                                 "statistic": f"paired_delta_{stat_name}"})
                if readout == "primary":
                    boot_store[(stratum, context)] = draws

    cor = pd.DataFrame(rows)
    fam = ((cor["readout_window"] == "primary") & (cor["arm"] != "native_minus_npad")
           & (cor["channel"] != "accessibility_mean") & (cor["context"].isin(CONTEXTS)))
    cor["bh_family"] = np.where(fam, "primary_readout_3channels_x_2arms_x_2strata_x_4contexts", "")
    q = np.full(cor.shape[0], np.nan)
    q[fam.to_numpy()] = bh(cor.loc[fam, "p_boot_two_sided"].tolist())
    cor["bh_q"] = q
    cor.to_csv(tables / "e3_correlations.tsv", sep="\t", index=False)

    def pick(stratum: str, arm: str, context: str = PRIMARY_CONTEXT,
             channel: str = PRIMARY_CHANNEL) -> pd.Series:
        m = ((cor["readout_window"] == PRIMARY_READOUT) & (cor["stratum"] == stratum)
             & (cor["context"] == context) & (cor["arm"] == arm)
             & (cor["channel"] == channel))
        return cor[m].iloc[0]

    preds = []
    e31_rows = []
    for stratum in ("random", "topeffect"):
        d = pick(stratum, "native_minus_npad")
        e31_rows.append((stratum, float(d["spearman"]), float(d["ci_lo"]), float(d["ci_hi"])))
    met31 = all(v > 0 for _, v, _, _ in e31_rows)
    preds.append({"prediction": "E3.1",
                  "statement": ("native-flank correlation with the measured reporter allele effect "
                                "exceeds the N-padded one, both strata (DNASE, HepG2 control, "
                                "primary readout)"),
                  "observed": "; ".join(f"{s}: delta rho {v:+.3f} [{lo:+.3f}, {hi:+.3f}]"
                                        for s, v, lo, hi in e31_rows),
                  "verdict": "met" if met31 else "not met"})

    for stratum in ("random", "topeffect"):
        n = pick(stratum, "native")
        val = float(n["spearman"])
        preds.append({"prediction": f"E3.2 ({stratum})",
                      "statement": "native-flank Spearman between 0.05 and 0.25",
                      "observed": f"rho {val:+.3f} [{float(n['ci_lo']):+.3f}, {float(n['ci_hi']):+.3f}]",
                      "verdict": "met" if 0.05 <= val <= 0.25 else "not met"})

    for stratum in ("random", "topeffect"):
        vals = {arm: float(pick(stratum, arm)["spearman"]) for arm in ARMS}
        macro = {arm: float(pick(stratum, arm, context="macro_hepg2")["spearman"]) for arm in ARMS}
        preds.append({"prediction": f"E3.3 ({stratum})",
                      "statement": f"both arms below the trained HyenaDNA difference head, {HYENADNA_REPORTER_RHO}",
                      "observed": ("; ".join(f"{a} rho {v:+.3f}" for a, v in vals.items())
                                   + " | same-definition macro Fisher-z over the two HepG2 contexts: "
                                   + "; ".join(f"{a} {v:+.3f}" for a, v in macro.items())),
                      "verdict": "met" if all(v < HYENADNA_REPORTER_RHO for v in vals.values())
                                 else "not met"})

    npad_ok = bool(np.isfinite(sc[[c for c in sc.columns if c.startswith("npad_")]]
                               .to_numpy(dtype=float)).any())
    preds.append({"prediction": "E3.4",
                  "statement": "the model API accepts N-padded sequences",
                  "observed": f"npad arm returned finite values for {int(sc.shape[0])} elements"
                              if npad_ok else "npad arm returned no finite values",
                  "verdict": "met" if npad_ok else "not met"})

    pd.DataFrame(preds).to_csv(tables / "e3_predictions.tsv", sep="\t", index=False)

    npad_var = {ch: float(np.nanstd(frames["primary"][f"npad_{ch}"])) for ch in ("atac", "dnase", "rna")}
    nat_var = {ch: float(np.nanstd(frames["primary"][f"native_{ch}"])) for ch in ("atac", "dnase", "rna")}
    pm = frames["primary"]
    arm_agreement = {}
    for ch in ("atac", "dnase", "rna"):
        a, b = pm[f"native_{ch}"].to_numpy(float), pm[f"npad_{ch}"].to_numpy(float)
        ok = np.isfinite(a) & np.isfinite(b)
        arm_agreement[ch] = {
            "spearman_native_vs_npad": spearman(a, b),
            "share_same_sign": float(np.mean(np.sign(a[ok]) == np.sign(b[ok]))) if ok.any() else float("nan"),
            "median_abs_native": float(np.nanmedian(np.abs(a))),
            "median_abs_npad": float(np.nanmedian(np.abs(b))),
        }
    summary = {
        "package": "e3-reporter-context",
        "what_ran": ("296 MPRA paired-SNV oligos scored through the hosted AlphaGenome model API at "
                     "16,384 bp, two padding arms x two alleles = 1,184 calls"),
        "n_elements_scored": int(sc.shape[0]),
        "n_random_stratum": int(sc["in_random_stratum"].sum()),
        "n_topeffect_stratum": int(sc["in_topeffect_stratum"].sum()),
        "n_long_range_blocks": int(sc["long_range_block"].nunique()),
        "primary": {"channel": PRIMARY_CHANNEL, "context": PRIMARY_CONTEXT,
                    "readout_window_bp": 128, "bootstrap_draws": N_BOOT, "bootstrap_seed": BOOT_SEED},
        "predicted_effect_sd_primary_readout": {"native": nat_var, "npad": npad_var},
        "arm_agreement_primary_readout": arm_agreement,
        "reference_number_traced": {
            "hyenadna_delta_ridge_macro_spearman": HYENADNA_REPORTER_RHO,
            "exact_value": 0.21622332457761556,
            "producing_file": ("Analysis/MASLD_Model_Benchmark/executions/"
                               "gse281364-hyenadna-window-ablation-20260905T192600Z/evaluation/arm_metrics.tsv"),
            "definition": ("supervised delta ridge, out of fold, 5 seeds, 1,033 elements, 239 blocks, "
                           "tanh of the mean Fisher-z Spearman over HepG2_control and HepG2_PAOA"),
            "not_like_for_like_with_this_arm": ("this arm is zero-shot, unfitted, 296 elements in 112 "
                                                "long-range blocks"),
        },
        "scale_caveats": [
            "a hosted-API log2 track ratio is not an Atlas quantile and is never pooled with one",
            "the model was trained on genomic, not episomal, context; neither arm simulates a plasmid",
            "a 126-bp insert carries no promoter, backbone or barcode in either arm",
        ],
    }
    json.dump(summary, (tables / "e3_summary.json").open("w"), indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float))
    print(pd.DataFrame(preds).to_string(index=False))


if __name__ == "__main__":
    main()
