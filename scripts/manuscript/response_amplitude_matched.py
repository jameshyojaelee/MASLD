#!/usr/bin/env python3
"""Does inherited liver-trait genetics mark transcriptionally BUFFERED genes?

Background. Colocalized genes are largely not differentially expressed, and the
overlap odds ratio falls monotonically as the |log2FC| floor rises (1.104 with
no floor to 0.395 at >1.00; see figS_deg_threshold_sensitivity). Two very
different things produce that shape:

  BUFFERING   colocalized genes genuinely respond less to disease, because
              inherited risk acts on constitutive, dosage-sensitive genes whose
              expression the network holds near-constant.

  OBSERVABILITY  colocalized genes are high-expression by construction (441/447
              are source-defined eGenes vs 427/1,261 established-state genes),
              and a high |log2FC| floor preferentially admits LOW-expression
              genes whose effect sizes are inflated conditional on passing FDR.
              This alone predicts a monotone decline.

This script separates them. Conditional on expression level, is the response
amplitude of colocalized genes shifted down relative to comparable genes?

Estimand
    E[|shrunk_logFC| | colocalized] - E[|shrunk_logFC| | matched control]
    over the 14,931-gene jointly testable universe, 447 colocalized genes.

Why shrunk_logFC. The ashr posterior shrinks each estimate in proportion to its
standard error, which removes most of the winner's-curse inflation that drives
the artifact. Raw |logFC| is reported as a sensitivity arm and is expected to be
more favourable to the buffering reading for exactly that reason -- it is the
weaker test, not the stronger one.

Primary specification matches on AveExpr ALONE. A second arm adds SE. That arm
is reported as conservative rather than better: SE is a function of the gene's
dispersion, so if buffering is real then SE partly MEDIATES the hypothesis and
matching on it adjusts away part of the signal. Do not read the SE-matched arm
as the corrected version of the primary.

Matching follows the contract in scripts/portal/gene_catalog_observability.py
(L2 propensity on standardized covariates, 0.2 SD logit caliper, deterministic
1:1 nearest neighbour, SMD balance gate). That function is hardwired to the
lncRNA-versus-protein-coding contrast and to a covariate set that does not apply
here, so the specification is reused rather than the code.

Decision rule, fixed before running:
    delta < 0, matched-null p < 0.05, SMD gate passes  -> buffering supported
    delta ~ 0                                          -> observability artifact
    SMD gate fails                                     -> untestable, claim neither

Read-only on all inputs. Writes one timestamped directory.
"""
from __future__ import annotations

import json
import subprocess
import sys
from bisect import bisect_left
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression

PROJ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEG = PROJ / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"
EVIDENCE = PROJ / "RNA-seq/results/manuscript_release/2026-07-15-r2/evidence_class_table.tsv"
SUSIE = (PROJ / "results/remediation/bg001/bg001-fragment-v211-gencode49-20260807T195243Z"
         / "frozen_sets/finite_susie_447.tsv")
OUTROOT = PROJ / "RNA-seq/results/response_amplitude"

# Contract constants, matching gene_catalog_observability.py.
MATCH_SEED = 20260811
CALIPER_SD = 0.2
MAX_ABS_SMD = 0.10
MIN_MATCH_RATE = 0.70
N_BOOT = 10_000
N_NULL_DRAWS = 10_000
N_DECILES = 10

ARMS = {
    "primary_aveexpr": ["AveExpr"],
    "sensitivity_aveexpr_se": ["AveExpr", "SE"],
}
OUTCOMES = {
    "shrunk": "abs_shrunk_logFC",   # primary
    "raw": "abs_logFC",             # sensitivity
}


def smd(a: np.ndarray, b: np.ndarray) -> float:
    """Standardized mean difference, pooled-variance denominator."""
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    denom = np.sqrt((va + vb) / 2.0)
    if not np.isfinite(denom) or denom == 0:
        return 0.0 if np.mean(a) == np.mean(b) else np.inf
    return float((np.mean(a) - np.mean(b)) / denom)


def load_universe() -> pd.DataFrame:
    d = pd.read_csv(DEG)
    d["g"] = d["gene"].astype(str).str.split(".").str[0]

    ev = pd.read_csv(EVIDENCE, sep="\t", usecols=["ensembl_bulk", "joint_testable"])
    ev["g"] = ev["ensembl_bulk"].astype(str).str.split(".").str[0]
    ev["jt"] = ev["joint_testable"].astype(str).str.upper() == "TRUE"
    jt_ids = set(ev.loc[ev["jt"], "g"])

    su = pd.read_csv(SUSIE, sep="\t")
    su["g"] = su["gene"].astype(str).str.split(".").str[0]

    d = d.loc[d["g"].isin(jt_ids)].copy()
    d["genetic"] = d["g"].isin(set(su["g"])).astype(int)
    d["abs_logFC"] = d["logFC"].abs()
    d["abs_shrunk_logFC"] = d["shrunk_logFC"].abs()

    need = ["AveExpr", "SE", "abs_logFC", "abs_shrunk_logFC"]
    before = len(d)
    d = d.dropna(subset=need).copy()
    if len(d) != before:
        print(f"  dropped {before - len(d)} genes with missing covariate/outcome")

    assert len(d) == 14931, f"expected 14931 jointly testable genes, got {len(d)}"
    assert int(d["genetic"].sum()) == 447, \
        f"expected 447 colocalized genes, got {int(d['genetic'].sum())}"
    return d


def match(frame: pd.DataFrame, covariates: list[str]) -> pd.DataFrame:
    """Deterministic 1:1 propensity nearest-neighbour within a logit caliper."""
    X = frame[covariates].to_numpy(dtype=float)
    mu, sd = X.mean(axis=0), X.std(axis=0, ddof=0)
    sd[sd == 0] = 1.0
    design = (X - mu) / sd
    y = frame["genetic"].to_numpy()

    model = LogisticRegression(penalty="l2", solver="lbfgs", max_iter=10_000,
                               tol=1e-10, random_state=MATCH_SEED)
    model.fit(design, y)
    assert int(model.n_iter_[0]) < model.max_iter, "propensity did not converge"

    p = np.clip(model.predict_proba(design)[:, 1], 1e-12, 1 - 1e-12)
    logit = np.log(p / (1 - p))
    f = frame.assign(propensity_logit=logit)
    caliper = CALIPER_SD * float(np.std(logit, ddof=1))

    treated = f.loc[f["genetic"] == 1]
    controls = f.loc[f["genetic"] == 0]

    # Match the hardest treated genes first: those with fewest eligible controls.
    sorted_ctrl = np.sort(controls["propensity_logit"].to_numpy())
    n_in_caliper = [
        int(np.searchsorted(sorted_ctrl, v + caliper, "right")
            - np.searchsorted(sorted_ctrl, v - caliper, "left"))
        for v in treated["propensity_logit"].to_numpy()
    ]
    treated = treated.assign(n_in_caliper=n_in_caliper).sort_values(
        ["n_in_caliper", "g"], kind="mergesort")

    available = sorted(zip(controls["propensity_logit"].to_numpy(),
                           controls["g"].astype(str)))
    pairs = []
    for row in treated.itertuples(index=False):
        if not available:
            break
        target = float(row.propensity_logit)
        pos = bisect_left(available, (target, ""))
        best, best_d = None, np.inf
        for j in (pos - 1, pos):           # nearest below and nearest at/above
            if 0 <= j < len(available):
                dist = abs(available[j][0] - target)
                if dist < best_d:
                    best, best_d = j, dist
        if best is None or best_d > caliper:
            continue
        pairs.append({"treated_g": row.g, "control_g": available[best][1],
                      "logit_distance": best_d})
        available.pop(best)

    return pd.DataFrame(pairs)


def analyse(frame: pd.DataFrame, pairs: pd.DataFrame, covariates: list[str],
            outcome: str, arm: str, tag: str) -> dict:
    idx = frame.set_index("g")
    t = idx.loc[pairs["treated_g"]]
    c = idx.loc[pairs["control_g"]]

    balance = [{"arm": arm, "outcome": tag, "covariate": cov,
                "treated_mean": float(t[cov].mean()),
                "control_mean": float(c[cov].mean()),
                "smd": smd(t[cov].to_numpy(), c[cov].to_numpy())}
               for cov in covariates]
    max_abs_smd = max(abs(b["smd"]) for b in balance)
    match_rate = len(pairs) / int(frame["genetic"].sum())

    diff = t[outcome].to_numpy() - c[outcome].to_numpy()
    delta = float(np.mean(diff))

    rng = np.random.default_rng(MATCH_SEED)
    boot = np.array([np.mean(rng.choice(diff, size=len(diff), replace=True))
                     for _ in range(N_BOOT)])
    ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5])

    wilcox_p = float(stats.wilcoxon(diff, alternative="two-sided").pvalue) \
        if np.any(diff != 0) else 1.0

    gate = "PASS" if (max_abs_smd < MAX_ABS_SMD and match_rate >= MIN_MATCH_RATE) \
        else "FAIL"
    return {
        "arm": arm, "outcome": tag, "n_pairs": len(pairs),
        "match_rate": match_rate, "max_abs_smd": max_abs_smd,
        "balance_gate": gate,
        "treated_mean": float(t[outcome].mean()),
        "control_mean": float(c[outcome].mean()),
        "delta": delta, "ci_low": float(ci_lo), "ci_high": float(ci_hi),
        "wilcoxon_p": wilcox_p,
        "_balance": balance,
    }


def matched_ladder(frame: pd.DataFrame, pairs: pd.DataFrame, arm: str) -> pd.DataFrame:
    """The floor ladder recomputed WITHIN matched pairs.

    The mean comparison above and the ladder are different estimands: two
    distributions can share a mean and differ in tail mass, and the ladder is
    entirely a statement about the tail. So test the tail directly. Each
    colocalized gene is paired with its expression-matched control, both are
    scored DE or not at each floor, and the discordant pairs are tested with
    McNemar -- the correct test for paired binary outcomes. If the unmatched
    ladder's decline is real biology it should survive here; if it is an
    expression artifact the matched odds ratios should sit at 1.
    """
    idx = frame.set_index("g")
    t = idx.loc[pairs["treated_g"]]
    c = idx.loc[pairs["control_g"]]
    rows = []
    for floor in [0.00, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.75, 1.00]:
        td = ((t["padj"] < 0.05) & (t["abs_logFC"] > floor)).to_numpy()
        cd = ((c["padj"] < 0.05) & (c["abs_logFC"] > floor)).to_numpy()
        b = int(np.sum(td & ~cd))     # treated DE, control not
        cc = int(np.sum(~td & cd))    # control DE, treated not
        # Exact binomial McNemar on the discordant pairs.
        p = float(stats.binomtest(b, b + cc, 0.5).pvalue) if (b + cc) else 1.0
        rows.append({"arm": arm, "lfc_floor": floor,
                     "n_treated_de": int(td.sum()), "n_control_de": int(cd.sum()),
                     "discordant_treated_only": b, "discordant_control_only": cc,
                     "matched_or": (b / cc) if cc else np.nan, "mcnemar_p": p})
    return pd.DataFrame(rows)


def decile_null(frame: pd.DataFrame, outcome: str) -> dict:
    """Cross-check: AveExpr-decile-matched random control sets.

    Follows the construction in program_observability_map/t3_vocabulary_lib.R --
    deciles are computed on expression alone, so the draw is outcome-blind.
    """
    f = frame.copy()
    f["decile"] = pd.qcut(f["AveExpr"], N_DECILES, labels=False, duplicates="drop")
    treated = f.loc[f["genetic"] == 1]
    pool = {d: sub[outcome].to_numpy()
            for d, sub in f.loc[f["genetic"] == 0].groupby("decile")}
    need = treated["decile"].value_counts().to_dict()
    for d, k in need.items():
        assert len(pool[d]) >= k, f"decile {d}: {len(pool[d])} controls for {k} needed"

    observed = float(treated[outcome].mean())
    rng = np.random.default_rng(MATCH_SEED)
    draws = np.array([
        float(np.mean(np.concatenate(
            [rng.choice(pool[d], size=k, replace=False) for d, k in need.items()])))
        for _ in range(N_NULL_DRAWS)])

    # Two-sided empirical p with the +1 correction.
    n_extreme = int(np.sum(np.abs(draws - draws.mean()) >= abs(observed - draws.mean())))
    return {
        "observed_mean": observed,
        "null_mean": float(draws.mean()),
        "null_sd": float(draws.std(ddof=1)),
        "delta_vs_null": observed - float(draws.mean()),
        "empirical_p_two_sided": (1 + n_extreme) / (1 + N_NULL_DRAWS),
        "n_draws": N_NULL_DRAWS,
    }


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    out = OUTROOT / stamp
    out.mkdir(parents=True, exist_ok=True)
    print(f"OUT {out}\n")

    frame = load_universe()
    print(f"universe {len(frame)}  colocalized {int(frame['genetic'].sum())}\n")

    # Unmatched contrast, for reference only -- this is the confounded comparison.
    raw_rows = []
    for tag, col in OUTCOMES.items():
        g = frame.loc[frame["genetic"] == 1, col]
        n = frame.loc[frame["genetic"] == 0, col]
        raw_rows.append({"outcome": tag, "genetic_mean": float(g.mean()),
                         "other_mean": float(n.mean()),
                         "delta": float(g.mean() - n.mean()),
                         "mannwhitney_p": float(
                             stats.mannwhitneyu(g, n, alternative="two-sided").pvalue)})
        # Expression itself: the confounder, quantified.
    conf = {"genetic_AveExpr_mean": float(frame.loc[frame["genetic"] == 1, "AveExpr"].mean()),
            "other_AveExpr_mean": float(frame.loc[frame["genetic"] == 0, "AveExpr"].mean()),
            "AveExpr_smd": smd(frame.loc[frame["genetic"] == 1, "AveExpr"].to_numpy(),
                               frame.loc[frame["genetic"] == 0, "AveExpr"].to_numpy())}
    pd.DataFrame(raw_rows).to_csv(out / "unmatched_contrast.tsv", sep="\t", index=False)
    print("UNMATCHED (confounded):")
    print(pd.DataFrame(raw_rows).to_string(index=False))
    print(f"\nconfounder: AveExpr genetic={conf['genetic_AveExpr_mean']:.3f} "
          f"other={conf['other_AveExpr_mean']:.3f}  SMD={conf['AveExpr_smd']:.3f}\n")

    results, balances, ladders = [], [], []
    for arm, covs in ARMS.items():
        pairs = match(frame, covs)
        pairs.to_csv(out / f"pairs_{arm}.tsv", sep="\t", index=False)
        for tag, col in OUTCOMES.items():
            r = analyse(frame, pairs, covs, col, arm, tag)
            balances.extend(r.pop("_balance"))
            results.append(r)
        ladders.append(matched_ladder(frame, pairs, arm))

    res = pd.DataFrame(results)
    res.to_csv(out / "matched_results.tsv", sep="\t", index=False)
    pd.DataFrame(balances).to_csv(out / "balance.tsv", sep="\t", index=False)
    lad = pd.concat(ladders, ignore_index=True)
    lad.to_csv(out / "matched_ladder.tsv", sep="\t", index=False)

    nulls = []
    for tag, col in OUTCOMES.items():
        n = decile_null(frame, col)
        n["outcome"] = tag
        nulls.append(n)
    pd.DataFrame(nulls).to_csv(out / "decile_matched_null.tsv", sep="\t", index=False)

    print("MATCHED (mean amplitude):")
    print(res.to_string(index=False))
    print("\nMATCHED LADDER (tail, the estimand the unmatched ladder actually makes):")
    print(lad.loc[lad["arm"] == "primary_aveexpr"].to_string(index=False))
    print("\nDECILE-MATCHED NULL (cross-check):")
    print(pd.DataFrame(nulls).to_string(index=False))

    primary = res[(res["arm"] == "primary_aveexpr") & (res["outcome"] == "shrunk")].iloc[0]
    if primary["balance_gate"] == "FAIL":
        verdict = "UNTESTABLE (balance gate failed) -- claim neither"
    elif primary["ci_low"] <= 0 <= primary["ci_high"]:
        verdict = ("NO AMPLITUDE DIFFERENCE -- the ladder is an observability "
                   "artifact; restate orthogonality as a statement about "
                   "expression level, not biology")
    elif primary["delta"] < 0:
        verdict = "BUFFERING SUPPORTED -- colocalized genes respond less than matched genes"
    else:
        verdict = "REVERSED -- colocalized genes respond MORE than matched genes"
    print(f"\nVERDICT ({primary['arm']}/{primary['outcome']}): {verdict}")
    print(f"  delta={primary['delta']:+.4f} "
          f"95% CI [{primary['ci_low']:+.4f}, {primary['ci_high']:+.4f}] "
          f"wilcoxon p={primary['wilcoxon_p']:.3g} "
          f"pairs={primary['n_pairs']} max|SMD|={primary['max_abs_smd']:.4f}")

    (out / "verdict.json").write_text(json.dumps(
        {"verdict": verdict, "primary": primary.to_dict(), "unmatched": raw_rows,
         "confounder": conf, "seed": MATCH_SEED,
         "decision_rule": "delta<0 & CI excludes 0 & SMD gate -> buffering",
         "generated": stamp}, indent=2, default=str))
    (out / "pip_freeze.txt").write_text(
        subprocess.run([sys.executable, "-m", "pip", "freeze"],
                       capture_output=True, text=True).stdout)
    print(f"\nWROTE {out}")


if __name__ == "__main__":
    main()
