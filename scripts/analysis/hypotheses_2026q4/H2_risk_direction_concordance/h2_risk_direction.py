"""H2: per anchored gene, does the disease-state direction match the risk-allele expression direction?

Usage: python -I h2_risk_direction.py <result_dir> <run_dir>
Follows Analysis/Hypotheses/H2-20261008T224908Z/PRESPEC_H2.md. All inputs are provisional candidates.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
H26 = ROOT / "Analysis/Hypotheses/H26-20261007T210504Z/run_22092011/h26_gene_direction_table.tsv"
CAND = ROOT / "RNA-seq/results/manuscript_release/candidates"
FITS = {  # label -> fit directory; arm primary_source_controls in each
    "u774": CAND / "bulk-corrected-controls-20261008T162308Z",
    "u738": CAND / "bulk-corrected-controls-20261008T162310Z",
    "r844": CAND / "bulk-corrected-controls-20260923T231939Z",
}
PRIMARY = "u774"
ARM = "primary_source_controls"
SEED = 20260923
B = 100_000
CONTROLS = ["PNPLA3", "HSD17B13", "MTARC1", "GPAM", "APOB", "CIDEB"]  # H26 calibration table only
SIGN_OF_CALL = {"inhibit": 1, "activate": -1}  # PRESPEC_H2 section 3


def load_fit(path):
    d = pd.read_csv(path / ARM / "deg_results.csv")
    d["ensembl"] = d["gene"].astype(str).str.replace(r"\.[0-9]+$", "", regex=True)
    dup = d["ensembl"][d["ensembl"].duplicated(keep=False)]
    dropped = sorted(set(dup))
    d = d[~d["ensembl"].isin(dropped)].copy()
    d["called"] = (d["padj"] < 0.05) & (d["logFC"].abs() > 0.5)
    groups = pd.read_csv(path / ARM / "cohort_group_counts.tsv", sep="\t")
    return d, dropped, groups


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan)
    ci = binomtest(k, n, 0.5).proportion_ci(confidence_level=0.95, method="wilson")
    return (float(ci.low), float(ci.high))


def pooled(flags, label):
    """flags: boolean array, risk_concordant for each gene with a disease call."""
    n = len(flags)
    k = int(flags.sum())
    lo, hi = wilson(k, n)
    p_binom = float(binomtest(k, n, 0.5, alternative="two-sided").pvalue) if n else np.nan
    return dict(fit=label, n_called=n, k_risk_concordant=k, n_risk_discordant=n - k,
                fraction=k / n if n else np.nan, wilson95_low=lo, wilson95_high=hi, binom_p_two_sided=p_binom)


def permutation(G, D, called_mask):
    """Shuffle G across the whole pool; count concordance among called genes. PRESPEC section 5."""
    rng = np.random.default_rng(SEED)
    idx = np.flatnonzero(called_mask)
    Dc = D[idx]
    k_obs = int((G[idx] == Dc).sum())
    ks = np.empty(B, dtype=np.int64)
    for b in range(B):
        Gp = rng.permutation(G)
        ks[b] = int((Gp[idx] == Dc).sum())
    m = ks.mean()
    p = (1 + int((np.abs(ks - m) >= abs(k_obs - m) - 1e-9).sum())) / (B + 1)
    n = len(idx)
    return dict(null_mean_fraction=m / n if n else np.nan,
                null_q025_fraction=float(np.quantile(ks, 0.025)) / n if n else np.nan,
                null_q975_fraction=float(np.quantile(ks, 0.975)) / n if n else np.nan,
                perm_p_two_sided=p, B=B, seed=SEED), ks


def main(result_dir, run_dir):
    result_dir, run_dir = Path(result_dir), Path(run_dir)
    assert (result_dir / "PRESPEC_H2.md").is_file()
    summary = {"prespec": str(result_dir / "PRESPEC_H2.md"), "provisional": True,
               "genetic_input": str(H26), "fits": {k: str(v) for k, v in FITS.items()}}

    t = pd.read_csv(H26, sep="\t", low_memory=False)
    pool = t[t["lineA_status"] == "call"].copy()
    # encoding checks (PRESPEC section 3)
    bad = pool[(pool["lineA_sign"] == 0) | (pool["lineA_call"].map(SIGN_OF_CALL) != pool["lineA_sign"])]
    if len(bad):
        raise SystemExit(f"lineA_call / lineA_sign inconsistent on {len(bad)} rows: {bad['symbol'].tolist()}")
    assert pool["ensembl"].is_unique
    pool["G"] = pool["lineA_sign"].astype(int)
    pool["lineA_ggt_only"] = pool["lineA_ggt_only"].astype(str).str.lower().isin(["true", "1"])
    pool["published_coding_mechanism"] = pool["published_coding_mechanism"].fillna("").astype(str)
    keep = ["ensembl", "symbol", "primary_evidence_class", "lineA_call", "G", "lineA_defining_trait",
            "lineA_directional_traits", "lineA_ggt_only", "published_coding_mechanism", "susie_driving_gwas"]
    out = pool[keep].reset_index(drop=True)
    summary["pool"] = {"n_lineA_call": len(out),
                       "by_class": out["primary_evidence_class"].value_counts().to_dict(),
                       "G_up_inhibit": int((out["G"] == 1).sum()), "G_down_activate": int((out["G"] == -1).sum())}

    rows, null_rows = [], []
    for label, path in FITS.items():
        d, dropped, groups = load_fit(path)
        summary.setdefault("fit_info", {})[label] = {
            "n_genes_tested": int(len(d)), "n_calls_all_genes": int(d["called"].sum()),
            "dropped_duplicate_unversioned_ids": dropped,
            "group_counts": groups.to_dict("records")}
        m = out[["ensembl"]].merge(d[["ensembl", "gene", "logFC", "padj", "called"]], on="ensembl", how="left")
        state = np.where(m["gene"].isna(), "not_tested", np.where(m["called"] == True, "called", "tested_no_call"))  # noqa: E712
        out[f"{label}_gene_versioned"] = m["gene"].values
        out[f"{label}_logFC"] = m["logFC"].values
        out[f"{label}_padj"] = m["padj"].values
        out[f"{label}_state"] = state
        called_mask = state == "called"
        D = np.sign(np.nan_to_num(m["logFC"].values)).astype(int)
        flag = np.where(~called_mask, "", np.where(out["G"].values == D, "risk_concordant", "risk_discordant"))
        out[f"{label}_flag"] = flag
        res = pooled(flag[called_mask] == "risk_concordant", label)
        res["n_pool_tested"] = int((state != "not_tested").sum())
        res["n_pool_not_tested"] = int((state == "not_tested").sum())
        res["n_called_disease_up"] = int((called_mask & (D == 1)).sum())
        res["n_called_disease_down"] = int((called_mask & (D == -1)).sum())
        res["n_called_G_up"] = int((called_mask & (out["G"].values == 1)).sum())
        res["n_called_G_down"] = int((called_mask & (out["G"].values == -1)).sum())
        perm, ks = permutation(out["G"].values, D, called_mask)
        res.update(perm)
        res["analysis"] = "primary" if label == PRIMARY else "robustness"
        rows.append(res)
        null_rows.append(pd.DataFrame({"fit": label, "k_concordant_null": ks}).groupby(["fit", "k_concordant_null"]).size()
                         .rename("n_permutations").reset_index())

    # sensitivity in the primary fit (PRESPEC section 6)
    called = out[f"{PRIMARY}_state"] == "called"
    sens = called & (out["published_coding_mechanism"] == "") & (~out["lineA_ggt_only"])
    ks_ = int((out.loc[sens, f"{PRIMARY}_flag"] == "risk_concordant").sum())
    ns_ = int(sens.sum())
    lo, hi = wilson(ks_, ns_)
    rows.append(dict(fit=PRIMARY, analysis="sensitivity_no_coding_no_ggt_only", n_called=ns_, k_risk_concordant=ks_,
                     n_risk_discordant=ns_ - ks_, fraction=ks_ / ns_ if ns_ else np.nan, wilson95_low=lo,
                     wilson95_high=hi,
                     binom_p_two_sided=float(binomtest(ks_, ns_, 0.5).pvalue) if ns_ else np.nan))
    pooled_tab = pd.DataFrame(rows)
    pooled_tab.to_csv(run_dir / "h2_pooled.tsv", sep="\t", index=False)
    pd.concat(null_rows).to_csv(run_dir / "h2_permutation_null.tsv", sep="\t", index=False)

    out = out.sort_values(["primary_evidence_class", "symbol"]).reset_index(drop=True)
    out.to_csv(run_dir / "h2_gene_flags.tsv", sep="\t", index=False)
    cols = ["symbol", "ensembl", "primary_evidence_class", "lineA_call", "G", "lineA_defining_trait",
            "lineA_ggt_only", "published_coding_mechanism", f"{PRIMARY}_logFC", f"{PRIMARY}_padj",
            "u738_state", "u738_flag", "u738_logFC", "r844_state", "r844_flag", "r844_logFC"]
    for flag in ("risk_concordant", "risk_discordant"):
        sub = out[out[f"{PRIMARY}_flag"] == flag][cols].sort_values(f"{PRIMARY}_padj")
        sub.to_csv(run_dir / f"h2_{flag}_primary.tsv", sep="\t", index=False)
        summary[f"{flag}_primary_symbols"] = sub["symbol"].tolist()

    # flag stability across fits for genes called in the primary fit
    called = out[f"{PRIMARY}_state"] == "called"  # recomputed after the sort above
    stab = out[called][["u738_state", "u738_flag", "r844_state", "r844_flag", f"{PRIMARY}_flag"]].copy()
    summary["stability"] = {
        other: {f"{a}->{b if b else s}": int(n) for (a, b, s), n in
                stab.groupby([f"{PRIMARY}_flag", f"{other}_flag", f"{other}_state"]).size().items()}
        for other in ("u738", "r844")}

    # controls (PRESPEC section 7): recorded, not scored
    ctl = []
    for g in CONTROLS:
        hit = t[t["symbol"] == g]
        r = {"symbol": g, "in_h26_table": len(hit) > 0, "n_h26_rows": len(hit)}
        if len(hit):
            tr = hit.iloc[0]
            r.update(lineA_status=tr["lineA_status"], lineA_call=tr["lineA_call"],
                     published_coding_mechanism=tr["published_coding_mechanism"])
        if g in set(out["symbol"]):
            o = out[out["symbol"] == g].iloc[0]
            for label in FITS:
                r.update({f"{label}_state": o[f"{label}_state"], f"{label}_flag": o[f"{label}_flag"],
                          f"{label}_logFC": o[f"{label}_logFC"], f"{label}_padj": o[f"{label}_padj"]})
        r["expected_H2_flag"] = "none set (PRESPEC section 7)"
        ctl.append(r)
    pd.DataFrame(ctl).to_csv(run_dir / "h2_controls.tsv", sep="\t", index=False)

    summary["pooled"] = pooled_tab.to_dict("records")
    with open(run_dir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print(pooled_tab.to_string())


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
