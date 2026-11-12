#!/usr/bin/env python3
"""AG2X Stage C analysis (PRESPEC_AG2X_stageC.md section 5). No hosted call, no fitted model, no p value.

Per variant v and donor set s: a_cp, a_oh = SC_ATAC readout (mean over the usable reference-protocol tracks)
under ccl20pos and otherhep conditioning; d_vs = |a_cp| - |a_oh|; D_v = mean over the three sets.
State-dependent (primary, CCL20+ stronger): d_vs > 0 in all three sets and the share of background D strictly
below D_v >= 0.99. Mirror (descriptive, other hepatocytes stronger): d_vs < 0 in all sets and the share of
background (-D) strictly below -D_v >= 0.99. Expected count among members = (share of background variants
meeting the same rule against the background pool) x number of members.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import c01_stageC_lib as S  # noqa: E402


def share_below(x: np.ndarray, pool: np.ndarray) -> np.ndarray:
    pool = np.sort(pool)
    return np.searchsorted(pool, x, side="left") / len(pool)


def main(all_tracks: bool = False):
    if (S.C / "NOT_TESTABLE.json").exists():
        raise SystemExit("Stage C is not testable (NOT_TESTABLE.json); no contrast is computed")
    for f in ["control_a_encode_cells.tsv", "control_b_reference.tsv"]:
        t = pd.read_csv(S.TAB / f, sep="\t")
        assert t["pass"].all(), f"{f}: control failed; no contrast is read"
    tag = "alltracks" if all_tracks else "primary"
    v = pd.read_csv(S.C / "inputs/variants_stageC.tsv", sep="\t").set_index("uid")
    A = {}
    for s in S.SETS:
        for st in S.SCORED_STATES:
            nm = S.cond_name(s, st)
            A[(s, st)] = S.readout(S.read_dir(S.RAW / nm), nm, all_tracks=all_tracks)
    W = pd.DataFrame({f"a_{st}_{s}": A[(s, st)] for s in S.SETS for st in S.SCORED_STATES}).reindex(v.index)
    n_missing = int(W.isna().any(axis=1).sum())
    W = W.dropna()
    for s in S.SETS:
        W[f"d_{s}"] = W[f"a_ccl20pos_{s}"].abs() - W[f"a_otherhep_{s}"].abs()
    d = W[[f"d_{s}" for s in S.SETS]].to_numpy()
    W["D"] = d.mean(1)
    W["all_pos"] = (d > 0).all(1)
    W["all_neg"] = (d < 0).all(1)
    W["role"] = v.loc[W.index, "stageC_role"]
    bg = W[W.role == "background"]
    W["pct_pos"] = share_below(W.D.to_numpy(), bg.D.to_numpy())
    W["pct_neg"] = share_below(-W.D.to_numpy(), -bg.D.to_numpy())
    W["state_dependent_ccl20pos_stronger"] = W.all_pos & (W.pct_pos >= S.CALL_PCT)
    W["state_dependent_otherhep_stronger"] = W.all_neg & (W.pct_neg >= S.CALL_PCT)
    W["label"] = "predicted; imputed; candidate_not_adopted"
    mem, bgw = W[W.role == "member"], W[W.role == "background"]
    res = {"readout": tag, "n_members": int(len(mem)), "n_background": int(len(bgw)), "n_missing_dropped": n_missing}
    for k in ["state_dependent_ccl20pos_stronger", "state_dependent_otherhep_stronger"]:
        f_bg = float(bgw[k].mean())
        res[k] = {"members_called": int(mem[k].sum()), "background_called": int(bgw[k].sum()),
                  "background_share": f_bg, "expected_members_from_background": f_bg * len(mem)}
    res["sign_consistency_share_members_all_pos"] = float(mem.all_pos.mean())
    res["sign_consistency_share_members_all_neg"] = float(mem.all_neg.mean())
    res["sign_consistency_share_background_all_pos"] = float(bgw.all_pos.mean())
    res["sign_consistency_share_background_all_neg"] = float(bgw.all_neg.mean())
    desc = []
    for s in S.SETS:
        cp, oh = W[f"a_ccl20pos_{s}"], W[f"a_otherhep_{s}"]
        nz = (cp != 0) & (oh != 0)
        desc.append({"set": s, "spearman_ccl20pos_vs_otherhep": spearmanr(cp, oh)[0],
                     "alt_minus_ref_sign_agreement": float((np.sign(cp[nz]) == np.sign(oh[nz])).mean()),
                     "median_abs_d": float(W[f"d_{s}"].abs().median()), "mean_d": float(W[f"d_{s}"].mean())})
    pd.DataFrame(desc).to_csv(S.TAB / f"stageC_descriptive_{tag}.tsv", sep="\t", index=False)
    roles = pd.read_csv(S.PKG / "stageB/inputs/variant_roles.tsv", sep="\t").drop_duplicates("uid").set_index("uid")
    W = W.join(roles[["loci", "max_pip", "coding_pip05"]], how="left")
    W.reset_index().rename(columns={"index": "uid"}).to_csv(S.TAB / f"stageC_variant_table_{tag}.tsv", sep="\t", index=False)
    W[W.state_dependent_ccl20pos_stronger | W.state_dependent_otherhep_stronger].reset_index().rename(
        columns={"index": "uid"}).to_csv(S.TAB / f"stageC_calls_{tag}.tsv", sep="\t", index=False)
    json.dump(res, open(S.TAB / f"stageC_summary_{tag}.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(all_tracks="--all-tracks" in sys.argv)
