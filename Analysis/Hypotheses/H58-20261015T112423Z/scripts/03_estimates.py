"""H58 step 03: families C and X, descriptive outputs and RESULT_H58.md.

Reads GATE_H58.tsv (written by 01) and attaches each member's call status.
Call rules and wording follow PRESPEC_H58.md ("Testing families", "Reading").

--debug-permute: debug mode for code development only. Trait-class labels are
permuted within ancestry (seeded) and B is cut to 200; output goes to the
given --out directory and the result file is named RESULT_H58_DEBUG_PERMUTED.md.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import h58_common as h  # noqa: E402
from h58_common import D, E, m  # noqa: E402

X_COLS = {"X1": "promoter_mass", "X2": "abc_mass", "X3": "accessible_mass",
          "X4": "resolved_no_context_mass", "X5": "unresolved_context_mass"}
X_LABEL = {"X1": "promoter-proximal", "X2": "ABC enhancer", "X3": "lineage-accessible",
           "X4": "resolved, no context", "X5": "unresolved"}
W_C = {"EUR": 16 / 26, "AMR": 10 / 26}


def f6(x: float) -> str:
    return "NA" if x is None or not np.isfinite(x) else f"{x:.6f}"


def iv(x: np.ndarray, level: str) -> tuple[float, float]:
    return h.pct(x, level)


def permute_labels(frame: pd.DataFrame, seed_id: str) -> pd.DataFrame:
    g = h.rng(seed_id)
    f = frame.copy()
    for anc, idx in f.groupby("ancestry").groups.items():
        f.loc[idx, "trait_scope"] = g.permutation(f.loc[idx, "trait_scope"].to_numpy())
    return f


def desc_analyse(sub: pd.DataFrame, label: str, n_boot: int) -> dict:
    t = m.cluster_table(m.assign_clusters(sub.copy()), "coding_mass", "arch_total")
    r = m.analyse(t, h.rng(f"H58:desc:{label}"), n_boot)
    r["label"] = label
    return r


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--debug-permute", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    debug = args.debug_permute
    if not debug:
        h.check_prespec(out)
    n_boot = 200 if debug else h.N_BOOT
    gate = pd.read_csv(out / "GATE_H58.tsv", sep="\t").set_index("member")
    extra = json.loads((out / "gate_extra.json").read_text())

    a = h.load_coding()
    mem = h.load_gated_members()
    sites = h.load_sites(out)
    rep = h.repaired_members(mem, sites)
    ctx = h.context_per_set(mem)
    ctx_rep = h.context_per_set(rep)
    if debug:
        a = permute_labels(a, "H58:debug:coding")
        ctx_rep = permute_labels(ctx_rep, "H58:debug:ctx")
        ctx = ctx_rep.copy()

    # ------------------------------------------------------------ family C
    cdes = h.Design(a, ["coding_mass"], "arch_total", h.STRATA)
    W = h.count_matrix(h.rng("H58:C:boot"), len(cdes.units), n_boot)
    fr, _ = cdes.fractions(W)
    ob, _ = cdes.fractions(None)
    fc, oc = fr["coding_mass"], ob["coding_mass"][0]
    k = cdes.k
    f_obs = {(s, x): oc[k(s, x)] for s in h.STRATA for x in (D, E)}
    f_bt = {(s, x): fc[:, k(s, x)] for s in h.STRATA for x in (D, E)}
    d_obs = {s: f_obs[(s, D)] - f_obs[(s, E)] for s in h.STRATA}
    d_bt = {s: f_bt[(s, D)] - f_bt[(s, E)] for s in h.STRATA}
    std_obs = {x: sum(W_C[s] * f_obs[(s, x)] for s in h.STRATA) for x in (D, E)}
    std_bt = {x: sum(W_C[s] * f_bt[(s, x)] for s in h.STRATA) for x in (D, E)}
    C = {
        "C1": dict(point=sum(W_C[s] * d_obs[s] for s in h.STRATA), draws=sum(W_C[s] * d_bt[s] for s in h.STRATA),
                   fd=std_obs[D], fe=std_obs[E], fd_bt=std_bt[D], fe_bt=std_bt[E]),
        "C2": dict(point=d_obs["EUR"], draws=d_bt["EUR"], fd=f_obs[("EUR", D)], fe=f_obs[("EUR", E)],
                   fd_bt=f_bt[("EUR", D)], fe_bt=f_bt[("EUR", E)]),
        "C3": dict(point=d_obs["AMR"], draws=d_bt["AMR"], fd=f_obs[("AMR", D)], fe=f_obs[("AMR", E)],
                   fd_bt=f_bt[("AMR", D)], fe_bt=f_bt[("AMR", E)]),
    }
    crow = []
    for mem_id, c in C.items():
        dr = c["draws"]
        ok = np.isfinite(dr)
        fin = dr[ok]
        g = gate.loc[mem_id]
        sd = float(np.std(fin, ddof=1))
        if not debug and abs(sd - g.boot_sd) > 1e-12:
            raise AssertionError(f"{mem_id}: bootstrap SD {sd} differs from gate {g.boot_sd}")
        lo95, hi95 = iv(fin, "95")
        lo98, hi98 = iv(fin, "98.333")
        with np.errstate(divide="ignore", invalid="ignore"):
            rb = (c["fd_bt"] / c["fe_bt"])[ok]
        rb = rb[np.isfinite(rb)]
        if bool(g.call_eligible):
            if lo98 > 0:
                call = "difference supported"
            elif lo98 <= 0 <= hi98:
                call = "no difference detected"
            else:
                call = "not supported: 98.333% interval below 0"
        else:
            call = "no call: " + str(g.no_call_reason)
        crow.append(dict(member=mem_id, strata=g.strata, fraction_direct=c["fd"], fraction_enzyme=c["fe"],
                         difference=c["point"], lo95=lo95, hi95=hi95, lo98_333=lo98, hi98_333=hi98,
                         ratio=c["fd"] / c["fe"] if c["fe"] > 0 else np.nan,
                         ratio_lo95=float(np.percentile(rb, 2.5)), ratio_hi95=float(np.percentile(rb, 97.5)),
                         n_finite=int(ok.sum()), boot_sd=sd, mde=h.MDE_FACTOR * sd, call=call))
    ctab = pd.DataFrame(crow)
    # descriptive from shared draws
    het = d_bt["EUR"] - d_bt["AMR"]
    het = het[np.isfinite(het)]
    eqw = 0.5 * d_bt["EUR"] + 0.5 * d_bt["AMR"]
    eqw = eqw[np.isfinite(eqw)]
    wcl = {}
    for x in (D, E):
        v = f_bt[("EUR", x)] - f_bt[("AMR", x)]
        v = v[np.isfinite(v)]
        wcl[x] = (f_obs[("EUR", x)] - f_obs[("AMR", x)], *iv(v, "95"))
    shared_desc = pd.DataFrame([
        dict(label="D(EUR) - D(AMR)", point=d_obs["EUR"] - d_obs["AMR"], lo95=iv(het, "95")[0], hi95=iv(het, "95")[1],
             lo98_333=np.nan, hi98_333=np.nan),
        dict(label="equal-weight 0.5 D(EUR) + 0.5 D(AMR)", point=0.5 * d_obs["EUR"] + 0.5 * d_obs["AMR"],
             lo95=iv(eqw, "95")[0], hi95=iv(eqw, "95")[1], lo98_333=iv(eqw, "98.333")[0], hi98_333=iv(eqw, "98.333")[1]),
        dict(label="enzyme EUR - enzyme AMR", point=wcl[E][0], lo95=wcl[E][1], hi95=wcl[E][2],
             lo98_333=np.nan, hi98_333=np.nan),
        dict(label="direct EUR - direct AMR", point=wcl[D][0], lo95=wcl[D][1], hi95=wcl[D][2],
             lo98_333=np.nan, hi98_333=np.nan),
    ])
    # descriptive refits with the producer estimator
    mvp = a.study.str.startswith("MVP_")
    subsets = {
        "nonEUR_pooled": a[a.ancestry != "EUR"],
        "MVP_only": a[mvp],
        "direct_EUR_vs_enzyme_nonEUR": a[((a.trait_scope == D) & (a.ancestry == "EUR")) | ((a.trait_scope == E) & (a.ancestry != "EUR"))],
        "direct_nonEUR_vs_enzyme_EUR": a[((a.trait_scope == D) & (a.ancestry != "EUR")) | ((a.trait_scope == E) & (a.ancestry == "EUR"))],
    }
    drows = []
    for lab, sub in subsets.items():
        r = desc_analyse(sub, lab, n_boot)
        drows.append({k_: r.get(k_) for k_ in ["label", "n_clusters_direct", "n_clusters_enzyme", "n_physical_clusters",
                                                "fraction_direct", "fraction_enzyme", "difference", "diff_boot_lo",
                                                "diff_boot_hi", "ratio", "ratio_boot_lo", "ratio_boot_hi"]})
    for anc in ["AFR", "EAS", "SAS"]:
        sub = m.assign_clusters(a[(a.ancestry == anc) & (a.trait_scope == E)].copy())
        t = m.cluster_table(sub, "coding_mass", "arch_total")
        gq = h.rng(f"H58:desc:enzyme_only_{anc}")
        idx = gq.integers(0, len(t), size=(n_boot, len(t)))
        bt = t.num.to_numpy()[idx].sum(1) / t.den.to_numpy()[idx].sum(1)
        drows.append(dict(label=f"enzyme_only_{anc}", n_clusters_enzyme=len(t),
                          fraction_enzyme=t.num.sum() / t.den.sum(),
                          diff_boot_lo=np.nan, diff_boot_hi=np.nan,
                          ratio_boot_lo=float(np.percentile(bt, 2.5)), ratio_boot_hi=float(np.percentile(bt, 97.5))))
    desc = pd.DataFrame(drows)
    desc = desc.rename(columns={"ratio_boot_lo": "ratio_or_fraction_lo95", "ratio_boot_hi": "ratio_or_fraction_hi95"})
    # EUR leave-one-study-out
    eur = a[a.ancestry == "EUR"]
    loso = []
    for study, grp in eur.groupby("study"):
        scope = ",".join(sorted(grp.trait_scope.unique()))
        sub = eur[eur.study != study]
        t = m.cluster_table(m.assign_clusters(sub.copy()), "coding_mass", "arch_total")
        r = m.analyse(t, h.rng(f"H58:LOSO:{study}"), n_boot)
        loso.append(dict(left_out=study, scope=scope, n_clusters_direct=r["n_clusters_direct"],
                         n_clusters_enzyme=r["n_clusters_enzyme"], difference=r.get("difference"),
                         diff_lo95=r.get("diff_boot_lo"), diff_hi95=r.get("diff_boot_hi"), ratio=r.get("ratio")))
    loso = pd.DataFrame(loso)

    # ------------------------------------------------------------ family X
    xcols = list(X_COLS.values())
    xdes = h.Design(ctx_rep, xcols, "noncoding_pip_mass", h.STRATA)
    Wx = h.count_matrix(h.rng("H58:X:boot"), len(xdes.units), n_boot)
    frx, _ = xdes.fractions(Wx)
    obx, _ = xdes.fractions(None)
    wx = {s: xdes.DEN[:, xdes.k(s, D)].sum() for s in h.STRATA}
    tw = sum(wx.values())
    wx = {s: v / tw for s, v in wx.items()}
    xrow, xdesc = [], []
    for xm, col in X_COLS.items():
        kd, ke = xdes.k("AMR", D), xdes.k("AMR", E)
        dr = frx[col][:, kd] - frx[col][:, ke]
        fin = dr[np.isfinite(dr)]
        pt = obx[col][0, kd] - obx[col][0, ke]
        g = gate.loc[xm]
        sd = float(np.std(fin, ddof=1))
        if not debug and abs(sd - g.boot_sd) > 1e-12:
            raise AssertionError(f"{xm}: bootstrap SD {sd} differs from gate {g.boot_sd}")
        lo95, hi95 = iv(fin, "95")
        lo99, hi99 = iv(fin, "99")
        if bool(g.call_eligible):
            if lo99 > 0:
                call = "difference supported (direct higher)"
            elif hi99 < 0:
                call = "difference supported (enzyme higher)"
            else:
                call = "no difference detected"
        else:
            call = "no call: " + str(g.no_call_reason)
        xrow.append(dict(member=xm, context=X_LABEL[xm], fraction_direct_AMR=obx[col][0, kd],
                         fraction_enzyme_AMR=obx[col][0, ke], difference=pt, lo95=lo95, hi95=hi95,
                         lo99=lo99, hi99=hi99, n_finite=len(fin), boot_sd=sd, call=call))
        # EUR and EUR+AMR standardized (descriptive)
        de = frx[col][:, xdes.k("EUR", D)] - frx[col][:, xdes.k("EUR", E)]
        pe = obx[col][0, xdes.k("EUR", D)] - obx[col][0, xdes.k("EUR", E)]
        st = wx["EUR"] * de + wx["AMR"] * dr
        ps = wx["EUR"] * pe + wx["AMR"] * pt
        for lab, p_, v in [("EUR repaired", pe, de), ("EUR+AMR standardized repaired", ps, st)]:
            v = v[np.isfinite(v)]
            xdesc.append(dict(context=X_LABEL[xm], label=lab, difference=p_, lo95=iv(v, "95")[0], hi95=iv(v, "95")[1]))
    xtab = pd.DataFrame(xrow)
    # pooled repaired (producer estimator) and adopted reproduced
    cpr = m.assign_clusters(ctx_rep.copy())
    cad = m.assign_clusters(ctx.copy())
    for col in xcols:
        r = m.analyse(m.cluster_table(cpr, col, "noncoding_pip_mass"), h.rng(f"H58:ctx_pooled_repaired:{col}"), n_boot)
        ra = m.analyse(m.cluster_table(cad, col, "noncoding_pip_mass"), m.fit_rng(col), n_boot)
        lab = [k_ for k_, v in X_COLS.items() if v == col][0]
        xdesc.append(dict(context=X_LABEL[lab], label="pooled adopted (reproduced)", difference=ra["difference"],
                          lo95=ra["diff_boot_lo"], hi95=ra["diff_boot_hi"],
                          fraction_direct=ra["fraction_direct"], fraction_enzyme=ra["fraction_enzyme"]))
        xdesc.append(dict(context=X_LABEL[lab], label="pooled repaired", difference=r["difference"],
                          lo95=r["diff_boot_lo"], hi95=r["diff_boot_hi"],
                          fraction_direct=r["fraction_direct"], fraction_enzyme=r["fraction_enzyme"]))
    # credible-set-weighted (companion panel) before and after repair
    for lab, frame in [("credible-set-weighted adopted", ctx), ("credible-set-weighted repaired", ctx_rep)]:
        s = frame.groupby("trait_scope")[xcols + ["noncoding_pip_mass"]].sum()
        for xm, col in X_COLS.items():
            xdesc.append(dict(context=X_LABEL[xm], label=lab,
                              fraction_direct=s.loc[D, col] / s.loc[D, "noncoding_pip_mass"],
                              fraction_enzyme=s.loc[E, col] / s.loc[E, "noncoding_pip_mass"]))
    xdesc = pd.DataFrame(xdesc)

    # ------------------------------------------------------------ write tables
    sfx = "_DEBUG_PERMUTED" if debug else ""
    ctab.to_csv(out / f"family_C{sfx}.tsv", sep="\t", index=False)
    xtab.to_csv(out / f"family_X{sfx}.tsv", sep="\t", index=False)
    shared_desc.to_csv(out / f"descriptive_shared_draws{sfx}.tsv", sep="\t", index=False)
    desc.to_csv(out / f"descriptive_refits{sfx}.tsv", sep="\t", index=False)
    loso.to_csv(out / f"EUR_leave_one_study_out{sfx}.tsv", sep="\t", index=False)
    xdesc.to_csv(out / f"context_descriptive{sfx}.tsv", sep="\t", index=False)

    # ------------------------------------------------------------ verdict (prespec "Reading")
    c3 = ctab.set_index("member").loc["C3"]
    g5, g6 = bool(extra["g5"]), bool(extra["g6"])
    het_excl = not (shared_desc.loc[0, "lo95"] <= 0 <= shared_desc.loc[0, "hi95"])
    if c3.call == "difference supported":
        coding_reading = ("C3 SUPPORTED: in MVP AMR (same cohort, LD panel and fine-mapping pipeline, comparable "
                          "coverage) direct-trait credible sets carry more coding plus splice posterior than "
                          "enzyme-trait credible sets. Print D(AMR) with its 98.333% interval as the within-ancestry "
                          "estimate; the pooled 23.34% vs 3.87% only describes the two arms as sampled.")
        if het_excl:
            coding_reading += (" D(EUR) - D(AMR) 95% interval excludes 0: the text may not say the excess holds "
                               "within each ancestry; report the AMR result only.")
        verdict = "C3 difference supported"
    elif c3.call == "no difference detected":
        coding_reading = ("C3 NO DIFFERENCE DETECTED: in MVP AMR no coding excess was detected although a difference "
                          "of the pooled size (0.195) was detectable. The pooled 6.03 is attributed to ancestry, LD "
                          "panel and fine-mapping coverage composition; the coding-architecture sentences are removed "
                          "or rewritten as composition statements.")
        verdict = "C3 no difference detected"
    else:
        coding_reading = ("C3 UNTESTABLE (" + c3.call + "): the available fine-mapping cannot separate trait class "
                          "from ancestry and coverage. Report D*, D(EUR) and D(AMR) as descriptive estimates, drop "
                          "trait-class wording for the 6.03 ratio, and name the missing analysis (SuSiE fine-mapping "
                          "of all 767 UKBB ALT/AST/GGT lead windows and of MVP EUR ALT, AST and NAFLD).")
        verdict = "C3 untestable"
    if g6:
        ctx_reading = ("G6 passed: the sentence that 36.49% of direct-trait and 77.85% of liver-enzyme noncoding "
                       "posterior lacks resolved deposited regulatory context, and its assay-coverage reading, are "
                       "withdrawn. X1-X4 (AMR, repaired) replace the adopted context differences, each with its 99% "
                       "interval and call.")
    else:
        ctx_reading = ("G6 failed: family X stopped; the adopted context differences stay withdrawn and only the "
                       "cause split is reported.")
    cause = extra["cause_split"]

    def tbl(df: pd.DataFrame) -> str:
        return df.to_markdown(index=False, floatfmt=".6f") if hasattr(df, "to_markdown") else df.to_string(index=False)

    try:
        import tabulate  # noqa: F401
        fmt = tbl
    except ImportError:
        fmt = lambda df: "```\n" + df.to_string(index=False) + "\n```"  # noqa: E731
    lines = [
        f"# RESULT_H58{sfx}", "",
        "DEBUG RUN ON PERMUTED LABELS WITH B=200. NOT A RESULT." if debug else
        "Prespecification: PRESPEC_H58.md (sha256 " + h.PRESPEC_SHA256 + ").", "",
        f"## Verdict: {verdict}", "", coding_reading, "", ctx_reading, "",
        "## Gates", "",
        f"- G5 (coding annotation, P2): {'pass' if g5 else 'FAIL'}; G6 (repair fidelity, P3+N2): {'pass' if g6 else 'FAIL'}",
        f"- G7 selection symmetry: EUR {'pass' if extra['g7_pass']['EUR'] else 'FAIL'}, AMR {'pass' if extra['g7_pass']['AMR'] else 'FAIL'}",
        "", fmt(gate.reset_index()[["member", "strata", "g7_pass", "boot_sd", "mde", "g2_pass", "g3_control",
                                    "g3_cover95", "g3_cover_corr", "g3_pass", "n_dropped", "call_eligible",
                                    "no_call_reason"]]), "",
        "## Family C (coding plus canonical splice, direct minus enzyme; Bonferroni over 3, 98.333% percentile intervals; "
        "10,000 shared draws over 81 pooled clusters, fit id H58:C:boot)", "",
        fmt(ctab), "",
        "Weights fixed at w(EUR)=16/26, w(AMR)=10/26. No p-value (prespec: labels not exchangeable).", "",
        "## Family X (repaired noncoding context, AMR, direct minus enzyme; Bonferroni over 5, 99% intervals; "
        "79 pooled clusters, fit id H58:X:boot)", "", fmt(xtab), "",
        "## Descriptive (outside the families)", "", fmt(shared_desc), "", fmt(desc), "",
        "EUR leave-one-study-out:", "", fmt(loso), "",
        f"Context weights (diagnostic): w(EUR)={wx['EUR']:.6f}, w(AMR)={wx['AMR']:.6f}", "",
        fmt(xdesc), "",
        "## Cause split of the adopted unresolved mass", "",
        "```", json.dumps(cause, indent=1, default=float), "```", "",
        "See DIAGNOSTICS_H58.txt, N1_coverage.tsv, S1_power.tsv, G7_coverage.tsv, cause_split_*.tsv.",
    ]
    (out / f"RESULT_H58{sfx}.md").write_text("\n".join(str(x) for x in lines) + "\n")
    print("\n".join(str(x) for x in lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
