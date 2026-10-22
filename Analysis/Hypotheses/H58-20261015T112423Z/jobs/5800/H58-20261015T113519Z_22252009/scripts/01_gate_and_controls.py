"""H58 step 01: gates G1-G7, positive controls P1-P3, negative controls N1a/b/c,
N1x, N2, the S1 power check and the prespecified diagnostics.

Prints bootstrap SDs and coverage counts only. Computes no point estimate,
interval or ratio for family members C1-C3 or X1-X5 (PRESPEC_H58.md, "Gate
computed before the outcome"). Exit 2 on a G1, G4 or G7-count stop.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import h58_common as h  # noqa: E402
from h58_common import D, E, m  # noqa: E402

LOG: list[str] = []


def say(*parts) -> None:
    line = " ".join(str(p) for p in parts)
    print(line, flush=True)
    LOG.append(line)


def stop(out: Path, check: str, detail: str) -> None:
    (out / "DIAGNOSTICS_H58.txt").write_text("\n".join(LOG) + "\n")
    (out / "RESULT_H58.md").write_text(
        "# RESULT_H58\n\nVerdict: STOPPED BY RULE before any estimate.\n\n"
        f"Failing check: {check}\n\n{detail}\n\nNothing further was computed "
        "(PRESPEC_H58.md, stop rules).\n")
    say(f"STOP {check}: {detail}")
    sys.exit(2)


def tf(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower() == "true"


# ------------------------------------------------------------------ expected counts (prespec)
EXP_CODING_INST = {(D, "EUR"): 106, (D, "AMR"): 95, (D, "AFR"): 5, (D, "EAS"): 1, (D, "SAS"): 0,
                   (E, "EUR"): 82, (E, "AMR"): 715, (E, "AFR"): 1664, (E, "EAS"): 3694,
                   (E, "SAS"): 53}
EXP_CODING_CLU = {"pooled": (19, 228, 233, 14), "EUR": (16, 14, 29, 1), "AMR": (10, 63, 67, 6),
                  "AFR": (1, 72, None, None), "EAS": (1, 139, None, None), "SAS": (0, 19, None, None),
                  "nonEUR": (10, 224, 225, 9), "MVP": (10, 115, 118, 7)}
EXP_CTX = {"pooled": (190, 5968, 18, 223, 227), "EUR": (94, 76, 15, 14, None),
           "AMR": (95, 676, 10, 62, None), "AFR": (0, 1628, 0, 70, None),
           "EAS": (1, 3535, 1, 136, None), "SAS": (0, 53, 0, 19, None)}
EXP_UNITS = {"coding": 81, "context": 79}
EXP_MEMBER_ROWS, EXP_SITES, EXP_RESOLVED_ROWS = 37112, 16324, 10670
EXP_G7 = {("EUR", D): (47, 47), ("EUR", E): (8, 767), ("AMR", D): (10, 37), ("AMR", E): (84, 271)}

X_COLS = {"X1": "promoter_mass", "X2": "abc_mass", "X3": "accessible_mass",
          "X4": "resolved_no_context_mass", "X5": "unresolved_context_mass"}


def census(frame: pd.DataFrame) -> tuple[int, int, int, int]:
    c = m.assign_clusters(frame.copy())
    nd = c.loc[c.trait_scope == D, "cluster_id"].nunique()
    ne = c.loc[c.trait_scope == E, "cluster_id"].nunique()
    nphys = c.cluster_id.nunique()
    return nd, ne, nphys, nd + ne - nphys


# ------------------------------------------------------------------ G7
def norm_key(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.rstrip("0").str.rstrip(".")


def g7_table(a: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    reg = pd.read_csv(h.REGISTRY, sep="\t")
    cs = pd.read_csv(h.CREDSETS, usecols=["study", "locus", "chromosome", "position"],
                     dtype={"locus": str})
    cs["k"] = norm_key(cs.locus)
    rows = []
    for (anc, scope, study), _ in a.groupby(["ancestry", "trait_scope", "study"]):
        rr = reg[reg.study_name == study]
        if len(rr) != 1:
            raise AssertionError(f"registry rows for {study}: {len(rr)}")
        lp = h.ROOT / "GWAS/finemapping" / str(rr.iloc[0].leadsnps_path)
        lead = pd.read_csv(lp, sep="\t", dtype={"locus": str}) if lp.exists() else pd.DataFrame()
        cst = cs[cs.study == study]
        keys = set(cst.k)
        n_lead = len(lead)
        nw = nf = nf_pos = 0
        kin = 0
        if n_lead:
            lead["k"] = norm_key(lead.locus)
            kin = len(keys & set(lead.k))
            for ch, b in lead.sort_values(["CHR", "BP"]).groupby("CHR"):
                bp = b.BP.to_numpy()
                wid = np.cumsum(np.concatenate([[True], np.diff(bp) > 1_000_000]))
                cpos = cst.loc[cst.chromosome.astype(str) == str(ch), "position"].to_numpy()
                kk = b.k.to_numpy()
                for w in np.unique(wid):
                    sel = wid == w
                    nw += 1
                    nf += int(bool(set(kk[sel]) & keys))
                    lo, hi = bp[sel].min() - 500_000, bp[sel].max() + 500_000
                    nf_pos += int(((cpos >= lo) & (cpos <= hi)).any())
        rows.append(dict(ancestry=anc, trait_scope=scope, study=study, lead_rows=n_lead,
                         lead_windows=nw, windows_with_cs=nf, windows_with_cs_pos500kb=nf_pos,
                         cs_locus_keys=len(keys), cs_keys_in_lead_file=kin))
    per = pd.DataFrame(rows)
    nonempty = per[per.lead_rows > 0]
    g = nonempty.groupby(["ancestry", "trait_scope"])[
        ["lead_windows", "windows_with_cs", "windows_with_cs_pos500kb"]].sum().reset_index()
    g["coverage"] = g.windows_with_cs / g.lead_windows
    g["coverage_pos500kb"] = g.windows_with_cs_pos500kb / g.lead_windows
    return per, g


# ------------------------------------------------------------------ split controls
def arm_clusters(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    c = m.assign_clusters(frame.copy())
    return c, sorted(c.cluster_id.unique())


def split_frame(c: pd.DataFrame, ids: list[str], perm: np.ndarray, n_d: int, n_e: int,
                stratum: str) -> pd.DataFrame:
    dset = {ids[j] for j in perm[:n_d]}
    eset = {ids[j] for j in perm[n_d:n_d + n_e]}
    sub = c[c.cluster_id.isin(dset | eset)].drop(columns="cluster_id").copy()
    sub["trait_scope"] = np.where(c.loc[sub.index, "cluster_id"].isin(dset), D, E)
    sub["pseudo_stratum"] = stratum
    sub["arm_cluster"] = c.loc[sub.index, "cluster_id"].to_numpy()
    return sub


def boot_diffs(design: h.Design, g: np.random.Generator, weights: dict[str, float],
               cols: list[str]) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    W = h.count_matrix(g, len(design.units))
    fr, _ = design.fractions(W)
    obs, _ = design.fractions(None)
    out, point = {}, {}
    for c in cols:
        d = sum(w * (fr[c][:, design.k(s, D)] - fr[c][:, design.k(s, E)]) for s, w in weights.items())
        p = sum(w * (obs[c][0, design.k(s, D)] - obs[c][0, design.k(s, E)]) for s, w in weights.items())
        out[c] = d[np.isfinite(d)]
        point[c] = float(p)
    return out, point


def run_split_controls(coding: pd.DataFrame, ctx_rep: pd.DataFrame, out: Path) -> pd.DataFrame:
    eas_c, eas_ids = arm_clusters(coding[(coding.ancestry == "EAS") & (coding.trait_scope == E)])
    afr_c, afr_ids = arm_clusters(coding[(coding.ancestry == "AFR") & (coding.trait_scope == E)])
    eas_x, eas_xids = arm_clusters(ctx_rep[(ctx_rep.ancestry == "EAS") & (ctx_rep.trait_scope == E)])
    say(f"N1 arms: EAS enzyme coding clusters {len(eas_ids)} (expected 139); AFR enzyme coding "
        f"clusters {len(afr_ids)} (expected 72); EAS enzyme context clusters {len(eas_xids)} "
        "(expected 136)")
    assert (len(eas_ids), len(afr_ids), len(eas_xids)) == (139, 72, 136)
    rows = []
    s1_rows = []
    for i in range(h.N_SPLIT):
        perm_a = h.rng(f"H58:N1a:split:{i}").permutation(len(eas_ids))
        # N1a mirrors C2 (16 vs 14)
        fa = split_frame(eas_c, eas_ids, perm_a, 16, 14, "EAS")
        des = h.Design(fa, ["coding_mass"], "arch_total", ["EAS"], stratum_col="pseudo_stratum")
        dr, pt = boot_diffs(des, h.rng(f"H58:N1a:boot:{i}"), {"EAS": 1.0}, ["coding_mass"])
        rows.append(dict(control="N1a", split=i, column="coding_mass", point=pt["coding_mass"],
                         n_finite=len(dr["coding_mass"]),
                         cover95=h.covers0(dr["coding_mass"], "95"),
                         cover_corr=h.covers0(dr["coding_mass"], "98.333")))
        # S1-C2 power on the N1a split
        for k in range(1, 6):
            rows_s1 = s1_power(fa, k, f"H58:S1:C2:{k}:{i}", {"EAS": 1.0})
            s1_rows.append(dict(arm="C2", k=k, split=i, **rows_s1))
        # N1b mirrors C3 (10 vs 63)
        perm_b = h.rng(f"H58:N1b:split:{i}").permutation(len(eas_ids))
        fb = split_frame(eas_c, eas_ids, perm_b, 10, 63, "EAS")
        des = h.Design(fb, ["coding_mass"], "arch_total", ["EAS"], stratum_col="pseudo_stratum")
        dr, pt = boot_diffs(des, h.rng(f"H58:N1b:boot:{i}"), {"EAS": 1.0}, ["coding_mass"])
        rows.append(dict(control="N1b", split=i, column="coding_mass", point=pt["coding_mass"],
                         n_finite=len(dr["coding_mass"]),
                         cover95=h.covers0(dr["coding_mass"], "95"),
                         cover_corr=h.covers0(dr["coding_mass"], "98.333")))
        for k in range(1, 4):
            rows_s1 = s1_power(fb, k, f"H58:S1:C3:{k}:{i}", {"EAS": 1.0})
            s1_rows.append(dict(arm="C3", k=k, split=i, **rows_s1))
        # N1c mirrors C1: N1a split i (EAS 16 vs 14) + AFR 10 vs 62, weights 16/26, 10/26
        perm_c = h.rng(f"H58:N1c:splitAFR:{i}").permutation(len(afr_ids))
        fc_afr = split_frame(afr_c, afr_ids, perm_c, 10, 62, "AFR")
        fc = pd.concat([fa, fc_afr], ignore_index=True)
        des = h.Design(fc, ["coding_mass"], "arch_total", ["EAS", "AFR"], stratum_col="pseudo_stratum")
        dr, pt = boot_diffs(des, h.rng(f"H58:N1c:boot:{i}"), {"EAS": 16 / 26, "AFR": 10 / 26},
                            ["coding_mass"])
        rows.append(dict(control="N1c", split=i, column="coding_mass", point=pt["coding_mass"],
                         n_finite=len(dr["coding_mass"]),
                         cover95=h.covers0(dr["coding_mass"], "95"),
                         cover_corr=h.covers0(dr["coding_mass"], "98.333")))
        # N1x mirrors family X (10 vs 62 EAS enzyme context clusters), repaired flags
        perm_x = h.rng(f"H58:N1x:split:{i}").permutation(len(eas_xids))
        fx = split_frame(eas_x, eas_xids, perm_x, 10, 62, "EAS")
        xcols = [X_COLS[x] for x in ["X1", "X2", "X3", "X4"]]
        des = h.Design(fx, xcols, "noncoding_pip_mass", ["EAS"], stratum_col="pseudo_stratum")
        dr, pt = boot_diffs(des, h.rng(f"H58:N1x:boot:{i}"), {"EAS": 1.0}, xcols)
        for c in xcols:
            rows.append(dict(control="N1x", split=i, column=c, point=pt[c], n_finite=len(dr[c]),
                             cover95=h.covers0(dr[c], "95"), cover_corr=h.covers0(dr[c], "99")))
        if i % 20 == 0:
            say(f"  split controls: {i + 1}/{h.N_SPLIT} done")
    res = pd.DataFrame(rows)
    res.to_csv(out / "N1_split_controls.tsv", sep="\t", index=False)
    s1 = pd.DataFrame(s1_rows)
    s1.to_csv(out / "S1_power_splits.tsv", sep="\t", index=False)
    return res, s1


def s1_power(frame: pd.DataFrame, k: int, fit_id: str, weights: dict) -> dict:
    """Convert k pseudo-direct clusters to all coding mass; share with 98.333% lower bound > 0."""
    g = h.rng(fit_id)
    dclu = sorted(frame.loc[frame.trait_scope == D, "arm_cluster"].unique())
    pick = set(np.asarray(dclu)[g.choice(len(dclu), size=k, replace=False)])
    f = frame.copy()
    sel = f.arm_cluster.isin(pick)
    f.loc[sel, "coding_mass"] = f.loc[sel, "arch_total"]
    des = h.Design(f, ["coding_mass"], "arch_total", list(weights), stratum_col="pseudo_stratum")
    dr, pt = boot_diffs(des, g, weights, ["coding_mass"])
    lo, _ = h.pct(dr["coding_mass"], "98.333")
    return dict(point=pt["coding_mass"], lo98=lo, lower_above0=bool(lo > 0))


# ------------------------------------------------------------------ main
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    out = Path(ap.parse_args().out)
    h.check_prespec(out)
    gate_rows = []

    a = h.load_coding()
    mem = h.load_gated_members()
    ctx = h.context_per_set(mem)
    sites = h.load_sites(out)
    rep = h.repaired_members(mem, sites)
    ctx_rep = h.context_per_set(rep)

    # ---- diagnostic 1 written by the sbatch wrapper (input_check.txt)
    say("== Diagnostic 1: input sha256 check: see input_check.txt (G0 passed if this runs)")

    # ---- diagnostic 3: ancestry x LD panel
    reg = pd.read_csv(h.REGISTRY, sep="\t")
    al = a.merge(reg[["study_name", "ld_panel"]].rename(columns={"study_name": "study"}),
                 on="study", how="left")
    say("== Diagnostic 3: ancestry x LD panel (instances)")
    say(pd.crosstab(al.ancestry, al.ld_panel).to_string())
    aliased = (pd.crosstab(al.ancestry, al.ld_panel) > 0).sum(axis=1).max() == 1

    # ---- G1
    say("== G1 / Diagnostic 2: counts")
    g1_fail = []
    inst = a.groupby(["trait_scope", "ancestry"]).size().to_dict()
    for key, exp in EXP_CODING_INST.items():
        got = inst.get(key, 0)
        if got != exp:
            g1_fail.append(f"coding instances {key}: {got} != {exp}")
    say("coding instances", {f"{k[0][:6]}|{k[1]}": v for k, v in sorted(inst.items())})
    subsets = {"pooled": a, "EUR": a[a.ancestry == "EUR"], "AMR": a[a.ancestry == "AMR"],
               "AFR": a[a.ancestry == "AFR"], "EAS": a[a.ancestry == "EAS"],
               "SAS": a[a.ancestry == "SAS"], "nonEUR": a[a.ancestry != "EUR"],
               "MVP": a[a.study.str.startswith("MVP_")]}
    for lab, sub in subsets.items():
        c = census(sub)
        exp = EXP_CODING_CLU[lab]
        say(f"coding clusters {lab:<7} direct/enzyme/physical/shared = {c}  expected {exp}")
        if any(e is not None and e != g for e, g in zip(exp, c)):
            g1_fail.append(f"coding clusters {lab}: {c} != {exp}")
    ctx_sub = {"pooled": ctx, **{s: ctx[ctx.ancestry == s] for s in ["EUR", "AMR", "AFR", "EAS", "SAS"]}}
    for lab, sub in ctx_sub.items():
        nd_i, ne_i = int((sub.trait_scope == D).sum()), int((sub.trait_scope == E).sum())
        c = census(sub)
        got = (nd_i, ne_i, c[0], c[1], c[2])
        exp = EXP_CTX[lab]
        say(f"context {lab:<7} inst d/e, clusters d/e, physical = {got}  expected {exp}")
        if any(e is not None and e != g for e, g in zip(exp, got)):
            g1_fail.append(f"context {lab}: {got} != {exp}")
    n_rows, n_sites = len(mem), mem[["chromosome", "position"]].drop_duplicates().shape[0]
    n_resolved = int((mem.promoter_context_status == "resolved").sum())
    say(f"gated member rows {n_rows} (exp {EXP_MEMBER_ROWS}); unique sites {n_sites} "
        f"(exp {EXP_SITES}); resolved rows {n_resolved} (exp {EXP_RESOLVED_ROWS})")
    for lab, got, exp in [("member rows", n_rows, EXP_MEMBER_ROWS), ("sites", n_sites, EXP_SITES),
                          ("resolved rows", n_resolved, EXP_RESOLVED_ROWS)]:
        if got != exp:
            g1_fail.append(f"{lab}: {got} != {exp}")
    # designs (also give the 81 / 79 unit counts and the nesting check)
    cdes = h.Design(a, ["coding_mass"], "arch_total", h.STRATA)
    xdes = h.Design(ctx_rep, list(X_COLS.values()), "noncoding_pip_mass", h.STRATA)
    xdes_adopted = h.Design(ctx, list(X_COLS.values()), "noncoding_pip_mass", h.STRATA)
    say(f"resampled pooled clusters: coding {len(cdes.units)} of {cdes.n_pooled_total} "
        f"(exp 81 of 233); context {len(xdes.units)} of {xdes.n_pooled_total} (exp 79 of 227); "
        "nesting check passed (every within-stratum cluster in one pooled cluster)")
    if len(cdes.units) != EXP_UNITS["coding"]:
        g1_fail.append(f"coding units {len(cdes.units)} != 81")
    if len(xdes.units) != EXP_UNITS["context"] or len(xdes_adopted.units) != EXP_UNITS["context"]:
        g1_fail.append(f"context units {len(xdes.units)} != 79")
    for des, lab in [(cdes, "coding"), (xdes, "context")]:
        for s in h.STRATA:
            for x in (D, E):
                if des.n_clusters(s, x) < 5:
                    g1_fail.append(f"{lab} {s} {x[:6]} clusters {des.n_clusters(s, x)} < 5")
    # repaired context instance set equals adopted
    if set(ctx_rep.credible_set_uid) != set(ctx.credible_set_uid):
        g1_fail.append("repaired context instance set differs from adopted")
    if not aliased:
        say("WARNING: ancestry and LD panel not fully aliased")
    if g1_fail:
        stop(out, "G1", "; ".join(g1_fail))
    say("G1 PASS")

    # ---- diagnostic 4: weights
    w_c = {s: cdes.DEN[:, cdes.k(s, D)].sum() for s in h.STRATA}
    tot = sum(w_c.values())
    w_c = {s: v / tot for s, v in w_c.items()}
    w_x = {s: xdes.DEN[:, xdes.k(s, D)].sum() for s in h.STRATA}
    totx = sum(w_x.values())
    w_x = {s: v / totx for s, v in w_x.items()}
    say(f"== Diagnostic 4: weights coding w(EUR)={w_c['EUR']:.6f} w(AMR)={w_c['AMR']:.6f} "
        f"(prespec 16/26={16/26:.6f}); context (diagnostic) w(EUR)={w_x['EUR']:.6f} "
        f"w(AMR)={w_x['AMR']:.6f}")
    if abs(w_c["EUR"] - 16 / 26) > 1e-9:
        say("NOTE: coding weight differs from 16/26 by more than 1e-9; fixed prespec weights 16/26, 10/26 used")

    # ---- diagnostic 5, 6
    say("== Diagnostic 5: largest single-study share per arm (instances; clusters holding the study)")
    for s in h.STRATA:
        for x in (D, E):
            sub = cdes.within[s][cdes.within[s].trait_scope == x]
            vc = sub.study.value_counts()
            ncl = sub.cluster_id.nunique()
            top = vc.index[0]
            ncl_top = sub[sub.study == top].cluster_id.nunique()
            say(f"  coding {s} {x[:6]}: {top} {vc.iloc[0]} of {len(sub)} instances; "
                f"in {ncl_top} of {ncl} clusters")
    say("== Diagnostic 6: clusters shared between classes")
    for s in h.STRATA:
        c = census(a[a.ancestry == s])
        say(f"  coding {s}: {c[3]} shared of {c[2]} physical")

    # ---- G7 (diagnostic 7)
    say("== G7 / Diagnostic 7: selection symmetry")
    per, cov = g7_table(a)
    per.to_csv(out / "G7_per_study.tsv", sep="\t", index=False)
    cov.to_csv(out / "G7_coverage.tsv", sep="\t", index=False)
    say(per.to_string(index=False))
    say(cov.to_string(index=False))
    g7_count_fail, g7_pass, g7_text = [], {}, {}
    for (s, x), exp in EXP_G7.items():
        r = cov[(cov.ancestry == s) & (cov.trait_scope == x)]
        got = (int(r.windows_with_cs.iloc[0]), int(r.lead_windows.iloc[0])) if len(r) else (0, 0)
        if got != exp:
            g7_count_fail.append(f"{s} {x[:6]} {got} != {exp}")
    if g7_count_fail:
        stop(out, "G7 count reproduction", "; ".join(g7_count_fail))
    for s in h.STRATA:
        cv = {x: float(cov[(cov.ancestry == s) & (cov.trait_scope == x)].coverage.iloc[0]) for x in (D, E)}
        lo, hi = min(cv.values()), max(cv.values())
        ratio = hi / lo if lo > 0 else np.inf
        g7_pass[s] = bool(lo > 0 and hi <= 2 * lo)
        rd = cov[(cov.ancestry == s) & (cov.trait_scope == D)].iloc[0]
        re_ = cov[(cov.ancestry == s) & (cov.trait_scope == E)].iloc[0]
        g7_text[s] = (f"{s} (enzyme arm fine-mapped at {int(re_.windows_with_cs)} of {int(re_.lead_windows)} "
                      f"lead windows; direct arm {int(rd.windows_with_cs)} of {int(rd.lead_windows)})")
        say(f"G7 {s}: coverage direct {cv[D]:.3f} enzyme {cv[E]:.3f} ratio {ratio:.2f} -> "
            f"{'PASS' if g7_pass[s] else 'FAIL'}")
    cs = pd.read_csv(h.CREDSETS, usecols=["study", "chromosome", "position"])
    for name, (ch, p) in {"PNPLA3": (22, 44324727), "TM6SF2": (19, 19379549)}.items():
        for st in ["UKBB_ALT", "UKBB_AST", "UKBB_GGT"]:
            lp = h.ROOT / "GWAS/finemapping" / reg.loc[reg.study_name == st, "leadsnps_path"].iloc[0]
            ld = pd.read_csv(lp, sep="\t")
            nl = int(((ld.CHR.astype(str) == str(ch)) & ((ld.BP - p).abs() < 500_000)).sum())
            nc = int(((cs.study == st) & (cs.chromosome.astype(str) == str(ch))
                      & ((cs.position - p).abs() < 500_000)).sum())
            say(f"  {name} {ch}:{p} {st}: lead rows within 500 kb {nl}; credible-set rows within 500 kb {nc}")

    # ---- P1 / G4
    say("== P1 (G4): reproduction")
    fam = pd.read_csv(h.FAMILY_TABLE, sep="\t").set_index("test")
    p1_fail = []
    ac = m.assign_clusters(a.copy())
    r = m.analyse(m.cluster_table(ac, "coding_mass", "arch_total"), m.fit_rng("coding_mass"))
    for key, exp in [("fraction_direct", 0.23342751888576008), ("fraction_enzyme", 0.038732412743244765),
                     ("ratio", 6.026671264533389)]:
        say(f"  coding {key} {r[key]!r} vs {exp!r}")
        if abs(r[key] - exp) > 1e-9:
            p1_fail.append(f"coding {key}")
    for key in ["diff_boot_lo", "diff_boot_hi", "ratio_boot_lo", "ratio_boot_hi", "difference"]:
        exp = float(fam.loc["coding_mass", key])
        say(f"  coding {key} {r[key]!r} vs table {exp!r}")
        if abs(r[key] - exp) > 1e-9:
            p1_fail.append(f"coding {key}")
    cens = (r["n_clusters_direct"], r["n_clusters_enzyme"], r["n_physical_clusters"], r["n_shared_clusters"])
    if cens != (19, 228, 233, 14):
        p1_fail.append(f"census {cens}")
    cx = m.assign_clusters(ctx.copy())
    for c in m.CTX_COLS:
        rc = m.analyse(m.cluster_table(cx, c, "noncoding_pip_mass"), m.fit_rng(c))
        exp = float(fam.loc[c, "difference"])
        say(f"  context {c} difference {rc['difference']!r} vs table {exp!r}")
        if abs(rc["difference"] - exp) > 1e-9:
            p1_fail.append(f"context {c}")
    for lab, sub, exp in [("EUR", a[a.ancestry == "EUR"], 4.319421760480571),
                          ("nonEUR", a[a.ancestry != "EUR"], 1.9004120095112411)]:
        t = m.cluster_table(m.assign_clusters(sub.copy()), "coding_mass", "arch_total")
        fd = t[t.trait_scope == D]
        fe = t[t.trait_scope == E]
        ratio = (fd.num.sum() / fd.den.sum()) / (fe.num.sum() / fe.den.sum())
        say(f"  within-subset {lab} ratio {ratio!r} vs review {exp!r}")
        if abs(ratio - exp) > 1e-9:
            p1_fail.append(f"{lab} ratio")
    if p1_fail:
        stop(out, "G4 (P1 reproduction)", "; ".join(p1_fail))
    say("G4 PASS")

    # ---- P2 / G5
    say("== P2 (G5): coding annotation at PNPLA3 / TM6SF2")
    p2 = mem[((mem.chromosome == "22") & (mem.position == 44324727))
             | ((mem.chromosome == "19") & (mem.position == 19379549))]
    ok = p2.consequence.astype(str).str.contains("missense_variant") & (p2.consequence_category == "protein_altering")
    say(f"  rows {len(p2)} (exp 98); passing {int(ok.sum())}")
    say(p2.groupby(["chromosome", "trait_scope", "ancestry"]).size().to_string())
    g5 = bool(len(p2) == 98 and ok.all())
    say(f"G5 {'PASS' if g5 else 'FAIL'}")

    # ---- P3, N2 / G6 ; diagnostic 10
    say("== Diagnostic 10: liftOver map counts over unique sites")
    say(f"  0: {int((sites.n_map == 0).sum())}  1: {int((sites.n_map == 1).sum())}  "
        f">1: {int((sites.n_map > 1).sum())}  total {len(sites)}")
    mm = mem.merge(sites, on=["chromosome", "position"], how="left", validate="m:1")
    mm.index = mem.index
    chg = pd.DataFrame({
        "ancestry": rep.ancestry, "trait_scope": rep.trait_scope,
        "promoter_changed": tf(mem.promoter_proximal).to_numpy() != tf(rep.promoter_proximal).to_numpy(),
        "accessible_changed": tf(mem.lineage_accessible).to_numpy() != tf(rep.lineage_accessible).to_numpy(),
        "abc_changed": tf(mem.abc_enhancer_overlap).to_numpy() != tf(rep.abc_enhancer_overlap).to_numpy(),
        "status_changed": (mem.promoter_context_status.to_numpy() != rep.promoter_context_status.to_numpy())})
    say("  rows whose flag changes after repair, by ancestry and class:")
    say(chg.groupby(["ancestry", "trait_scope"])[["promoter_changed", "accessible_changed", "abc_changed",
                                                   "status_changed"]].sum().to_string())
    resolved_no_map = int(((mem.promoter_context_status == "resolved").to_numpy() & (mm.n_map != 1).to_numpy()).sum())
    say(f"  rows recorded resolved whose site lacks a 1:1 map (kept as recorded): {resolved_no_map}")

    say("== P3: repair fidelity on rows recorded as resolved")
    res = mm[mm.promoter_context_status == "resolved"]
    p3_ok = True
    p3 = {}
    for lab, rec_col, new_col in [("promoter", "promoter_proximal", "promoter_new"),
                                  ("accessible", "lineage_accessible", "accessible_new"),
                                  ("abc", "abc_enhancer_overlap", "abc_new")]:
        rec = tf(res[rec_col]).to_numpy()
        new = (res[new_col].astype(bool) & (res.n_map == 1 if lab != "abc" else True)).to_numpy()
        tab = pd.crosstab(pd.Series(rec, name=f"recorded_{lab}"), pd.Series(new, name=f"reannotated_{lab}"))
        say(tab.to_string())
        agree = float((rec == new).mean())
        recall = float(new[rec].mean()) if rec.any() else np.nan
        prec = float(rec[new].mean()) if new.any() else np.nan
        spec = float((~new[~rec]).mean())
        p3[lab] = dict(agreement=agree, recall=recall, precision=prec, false_kept_false=spec,
                       n=int(len(rec)), n_true_rec=int(rec.sum()), n_true_new=int(new.sum()))
        if lab == "abc":
            okk = (recall >= 0.99 or not rec.any()) and spec >= 0.999
        else:
            okk = agree >= 0.995 and (recall >= 0.99 or not rec.any()) and (prec >= 0.99 or not new.any())
        p3_ok &= bool(okk)
        say(f"  {lab}: n {len(rec)} agreement {agree:.6f} recall {recall:.6f} precision {prec:.6f} "
            f"recorded-false kept false {spec:.6f} -> {'ok' if okk else 'FAIL'}")
    say("== N2: repair leaves resolved EUR rows unchanged")
    n2_ok = True
    eur_res_ids = mem.index[(mem.ancestry == "EUR") & (mem.promoter_context_status == "resolved")]
    say(f"  EUR rows resolved before repair: {len(eur_res_ids)} (exp 1548)")

    def eur_ctx_diffs(members: pd.DataFrame) -> dict:
        ps = m.parse_locus(m.build_context_masses(members), mem)
        ancm = mem.drop_duplicates("credible_set_uid").set_index("credible_set_uid").ancestry
        ps = ps[ps.credible_set_uid.map(ancm) == "EUR"]
        ps = m.assign_clusters(ps)
        outd = {}
        for c in m.CTX_COLS:
            t = m.cluster_table(ps, c, "noncoding_pip_mass")
            fd, fe = t[t.trait_scope == D], t[t.trait_scope == E]
            outd[c] = fd.num.sum() / fd.den.sum() - fe.num.sum() / fe.den.sum()
        return outd
    before = eur_ctx_diffs(mem.loc[eur_res_ids])
    after = eur_ctx_diffs(rep.loc[eur_res_ids])
    n2 = {}
    for c in ["promoter_mass", "abc_mass", "accessible_mass", "resolved_no_context_mass"]:
        d = abs(after[c] - before[c])
        n2[c] = d
        n2_ok &= bool(d <= 0.001)
        say(f"  {c}: |change| {d:.3e} (tolerance 0.001)")
    if len(eur_res_ids) != 1548:
        n2_ok = False
    ub = eur_ctx_diffs(mem[mem.ancestry == "EUR"])
    ua = eur_ctx_diffs(rep[rep.ancestry == "EUR"])
    say("  unrestricted EUR change (diagnostic, not gated): "
        + ", ".join(f"{c} {ua[c] - ub[c]:+.6f}" for c in m.CTX_COLS))
    eu_un = mem[(mem.ancestry == "EUR") & (mem.promoter_context_status == "unresolved_liftover")]
    say(f"  EUR unresolved rows {len(eu_un)} (exp 11), normalized PIP sum {eu_un.normalized_susie_pip.sum():.6f}")
    g6 = bool(p3_ok and n2_ok)
    say(f"G6 {'PASS' if g6 else 'FAIL'} (P3 {'ok' if p3_ok else 'fail'}, N2 {'ok' if n2_ok else 'fail'})")

    # ---- diagnostic 11: cause split
    say("== Diagnostic 11: adopted unresolved share and cause split")
    lk = pd.read_csv(h.LOOKUP, dtype=str, keep_default_na=False)
    lk["site"] = lk.chromosome + ":" + lk.position
    lk["vid"] = lk.site + ":" + lk.allele1 + ":" + lk.allele2
    lk["has_hg38"] = (lk.chr_hg38 != "") & (lk.pos_hg38 != "")
    hg = lk.groupby("vid").has_hg38.any()
    sites_in = set(lk.site)
    um = mem.copy()
    cons = um.consequence.astype(str)
    is_utr = cons.str.contains(r"(?:^|,)[35]_prime_UTR_variant(?:$|,)", regex=True)
    is_syn = cons.str.contains(r"(?:^|,)synonymous_variant(?:$|,)", regex=True)
    um["noncoding_dna"] = (um.consequence_category == "other_noncoding") | (
        (um.consequence_category == "synonymous_or_utr") & is_utr & ~is_syn)
    um["unres"] = um.promoter_context_status == "unresolved_liftover"
    say(f"  promoter and lineage status identical on all rows: "
        f"{bool((mem.promoter_context_status == mem.lineage_accessibility_status).all())}")
    sh = um.groupby(["trait_scope", "ancestry"]).apply(lambda d: pd.Series({
        "rows": len(d), "unres_rows": int(d.unres.sum()),
        "unres_pip_share_all": d.loc[d.unres, "normalized_susie_pip"].sum() / d.normalized_susie_pip.sum(),
        "unres_pip_share_noncoding": (d.loc[d.unres & d.noncoding_dna, "normalized_susie_pip"].sum()
                                      / max(d.loc[d.noncoding_dna, "normalized_susie_pip"].sum(), 1e-300))}))
    say(sh.to_string())
    u = um[um.unres].copy()
    u["site"] = u.chromosome + ":" + u.position.astype(str)
    u["vid"] = u.site + ":" + u.allele1.astype(str) + ":" + u.allele2.astype(str)
    site_in = u.site.isin(sites_in)
    vid_in = u.vid.isin(set(lk.vid))
    vid_hg = u.vid.map(hg).fillna(False).astype(bool)
    u["cause"] = np.select([~site_in, site_in & ~vid_in, vid_in & ~vid_hg],
                           ["a_site_absent", "b_other_allele_order", "c_present_empty_hg38"],
                           default="e_present_with_hg38")
    u = u.merge(sites[["chromosome", "position", "n_map"]], on=["chromosome", "position"], how="left")
    u["d_no_1to1_map"] = u.n_map != 1
    cs_tab = u.groupby(["trait_scope", "ancestry", "cause"]).agg(
        rows=("vid", "size"), pip=("normalized_susie_pip", "sum")).reset_index()
    cs_nc = u[u.noncoding_dna].groupby(["trait_scope", "ancestry", "cause"]).agg(
        rows=("vid", "size"), pip=("normalized_susie_pip", "sum")).reset_index()
    cs_tab.to_csv(out / "cause_split_all_rows.tsv", sep="\t", index=False)
    cs_nc.to_csv(out / "cause_split_noncoding_rows.tsv", sep="\t", index=False)
    say(cs_tab.to_string(index=False))
    tot = u.groupby("cause").normalized_susie_pip.sum()
    tot_nc = u[u.noncoding_dna].groupby("cause").normalized_susie_pip.sum()
    def share(t):
        ab = t.get("a_site_absent", 0) + t.get("b_other_allele_order", 0)
        return ab / (ab + t.get("c_present_empty_hg38", 0))
    cause = dict(pip_by_cause_all=tot.to_dict(), pip_by_cause_noncoding=tot_nc.to_dict(),
                 rows_by_cause=u.cause.value_counts().to_dict(),
                 share_failed_lookup_all=float(share(tot)),
                 share_failed_lookup_noncoding=float(share(tot_nc)),
                 d_rows_without_1to1_map=int(u.d_no_1to1_map.sum()))
    say(f"  cause split: {json.dumps(cause, default=float)}")

    # ---- shared bootstraps: SD and MDE only (G2), finite counts (diagnostic 8)
    say("== G2 / Diagnostic 8: shared bootstrap SD and MDE")
    Wc = h.count_matrix(h.rng("H58:C:boot"), len(cdes.units))
    fr, den = cdes.fractions(Wc)
    fc = fr["coding_mass"]
    dE = fc[:, cdes.k("EUR", D)] - fc[:, cdes.k("EUR", E)]
    dA = fc[:, cdes.k("AMR", D)] - fc[:, cdes.k("AMR", E)]
    members_c = {"C1": 16 / 26 * dE + 10 / 26 * dA, "C2": dE, "C3": dA}
    thr = float(pd.read_csv(h.FAMILY_TABLE, sep="\t").set_index("test").loc["coding_mass", "difference"])
    strata_of = {"C1": ["EUR", "AMR"], "C2": ["EUR"], "C3": ["AMR"]}
    for k_, v in members_c.items():
        fin = v[np.isfinite(v)]
        sd = float(np.std(fin, ddof=1))
        mde = h.MDE_FACTOR * sd
        gate_rows.append(dict(member=k_, family="C", strata="+".join(strata_of[k_]),
                              g7_pass=all(g7_pass[s] for s in strata_of[k_]),
                              n_finite=len(fin), n_dropped=h.N_BOOT - len(fin), boot_sd=sd, mde=mde,
                              mde_threshold=thr, g2_pass=bool(mde <= thr)))
        say(f"  {k_}: finite {len(fin)} of {h.N_BOOT}; SD {sd:.6f}; MDE {mde:.6f} "
            f"(threshold {thr:.6f}, factor {h.MDE_FACTOR:.4f}) -> G2 {'PASS' if mde <= thr else 'FAIL'}")
    Wx = h.count_matrix(h.rng("H58:X:boot"), len(xdes.units))
    frx, _ = xdes.fractions(Wx)
    for xm, col in X_COLS.items():
        v = frx[col][:, xdes.k("AMR", D)] - frx[col][:, xdes.k("AMR", E)]
        fin = v[np.isfinite(v)]
        sd = float(np.std(fin, ddof=1))
        gate_rows.append(dict(member=xm, family="X", strata="AMR", g7_pass=g7_pass["AMR"],
                              n_finite=len(fin), n_dropped=h.N_BOOT - len(fin), boot_sd=sd,
                              mde=h.MDE_FACTOR * sd, mde_threshold=np.nan, g2_pass=np.nan,
                              all_draws_equal=bool(np.ptp(fin) == 0) if len(fin) else True))
        say(f"  {xm} ({col}): finite {len(fin)}; SD {sd:.6g}; MDE {h.MDE_FACTOR * sd:.6g} (printed, not gated)")

    # ---- N1 split controls, S1 (diagnostic 9)
    say("== G3 / Diagnostic 9: random-split negative controls (200 splits each)")
    n1, s1 = run_split_controls(a, ctx_rep, out)
    cov_tab = n1.groupby(["control", "column"]).agg(
        median_point=("point", "median"), cover95=("cover95", "sum"), cover_corr=("cover_corr", "sum"),
        min_finite=("n_finite", "min")).reset_index()
    cov_tab["g3_pass"] = (cov_tab.cover95 >= 180) & (cov_tab.cover_corr >= 190)
    cov_tab["median_within_0.01"] = cov_tab.median_point.abs() <= 0.01
    cov_tab.to_csv(out / "N1_coverage.tsv", sep="\t", index=False)
    say(cov_tab.to_string(index=False))
    s1_tab = s1.groupby(["arm", "k"]).agg(share_lower_above0=("lower_above0", "mean"),
                                          median_shift=("point", "median")).reset_index()
    s1_tab.to_csv(out / "S1_power.tsv", sep="\t", index=False)
    say("S1 power (share of 200 splits with 98.333% lower bound above 0):")
    say(s1_tab.to_string(index=False))
    gov = {"C1": ("N1c", "coding_mass"), "C2": ("N1a", "coding_mass"), "C3": ("N1b", "coding_mass"),
           "X1": ("N1x", "promoter_mass"), "X2": ("N1x", "abc_mass"), "X3": ("N1x", "accessible_mass"),
           "X4": ("N1x", "resolved_no_context_mass")}
    for gr in gate_rows:
        if gr["member"] in gov:
            ctl, col = gov[gr["member"]]
            row = cov_tab[(cov_tab.control == ctl) & (cov_tab.column == col)].iloc[0]
            gr.update(g3_control=ctl, g3_cover95=int(row.cover95), g3_cover_corr=int(row.cover_corr),
                      g3_pass=bool(row.g3_pass))
        else:
            gr.update(g3_control="none (N1x covers X1-X4 only)", g3_cover95=np.nan,
                      g3_cover_corr=np.nan, g3_pass=False)

    # ---- call eligibility
    for gr in gate_rows:
        reasons = []
        if gr["family"] == "C":
            if not g5:
                reasons.append("G5 failed: coding family stopped")
            if not gr["g7_pass"]:
                fails = [s for s in strata_of[gr["member"]] if not g7_pass[s]]
                reasons.append("G7 failed in " + ",".join(g7_text[s] for s in fails))
            if not gr["g2_pass"]:
                reasons.append(f"G2 failed: MDE {gr['mde']:.4f} > {thr:.4f}")
        else:
            if not g6:
                reasons.append("G6 failed: family X stopped")
            if not gr["g7_pass"]:
                reasons.append("G7 failed in " + g7_text["AMR"])
            if gr.get("all_draws_equal"):
                reasons.append("all bootstrap draws equal")
        if not gr["g3_pass"]:
            reasons.append(f"G3 failed or not governed ({gr['g3_control']})")
        if gr["n_dropped"] > 100:
            reasons.append(f"{gr['n_dropped']} dropped draws > 100")
        gr["call_eligible"] = not reasons
        gr["no_call_reason"] = "; ".join(reasons)
    gate = pd.DataFrame(gate_rows)
    gate["G5_pass"] = g5
    gate["G6_pass"] = g6
    gate["G7_EUR_pass"] = g7_pass["EUR"]
    gate["G7_AMR_pass"] = g7_pass["AMR"]
    gate.to_csv(out / "GATE_H58.tsv", sep="\t", index=False)
    say("== GATE_H58.tsv")
    say(gate.drop(columns=["mde_threshold"]).to_string(index=False))
    extra = dict(weights_coding=w_c, weights_context=w_x, cause_split=cause, p3=p3,
                 n2=n2, g5=g5, g6=g6, g7_pass=g7_pass, aliased=bool(aliased))
    (out / "gate_extra.json").write_text(json.dumps(extra, indent=1, default=float))
    (out / "DIAGNOSTICS_H58.txt").write_text("\n".join(LOG) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
