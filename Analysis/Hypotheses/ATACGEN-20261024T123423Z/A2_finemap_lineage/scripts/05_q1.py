#!/usr/bin/env python3
"""A2 Q1: share of fine-mapped posterior in lineage chromatin vs matched expectation.

Unit: credible set nested in LD block; each block counts once within a trait class.
Families: Q1-main (5 classes x 5 lineages x {specific, shared}), Q1-contrast (5 x 4), BH.
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from a2common import (LIN, NONHEP, RES, ROOT, SEED, W, add_bins, bh, boot_p, boot_ratio, ci,
                      ld_tertiles, matched_expectation, peak_classes)

MIN_POS = 5  # Amendment A5: units with any observed posterior needed for a bootstrap test

CLASSES_ORDER = ["direct", "liver_enzyme", "sequela", "calib_albumin", "calib_platelet"]
RES.mkdir(exist_ok=True)

tab = pd.read_csv(W / "cs_table.tsv", sep="\t", dtype={"locus": str})
mem = pd.read_parquet(W / "cs_members.parquet")
win = pd.read_parquet(W / "windows.parquet")
pos = pd.read_parquet(W / "pos_annot.parquet")
anc = pd.read_parquet(W / "anc_annot.parquet")

for d in (mem, win):
    d["allele1"] = d.allele1.astype(str).str.upper()
    d["allele2"] = d.allele2.astype(str).str.upper()
    d["chromosome"] = d.chromosome.astype(np.int64)
    d["position"] = d.position.astype(np.int64)
study_anc = mem.drop_duplicates("study").set_index("study").ancestry
win["study"] = win.study.astype(str)
win["locus"] = win.locus.astype(str)
win["ancestry"] = win.study.map(study_anc)
win["wid"] = win.study + "|" + win.locus
mem["wid"] = mem.study + "|" + mem.locus.astype(str)
n_win_raw = len(win)
win = win.merge(pos, on=["chromosome", "position"], how="inner")
n_win_mapped = len(win)
win = win.merge(anc[["ancestry", "chromosome", "position", "allele1", "allele2", "maf", "ld"]],
                on=["ancestry", "chromosome", "position", "allele1", "allele2"], how="left")
classes = peak_classes(win, "in_") + peak_classes(win, "inb_")
edges = ld_tertiles(win)
pd.DataFrame([(a, *e) for a, e in edges.items()], columns=["ancestry", "ld_q1", "ld_q2"]).to_csv(
    W / "ld_tertile_edges.tsv", sep="\t", index=False)
win = add_bins(win, edges)
win = win.reset_index(drop=True)

# members -> window rows
m = mem.merge(win, on=["wid", "chromosome", "position", "allele1", "allele2"], how="left",
              suffixes=("", "_w"), indicator=True)
m["mapped"] = m["_merge"] == "both"
cs_map = m.groupby("cs_uid").apply(lambda g: g.pip_norm[g.mapped].sum(), include_groups=False)
testable = cs_map[cs_map >= 0.95].index
lift_summary = pd.DataFrame({
    "quantity": ["window_rows", "window_rows_mapped_hg38", "cs_kept", "cs_testable",
                 "cs_untestable_lt95_mapped", "member_rows", "member_rows_mapped",
                 "mean_mapped_mass_testable"],
    "value": [n_win_raw, n_win_mapped, len(cs_map), len(testable), int((cs_map < 0.95).sum()),
              len(m), int(m.mapped.sum()), float(cs_map[testable].mean())]})
lift_summary.to_csv(RES / "q1_liftover_testability.tsv", sep="\t", index=False)
print(lift_summary.to_string(), flush=True)

m = m[m.mapped & m.cs_uid.isin(testable)].copy().reset_index(drop=True)
m["w"] = m.pip_norm / m.cs_uid.map(cs_map)
cols_keep = ["cs_uid", "wid", "w", "chromosome", "position", "pos38", "bin_full", "bin_tm",
             "bin_t"] + classes
mm = m[cols_keep]

full = matched_expectation(win, mm, "wid", classes)
tm = matched_expectation(win, mm, "wid", classes, levels=("bin_tm", "bin_t"))
unm = matched_expectation(win, mm, "wid", classes, levels=())
lev = full.level_used.value_counts(normalize=True).rename("fraction").reset_index()
lev.to_csv(RES / "q1_matching_level_used.tsv", sep="\t", index=False)
print(lev.to_string(), flush=True)

blocks = pd.read_csv(ROOT / "GWAS/finemapping/data/ld_ref/1kg_eur/approx_LD_blocks.txt",
                     sep=r"\s+")
tab = tab[tab.cs_uid.isin(testable)].copy()


def block_of(chrom, p):
    b = blocks[blocks.chr == chrom].sort_values("start")
    j = np.searchsorted(b.start.to_numpy(), p, side="right") - 1
    j = np.clip(j, 0, len(b) - 1)
    return f"{chrom}:{b.start.to_numpy()[j]}-{b.stop.to_numpy()[j]}"


tab["block"] = [block_of(c, p) for c, p in zip(tab.lead_chr, tab.lead_pos)]
tab["eur"] = np.where(tab.ancestry == "EUR", "EUR", "nonEUR")
# Amendment A3: coding/splice posterior per CS (tier 1/2 only) for a noncoding-only stratum
arch = pd.read_csv(ROOT / ("GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/"
                           "run_19783321/credible_set_pip_architecture.tsv"), sep="\t",
                   dtype={"locus": str})
arch["cs_uid"] = arch.study + "|" + arch.locus + "|" + arch.susie_cs.astype(float).astype(str)
arch["coding_splice_mass"] = arch.protein_altering_pip_mass + arch.canonical_splice_pip_mass
tab["coding_splice_mass"] = tab.cs_uid.map(arch.set_index("cs_uid").coding_splice_mass)
print("coding annotation joined for", int(tab.coding_splice_mass.notna().sum()), "of",
      int(tab.trait_class.isin(["direct", "liver_enzyme"]).sum()), "tier 1/2 CSs", flush=True)
tab["coding"] = np.select([tab.coding_splice_mass >= 0.5, tab.coding_splice_mass < 0.5],
                          ["coding", "noncoding"], "unannotated")
m["block"] = m.cs_uid.map(tab.set_index("cs_uid").block)
dcols = [f"da_{l}" for l in ("hepatocyte", "stellate")] + \
    [f"dadir_{l}" for l in ("hepatocyte", "stellate")]
m[["cs_uid", "wid", "w", "study", "trait_class", "block", "chromosome", "position", "allele1",
   "allele2", "pos38", "susie_pip", "pip_norm", "tss_dist", "tss_bin", "nearest_gene"] + dcols +
  [f"pk_{l}" for l in LIN] + classes].to_parquet(W / "q1_members_annot.parquet", index=False)
win[["wid", "ancestry", "chromosome", "position", "pos38", "tss_bin"] + dcols + classes].to_parquet(
    W / "q1_windows_annot.parquet", index=False)


def per_cs(df, tagname):
    ob = df[classes].mul(df.w, axis=0).groupby(df.cs_uid).sum().add_prefix("obs_")
    ex = df[[f"exp_{c}" for c in classes]].mul(df.w, axis=0).groupby(df.cs_uid).sum()
    o = ob.join(ex)
    o["method"] = tagname
    return o


cs_full = per_cs(full, "matched_tss_maf_ld")
cs_tm = per_cs(tm, "matched_tss_maf")
cs_unm = per_cs(unm, "unmatched_window")
cs_all = pd.concat([cs_full, cs_tm, cs_unm]).reset_index()
cs_all = cs_all.merge(tab[["cs_uid", "study", "trait", "trait_class", "ancestry", "eur", "block",
                           "lead_variant", "lead_pip", "n_members", "coding_splice_mass", "coding"]],
                      on="cs_uid")
cs_all.to_csv(RES / "q1_cs_level.tsv.gz", sep="\t", index=False)

MAIN = [f"spec_{l}" for l in LIN] + [f"shared_{l}" for l in LIN]
REPORT = MAIN + [f"any_{l}" for l in LIN] + ["Honly", "Nonly", "HN", "anypeak"] + \
    [f"specb_{l}" for l in LIN] + [f"anyb_{l}" for l in LIN]
rows, contr = [], []
for method, cdf in cs_all.groupby("method"):
    for tclass in CLASSES_ORDER:
        for strat, sdf in [("all", cdf[cdf.trait_class == tclass])] + (
                [(e, cdf[(cdf.trait_class == tclass) & (cdf.eur == e)]) for e in ("EUR", "nonEUR")]
                + [("noncoding_cs", cdf[(cdf.trait_class == tclass) & (cdf.coding == "noncoding")])]
                if tclass in ("direct", "liver_enzyme") else []):
            if sdf.empty:
                continue
            units = sdf.groupby("block")[[f"obs_{c}" for c in REPORT] +
                                         [f"exp_{c}" for c in REPORT]].mean()
            r = boot_ratio(units, REPORT, seed=SEED)
            npos = {c: int((units[f"obs_{c}"] > 0).sum()) for c in REPORT}
            for i, c in enumerate(REPORT):
                with np.errstate(divide="ignore"):
                    lr = np.log(r["boot"][:, i])
                lo, hi = ci(r["boot"][:, i])
                rows.append({"method": method, "trait_class": tclass, "stratum": strat,
                             "peak_class": c, "n_cs": len(sdf), "n_blocks": len(units),
                             "n_studies": sdf.study.nunique(),
                             "share_obs": r["obs"][i], "share_exp": r["exp"][i],
                             "ratio": r["ratio"][i], "ratio_lo": lo, "ratio_hi": hi,
                             "n_blocks_obs_pos": npos[c],
                             "p_boot": boot_p(lr) if npos[c] >= MIN_POS else np.nan})
            if strat == "all":
                ih = REPORT.index("spec_hepatocyte")
                for l in NONHEP:
                    il = REPORT.index(f"spec_{l}")
                    with np.errstate(divide="ignore", invalid="ignore"):
                        rr = r["boot"][:, ih] / r["boot"][:, il]
                        lr = np.log(rr)
                        point = r["ratio"][ih] / r["ratio"][il]
                    lo, hi = ci(rr)
                    ok = min(npos["spec_hepatocyte"], npos[f"spec_{l}"]) >= MIN_POS
                    contr.append({"method": method, "trait_class": tclass,
                                  "contrast": f"spec_hepatocyte_vs_spec_{l}",
                                  "ratio_of_ratios": point, "lo": lo, "hi": hi,
                                  "n_blocks_obs_pos_hep": npos["spec_hepatocyte"],
                                  "n_blocks_obs_pos_other": npos[f"spec_{l}"],
                                  "p_boot": boot_p(lr) if ok else np.nan})
res = pd.DataFrame(rows)
res["family"] = ""
res["q_bh"] = np.nan
fam = (res.method == "matched_tss_maf_ld") & (res.stratum == "all") & res.peak_class.isin(MAIN)
res.loc[fam, "family"] = "Q1-main (50 prespecified; untestable rows have p NA)"
# untestable rows enter BH as p = 1 so the family size stays as prespecified
res.loc[fam, "q_bh"] = bh(res.loc[fam, "p_boot"].fillna(1.0))
res.loc[res.p_boot.isna(), "q_bh"] = np.nan
assert fam.sum() == 50, fam.sum()
res.to_csv(RES / "q1_oe_by_class_lineage.tsv", sep="\t", index=False)
con = pd.DataFrame(contr)
con["family"] = ""
con["q_bh"] = np.nan
f2 = con.method == "matched_tss_maf_ld"
con.loc[f2, "family"] = "Q1-contrast (20 prespecified; untestable rows have p NA)"
con.loc[f2, "q_bh"] = bh(con.loc[f2, "p_boot"].fillna(1.0))
con.loc[con.p_boot.isna(), "q_bh"] = np.nan
con.to_csv(RES / "q1_contrasts.tsv", sep="\t", index=False)
show = res[(res.method == "matched_tss_maf_ld") & (res.stratum == "all") &
           res.peak_class.isin(MAIN + ["anypeak"])]
pd.set_option("display.width", 250)
print(show[["trait_class", "peak_class", "n_cs", "n_blocks", "share_obs", "share_exp", "ratio",
            "ratio_lo", "ratio_hi", "q_bh"]].round(4).to_string(), flush=True)
print(con[f2].round(4).to_string(), flush=True)
