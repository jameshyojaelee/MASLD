#!/usr/bin/env python3
"""Pazoki 2022 PDFF (GCST90267352): are the effect/other allele labels reversed?

Reads the as-downloaded file A and its three derived copies (B hg19 table, C
fine-mapping table, D LDSC munged table) and compares A with two independently
oriented PDFF GWAS (2021_34128465 = R1, 2021_34957434 = R2) on chromosome,
position and unordered alleles. A, C and both references share GRCh37 positions
(A equals the LDSC 1000G reference at rs1260326 and rs641738). Nothing is
modified; all outputs go to the directory given as argv[1].

Usage:
  pazoki_orientation_check.py out_dir                       # as-downloaded files
  pazoki_orientation_check.py out_dir corrected corrected_dir  # corrected copies of A-D
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(sys.argv[1])
assert OUT.is_dir(), "out_dir must exist (created by the sbatch wrapper)"
MODE = sys.argv[2] if len(sys.argv) > 2 else "original"
assert MODE in ("original", "corrected")
G = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS")
P = {
    "A": G / "MR_Data/GCST90267352_PDFF_Pazoki2022.tsv.gz",
    "B": G / "MR_Data/hg19/GCST90267352_PDFF_Pazoki2022_hg19.tsv.gz",
    "C": G / "finemapping/data/sumstats/2022_36402844_PDFF_EUR_preprocessed.tsv",
    "D": G / "ldsc/munged/pdff.sumstats.gz",
    "R1": G / "finemapping/data/sumstats/2021_34128465_PDFF_EUR_preprocessed.tsv",
    "R2": G / "finemapping/data/sumstats/2021_34957434_PDFF_EUR_preprocessed.tsv",
}
if MODE == "corrected":
    CD = Path(sys.argv[3])
    for k in ("A", "B", "C", "D"):
        P[k] = CD / P[k].name
assert all(p.exists() for p in P.values()), [str(p) for p in P.values() if not p.exists()]
THR = [5e-8, 1e-6, 1e-4, 0.05]
LOCI = [("PNPLA3", "rs738409"), ("TM6SF2", "rs58542926"), ("GCKR", "rs1260326"), ("MBOAT7", "rs641738")]
np.random.seed(20260923)


def log(msg):
    print(msg, flush=True)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 23), b""):
            h.update(chunk)
    return h.hexdigest()


with open(OUT / ("md5_before.txt" if MODE == "original" else "md5_corrected.txt"), "w") as fh:
    for k, p in P.items():
        fh.write(f"{md5(p)}  {k}  {p}\n")
log(f"mode {MODE}; md5 written")

# ---- read ------------------------------------------------------------------
STR = {"variant_id": str, "chromosome": str, "effect_allele": str, "other_allele": str}
A = pd.read_csv(P["A"], sep="\t", dtype=STR, low_memory=False)
B = pd.read_csv(P["B"], sep="\t", dtype={**STR, "hg19_chr": str}, low_memory=False)
C = pd.read_csv(P["C"], sep="\t", dtype={"chromosome": str, "allele1": str, "allele2": str}, low_memory=False)
D = pd.read_csv(P["D"], sep="\t", dtype={"SNP": str, "A1": str, "A2": str}, low_memory=False)
R = {k: pd.read_csv(P[k], sep="\t", dtype={"chromosome": str, "allele1": str, "allele2": str}, low_memory=False)
     for k in ("R1", "R2")}
log(f"rows A={len(A)} B={len(B)} C={len(C)} D={len(D)} R1={len(R['R1'])} R2={len(R['R2'])}")
cols = {k: list(v.columns) for k, v in {"A": A, "B": B, "C": C, "D": D, "R1": R["R1"], "R2": R["R2"]}.items()}
eaf_cols = {k: [c for c in v if any(s in c.lower() for s in ("freq", "eaf", "maf", "_af", "af_"))] for k, v in cols.items()}


# ---- standardised tables for alignment --------------------------------------
def finish(d):
    d = d.copy()
    first, second = d["first"], d["second"]
    lo = first.where(first < second, second)
    hi = second.where(first < second, first)
    d["key"] = d["chrom"].astype(str) + ":" + d["pos"].astype("int64").astype(str) + ":" + lo + ":" + hi
    d["pal"] = (first + second).isin(["AT", "TA", "CG", "GC"])
    d = d[~d["key"].duplicated(keep=False)]
    return d[np.isfinite(d["z"])]


def std_a(df):
    d = df.loc[df["p_value"] < 0.05, ["variant_id", "chromosome", "base_pair_location", "effect_allele",
                                      "other_allele", "z_value", "p_value"]].copy()
    d.columns = ["id", "chrom", "pos", "first", "second", "z", "p"]
    d["first"], d["second"] = d["first"].str.upper(), d["second"].str.upper()
    return finish(d)


def std_r(df):
    d = df.loc[df["pval"] < 0.05, ["chromosome", "position", "allele1", "allele2", "beta", "se", "pval"]].copy()
    d.columns = ["chrom", "pos", "first", "second", "beta", "se", "p"]
    d["first"], d["second"] = d["first"].str.upper(), d["second"].str.upper()
    d["z"] = d["beta"] / d["se"]
    d["id"] = ""
    return finish(d[["id", "chrom", "pos", "first", "second", "z", "p"]])


def align(x, y):
    m = x.merge(y[["key", "first", "z", "p"]], on="key", suffixes=("", "_y"))
    m["same"] = m["first"] == m["first_y"]
    m["z_y_in_x"] = np.where(m["same"], m["z_y"], -m["z_y"])
    m["discordant"] = np.sign(m["z"]) != np.sign(m["z_y_in_x"])
    return m


def summarise(s):
    n = len(s)
    out = {"n": n, "n_discordant": int(s["discordant"].sum()) if n else 0}
    if n:
        out["frac_discordant"] = out["n_discordant"] / n
        ci = stats.binomtest(out["n_discordant"], n).proportion_ci(confidence_level=0.95, method="exact")
        out["ci_lo"], out["ci_hi"] = ci.low, ci.high
        if n >= 3:
            lr = stats.linregress(s["z_y_in_x"], s["z"])
            out["slope_z_x_on_z_y"], out["r"] = lr.slope, lr.rvalue
    return out


As = std_a(A)
Rs = {k: std_r(v) for k, v in R.items()}
log(f"standardised p<0.05 rows: A={len(As)} R1={len(Rs['R1'])} R2={len(Rs['R2'])}")
pairs = {"A_vs_R1": align(As, Rs["R1"]), "A_vs_R2": align(As, Rs["R2"]), "R1_vs_R2": align(Rs["R1"], Rs["R2"])}

tiers, by_orient, by_chr, loc_lvl = [], [], [], []
for name, m in pairs.items():
    for thr in THR:
        for pal in (False, True):
            s = m[(m["p"] < thr) & (m["p_y"] < thr) & (m["pal"] == pal)]
            tiers.append({"pair": name, "p_threshold_both": thr, "palindromic": pal, **summarise(s)})
    for thr in (5e-8, 1e-4):
        base = m[(m["p"] < thr) & (m["p_y"] < thr) & (~m["pal"])]
        for same, s in base.groupby("same"):
            by_orient.append({"pair": name, "p_threshold_both": thr, "same_allele_order": bool(same), **summarise(s)})
    base = m[(m["p"] < 1e-4) & (m["p_y"] < 1e-4) & (~m["pal"])]
    for chrom, s in base.groupby("chrom"):
        by_chr.append({"pair": name, "chrom": chrom, **summarise(s)})
    for thr in (5e-8, 1e-6):
        s = m[(m["p"] < thr) & (m["p_y"] < thr) & (~m["pal"])].copy()
        if len(s):
            s["bin"] = s["chrom"].astype(str) + ":" + (s["pos"].astype("int64") // 1_000_000).astype(str)
            lead = s.loc[s.groupby("bin")["p"].idxmin()]
            loc_lvl.append({"pair": name, "p_threshold_both": thr, "unit": "lead variant per 1-Mb bin", **summarise(lead)})
    m[(m["p"] < 5e-8) & (m["p_y"] < 5e-8)][["id", "chrom", "pos", "first", "second", "z", "p", "first_y", "z_y", "p_y",
                                           "same", "z_y_in_x", "pal", "discordant"]] \
        .sort_values(["chrom", "pos"]).to_csv(OUT / f"sig_aligned_variants_5e-8_{name}.tsv", sep="\t", index=False)

pd.DataFrame(tiers).to_csv(OUT / "concordance_by_tier.tsv", sep="\t", index=False)
pd.DataFrame(by_orient).to_csv(OUT / "concordance_by_allele_order.tsv", sep="\t", index=False)
pd.DataFrame(by_chr).to_csv(OUT / "concordance_by_chromosome_p1e-4.tsv", sep="\t", index=False)
pd.DataFrame(loc_lvl).to_csv(OUT / "concordance_locus_level_1Mb.tsv", sep="\t", index=False)
log("concordance tables written")

# ---- four landmark loci in every file -------------------------------------
recs = []


def add(src, gene, rsid, chrom, pos, first, second, effect, kind, p=np.nan):
    recs.append({"source": src, "gene": gene, "rsid": rsid, "chrom": chrom, "pos": pos, "first_allele": first,
                 "second_allele": second, "effect_for_first": effect, "effect_kind": kind, "p": p,
                 "implied_raising_allele": (first if effect > 0 else second) if np.isfinite(effect) else None})


for gene, rsid in LOCI:
    for _, a in A.loc[A["variant_id"] == rsid].iterrows():
        chrom, pos = a["chromosome"], int(a["base_pair_location"])
        aset = {a["effect_allele"].upper(), a["other_allele"].upper()}
        add("A", gene, rsid, chrom, pos, a["effect_allele"], a["other_allele"], a["beta"], "beta", a["p_value"])
        for _, b in B.loc[B["variant_id"] == rsid].iterrows():
            add("B", gene, rsid, chrom, pos, b["effect_allele"], b["other_allele"], b["beta"], "beta", b["p_value"])
        for _, c in C.loc[(C["chromosome"] == chrom) & (C["position"] == pos)].iterrows():
            if {c["allele1"].upper(), c["allele2"].upper()} == aset:
                add("C", gene, rsid, chrom, pos, c["allele1"], c["allele2"], c["beta"], "beta", c["pval"])
        for _, d in D.loc[D["SNP"] == rsid].iterrows():
            if isinstance(d["A1"], str):
                add("D", gene, rsid, chrom, pos, d["A1"], d["A2"], d["Z"], "Z", np.nan)
        for k, rdf in R.items():
            for _, r in rdf.loc[(rdf["chromosome"] == chrom) & (rdf["position"] == pos)].iterrows():
                if {r["allele1"].upper(), r["allele2"].upper()} == aset:
                    add(k, gene, rsid, chrom, pos, r["allele1"], r["allele2"], r["beta"], "beta", r["pval"])
loci = pd.DataFrame(recs)
loci.to_csv(OUT / "four_loci_long.tsv", sep="\t", index=False)
wide = loci.pivot_table(index=["gene", "rsid", "chrom", "pos"], columns="source", values="implied_raising_allele",
                        aggfunc="first").reset_index()
wide.to_csv(OUT / "four_loci_raising_allele.tsv", sep="\t", index=False)
log("four-locus tables written")

# ---- derived files vs A (key based: B drops rows and C is not in A's row order) ------
cons = []


def add_cons(f, check, value):
    cons.append({"file": f, "check": check, "value": value})


def keyed(df, chrom, pos, a1, a2):
    u1, u2 = df[a1].str.upper(), df[a2].str.upper()
    lo, hi = u1.where(u1 < u2, u2), u2.where(u1 < u2, u1)
    return df[chrom].astype(str) + ":" + df[pos].astype("int64").astype(str) + ":" + lo + ":" + hi


ka = keyed(A, "chromosome", "base_pair_location", "effect_allele", "other_allele")
a_ok = ~ka.duplicated(keep=False)
a_small = pd.DataFrame({"key": ka[a_ok].values, "ea": A.loc[a_ok, "effect_allele"].str.upper().values,
                        "beta_a": A.loc[a_ok, "beta"].values})
add_cons("A", "rows", len(A))
add_cons("A", "rows with a unique chromosome:position:allele-pair key", int(a_ok.sum()))
# B: every B row should be an A row with identical first ten columns
bu = B[~B["variant_id"].duplicated(keep=False)]
au = A[~A["variant_id"].duplicated(keep=False)]
mb = bu.merge(au, on="variant_id", suffixes=("_B", "_A"))
add_cons("B", "rows", len(B))
add_cons("B", "rows with a unique variant_id", len(bu))
add_cons("B", "rows matched to A by variant_id", len(mb))
add_cons("B", "A rows whose variant_id is absent from B", int((~A["variant_id"].isin(B["variant_id"])).sum()))
add_cons("B", "B rows whose variant_id is absent from A", int((~B["variant_id"].isin(A["variant_id"])).sum()))
for c in [c for c in A.columns if c != "variant_id"]:
    x, y = mb[c + "_B"].values, mb[c + "_A"].values
    same = np.array_equal(x, y, equal_nan=True) if A[c].dtype.kind in "fc" else bool((x == y).all())
    add_cons("B", f"column {c} identical to A on matched rows", bool(same))
# C: allele1 should be A's effect allele and beta A's beta
kc = keyed(C, "chromosome", "position", "allele1", "allele2")
c_ok = ~kc.duplicated(keep=False)
c_small = pd.DataFrame({"key": kc[c_ok].values, "a1": C.loc[c_ok, "allele1"].str.upper().values,
                        "beta_c": C.loc[c_ok, "beta"].values})
mc = a_small.merge(c_small, on="key")
add_cons("C", "rows", len(C))
add_cons("C", "rows with a unique key", int(c_ok.sum()))
add_cons("C", "unique-key rows matched to A", len(mc))
add_cons("C", "same chromosome and position order as A, row by row",
         bool(len(C) == len(A) and (C["position"].values == A["base_pair_location"].values).all()
              and (C["chromosome"].values == A["chromosome"].values).all()))
add_cons("C", "fraction of matched rows with allele1 == A effect_allele", float((mc["a1"] == mc["ea"]).mean()))
add_cons("C", "fraction of matched rows with beta equal to A beta (rtol 1e-9)",
         float(np.isclose(mc["beta_c"], mc["beta_a"], rtol=1e-9).mean()))
# D: A1 should be A's effect allele and Z the sign of A's z
am = au[["variant_id", "effect_allele", "other_allele", "z_value"]]
dm = D.dropna(subset=["A1"]).merge(am, left_on="SNP", right_on="variant_id")
add_cons("D", "rows", len(D))
add_cons("D", "rows with empty A1/A2/Z/N", int(D["A1"].isna().sum()))
add_cons("D", "rows joined to A by rsid", len(dm))
add_cons("D", "fraction A1 == A effect_allele", float((dm["A1"].str.upper() == dm["effect_allele"].str.upper()).mean()))
add_cons("D", "fraction A2 == A other_allele", float((dm["A2"].str.upper() == dm["other_allele"].str.upper()).mean()))
nz = dm[(dm["Z"] != 0) & (dm["z_value"] != 0)]
add_cons("D", "fraction sign(Z) == sign(A z_value)", float((np.sign(nz["Z"]) == np.sign(nz["z_value"])).mean()))
cons_df = pd.DataFrame(cons)
cons_df.to_csv(OUT / "derived_file_consistency.tsv", sep="\t", index=False)
log("derived-file checks written")

# ---- hg19 columns of B (report only) ---------------------------------------
b = B.copy()
b["offset"] = pd.to_numeric(b["hg19_pos"], errors="coerce") - b["base_pair_location"]
hg = {"rows": len(b), "hg19_pos_missing": int(b["offset"].isna().sum()), "offset_zero": int((b["offset"] == 0).sum()),
      "offset_nonzero": int((b["offset"].fillna(0) != 0).sum()),
      "hg19_chr_equals_chromosome": float((b["hg19_chr"] == b["chromosome"]).mean())}
(OUT / "hg19_summary.json").write_text(json.dumps(hg, indent=2))
b.groupby("chromosome")["offset"].agg(["count", "median", "min", "max"]).reset_index() \
    .to_csv(OUT / "hg19_offsets_by_chr.tsv", sep="\t", index=False)
ref37 = {}
for chrom in ("2", "19", "22"):
    t = pd.read_csv(G / f"ldsc_ref/eur_w_ld_chr/{chrom}.l2.ldscore.gz", sep="\t", usecols=["SNP", "BP"])
    ref37.update(dict(zip(t["SNP"], t["BP"])))
rows = []
for gene, rsid in LOCI:
    a = A.loc[A["variant_id"] == rsid].iloc[0]
    bb = B.loc[B["variant_id"] == rsid].iloc[0]
    n_ref = {k: int(((v["chromosome"] == a["chromosome"]) & (v["position"] == a["base_pair_location"])).sum()) for k, v in R.items()}
    rows.append({"gene": gene, "rsid": rsid, "A_base_pair_location": int(a["base_pair_location"]),
                 "B_hg19_pos": float(pd.to_numeric(bb["hg19_pos"])),
                 "B_offset": float(pd.to_numeric(bb["hg19_pos"])) - float(a["base_pair_location"]),
                 "ldsc_1000G_GRCh37_BP": ref37.get(rsid, np.nan), "rows_in_R1_at_A_pos": n_ref["R1"],
                 "rows_in_R2_at_A_pos": n_ref["R2"]})
pd.DataFrame(rows).to_csv(OUT / "hg19_four_loci.tsv", sep="\t", index=False)
log("hg19 checks written")

# ---- decision ---------------------------------------------------------------
t = pd.DataFrame(tiers)


def frac(pair, thr=5e-8, pal=False):
    r = t[(t["pair"] == pair) & (t["p_threshold_both"] == thr) & (t["palindromic"] == pal)].iloc[0]
    return float(r.get("frac_discordant", np.nan)), int(r["n"])


def cval(file, check):
    v = cons_df[(cons_df["file"] == file) & (cons_df["check"] == check)]["value"]
    return v.iloc[0] if len(v) else None


f1, n1 = frac("A_vs_R1")
f2, n2 = frac("A_vs_R2")
f12, n12 = frac("R1_vs_R2")
ref_raising = loci[loci["source"].isin(["R1", "R2"])].groupby(["gene", "source"])["implied_raising_allele"].first().unstack()
refs_agree_loci = bool(all(ref_raising.loc[g, "R1"] == ref_raising.loc[g, "R2"] for g, _ in LOCI))
src_raising = {s: loci[loci["source"] == s].groupby("gene")["implied_raising_allele"].first() for s in ("A", "B", "C", "D")}
loci_match_refs = {s: {g: bool(v[g] == ref_raising.loc[g, "R1"]) for g in v.index} for s, v in src_raising.items()}
b_col_checks = [bool(v) for v in cons_df[(cons_df["file"] == "B") & cons_df["check"].str.startswith("column ")]["value"]]
c_match_frac = cval("C", "unique-key rows matched to A") / max(cval("C", "rows with a unique key"), 1)
labels_consistent = bool(
    all(b_col_checks) and (cval("C", "fraction of matched rows with allele1 == A effect_allele") or 0) >= 0.999
    and c_match_frac >= 0.99 and (cval("D", "fraction A1 == A effect_allele") or 0) >= 0.999
    and (cval("D", "fraction sign(Z) == sign(A z_value)") or 0) >= 0.999)
decision = {
    "mode": MODE,
    "frac_discordant_A_vs_R1_5e-8_nonpalindromic": f1, "n_A_vs_R1": n1,
    "frac_discordant_A_vs_R2_5e-8_nonpalindromic": f2, "n_A_vs_R2": n2,
    "frac_discordant_R1_vs_R2_5e-8_nonpalindromic": f12, "n_R1_vs_R2": n12,
    "references_agree": bool(f12 <= 0.05 and refs_agree_loci),
    "derived_copies_B_C_D_share_A_labels": labels_consistent,
    "C_matched_fraction_of_unique_keys": c_match_frac,
    "four_loci_source_matches_reference_raising_allele": loci_match_refs,
    "eaf_like_columns_present": {k: v for k, v in eaf_cols.items() if v},
    "columns": cols,
    "row_counts": {"A": len(A), "B": len(B), "C": len(C), "D": len(D), "R1": len(R["R1"]), "R2": len(R["R2"])},
}
if MODE == "original":
    decision["reversal_global"] = bool(f1 >= 0.95 and f2 >= 0.95)
    decision["four_loci_all_reversed"] = bool(all(not v for s in ("A", "B", "C") for v in loci_match_refs[s].values())
                                              and all(not v for v in loci_match_refs["D"].values()))
    proceed = (decision["reversal_global"] and decision["references_agree"] and decision["four_loci_all_reversed"]
               and decision["derived_copies_B_C_D_share_A_labels"])
    decision["fix"] = "swap_effect_other_labels_keep_beta_z_se_p" if proceed else "STOP_no_fix_chosen"
    decision["fix_note"] = ("No allele-frequency column exists in A-D, so swapping the labels and negating beta/z are the "
                            "same statement about the data; the label swap is chosen because the effect allele then "
                            "equals the PDFF-raising allele at the four loci, as in both references.")
else:
    decision["corrected_ok"] = bool(
        f1 <= 0.01 and f2 <= 0.01 and decision["references_agree"] and decision["derived_copies_B_C_D_share_A_labels"]
        and all(v for s in ("A", "B", "C", "D") for v in loci_match_refs[s].values()))
(OUT / "decision.json").write_text(json.dumps(decision, indent=2, default=str))
(OUT / "pip_freeze.txt").write_text(subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True).stdout)
(OUT / "environment.txt").write_text(sys.version + "\n")
log(json.dumps({k: decision[k] for k in decision if k not in ("columns",)}, indent=2, default=str))
log("done")
