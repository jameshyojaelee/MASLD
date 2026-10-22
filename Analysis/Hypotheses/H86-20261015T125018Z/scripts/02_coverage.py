#!/usr/bin/env python3
"""H86 step 2 (prespec PRESPEC_H86.md, sha256 7e96f97f...).

Joins the PXD051911 adjusted liver protein DE table to Gene Catalog r8 by gene symbol, applies the
prespecified symbol-resolution order, runs gates G1-G4, controls PC-1..PC-3, NC-1, NC-2, diagnostics
D1-D9 and stop rules S1-S8 in order, and writes T1, T1b and RESULT_H86.md with the verdict computed here.
There is no test and no new p value. No step is random; the seed is set because the prespec asks for it.

--debug-permute-classes: for code checks only. Gates and controls run on the real Catalog (their values
were known or disclosed before writing); then evidence classes are permuted across Catalog genes before
T1, T1b and D7 are computed, and D5 is computed on a random 300-symbol subset. Outputs are prefixed DEBUG_.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
H86 = Path(__file__).resolve().parent.parent
PRESPEC_SHA = "7e96f97f0185657daa34eb561898e44358d1c03a8a2a9e6bb399e9678b2ebbc1"
R4 = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z"

CLASSES = ["concordant", "established_state_associated", "genetically_anchored", "unresolved", "untested"]
CLASSED = ["concordant", "established_state_associated", "genetically_anchored"]
EXPECTED_CLASS_N = {"untested": 13551, "unresolved": 13127, "established_state_associated": 1129,
                    "genetically_anchored": 400, "concordant": 26}
PC1_EXPECTED = {"concordant": (26, 11, 6), "established_state_associated": (1129, 170, 40),
                "genetically_anchored": (400, 199, 9), "unresolved": (13127, 5163, 218),
                "untested": (13551, 409, 10)}
PC2_GENES = ["ALB", "APOA1", "APOB", "APOE", "CYP2E1", "CYP3A4", "FASN", "HSD17B13", "GCKR", "PLIN2",
             "COL1A1", "ACTA2", "LUM", "MTARC1", "PNPLA3", "TM6SF2"]
PC3_PAIRS = [("MARC1", "MTARC1"), ("MARC2", "MTARC2"), ("AARS", "AARS1"), ("CD97", "ADGRE5"),
             ("CERT", "CERT1"), ("ATP5MD", "ATP5MK")]
NC1_GENES = ["INS", "GCG", "NEUROD1", "SYP", "OLIG2", "TNNT2", "CD8A", "MYH7", "TTN", "KRT1"]
SIG_REMAP_ALLOWED = {"ATP5MD", "H2AFY", "H2AFY2", "HIST1H1B", "LARS", "MARC2", "METTL7A"}
D4_EXPECTED = {"n_unmatched": 119, "single": 114, "ambiguous": 5, "none": 0, "single_in_catalog": 90,
               "target_has_direct_row": 0, "agree": 109, "both_unique": 113}
D4_AMBIGUOUS_EXPECTED = ["LOR", "MPP6", "QARS", "SARS", "SMAP"]
D9_EXPECTED = (117, 283, ["PROGRAM_MEMBERSHIP_UNMAPPED_SOURCE"])
EXTRA_ROWS = ["identity quarantine (no Catalog gene)", "not in Catalog", "ambiguous or colliding symbol"]
DIRECT_ONLY_ROW = "no direct Catalog symbol (outdated protein symbols not remapped)"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bh(p):
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


class Run:
    def __init__(self, out, catalog, seed, debug):
        self.out = Path(out)
        self.catalog = catalog
        self.seed = seed
        self.debug = debug
        self.prefix = "DEBUG_" if debug else ""
        self.lines = []
        self.gates = []
        self.result = {"prespec_sha256": PRESPEC_SHA, "debug_permuted_classes": debug, "seed": seed}

    def say(self, *a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        self.lines.append(s)

    def gate(self, gid, what, predicted, observed, ok, rule=""):
        self.gates.append({"id": gid, "check": what, "predicted": str(predicted), "observed": str(observed),
                           "rule": rule, "pass": bool(ok)})
        self.say(f"[{gid}] {'PASS' if ok else 'FAIL'} | {what} | predicted {predicted} | observed {observed}")
        return bool(ok)

    def w(self, df, name):
        df.to_csv(self.out / f"{self.prefix}{name}", sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--debug-permute-classes", action="store_true")
    a = ap.parse_args()
    np.random.seed(a.seed)
    rng = np.random.default_rng(a.seed)
    R = Run(a.out, a.catalog, a.seed, a.debug_permute_classes)
    R.say(f"[02] H86 coverage; out={R.out}; catalog={a.catalog}; seed={a.seed}; debug={R.debug}")

    # ---------------- S1 / D1: hashes --------------------------------------------------------------
    cat_dir = ROOT / a.catalog
    expected = [
        ("prespec_copy", H86 / "PRESPEC_H86.md", None, PRESPEC_SHA),
        ("prespec_original", ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H86.md", None, PRESPEC_SHA),
        ("catalog_gene_index", cat_dir / "passport_gene_index.parquet", 468456, "ffd3c4650c77a899e9ebed7327697398c4b6462cfffcd7f0e635cd0d023db06b"),
        ("catalog_quarantine", cat_dir / "passport_identity_quarantine.tsv", 243387, "f696105dee84a5837f080feadf8dc8f44164b2febaafce3f676883c758e50e6d"),
        ("catalog_domain_coverage", cat_dir / "passport_domain_coverage.tsv", 3868, "d5c541b71700f8ec2365c0ccb9b9875ffbf934cf7f8c3e24fa2416e75e3f194e"),
        ("protein_de", ROOT / "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv", 1361989, "4cab14f405049964a441d45a298e6209f5067fb48ec233873b3c839800cb6df3"),
        ("deposit_liver", ROOT / "data/PXD051911/liver_protein_quant.txt", 5599042, "c4fafdfb2339d54d17a91779098012e08448aa5ebcc8a865988dcc49f92c969e"),
        ("deposit_meta", ROOT / "data/PXD051911/meta_data.txt", 35867, "340874d4f720f909fcfca428fed6da16ba3f7cfdf4d9fd53dd6e0378198ba2bc"),
        ("projection_script", ROOT / "Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R", 29975, "20a3c06cfa4e1a726339cedc4b701ab0d9d276585791d1a3c7202045dc19dd6b"),
        ("gate_cov_py", ROOT / R4 / "gate/R11-09/cov.py", 1059, "ce6af28bf2431a8fe13a388d1337ec566cb8df8f875278847f684b5d6924e0f0"),
        ("writer_remap_R", ROOT / R4 / "write/H86/remap_feasibility.R", 2092, "d1332e26541e4c779eae34c881ae79f3da7596dd056e62b01abcd6fff7c4e528"),
    ]
    rows = []
    for role, path, size, sha in expected:
        ex = path.exists()
        obs = sha256(path) if ex else None
        sz = path.stat().st_size if ex else None
        ok = ex and obs == sha and (size is None or sz == size)
        rows.append({"role": role, "path": str(path), "size": sz, "sha256_expected": sha, "sha256_observed": obs, "match": ok})
    hashes = pd.DataFrame(rows)
    R.w(hashes, "D1_inputs_sha256.tsv")
    step01 = (R.out / "STEP01_STATUS.txt").read_text().strip() if (R.out / "STEP01_STATUS.txt").exists() else "MISSING"
    ok_s1 = R.gate("S1/D1", "every input sha256 and size match the prespec; step 01 status OK",
                   "11 of 11 match; OK", f"{int(hashes.match.sum())} of {len(hashes)} match; step01 {step01}",
                   bool(hashes.match.all()) and step01 == "OK", "any mismatch stops")
    if not ok_s1:
        bad = hashes.loc[~hashes.match, "role"].tolist()
        return finish(R, "S1", f"input hash mismatch: {bad}; step01={step01}")

    versions = pd.read_csv(R.out / "versions.tsv", sep="\t", keep_default_na=False)
    vd = dict(zip(versions["item"], versions["value"]))
    R.result["versions"] = vd
    R.say("[D1] versions:", json.dumps(vd))

    # ---------------- G1 / D2 / D3 -----------------------------------------------------------------
    cat = pd.read_parquet(cat_dir / "passport_gene_index.parquet",
                          columns=["ensembl_id", "symbol", "symbol_collision", "primary_evidence_class"])
    cls_n = cat.primary_evidence_class.value_counts(dropna=False).to_dict()
    d3_ok = (len(cat) == 28233 and cat.symbol.nunique() == 28233 and cat.ensembl_id.nunique() == 28233
             and int(cat.symbol_collision.sum()) == 0 and cls_n == EXPECTED_CLASS_N)
    R.gate("G1/D3", "Catalog rows, unique symbols, unique Ensembl IDs, collisions, class totals",
           f"28233/28233/28233/0 {EXPECTED_CLASS_N}",
           f"{len(cat)}/{cat.symbol.nunique()}/{cat.ensembl_id.nunique()}/{int(cat.symbol_collision.sum())} {cls_n}", d3_ok)

    de = pd.read_csv(ROOT / "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv",
                     sep="\t", keep_default_na=False, na_values=[""])
    n_sep = int(de.gene.str.contains(r"[;,|]").sum())
    bh_diff = float(np.max(np.abs(bh(de.pvalue.values) - de.padj.values)))
    n_sig = int((de.padj < 0.05).sum())
    d2_ok = len(de) == 6071 and de.gene.nunique() == 6071 and n_sep == 0 and bh_diff < 1e-12 and n_sig == 290
    R.gate("G1/D2", "DE table rows / unique / separators / BH recompute max diff / padj<0.05",
           "6071/6071/0/<1e-12/290", f"{len(de)}/{de.gene.nunique()}/{n_sep}/{bh_diff:.3g}/{n_sig}", d2_ok)
    R.result["bh_max_diff"] = bh_diff

    # sample counts (descriptive, Unit of replication section)
    dep_cols = pd.read_csv(ROOT / "data/PXD051911/liver_protein_quant.txt", sep="\t", nrows=0).columns.tolist()
    sample_cols = [c for c in dep_cols if c not in ("ProteinAccessions", "Genes", "ProteinDescriptions")]
    meta = pd.read_csv(ROOT / "data/PXD051911/meta_data.txt", sep="\t")
    meta = meta[meta.liver_proteomics_filename.notna() & (meta.liver_proteomics_filename != "")
                & meta.liver_proteomics_filename.isin(sample_cols)].drop_duplicates("liver_proteomics_filename")
    meta["grp"] = np.where(meta.saf_diagnosis == "No_MASLD", "control", "MASLD")
    meta["batch"] = np.where(meta.liver_proteomics_filename.str.startswith("2019"), "2019", "2020")
    sample_tab = meta.groupby(["batch", "grp"]).size().to_dict()
    n_missing_cov = int(meta[["alder", "bmi", "gender"]].isna().sum().sum())
    s_obs = (f"{len(meta)} samples, {meta.patient_name.nunique()} patients, "
             f"{int((meta.grp == 'control').sum())} control / {int((meta.grp == 'MASLD').sum())} MASLD, "
             f"2019 {sample_tab.get(('2019', 'control'), 0)}/{sample_tab.get(('2019', 'MASLD'), 0)}, "
             f"2020 {sample_tab.get(('2020', 'control'), 0)}/{sample_tab.get(('2020', 'MASLD'), 0)}, "
             f"missing age/BMI/sex {n_missing_cov}")
    R.gate("Unit", "PXD051911 liver samples (descriptive; not a stop rule)",
           "58 samples, 58 patients, 12 control / 46 MASLD, 2019 5/28, 2020 7/18, missing age/BMI/sex 0", s_obs,
           s_obs == "58 samples, 58 patients, 12 control / 46 MASLD, 2019 5/28, 2020 7/18, missing age/BMI/sex 0",
           "descriptive")
    R.result["samples"] = s_obs
    if not (d3_ok and d2_ok):
        return finish(R, "G1", "Catalog or DE-table counts differ from the prespec although hashes match")

    cat_sym = set(cat.symbol)
    de_sym = set(de.gene)
    sym2class = dict(zip(cat.symbol, cat.primary_evidence_class))
    sym2ens = dict(zip(cat.symbol, cat.ensembl_id))
    de_idx = de.set_index("gene")

    # ---------------- PC-1 / S2: direct join reproduction ------------------------------------------
    sig_sym = set(de.gene[de.padj < 0.05])
    obs_pc1 = {}
    for c in CLASSES:
        s = cat.symbol[cat.primary_evidence_class == c]
        obs_pc1[c] = (len(s), int(s.isin(de_sym).sum()), int(s.isin(sig_sym).sum()))
    tot_q = sum(v[1] for v in obs_pc1.values())
    tot_s = sum(v[2] for v in obs_pc1.values())
    ok_pc1 = R.gate("PC-1/G2", "direct-join triples total/quantified/significant per class; totals",
                    f"{PC1_EXPECTED}; 5952 / 283", f"{obs_pc1}; {tot_q} / {tot_s}",
                    obs_pc1 == PC1_EXPECTED and tot_q == 5952 and tot_s == 283, "any difference stops (S2)")
    R.result["pc1"] = {k: list(v) for k, v in obs_pc1.items()}
    if not ok_pc1:
        return finish(R, "S2", "direct-join reproduction failed")

    # ---------------- D4 / S7: remap of the unmatched symbols ---------------------------------------
    um = pd.read_csv(R.out / "unmatched_119.tsv", sep="\t", keep_default_na=False, na_values=[""])
    hall = pd.read_csv(R.out / "hgnc_all_symbols.tsv", sep="\t", keep_default_na=False, na_values=[""])
    hall_idx = hall.set_index("x")
    um_py = sorted(de_sym - cat_sym)
    consistent = sorted(um.x) == um_py and all(
        (pd.isna(hall_idx.at[x, "Suggested.Symbol"]) and pd.isna(s)) or hall_idx.at[x, "Suggested.Symbol"] == s
        for x, s in zip(um.x, um["Suggested.Symbol"]))
    if not consistent:
        raise RuntimeError("step 01 unmatched set or HGNChelper calls disagree with step 02; code error")
    um["single"] = um.n_sugg == 1
    um["single_in_catalog"] = um.single & um["Suggested.Symbol"].isin(cat_sym)
    um["target_has_direct_row"] = um.single_in_catalog & um["Suggested.Symbol"].isin(de_sym)
    both = um[um.single & um.org_symbol.notna()]
    agree = int((both["Suggested.Symbol"] == both.org_symbol).sum())
    disagree = both[both["Suggested.Symbol"] != both.org_symbol][["x", "Suggested.Symbol", "org_symbol"]]
    amb = sorted(um.x[um.n_sugg > 1])
    tgt = um.loc[um.single_in_catalog, "Suggested.Symbol"]
    dup_targets = sorted(tgt[tgt.duplicated(keep=False)].unique())
    d4 = {"n_unmatched": len(um), "single": int(um.single.sum()), "ambiguous": int((um.n_sugg > 1).sum()),
          "none": int((um.n_sugg == 0).sum()), "single_in_catalog": int(um.single_in_catalog.sum()),
          "target_has_direct_row": int(um.target_has_direct_row.sum()), "agree": agree, "both_unique": len(both)}
    R.result["d4"] = d4
    R.result["d4_ambiguous"] = amb
    R.result["d4_disagreements"] = disagree.to_dict("records")
    R.result["d4_duplicate_targets"] = dup_targets
    R.say("[D4] ambiguous:", amb, "| HGNChelper vs org.Hs.eg.db disagreements (HGNChelper used):",
          disagree.to_dict("records"), "| duplicate remap targets:", dup_targets)
    um.to_csv(R.out / f"{R.prefix}D4_unmatched_119_resolution.tsv", sep="\t", index=False)
    ver_ok = (vd.get("HGNChelper") == "0.8.15" and "Nov 16" in vd.get("HGNChelper_map_date", "")
              and "2024" in vd.get("HGNChelper_map_date", "") and vd.get("org.Hs.eg.db") == "3.20.0")
    ok_s7 = R.gate("D4/S7", "package versions, map date, remap counts single/ambiguous/none, in Catalog, "
                   "targets with direct row, HGNChelper-org.Hs.eg.db agreement",
                   f"HGNChelper 0.8.15 map 2024-11-16, org.Hs.eg.db 3.20.0; {D4_EXPECTED}; ambiguous {D4_AMBIGUOUS_EXPECTED}",
                   f"HGNChelper {vd.get('HGNChelper')} map {vd.get('HGNChelper_map_date')}, org.Hs.eg.db {vd.get('org.Hs.eg.db')}; {d4}; ambiguous {amb}",
                   ver_ok and d4 == D4_EXPECTED, "any difference stops before T1 (S7)")
    if not ok_s7:
        return finish(R, "S7", "package version, map date or D4 remap counts differ from the prespec")

    # ---------------- D9 / S8 ----------------------------------------------------------------------
    qt = pd.read_csv(cat_dir / "passport_identity_quarantine.tsv", sep="\t", keep_default_na=False)
    q_src = set(qt.source_symbol) - {""}
    ov = de_sym & q_src & cat_sym
    ov_rows = qt[qt.source_symbol.isin(ov)]
    d9 = (len(ov), len(ov_rows), sorted(ov_rows.reason_code.unique().tolist()))
    ok_s8 = R.gate("D9/S8", "direct protein symbols that are also quarantine source_symbol: symbols, rows, reason codes",
                   D9_EXPECTED, d9, d9 == D9_EXPECTED, "any difference stops before T1 (S8)")
    R.result["d9"] = list(d9)
    if not ok_s8:
        return finish(R, "S8", "D9 differs from 117 symbols / 283 rows / one reason code")

    # ---------------- G3 / S3 ----------------------------------------------------------------------
    ok_g3 = R.gate("G3", "unmatched symbols with a single HGNChelper suggestion", ">=100 of 119 (writer 114)",
                   f"{d4['single']} of {d4['n_unmatched']}", d4["single"] >= 100, "below 100: direct join only (S3)")

    # resolution category for each unmatched DE symbol (Definitions steps 2-5)
    def category(r):
        if r.n_sugg > 1:
            return "ambiguous or colliding symbol"
        if r.single_in_catalog:
            if r.target_has_direct_row or r["Suggested.Symbol"] in dup_targets:
                return "ambiguous or colliding symbol"
            return "remapped"
        sug = r["Suggested.Symbol"] if r.single else None
        if r.x in q_src or (sug is not None and sug in q_src):
            return "identity quarantine (no Catalog gene)"
        return "not in Catalog"
    um["category"] = um.apply(category, axis=1)
    remap = um[um.category == "remapped"][["x", "Suggested.Symbol"]].rename(columns={"x": "old", "Suggested.Symbol": "target"})

    # ---------------- PC-3 -------------------------------------------------------------------------
    rm_map = dict(zip(remap.old, remap.target))
    pc3_rows = []
    for old, new in PC3_PAIRS:
        r = um[um.x == old]
        single_ok = len(r) == 1 and bool(r.single.iloc[0]) and r["Suggested.Symbol"].iloc[0] == new and new in cat_sym
        pc3_rows.append({"old": old, "target": new, "single_and_catalog": single_ok, "remapped": rm_map.get(old) == new,
                         "class_r8": sym2class.get(new), "padj_old": float(de_idx.at[old, "padj"]) if old in de_idx.index else None})
    pc3 = pd.DataFrame(pc3_rows)

    def cls_ok(new, want_class, want_sig):
        r = pc3[pc3.target == new].iloc[0]
        return bool(r.remapped and r.class_r8 == want_class and ((r.padj_old < 0.05) == want_sig))
    pc3_class_ok = (cls_ok("MTARC1", "unresolved", False) and cls_ok("MTARC2", "genetically_anchored", True)
                    and cls_ok("ATP5MK", "unresolved", True) and cls_ok("AARS1", "unresolved", False)
                    and cls_ok("ADGRE5", "unresolved", False) and cls_ok("CERT1", "unresolved", False))
    R.w(pc3, "PC3_remap_known_answers.tsv")
    ok_pc3 = R.gate("PC-3", "six remap known answers single + Catalog; MTARC1 quantified unresolved; MTARC2 sig GA; "
                    "ATP5MK sig unresolved; AARS1/ADGRE5/CERT1 quantified not sig unresolved",
                    "all 6 hold, classes as written", pc3.to_dict("records"),
                    bool(pc3.single_and_catalog.all()) and bool(pc3.remapped.all()) and pc3_class_ok,
                    "all must hold, else remap not used (S4)")

    # ---------------- NC-2 -------------------------------------------------------------------------
    ens = pd.read_csv(R.out / "nc2_alias_ensembl.tsv", sep="\t", keep_default_na=False, na_values=[""])
    ens = ens[ens.ENSEMBL.notna()]
    old2ens = ens.groupby("ALIAS").ENSEMBL.apply(set).to_dict()
    pairs = remap.sort_values("old").reset_index(drop=True)
    pairs["target_ensembl"] = pairs.target.map(sym2ens)
    pairs["shift_target"] = np.roll(pairs.target.values, -1)
    pairs["shift_target_ensembl"] = pairs.shift_target.map(sym2ens)
    pairs["has_ids"] = pairs.old.map(lambda s: s in old2ens)
    pairs["agree"] = [t in old2ens.get(o, set()) for o, t in zip(pairs.old, pairs.target_ensembl)]
    pairs["agree_shift"] = [t in old2ens.get(o, set()) for o, t in zip(pairs.old, pairs.shift_target_ensembl)]
    n_ids = int(pairs.has_ids.sum())
    n_agree = int(pairs.agree[pairs.has_ids].sum())
    n_shift = int(pairs.agree_shift[pairs.has_ids].sum())
    R.w(pairs, "NC2_identifier_check.tsv")
    ok_nc2 = R.gate("NC-2", "org.Hs.eg.db ALIAS->ENSEMBL of old symbol contains Catalog ensembl_id of HGNChelper target; "
                    "shifted pairing (Python code-point sort of old symbols)",
                    "89 of 90 with IDs; >=85 of 89 agree [Estimate]; shifted <10",
                    f"{n_ids} of {len(pairs)} with IDs; {n_agree} agree; shifted {n_shift}",
                    n_agree >= 80 and n_shift < 10, "fewer than 80 agree, or 10+ shifted agreements: remap not used (S4)")
    R.result["nc2"] = {"n_pairs": len(pairs), "n_with_ids": n_ids, "agree": n_agree, "agree_shifted": n_shift}

    remap_used = True
    remap_off = None
    if not ok_g3:
        remap_used, remap_off = False, "S3"
    elif not (ok_pc3 and ok_nc2):
        remap_used, remap_off = False, "S4"
    R.result["remap_used"] = remap_used
    R.say(f"[remap] used={remap_used} reason_off={remap_off}")

    # ---------------- assignment of every protein row ----------------------------------------------
    assign = pd.DataFrame({"protein_symbol": de.gene, "logFC": de.logFC, "padj": de.padj})
    assign["catalog_symbol"] = np.where(assign.protein_symbol.isin(cat_sym), assign.protein_symbol, None)
    assign["route"] = np.where(assign.catalog_symbol.notna(), "direct", None)
    um_cat = dict(zip(um.x, um.category))
    for i, s in assign.protein_symbol.items():
        if assign.at[i, "route"] == "direct":
            continue
        if remap_used:
            c = um_cat[s]
            if c == "remapped":
                assign.at[i, "catalog_symbol"] = rm_map[s]
                assign.at[i, "route"] = "remapped"
            else:
                assign.at[i, "route"] = c
        else:
            assign.at[i, "route"] = DIRECT_ONLY_ROW
    assign["sig"] = assign.padj < 0.05

    # ---------------- PC-2 / NC-1 / S5 -------------------------------------------------------------
    quant_sym = set(assign.catalog_symbol.dropna())
    pc2 = {g: g in quant_sym for g in PC2_GENES}
    nc1 = {g: g in quant_sym for g in NC1_GENES}
    n_pc2, n_nc1 = sum(pc2.values()), sum(nc1.values())
    pc2_pred_match = (not pc2["PNPLA3"]) and (not pc2["TM6SF2"]) and \
        (assign.loc[assign.catalog_symbol == "MTARC1", "route"].tolist() == (["remapped"] if remap_used else []))
    ok_pc2 = R.gate("PC-2", f"known liver proteins quantified ({'after remap' if remap_used else 'direct join'})",
                     "14 of 16; PNPLA3, TM6SF2 absent; MTARC1 only via MARC1", f"{n_pc2} of 16 {pc2}; detail prediction met: {pc2_pred_match}",
                     n_pc2 >= 12, "below 12 stops (S5)")
    ok_nc1 = R.gate("NC-1", "tissue-restricted non-liver proteins quantified",
                    "MYH7, TTN, KRT1 quantified; other 7 absent (3 of 10)", f"{n_nc1} of 10 {nc1}; "
                    f"all 10 Catalog symbols: {all(g in cat_sym for g in NC1_GENES)}",
                    n_nc1 <= 3, "more than 3 stops (S5)")
    R.result["pc2"] = {"n": n_pc2, "detail": pc2, "detail_prediction_met": pc2_pred_match}
    R.result["nc1"] = {"n": n_nc1, "detail": nc1}
    if not (ok_pc2 and ok_nc1):
        return finish(R, "S5", "known-answer controls: the protein table does not behave like a bulk liver proteome")

    # ---------------- outcome stage (blind in the final run) -------------------------------------
    if R.debug:
        cat = cat.copy()
        cat["primary_evidence_class"] = rng.permutation(cat.primary_evidence_class.values)
        sym2class = dict(zip(cat.symbol, cat.primary_evidence_class))
        R.say("[DEBUG] evidence classes permuted across Catalog genes; T1/T1b/D7 below are not real")

    assign["class"] = assign.catalog_symbol.map(sym2class)
    R.w(assign, "protein_row_assignment.tsv")

    # ---------------- D8 / S6 conservation ----------------------------------------------------------
    per_gene = assign[assign.catalog_symbol.notna()]
    dup_gene = int(per_gene.catalog_symbol.duplicated().sum())
    sig_added = set(assign.protein_symbol[(assign.route == "remapped") & assign.sig])

    # ---------------- T1 ---------------------------------------------------------------------------
    def row(label, n_genes, sub):
        nq = len(sub)
        ns = int(sub.sig.sum())
        return {"row": label, "n_genes": n_genes, "n_quantified": nq,
                "pct_quantified": round(100 * nq / n_genes, 2) if n_genes else None,
                "n_protein_sig": ns, "n_sig_up": int((sub.sig & (sub.logFC > 0)).sum()),
                "n_sig_down": int((sub.sig & (sub.logFC < 0)).sum()),
                "pct_sig_of_quantified": round(100 * ns / nq, 2) if nq else None,
                "n_quantified_direct": int((sub.route == "direct").sum()),
                "n_quantified_remapped": int((sub.route == "remapped").sum())}
    class_n = cat.primary_evidence_class.value_counts().to_dict()
    t1 = [row(c, class_n[c], per_gene[per_gene["class"] == c]) for c in CLASSES]
    extra = EXTRA_ROWS if remap_used else [DIRECT_ONLY_ROW]
    for lab in extra:
        r = row(lab, None, assign[assign.route == lab])
        r["n_genes"] = None
        t1.append(r)
    t1_total = row("all protein rows", None, assign)
    t1_total["n_genes"] = None
    t1_out = pd.DataFrame(t1 + [t1_total])
    t1 = pd.DataFrame(t1)
    sum_q, sum_s = int(t1.n_quantified.sum()), int(t1.n_protein_sig.sum())
    ok_d8 = R.gate("D8/S6", "T1 rows sum to all protein rows and significant rows; one protein row per Catalog gene; "
                   "remap-added significant genes within the 7 listed",
                   f"6071 / 290 / 0 duplicated genes / subset of {sorted(SIG_REMAP_ALLOWED)}",
                   f"{sum_q} / {sum_s} / {dup_gene} duplicated genes / {sorted(sig_added)}",
                   sum_q == 6071 and sum_s == 290 and dup_gene == 0 and sig_added <= SIG_REMAP_ALLOWED
                   and assign.route.notna().all(), "any failure stops (S6)")
    if not ok_d8:
        return finish(R, "S6", "conservation failed: the join double counts or drops protein rows")
    R.w(t1_out, "T1_protein_coverage_by_class.tsv")

    sub_cl = per_gene[per_gene["class"].isin(CLASSED)]
    t1b = pd.DataFrame([r for r in t1.to_dict("records") if r["row"] in CLASSED] +
                       [row("classed genes (concordant + established + genetically anchored)",
                            sum(class_n[c] for c in CLASSED), sub_cl)])
    R.w(t1b, "T1b_classed_genes.tsv")

    # ---------------- D5 (blind): direct matches HGNChelper marks not approved ---------------------
    direct_syms = sorted(assign.protein_symbol[assign.route == "direct"])
    if R.debug:
        direct_syms = sorted(rng.choice(direct_syms, size=300, replace=False).tolist())
    d5 = hall_idx.loc[direct_syms, ["Approved", "Suggested.Symbol", "n_sugg"]].reset_index()
    d5 = d5[~d5.Approved.astype(str).str.upper().eq("TRUE")].copy()
    d5["suggestion_is_other_catalog_symbol"] = [(n == 1 and s != x and s in cat_sym)
                                                 for x, s, n in zip(d5.x, d5["Suggested.Symbol"], d5.n_sugg)]
    R.w(d5, "D5_direct_not_approved.tsv")
    R.result["d5"] = {"n_direct_checked": len(direct_syms), "n_not_approved": len(d5),
                      "n_suggestion_other_catalog_symbol": int(d5.suggestion_is_other_catalog_symbol.sum()),
                      "above_20_limit_of_join": len(d5) > 20}
    R.say(f"[D5] direct matches checked {len(direct_syms)}; not approved by HGNChelper {len(d5)}; "
          f"of which single suggestion is another Catalog symbol {int(d5.suggestion_is_other_catalog_symbol.sum())}; "
          f"prediction [Estimate] fewer than 20")

    # ---------------- D6: single suggestions that are not Catalog symbols ---------------------------
    d6 = um[um.single & ~um.single_in_catalog][["x", "Suggested.Symbol"]].copy()
    d6["in_quarantine"] = [(x in q_src) or (s in q_src) for x, s in zip(d6.x, d6["Suggested.Symbol"])]
    R.w(d6, "D6_single_not_in_catalog.tsv")
    R.result["d6"] = {"n": len(d6), "in_quarantine": int(d6.in_quarantine.sum()),
                      "absent": int((~d6.in_quarantine).sum()), "CAST_in_quarantine": bool(d6.loc[d6.x == "CAST", "in_quarantine"].any())}
    R.say(f"[D6] {R.result['d6']}; quarantine: {d6.x[d6.in_quarantine].tolist()}; absent: {d6.x[~d6.in_quarantine].tolist()}")

    # ---------------- D7 (blind, descriptive): deposit detection ---------------------------------
    dep = pd.read_csv(R.out / "deposit_symbols.tsv", sep="\t", keep_default_na=False, na_values=[""])
    dep["catalog_symbol"] = np.where(dep.symbol.isin(cat_sym), dep.symbol, None)
    if remap_used:
        sug = hall_idx["Suggested.Symbol"]
        nsug = hall_idx["n_sugg"]
        for i, s in dep.symbol.items():
            if dep.at[i, "catalog_symbol"] is None and nsug.get(s, 0) == 1 and sug.get(s) in cat_sym:
                dep.at[i, "catalog_symbol"] = sug.get(s)
    det_any = set(dep.catalog_symbol.dropna())
    det_single = set(dep.catalog_symbol[~dep.multi_gene_group.astype(str).str.upper().eq("TRUE")].dropna())
    sig_genes = set(per_gene.catalog_symbol[per_gene.sig])
    d7 = []
    for c in CLASSES:
        g = set(cat.symbol[cat.primary_evidence_class == c])
        d7.append({"class": c, "n_genes": len(g), "n_detected_in_deposit": len(g & det_any),
                   "pct_detected": round(100 * len(g & det_any) / len(g), 2),
                   "n_detected_only_in_multigene_group": len((g & det_any) - det_single),
                   "n_quantified": len(g & quant_sym), "n_detected_not_quantified": len((g & det_any) - quant_sym),
                   "n_quantified_not_detected": len((g & quant_sym) - det_any), "n_protein_sig": len(g & sig_genes)})
    d7 = pd.DataFrame(d7)
    R.w(d7, "D7_deposit_detection_by_class.tsv")
    R.say("[D7]\n" + d7.to_string(index=False))
    if int(d7.n_quantified_not_detected.sum()) != 0:
        R.say("[D7] WARNING: quantified genes not found in the deposit Genes field:", int(d7.n_quantified_not_detected.sum()))

    verdict = "T1_REMAPPED" if remap_used else f"DIRECT_JOIN_ONLY_{remap_off}"
    return finish(R, verdict, None, t1_out=t1_out, t1b=t1b, d7=d7, remap_off=remap_off)


def md_table(df):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if v is None or (isinstance(v, float) and np.isnan(v)):
                cells.append("")
            elif isinstance(v, (bool, np.bool_)):
                cells.append("yes" if v else "no")
            elif isinstance(v, (int, np.integer)) or (isinstance(v, float) and float(v).is_integer() and c.startswith("n_")):
                cells.append(f"{int(v):,}")
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def finish(R, code, reason, t1_out=None, t1b=None, d7=None, remap_off=None):
    stops = {"S1": "untestable: input hash mismatch (S1); no number is used",
             "G1": "untestable: input counts differ from the prespec (G1); no number is used",
             "S2": "untestable: direct-join reproduction failed (S2); no number is used",
             "S5": "untestable: known-answer controls failed (S5); no coverage number is reported",
             "S6": "untestable: conservation failed (S6); the join double counts or drops rows",
             "S7": "untestable: versions or D4 remap counts differ (S7); amendment needed before T1",
             "S8": "untestable: D9 differs (S8); stop before T1"}
    if code in stops:
        verdict = stops[code]
        status = "stopped_by_rule"
    elif code == "T1_REMAPPED":
        verdict = ("completed: all gates and controls passed; T1 reported with the HGNChelper alias remap. "
                   "Resource annotation, no test.")
        status = "completed"
    else:
        verdict = (f"completed with fallback ({remap_off}): the remap is not used; T1 is the direct join, labelled "
                   "'outdated protein symbols not remapped'.")
        status = "completed_direct_join_only"
    R.result.update({"verdict_code": code, "verdict": verdict, "status": status, "stop_reason": reason,
                     "gates": R.gates})
    if t1_out is not None:
        R.result["T1"] = json.loads(t1_out.to_json(orient="records"))
        R.result["T1b"] = json.loads(t1b.to_json(orient="records"))
        R.result["D7"] = json.loads(d7.to_json(orient="records"))
    (R.out / f"{R.prefix}result_H86.json").write_text(json.dumps(R.result, indent=1, default=str))

    L = []
    title = "RESULT H86" + (" (DEBUG: classes permuted, not a result)" if R.debug else "")
    L.append(f"# {title}\n")
    L.append("Gene Catalog r8 evidence classes against PXD051911 liver protein quantification and protein BH q < 0.05.")
    L.append(f"Prespec `PRESPEC_H86.md` sha256 `{PRESPEC_SHA}`. Seed {R.seed}; no step is random. "
             "Catalog r8 is a candidate release (corrected COLOC rerun pending); every number here is a candidate annotation.\n")
    L.append(f"**Verdict (computed by the prespec rule): {verdict}**\n")
    if reason:
        L.append(f"Stop reason: {reason}\n")
    L.append("## Gates, controls and diagnostics in prespec order\n")
    L.append(md_table(pd.DataFrame(R.gates)[["id", "check", "predicted", "observed", "rule", "pass"]]))
    L.append("")
    if t1_out is not None:
        L.append("## T1. Protein coverage by evidence class\n")
        L.append("Unit: Catalog gene (counting only; genes are not replicates). Protein-significant = quantified and "
                 "`padj < 0.05`, the inherited Benjamini-Hochberg family over the 6,071 proteins of the PXD051911 fit "
                 "(58 patients, 12 control / 46 MASLD; covariates acquisition batch, age, BMI, sex). H86 adds no test "
                 "and no new p value. Non-class rows count protein rows, not Catalog genes.\n")
        L.append(md_table(t1_out))
        L.append("\n## T1b. Classed genes only\n")
        L.append(md_table(t1b))
        L.append("\n## D7. Deposit detection (descriptive)\n")
        L.append(md_table(d7))
        L.append("")
        d = R.result
        L.append("## Diagnostics\n")
        L.append(f"- D4: {d['d4']}; ambiguous {d['d4_ambiguous']}; disagreements (HGNChelper used) {d['d4_disagreements']}.")
        L.append(f"- D5: {d['d5']} (prediction [Estimate]: fewer than 20).")
        L.append(f"- D6: {d['d6']}.")
        L.append(f"- D9: {d['d9']}.")
        L.append(f"- NC-2: {d['nc2']}.")
        L.append("")
        t = t1_out.set_index("row")
        e, gq = t.loc["established_state_associated"], t.loc["genetically_anchored"]
        if code == "T1_REMAPPED":
            sent = (f"Protein coverage depends on the evidence class. Of the 1,129 established state-associated genes, "
                    f"{int(e.n_quantified):,} were quantified in the PXD051911 liver proteome and {int(e.n_protein_sig):,} "
                    f"of those reached BH q < 0.05; for the 400 genetically anchored genes the counts are "
                    f"{int(gq.n_quantified):,} and {int(gq.n_protein_sig):,}. A gene without a protein row was not "
                    "quantified and has no protein evidence in either direction.")
        else:
            sent = ("119 protein symbols do not match a Catalog symbol: 95 are previous or alias gene names (for example "
                    "MARC1, now MTARC1) and 24 are current names of genes absent from the Catalog. "
                    f"Of the 1,129 established state-associated genes, {int(e.n_quantified):,} were quantified "
                    f"(direct symbol join) and {int(e.n_protein_sig):,} reached BH q < 0.05; for the 400 genetically "
                    f"anchored genes the counts are {int(gq.n_quantified):,} and {int(gq.n_protein_sig):,}.")
        L.append("## Proposed sentence (filled from T1 by code; not installed)\n")
        L.append(sent)
        R.result["proposed_sentence"] = sent
        (R.out / f"{R.prefix}result_H86.json").write_text(json.dumps(R.result, indent=1, default=str))
        L.append("\nClasses are not compared statistically: coverage differences mostly reflect expression level and "
                 "protein abundance [Inference], not class biology. With 12 controls, protein nonsignificance is not "
                 "evidence of no change. Significance columns go to Supplementary T1 and the Figure 5B legend only, "
                 "never into the Catalog summary (proteomics is program_context_only).")
    else:
        L.append("No coverage number is reported (stop rule).")
    (R.out / f"{R.prefix}RESULT_H86.md").write_text("\n".join(L) + "\n")
    (R.out / f"{R.prefix}diagnostics.log").write_text("\n".join(R.lines) + "\n")
    R.say(f"[02] verdict: {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
