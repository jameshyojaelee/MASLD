#!/usr/bin/env python3
"""H66: per-gene follow-up assay observability for every Gene Catalog r8 gene.

Implements PRESPEC_H66.md (sha256 66be821e9195362a44c837835ed2b2e2c730da96e83e5ce6c3040b5594fd6b6c)
exactly: gates and stop rules in order (S1, S6, G1, G2, G4a, G4), four fields (liver proteome,
plasma proteome, CosMx panel, donor-collapsed single-cell lineage), acceptance checks A1 to A4,
NC-1 (A2) and NC-2, the descriptive secondary class models with Holm over 4, common-support
sensitivity, D1 to D7 diagnostics, and RESULT_H66.md with the verdict computed in code.

Modes
  full  (default; the prespec command line): everything, writes RESULT_H66.md.
  gates : S1, S6, G1, G2, G4a, G4 and the mapping diagnostics only; no field is built.
  null  : debug only. Gates on real data, then proteome rows, pseudobulk row labels and Catalog
          classes are permuted before any field is built; small B; writes RESULT_NULLDEBUG.md.

The unit of the annotation is the gene. Inside each field the unit is the participant or donor.
"""
import argparse
import gzip
import hashlib
import json
import os
import re
import sys
import time
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from statsmodels.tools.sm_exceptions import ConvergenceWarning

try:
    from statsmodels.tools.sm_exceptions import PerfectSeparationWarning
except ImportError:  # older statsmodels
    PerfectSeparationWarning = None

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
H66 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRESPEC = os.path.join(H66, "PRESPEC_H66.md")
PRESPEC_SHA = "66be821e9195362a44c837835ed2b2e2c730da96e83e5ce6c3040b5594fd6b6c"
ROUND4 = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z"
PREFIX = {
    ".../pseudobulk/": "Analysis/SingleCell/results_gpu_v2/pseudobulk/",
    ".../gate/NG05/": f"{ROUND4}/gate/NG05/",
    ".../write/H66/": f"{ROUND4}/write/H66/",
}
LINEAGES = ["B_cells", "Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Hepatocytes",
            "Macrophages", "Plasma_cells", "Resident_NK", "T_cells"]
CC_NOT_USED = ["Basophils", "Circulating_NK/NKT", "Mono+mono_derived_cells", "Neutrophils",
               "cDC1s", "cDC2s", "pDCs"]
ANCH, EST = "genetically_anchored", "established_state_associated"
OR_PATTERN = r"^OR\d+[A-Z]+\d+[A-Z]?\d*$"
OR_PATTERN_SHORT = r"^OR\d+[A-Z]+\d+"

# Fixed CosMx shared-probe expansion (prespec, "Field definitions")
COSMX_EXP = {
    "CCL3/L1/L3": ["CCL3", "CCL3L1", "CCL3L3"], "CCL4/L1/L2": ["CCL4", "CCL4L1", "CCL4L2"],
    "CXCL1/2/3": ["CXCL1", "CXCL2", "CXCL3"], "EIF5A/L1": ["EIF5A", "EIF5AL1"],
    "FCGR3A/B": ["FCGR3A", "FCGR3B"], "HBA1/2": ["HBA1", "HBA2"], "HCAR2/3": ["HCAR2", "HCAR3"],
    "HLA-DQB1/2": ["HLA-DQB1", "HLA-DQB2"], "HSPA1A/B": ["HSPA1A", "HSPA1B"],
    "IFNA1/13": ["IFNA1", "IFNA13"], "IFNL2/3": ["IFNL2", "IFNL3"],
    "KRT6A/B/C": ["KRT6A", "KRT6B", "KRT6C"], "MAP1LC3B/2": ["MAP1LC3B", "MAP1LC3B2"],
    "MZT2A/B": ["MZT2A", "MZT2B"], "PF4/V1": ["PF4", "PF4V1"], "SAA1/2": ["SAA1", "SAA2"],
    "TNXA/B": ["TNXA", "TNXB"], "TPSAB1/B2": ["TPSAB1", "TPSB2"], "XCL1/2": ["XCL1", "XCL2"],
}
COSMX_HEADER_MD5 = "80ef1b7875ecc4eedd3bc8edf48404b7"

KA_LIVER = ["ALB", "APOB", "CYP3A4", "PNPLA3", "HSD17B13", "SERPINA1", "MTARC1"]
KA_LIVER_REPORT = ["TM6SF2"]
KA_PLASMA = ["ALB", "APOA1", "APOB", "TF", "HP", "FGA", "CRP", "SERPINA1"]
KA_SC = [("KRT19", "Cholangiocytes"), ("EPCAM", "Cholangiocytes"), ("LUM", "Fibroblasts"),
         ("DCN", "Fibroblasts"), ("COL1A1", "Fibroblasts"), ("ALB", "Hepatocytes"),
         ("PECAM1", "Endothelial_cells"), ("CD68", "Macrophages"), ("CD3E", "T_cells"),
         ("MS4A1", "B_cells"), ("JCHAIN", "Plasma_cells")]
KA_SC_REQUIRED = {("KRT19", "Cholangiocytes"), ("EPCAM", "Cholangiocytes")}
KA_ALL = ["ALB", "APOB", "CYP3A4", "PNPLA3", "HSD17B13", "SERPINA1", "TM6SF2", "MTARC1", "APOA1",
          "TF", "HP", "FGA", "CRP", "KRT19", "EPCAM", "LUM", "DCN", "COL1A1", "PECAM1", "CD68",
          "CD3E", "MS4A1", "JCHAIN"]

# Thresholds (prespec; never tuned)
SC_MIN_CELLS, SC_MIN_LIB, SC_MIN_COUNT, SC_MIN_CPM, SC_DONOR_FRAC, SC_MIN_ELIG = 20, 10_000, 5, 10.0, 0.5, 10
SC_COL_MIN_CELLS = 10  # G4a: a pseudobulk column is a run with >= 10 cells of the lineage
G4_MIN_RHO = 0.5
A2_MAX_FRAC = 0.02
A3_MIN_DIFF, A3_MIN_RHO = 0.30, 0.80
NC2_OR_LO, NC2_OR_HI, NC2_P_LO, NC2_P_HI = 0.85, 1.15, 0.02, 0.10
MIN_EVENTS = 10

# G1 targets: every number in "n by group" (prespec), from the audit outputs
G1 = {
    "catalog_n": 28233, "catalog_unique_ens": 28233, "catalog_unique_sym": 28233,
    "class_untested": 13551, "class_unresolved": 13127, "class_established_state_associated": 1129,
    "class_genetically_anchored": 400, "class_concordant": 26,
    "catalog_ens_in_hgnc": 23136,
    "secondary_n": 1529, "secondary_anchored": 400, "secondary_established": 1129,
    "secondary_aveexpr_missing": 0, "secondary_pc_anchored": 389, "secondary_pc_established": 1083,
    "liver_groups": 7096, "liver_route_accession": 7076, "liver_route_current": 14,
    "liver_route_prev": 2, "liver_prev_symbols": ["CYHR1", "VNN3"], "liver_route_unmapped": 4,
    "liver_unmapped_no_gene_name": 2, "liver_unmapped_named": ["IGHV3-43D", "PALM2"],
    "liver_ens_genes": 7187, "liver_ens_in_catalog": 7133, "liver_multi_gene_groups": 64,
    "liver_sample_cols": 58, "liver_participants": 58, "liver_sample_group": {"initial_sample": 58},
    "liver_saf": {"MASH": 17, "MASL": 29, "No_MASLD": 12},
    "plasma_groups": 578, "plasma_route_accession": 573, "plasma_route_current": 1,
    "plasma_route_prev": 0, "plasma_route_unmapped": 4, "plasma_unmapped_all_IGHV": True,
    "plasma_ens_genes": 597, "plasma_ens_in_catalog": 578, "plasma_multi_gene_groups": 19,
    "plasma_sample_cols": 184, "plasma_baseline_cols": 143, "plasma_baseline_participants": 143,
    "plasma_baseline_sample_group": {"at_bariatric_surgery_sample": 7, "initial_sample": 136},
    "plasma_followup_cols": 41, "plasma_participants": 143, "plasma_followup_in_baseline": 41,
    "plasma_cols_unassigned": 0,
    "cosmx_header_cols": 1178, "cosmx_fov_cellid": 2, "cosmx_systemcontrol": 197, "cosmx_negative": 11,
    "cosmx_targets": 968, "cosmx_single_current": 944, "cosmx_single_prev": ["AATK", "DDX58", "IL1RL1"],
    "cosmx_single_prev_to": ["IL33R", "LMTK1", "RIGI"],
    "cosmx_shared_probes": 19, "cosmx_shared_genes": 42, "cosmx_shared_genes_in_catalog": 40,
    "cosmx_shared_genes_not_in_catalog": ["CCL3L1", "CCL4L1"], "cosmx_no_gene": ["HLA-DRB", "MHC I"],
    "cosmx_header_md5_all_four": COSMX_HEADER_MD5,
    "sc_runs": 275, "sc_donors": 104,
    "sc_donors_by_dataset": {"GSE136103": 10, "GSE174748": 4, "GSE185477": 3, "GSE189600": 4,
                             "GSE202379": 46, "GSE244832": 18, "Liver_Atlas": 19},
    "sc_gse289173_present": False,
    "sc_pb_cols": {"B_cells": 59, "Cholangiocytes": 214, "Endothelial_cells": 245, "Fibroblasts": 218,
                   "Hepatocytes": 240, "Macrophages": 196, "Plasma_cells": 122, "Resident_NK": 98,
                   "T_cells": 161},
    "sc_pb_donors": {"B_cells": 27, "Cholangiocytes": 81, "Endothelial_cells": 96, "Fibroblasts": 88,
                     "Hepatocytes": 93, "Macrophages": 80, "Plasma_cells": 57, "Resident_NK": 54,
                     "T_cells": 76},
    "sc_eligible": {"B_cells": 25, "Cholangiocytes": 73, "Endothelial_cells": 95, "Fibroblasts": 80,
                    "Hepatocytes": 85, "Macrophages": 70, "Plasma_cells": 42, "Resident_NK": 43,
                    "T_cells": 75},
    "sc_eligible_50cells": {"B_cells": 21, "Cholangiocytes": 60, "Endothelial_cells": 88,
                            "Fibroblasts": 65, "Hepatocytes": 80, "Macrophages": 59, "Plasma_cells": 25,
                            "Resident_NK": 33, "T_cells": 64},
    "sc_lineages_used": 9,
    "sc_all_integer": True, "sc_all_cols_in_cellcount": True, "sc_cols_gse289173": 0,
    "pb_rows": 37533, "pb_rows_identical_9": True, "pb_rows_ensembl": 12468, "pb_rows_symbol": 25065,
    "pb_catalog_genes_with_row": 23593, "pb_rows_route_ensembl_in_catalog": 2766,
    "pb_rows_route_current_in_catalog": 20594, "pb_rows_route_prev_in_catalog": 291,
    "pb_rows_route_catsym_in_catalog": 16, "pb_catalog_genes_multi_row": 70,
    "pb_catalog_genes_symbol_only_rule": 20811, "pb_catalog_genes_no_row": 4640,
    "ka_in_catalog": 23,
    "ka_classes": {"established_state_associated": 4, "unresolved": 18, "untested": 1},
    "ka_established": ["COL1A1", "CYP3A4", "EPCAM", "LUM"], "ka_untested": ["CD68"],
    "ka_sc_symbol_row": 11,
    "or_genes": 387, "or_classes": {"established_state_associated": 1, "unresolved": 4, "untested": 382},
    "or_genes_with_pb_row": 292, "or_genes_short_pattern": 388,
    "nc2_pool_total": 13553, "nc2_pool_inside_range": 13522, "nc2_pool_min_per_decile": 447,
    "nc2_need_max_per_decile": 153,
}


# ----------------------------------------------------------------------------------------------
class Diag:
    STEPS = {"1": "input sha256", "2": "G1 counts and G2 PC-A", "3": "D1 mapping", "4": "D2 symbol-only",
             "5": "D3 CosMx targets", "6": "D4 CosMx headers", "7": "D5 single-cell donors, G4a, G4",
             "8": "D6 single-cell sensitivity", "9": "A1 A2 A3 A4", "10": "NC-2",
             "11": "coverage by class, secondary models, common support", "12": "D7 bootstrap"}

    def __init__(self, out):
        self.out = out
        self.rows = []
        self.fh = open(os.path.join(out, "h66.log"), "a")

    def msg(self, s):
        line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {s}"
        print(line, flush=True)
        self.fh.write(line + "\n")
        self.fh.flush()

    def add(self, step, check, value, expected=None, ok=None, note=""):
        enc = lambda v: v if isinstance(v, (str, int, float, bool, type(None))) else json.dumps(v, sort_keys=True, default=str)
        self.rows.append(dict(step=int(step), step_name=self.STEPS[str(step)], check=check, value=enc(value),
                              expected=enc(expected), passed=ok, note=note))
        tail = "" if expected is None else f" | expected {enc(expected)}"
        tag = "" if ok is None else (" | PASS" if ok else " | FAIL")
        self.msg(f"[D{step}] {check}: {enc(value)}{tail}{tag}{(' | ' + note) if note else ''}")

    def flush(self):
        d = pd.DataFrame(self.rows)
        if len(d):
            d = d.sort_values("step", kind="stable")
            d.to_csv(os.path.join(self.out, "diagnostics.tsv"), sep="\t", index=False)
            self.msg("==== diagnostics in prespecified order (also diagnostics.tsv) ====")
            for r in d.itertuples():
                self.fh.write(f"  {r.step}\t{r.step_name}\t{r.check}\t{r.value}\t{r.expected}\t{r.passed}\t{r.note}\n")
            self.fh.flush()


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_input_table(text):
    rows, in_tab = [], False
    for line in text.splitlines():
        if line.startswith("| Role | Path | Bytes | sha256 |"):
            in_tab = True
            continue
        if in_tab:
            if not line.startswith("|"):
                break
            if line.startswith("|---"):
                continue
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if len(cells) != 4:
                raise RuntimeError(f"bad input-table row: {line}")
            role, path, nb, sha = cells
            path, sha = path.strip("`"), sha.strip("`")
            for k, v in PREFIX.items():
                if path.startswith(k):
                    path = v + path[len(k):]
            if path.startswith("..."):
                raise RuntimeError(f"unexpanded path {path}")
            rows.append(dict(role=role, path=path, bytes=int(nb.replace(",", "")) if nb else None, sha256=sha))
    return pd.DataFrame(rows)


def write_result(out, mode, title_status, verdict, body_lines):
    name = "RESULT_H66.md" if mode == "full" else ("RESULT_NULLDEBUG.md" if mode == "null" else "RESULT_GATES_DEBUG.md")
    with open(os.path.join(out, name), "w") as fh:
        fh.write(f"# H66 result{'' if mode == 'full' else ' (' + mode + ' debug run; not a result)'}\n\n")
        fh.write(f"Status: {title_status}\n\nVerdict (computed in code by the prespec rule): {verdict}\n\n")
        fh.write("\n".join(body_lines) + "\n")
    return name


def fit_logit(df, formula, term="anch"):
    """Return dict with OR, CI, p, converged flag; None values when singular or not converged."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            m = smf.logit(formula, data=df).fit(disp=0, maxiter=200)
        except Exception as e:  # singular matrix, perfect separation raised as error
            return dict(ok=False, reason=type(e).__name__ + ": " + str(e)[:120])
    bad = [str(x.message)[:80] for x in w if issubclass(x.category, ConvergenceWarning)
           or (PerfectSeparationWarning is not None and issubclass(x.category, PerfectSeparationWarning))]
    conv = bool(m.mle_retvals.get("converged", False))
    b, se = float(m.params[term]), float(m.bse[term])
    if not conv or bad or not np.isfinite(b) or not np.isfinite(se) or se <= 0:
        return dict(ok=False, reason="not converged or separation: " + ";".join(bad)[:200], b=b, se=se)
    return dict(ok=True, OR=float(np.exp(b)), lo=float(np.exp(b - 1.959963984540054 * se)),
                hi=float(np.exp(b + 1.959963984540054 * se)), p=float(m.pvalues[term]), b=b, se=se)


# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20260923)
    ap.add_argument("--mode", choices=["full", "gates", "null"], default="full")
    ap.add_argument("--nc2-draws", type=int, default=1000)
    ap.add_argument("--boot", type=int, default=1000)
    a = ap.parse_args()
    out = a.out
    if not os.path.isdir(out):
        sys.exit(f"output directory missing: {out}")
    if a.mode == "full" and (a.nc2_draws != 1000 or a.boot != 1000 or a.seed != 20260923):
        sys.exit("full mode must use the prespecified seed and B (1,000 NC-2 draws, 1,000 bootstrap resamples)")
    D = Diag(out)
    D.msg(f"H66 start mode={a.mode} seed={a.seed} nc2={a.nc2_draws} boot={a.boot} out={out}")
    stops = []  # (rule, text)

    # ------------------------------------------------------------------ S1 prespec and inputs
    got = sha256_file(PRESPEC)
    D.add(1, "PRESPEC_H66.md sha256", got, PRESPEC_SHA, got == PRESPEC_SHA)
    if got != PRESPEC_SHA:
        stops.append(("S1", "PRESPEC_H66.md sha256 differs"))
        return finish_stop(a, D, out, stops)
    ptxt = open(PRESPEC).read()
    tab = parse_input_table(ptxt)
    D.add(1, "input-table rows parsed", len(tab), 35, len(tab) == 35)
    if len(tab) != 35:
        stops.append(("S1", "input table did not parse to 35 rows"))
        return finish_stop(a, D, out, stops)
    res = []
    for r in tab.itertuples():
        p = os.path.join(ROOT, r.path)
        if not os.path.exists(p):
            res.append(dict(role=r.role, path=r.path, bytes_expected=r.bytes, bytes=None, sha256_expected=r.sha256, sha256=None, ok=False))
            continue
        nb = os.path.getsize(p)
        s = sha256_file(p)
        ok = (s == r.sha256) and (pd.isna(r.bytes) or nb == int(r.bytes))
        res.append(dict(role=r.role, path=r.path, bytes_expected=r.bytes, bytes=nb, sha256_expected=r.sha256, sha256=s, ok=ok))
    res = pd.DataFrame(res)
    res.to_csv(os.path.join(out, "input_sha256.tsv"), sep="\t", index=False)
    for r in res.itertuples():
        D.add(1, f"sha256 {r.path}", r.sha256, r.sha256_expected, bool(r.ok))
    if not res.ok.all():
        stops.append(("S1", "sha256 differs for: " + ", ".join(res.loc[~res.ok, "path"])))
        return finish_stop(a, D, out, stops)
    P = {r.role: os.path.join(ROOT, r.path) for r in tab.itertuples()}
    PB = {L: os.path.join(ROOT, f"Analysis/SingleCell/results_gpu_v2/pseudobulk/{L}_pseudobulk.csv") for L in LINEAGES}
    for L in LINEAGES:
        if PB[L] not in set(P.values()):
            raise RuntimeError(f"pseudobulk {L} not in hashed input table")
    COSMX = [os.path.join(ROOT, "data/external/govaere2026_natgenetics/cosmx", f) for f in
             ["GSM9351704_Leuven_1_exprMat_file.csv.gz", "GSM9351705_Leuven_2_exprMat_file.csv.gz",
              "GSM9351706_Leuven_3_exprMat_file.csv.gz", "GSM9351707_Leuven_4_exprMat_file.csv.gz"]]
    for f in COSMX:
        assert f in set(P.values())
    run_map_path = os.path.join(out, "run_to_donor.tsv")
    if not os.path.exists(run_map_path):
        sys.exit("run_to_donor.tsv missing: run 00_donor_map.R first")
    r_in = pd.read_csv(os.path.join(out, "00_donor_map_input_sha256.tsv"), sep="\t")
    D.add(1, "00_donor_map.R input sha256 all match", bool(r_in.ok.all()), True, bool(r_in.ok.all()))

    # ------------------------------------------------------------------ load
    cat = pd.read_parquet(P["Catalog r8 gene index"], columns=["ensembl_id", "symbol", "primary_evidence_class"])
    cat["cls"] = cat["primary_evidence_class"].astype(str)
    h = pd.read_csv(P["HGNC complete set (2026-08-27 snapshot)"], sep="\t", dtype=str, low_memory=False)
    h = h[h.status == "Approved"]
    hens = h.dropna(subset=["ensembl_gene_id"])
    hgnc_ens = set(hens.ensembl_gene_id)
    acc2ens = {}
    for e, u in zip(hens.ensembl_gene_id, hens.uniprot_ids):
        if isinstance(u, str):
            for x in u.split("|"):
                acc2ens.setdefault(x.strip(), set()).add(e)
    cur = dict(zip(hens.symbol, hens.ensembl_gene_id))
    prev = {}
    for e, v in zip(hens.ensembl_gene_id, hens.prev_symbol):
        if isinstance(v, str):
            for s in v.split("|"):
                prev.setdefault(s.strip(), set()).add(e)
    ens2hgncsym = dict(zip(hens.ensembl_gene_id, hens.symbol))
    gen = pd.read_csv(P["GENCODE v49 biotype"], sep="\t")
    deg = pd.read_csv(P["Bulk AveExpr (adopted refit)"])
    deg["ens"] = deg["gene"].str.split(".").str[0]
    c2 = cat.merge(gen[["ensembl_base", "gene_biotype"]], left_on="ensembl_id", right_on="ensembl_base", how="left") \
            .merge(deg[["ens", "AveExpr"]], left_on="ensembl_id", right_on="ens", how="left")
    if len(c2) != len(cat):
        raise RuntimeError("biotype/AveExpr merge changed the Catalog row count")
    c2 = c2.drop(columns=["ensembl_base", "ens"])
    c2["pc"] = (c2.gene_biotype == "protein_coding").astype(int)
    cens = set(cat.ensembl_id)
    cat_nohgnc = c2[~c2.ensembl_id.isin(hgnc_ens)]
    nohgnc_sym2ens = dict(zip(cat_nohgnc.symbol, cat_nohgnc.ensembl_id))
    sym2cat = dict(zip(cat.symbol, cat.ensembl_id))
    meta = pd.read_csv(P["PXD051911 sample metadata"], sep="\t", dtype=str)
    runmap = pd.read_csv(run_map_path, sep="\t")
    cc = pd.read_csv(P["Run-level lineage cell counts"], sep="\t")
    g1 = {}

    # ------------------------------------------------------------------ S6 GSE289173
    run2ds = dict(zip(runmap.run, runmap.dataset))
    run2donor = dict(zip(runmap.run, runmap.donor))
    s6 = bool(runmap.dataset.astype(str).str.contains("GSE289173").any()) or any("GSE289173" in str(x) for x in cc["sample"])
    g1["sc_gse289173_present"] = s6

    # ------------------------------------------------------------------ G1 Catalog and secondary
    g1["catalog_n"] = len(cat)
    g1["catalog_unique_ens"] = int(cat.ensembl_id.nunique())
    g1["catalog_unique_sym"] = int(cat.symbol.nunique())
    vc = cat.cls.value_counts()
    for k in ["untested", "unresolved", "established_state_associated", "genetically_anchored", "concordant"]:
        g1[f"class_{k}"] = int(vc.get(k, 0))
    g1["catalog_ens_in_hgnc"] = len(cens & hgnc_ens)
    sec = c2[c2.cls.isin([ANCH, EST])].copy()
    g1["secondary_n"] = len(sec)
    g1["secondary_anchored"] = int((sec.cls == ANCH).sum())
    g1["secondary_established"] = int((sec.cls == EST).sum())
    g1["secondary_aveexpr_missing"] = int(sec.AveExpr.isna().sum())
    g1["secondary_pc_anchored"] = int(sec[(sec.cls == ANCH)].pc.sum())
    g1["secondary_pc_established"] = int(sec[(sec.cls == EST)].pc.sum())

    # ------------------------------------------------------------------ proteome mapping
    def map_proteome(t):
        recs, drop_amb = [], []
        for i, (acc, g) in enumerate(zip(t.ProteinAccessions, t.Genes)):
            accs = [re.sub(r"-\d+$", "", x.strip()) for x in str(acc).split(";")]
            es = set().union(*[acc2ens.get(x, set()) for x in accs])
            syms = [s.strip() for s in str(g).split(";") if s.strip() and s.strip() != "nan"] if isinstance(g, str) else []
            route = "accession"
            if not es:
                es = set().union(*[{cur[s]} for s in syms if s in cur]) if syms else set()
                route = "current_symbol"
                if not es:
                    for s in syms:
                        pp = prev.get(s, set())
                        if len(pp) == 1:
                            es |= pp
                        elif len(pp) > 1:
                            drop_amb.append(dict(group=i, symbol=s, n_genes=len(pp), genes=";".join(sorted(pp))))
                    route = "prev_symbol" if es else "unmapped"
            d2 = {nohgnc_sym2ens[s] for s in syms if s in nohgnc_sym2ens}
            recs.append(dict(group=i, ProteinAccessions=acc, Genes=g if isinstance(g, str) else "", route=route,
                             hgnc_ens=";".join(sorted(es)), n_hgnc_ens=len(es),
                             d2_catalog_symbol_ens=";".join(sorted(d2)), all_ens=";".join(sorted(es | d2)),
                             n_all_ens=len(es | d2)))
        return pd.DataFrame(recs), pd.DataFrame(drop_amb, columns=["group", "symbol", "n_genes", "genes"])

    prot = {}
    for tag, role in [("liver", "Liver proteome"), ("plasma", "Plasma proteome (unfiltered)")]:
        t = pd.read_csv(P[role], sep="\t")
        samp = list(t.columns[3:])
        if list(t.columns[:3]) != ["ProteinAccessions", "Genes", "ProteinDescriptions"]:
            raise RuntimeError(f"{tag}: unexpected first columns")
        mp, amb = map_proteome(t)
        prot[tag] = dict(t=t, samp=samp, map=mp, amb=amb)
        g1[f"{tag}_groups"] = len(t)
        rc = mp.route.value_counts()
        g1[f"{tag}_route_accession"] = int(rc.get("accession", 0))
        g1[f"{tag}_route_current"] = int(rc.get("current_symbol", 0))
        g1[f"{tag}_route_prev"] = int(rc.get("prev_symbol", 0))
        g1[f"{tag}_route_unmapped"] = int(rc.get("unmapped", 0))
        ens_all = set(";".join(mp.hgnc_ens[mp.hgnc_ens != ""]).split(";")) - {""}
        g1[f"{tag}_ens_genes"] = len(ens_all)
        g1[f"{tag}_ens_in_catalog"] = len(ens_all & cens)
        g1[f"{tag}_multi_gene_groups"] = int((mp.n_hgnc_ens > 1).sum())
        g1[f"{tag}_sample_cols"] = len(samp)
        um = mp[mp.route == "unmapped"]
        if tag == "liver":
            g1["liver_prev_symbols"] = sorted(mp.loc[mp.route == "prev_symbol", "Genes"])
            g1["liver_unmapped_no_gene_name"] = int((um.Genes == "").sum())
            g1["liver_unmapped_named"] = sorted(um.Genes[um.Genes != ""])
            m = meta[meta.liver_proteomics_filename.notna()]
            if set(m.liver_proteomics_filename) != set(samp) or len(m) != len(samp):
                raise RuntimeError("liver: metadata files do not equal quant columns")
            g1["liver_participants"] = int(m.patient_name.nunique())
            g1["liver_sample_group"] = {k: int(v) for k, v in m.sample_group.value_counts().items()}
            g1["liver_saf"] = {k: int(v) for k, v in m.saf_diagnosis.value_counts().items()}
            prot[tag]["base"] = samp
        else:
            g1["plasma_unmapped_all_IGHV"] = bool(len(um) > 0 and um.Genes.str.startswith("IGHV").all())
            b = meta[meta.plasma_proteomics_filename.notna()]
            fu = meta[meta.plasma_follow_up_proteomics_filename.notna()]
            base = [c for c in samp if c in set(b.plasma_proteomics_filename)]
            fucols = [c for c in samp if c in set(fu.plasma_follow_up_proteomics_filename)]
            g1["plasma_baseline_cols"] = len(base)
            g1["plasma_baseline_participants"] = int(b[b.plasma_proteomics_filename.isin(samp)].patient_name.nunique())
            g1["plasma_baseline_sample_group"] = {k: int(v) for k, v in b.sample_group.value_counts().items()}
            g1["plasma_followup_cols"] = len(fucols)
            g1["plasma_participants"] = int(pd.concat([b.patient_name, fu.patient_name]).nunique())
            g1["plasma_followup_in_baseline"] = int(fu.patient_name.isin(set(b.patient_name)).sum())
            g1["plasma_cols_unassigned"] = len(set(samp) - set(base) - set(fucols))
            if len(base) != len(b) or b.patient_name.nunique() != len(b):
                raise RuntimeError("plasma: baseline is not one sample per participant")
            prot[tag]["base"] = base

    # ------------------------------------------------------------------ CosMx parsing
    heads = []
    for f in COSMX:
        with gzip.open(f, "rb") as fh:
            line = fh.readline()
        heads.append(dict(file=os.path.basename(f), md5=hashlib.md5(line).hexdigest(),
                          cols=line.decode().rstrip("\r\n").split(",")))
    cols = heads[0]["cols"]
    sysc = [c for c in cols if re.match(r"^SystemControl\d+$", c)]
    neg = [c for c in cols if re.match(r"^(Negative\d+|NegativeAdd)$", c)]
    fovc = [c for c in cols if c in ("fov", "cell_ID")]
    targets = [c for c in cols if c not in fovc and c not in sysc and c not in neg]

    def cosmx_targets(cols_):
        s_ = [c for c in cols_ if not re.match(r"^SystemControl\d+$", c) and not re.match(r"^(Negative\d+|NegativeAdd)$", c)
              and c not in ("fov", "cell_ID")]
        return s_
    tsets = [set(cosmx_targets(hh["cols"])) for hh in heads]
    headers_identical = len({hh["md5"] for hh in heads}) == 1
    panel_targets = set.intersection(*tsets)  # A1 CosMx rule: intersection if headers differ
    tgt_rows = []
    for c in targets:
        if c not in panel_targets:
            continue
        if "/" in c:
            genes = COSMX_EXP.get(c)
            if genes is None:
                tgt_rows.append(dict(probe=c, kind="shared", genes="", ensembl="", route="unexpanded"))
            else:
                ens = [cur.get(g, "") for g in genes]
                tgt_rows.append(dict(probe=c, kind="shared", genes=";".join(genes), ensembl=";".join(ens),
                                     route="fixed_expansion_current_symbol"))
        elif c in cur:
            tgt_rows.append(dict(probe=c, kind="single", genes=c, ensembl=cur[c], route="current_symbol"))
        elif len(prev.get(c, set())) == 1:
            e = next(iter(prev[c]))
            tgt_rows.append(dict(probe=c, kind="single", genes=ens2hgncsym.get(e, ""), ensembl=e, route="prev_symbol"))
        else:
            tgt_rows.append(dict(probe=c, kind="single", genes="", ensembl="", route="no_single_gene"))
    tg = pd.DataFrame(tgt_rows)
    single_ens = set(tg.loc[(tg.kind == "single") & (tg.ensembl != ""), "ensembl"])
    shared_ens = set(e for s in tg.loc[tg.kind == "shared", "ensembl"] for e in s.split(";") if e)
    shared_syms = [g for s in tg.loc[tg.kind == "shared", "genes"] for g in s.split(";") if g]
    g1["cosmx_header_cols"] = len(cols)
    g1["cosmx_fov_cellid"] = len(fovc)
    g1["cosmx_systemcontrol"] = len(sysc)
    g1["cosmx_negative"] = len(neg)
    g1["cosmx_targets"] = len(targets)
    g1["cosmx_single_current"] = int((tg.route == "current_symbol").sum())
    g1["cosmx_single_prev"] = sorted(tg.loc[tg.route == "prev_symbol", "probe"])
    g1["cosmx_single_prev_to"] = sorted(tg.loc[tg.route == "prev_symbol", "genes"])
    g1["cosmx_shared_probes"] = int((tg.kind == "shared").sum())
    g1["cosmx_shared_genes"] = len(shared_syms)
    g1["cosmx_shared_genes_in_catalog"] = int(sum(cur.get(g) in cens for g in shared_syms))
    g1["cosmx_shared_genes_not_in_catalog"] = sorted(g for g in shared_syms if cur.get(g) not in cens)
    g1["cosmx_no_gene"] = sorted(tg.loc[tg.route == "no_single_gene", "probe"])
    g1["cosmx_header_md5_all_four"] = heads[0]["md5"] if headers_identical else "DIFFERENT"

    # ------------------------------------------------------------------ single-cell inputs
    pbm, pb_rows0, pb_meta = {}, None, {}
    rows_identical, all_int, all_in_cc, cols_289 = True, True, True, 0
    cc_runs = set(cc["sample"])
    for L in LINEAGES:
        m = pd.read_csv(PB[L], index_col=0)
        m.index = m.index.astype(str)
        if pb_rows0 is None:
            pb_rows0 = list(m.index)
        else:
            rows_identical &= list(m.index) == pb_rows0
        v = m.to_numpy(dtype=np.float64)
        all_int &= bool(np.all(np.abs(v - np.round(v)) < 1e-6))
        all_in_cc &= bool(all(c in cc_runs for c in m.columns))
        cols_289 += int(sum(("GSE289173" in str(run2ds.get(c, ""))) or ("GSE289173" in c) for c in m.columns))
        pbm[L] = m
    g1["sc_all_integer"] = all_int
    g1["sc_all_cols_in_cellcount"] = all_in_cc
    g1["sc_cols_gse289173"] = cols_289
    if cols_289 > 0:
        s6 = True
        g1["sc_gse289173_present"] = True
    g1["sc_runs"] = int(cc["sample"].nunique())
    g1["sc_donors"] = int(runmap.donor.nunique())
    g1["sc_donors_by_dataset"] = {k: int(v) for k, v in runmap.groupby("dataset").donor.nunique().items()}

    # per-lineage donor tables (cells only from runs that are columns of that lineage's pseudobulk)
    lin_info, g4a_rows = {}, []
    for L in LINEAGES:
        m = pbm[L]
        runs = list(m.columns)
        donors = [run2donor.get(r, r) for r in runs]
        dsum = m.T.groupby(np.array(donors)).sum().T  # genes x donors, raw counts summed over runs
        lib = dsum.sum(axis=0)
        ccl = cc.set_index("sample")[f"n_{L}"]
        cells = pd.Series([ccl.get(r, 0) for r in runs], index=runs).groupby(np.array(donors)).sum()
        cells = cells.reindex(dsum.columns)
        dsets = pd.Series({d: run2ds.get(r) for r, d in zip(runs, donors)}).reindex(dsum.columns)
        elig = (cells >= SC_MIN_CELLS) & (lib >= SC_MIN_LIB)
        elig50 = (cells >= 50) & (lib >= SC_MIN_LIB)
        # G4a: pseudobulk columns == runs with >= 10 cells of the lineage
        exp_cols = set(cc.loc[cc[f"n_{L}"] >= SC_COL_MIN_CELLS, "sample"])
        g4a_diff = sorted(set(runs) ^ exp_cols)
        rho = float(spearmanr(cells.values, lib.values).correlation)
        lin_info[L] = dict(dsum=dsum, lib=lib, cells=cells, dsets=dsets, elig=elig, elig50=elig50,
                           n_cols=len(runs), n_donors=dsum.shape[1], n_elig=int(elig.sum()), n_elig50=int(elig50.sum()),
                           g4a_diff=g4a_diff, rho=rho, cc_total=int(cc[f"n_{L}"].sum()))
        g4a_rows.append(dict(lineage=L, pb_cols=len(runs), runs_ge10_cells=len(exp_cols), symmetric_difference=len(g4a_diff),
                             g4_spearman_cells_vs_lib=rho, donors_with_col=dsum.shape[1], eligible=int(elig.sum()),
                             eligible_50cells=int(elig50.sum()), min_lib_eligible=float(lib[elig].min()) if elig.any() else np.nan))
    del pbm
    g1["sc_pb_cols"] = {L: lin_info[L]["n_cols"] for L in LINEAGES}
    g1["sc_pb_donors"] = {L: lin_info[L]["n_donors"] for L in LINEAGES}
    g1["sc_eligible"] = {L: lin_info[L]["n_elig"] for L in LINEAGES}
    g1["sc_eligible_50cells"] = {L: lin_info[L]["n_elig50"] for L in LINEAGES}
    used = [L for L in LINEAGES if lin_info[L]["cc_total"] > 0 and lin_info[L]["n_elig"] >= SC_MIN_ELIG]
    g1["sc_lineages_used"] = len(used)

    # pseudobulk row mapping (priority: ensembl > current symbol > unique prev symbol > Catalog symbol of no-HGNC gene)
    PRI = {"ensembl": 1, "hgnc_current": 2, "hgnc_prev_unique": 3, "catalog_symbol_no_hgnc": 4, "unmapped": 9}
    rmap = []
    for i, s in enumerate(pb_rows0):
        mm = re.match(r"^(ENSG\d+)(\.\d+)?$", s)
        if mm:
            e, r = mm.group(1), "ensembl"
        elif s in cur:
            e, r = cur[s], "hgnc_current"
        elif len(prev.get(s, set())) == 1:
            e, r = next(iter(prev[s])), "hgnc_prev_unique"
        elif s in nohgnc_sym2ens:
            e, r = nohgnc_sym2ens[s], "catalog_symbol_no_hgnc"
        else:
            e, r = None, "unmapped"
        rmap.append((i, s, e, r))
    rmap = pd.DataFrame(rmap, columns=["row_idx", "row", "ens", "route"])
    rmap["in_catalog"] = rmap.ens.isin(cens)
    rmap["priority"] = rmap.route.map(PRI)
    mc = rmap[rmap.in_catalog].copy()
    best = mc.groupby("ens").priority.transform("min")
    mc["used"] = mc.priority == best
    rmap = rmap.merge(mc[["row_idx", "used"]], on="row_idx", how="left")
    rmap["used"] = rmap["used"].eq(True)
    multi = mc.ens.value_counts()
    multi_genes = multi[multi > 1].index
    g1["pb_rows"] = len(pb_rows0)
    g1["pb_rows_identical_9"] = rows_identical
    g1["pb_rows_ensembl"] = int((rmap.route == "ensembl").sum())
    g1["pb_rows_symbol"] = int((rmap.route != "ensembl").sum())
    g1["pb_catalog_genes_with_row"] = int(mc.ens.nunique())
    for k, r in [("ensembl", "ensembl"), ("current", "hgnc_current"), ("prev", "hgnc_prev_unique"), ("catsym", "catalog_symbol_no_hgnc")]:
        g1[f"pb_rows_route_{k}_in_catalog"] = int(((mc.route == r)).sum())
    g1["pb_catalog_genes_multi_row"] = len(multi_genes)
    sym_only = set()
    for s in rmap.loc[rmap.route != "ensembl", "row"]:
        if s in cur:
            e = {cur[s]}
        else:
            pp = prev.get(s, set())
            e = pp if len(pp) == 1 else set()
        sym_only |= (e & cens)
    g1["pb_catalog_genes_symbol_only_rule"] = len(sym_only)
    g1["pb_catalog_genes_no_row"] = len(cens) - int(mc.ens.nunique())

    # known-answer and OR sets
    ka = c2[c2.symbol.isin(KA_ALL)]
    g1["ka_in_catalog"] = len(ka)
    g1["ka_classes"] = {k: int(v) for k, v in ka.cls.value_counts().items()}
    g1["ka_established"] = sorted(ka.loc[ka.cls == EST, "symbol"])
    g1["ka_untested"] = sorted(ka.loc[ka.cls == "untested", "symbol"])
    g1["ka_sc_symbol_row"] = int(sum(g in set(pb_rows0) for g, _ in KA_SC))
    org = c2[c2.symbol.str.match(OR_PATTERN, na=False)]
    g1["or_genes"] = len(org)
    g1["or_classes"] = {k: int(v) for k, v in org.cls.value_counts().items()}
    g1["or_genes_with_pb_row"] = int(org.ensembl_id.isin(set(mc.ens)).sum())
    g1["or_genes_short_pattern"] = int(c2.symbol.str.match(OR_PATTERN_SHORT, na=False).sum())

    # NC-2 pool
    sec_dec, edges = pd.qcut(sec.AveExpr, 10, labels=False, retbins=True)
    sec["dec"] = sec_dec.astype(int)
    need = sec.groupby(["dec", "cls"]).size().unstack(fill_value=0)
    pool = c2[c2.cls.isin(["unresolved", "untested"]) & (c2.pc == 1) & c2.AveExpr.notna() & ~c2.ensembl_id.isin(set(sec.ensembl_id))].copy()
    g1["nc2_pool_total"] = len(pool)
    pool = pool[(pool.AveExpr >= edges[0]) & (pool.AveExpr <= edges[-1])].copy()
    g1["nc2_pool_inside_range"] = len(pool)
    pool["dec"] = pd.cut(pool.AveExpr, bins=edges, labels=False, include_lowest=True).astype(int)
    pool_n = pool.dec.value_counts().sort_index()
    g1["nc2_pool_min_per_decile"] = int(pool_n.min())
    g1["nc2_need_max_per_decile"] = int(need.sum(axis=1).max())  # pseudo-anchored + pseudo-established per decile

    # ------------------------------------------------------------------ G1 compare
    g1_fail = []
    for k, exp in G1.items():
        obs = g1.get(k, "MISSING")
        ok = (obs == exp)
        if not ok:
            g1_fail.append(k)
        D.add(2, f"G1 {k}", obs, exp, ok)
    extra = sorted(set(g1) - set(G1))
    if extra:
        raise RuntimeError(f"internal: G1 keys without target {extra}")
    json.dump(dict(observed=g1, expected=G1, failed=g1_fail), open(os.path.join(out, "g1_counts.json"), "w"), indent=1, default=str)

    # ------------------------------------------------------------------ G2 PC-A gate reproduction (gate's symbol rule, copied)
    def symset(path, col="Genes"):
        t_ = pd.read_csv(path, sep="\t", usecols=[col])
        s_ = set()
        for v_ in t_[col].dropna().astype(str):
            s_.update(x.strip() for x in v_.split(";"))
        return s_

    def symset_frac(path, thr=0.5):
        t_ = pd.read_csv(path, sep="\t")
        q_ = t_.iloc[:, 3:]
        frac_ = q_.notna().mean(axis=1)
        s_ = set()
        for v_, f_ in zip(t_["Genes"].astype(str), frac_):
            if f_ >= thr:
                s_.update(x.strip() for x in v_.split(";"))
        return s_, q_.shape[1]
    gl = symset(P["Liver proteome"])
    gp = symset(P["Plasma proteome (unfiltered)"])
    gl50, _ = symset_frac(P["Liver proteome"])
    gp50, _ = symset_frac(P["Plasma proteome (unfiltered)"])
    gcos = {c for c in cols if not any(k in c for k in ("NegPrb", "FalseCode", "SystemControl", "fov", "cell_ID", "cell"))}
    rep = cat[["ensembl_id", "symbol", "cls"]].copy()
    rep["cosmx"] = rep.symbol.isin(gcos).astype(int)
    rep["liver_prot"] = rep.symbol.isin(gl).astype(int)
    rep["liver_prot50"] = rep.symbol.isin(gl50).astype(int)
    rep["plasma_prot"] = rep.symbol.isin(gp).astype(int)
    rep["plasma_prot50"] = rep.symbol.isin(gp50).astype(int)
    gate = pd.read_csv(P["Gate per-gene file (PC-A reproduction only)"], sep="\t")
    gm = rep.merge(gate[["ensembl_id", "cosmx", "liver_prot", "liver_prot50", "plasma_prot", "plasma_prot50", "sc_lineage"]],
                   on="ensembl_id", how="outer", suffixes=("", "_gate"), indicator=True)
    g2_ok = bool((gm._merge == "both").all()) and len(gm) == 28233
    D.add(2, "G2 gate file rows joined on ensembl_id", int((gm._merge == "both").sum()), 28233, g2_ok)
    for col in ["liver_prot", "liver_prot50", "plasma_prot", "plasma_prot50", "cosmx"]:
        nd = int((gm[col] != gm[col + "_gate"]).sum())
        D.add(2, f"G2 differences {col}", nd, 0, nd == 0)
        g2_ok &= nd == 0
    for col, ea, ee in [("liver_prot", 224, 226), ("plasma_prot", 15, 50), ("cosmx", 25, 157)]:
        sa, se_ = int(rep.loc[rep.cls == ANCH, col].sum()), int(rep.loc[rep.cls == EST, col].sum())
        ok = (sa, se_) == (ea, ee)
        D.add(2, f"G2 class sums {col} anchored/established", f"{sa}/{se_}", f"{ea}/{ee}", ok)
        g2_ok &= ok

    # ------------------------------------------------------------------ D1 to D5 (mapping, panel, single-cell donors)
    for tag in ["liver", "plasma"]:
        mp = prot[tag]["map"]
        D.add(3, f"D1 {tag} protein groups by route", {k: int(v) for k, v in mp.route.value_counts().items()})
        D.add(3, f"D1 {tag} groups with >1 HGNC gene", int((mp.n_hgnc_ens > 1).sum()))
        D.add(3, f"D1 {tag} dropped ambiguous previous symbols", prot[tag]["amb"].symbol.tolist())
        mp.to_csv(os.path.join(out, f"d1_protein_group_map_{tag}.tsv"), sep="\t", index=False)
        prot[tag]["amb"].to_csv(os.path.join(out, f"d1_dropped_ambiguous_symbols_{tag}.tsv"), sep="\t", index=False)
        d2 = mp[mp.d2_catalog_symbol_ens != ""]
        d2g = set(";".join(d2.d2_catalog_symbol_ens).split(";")) - {""}
        D.add(4, f"D2 {tag}: Catalog genes without HGNC Ensembl matched by Catalog symbol", len(d2g),
              note=";".join(sorted(cat.set_index("ensembl_id").loc[sorted(d2g), "symbol"])) if d2g else "")
    D.add(3, "D1 pseudobulk rows by route (all rows)", {k: int(v) for k, v in rmap.route.value_counts().items()})
    D.add(3, "D1 pseudobulk rows into Catalog by route", {k: int(v) for k, v in mc.route.value_counts().items()})
    D.add(3, "D1 Catalog genes reached by more than one row", len(multi_genes))
    rmap.to_csv(os.path.join(out, "d1_pseudobulk_row_map.tsv.gz"), sep="\t", index=False)
    mcm = mc[mc.ens.isin(multi_genes)].sort_values(["ens", "priority"])
    mcm.to_csv(os.path.join(out, "d1_pseudobulk_multi_row_genes.tsv"), sep="\t", index=False)
    D.add(4, "D2 single-cell: Catalog genes reached only through the Catalog symbol (no HGNC Ensembl)",
          int(mc[mc.used & (mc.route == "catalog_symbol_no_hgnc")].ens.nunique()))
    d2cos = sorted(set(targets) & set(nohgnc_sym2ens))
    D.add(4, "D2 CosMx: target names equal to a no-HGNC Catalog symbol (not used; CosMx maps by HGNC only)", len(d2cos), note=";".join(d2cos))
    tg.to_csv(os.path.join(out, "d3_cosmx_targets.tsv"), sep="\t", index=False)
    D.add(5, "D3 CosMx targets by route", {k: int(v) for k, v in tg.route.value_counts().items()})
    D.add(5, "D3 CosMx target table rows", len(tg), 968, len(tg) == 968)
    tg_nocat = tg[(tg.ensembl != "") & tg.ensembl.apply(lambda s: not any(e in cens for e in s.split(";") if e))]
    D.add(5, "D3 CosMx probes whose HGNC Ensembl IDs reach no Catalog key (reported; the CosMx rule has no Catalog-symbol route)",
          len(tg_nocat), note=";".join(tg_nocat.probe))
    pd.DataFrame([{k: v for k, v in hh.items() if k != "cols"} | {"n_cols": len(hh["cols"])} for hh in heads]) \
        .to_csv(os.path.join(out, "d4_cosmx_headers.tsv"), sep="\t", index=False)
    D.add(6, "D4 CosMx header md5 (Leuven_1..4)", [hh["md5"] for hh in heads])
    D.add(6, "D4 CosMx headers identical across the four exprMat files", headers_identical, True, headers_identical,
          note="" if headers_identical else f"panel = intersection; differing targets: {sorted(set.union(*tsets) - panel_targets)}")
    g4 = pd.DataFrame(g4a_rows)
    g4.to_csv(os.path.join(out, "d5_sc_lineage_donors_g4.tsv"), sep="\t", index=False)
    s5 = []
    for r in g4.itertuples():
        D.add(7, f"G4a {r.lineage}: pseudobulk columns vs runs with >=10 cells, symmetric difference", r.symmetric_difference, 0, r.symmetric_difference == 0)
        D.add(7, f"G4 {r.lineage}: Spearman donor cells vs donor library size", round(r.g4_spearman_cells_vs_lib, 4), ">=0.5",
              bool(r.g4_spearman_cells_vs_lib >= G4_MIN_RHO))
        if r.lineage in used and (r.symmetric_difference != 0 or not r.g4_spearman_cells_vs_lib >= G4_MIN_RHO):
            s5.append(r.lineage)
    ebd = []
    for L in LINEAGES:
        li = lin_info[L]
        e = li["dsets"][li["elig"]].value_counts().sort_index()
        e_all = li["dsets"].value_counts().sort_index()
        for ds in sorted(set(e_all.index)):
            ebd.append(dict(lineage=L, dataset=ds, donors_with_col=int(e_all.get(ds, 0)), eligible=int(e.get(ds, 0)),
                            eligible_50cells=int(li["dsets"][li["elig50"]].value_counts().get(ds, 0))))
        D.add(7, f"D5 {L} eligible donors by dataset", {k: int(v) for k, v in e.items()})
    pd.DataFrame(ebd).to_csv(os.path.join(out, "d5_sc_eligible_by_dataset.tsv"), sep="\t", index=False)
    D.add(7, "D5 lineages used (cell counts present and >=10 eligible donors)", used)
    D.add(7, "D5 lineages in the cell-count table without counts (not used)",
          {x: int(cc[f"n_{x}"].sum()) for x in CC_NOT_USED})
    D.add(7, "D5 donors per dataset (all runs)", g1["sc_donors_by_dataset"])

    # ------------------------------------------------------------------ stops S1/S2/S6
    if s6:
        stops.append(("S6", "GSE289173 appears in a single-cell input"))
    if g1_fail or not g2_ok:
        stops.append(("S2", f"G1 differences: {g1_fail}; G2 pass: {g2_ok}"))
    if stops:
        return finish_stop(a, D, out, stops)
    if a.mode == "gates":
        D.msg("gates mode: S1, S6, G1, G2, G4a, G4 and mapping diagnostics done; no field built")
        D.flush()
        write_result(out, "gates", "gates passed (debug)", "gates only; no field built",
                     [f"S5 lineages flagged: {s5}"])
        return 0

    # ------------------------------------------------------------------ null-debug permutations (never in full mode)
    if a.mode == "null":
        prng = np.random.default_rng(987654321)
        for tag in ["liver", "plasma"]:
            t = prot[tag]["t"]
            q = t[prot[tag]["samp"]].to_numpy()
            t.loc[:, prot[tag]["samp"]] = q[prng.permutation(len(t))]
        perm_rows = prng.permutation(len(pb_rows0))
        okx = c2.AveExpr.notna().to_numpy()  # permute classes among genes with AveExpr so the class models run
        cl = c2.cls.to_numpy().copy()
        cl[okx] = cl[okx][prng.permutation(int(okx.sum()))]
        c2["cls"] = cl
        D.msg("NULL DEBUG: proteome rows, pseudobulk rows and Catalog classes permuted; outputs are not results")
    else:
        perm_rows = None

    # ------------------------------------------------------------------ build fields
    F = c2[["ensembl_id", "symbol"]].copy().set_index("ensembl_id")
    for tag in ["liver", "plasma"]:
        t, mp = prot[tag]["t"], prot[tag]["map"]
        q_all = t[prot[tag]["samp"]].notna()
        q_base = t[prot[tag]["base"]].notna()
        grp = pd.DataFrame({"group": np.arange(len(t)), "any": (q_all.sum(axis=1) >= 1).to_numpy().astype(int),
                            "frac": q_base.mean(axis=1).to_numpy(), "shared": (mp.n_all_ens > 1).astype(int).to_numpy(),
                            "route": mp.route.to_numpy(), "all_ens": mp.all_ens.to_numpy(),
                            "d2": mp.d2_catalog_symbol_ens.to_numpy()})
        ex = grp.assign(ens=grp.all_ens.str.split(";")).explode("ens")
        ex = ex[ex.ens.isin(cens)]
        ex["via_catalog_symbol"] = [e in d.split(";") for e, d in zip(ex.ens, ex.d2)]
        agg = ex.groupby("ens").agg(any_=("any", "max"), frac=("frac", "max"), shared=("shared", "max"),
                                    n_groups=("group", "nunique"), routes=("route", lambda s: ";".join(sorted(set(s)))),
                                    only_catsym=("via_catalog_symbol", "min"))
        F[f"{tag}_proteome_any"] = agg.any_.reindex(F.index).fillna(0).astype(int)
        F[f"{tag}_proteome_frac"] = agg.frac.reindex(F.index).fillna(0.0).round(6)
        F[f"{tag}_proteome_50"] = (F[f"{tag}_proteome_frac"] >= 0.5).astype(int)
        F[f"{tag}_proteome_shared_group"] = agg.shared.reindex(F.index).fillna(0).astype(int)
        F[f"{tag}_proteome_n_groups"] = agg.n_groups.reindex(F.index).fillna(0).astype(int)
        F[f"{tag}_proteome_route"] = agg.routes.reindex(F.index).fillna("not_in_proteome")
        F[f"{tag}_proteome_symbol_only_match"] = agg.only_catsym.reindex(F.index).eq(True).astype(int)
        D.add(4, f"D2 {tag}: Catalog genes whose only match is the Catalog symbol, detected (any)",
              int(((F[f"{tag}_proteome_symbol_only_match"] == 1) & (F[f"{tag}_proteome_any"] == 1)).sum()))
    F["cosmx_panel"] = np.where(F.index.isin(single_ens), "on_panel",
                                np.where(F.index.isin(shared_ens), "on_panel_shared_probe", "not_on_panel"))
    F["cosmx_any"] = (F.cosmx_panel != "not_on_panel").astype(int)

    # single-cell detection per lineage at the donor level
    used_rows = rmap[rmap.used].copy()
    if perm_rows is not None:
        used_rows["row_idx"] = perm_rows[used_rows.row_idx.to_numpy()]
    gene_of = used_rows.ens.to_numpy()
    ridx = used_rows.row_idx.to_numpy()
    fr = {}            # primary fraction per lineage (gene index)
    sens = {k: {} for k in ["cells50", "dataset_balanced"]}
    for L in used:
        li = lin_info[L]
        dsum, lib, elig = li["dsum"], li["lib"], li["elig"]
        X = dsum.loc[:, elig].to_numpy()
        libe = lib[elig].to_numpy()
        det = (X >= SC_MIN_COUNT) & ((X / libe[None, :]) * 1e6 >= SC_MIN_CPM)  # rows x eligible donors
        gdet = pd.DataFrame(det[ridx], index=gene_of).groupby(level=0).max()  # gene detected if any row passes
        donors_e = list(dsum.columns[elig.to_numpy()])
        gdet.columns = donors_e
        fr[L] = gdet.mean(axis=1)
        e50 = li["elig50"][elig].to_numpy()
        sens["cells50"][L] = gdet.loc[:, e50].mean(axis=1)
        ds_e = li["dsets"][elig]
        keep = [ds for ds, n in ds_e.value_counts().items() if n >= 3]
        sens["dataset_balanced"][L] = pd.concat([gdet.loc[:, (ds_e == ds).to_numpy()].mean(axis=1) for ds in sorted(keep)], axis=1).mean(axis=1)
        D.add(8, f"D6 {L}: datasets with >=3 eligible donors used in the dataset-balanced fraction", sorted(keep))
    FR = pd.DataFrame(fr).reindex(columns=used)
    has_row = F.index.isin(FR.index)
    FRg = FR.reindex(F.index)
    det_lin = FRg >= SC_DONOR_FRAC
    F["sc_detected"] = np.where(has_row, det_lin.any(axis=1).astype(int), np.nan)
    F["sc_lineages"] = np.where(has_row, det_lin.apply(lambda r: ";".join([L for L in used if r[L]]), axis=1), None)
    bestL = FRg.fillna(-1).to_numpy().argmax(axis=1)
    F["sc_best_lineage"] = np.where(has_row, np.array(used)[bestL], None)
    F["sc_best_donor_frac"] = np.where(has_row, FRg.max(axis=1).round(6), np.nan)
    F["sc_na_reason"] = np.where(has_row, "", "no pseudobulk row under the row-mapping rule")
    for L in used:
        F[f"sc_frac_{L}"] = FRg[L].round(6)

    # D6 sensitivity (reported only)
    d6 = []
    prim = F.sc_detected
    for name, thr, src in [("primary", 0.5, fr), ("cells50_floor", 0.5, sens["cells50"]), ("donor_frac_0.25", 0.25, fr),
                           ("donor_frac_0.75", 0.75, fr), ("dataset_balanced", 0.5, sens["dataset_balanced"])]:
        M = pd.DataFrame(src).reindex(columns=used).reindex(F.index)
        dd = np.where(has_row, (M >= thr).any(axis=1).astype(int), np.nan)
        ka_pass = []
        for g, L in KA_SC:
            e = sym2cat[g]
            ka_pass.append(bool(M.loc[e, L] >= thr) if L in used else False)
        d6.append(dict(variant=name, n_detected=int(np.nansum(dd)), n_measured=int(has_row.sum()),
                       agree_with_primary=int(np.nansum(dd == prim)), ka_pass=int(sum(ka_pass)),
                       ka_required_pass=bool(ka_pass[0] and ka_pass[1])))
    d6 = pd.DataFrame(d6)
    d6.to_csv(os.path.join(out, "d6_sc_sensitivity.tsv"), sep="\t", index=False)
    for r in d6.itertuples():
        D.add(8, f"D6 {r.variant}", dict(n_detected=r.n_detected, n_measured=r.n_measured, agree=r.agree_with_primary,
                                          ka_pass_of_11=r.ka_pass, ka_required=r.ka_required_pass))

    # ------------------------------------------------------------------ G3 / A4 coverage
    a4 = {}
    for fld, cols_ in [("liver_proteome", ["liver_proteome_any", "liver_proteome_frac", "liver_proteome_50"]),
                       ("plasma_proteome", ["plasma_proteome_any", "plasma_proteome_frac", "plasma_proteome_50"]),
                       ("cosmx_panel", ["cosmx_panel"])]:
        ok = all(F[c].notna().all() for c in cols_) and len(F) == 28233
        if fld == "cosmx_panel":
            ok &= set(F.cosmx_panel.unique()) <= {"on_panel", "on_panel_shared_probe", "not_on_panel"}
        a4[fld] = bool(ok)
    sc_na = F.sc_detected.isna()
    a4["sc_lineage"] = bool(len(F) == 28233 and F.loc[~sc_na, "sc_detected"].isin([0, 1]).all()
                            and (F.loc[sc_na, "sc_na_reason"] != "").all() and int(sc_na.sum()) == G1["pb_catalog_genes_no_row"])
    for k, v in a4.items():
        D.add(9, f"A4/G3 complete coverage {k}", v, True, v)
    D.add(9, "A4/G3 sc_detected NA (no pseudobulk row)", int(sc_na.sum()), G1["pb_catalog_genes_no_row"])

    # ------------------------------------------------------------------ A1 known answers
    sym = F.reset_index().set_index("symbol")
    a1_rows = []
    for g in KA_LIVER + KA_LIVER_REPORT:
        a1_rows.append(dict(field="liver_proteome", gene=g, lineage="", required=g in KA_LIVER,
                            value=int(sym.loc[g, "liver_proteome_any"]), frac=float(sym.loc[g, "liver_proteome_frac"])))
    for g in KA_PLASMA:
        a1_rows.append(dict(field="plasma_proteome", gene=g, lineage="", required=True,
                            value=int(sym.loc[g, "plasma_proteome_any"]), frac=float(sym.loc[g, "plasma_proteome_frac"]),
                            plasma_50=int(sym.loc[g, "plasma_proteome_50"])))
    for g, L in KA_SC:
        v = float(sym.loc[g, f"sc_frac_{L}"]) if L in used else np.nan
        a1_rows.append(dict(field="sc_lineage", gene=g, lineage=L, required=(g, L) in KA_SC_REQUIRED,
                            value=int(v >= SC_DONOR_FRAC) if np.isfinite(v) else 0, frac=v))
    a1 = pd.DataFrame(a1_rows)
    a1.to_csv(os.path.join(out, "a1_known_answer.tsv"), sep="\t", index=False)
    A1 = {}
    lv = a1[(a1.field == "liver_proteome") & a1.required]
    A1["liver_proteome"] = bool((lv.value == 1).all())
    D.add(9, "A1 liver proteome known answers detected (rule: all 7)", f"{int(lv.value.sum())}/7", "7/7", A1["liver_proteome"],
          note=f"TM6SF2 (reported, not required) = {int(sym.loc['TM6SF2', 'liver_proteome_any'])}")
    pl = a1[a1.field == "plasma_proteome"]
    A1["plasma_proteome"] = bool((pl.value == 1).all())
    D.add(9, "A1 plasma proteome known answers detected (rule: all 8)", f"{int(pl.value.sum())}/8", "8/8", A1["plasma_proteome"],
          note=f"plasma_proteome_50: {int(pl.plasma_50.sum())}/8 (predicted >=7; reported, not a rule)")
    sc_a1 = a1[a1.field == "sc_lineage"]
    req_ok = bool(sc_a1[sc_a1.required].value.eq(1).all())
    n_ok = int(sc_a1.value.sum())
    A1["sc_lineage"] = req_ok and n_ok >= 9
    D.add(9, "A1 single-cell known answers (rule: KRT19 and EPCAM in Cholangiocytes, and >=9 of 11)", f"{n_ok}/11; required pair {req_ok}",
          ">=9/11 and required pair", A1["sc_lineage"], note="; ".join(f"{r.gene}@{r.lineage}={r.frac:.3f}" for r in sc_a1.itertuples()))
    neg_in_field = sorted(set(neg) & set(tg.probe))
    exp_genes_ok = (g1["cosmx_shared_probes"] == 19 and g1["cosmx_shared_genes"] == 42 and (tg.route != "unexpanded").all())
    a1_cos = (len(panel_targets & set(targets)) == 968 if headers_identical else True) and exp_genes_ok and not neg_in_field
    A1["cosmx_panel"] = bool(a1_cos)
    D.add(9, "A1 CosMx structural (968 named targets, 19 shared probes to 42 genes, no Negative probe in the field, D4)",
          dict(targets=len(tg), shared=g1["cosmx_shared_probes"], shared_genes=g1["cosmx_shared_genes"],
               negative_in_field=len(neg_in_field), headers_identical=headers_identical,
               catalog_on_panel=int((F.cosmx_panel == "on_panel").sum()),
               catalog_on_shared_probe=int((F.cosmx_panel == "on_panel_shared_probe").sum())),
          ok=A1["cosmx_panel"])

    # ------------------------------------------------------------------ A2 / NC-1 olfactory receptors
    orF = F[F.symbol.str.match(OR_PATTERN, na=False)]
    A2, a2_rows = {}, []
    for fld, col, denom in [("liver_proteome", "liver_proteome_any", len(orF)), ("plasma_proteome", "plasma_proteome_any", len(orF)),
                            ("cosmx_panel", "cosmx_any", len(orF)), ("sc_lineage", "sc_detected", int(orF.sc_detected.notna().sum()))]:
        n = int(np.nansum(orF[col]))
        ok = n <= A2_MAX_FRAC * denom
        A2[fld] = bool(ok)
        a2_rows.append(dict(field=fld, or_detected=n, measurable=denom, limit=A2_MAX_FRAC * denom, passed=ok,
                            genes=";".join(sorted(orF.loc[orF[col] == 1, "symbol"]))))
        D.add(9, f"A2/NC-1 olfactory receptor genes called present: {fld}", f"{n}/{denom}", f"<= {A2_MAX_FRAC * denom:.2f}", ok,
              note=";".join(sorted(orF.loc[orF[col] == 1, "symbol"]))[:300])
    pd.DataFrame(a2_rows).to_csv(os.path.join(out, "a2_or_negative_control.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------------ A3 abundance control
    u = c2[(c2.cls == "unresolved") & (c2.pc == 1) & c2.AveExpr.notna()][["ensembl_id", "AveExpr"]].copy()
    D.add(9, "A3 unresolved protein-coding genes with AveExpr", len(u), 12679, len(u) == 12679)
    u["dec"] = pd.qcut(u.AveExpr, 10, labels=False).astype(int)
    u = u.join(F, on="ensembl_id")
    A3, a3_rows = {}, []
    for fld, col in [("liver_proteome", "liver_proteome_any"), ("plasma_proteome", "plasma_proteome_any"),
                     ("sc_lineage", "sc_detected"), ("liver_proteome_50 (report)", "liver_proteome_50"),
                     ("plasma_proteome_50 (report)", "plasma_proteome_50"), ("cosmx_panel (not applicable)", "cosmx_any")]:
        uu = u[u[col].notna()]
        rate = uu.groupby("dec")[col].mean().reindex(range(10))
        n = uu.groupby("dec")[col].size().reindex(range(10))
        diff = float(rate[9] - rate[0])
        rho = float(spearmanr(np.arange(10), rate.values).correlation)
        ok = diff >= A3_MIN_DIFF and rho >= A3_MIN_RHO
        if fld in ("liver_proteome", "plasma_proteome", "sc_lineage"):
            A3[fld] = bool(ok)
        for d_ in range(10):
            a3_rows.append(dict(field=fld, decile=d_, n=int(n[d_]), rate=float(rate[d_])))
        D.add(9, f"A3 abundance control {fld}: decile 9 minus decile 0; Spearman", f"{diff:.3f}; {rho:.3f} (n={len(uu)})",
              ">=0.30; >=0.80", ok if fld in A3 else None)
    pd.DataFrame(a3_rows).to_csv(os.path.join(out, "a3_abundance_by_decile.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------------ ship decisions
    ship, why = {}, {}
    for fld in ["liver_proteome", "plasma_proteome", "cosmx_panel", "sc_lineage"]:
        r = []
        if not A1[fld]:
            r.append("S3: failed A1")
        if not a4[fld]:
            r.append("S3: failed A4")
        if not A2[fld]:
            r.append("S4: failed A2")
        if fld == "sc_lineage" and s5:
            r.append(f"S5: G4/G4a failed in {s5}")
        ship[fld] = not r
        note = [] if fld == "cosmx_panel" else ([] if A3[fld] else ["A3 failed (shipped with the failure noted)"])
        why[fld] = r + note
        D.add(9, f"ship {fld}", ship[fld], note="; ".join(why[fld]))
    pd.DataFrame([dict(field=k, A1=A1[k], A2=A2[k], A3=A3.get(k, "not applicable"), A4=a4[k], shipped=ship[k],
                       reasons="; ".join(why[k])) for k in ship]).to_csv(os.path.join(out, "acceptance.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------------ write per-gene fields
    F_out = F.reset_index()
    F_out.to_csv(os.path.join(out, "h66_fields_all_including_unshipped.tsv.gz"), sep="\t", index=False)
    keep = ["ensembl_id", "symbol"]
    lim = pd.Series([[] for _ in range(len(F_out))])
    if ship["liver_proteome"]:
        keep += ["liver_proteome_any", "liver_proteome_frac", "liver_proteome_50", "liver_proteome_shared_group"]
    if ship["plasma_proteome"]:
        keep += ["plasma_proteome_any", "plasma_proteome_frac", "plasma_proteome_50", "plasma_proteome_shared_group"]
    if ship["cosmx_panel"]:
        keep += ["cosmx_panel"]
    if ship["sc_lineage"]:
        keep += ["sc_detected", "sc_lineages", "sc_best_lineage", "sc_best_donor_frac"]
    shipped = F_out[keep].copy()
    limtxt = []
    for r in F_out.itertuples():
        s_ = []
        if ship["liver_proteome"] and r.liver_proteome_any == 0:
            s_.append("not observed in PXD051911 liver proteome; absence is untestable, not negative")
        if ship["plasma_proteome"] and r.plasma_proteome_any == 0:
            s_.append("not observed in PXD051911 plasma proteome; absence is untestable, not negative")
        if ship["cosmx_panel"] and r.cosmx_panel == "not_on_panel":
            s_.append("not on the Govaere CosMx panel; absence is untestable, not negative")
        if ship["sc_lineage"] and r.sc_detected == 0:
            s_.append("not observed in single-cell lineage pseudobulk; absence is untestable, not negative")
        limtxt.append("; ".join(s_))
    shipped["h66_limitation"] = limtxt
    shipped.to_csv(os.path.join(out, "h66_catalog_fields_shipped.tsv.gz"), sep="\t", index=False)

    # D1 changes vs the gate's symbol rule
    chg = []
    gmi = gate.set_index("ensembl_id")
    for new, old in [("liver_proteome_any", "liver_prot"), ("liver_proteome_50", "liver_prot50"),
                     ("plasma_proteome_any", "plasma_prot"), ("plasma_proteome_50", "plasma_prot50"),
                     ("cosmx_any", "cosmx"), ("sc_detected", "sc_lineage")]:
        o = gmi[old].reindex(F.index)
        n = F[new]
        dmask = (o != n) & ~(n.isna())
        dna = n.isna() & (o == 1)
        for e in F.index[dmask | dna]:
            chg.append(dict(ensembl_id=e, symbol=F.at[e, "symbol"], field=new, gate_field=old, old=int(o[e]),
                            new=("NA" if pd.isna(n[e]) else int(n[e]))))
    chg = pd.DataFrame(chg, columns=["ensembl_id", "symbol", "field", "gate_field", "old", "new"])
    chg.to_csv(os.path.join(out, "d1_changes_vs_gate.tsv"), sep="\t", index=False)
    if len(chg):
        D.add(3, "D1 genes whose call changes from the gate's symbol rule (field: old->new counts)",
              {f"{k[0]} {k[1]}->{k[2]}": int(v) for k, v in chg.groupby(["field", "old", "new"]).size().items()})
    mt = chg[chg.symbol == "MTARC1"]
    D.add(3, "D1 MTARC1 changes", mt[["field", "old", "new"]].to_dict("records"))

    # ------------------------------------------------------------------ secondary description (class; S7 applies)
    ss = np.random.SeedSequence(a.seed)
    child = ss.spawn(2)
    rng_nc2, rng_boot = np.random.default_rng(child[0]), np.random.default_rng(child[1])
    FIELDS = [("liver_proteome_any", "liver_proteome"), ("plasma_proteome_any", "plasma_proteome"),
              ("cosmx_any", "cosmx_panel"), ("sc_detected", "sc_lineage")]
    Fd = F[[f for f, _ in FIELDS]]

    # NC-2 abundance-matched pseudo-classes
    pool = pool.join(Fd, on="ensembl_id").sort_values("ensembl_id")
    pools = {d_: pool[pool.dec == d_].reset_index(drop=True) for d_ in range(10)}
    na_, ne_ = need.get(ANCH, pd.Series(0, index=range(10))), need.get(EST, pd.Series(0, index=range(10)))
    short = {d_: len(pools[d_]) < int(na_[d_] + ne_[d_]) for d_ in range(10)}
    D.add(10, "NC-2 deciles drawn with replacement (pool short)", [d_ for d_, v in short.items() if v])
    nc2 = []
    t0 = time.time()
    for b in range(a.nc2_draws):
        parts = []
        for d_ in range(10):
            pd_ = pools[d_]
            k_a, k_e = int(na_[d_]), int(ne_[d_])
            idx = rng_nc2.choice(len(pd_), size=k_a + k_e, replace=short[d_])
            s_ = pd_.iloc[idx].copy()
            s_["anch"] = [1] * k_a + [0] * k_e
            parts.append(s_)
        dr = pd.concat(parts, ignore_index=True)
        for col, fld in FIELDS:
            dd = dr[dr[col].notna()]
            r = fit_logit(dd.assign(y=dd[col].astype(int)), "y ~ anch + bs(AveExpr, df=4)")
            nc2.append(dict(draw=b, field=fld, ok=r["ok"], OR=r.get("OR", np.nan), p=r.get("p", np.nan),
                            reason=r.get("reason", "")))
        if b % 100 == 0:
            D.msg(f"NC-2 draw {b} ({time.time() - t0:.0f} s)")
    nc2 = pd.DataFrame(nc2)
    nc2.to_csv(os.path.join(out, "nc2_draws.tsv.gz"), sep="\t", index=False)
    NC2, nc2s = {}, []
    for _, fld in FIELDS:
        x = nc2[(nc2.field == fld) & nc2.ok]
        med = float(np.median(x.OR)) if len(x) else np.nan
        fp = float((x.p < 0.05).mean()) if len(x) else np.nan
        ok = bool(len(x) > 0 and NC2_OR_LO <= med <= NC2_OR_HI and NC2_P_LO <= fp <= NC2_P_HI)
        NC2[fld] = ok
        nc2s.append(dict(field=fld, draws=a.nc2_draws, valid=len(x), excluded=a.nc2_draws - len(x), median_OR=med,
                         frac_p_lt_0_05=fp, passed=ok))
        D.add(10, f"NC-2 {fld}: valid draws; median OR; fraction p<0.05", f"{len(x)}; {med:.3f}; {fp:.3f}",
              "median 0.85-1.15; fraction 0.02-0.10", ok)
    pd.DataFrame(nc2s).to_csv(os.path.join(out, "nc2_summary.tsv"), sep="\t", index=False)

    # coverage by class (all five classes)
    cls_ = c2.set_index("ensembl_id").cls.reindex(F.index)
    cov = []
    for c_ in ["genetically_anchored", "concordant", "established_state_associated", "unresolved", "untested"]:
        sub = F[cls_ == c_]
        for col in ["liver_proteome_any", "liver_proteome_50", "plasma_proteome_any", "plasma_proteome_50", "cosmx_any", "sc_detected"]:
            v = sub[col].dropna()
            cov.append(dict(cls=c_, field=col, n_class=len(sub), n_measured=len(v), n_detected=int(v.sum()),
                            pct=100 * float(v.mean()) if len(v) else np.nan))
    cov = pd.DataFrame(cov)
    cov.to_csv(os.path.join(out, "class_coverage_S7_restricted.tsv"), sep="\t", index=False)
    for r in cov.itertuples():
        D.add(11, f"coverage {r.cls} {r.field}", f"{r.n_detected}/{r.n_measured} ({r.pct:.1f}%)")

    # secondary models
    S = c2[c2.cls.isin([ANCH, EST])][["ensembl_id", "cls", "AveExpr", "pc"]].copy()
    S["anch"] = (S.cls == ANCH).astype(int)
    S = S.join(Fd, on="ensembl_id")
    q05, q95 = np.quantile(S.loc[S.anch == 1, "AveExpr"], [0.05, 0.95])
    sec_rows = []
    for col, fld in FIELDS:
        dd = S[S[col].notna()].assign(y=lambda d: d[col].astype(int))
        ev = dd.groupby("anch").y.sum()
        nn = dd.groupby("anch").y.size()
        pa, pe = 100 * ev.get(1, 0) / nn.get(1, 1), 100 * ev.get(0, 0) / nn.get(0, 1)
        r0 = fit_logit(dd, "y ~ anch")
        if r0["ok"]:
            A_ = r0["OR"]
        else:
            a_, b_, c_, d__ = ev.get(1, 0), nn.get(1, 0) - ev.get(1, 0), ev.get(0, 0), nn.get(0, 0) - ev.get(0, 0)
            A_ = (a_ * d__) / (b_ * c_) if b_ * c_ > 0 else np.nan
        few = int(min(ev.get(1, 0), ev.get(0, 0))) < MIN_EVENTS
        r1 = fit_logit(dd, "y ~ anch + bs(AveExpr, df=4) + pc") if not few else dict(ok=False, reason=f"a class has fewer than {MIN_EVENTS} detected genes")
        dcs = dd[(dd.AveExpr >= q05) & (dd.AveExpr <= q95)]
        evcs = dcs.groupby("anch").y.sum()
        fewcs = int(min(evcs.get(1, 0), evcs.get(0, 0))) < MIN_EVENTS
        r2 = fit_logit(dcs, "y ~ anch + bs(AveExpr, df=4) + pc") if not fewcs else dict(ok=False, reason="few events")
        sec_rows.append(dict(field=fld, column=col, n_anch=int(nn.get(1, 0)), n_est=int(nn.get(0, 0)),
                             det_anch=int(ev.get(1, 0)), det_est=int(ev.get(0, 0)), pct_anch=pa, pct_est=pe, OR_unadj=A_,
                             model_ok=r1["ok"], OR_adj=r1.get("OR", np.nan), lo=r1.get("lo", np.nan), hi=r1.get("hi", np.nan),
                             p_nominal=r1.get("p", np.nan), model_note=r1.get("reason", ""),
                             cs_n_anch=int((dcs.anch == 1).sum()), cs_n_est=int((dcs.anch == 0).sum()),
                             cs_ok=r2["ok"], cs_OR_adj=r2.get("OR", np.nan), cs_lo=r2.get("lo", np.nan), cs_hi=r2.get("hi", np.nan),
                             cs_p=r2.get("p", np.nan), cs_note=r2.get("reason", "")))
    sm = pd.DataFrame(sec_rows)
    pv = sm.p_nominal.fillna(1.0).to_numpy()  # Holm over the 4 fields; a field with proportions only enters as p = 1
    sm["p_holm4"] = multipletests(pv, method="holm")[1]
    sm.loc[~sm.model_ok, "p_holm4"] = np.nan
    sm["nc2_passed"] = sm.field.map(NC2)
    sm["label"] = np.where(sm.nc2_passed, "descriptive", "descriptive; abundance not controlled by the spline")
    sm.to_csv(os.path.join(out, "secondary_models_S7_restricted.tsv"), sep="\t", index=False)
    for r in sm.itertuples():
        D.add(11, f"secondary {r.field}", f"anch {r.det_anch}/{r.n_anch} ({r.pct_anch:.1f}%), est {r.det_est}/{r.n_est} ({r.pct_est:.1f}%); "
              f"OR unadj {r.OR_unadj:.3f}; adj {r.OR_adj:.3f} ({r.lo:.3f} to {r.hi:.3f}), p {r.p_nominal:.3g}, Holm {r.p_holm4:.3g}; {r.label}",
              note=r.model_note)
        D.add(11, f"common support {r.field} (anchored 5-95% AveExpr {q05:.2f} to {q95:.2f})",
              f"n {r.cs_n_anch}/{r.cs_n_est}; OR adj {r.cs_OR_adj:.3f} ({r.cs_lo:.3f} to {r.cs_hi:.3f}), p {r.cs_p:.3g}", note=r.cs_note)

    # D7 bootstrap
    d7 = []
    for col, fld in FIELDS:
        dd = S[S[col].notna()].sort_values("ensembl_id")
        ya = dd.loc[dd.anch == 1, col].to_numpy(dtype=float)
        ye = dd.loc[dd.anch == 0, col].to_numpy(dtype=float)
        obs = 100 * (ya.mean() - ye.mean())
        bs_ = np.empty(a.boot)
        for b in range(a.boot):
            bs_[b] = 100 * (ya[rng_boot.integers(0, len(ya), len(ya))].mean() - ye[rng_boot.integers(0, len(ye), len(ye))].mean())
        lo, hi = np.percentile(bs_, [2.5, 97.5])
        d7.append(dict(field=fld, diff_pct_points=obs, boot_lo=lo, boot_hi=hi, B=a.boot))
        D.add(12, f"D7 {fld} anchored minus established, % detected", f"{obs:.2f} ({lo:.2f} to {hi:.2f})")
    pd.DataFrame(d7).to_csv(os.path.join(out, "d7_bootstrap_S7_restricted.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------------ reading templates
    ASSAY = {"liver_proteome": ("the PXD051911 liver proteome", "detected"),
             "plasma_proteome": ("the PXD051911 plasma proteome", "detected"),
             "cosmx_panel": ("the Govaere CosMx panel", "on the panel"),
             "sc_lineage": ("single-cell lineage pseudobulk", "detected")}
    sentences = []
    for r in sm.itertuples():
        assay, verb = ASSAY[r.field]
        if not r.model_ok:
            s_ = (f"In {assay}, {r.pct_anch:.1f}% of genetically anchored and {r.pct_est:.1f}% of established-state genes were {verb} "
                  f"(proportions only: {r.model_note}).")
        else:
            s_ = (f"In {assay}, {r.pct_anch:.1f}% of genetically anchored and {r.pct_est:.1f}% of established-state genes were {verb} "
                  f"(OR unadjusted {r.OR_unadj:.2f}; adjusted for bulk RNA abundance and biotype {r.OR_adj:.2f}, 95% CI {r.lo:.2f} to {r.hi:.2f}).")
            if r.lo <= 1 <= r.hi:
                s_ += f" No class difference was detected after adjustment (95% CI {r.lo:.2f} to {r.hi:.2f})."
            la, lb = np.log(r.OR_unadj), np.log(r.OR_adj)
            if np.sign(la) == np.sign(lb) and abs(lb) < 0.5 * abs(la) and r.nc2_passed:
                s_ += " Most of the difference follows their higher RNA abundance."
            if not r.nc2_passed:
                s_ += " [abundance not controlled by the spline: NC-2 failed]"
        if r.field == "cosmx_panel":
            s_ += " Panel membership reflects panel design, not expression."
        sentences.append(dict(field=r.field, sentence=s_))
    pd.DataFrame(sentences).to_csv(os.path.join(out, "secondary_sentences_S7_restricted.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------------ verdict and RESULT
    shipped_f = [k for k, v in ship.items() if v]
    not_f = [k for k, v in ship.items() if not v]
    parts_v = ["Fields shipped: " + (", ".join(shipped_f) if shipped_f else "none")]
    if not_f:
        parts_v.append("not shipped: " + ", ".join(k + " (" + "; ".join(why[k]) + ")" for k in not_f))
    else:
        parts_v.append("not shipped: none")
    a3_noted = [k for k in A3 if (not A3[k]) and ship[k]]
    if a3_noted:
        parts_v.append("A3 failed, shipped with the failure noted: " + ", ".join(a3_noted))
    verdict = "; ".join(parts_v)
    n_liver = int(F.liver_proteome_any.sum()); n_plasma = int(F.plasma_proteome_any.sum())
    n_cos = int(F.cosmx_any.sum()); n_sc = int(np.nansum(F.sc_detected)); n_sc_m = int(F.sc_detected.notna().sum())
    body = [
        "## Gates", "",
        "S1 input sha256: pass (35 inputs). S6 GSE289173: absent. G1: all counts equal the prespec. G2 PC-A: 0 differences.",
        f"G4a: {'pass' if all(len(lin_info[L]['g4a_diff']) == 0 for L in used) else 'FAIL'} for all 9 lineages. "
        f"G4 Spearman (donor cells vs library size): " + ", ".join(f"{L} {lin_info[L]['rho']:.2f}" for L in used) + ".",
        f"S5 lineages flagged: {s5 if s5 else 'none'}.", "",
        "## Acceptance per field", "",
        "| Field | A1 | A2 (OR genes called present) | A3 (d9-d0; rho) | A4 | Shipped | Notes |", "|---|---|---|---|---|---|---|",
    ]
    a2d = pd.DataFrame(a2_rows).set_index("field")
    a3d = pd.DataFrame(a3_rows)
    for k in ["liver_proteome", "plasma_proteome", "cosmx_panel", "sc_lineage"]:
        if k in A3:
            rr = a3d[a3d.field == k].set_index("decile").rate
            a3s = f"{'pass' if A3[k] else 'fail'} ({rr[9] - rr[0]:.3f}; {spearmanr(np.arange(10), rr.values).correlation:.3f})"
        else:
            a3s = "not applicable"
        body.append(f"| {k} | {'pass' if A1[k] else 'FAIL'} | {'pass' if A2[k] else 'FAIL'} ({a2d.at[k, 'or_detected']}/{a2d.at[k, 'measurable']}) | {a3s} | "
                    f"{'pass' if a4[k] else 'FAIL'} | {'yes' if ship[k] else 'no'} | {'; '.join(why[k])} |")
    body += [
        "", "## Field totals (not by class)", "",
        f"Liver proteome detected (any of 58 participants): {n_liver} of 28,233 Catalog genes; detected in at least half: {int(F.liver_proteome_50.sum())}.",
        f"Plasma proteome detected (any of 184 samples): {n_plasma}; detected in at least half of 143 baseline participants: {int(F.plasma_proteome_50.sum())}.",
        f"CosMx panel: {int((F.cosmx_panel == 'on_panel').sum())} genes on a single-gene probe, {int((F.cosmx_panel == 'on_panel_shared_probe').sum())} on a shared probe.",
        f"Single-cell: {n_sc} of {n_sc_m} measurable genes detected in at least one of {len(used)} lineages; {int(F.sc_detected.isna().sum())} genes have no pseudobulk row (NA).",
        "", "Known answers: " + "; ".join(f"{r.field}:{r.gene}{('@' + r.lineage) if r.lineage else ''}={r.value}" for r in a1.itertuples()),
        "",
        "## Secondary class description (S7: not to be printed outside this directory until recomputed on the adopted post-COLOC Catalog)", "",
        "Family: 4 abundance-adjusted ORs (liver_proteome_any, plasma_proteome_any, CosMx on panel, sc_detected), Holm. Descriptive only; genes are not independent.", "",
    ]
    for r in sm.itertuples():
        body.append(f"- {r.field}: anchored {r.det_anch}/{r.n_anch} ({r.pct_anch:.1f}%), established {r.det_est}/{r.n_est} ({r.pct_est:.1f}%); "
                    f"OR unadjusted {r.OR_unadj:.2f}; adjusted {r.OR_adj:.2f} (95% CI {r.lo:.2f} to {r.hi:.2f}), nominal p {r.p_nominal:.3g}, Holm p {r.p_holm4:.3g}; "
                    f"NC-2 {'pass' if r.nc2_passed else 'fail'}; common support OR {r.cs_OR_adj:.2f} ({r.cs_lo:.2f} to {r.cs_hi:.2f}). {r.model_note}")
    body += ["", "NC-2: " + "; ".join(f"{x['field']} median OR {x['median_OR']:.3f}, p<0.05 in {x['frac_p_lt_0_05']:.3f} of {x['valid']} valid draws ({'pass' if x['passed'] else 'fail'})" for x in nc2s),
             "", "D7 bootstrap (anchored minus established, percentage points): " + "; ".join(f"{x['field']} {x['diff_pct_points']:.1f} ({x['boot_lo']:.1f} to {x['boot_hi']:.1f})" for x in d7),
             "", "Allowed sentences (templates filled by code):", ""] + [f"- {x['sentence']}" for x in sentences]
    body += ["", "## Files", "", "h66_catalog_fields_shipped.tsv.gz (candidate Catalog columns, keyed on Ensembl ID); h66_fields_all_including_unshipped.tsv.gz; "
             "acceptance.tsv; a1_known_answer.tsv; a2_or_negative_control.tsv; a3_abundance_by_decile.tsv; diagnostics.tsv; d1_* to d7_*; nc2_*; "
             "*_S7_restricted.tsv (class description)."]
    D.flush()
    name = write_result(out, a.mode, "completed", verdict, body)
    json.dump(dict(verdict=verdict, ship=ship, why=why, A1=A1, A2=A2, A3=A3, A4=a4, NC2=NC2, s5=s5,
                   totals=dict(liver_any=n_liver, plasma_any=n_plasma, cosmx_any=n_cos, sc_detected=n_sc, sc_measured=n_sc_m)),
              open(os.path.join(out, "result_summary.json"), "w"), indent=1, default=str)
    D.msg(f"wrote {name}; verdict: {verdict}")
    return 0


def finish_stop(a, D, out, stops):
    D.flush()
    verdict = "STOPPED by " + "; ".join(f"{r} ({t})" for r, t in stops) + ". Nothing ships."
    write_result(out, a.mode, "stopped_by_rule", verdict, ["No field is shipped (S1, S2 or S6)."])
    json.dump(dict(verdict=verdict, stops=stops), open(os.path.join(out, "result_summary.json"), "w"), indent=1)
    D.msg(verdict)
    return 0


if __name__ == "__main__":
    sys.exit(main())
