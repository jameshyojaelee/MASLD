#!/usr/bin/env python
"""H34: excess Y-negative immune cells in male liver donors versus stage.

Implements Analysis/Hypotheses/H34-<ts>/PRESPEC_H34.md.
Usage: h34_y_loss.py <repo_root> <out_dir>
"""
import gzip
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

SEED = 34
REPO = Path(sys.argv[1])
OUT = Path(sys.argv[2])
# Optional third argument: read cached pass1/pass2 arrays from an earlier run's
# intermediate directory (read-only reuse; nothing there is rewritten).
INTER = Path(sys.argv[3]) if len(sys.argv) > 3 else OUT / "intermediate"
INTER.mkdir(exist_ok=True)
ATLAS = REPO / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
ROSTER = REPO / ("Analysis/SingleCell/candidates/ambient-program-recalibration-all117-"
                 "candidate-2026-08-15-v2/results/current_source_diagnosis_donor_roster.tsv")
GTF = Path("/gpfs/commons/home/jameslee/reference_genome/cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz")
GEO = OUT / "geo_sample_characteristics"

Y_CAND = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "ZFY", "TXLNGY",
          "NLGN4Y", "TTTY14", "LINC00278"]
PROFIB = ["SPP1", "TREM2", "CD9", "GPNMB", "TGFB1", "PDGFB", "LGALS3", "FABP5", "CD63"]
MYELOID = ["Macrophages", "Mono+mono derived cells"]
LYMPHOID = ["T cells", "Circulating NK/NKT", "Resident NK", "B cells", "Plasma cells"]
IMMUNE_FAMILY = ["Macrophages", "Mono+mono derived cells", "T cells", "Circulating NK/NKT",
                 "Resident NK", "B cells", "Plasma cells", "cDC2s", "Neutrophils"]
NONIMMUNE = ["Hepatocytes", "Endothelial cells", "Fibroblasts", "Cholangiocytes"]
STAGE_CODE = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2, "Cirrhosis": 3}
N_CTRL = 50
N_NEIGH = 50
CHUNK = 100_000

log = lambda *a: print(*a, flush=True)


def decode(arr):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in arr])


# ---------------------------------------------------------------- metadata
def library_to_donor():
    m = {}
    for ds in ["GSE244832", "GSE185477", "GSE202379", "GSE136103"]:
        dp = pd.read_csv(REPO / f"data/{ds}/metadata/donor_pairing.csv", dtype=str)
        for _, r in dp.iterrows():
            for s in r["rna_srrs"].split(";"):
                assert s not in m
                m[s] = f"{ds}_{r['donor_id']}"
    runs = pd.read_csv(REPO / "Liver_Atlas/metadata/SraRunTable.csv", dtype=str)
    smp = pd.read_csv(REPO / "Liver_Atlas/metadata/GSE192740_sampleInfo_scRNAseq.tsv", sep="\t", dtype=str)
    j = runs[["Run", "shortfilename"]].merge(
        smp[["characteristics: shortFileName", "title"]].rename(
            columns={"characteristics: shortFileName": "shortfilename"}), on="shortfilename")
    j["donor_short"] = j["title"].map(lambda t: (re.search(r"H[0-9]+", t) or [None])[0])
    j = j.dropna(subset=["donor_short"])
    assert j.Run.nunique() == 48 and j.donor_short.nunique() == 19, "Liver Atlas donor map drift"
    for r, d in zip(j.Run, j.donor_short):
        m[r] = f"Liver_Atlas_{d}"
    return m


def parse_soft(path):
    recs, cur = [], None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("^SAMPLE"):
            cur = {"gsm": line.split(" = ")[1]}
            recs.append(cur)
        elif cur is not None and line.startswith("!Sample_title"):
            cur["title"] = line.split(" = ", 1)[1]
        elif cur is not None and line.startswith("!Sample_characteristics_ch1"):
            kv = line.split(" = ", 1)[1]
            if ": " in kv:
                k, v = kv.split(": ", 1)
                cur[k.strip().lower()] = v.strip()
    return pd.DataFrame(recs)


def geo_metadata():
    rows = []
    g = parse_soft(GEO / "GSE202379.soft")
    for pid, sub in g.groupby("patient id"):
        saf = sub["saf score"].iloc[0]
        if saf.lower() == "healthy control":
            F = 0
        elif saf.lower() == "end stage":
            F = 4
        else:
            mF = re.search(r"F(\d)", saf)
            F = int(mF.group(1)) if mF else np.nan
        rows.append(dict(donor=f"GSE202379_{pid}", meta_sex={"M": "Male", "F": "Female"}.get(
            sub["gender"].iloc[0], np.nan), age=pd.to_numeric(sub["age"].iloc[0], errors="coerce"),
            saf=saf, F_saf=F, geo_disease=sub["disease status"].iloc[0],
            n_geo_samples=len(sub), consistent=int(sub[["gender", "age", "saf score"]].nunique().max() == 1)))
    g = parse_soft(GEO / "GSE136103.soft")
    g = g[g.get("sex").notna()] if "sex" in g else g
    for _, r in g.iterrows():
        d = "GSE136103_" + r["title"].split("_")[0].lower()
        rows.append(dict(donor=d, meta_sex=r["sex"], age=np.nan))
    g = parse_soft(GEO / "GSE185477.soft")
    for _, r in g.iterrows():
        mm = re.search(r"Donor (\d+)", r["title"])
        if mm:
            rows.append(dict(donor=f"GSE185477_D{int(mm.group(1)):02d}",
                             meta_sex={"M": "Male", "F": "Female"}.get(r["sex"], r["sex"]), age=np.nan))
    md = pd.DataFrame(rows).drop_duplicates("donor")
    return md


def gene_chrom():
    out = {}
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split("\t", 4)
            if p[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^"]+)"', line).group(1)
            out[gid] = p[0]
    return out


# ---------------------------------------------------------------- passes
def open_atlas():
    f = h5py.File(ATLAS, "r")
    X = f["raw/X"]
    indptr = X["indptr"][:].astype(np.int64)
    genes = decode(f["raw/var/_index"][:])
    gene_ids = decode(f["var/gene_ids"][:])
    obs = f["obs"]
    lib_cats = decode(obs["sample/categories"][:]); lib_codes = obs["sample/codes"][:]
    ct_cats = decode(obs["cell_type/categories"][:]); ct_codes = obs["cell_type/codes"][:]
    ds_cats = decode(obs["dataset/categories"][:]); ds_codes = obs["dataset/codes"][:]
    return f, X, indptr, genes, gene_ids, lib_cats, lib_codes, ct_cats, ct_codes, ds_cats, ds_codes


def iter_chunks(X, indptr, n_genes):
    n = len(indptr) - 1
    for s in range(0, n, CHUNK):
        e = min(s + CHUNK, n)
        a, b = indptr[s], indptr[e]
        data = X["data"][a:b].astype(np.float32)
        idx = X["indices"][a:b]
        ptr = indptr[s:e + 1] - a
        yield s, e, sparse.csr_matrix((data, idx, ptr), shape=(e - s, n_genes))


def pass1(X, indptr, n_genes, lib_codes, ct_codes, n_lib, n_ct):
    """Per (library x cell type) gene sums and per-cell totals."""
    ng = n_lib * n_ct
    acc = np.zeros((ng, n_genes), dtype=np.float64)
    tot = np.zeros(len(indptr) - 1, dtype=np.float64)
    for s, e, M in iter_chunks(X, indptr, n_genes):
        g = lib_codes[s:e].astype(np.int64) * n_ct + ct_codes[s:e]
        G = sparse.csr_matrix((np.ones(e - s, dtype=np.float32), (g, np.arange(e - s))), shape=(ng, e - s))
        acc += (G @ M).toarray()
        tot[s:e] = np.asarray(M.sum(axis=1)).ravel()
        log("pass1", e)
    return acc, tot


def pass2(X, indptr, n_genes, cols, keep_rows_mask):
    """Per-cell counts for selected columns, plus full rows for keep_rows_mask."""
    sub_parts, full_parts = [], []
    for s, e, M in iter_chunks(X, indptr, n_genes):
        sub_parts.append(M[:, cols].tocsr())
        km = keep_rows_mask[s:e]
        full_parts.append(M[km])
        log("pass2", e)
    return sparse.vstack(sub_parts).tocsr(), sparse.vstack(full_parts).tocsr()


# ---------------------------------------------------------------- statistics
def excess_zero(y, n, g, eligible_g, donor_of_g, n_groups):
    """Return per-group observed zeros, expected zeros and cell counts (eligible groups only)."""
    sy = np.bincount(g, weights=y, minlength=n_groups)
    sn = np.bincount(g, weights=n, minlength=n_groups)
    f = np.divide(sy, sn, out=np.zeros_like(sy), where=sn > 0)
    expz = np.exp(-n * f[g])
    O = np.bincount(g, weights=(y == 0).astype(float), minlength=n_groups)
    E = np.bincount(g, weights=expz, minlength=n_groups)
    return O * eligible_g, E * eligible_g


def ols(df, formula):
    # A dataset with one donor has leverage 1 under the dataset fixed effect,
    # which makes HC3 undefined and adds no within-dataset stage information.
    if "C(dataset)" in formula:
        k = df.groupby("dataset")["dataset"].transform("size")
        df = df[k >= 2]
    m = smf.ols(formula, data=df).fit(cov_type="HC3")
    return m


def main():
    np.random.seed(SEED)
    rng = np.random.default_rng(SEED)
    summary = {}
    f, X, indptr, genes, gene_ids, lib_cats, lib_codes, ct_cats, ct_codes, ds_cats, ds_codes = open_atlas()
    n_genes = len(genes); n_lib = len(lib_cats); n_ct = len(ct_cats)
    gidx = {g: i for i, g in enumerate(genes)}
    lib2don = library_to_donor()
    lib_donor = np.array([lib2don.get(l, l) for l in lib_cats])
    lib_ds = pd.Series(ds_cats[ds_codes]).groupby(lib_codes).agg(lambda s: s.mode().iloc[0])
    lib_ds = lib_ds.reindex(range(n_lib)).values
    ct_name = {i: c for i, c in enumerate(ct_cats)}

    p1 = INTER / "pass1.npz"
    if p1.exists():
        z = np.load(p1); acc, tot = z["acc"], z["tot"]
    else:
        acc, tot = pass1(X, indptr, n_genes, lib_codes, ct_codes, n_lib, n_ct)
        np.savez(p1, acc=acc, tot=tot)
    acc3 = acc.reshape(n_lib, n_ct, n_genes)
    lib_gene = acc3.sum(axis=1)  # library x gene

    # ---- sex inference (candidate Y set; Y-set retention needs sex first)
    ycand = [g for g in Y_CAND if g in gidx]
    missing = [g for g in Y_CAND if g not in gidx]
    yc_idx = [gidx[g] for g in ycand]
    libdf = pd.DataFrame({"library": lib_cats, "donor": lib_donor, "dataset": lib_ds,
                          "Y": lib_gene[:, yc_idx].sum(1), "XIST": lib_gene[:, gidx["XIST"]],
                          "total": lib_gene.sum(1), "n_cells": np.bincount(lib_codes, minlength=n_lib)})

    def call(y, x):
        if y + x < 50:
            return "uncallable"
        r = np.log2((y + 0.5) / (x + 0.5))
        return "Male" if r > 2 else ("Female" if r < -2 else "ambiguous")
    libdf["lib_sex"] = [call(a, b) for a, b in zip(libdf.Y, libdf.XIST)]
    don = libdf.groupby("donor").agg(dataset=("dataset", "first"), Y=("Y", "sum"), XIST=("XIST", "sum"),
                                     total=("total", "sum"), n_cells=("n_cells", "sum"),
                                     n_libraries=("library", "size")).reset_index()
    don["log2_Y_over_XIST"] = np.log2((don.Y + 0.5) / (don.XIST + 0.5))
    don["sex"] = [call(a, b) for a, b in zip(don.Y, don.XIST)]
    libdf = libdf.merge(don[["donor", "sex"]].rename(columns={"sex": "donor_sex"}), on="donor")
    libdf["discordant"] = (libdf.lib_sex != libdf.donor_sex) & libdf.lib_sex.isin(["Male", "Female"])
    md = geo_metadata()
    don = don.merge(md, on="donor", how="left")
    roster = pd.read_csv(ROSTER, sep="\t")
    don = don.merge(roster[["donor", "disease_stage_coarse", "exclude"]], on="donor", how="left")
    don["stage"] = don.disease_stage_coarse.map(STAGE_CODE)
    don.to_csv(OUT / "donor_sex_calls.tsv", sep="\t", index=False)
    libdf.to_csv(OUT / "library_sex_calls.tsv", sep="\t", index=False)
    conc = don.dropna(subset=["meta_sex"])
    summary["sex_meta_concordance"] = {"n_with_meta": int(len(conc)),
                                       "n_concordant": int((conc.sex == conc.meta_sex).sum()),
                                       "discordant": conc.loc[conc.sex != conc.meta_sex, ["donor", "sex", "meta_sex"]].values.tolist()}
    summary["donor_sex_counts_all"] = don.sex.value_counts().to_dict()
    summary["discordant_libraries"] = libdf.loc[libdf.discordant, "library"].tolist()
    summary["y_candidates_missing_from_matrix"] = missing

    good_lib = ~libdf.set_index("library").loc[lib_cats, "discordant"].values
    male_lib = (libdf.set_index("library").loc[lib_cats, "donor_sex"].values == "Male") & good_lib
    fem_lib = (libdf.set_index("library").loc[lib_cats, "donor_sex"].values == "Female") & good_lib

    # ---- Y gene retention
    mrate = lib_gene[male_lib].sum(0) / lib_gene[male_lib].sum()
    frate = lib_gene[fem_lib].sum(0) / lib_gene[fem_lib].sum()
    ret = pd.DataFrame({"gene": ycand, "male_umi": [lib_gene[male_lib, gidx[g]].sum() for g in ycand],
                        "female_umi": [lib_gene[fem_lib, gidx[g]].sum() for g in ycand],
                        "male_rate": [mrate[gidx[g]] for g in ycand], "female_rate": [frate[gidx[g]] for g in ycand]})
    ret["female_over_male"] = ret.female_rate / ret.male_rate
    ret["kept"] = (ret.female_over_male < 0.05) & (ret.male_umi >= 100)
    ret.to_csv(OUT / "y_gene_retention.tsv", sep="\t", index=False)
    yset = ret.loc[ret.kept, "gene"].tolist()
    summary["y_set"] = yset
    if len(yset) < 3:
        (OUT / "STOP.txt").write_text(f"Fewer than 3 Y genes retained: {yset}\n")
        json.dump(summary, open(OUT / "summary.json", "w"), indent=2, default=str); return

    # ---- control gene matching
    chrom = gene_chrom()
    def grp_rate(names):
        cti = [i for i, c in enumerate(ct_cats) if c in names]
        s = acc3[male_lib][:, cti, :].sum(axis=(0, 1))
        return s / s.sum()
    feats = np.vstack([np.log10(grp_rate(MYELOID) + 1e-9), np.log10(grp_rate(LYMPHOID) + 1e-9),
                       np.log10(grp_rate(["Hepatocytes"]) + 1e-9)]).T
    bad = re.compile(r"^(MT-|HLA-|IG[HKL][VDJC]|TR[ABDG][VDJC])")
    auto = np.array([chrom.get(gid, "") in {f"chr{i}" for i in range(1, 23)} for gid in gene_ids])
    pool = np.where(auto & np.array([not bad.match(g) for g in genes]) &
                    ~np.isin(genes, PROFIB) & (feats > -9 + 1e-6).all(1))[0]
    neigh = {}
    for g in yset:
        d = np.sqrt(((feats[pool] - feats[gidx[g]]) ** 2).sum(1))
        neigh[g] = pool[np.argsort(d)[:N_NEIGH]]
    ctrl_sets = []
    perms = {g: rng.permutation(N_NEIGH) for g in yset}
    for k in range(N_CTRL):
        used, members = set(), []
        for g in yset:
            for j in range(N_NEIGH):
                cand = neigh[g][perms[g][(k + j) % N_NEIGH]]
                if cand not in used:
                    used.add(cand); members.append(int(cand)); break
        ctrl_sets.append(members)
    cm = pd.DataFrame([(k, genes[i], yset[j]) for k, s in enumerate(ctrl_sets) for j, i in enumerate(s)],
                      columns=["set", "gene", "matched_to"])
    cm.to_csv(OUT / "control_gene_sets.tsv", sep="\t", index=False)
    mf = pd.DataFrame(feats[[gidx[g] for g in yset] + sorted({i for s in ctrl_sets for i in s})],
                      columns=["log10_myeloid", "log10_lymphoid", "log10_hep"])
    mf.insert(0, "gene", yset + [genes[i] for i in sorted({i for s in ctrl_sets for i in s})])
    mf.insert(1, "is_Y", [True] * len(yset) + [False] * (len(mf) - len(yset)))
    mf.to_csv(OUT / "matching_features.tsv", sep="\t", index=False)

    # ---- pass 2
    ctrl_genes = sorted({i for s in ctrl_sets for i in s})
    cols = [gidx[g] for g in yset] + ctrl_genes + [gidx[g] for g in PROFIB if g in gidx]
    colpos = {c: i for i, c in enumerate(cols)}
    cell_male = male_lib[lib_codes]
    mac_cti = [i for i, c in enumerate(ct_cats) if c in MYELOID]
    keep_full = cell_male & np.isin(ct_codes, mac_cti)
    p2 = INTER / "pass2.npz"
    if p2.exists():
        z = np.load(p2, allow_pickle=True)
        sub = sparse.csr_matrix((z["sd"], z["si"], z["sp"]), shape=tuple(z["ss"]))
        full = sparse.csr_matrix((z["fd"], z["fi"], z["fp"]), shape=tuple(z["fs"]))
    else:
        sub, full = pass2(X, indptr, n_genes, cols, keep_full)
        np.savez(p2, sd=sub.data, si=sub.indices, sp=sub.indptr, ss=sub.shape,
                 fd=full.data, fi=full.indices, fp=full.indptr, fs=full.shape)
    ycol = [colpos[gidx[g]] for g in yset]
    yc = np.asarray(sub[:, ycol].sum(1)).ravel()
    S = np.zeros((len(cols), N_CTRL), dtype=np.float32)
    for k, s in enumerate(ctrl_sets):
        for i in s:
            S[colpos[i], k] = 1
    cc = np.asarray((sub @ sparse.csr_matrix(S)).todense())  # cells x 50

    n = tot
    cells = pd.DataFrame({"lib": lib_codes, "ct": ct_codes, "n": n, "Y": yc})
    cells["donor"] = lib_donor[lib_codes]
    cells["sex"] = np.where(cell_male, "Male", np.where(fem_lib[lib_codes], "Female", "other"))
    cells["cell_type"] = ct_cats[ct_codes]

    # ---- cells per type
    cpt = cells.groupby(["cell_type", "sex"]).size().unstack(fill_value=0)
    cpt.to_csv(OUT / "cells_per_type_by_sex.tsv", sep="\t")

    # ---- female false Y-positive rate
    fc = cells[cells.sex == "Female"].copy()
    fc["depth_decile"] = pd.qcut(np.log10(fc.n), 10, labels=False)
    fem_ct = fc.groupby("cell_type").agg(n_cells=("Y", "size"), frac_Ypos=("Y", lambda v: (v > 0).mean()),
                                         median_umi=("n", "median")).reset_index()
    fem_dep = fc.groupby("depth_decile").agg(n_cells=("Y", "size"), frac_Ypos=("Y", lambda v: (v > 0).mean()),
                                             median_umi=("n", "median")).reset_index()
    fem_ct.to_csv(OUT / "female_false_Ypos_by_celltype.tsv", sep="\t", index=False)
    fem_dep.to_csv(OUT / "female_false_Ypos_by_depth_decile.tsv", sep="\t", index=False)
    mc = cells[cells.sex == "Male"]
    male_raw = mc.groupby("cell_type").agg(n_cells=("Y", "size"), frac_Yneg=("Y", lambda v: (v == 0).mean()),
                                           median_umi=("n", "median")).reset_index()
    male_raw.to_csv(OUT / "male_raw_Yneg_by_celltype.tsv", sep="\t", index=False)

    # ---- per-donor excess Y-negative statistic
    male_idx = np.where(cell_male)[0]
    g = (lib_codes[male_idx].astype(np.int64) * n_ct + ct_codes[male_idx])
    ngr = n_lib * n_ct
    nm = n[male_idx]
    ncell_g = np.bincount(g, minlength=ngr)
    yg = np.bincount(g, weights=yc[male_idx], minlength=ngr)
    elig = ((ncell_g >= 20) & (yg >= 10)).astype(float)
    gr_lib = np.arange(ngr) // n_ct; gr_ct = np.arange(ngr) % n_ct

    def donor_table(y):
        O, E = excess_zero(y, nm, g, elig, None, ngr)
        t = pd.DataFrame({"donor": lib_donor[gr_lib], "cell_type": ct_cats[gr_ct], "O": O, "E": E,
                          "N": ncell_g * elig})
        t = t.groupby(["donor", "cell_type"], as_index=False)[["O", "E", "N"]].sum()
        t = t[t.N >= 50]
        t["excess"] = (t.O - t.E) / t.N
        return t.set_index(["donor", "cell_type"])
    ty = donor_table(yc[male_idx])
    tk = [donor_table(cc[male_idx, k]) for k in range(N_CTRL)]
    Ek = pd.concat([t["excess"].rename(k) for k, t in enumerate(tk)], axis=1)
    est = ty.copy()
    est["E_Y"] = est.excess
    est["E_ctrl_median"] = Ek.median(axis=1)
    est["D_pp"] = 100 * (est.E_Y - est.E_ctrl_median)
    est["obs_Yneg_frac"] = est.O / est.N
    est["exp_Yneg_frac"] = est.E / est.N
    est = est.reset_index().merge(don[["donor", "dataset", "stage", "disease_stage_coarse", "exclude",
                                       "age", "F_saf"]], on="donor", how="left")
    hep = est[est.cell_type == "Hepatocytes"].set_index("donor").D_pp
    est["D_minus_hep_pp"] = est.D_pp - est.donor.map(hep)
    est.to_csv(OUT / "donor_celltype_excess_Yneg.tsv", sep="\t", index=False)
    Dk = {}
    for k in range(N_CTRL):
        others = Ek.drop(columns=k).median(axis=1)
        Dk[k] = 100 * (Ek[k] - others)
    Dk = pd.DataFrame(Dk)

    # ---- stop rule
    stage_ok = est[(est.exclude == False) & est.stage.notna()]
    n_mac = int((stage_ok.cell_type == "Macrophages").sum())
    summary["n_male_stage_donors_with_50_macrophages"] = n_mac
    summary["n_male_stage_donors_any"] = int(((don.sex == "Male") & (don.exclude == False) & don.stage.notna()).sum())
    if n_mac < 15:
        (OUT / "STOP.txt").write_text(f"Only {n_mac} male stage donors with >=50 macrophages\n")

    # ---- calibration summary
    cal = est.groupby("cell_type").agg(n_donors=("D_pp", "size"), median_obs_Yneg=("obs_Yneg_frac", "median"),
                                       median_exp_Yneg=("exp_Yneg_frac", "median"),
                                       median_E_Y=("E_Y", "median"), median_E_ctrl=("E_ctrl_median", "median"),
                                       median_D_pp=("D_pp", "median"), mean_D_pp=("D_pp", "mean")).reset_index()
    cal.to_csv(OUT / "calibration_by_celltype.tsv", sep="\t", index=False)

    # ---- stage models
    rows, ctrl_rows = [], []
    for ct in IMMUNE_FAMILY + NONIMMUNE:
        d = stage_ok[stage_ok.cell_type == ct].copy()
        fam = ct in IMMUNE_FAMILY
        r = dict(cell_type=ct, family=fam, n_donors=len(d), n_stage_levels=d.stage.nunique(),
                 n_datasets=d.dataset.nunique())
        if len(d) >= 15 and d.stage.nunique() >= 2:
            m = ols(d, "D_pp ~ stage + C(dataset)")
            ci = m.conf_int().loc["stage"]
            r.update(tested=True, n_model=int(m.nobs), beta_pp_per_stage=m.params["stage"], ci_lo=ci[0],
                     ci_hi=ci[1], p=m.pvalues["stage"])
            dh = d.dropna(subset=["D_minus_hep_pp"])
            if len(dh) >= 15 and dh.stage.nunique() >= 2:
                mh = ols(dh, "D_minus_hep_pp ~ stage + C(dataset)")
                r.update(hepref_n=len(dh), hepref_beta=mh.params["stage"], hepref_p=mh.pvalues["stage"])
            b = d[(d.dataset == "GSE202379")].dropna(subset=["F_saf", "age"])
            if len(b) >= 10 and b.F_saf.nunique() >= 2:
                mb = ols(b, "D_pp ~ F_saf + age")
                cib = mb.conf_int().loc["F_saf"]
                r.update(B_n=len(b), B_beta_pp_per_F=mb.params["F_saf"], B_ci_lo=cib[0], B_ci_hi=cib[1],
                         B_p=mb.pvalues["F_saf"], B_age_beta=mb.params["age"], B_age_p=mb.pvalues["age"])
            # negative control slopes
            slopes = []
            for k in range(N_CTRL):
                dk = d[["donor", "dataset", "stage"]].copy()
                dk["D"] = Dk.loc[list(zip(d.donor, [ct] * len(d))), k].values
                dk = dk.dropna()
                mk = ols(dk, "D ~ stage + C(dataset)")
                slopes.append(mk.params["stage"])
                ctrl_rows.append(dict(cell_type=ct, set=k, n=len(dk), beta=mk.params["stage"], p=mk.pvalues["stage"]))
            slopes = np.array(slopes)
            r.update(ctrl_slope_median=np.median(slopes), ctrl_slope_p95=np.percentile(slopes, 95),
                     ctrl_slope_p5=np.percentile(slopes, 5),
                     Y_slope_rank_frac=(np.sum(slopes >= r["beta_pp_per_stage"]) + 1) / (N_CTRL + 1))
        else:
            r.update(tested=False)
        rows.append(r)
    res = pd.DataFrame(rows)
    tf = res.family & res.tested
    res["q_BH_immune"] = np.nan
    if tf.any():
        res.loc[tf, "q_BH_immune"] = multipletests(res.loc[tf, "p"], method="fdr_bh")[1]
    res["positive_call"] = (res.q_BH_immune < 0.05) & (res.beta_pp_per_stage > 0) & \
        (np.sign(res.get("B_beta_pp_per_F", np.nan)) == 1) & (res.beta_pp_per_stage > res.get("ctrl_slope_p95", np.inf))
    res.to_csv(OUT / "stage_association.tsv", sep="\t", index=False)
    cr = pd.DataFrame(ctrl_rows)
    if len(cr):
        cr["q_BH_within_set"] = cr.groupby("set").p.transform(
            lambda p: multipletests(p, method="fdr_bh")[1] if len(p) else p)
    cr.to_csv(OUT / "negative_control_stage_slopes.tsv", sep="\t", index=False)

    # ---- within-donor Y-neg vs Y-pos macrophages
    full_rows = np.where(keep_full)[0]
    fl = pd.DataFrame({"row": np.arange(len(full_rows)), "lib": lib_codes[full_rows],
                       "ct": ct_cats[ct_codes[full_rows]], "n": n[full_rows], "Y": yc[full_rows]})
    fl["donor"] = lib_donor[fl.lib]
    fl["bin"] = np.floor(np.log10(np.maximum(fl.n, 1)) / 0.1).astype(int)
    pf_idx = [gidx[g] for g in PROFIB if g in gidx]

    def matched_diff(neg_mask, drop_cols, types, rs):
        sel = fl[fl.ct.isin(types)].copy()
        sel["neg"] = neg_mask[sel.row.values]
        picks = []
        for (lib, b), grp in sel.groupby(["lib", "bin"]):
            ng_, pg_ = grp[grp.neg], grp[~grp.neg]
            m_ = min(len(ng_), len(pg_))
            if m_ == 0:
                continue
            picks.append(ng_.sample(m_, random_state=rs).assign(state="neg"))
            picks.append(pg_.sample(m_, random_state=rs).assign(state="pos"))
        if not picks:
            return pd.DataFrame()
        pk = pd.concat(picks)
        keep_cols = np.setdiff1d(np.arange(n_genes), drop_cols)
        out = []
        for dnr, gd in pk.groupby("donor"):
            nn = (gd.state == "neg").sum()
            if nn < 20:
                continue
            rec = dict(donor=dnr, n_matched_per_state=int(nn))
            for st in ["neg", "pos"]:
                v = np.asarray(full[gd.loc[gd.state == st, "row"].values].sum(0)).ravel()
                lib_size = v[keep_cols].sum()
                cpm = v[pf_idx] / lib_size * 1e6
                rec[f"score_{st}"] = np.mean(np.log2(cpm + 1))
                rec[f"umi_{st}"] = float(gd.loc[gd.state == st, "n"].sum())
                for gname, cv in zip([g for g in PROFIB if g in gidx], cpm):
                    rec[f"{gname}_{st}"] = np.log2(cv + 1)
            rec["diff"] = rec["score_neg"] - rec["score_pos"]
            out.append(rec)
        return pd.DataFrame(out)

    yneg_full = yc[full_rows] == 0
    wd_rows, wd_ctrl = [], []
    for label, types in [("Macrophages", ["Macrophages"]), ("Myeloid_sensitivity", MYELOID)]:
        t = matched_diff(yneg_full, [gidx[g] for g in yset], types, SEED)
        t.to_csv(OUT / f"within_donor_{label}_Y.tsv", sep="\t", index=False)
        ctrl_means = []
        for k in range(N_CTRL):
            tk_ = matched_diff(cc[full_rows, k] == 0, ctrl_sets[k], types, SEED + k + 1)
            if len(tk_):
                ctrl_means.append(tk_["diff"].mean())
                wd_ctrl.append(dict(analysis=label, set=k, n_donors=len(tk_), mean_diff=tk_["diff"].mean(),
                                    median_diff=tk_["diff"].median()))
        r = dict(analysis=label, n_donors=len(t))
        if len(t) >= 5:
            r.update(mean_diff_log2=t["diff"].mean(), median_diff_log2=t["diff"].median(),
                     n_neg_higher=int((t["diff"] > 0).sum()),
                     wilcoxon_p=stats.wilcoxon(t["diff"]).pvalue,
                     ttest_p=stats.ttest_1samp(t["diff"], 0).pvalue,
                     ci_lo=t["diff"].mean() - stats.t.ppf(0.975, len(t) - 1) * t["diff"].std(ddof=1) / np.sqrt(len(t)),
                     ci_hi=t["diff"].mean() + stats.t.ppf(0.975, len(t) - 1) * t["diff"].std(ddof=1) / np.sqrt(len(t)),
                     ctrl_mean_diff_median=np.median(ctrl_means) if ctrl_means else np.nan,
                     ctrl_mean_diff_p5=np.percentile(ctrl_means, 5) if ctrl_means else np.nan,
                     ctrl_mean_diff_p95=np.percentile(ctrl_means, 95) if ctrl_means else np.nan,
                     Y_rank_frac_ge=(np.sum(np.array(ctrl_means) >= t["diff"].mean()) + 1) / (len(ctrl_means) + 1))
        wd_rows.append(r)
    pd.DataFrame(wd_rows).to_csv(OUT / "within_donor_summary.tsv", sep="\t", index=False)
    pd.DataFrame(wd_ctrl).to_csv(OUT / "within_donor_negative_control.tsv", sep="\t", index=False)

    summary["n_donors_by_sex_in_stage_roster"] = don[(don.exclude == False) & don.stage.notna()].sex.value_counts().to_dict()
    json.dump(summary, open(OUT / "summary.json", "w"), indent=2, default=str)
    log("done")


if __name__ == "__main__":
    main()
