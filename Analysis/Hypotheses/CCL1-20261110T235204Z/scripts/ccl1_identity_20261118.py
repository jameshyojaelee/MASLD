"""CCL1 (2026-11-18 deduplicated unique_counts rerun + absolute contamination sensitivity; N03).
Replaces scalesc_human_annotated_celltypist.h5ad and Aug-13 stage metadata with NUCLEUS-FIX
unique_counts.h5ad (650,083 unique cells) and donor_metadata.tsv (strictly filtering
exclude != True & exclude_stage_analysis != True; dropping P98; D11=MASL, D14=MASH),
reading obs['donor'] directly.
Adds absolute contamination-filtered sensitivity alongside F1 (cs_all <= q95):
  - PTPRC == 0 and non-hepatocyte marker UMI sum == 0 (full panel & cholangiocyte+PTPRC).
Usage: python -I ccl1_identity_20261118.py <new_out_dir>
"""
import sys
from itertools import product
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
HYP = ROOT / "Analysis/Hypotheses"
V5 = HYP / "V5-20261021T122723Z"
V3 = HYP / "V3-20261008T015239Z"
DIG1 = HYP / "DIG1-20261024T113951Z"
NF = HYP / "NUCLEUS-FIX-20261107/unique_counts"
FLEX_H5 = V5 / "derived/gse298719_qc_cells.h5ad"
FLEX_QC = V5 / "outputs/donor_qc.tsv"
SN_DIR = V3 / "downloads_untrusted"
SN_LAB = DIG1 / "out_snrna_labels/lineage_per_nucleus.tsv.gz"
SN_META = V3 / "pseudobulk/donor_metadata.tsv"
ATLAS = NF / "unique_counts.h5ad"
DONOR_META = NF / "metadata/donor_metadata.tsv"
G289 = ROOT / "data/public_r4/GSE289173/derived_H90-20261015T120337Z/gse289173_qc_nuclei.h5ad"

assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
SEED = 20261110
rng = np.random.default_rng(SEED)
K_MATCH, MIN_POS, MIN_NEG, MIN_CHOL, MIN_DONORS = 10, 10, 50, 10, 5
MARGIN_H, MARGIN_D = 1.0, 2.0
HEP = ["ALB", "APOA1", "SERPINA1", "TF", "CYP3A4", "CYP2E1"]
DUC = ["KRT19", "KRT7", "EPCAM", "SOX9", "CFTR", "SPP1", "TACSTD2"]
NONHEP_EXTRA = ["COL1A1", "COL3A1", "DCN", "ACTA2", "PECAM1", "STAB2", "VWF", "CDH5", "CD163", "CD68", "LYZ"]
NONHEP_ALL = sorted(set(DUC) | set(NONHEP_EXTRA))
SCORE = ["KRT8", "SQSTM1"]
CORE = ["CXCL1", "CXCL8", "BIRC3", "TNFAIP3", "ICAM1", "ZC3H12A", "RELB", "NFKBIA", "NFKB2", "MAP3K14", "UBD", "FCAMR"]
STRICT_GENES = ["CCL20"] + HEP + DUC + SCORE + CORE + ["PTPRC"]
GENES = STRICT_GENES + [g for g in NONHEP_EXTRA if g not in STRICT_GENES]
PRIMARY = ["GSE298719_flex", "GSE256398_snrna"]
SECONDARY = ["GSE244832", "GSE202379", "GSE289173"]


# ---------------------------------------------------------------- statistics
def flips(n):
    if n <= 20:
        return np.array(list(product([1, -1], repeat=n)), dtype=np.int8)
    return np.where(np.random.default_rng(SEED).random((200_000, n)) < 0.5, 1, -1).astype(np.int8)


def signflip_greater(v):
    """One-sided exact sign-flip: H1 mean(v) > 0."""
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 3:
        return np.nan
    return float(np.mean((flips(len(v)) * v).mean(1) >= v.mean() - 1e-12))


def signflip_two(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 3:
        return np.nan
    return float(np.mean(np.abs((flips(len(v)) * v).mean(1)) >= abs(v.mean()) - 1e-12))


def holm(p):
    p = np.asarray(p, float); out = np.full(len(p), np.nan); ok = np.where(np.isfinite(p))[0]
    if len(ok):
        q = p[ok]; o = np.argsort(q); adj = np.maximum.accumulate(q[o] * (len(q) - np.arange(len(q))))
        out[ok[o]] = np.minimum(1, adj)
    return out


def pct(x):
    return (stats.rankdata(x) - 0.5) / len(x)


def log2cpm(counts, lib):
    return np.log2(counts / lib * 1e6 + 1)


def flag_true(s):
    return s.astype(str).str.strip().str.lower().eq("true")


# ---------------------------------------------------------------- loading
def obs_str(o, k):
    g = o[k]
    if isinstance(g, h5py.Group) and "categories" in g:
        c = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]])
        codes = g["codes"][:]
        return np.where(codes >= 0, c[np.maximum(codes, 0)], "")
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[:]])


def obs_float(o, k):
    if k not in o:
        return None
    g = o[k]
    if isinstance(g, h5py.Dataset):
        return np.asarray(g[:], float)
    if isinstance(g, h5py.Group) and "categories" in g:
        c = np.array([float(x) for x in g["categories"][:]])
        codes = g["codes"][:]
        return np.where(codes >= 0, c[np.maximum(codes, 0)], np.nan)
    return None


def stream_csr(grp, rows, ncol, gcols, block):
    indptr = grp["indptr"][:]
    umi = np.zeros(len(rows)); ngene = np.zeros(len(rows)); G = np.zeros((len(rows), len(gcols)), np.float32)
    nrow = len(indptr) - 1; pos = 0
    for b0 in range(0, nrow, block):
        b1 = min(b0 + block, nrow)
        sel = rows[(rows >= b0) & (rows < b1)]
        if len(sel) == 0:
            continue
        lo, hi = indptr[b0], indptr[b1]
        B = sparse.csr_matrix((grp["data"][lo:hi].astype(np.float32), grp["indices"][lo:hi], indptr[b0:b1 + 1] - lo), shape=(b1 - b0, ncol))
        B = B[sel - b0]
        k = len(sel)
        umi[pos:pos + k] = np.asarray(B.sum(axis=1)).ravel(); ngene[pos:pos + k] = np.diff(B.indptr)
        G[pos:pos + k] = B[:, gcols].toarray()
        pos += k
        print("  block", b0, k, flush=True)
    assert pos == len(rows)
    assert np.allclose(G, np.round(G)), "not integer counts"
    return umi, ngene, G


def gene_cols(gnames):
    """Listed genes by symbol; strict symbols cannot be duplicated."""
    first, cnt = {}, {}
    for i, g in enumerate(gnames):
        first.setdefault(g, i); cnt[g] = cnt.get(g, 0) + 1
    dl = [g for g in STRICT_GENES if cnt.get(g, 0) > 1]
    assert not dl, f"listed genes with duplicate rows: {dl}"
    present = [g for g in GENES if g in first and cnt.get(g, 0) == 1]
    return present, [first[g] for g in present]


def load_flex():
    dq = pd.read_csv(FLEX_QC, sep="\t", dtype={"id": str}); dq = dq[dq.included]
    with h5py.File(FLEX_H5, "r") as f:
        o = f["obs"]; gsm = obs_str(o, "gsm"); lab = obs_str(o, "label")
        conf = obs_float(o, "conf")
        if conf is None:
            conf = np.full(len(gsm), np.nan)
        gn = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var"]["_index"][:]]
        rows = np.where(np.isin(gsm, dq.gsm.to_numpy()))[0]
        present, gc = gene_cols(gn)
        umi, ng, G = stream_csr(f["X"], rows, len(gn), gc, 20_000)
    lab = np.where(lab[rows] == "Hepatocytes", "hep", np.where(lab[rows] == "Cholangiocytes", "chol", "other"))
    dis = dict(zip(dq.gsm, dq.disease.astype(int)))
    return dict(donor=gsm[rows], lineage=lab, umi=umi, ngene=ng, G=G, genes=present, conf=np.asarray(conf, float)[rows], disease=dis)


def load_snrna():
    import scanpy as sc
    labs = pd.read_csv(SN_LAB, sep="\t", index_col=0)
    md = pd.read_csv(SN_META, sep="\t")
    files = sorted(SN_DIR.glob("GSM*_CB_raw_feature_bc_matrix_filtered.h5"))
    assert len(files) == 17
    meta_parts, G_parts, present0 = [], [], None
    for fpath in files:
        sample = fpath.name.split("_")[1]
        a = sc.read_10x_h5(fpath, gex_only=True)
        dup = pd.Index(a.var_names)[pd.Index(a.var_names).duplicated(keep=False)]
        dl = sorted(set(dup) & set(STRICT_GENES))
        assert not dl, f"listed gene has duplicate symbol rows: {dl}"
        present, gc = gene_cols(list(a.var_names))
        a.var_names_make_unique()
        a.obs_names = [f"{sample}_{b}" for b in a.obs_names]
        X = sparse.csr_matrix(a.X)
        mt = np.char.startswith(np.asarray(a.var_names, dtype=str), "MT-")
        umi = np.asarray(X.sum(axis=1)).ravel(); ng = np.diff(X.indptr)
        pmt = np.asarray(X[:, mt].sum(axis=1)).ravel() / np.maximum(umi, 1) * 100
        keep = (ng >= 200) & (pmt < 5) & np.isin(a.obs_names, labs.index)
        assert present0 is None or present == present0
        present0 = present
        G_parts.append(X[keep][:, gc].toarray().astype(np.float32))
        ob = np.asarray(a.obs_names)[keep]
        assert len(ob) == int((labs.donor == sample).sum()), f"{sample}: {len(ob)} vs {(labs.donor == sample).sum()}"
        meta_parts.append(pd.DataFrame(dict(cell=ob, donor=sample, umi=umi[keep], ngene=ng[keep])))
        print("  snrna", sample, len(ob), flush=True)
        del a, X
    P = pd.concat(meta_parts, ignore_index=True); G = np.vstack(G_parts); del G_parts
    assert np.allclose(G, np.round(G))
    lin = labs.loc[P.cell, "lineage"].to_numpy()
    lab = np.where(lin == "hepatocyte", "hep", np.where(lin == "cholangiocyte", "chol", "other"))
    return dict(donor=P.donor.to_numpy(), lineage=lab, umi=P.umi.to_numpy(float), ngene=P.ngene.to_numpy(float), G=G,
                genes=present0, conf=np.full(len(P), np.nan), disease=dict(zip(md.donor.astype(str), md.disease.astype(int))))


def load_atlas():
    meta = pd.read_csv(DONOR_META, sep="\t", dtype=str)
    eligible = meta[
        meta.dataset.isin(["GSE202379", "GSE244832"])
        & (~flag_true(meta.exclude))
        & (~flag_true(meta.exclude_stage_analysis))
        & meta.disease_stage_coarse.isin(["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"])
    ].copy()
    dis = {r.donor: int(r.disease_stage_coarse != "Healthy") for r in eligible.itertuples(index=False)}
    assert "GSE202379_P98" not in dis
    assert dis["GSE244832_D11"] == 1 and dis["GSE244832_D14"] == 1
    with h5py.File(ATLAS, "r") as f:
        o = f["obs"]
        ds_all = obs_str(o, "dataset"); donor_all = obs_str(o, "donor"); ct_all = obs_str(o, "cell_type"); prep = obs_str(o, "preparation_method")
        conf_all = obs_float(o, "cell_type_conf")
        if conf_all is None:
            conf_all = np.full(len(ds_all), np.nan)
        gn = [x.decode() if isinstance(x, bytes) else str(x) for x in f["raw"]["var"]["_index"][:]]
        rows = np.where(np.isin(ds_all, ["GSE202379", "GSE244832"]) & np.isin(donor_all, list(dis)))[0]
        present, gc = gene_cols(gn)
        umi, ng, G = stream_csr(f["raw"]["X"], rows, len(gn), gc, 20_000)
    print("atlas preparation:", pd.crosstab(ds_all[rows], prep[rows]).to_string(), flush=True)
    assert not pd.Series(prep[rows]).isin(["cd45_positive", "cd45_negative", "sorted"]).any()
    ct = ct_all[rows]; don_r = donor_all[rows]
    out = {}
    for ds in ["GSE244832", "GSE202379"]:
        m = ds_all[rows] == ds
        lab = np.where(ct[m] == "Hepatocytes", "hep", np.where(ct[m] == "Cholangiocytes", "chol", "other"))
        out[ds] = dict(donor=don_r[m], lineage=lab, umi=umi[m], ngene=ng[m], G=G[m], genes=present,
                       conf=np.asarray(conf_all, float)[rows][m], disease=dis)
    return out


def load_289173():
    with h5py.File(G289, "r") as f:
        o = f["obs"]; gsm = obs_str(o, "gsm"); grp = obs_str(o, "group"); lab = obs_str(o, "label")
        conf = obs_float(o, "conf")
        if conf is None:
            conf = np.full(len(gsm), np.nan)
        gn = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var"]["_index"][:]]
        rows = np.where(gsm != "GSM8785704")[0]  # exclude sex-discordant FL723
        present, gc = gene_cols(gn)
        umi, ng, G = stream_csr(f["X"], rows, len(gn), gc, 20_000)
    lab = np.where(lab[rows] == "Hepatocytes", "hep", np.where(lab[rows] == "Cholangiocytes", "chol", "other"))
    dis = {d: int(g != "ctrl") for d, g in zip(gsm[rows], grp[rows])}
    return dict(donor=gsm[rows], lineage=lab, umi=umi, ngene=ng, G=G, genes=present, conf=np.asarray(conf, float)[rows], disease=dis)


# ---------------------------------------------------------------- per cohort
def match(Fp, Fn):
    _, nn = cKDTree(Fn).query(Fp, k=K_MATCH)
    return np.bincount(nn.ravel(), minlength=len(Fn)).astype(float)


def analyse(tag, D):
    gi = {g: i for i, g in enumerate(D["genes"])}
    hepg = [g for g in HEP if g in gi]; ducg = [g for g in DUC if g in gi]
    nonhepg = [g for g in NONHEP_ALL if g in gi]
    G = D["G"]; umi = D["umi"]; ng = D["ngene"]
    cp = lambda gs, ix: np.log2(1 + 1e4 * G[np.ix_(ix, [gi[g] for g in gs])].sum(1) / umi[ix])
    det_rows, pair_rows, abs_pair_rows, gene_rows, ctrl_rows, geo_rows, share_rows, conc_rows, plot_rows = [], [], [], [], [], [], [], [], []
    for d in sorted(np.unique(D["donor"])):
        dis = D["disease"].get(d, np.nan)
        m_d = D["donor"] == d
        H = np.where(m_d & (D["lineage"] == "hep"))[0]; C = np.where(m_d & (D["lineage"] == "chol"))[0]
        ccl = G[:, gi["CCL20"]]
        pos = ccl[H] > 0
        det_rows.append(dict(cohort=tag, donor=d, disease=dis, n_hep=len(H), n_chol=len(C), n_ccl20_pos_hep=int(pos.sum()),
                             frac_ccl20_pos_hep=pos.mean() if len(H) else np.nan, ccl20_umi_hep=float(ccl[H].sum()),
                             ccl20_umi_chol=float(ccl[C].sum()), ccl20_umi_all=float(ccl[m_d].sum())))
        if len(H) < 20:
            continue
        hs_all, cs_all = cp(hepg, H), cp(ducg, H)
        thr95 = np.quantile(cs_all, 0.95)
        krt = (G[H, gi["KRT19"]] > 0) | (G[H, gi["KRT7"]] > 0)
        ptprc_sum = G[H, gi["PTPRC"]] if "PTPRC" in gi else np.zeros(len(H))
        duc_sum = G[np.ix_(H, [gi[g] for g in ducg])].sum(1) if ducg else np.zeros(len(H))
        nonhep_sum = G[np.ix_(H, [gi[g] for g in nonhepg])].sum(1) if nonhepg else np.zeros(len(H))
        clean_duc_ptprc = (ptprc_sum == 0) & (duc_sum == 0)
        clean_abs = (ptprc_sum == 0) & (nonhep_sum == 0)
        # hepatocyte share of CCL20 UMIs and filters (MASLD donors with any CCL20)
        tot = ccl[m_d].sum()
        if dis == 1 and tot > 0:
            hep_c = ccl[H]
            pos_hep_c = hep_c > 0
            share_rows.append(dict(
                cohort=tag, donor=d, ccl20_umi_all=float(tot), share_hep=hep_c.sum() / tot,
                share_hep_F1=hep_c[cs_all <= thr95].sum() / tot, share_hep_F2=hep_c[~krt].sum() / tot,
                share_hep_F_abs_duc_ptprc=hep_c[clean_duc_ptprc].sum() / tot,
                share_hep_F_abs=hep_c[clean_abs].sum() / tot,
                share_chol=ccl[C].sum() / tot,
                frac_hep_removed_F1=(cs_all > thr95).mean(), frac_hep_removed_F2=krt.mean(),
                frac_hep_removed_F_abs_duc_ptprc=(~clean_duc_ptprc).mean(),
                frac_hep_removed_F_abs=(~clean_abs).mean(),
                n_ccl20_pos_hep=int(pos_hep_c.sum()),
                n_ccl20_pos_F_abs=int((pos_hep_c & clean_abs).sum()),
                frac_ccl20_pos_kept_F_abs=(pos_hep_c & clean_abs).sum() / pos_hep_c.sum() if pos_hep_c.any() else np.nan,
                frac_ccl20_pos_kept_F1=(pos_hep_c & (cs_all <= thr95)).sum() / pos_hep_c.sum() if pos_hep_c.any() else np.nan))
        # concentration inside hepatocytes
        if dis == 1 and ccl[H].sum() > 0:
            c = np.sort(ccl[H])[::-1]; cum = np.cumsum(c) / c.sum()
            k60 = int(np.searchsorted(cum, 0.6 - 1e-12) + 1)
            core = [g for g in CORE if g in gi]
            L = np.log1p(G[np.ix_(H, [gi[g] for g in core])] * (1e4 / umi[H][:, None]))
            sd = L.std(0); sd[sd == 0] = np.nan; Z = (L - L.mean(0)) / sd
            sc_ = np.nan_to_num(np.nanmean(Z, 1))
            B = np.column_stack([np.ones(len(H)), stats.rankdata(np.log(umi[H])), stats.rankdata(np.log(ng[H]))])
            r = sc_ - B @ np.linalg.lstsq(B, sc_, rcond=None)[0]
            top = pct(r) >= 0.95
            top_raw = pct(sc_ + 1e-9 * stats.rankdata(np.log(umi[H]))) >= 0.95
            conc_rows.append(dict(cohort=tag, donor=d, n_hep=len(H), ccl20_umi_hep=float(ccl[H].sum()), frac_hep_ccl20_pos=pos.mean(),
                                  frac_hep_holding_60pct=k60 / len(H), top5_core_share_of_hep_ccl20=ccl[H][top].sum() / ccl[H].sum(),
                                  top5_core_share_of_hep_umi=umi[H][top].sum() / umi[H].sum(),
                                  top5_rawcore_share_of_hep_ccl20=ccl[H][top_raw].sum() / ccl[H].sum(),
                                  top5_rawcore_share_of_hep_umi=umi[H][top_raw].sum() / umi[H].sum(),
                                  hep_share=ccl[H].sum() / tot if tot > 0 else np.nan))
        # positive control: cholangiocytes vs depth-matched hepatocytes
        if dis == 1 and len(C) >= MIN_CHOL:
            HC = np.concatenate([H, C])
            F = np.column_stack([pct(np.log(umi[HC])), pct(np.log(ng[HC]))])
            w = match(F[len(H):], F[:len(H)])
            a = log2cpm(G[C].sum(0), umi[C].sum()); b = log2cpm(w @ G[H], w @ umi[H])
            diff = a - b
            ctrl_rows.append(dict(cohort=tag, donor=d, n_chol=len(C), dH=np.mean([diff[gi[g]] for g in hepg]),
                                  dD=np.mean([diff[gi[g]] for g in ducg]), **{f"d_{g}": diff[gi[g]] for g in hepg + ducg}))
        # CCL20+ vs matched CCL20- hepatocytes
        if dis != 1 or pos.sum() < MIN_POS or (~pos).sum() < MIN_NEG:
            continue
        P = H[pos]; N = H[~pos]
        inj_L = np.log1p(G[np.ix_(H, [gi[g] for g in SCORE])] * (1e4 / umi[H][:, None]))
        sd = inj_L.std(0); sd[sd == 0] = np.nan; inj = np.nan_to_num(np.nanmean((inj_L - inj_L.mean(0)) / sd, 1))
        B = np.column_stack([np.ones(len(H)), stats.rankdata(np.log(umi[H])), stats.rankdata(np.log(ng[H]))])
        inj_r = inj - B @ np.linalg.lstsq(B, inj, rcond=None)[0]
        feats = {"depth": np.column_stack([pct(np.log(umi[H])), pct(np.log(ng[H]))])}
        feats["depth+injury"] = np.column_stack([feats["depth"], pct(inj_r)])
        chol_pb = log2cpm(G[C].sum(0), umi[C].sum()) if len(C) >= MIN_CHOL else np.full(len(gi), np.nan)
        hep_pb = log2cpm(G[H].sum(0), umi[H].sum())
        for mode, F in feats.items():
            w = match(F[pos], F[~pos])
            pp = log2cpm(G[P].sum(0), umi[P].sum()); mm = log2cpm(w @ G[N], w @ umi[N])
            vs_m = pp - mm; vs_c = pp - chol_pb
            hs, cs = cp(hepg, H), cp(ducg, H)
            pair_rows.append(dict(cohort=tag, donor=d, mode=mode, n_pos=len(P), n_chol=len(C),
                                  logumi_ratio=np.log2(umi[P].mean() / ((w @ umi[N]) / w.sum())),
                                  dH=np.mean([vs_m[gi[g]] for g in hepg]), dD_vs_matched=np.mean([vs_m[gi[g]] for g in ducg]),
                                  gD=np.mean([vs_c[gi[g]] for g in ducg]) if len(C) >= MIN_CHOL else np.nan,
                                  gH_vs_chol=np.mean([vs_c[gi[g]] for g in hepg]) if len(C) >= MIN_CHOL else np.nan,
                                  chol_minus_hep_D=np.mean([(chol_pb - hep_pb)[gi[g]] for g in ducg]) if len(C) >= MIN_CHOL else np.nan,
                                  PTPRC_vs_matched=vs_m[gi["PTPRC"]] if "PTPRC" in gi else np.nan))
            for g in hepg + ducg + ["CCL20", "PTPRC"] + SCORE:
                if g in gi:
                    gene_rows.append(dict(cohort=tag, donor=d, mode=mode, gene=g, set="H" if g in HEP else "D" if g in DUC else "other",
                                          pos=pp[gi[g]], matched=mm[gi[g]], chol=chol_pb[gi[g]], hep_all=hep_pb[gi[g]],
                                          vs_matched=vs_m[gi[g]], vs_chol=vs_c[gi[g]]))
            if mode == "depth":
                Nw = w / w.sum()
                cf = D["conf"]
                geo_rows.append(dict(cohort=tag, donor=d, n_pos=len(P), n_chol=len(C),
                                     hep_score_pos=np.median(hs[pos]), hep_score_matched=np.average(hs[~pos], weights=w) if w.sum() else np.nan,
                                     hep_score_hep_median=np.median(hs), hep_score_chol=np.median(cp(hepg, C)) if len(C) else np.nan,
                                     chol_score_pos=np.median(cs[pos]), chol_score_matched=Nw @ cs[~pos],
                                     chol_score_hep_median=np.median(cs), chol_score_chol=np.median(cp(ducg, C)) if len(C) else np.nan,
                                     frac_above_hep95_pos=(cs[pos] > thr95).mean(), frac_above_hep95_matched=Nw @ (cs[~pos] > thr95),
                                     frac_krt_pos=krt[pos].mean(), frac_krt_matched=Nw @ krt[~pos],
                                     frac_abs_clean_pos=clean_abs[pos].mean(), frac_abs_clean_matched=Nw @ clean_abs[~pos],
                                     conf_pos=np.nanmedian(cf[P]) if np.isfinite(cf[P]).any() else np.nan,
                                     conf_matched=(Nw @ np.nan_to_num(cf[N])) if np.isfinite(cf[N]).any() else np.nan))
                geo_rows[-1]["t_position"] = ((geo_rows[-1]["chol_score_pos"] - geo_rows[-1]["chol_score_hep_median"]) /
                                              (geo_rows[-1]["chol_score_chol"] - geo_rows[-1]["chol_score_hep_median"])) if len(C) else np.nan
                plot_rows.append(pd.DataFrame(dict(donor=d, kind=np.where(pos, "ccl20_pos_hep", "hep"), hep_score=hs, chol_score=cs)))
                if len(C):
                    plot_rows.append(pd.DataFrame(dict(donor=d, kind="chol", hep_score=cp(hepg, C), chol_score=cp(ducg, C))))
        # absolute contamination-filtered identity test (PTPRC==0 & all non-hepatocyte markers==0)
        H_abs = H[clean_abs]; pos_abs = ccl[H_abs] > 0
        if pos_abs.sum() >= MIN_POS and (~pos_abs).sum() >= MIN_NEG:
            Pa = H_abs[pos_abs]; Na = H_abs[~pos_abs]
            Fa = np.column_stack([pct(np.log(umi[H_abs])), pct(np.log(ng[H_abs]))])
            wa = match(Fa[pos_abs], Fa[~pos_abs])
            ppa = log2cpm(G[Pa].sum(0), umi[Pa].sum()); mma = log2cpm(wa @ G[Na], wa @ umi[Na])
            vs_ma = ppa - mma; vs_ca = ppa - chol_pb
            abs_pair_rows.append(dict(
                cohort=tag, donor=d, n_hep_abs=len(H_abs), n_pos_abs=len(Pa), n_chol=len(C),
                dH=np.mean([vs_ma[gi[g]] for g in hepg]),
                dD_vs_matched=np.mean([vs_ma[gi[g]] for g in ducg]),
                gD=np.mean([vs_ca[gi[g]] for g in ducg]) if len(C) >= MIN_CHOL else np.nan))
        print(tag, d, len(H), len(C), int(pos.sum()), flush=True)
    f = lambda r: pd.DataFrame(r)
    print(tag, "genes present", D["genes"], "labels", pd.Series(D["lineage"]).value_counts().to_dict(), flush=True)
    return dict(det=f(det_rows), pair=f(pair_rows), abs_pair=f(abs_pair_rows), gene=f(gene_rows), ctrl=f(ctrl_rows),
                geo=f(geo_rows), share=f(share_rows), conc=f(conc_rows),
                plot=pd.concat(plot_rows, ignore_index=True) if plot_rows else pd.DataFrame())


# ---------------------------------------------------------------- run
R = {}
print("loading Flex", flush=True); R["GSE298719_flex"] = analyse("GSE298719_flex", load_flex())
print("loading snRNA", flush=True); R["GSE256398_snrna"] = analyse("GSE256398_snrna", load_snrna())
print("loading atlas (deduplicated unique_counts)", flush=True)
AT = load_atlas()
for ds in ["GSE244832", "GSE202379"]:
    R[ds] = analyse(ds, AT[ds])
del AT
print("loading GSE289173", flush=True); R["GSE289173"] = analyse("GSE289173", load_289173())

cat_ = lambda k: pd.concat([R[c][k] for c in R if len(R[c][k])], ignore_index=True)
DET, PAIR, ABS_PAIR, GENE, CTRL, GEO, SHARE, CONC = (
    cat_(k) for k in ["det", "pair", "abs_pair", "gene", "ctrl", "geo", "share", "conc"])
for name, T in [("detectability_by_donor", DET), ("pair_by_donor", PAIR), ("abs_clean_pair_by_donor", ABS_PAIR),
                ("gene_by_donor", GENE), ("control_by_donor", CTRL), ("geometry_by_donor", GEO),
                ("hep_share_filters_by_donor", SHARE), ("concentration_by_donor", CONC)]:
    T.to_csv(OUT / f"{name}.tsv", sep="\t", index=False)

# detectability per cohort
elig = (PAIR[PAIR["mode"] == "depth"].groupby("cohort").donor.nunique()).reindex(list(R)).fillna(0).astype(int)
abs_elig = (ABS_PAIR.groupby("cohort").donor.nunique()).reindex(list(R)).fillna(0).astype(int) if len(ABS_PAIR) else pd.Series(0, index=list(R))
DS = DET[DET.disease == 1].groupby("cohort").agg(masld_donors=("donor", "nunique"), median_frac_ccl20_pos_hep=("frac_ccl20_pos_hep", "median"),
                                                 max_frac_ccl20_pos_hep=("frac_ccl20_pos_hep", "max"),
                                                 donors_ge10_pos=("n_ccl20_pos_hep", lambda v: int((v >= MIN_POS).sum()))).reindex(list(R))
DS["eligible_donors"] = elig; DS["enters"] = DS.eligible_donors >= MIN_DONORS
DS["eligible_donors_abs_clean"] = abs_elig; DS["enters_abs_clean"] = DS.eligible_donors_abs_clean >= MIN_DONORS
DS.to_csv(OUT / "cohort_detectability.tsv", sep="\t")

# primary and secondary tests
rows = []
for c in R:
    for mode in ["depth", "depth+injury"]:
        o = PAIR[(PAIR.cohort == c) & (PAIR["mode"] == mode)]
        if not DS.loc[c, "enters"]:
            rows.append(dict(cohort=c, mode=mode, n_donors_H=len(o), status="CCL20 not detectable enough")); continue
        gd = o.gD.dropna()
        pH = signflip_greater(o.dH + MARGIN_H); pD = signflip_greater(-(gd + MARGIN_D)) if len(gd) >= MIN_DONORS else np.nan
        rows.append(dict(cohort=c, mode=mode, status="analysed", n_donors_H=len(o), mean_dH=o.dH.mean(), min_dH=o.dH.min(), max_dH=o.dH.max(),
                         donors_dH_gt_minus1=int((o.dH > -MARGIN_H).sum()), p_H=pH, n_donors_D=len(gd), mean_gD=gd.mean(), min_gD=gd.min(),
                         max_gD=gd.max(), donors_gD_lt_minus2=int((gd < -MARGIN_D).sum()), p_D=pD, p_IUT=max(pH, pD) if np.isfinite(pD) else np.nan,
                         min_attainable_p=2.0 ** -min(len(o), len(gd)) if len(gd) else np.nan,
                         mean_dD_vs_matched=o.dD_vs_matched.mean(), donors_dD_vs_matched_pos=int((o.dD_vs_matched > 0).sum()),
                         p_dD_vs_matched_two_sided=signflip_two(o.dD_vs_matched)))
T = pd.DataFrame(rows)
T["family"] = np.where(T.cohort.isin(PRIMARY) & (T["mode"] == "depth"), "primary", "secondary")
T["p_holm"] = np.nan
for fam_mask in [(T.family == "primary"), (T.family == "secondary") & (T["mode"] == "depth"), (T["mode"] == "depth+injury")]:
    idx = T.index[fam_mask & (T.status == "analysed")]
    T.loc[idx, "p_holm"] = holm(T.loc[idx, "p_IUT"].to_numpy())
T.to_csv(OUT / "tests.tsv", sep="\t", index=False)

# absolute clean sensitivity tests
abs_rows = []
for c in R:
    o = ABS_PAIR[ABS_PAIR.cohort == c] if len(ABS_PAIR) else pd.DataFrame()
    if len(o) < MIN_DONORS:
        abs_rows.append(dict(cohort=c, n_donors_H=len(o), status="fewer than 5 eligible donors after absolute contamination filter")); continue
    gd = o.gD.dropna()
    pH = signflip_greater(o.dH + MARGIN_H); pD = signflip_greater(-(gd + MARGIN_D)) if len(gd) >= MIN_DONORS else np.nan
    abs_rows.append(dict(cohort=c, status="analysed", n_donors_H=len(o), mean_dH=o.dH.mean(), min_dH=o.dH.min(), max_dH=o.dH.max(),
                         p_H=pH, n_donors_D=len(gd), mean_gD=gd.mean(), p_D=pD, p_IUT=max(pH, pD) if np.isfinite(pD) else np.nan,
                         min_attainable_p=2.0 ** -min(len(o), len(gd)) if len(gd) else np.nan))
Tabs = pd.DataFrame(abs_rows)
Tabs["family"] = np.where(Tabs.cohort.isin(PRIMARY), "primary", "secondary")
Tabs["p_holm"] = np.nan
for fam_mask in [(Tabs.family == "primary"), (Tabs.family == "secondary")]:
    idx = Tabs.index[fam_mask & (Tabs.status == "analysed")]
    Tabs.loc[idx, "p_holm"] = holm(Tabs.loc[idx, "p_IUT"].to_numpy())
Tabs.to_csv(OUT / "abs_clean_tests.tsv", sep="\t", index=False)

# per gene summaries (descriptive)
gs = []
for (c, mode, g), o in GENE.groupby(["cohort", "mode", "gene"]):
    for col in ["vs_matched", "vs_chol"]:
        v = o[col].dropna()
        if len(v):
            gs.append(dict(cohort=c, mode=mode, gene=g, set=o.set.iloc[0], contrast=col, n=len(v), mean=v.mean(), min=v.min(), max=v.max(),
                           n_pos=int((v > 0).sum()), p_two_sided=signflip_two(v)))
pd.DataFrame(gs).to_csv(OUT / "gene_summary.tsv", sep="\t", index=False)

# control summary
cs_ = []
for c, o in CTRL.groupby("cohort"):
    cs_.append(dict(cohort=c, n=len(o), mean_dH=o.dH.mean(), donors_dH_le_minus1=int((o.dH <= -1).sum()), mean_dD=o.dD.mean(),
                    donors_dD_ge_2=int((o.dD >= 2).sum()), KRT19_pos=int((o.d_KRT19 > 0).sum()) if "d_KRT19" in o else np.nan,
                    KRT7_pos=int((o.d_KRT7 > 0).sum()) if "d_KRT7" in o else np.nan, ALB_neg=int((o.d_ALB < 0).sum()),
                    mean_KRT19=o.d_KRT19.mean() if "d_KRT19" in o else np.nan, mean_KRT7=o.d_KRT7.mean() if "d_KRT7" in o else np.nan,
                    mean_ALB=o.d_ALB.mean(),
                    passes=bool(o.dH.mean() <= -1 and o.dD.mean() >= 2 and (o.dH <= -1).mean() >= 0.8 and (o.dD >= 2).mean() >= 0.8)))
pd.DataFrame(cs_).to_csv(OUT / "control_summary.tsv", sep="\t", index=False)


# geometry, share and concentration summaries
def rng_summary(T_, cols):
    out = []
    for c, o in T_.groupby("cohort"):
        r = dict(cohort=c, n=len(o))
        for k in cols:
            v = o[k].dropna()
            r[f"{k}_median"] = v.median() if len(v) else np.nan
            r[f"{k}_min"] = v.min() if len(v) else np.nan
            r[f"{k}_max"] = v.max() if len(v) else np.nan
        out.append(r)
    return pd.DataFrame(out)


rng_summary(GEO, ["chol_score_pos", "chol_score_matched", "chol_score_chol", "t_position", "frac_above_hep95_pos",
                  "frac_above_hep95_matched", "frac_krt_pos", "frac_krt_matched", "frac_abs_clean_pos",
                  "frac_abs_clean_matched", "conf_pos", "conf_matched"]).to_csv(OUT / "geometry_summary.tsv", sep="\t", index=False)
rng_summary(SHARE, ["share_hep", "share_hep_F1", "share_hep_F2", "share_hep_F_abs_duc_ptprc", "share_hep_F_abs",
                    "share_chol", "frac_hep_removed_F1", "frac_hep_removed_F2", "frac_hep_removed_F_abs_duc_ptprc",
                    "frac_hep_removed_F_abs", "frac_ccl20_pos_kept_F1", "frac_ccl20_pos_kept_F_abs"]).to_csv(
    OUT / "hep_share_filters_summary.tsv", sep="\t", index=False)
CONC["top5_share_of_all_ccl20"] = CONC.top5_core_share_of_hep_ccl20 * CONC.hep_share
CONC["top5_raw_share_of_all_ccl20"] = CONC.top5_rawcore_share_of_hep_ccl20 * CONC.hep_share
rng_summary(CONC, ["frac_hep_ccl20_pos", "frac_hep_holding_60pct", "top5_core_share_of_hep_ccl20", "top5_core_share_of_hep_umi",
                   "top5_rawcore_share_of_hep_ccl20", "top5_rawcore_share_of_hep_umi", "top5_raw_share_of_all_ccl20",
                   "hep_share", "top5_share_of_all_ccl20"]).to_csv(OUT / "concentration_summary.tsv", sep="\t", index=False)

# score plot per cohort
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
for c in R:
    Pl = R[c]["plot"]
    if Pl.empty:
        continue
    fig, ax = plt.subplots(figsize=(3.2, 3.0))
    hep = Pl[Pl.kind == "hep"]
    if len(hep) > 30000:
        hep = hep.iloc[rng.choice(len(hep), 30000, replace=False)]
    for sub, col, s, lab in [(hep, "#9E9E9E", 1, "Hepatocytes"), (Pl[Pl.kind == "chol"], "#2166AC", 3, "Cholangiocytes"),
                             (Pl[Pl.kind == "ccl20_pos_hep"], "#C51B7D", 3, "CCL20+ hepatocytes")]:
        ax.scatter(sub.hep_score, sub.chol_score, s=s, c=col, lw=0, alpha=0.6, rasterized=True, label=lab)
    ax.set_xlabel("Hepatocyte score, log2(1 + CP10k)"); ax.set_ylabel("Cholangiocyte score, log2(1 + CP10k)")
    ax.legend(frameon=False, fontsize=6, markerscale=3)
    for sp_ in ["top", "right"]:
        ax.spines[sp_].set_visible(False)
    fig.tight_layout(); fig.savefig(OUT / f"score_plot_{c}.pdf"); plt.close(fig)

with open(OUT / "pip_freeze.txt", "w") as fh:
    import importlib.metadata as im
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500); pd.set_option("display.max_columns", 50)
for n in ["cohort_detectability", "tests", "abs_clean_tests", "control_summary", "geometry_summary", "hep_share_filters_summary", "concentration_summary"]:
    print("\n==", n); print(pd.read_csv(OUT / f"{n}.tsv", sep="\t").round(3).to_string())
print("\n== gene summary depth"); G_ = pd.read_csv(OUT / "gene_summary.tsv", sep="\t"); print(G_[G_["mode"] == "depth"].round(3).to_string())
