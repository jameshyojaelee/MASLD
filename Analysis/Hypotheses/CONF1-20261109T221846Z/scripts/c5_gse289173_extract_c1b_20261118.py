#!/usr/bin/env python3
"""CONF1 C5 & C1b (2026-11-18): extract GSE289173 hepatocyte pseudobulk + H90 sex QC + series labels."""
from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
H5AD = ROOT / "data/public_r4/GSE289173/derived_H90-20261015T120337Z/gse289173_qc_nuclei.h5ad"
SEX_TSV = ROOT / "Analysis/Hypotheses/H90-20261015T120337Z/outputs/D5_PC1_sex.tsv"
SETS_ORIG = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z/out_atlas_state_replication/gene_sets.tsv"
SETS_20261118 = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z/out_atlas_state_replication_20261118/gene_sets.tsv"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c5_extract_20261118"
OUT.mkdir(parents=True, exist_ok=True)

STAGE_ORD = {"ctrl": 0, "MASL": 1, "eMASH": 2, "aMASH": 3}
NONHEP_MARKERS = ["COL1A1", "COL3A1", "DCN", "PECAM1", "STAB2", "PTPRC", "CD163", "LYZ"]
C1B_ANCHOR = "CCL20"
SCORE_GENES = ["ALB", "APOB", "CYP3A4", "KRT7", "KRT19", "SPP1", "EPCAM"]
ZON_PORTAL = ["HAL", "ASS1", "ASL", "CPS1", "PCK1"]
ZON_CENTRAL = ["GLUL", "CYP2E1", "CYP1A2", "OAT", "SLC1A2"]
MIN_POS, K_MATCH, N_RAND = 10, 10, 20
SEED_ANCHOR = 20261109
SEED_SIGNFLIP = 20261008
DET_GATE_MIN_DONORS = 6


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def obs_str(series: pd.Series) -> pd.Series:
    return series.astype(str)


def pct(v: np.ndarray) -> np.ndarray:
    return stats.rankdata(v, method="average") / (len(v) + 1.0)


def signflip_exact_or_mc(v: np.ndarray, seed: int = SEED_SIGNFLIP, n_mc: int = 200000) -> float:
    v = np.asarray(v, dtype=np.float64)
    v = v[np.isfinite(v)]
    n = len(v)
    if n == 0:
        return np.nan
    obs = abs(np.mean(v))
    if n <= 20:
        signs = (((np.arange(1 << n)[:, None] >> np.arange(n)) & 1) * 2 - 1).astype(np.float64)
        null = np.abs((signs @ v) / float(n))
        return float(np.mean(null >= obs - 1e-12))
    rng = np.random.default_rng(seed)
    ge = 0
    chunk = 50000
    for s in range(0, n_mc, chunk):
        m = min(chunk, n_mc - s)
        sg = rng.choice([-1.0, 1.0], size=(m, n))
        ge += int(np.sum(np.abs((sg @ v) / float(n)) >= obs - 1e-12))
    return float((1 + ge) / (n_mc + 1))


def main():
    sets_path = SETS_20261118 if SETS_20261118.exists() else SETS_ORIG
    sets_df = pd.read_csv(sets_path, sep="\t")
    if SETS_ORIG.exists() and SETS_20261118.exists():
        assert sha256_file(SETS_ORIG) == sha256_file(SETS_20261118), "gene_sets.tsv SHA256 mismatch!"
    sets_df.to_csv(OUT / "gene_sets_used.tsv", sep="\t", index=False)
    nfkb_genes = [g.strip() for g in sets_df.loc[sets_df["set"] == "NFKB", "genes"].iloc[0].split(",") if g.strip()]
    sets_dict = {
        r["set"]: [g.strip() for g in str(r["genes"]).split(",") if g.strip()]
        for _, r in sets_df.iterrows()
    }

    sex_df = pd.read_csv(SEX_TSV, sep="\t")
    sex_df["gsm"] = sex_df["gsm"].astype(str)
    sex_df["concordant_bool"] = sex_df["concordant"].astype(str).str.strip().str.lower().isin({"true", "1", "yes"})

    A = ad.read_h5ad(H5AD)
    obs = A.obs.copy()
    obs["gsm"] = obs_str(obs["gsm"])
    obs["title"] = obs_str(obs["title"])
    obs["series_raw"] = obs_str(obs["series"])
    obs["group"] = obs_str(obs["group"])
    obs["sex"] = obs_str(obs["sex"])
    obs["label"] = obs_str(obs["label"])

    hep_mask = (obs["label"] == "Hepatocytes").to_numpy()
    Ah = A[hep_mask].copy()
    oh = Ah.obs.copy().reset_index(drop=True)
    oh["gsm"] = obs_str(oh["gsm"])
    oh["title"] = obs_str(oh["title"])
    oh["series_raw"] = obs_str(oh["series"])
    oh["group"] = obs_str(oh["group"])
    oh["sex"] = obs_str(oh["sex"])

    X = Ah.X.tocsr()
    # Assert raw integer UMI counts
    assert np.allclose(X.data, np.round(X.data)), "Non-integer counts in gse289173_qc_nuclei.h5ad!"
    genes = pd.Index(Ah.var_names.astype(str))
    assert not genes.duplicated().any(), "Duplicate symbols in gse289173_qc_nuclei.h5ad!"

    donors_all = sorted(oh["gsm"].unique())
    pb_cols = {}
    donor_rows = []
    sex_lookup = sex_df.set_index("gsm")

    for gsm in donors_all:
        m = (oh["gsm"] == gsm).to_numpy()
        sub = oh.loc[m]
        title = sub["title"].iloc[0]
        grp = sub["group"].iloc[0]
        sex_rep = sub["sex"].iloc[0]
        ser_raw = sub["series_raw"].iloc[0]
        # Explicit assignment matching: ifelse(grepl("^ctrl_", title), "ctrl_series", "FL_series")
        series = "ctrl_series" if title.startswith("ctrl_") else "FL_series"
        if gsm not in sex_lookup.index:
            raise KeyError(f"Missing {gsm} in D5_PC1_sex.tsv")
        sx = sex_lookup.loc[gsm]
        sex_call = str(sx["sex_call"]).strip().lower()
        sex_conc = bool(sx["concordant_bool"])
        xist_cpm = float(sx["xist_cpm"])
        y_cpm_sum = float(sx["y_cpm_sum"])
        female_rep = 1 if sex_rep.lower() == "female" else 0
        female_corr = 1 if sex_call == "female" else 0

        pb_cols[gsm] = np.asarray(X[m].sum(axis=0)).ravel().astype(np.int64)
        donor_rows.append({
            "donor": gsm,
            "title": title,
            "series_raw": ser_raw,
            "series": series,
            "group": grp,
            "female": female_rep,
            "female_reported": female_rep,
            "female_corrected": female_corr,
            "sex_reported": sex_rep,
            "sex_call": sex_call,
            "xist_cpm": xist_cpm,
            "y_cpm_sum": y_cpm_sum,
            "sex_concordant": sex_conc,
            "exclude_sex_discordant": (not sex_conc),
            "n_hep": int(m.sum()),
            "stage_ord": STAGE_ORD[grp],
        })

    donors_df = pd.DataFrame(donor_rows)
    assert len(donors_df) == 48
    donors_df.to_csv(OUT / "donors.tsv", sep="\t", index=False)
    donors_df[~donors_df["exclude_sex_discordant"]].to_csv(OUT / "donors_sex_concordant.tsv", sep="\t", index=False)
    donors_df.loc[
        donors_df["exclude_sex_discordant"],
        ["donor", "title", "series", "group", "sex_reported", "sex_call", "xist_cpm", "y_cpm_sum", "n_hep", "sex_concordant"],
    ].to_csv(OUT / "sex_qc_summary.tsv", sep="\t", index=False)

    pb_df = pd.DataFrame(pb_cols, index=genes).rename_axis("symbol")
    pb_df.to_csv(OUT / "pseudobulk_hep.tsv.gz", sep="\t", compression="gzip")

    # ---- C1b: CCL20 detection gate & sign-flip check (excluding sex-discordant donor) ----
    gi = {g: i for i, g in enumerate(genes)}
    ccl20_idx = gi[C1B_ANCHOR]
    conc_donors = set(donors_df.loc[~donors_df["exclude_sex_discordant"], "donor"])
    det_rows = []
    for _, r in donors_df.iterrows():
        d = r["donor"]
        m = (oh["gsm"] == d).to_numpy()
        Xd = X[m]
        n_pos = int((Xd[:, ccl20_idx] > 0).sum())
        n_neg = int(len(sub := oh.loc[m]) - n_pos)
        elig = (d in conc_donors) and (n_pos >= MIN_POS) and (n_neg >= 5 * MIN_POS)
        det_rows.append({
            "donor": d,
            "title": r["title"],
            "series": r["series"],
            "group": r["group"],
            "stage": r["stage_ord"],
            "sex_concordant": r["sex_concordant"],
            "n_hep": int(r["n_hep"]),
            "n_ccl20_pos": n_pos,
            "frac_CCL20": n_pos / float(max(int(r["n_hep"]), 1)),
            "eligible": elig,
        })
    det_df = pd.DataFrame(det_rows)
    det_df.to_csv(OUT / "c1b_detection_gate.tsv", sep="\t", index=False)

    elig_dons = det_df.loc[det_df["eligible"], "donor"].tolist()
    n_elig = len(elig_dons)
    if n_elig < DET_GATE_MIN_DONORS:
        pd.DataFrame([{
            "test": "C1b",
            "n_eligible": n_elig,
            "gate_min": DET_GATE_MIN_DONORS,
            "decision": "not_testable_detection_gate",
        }]).to_csv(OUT / "c1b_primary.tsv", sep="\t", index=False)
        pd.DataFrame(columns=["anchor", "random_anchors"]).to_csv(OUT / "random_anchors.tsv", sep="\t", index=False)
    else:
        # Compute DIG1-56 NN excess if testable
        rng = np.random.default_rng(SEED_ANCHOR)
        conc_mask = oh["gsm"].isin(conc_donors).to_numpy()
        Xc = X[conc_mask]
        ohc = oh.loc[conc_mask].reset_index(drop=True)
        det_all = np.asarray((Xc > 0).mean(axis=0)).ravel()
        da = det_all[ccl20_idx]
        setg = set(g for gs in sets_dict.values() for g in gs)
        excl = set(SCORE_GENES) | set(ZON_PORTAL) | set(ZON_CENTRAL) | setg | {C1B_ANCHOR}
        cand = sorted(g for g, dv in zip(genes, det_all) if 0.67 * da <= dv <= 1.5 * da and g not in excl and not g.startswith("ENSG"))
        rands = list(rng.choice(cand, size=min(N_RAND, len(cand)), replace=False))
        pd.DataFrame([{"anchor": C1B_ANCHOR, "random_anchors": ",".join(rands)}]).to_csv(OUT / "random_anchors.tsv", sep="\t", index=False)

        side = sorted(set(SCORE_GENES) | set(ZON_PORTAL) | set(ZON_CENTRAL))
        side = [g for g in side if g in gi]
        gsel_g = sorted(setg & set(genes))
        gsel_idx = np.array([gi[g] for g in gsel_g])
        gi_sub = {g: i for i, g in enumerate(gsel_g)}

        ex_list = []
        for d in elig_dons:
            md = (ohc["gsm"] == d).to_numpy()
            Xd = Xc[md]
            umi = np.asarray(Xd.sum(axis=1)).ravel()
            ngene = np.asarray((Xd > 0).sum(axis=1)).ravel()
            Ls = np.log1p(Xd[:, [gi[g] for g in side]].toarray() * (1e4 / umi[:, None]))
            sd = Ls.std(0)
            sd[sd == 0] = np.nan
            Zs = (Ls - Ls.mean(0)) / sd
            smap = {g: i for i, g in enumerate(side)}
            inj = np.nan_to_num(np.nanmean(Zs[:, [smap[g] for g in SCORE_GENES if g in smap]], 1))
            B = np.column_stack([np.ones(len(umi)), stats.rankdata(np.log(umi)), stats.rankdata(np.log(ngene))])
            inj_r = inj - B @ np.linalg.lstsq(B, inj, rcond=None)[0]
            F = np.column_stack([pct(np.log(umi)), pct(np.log(ngene)), pct(inj_r)])
            Xg = Xd[:, gsel_idx].tocsr()
            tot_d = umi

            lfc_map = {}
            for anc in [C1B_ANCHOR] + rands:
                pos = np.asarray((Xd[:, gi[anc]] > 0).todense()).ravel()
                if pos.sum() < MIN_POS or (~pos).sum() < 5 * MIN_POS:
                    continue
                P = np.where(pos)[0]
                N = np.where(~pos)[0]
                tree = cKDTree(F[N])
                _, nn = tree.query(F[P], k=K_MATCH)
                w = np.bincount(N[nn.ravel()], minlength=len(umi)).astype(float)
                sp = np.asarray(Xg[P].sum(axis=0)).ravel()
                sn = Xg.T @ w
                lfc_map[anc] = np.log2(sp / tot_d[P].sum() * 1e6 + 1) - np.log2(sn / (w @ tot_d) * 1e6 + 1)
            null_vecs = [lfc_map[r] for r in rands if r in lfc_map]
            ex = lfc_map[C1B_ANCHOR] - np.mean(np.vstack(null_vecs), axis=0)
            ex_list.append(ex)

        Ex = np.vstack(ex_list)
        nf_idx = [gi_sub[g] for g in nfkb_genes if g in gi_sub]
        donor_nfkb = np.nanmean(Ex[:, nf_idx], axis=1)
        p_c1b = signflip_exact_or_mc(donor_nfkb)
        pd.DataFrame([{
            "test": "C1b",
            "n_eligible": n_elig,
            "mean_excess_nfkb": float(np.mean(donor_nfkb)),
            "p_signflip": p_c1b,
            "decision": "replicated" if (np.mean(donor_nfkb) > 0 and p_c1b < 0.05) else "not_replicated",
        }]).to_csv(OUT / "c1b_primary.tsv", sep="\t", index=False)

    pip_out = subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True)
    (OUT / "pip_freeze.txt").write_text(pip_out)
    print(f"DONE c5_extract: 48 total ({len(conc_donors)} sex-concordant), C1b eligible={n_elig}", flush=True)


if __name__ == "__main__":
    main()
