#!/usr/bin/env python3
"""V4: Stromal ECM spatial organization in HMSMA (HRA007511) and Vu Visium.

Reuses the byte-pinned paper engine (04_spatial_projection.py) for scoring,
residualization, tissue-island graphs, Moran's I and matched-gene null sets.
Adds per-array null tracking so the statistic can be summarized per donor and
tested across donors. See PRESPEC_V4.md in the output directory.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
ENGINE = ROOT / "Analysis/Multimodal_Program_Projection/scripts/04_spatial_projection.py"
ENGINE_SHA256 = "d198583734ebc02af62973e644963a4d50b87d6de37faa7ef8084b854b94fa94"
HOTSPOT = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
REGISTRY_SHA256 = "b134b5f8e46a716e857af22271aa54e1294432b9f9ab735559d7f603cdd29fe7"
PROGRAM_UID = "hotspot_hepatocytes_f05c535ae5bbc0b9"
GSA_XLSX = ROOT / "data/HRA007511/metadata/HRA007511.xlsx"
HMSMA_H5AD = ROOT / "Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled/HRA007511_merged.h5ad"
C2L = {
    "GSE192741": ROOT / "Analysis/Spatial/results/cell2location/spatial_model/spatial_deconvolved.h5ad",
    "Vu_et_al_2025": ROOT / "Analysis/Spatial/results/cell2location/spatial_model_vu/spatial_deconvolved_vu.h5ad",
}
PAPER_OBSERVED = {"GSE192741": 0.19266400657890104, "Vu_et_al_2025": 0.2768127207447223}
HEP_MARKERS = ["ALB", "APOA1", "APOA2", "APOC3", "APOH", "AHSG", "FGA", "FGB", "FGG",
               "HP", "SERPINA1", "TF", "TTR", "HPX", "RBP4"]
INFERRED_CONTROLS = {"27", "28", "29"}
N_NULL = 9_999
N_SENS = 999
N_SPOT_PERM = 999
N_FLIP = 99_999
N_BOOT = 9_999
SEED_SPOT, SEED_BOOT, SEED_FLIP = 7, 20261007, 20261008
DATASET_INDEX = {"GSE192741": 0, "Vu_et_al_2025": 1, "HRA007511_HMSMA": 2}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(8 << 20):
            h.update(chunk)
    return h.hexdigest()


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_engine():
    if sha256(ENGINE) != ENGINE_SHA256:
        raise RuntimeError("engine hash drift")
    for key in [k for k in os.environ if k.startswith("FIG4_SPATIAL")]:
        os.environ.pop(key)
    os.environ["MASLD_PROJECT_ROOT"] = str(ROOT)
    out_dir = ROOT / "Analysis/Multimodal_Program_Projection/results/spatial"
    if not out_dir.is_dir():
        raise RuntimeError("engine output dir must pre-exist (import-time mkdir must be a no-op)")
    spec = importlib.util.spec_from_file_location("pinned_engine", ENGINE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.N_NULL, mod.N_SENSITIVITY_NULL, mod.SEED = N_NULL, N_SENS, 42
    mod.GRAPH_DISTANCE_MULTIPLIER, mod.MIN_GRAPH_SPOTS, mod.SMOKE = 2.5, 8, False
    return mod


def load_program():
    reg_path = HOTSPOT / "program_registry_v2.tsv"
    if sha256(reg_path) != REGISTRY_SHA256:
        raise RuntimeError("registry hash drift")
    reg = pd.read_csv(reg_path, sep="\t", dtype=str)
    mem = pd.read_csv(HOTSPOT / "program_membership_v2.tsv", sep="\t", dtype=str)
    mem = mem[mem["mapped_symbol_status"] == "gencode_v49_unique_symbol_confirmed"].copy()
    mem["original_l1_weight"] = mem["original_l1_weight"].astype(float)
    robust = set(reg.loc[reg["robust_display"] == "TRUE", "program_uid"])
    union = set(mem.loc[mem["program_uid"].isin(robust), "mapped_symbol"].astype(str))
    prog = mem[mem["program_uid"] == PROGRAM_UID]
    weights = prog.groupby("mapped_symbol")["original_l1_weight"].sum()
    if PROGRAM_UID not in robust:
        raise RuntimeError("Stromal ECM not in robust_display universe")
    return weights, union


def gsa_donor_map() -> pd.DataFrame:
    x = pd.ExcelFile(GSA_XLSX)
    ind, smp, exp, run = (x.parse(s) for s in ("Individual", "Sample", "Experiment", "Run"))
    m = (run[["Accession", "Experiment accession"]].rename(columns={"Accession": "run"})
         .merge(exp[["Accession", "BioSample accession", "Library Construction / Experimental Design"]]
                .rename(columns={"Accession": "Experiment accession", "Library Construction / Experimental Design": "design"}),
                on="Experiment accession")
         .merge(smp[["Accession", "Sample Name", "Individual Accession"]].rename(columns={"Accession": "BioSample accession"}),
                on="BioSample accession")
         .merge(ind[["Accession", "Individual Identifier", "Gender", "Disease Stage"]]
                .rename(columns={"Accession": "Individual Accession"}), on="Individual Accession"))
    st = m["design"].astype(str).str.contains("Visium")
    scrna_ind = set(m.loc[~st, "Individual Identifier"].astype(str))
    st_tab = m[st].drop_duplicates("Sample Name").copy()
    st_tab["array_id"] = "HRA_" + st_tab["Sample Name"].astype(str).str.replace("-ST", "", regex=False)
    st_tab["individual"] = st_tab["Individual Identifier"].astype(str)
    st_tab["has_scrna_in_gsa"] = st_tab["individual"].isin(scrna_ind)
    st_tab["group"] = np.where(st_tab["individual"].isin(INFERRED_CONTROLS), "inferred_control", "primary")
    out = st_tab[["array_id", "Sample Name", "individual", "Gender", "has_scrna_in_gsa", "group"]].rename(
        columns={"Sample Name": "gsa_sample_name", "Gender": "sex"})
    if len(out) != 35 or out["individual"].nunique() != 35 or out["array_id"].nunique() != 35:
        raise RuntimeError("GSA map is not 35 arrays to 35 individuals")
    if set(out.loc[out["has_scrna_in_gsa"], "individual"]) != INFERRED_CONTROLS:
        raise RuntimeError("scRNA-overlap individuals differ from 27/28/29")
    return out.sort_values("array_id").reset_index(drop=True)


def load_c2l(eng, dataset):
    import anndata as ad
    adata = ad.read_h5ad(C2L[dataset])
    adata.var_names = adata.var_names.astype(str)
    obs = adata.obs.copy()
    obs["sample_id"] = obs["sample_id"].astype(str)
    obs["individual"] = obs["individual"].astype(str)
    counts = eng.dense_counts(adata)
    names = [str(x).replace("means_per_cluster_mu_fg_", "") for x in adata.uns["mod"]["factor_names"]]
    hep = np.asarray(adata.obsm["q05_cell_abundance_w_sf"])[:, names.index("Hepatocytes")]
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    return counts, obs, coords, pd.Index(adata.var_names), hep


def load_hmsma(donors: pd.DataFrame):
    import anndata as ad
    adata = ad.read_h5ad(HMSMA_H5AD)
    obs = adata.obs.copy()
    obs["sample_id"] = obs["sample_id"].astype(str)
    dmap = dict(zip(donors["array_id"], donors["individual"]))
    missing = sorted(set(obs["sample_id"]) - set(dmap))
    if missing:
        raise RuntimeError(f"HMSMA arrays without GSA individual: {missing}")
    obs["individual"] = obs["sample_id"].map(dmap)
    obs["condition"] = obs["individual"].map(lambda i: "inferred_control" if i in INFERRED_CONTROLS else "primary")
    x = adata.X.tocsr().astype(np.float32)
    lib = np.asarray(x.sum(axis=1)).ravel()
    ngenes = np.diff(x.indptr)
    obs_lib = pd.to_numeric(obs["total_counts"]).to_numpy(float)
    audit = {"max_abs_total_counts_diff": float(np.max(np.abs(lib - obs_lib))),
             "max_abs_n_genes_diff": float(np.max(np.abs(ngenes - pd.to_numeric(obs["n_genes"]).to_numpy(float))))}
    obs["total_counts"] = lib
    obs["n_genes_by_counts"] = ngenes
    sym = adata.var["gene_symbol"].astype(str).to_numpy()
    uniq, inv = np.unique(sym, return_inverse=True)
    agg = sparse.csr_matrix((np.ones(len(sym), dtype=np.float32), (np.arange(len(sym)), inv)),
                            shape=(len(sym), len(uniq)))
    counts = (x @ agg).tocsr()
    det = np.asarray((counts > 0).mean(axis=0)).ravel()
    keep = det >= 0.01
    counts = counts[:, keep].tocsr()
    audit.update({"n_ensembl": len(sym), "n_symbols": len(uniq), "n_symbols_detected_1pct": int(keep.sum()),
                  "n_spots": int(x.shape[0])})
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    del adata, x
    return counts, obs, coords, pd.Index(uniq[keep]), audit


def per_array_null(eng, norm, sets, weights, design, graphs, arrays, chunk=64):
    out = np.empty((sets.shape[0], len(arrays)), dtype=np.float64)
    for start in range(0, sets.shape[0], chunk):
        block = sets[start:start + chunk]
        uniq, inverse = np.unique(block, return_inverse=True)
        z = eng.extract_z(norm, uniq)
        mapped = inverse.reshape(block.shape)
        scores = np.zeros((norm.shape[0], len(block)), dtype=np.float32)
        for j, w in enumerate(weights):
            scores += z[:, mapped[:, j]] * np.float32(w)
        by = eng.moran_columns(eng.residualize(scores, design), graphs)
        out[start:start + len(block)] = np.column_stack([by[a] for a in arrays])
    return out


def spot_permutation(values, graphs, arrays, rng, n_perm):
    """Within-island spot-label permutation of the residual score."""
    rows = []
    for a in arrays:
        vals, wts, perm_vals = [], [], []
        for idx, left, right in graphs[a]:
            v = values[idx].astype(float)
            v = v - v.mean()
            denom = float(v @ v)
            obs_i = len(idx) * float(np.sum(v[left] * v[right])) / (len(left) * denom)
            P = np.stack([rng.permutation(v) for _ in range(n_perm)], axis=1)
            perm_i = len(idx) * np.sum(P[left] * P[right], axis=0) / (len(left) * denom)
            vals.append(obs_i); wts.append(len(idx)); perm_vals.append(perm_i)
        w = np.asarray(wts, float)
        obs_a = float(np.average(vals, weights=w))
        perm_a = np.average(np.vstack(perm_vals), axis=0, weights=w)
        rows.append({"array_id": a, "spotperm_observed": obs_a, "spotperm_null_mean": float(perm_a.mean()),
                     "spotperm_null_q95": float(np.quantile(perm_a, 0.95)),
                     "spotperm_p": (1 + int(np.sum(perm_a >= obs_a))) / (n_perm + 1)})
    return pd.DataFrame(rows)


def run_dataset(eng, dataset, mode, counts, obs, coords, var_names, c2l_hep, weights_all, union, biotype, outdir):
    tag = f"{dataset}__{mode}"
    log(f"{tag}: spots={counts.shape[0]} genes={counts.shape[1]}")
    library = pd.to_numeric(obs["total_counts"]).to_numpy(float)
    norm = eng.log_normalize(counts, library).tocsc()
    mean = np.asarray(norm.mean(axis=0)).ravel()
    detection = np.asarray((counts > 0).mean(axis=0)).ravel()
    g2i = {g: i for i, g in enumerate(var_names)}

    markers = [g for g in HEP_MARKERS if g in g2i and detection[g2i[g]] >= 0.01 and g not in weights_all.index]
    if mode == "c2l":
        covariate = np.asarray(c2l_hep, float)
    else:
        if len(markers) < 10:
            raise RuntimeError(f"S2 stop: only {len(markers)} hepatocyte markers measured in {dataset}")
        covariate = eng.extract_z(norm, np.array([g2i[g] for g in markers])).mean(axis=1).astype(float)

    pool = eng.matched_pool_table(var_names, mean, detection, biotype)
    pool["lineage_corr"] = eng.abundance_correlation(norm, covariate)
    pool["lineage_corr_bin"] = pd.qcut(pool["lineage_corr"].rank(method="first"), 10, labels=False, duplicates="drop")

    m = weights_all.rename("w").reset_index().rename(columns={"mapped_symbol": "gene_symbol"})
    m["present"] = m["gene_symbol"].isin(g2i)
    m["detected"] = m["gene_symbol"].map(lambda g: detection[g2i[g]] >= 0.01 if g in g2i else False)
    measured = m[m["present"] & m["detected"]].copy()
    retained = float(measured["w"].sum())
    testable = len(measured) >= 8 and retained >= 0.20
    info = {"dataset": dataset, "mode": mode, "n_measured": len(measured), "retained_l1_weight": retained,
            "testable": testable, "n_markers": len(markers), "markers": ",".join(markers)}
    if not testable:
        return info, None, None
    measured["weight"] = measured["w"] / retained
    idx = np.array([g2i[g] for g in measured["gene_symbol"]], dtype=int)
    w = measured["weight"].to_numpy(float)
    z = eng.extract_z(norm, idx)
    score = z @ w
    equal_score = z.mean(axis=1)
    top = int(np.argmax(w))
    lw = np.delete(w, top); lw = lw / lw.sum()
    leave_score = np.delete(z, top, axis=1) @ lw

    design = eng.design_matrix(obs, covariate)
    residual = eng.residualize(score, design)
    graphs, gaudit = eng.section_graphs(obs, coords, k=6)
    gaudit.insert(0, "mode", mode); gaudit.insert(0, "dataset", dataset)
    arrays = sorted(graphs)
    res_by = eng.moran_columns(residual, graphs)
    eq_by = eng.moran_columns(eng.residualize(equal_score, design), graphs)
    lv_by = eng.moran_columns(eng.residualize(leave_score, design), graphs)
    info["observed_donor_collapsed_residual_I"] = float(eng.donor_collapse(res_by, obs, dataset)[0])
    if dataset in PAPER_OBSERVED and mode == "c2l":
        info["paper_observed_residual_I"] = PAPER_OBSERVED[dataset]
        info["reproduces_paper"] = abs(info["observed_donor_collapsed_residual_I"] - PAPER_OBSERVED[dataset]) < 1e-5

    rng = np.random.default_rng(42 + 1000 * DATASET_INDEX[dataset])
    sets, _, _ = eng.build_matched_sets(measured, pool, union, N_NULL, rng)
    info["matched_set_sha256"] = hashlib.sha256(sets.tobytes()).hexdigest()
    log(f"{tag}: observed={info['observed_donor_collapsed_residual_I']:.6f}; running {N_NULL} null sets")
    null = per_array_null(eng, norm, sets, w, design, graphs, arrays)
    sens = sets[:N_SENS]
    null_eq = per_array_null(eng, norm, sens, np.repeat(1.0 / len(w), len(w)), design, graphs, arrays)
    null_lv = per_array_null(eng, norm, np.delete(sens, top, axis=1), lw, design, graphs, arrays)
    np.save(outdir / f"null_per_array__{tag}.npy", null)

    perm = spot_permutation(residual, graphs, arrays, np.random.default_rng(SEED_SPOT), N_SPOT_PERM)
    a2d = obs[["sample_id", "individual", "condition"]].drop_duplicates().set_index("sample_id")
    rows = []
    for j, a in enumerate(arrays):
        I = float(res_by[a][0]); nl = null[:, j]
        rows.append({"dataset": dataset, "mode": mode, "array_id": a,
                     "individual": a2d.loc[a, "individual"], "condition": str(a2d.loc[a, "condition"]),
                     "residual_I": I, "null_mean": float(nl.mean()), "null_sd": float(nl.std(ddof=1)),
                     "null_q95": float(np.quantile(nl, 0.95)), "centered_d": I - float(nl.mean()),
                     "matched_null_p": (1 + int(np.sum(nl >= I))) / (N_NULL + 1),
                     "equal_centered": float(eq_by[a][0]) - float(null_eq[:, j].mean()),
                     "leave_top_centered": float(lv_by[a][0]) - float(null_lv[:, j].mean())})
    arr = pd.DataFrame(rows).merge(perm, on="array_id")
    gene_tab = measured[["gene_symbol", "w", "weight"]].assign(dataset=dataset, mode=mode)
    return info, arr, (null, arrays, gaudit, gene_tab)


def donor_tests(arr: pd.DataFrame, null: np.ndarray, arrays: list, keep_individuals: set, label: str):
    sub = arr[arr["individual"].isin(keep_individuals)]
    cols = [arrays.index(a) for a in sub["array_id"]]
    donor_ids = sub["individual"].to_numpy()
    uniq = sorted(set(donor_ids))
    D = np.zeros((len(uniq), len(cols)))
    for i, d in enumerate(uniq):
        D[i, donor_ids == d] = 1.0 / np.sum(donor_ids == d)
    obs_I = D @ sub["residual_I"].to_numpy()
    null_I = null[:, cols] @ D.T  # n_null x donors
    d = obs_I - null_I.mean(axis=0)
    n = len(d)
    flips = np.random.default_rng(SEED_FLIP).choice([-1.0, 1.0], size=(N_FLIP, n))
    flip_means = flips @ d / n
    boot = np.random.default_rng(SEED_BOOT).integers(0, n, size=(N_BOOT, n))
    boot_means = d[boot].mean(axis=1)
    cohort_obs = float(obs_I.mean()); cohort_null = null_I.mean(axis=1)
    return {"analysis_set": label, "n_donors": n, "n_arrays": len(cols),
            "mean_centered_d": float(d.mean()), "median_centered_d": float(np.median(d)),
            "boot95_low": float(np.quantile(boot_means, 0.025)), "boot95_high": float(np.quantile(boot_means, 0.975)),
            "n_donors_d_positive": int(np.sum(d > 0)),
            "n_donors_above_own_null_q95": int(np.sum(obs_I > np.quantile(null_I, 0.95, axis=0))),
            "signflip_p_one_sided": (1 + int(np.sum(flip_means >= d.sum() / n - 1e-12))) / (N_FLIP + 1),
            "cohort_mean_I": cohort_obs, "cohort_null_mean": float(cohort_null.mean()),
            "cohort_null_sd": float(cohort_null.std(ddof=1)),
            "cohort_centered": cohort_obs - float(cohort_null.mean()),
            "cohort_matched_null_p": (1 + int(np.sum(cohort_null >= cohort_obs))) / (N_NULL + 1),
            "cohort_percentile_of_observed": float(np.mean(cohort_null < cohort_obs))}, pd.DataFrame(
        {"analysis_set": label, "individual": uniq, "donor_residual_I": obs_I, "donor_null_mean": null_I.mean(axis=0),
         "donor_centered_d": d})


def holm(p):
    p = np.asarray(p, float); order = np.argsort(p); m = len(p)
    adj = np.empty(m); running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * p[i])); adj[i] = running
    return adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--smoke", action="store_true", help="tiny null sizes; code check only")
    args = ap.parse_args()
    out = args.outdir
    if args.smoke:
        global N_NULL, N_SENS, N_SPOT_PERM, N_FLIP, N_BOOT
        N_NULL, N_SENS, N_SPOT_PERM, N_FLIP, N_BOOT = 49, 19, 19, 999, 199
        (out / "PRESPEC_V4.md").write_text("smoke\n")
    if not (out / "PRESPEC_V4.md").is_file():
        raise RuntimeError("PRESPEC_V4.md must exist before running")
    if (out / "summary_donor_tests.tsv").exists():
        raise RuntimeError("outputs already exist; refusing to overwrite")
    eng = load_engine()
    weights, union = load_program()
    gm = pd.read_csv(eng.GENE_META_FILE, sep="\t", usecols=["gene_name", "gene_biotype"], compression="gzip")
    gm = gm.drop_duplicates("gene_name")
    biotype = dict(zip(gm["gene_name"].astype(str), gm["gene_biotype"].astype(str)))
    donors = gsa_donor_map()
    donors.to_csv(out / "hmsma_gsa_donor_map.tsv", sep="\t", index=False)

    infos, arrays_all, audits, genes = [], [], [], []
    store = {}
    for dataset in ("GSE192741", "Vu_et_al_2025"):
        counts, obs, coords, var, hep = load_c2l(eng, dataset)
        for mode in ("c2l", "proxy"):
            info, arr, extra = run_dataset(eng, dataset, mode, counts, obs, coords, var, hep, weights, union, biotype, out)
            infos.append(info)
            if arr is not None:
                arrays_all.append(arr); audits.append(extra[2]); genes.append(extra[3]); store[(dataset, mode)] = (arr, extra)
        del counts
    pd.DataFrame(infos).to_csv(out / "run_info.tsv", sep="\t", index=False)
    s1_ok = all(i.get("reproduces_paper", False) for i in infos if i["mode"] == "c2l")
    if not s1_ok:
        (out / "STOP_S1.txt").write_text("Engine did not reproduce the paper observed residual Moran's I; HMSMA not run.\n")
        log("S1 STOP"); return

    counts, obs, coords, var, haudit = load_hmsma(donors)
    json.dump(haudit, open(out / "hmsma_load_audit.json", "w"), indent=2)
    info, arr, extra = run_dataset(eng, "HRA007511_HMSMA", "proxy", counts, obs, coords, var, None, weights, union, biotype, out)
    infos.append(info)
    pd.DataFrame(infos).to_csv(out / "run_info.tsv", sep="\t", index=False)
    if arr is None:
        (out / "STOP_S2.txt").write_text("Stromal ECM untestable in HMSMA.\n"); return
    arrays_all.append(arr); audits.append(extra[2]); genes.append(extra[3]); store[("HRA007511_HMSMA", "proxy")] = (arr, extra)

    pd.concat(arrays_all).to_csv(out / "per_array_results.tsv", sep="\t", index=False)
    pd.concat(audits).to_csv(out / "graph_audit.tsv", sep="\t", index=False)
    pd.concat(genes).to_csv(out / "measured_program_genes.tsv", sep="\t", index=False)

    summaries, donor_rows = [], []
    for (dataset, mode), (arr, extra) in store.items():
        null, arrays = extra[0], extra[1]
        sets = {"all": set(arr["individual"])}
        if dataset == "HRA007511_HMSMA":
            prim = set(arr.loc[arr["condition"] == "primary", "individual"])
            sets = {"primary_32": prim, "all_35": set(arr["individual"])}
        for label, keep in sets.items():
            s, drows = donor_tests(arr, null, arrays, keep, label)
            s.update({"dataset": dataset, "mode": mode}); summaries.append(s)
            drows.insert(0, "mode", mode); drows.insert(0, "dataset", dataset); donor_rows.append(drows)
    # worst-case drop-4 for HMSMA primary
    arr, extra = store[("HRA007511_HMSMA", "proxy")]
    prim = arr[arr["condition"] == "primary"]
    if len(prim) < 20:
        (out / "STOP_S3.txt").write_text("Fewer than 20 primary HMSMA donors with a graph.\n")
    drop = set(prim.sort_values("centered_d", ascending=False)["individual"].head(4))
    s, drows = donor_tests(arr, extra[0], extra[1], set(prim["individual"]) - drop, "primary_drop_top4")
    s.update({"dataset": "HRA007511_HMSMA", "mode": "proxy", "dropped": ",".join(sorted(drop))}); summaries.append(s)
    summ = pd.DataFrame(summaries)
    fam = summ[((summ["dataset"] == "HRA007511_HMSMA") & (summ["analysis_set"] == "primary_32"))
               | ((summ["dataset"] == "Vu_et_al_2025") & (summ["mode"] == "c2l"))].copy()
    fam["holm_p"] = holm(fam["signflip_p_one_sided"].to_numpy())
    summ = summ.merge(fam[["dataset", "mode", "analysis_set", "holm_p"]], how="left", on=["dataset", "mode", "analysis_set"])
    summ["family_member"] = summ["holm_p"].notna()
    summ.to_csv(out / "summary_donor_tests.tsv", sep="\t", index=False)
    pd.concat(donor_rows).to_csv(out / "per_donor_results.tsv", sep="\t", index=False)
    log("done")


if __name__ == "__main__":
    main()
