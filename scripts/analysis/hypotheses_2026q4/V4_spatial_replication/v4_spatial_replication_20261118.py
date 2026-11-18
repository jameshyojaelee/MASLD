#!/usr/bin/env python3
"""V4 (20261118 candidate, N08): Stromal ECM spatial replication across GSE192741,
Vu_et_al_2025, and HRA007511 (HMSMA) using:
1. Within-section nonlinear covariate residualization ([1, log1p(cov_s), cov_s,
   log1p(total_counts_s), log1p(n_genes_s)], ridge lambda=1e-6, within-section std).
2. Weighted single-gene mean Moran's I across 9,999 matched-gene draws.
3. Exact 2^n sign-flip enumeration whenever n <= 15 (GSE192741 n=4 -> 16,
   Vu n=10 -> 1,024); MC 99,999 sign-flips when n > 15 (HMSMA n=32, 35).
4. Explicit HMSMA control-vs-disease array statistics comparing the 3 inferred
   controls (donors 27, 28, 29) against the 32 primary disease donors via
   exact permutation over all C(35, 3) = 6,545 splits + Welch t + Mann-Whitney.
"""

from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import os
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats


BASE = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
ENGINE_PATH = BASE / "Analysis/Multimodal_Program_Projection/scripts/04_spatial_projection_20261118.py"
PROGRAM_UID = "hotspot_hepatocytes_f05c535ae5bbc0b9"  # Stromal ECM
GSA_XLSX = BASE / "data/HRA007511/metadata/HRA007511.xlsx"
HMSMA_H5AD = BASE / "Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled/HRA007511_merged.h5ad"
OUT = BASE / "Analysis/Hypotheses/V4-20261118"

HEP_MARKERS = [
    "ALB", "APOA1", "APOA2", "APOC3", "APOH", "AHSG", "FGA", "FGB", "FGG",
    "HP", "SERPINA1", "TF", "TTR", "HPX", "RBP4",
]
INFERRED_CONTROLS = {"27", "28", "29"}
N_NULL = 9_999
N_SENS = 999
N_FLIP_MC = 99_999
N_BOOT = 9_999
SEED_BOOT, SEED_FLIP = 20261118, 20261119
DATASET_INDEX = {"GSE192741": 0, "Vu_et_al_2025": 1, "HRA007511_HMSMA": 2}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_engine():
    spec = importlib.util.spec_from_file_location("spatial_engine_20261118", ENGINE_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_program(eng):
    reg = pd.read_csv(eng.V2_REGISTRY_FILE, sep="\t", dtype=str)
    mem = pd.read_csv(eng.V2_MEMBERSHIP_FILE, sep="\t", dtype=str)
    gene_meta = pd.read_csv(eng.GENE_META_FILE, sep="\t", usecols=["gene_name", "chromosome"], compression="gzip")
    primary_chr = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}
    prim_counts = gene_meta[gene_meta["chromosome"].isin(primary_chr)].drop_duplicates("gene_name").groupby("gene_name").size().to_dict()
    sym = mem["mapped_symbol"].fillna("").astype(str).str.strip()
    canon = mem["canonical_gene"].fillna("").astype(str).str.strip()
    mem["gene_symbol"] = np.where(sym != "", sym, canon)
    in_prim = (
        (mem["mapped_symbol_status"] == "gencode_v49_unique_symbol_confirmed")
        | (mem["mapped_symbol_status"] == "gencode_v49_unambiguous_ensembl_to_symbol")
        | ((mem["mapped_symbol_status"] == "unmapped_or_ambiguous") & (mem["gene_symbol"].map(lambda g: prim_counts.get(g, 0) == 1)))
    )
    mem = mem[in_prim & (mem["gene_symbol"] != "")].copy()
    mem["original_l1_weight"] = mem["original_l1_weight"].astype(float)
    robust = set(reg.loc[reg["robust_display"] == "TRUE", "program_uid"])
    union = set(mem.loc[mem["program_uid"].isin(robust), "gene_symbol"].astype(str))
    prog = mem[mem["program_uid"] == PROGRAM_UID]
    weights = prog.groupby("gene_symbol")["original_l1_weight"].sum()
    return weights, union


def gsa_donor_map() -> pd.DataFrame:
    x = pd.ExcelFile(GSA_XLSX)
    ind, smp, exp, run = (x.parse(s) for s in ("Individual", "Sample", "Experiment", "Run"))
    m = (
        run[["Accession", "Experiment accession"]].rename(columns={"Accession": "run"})
        .merge(exp[["Accession", "BioSample accession", "Library Construction / Experimental Design"]]
               .rename(columns={"Accession": "Experiment accession", "Library Construction / Experimental Design": "design"}),
               on="Experiment accession")
        .merge(smp[["Accession", "Sample Name", "Individual Accession"]].rename(columns={"Accession": "BioSample accession"}),
               on="BioSample accession")
        .merge(ind[["Accession", "Individual Identifier", "Gender", "Disease Stage"]]
               .rename(columns={"Accession": "Individual Accession"}), on="Individual Accession")
    )
    st = m["design"].astype(str).str.contains("Visium")
    scrna_ind = set(m.loc[~st, "Individual Identifier"].astype(str))
    st_tab = m[st].drop_duplicates("Sample Name").copy()
    st_tab["array_id"] = "HRA_" + st_tab["Sample Name"].astype(str).str.replace("-ST", "", regex=False)
    st_tab["individual"] = st_tab["Individual Identifier"].astype(str)
    st_tab["has_scrna_in_gsa"] = st_tab["individual"].isin(scrna_ind)
    st_tab["group"] = np.where(st_tab["individual"].isin(INFERRED_CONTROLS), "inferred_control", "primary")
    out = st_tab[["array_id", "Sample Name", "individual", "Gender", "has_scrna_in_gsa", "group"]].rename(
        columns={"Sample Name": "gsa_sample_name", "Gender": "sex"}
    )
    return out.sort_values("array_id").reset_index(drop=True)


def load_c2l(eng, dataset):
    adata = ad.read_h5ad(eng.DATASETS[dataset])
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
    adata = ad.read_h5ad(HMSMA_H5AD)
    obs = adata.obs.copy()
    obs["sample_id"] = obs["sample_id"].astype(str)
    dmap = dict(zip(donors["array_id"], donors["individual"]))
    obs["individual"] = obs["sample_id"].map(dmap)
    obs["condition"] = obs["individual"].map(lambda i: "inferred_control" if i in INFERRED_CONTROLS else "primary")
    x = adata.X.tocsr().astype(np.float32)
    lib = np.asarray(x.sum(axis=1)).ravel()
    ngenes = np.diff(x.indptr)
    obs["total_counts"] = lib
    obs["n_genes_by_counts"] = ngenes
    sym = adata.var["gene_symbol"].astype(str).to_numpy()
    uniq, inv = np.unique(sym, return_inverse=True)
    agg = sparse.csr_matrix((np.ones(len(sym), dtype=np.float32), (np.arange(len(sym)), inv)), shape=(len(sym), len(uniq)))
    counts = (x @ agg).tocsr()
    det = np.asarray((counts > 0).mean(axis=0)).ravel()
    keep = det >= 0.01
    counts = counts[:, keep].tocsr()
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    audit = {"n_spots": int(x.shape[0]), "n_symbols_detected_1pct": int(keep.sum())}
    del adata, x
    return counts, obs, coords, pd.Index(uniq[keep]), audit


def run_dataset(eng, dataset, mode, counts, obs, coords, var_names, c2l_hep, weights_all, union, biotype):
    tag = f"{dataset}__{mode}"
    log(f"{tag}: spots={counts.shape[0]} genes={counts.shape[1]}")
    library = pd.to_numeric(obs["total_counts"]).to_numpy(float)
    norm = eng.log_normalize(counts, library).tocsc()
    mean = np.asarray(norm.mean(axis=0)).ravel()
    detection = np.asarray((counts > 0).mean(axis=0)).ravel()
    g2i = {g: i for i, g in enumerate(var_names)}

    markers = [g for g in HEP_MARKERS if g in g2i and detection[g2i[g]] >= 0.01 and g not in weights_all.index]
    if mode == "c2l":
        covariate = np.maximum(np.asarray(c2l_hep, float), 0.0)
    else:
        # Mean log1p CPM10k across hepatocyte markers (nonnegative for [1, log1p(cov), cov, ...])
        m_idx = np.array([g2i[g] for g in markers], dtype=int)
        covariate = np.maximum(np.asarray(norm[:, m_idx].mean(axis=1)).ravel().astype(float), 0.0)

    pool = eng.matched_pool_table(var_names, mean, detection, biotype)
    pool["lineage_corr"] = eng.abundance_correlation(norm, covariate)
    pool["lineage_corr_bin"] = pd.qcut(pool["lineage_corr"].rank(method="first"), 10, labels=False, duplicates="drop")

    graphs, gaudit = eng.section_graphs(obs, coords, k=6)
    gaudit.insert(0, "mode", mode); gaudit.insert(0, "dataset", dataset)
    res_sg = eng.compute_within_section_single_gene_moran(norm, obs, graphs, covariate, compute_raw=True)
    bar_res = res_sg["bar_I_resid"]
    finite_m = np.isfinite(bar_res)

    m = weights_all.rename("w").reset_index().rename(columns={"gene_symbol": "gene_symbol"})
    m["present"] = m["gene_symbol"].isin(g2i)
    m["detected"] = m["gene_symbol"].map(lambda g: (detection[g2i[g]] >= 0.01 and bool(finite_m[g2i[g]])) if g in g2i else False)
    measured = m[m["present"] & m["detected"]].copy()
    retained = float(measured["w"].sum())
    testable = len(measured) >= 8 and retained >= 0.20
    info = {"dataset": dataset, "mode": mode, "n_measured": len(measured), "retained_l1_weight": retained, "testable": testable}
    if not testable:
        return info, None, None

    measured["weight"] = measured["w"] / retained
    idx = np.array([g2i[g] for g in measured["gene_symbol"]], dtype=int)
    w = measured["weight"].to_numpy(np.float64)
    rng = np.random.default_rng(42 + 1000 * DATASET_INDEX[dataset])
    sets, _, _ = eng.build_matched_sets(measured, pool, union, N_NULL, rng, finite_mask=finite_m)
    info["matched_set_sha256"] = hashlib.sha256(sets.tobytes()).hexdigest()

    arrays = res_sg["arrays"]
    obs_arrays = np.array([eng.nan_weighted_dot(res_sg["sec_resid"][a][idx], w) for a in arrays], dtype=np.float64)
    null_per_array = np.column_stack([eng.nan_weighted_mat(res_sg["sec_resid"][a][sets], w) for a in arrays])  # (N_NULL, n_arrays)

    # Sensitivity schemes (cohort mean across donors)
    top = int(np.argmax(w))
    eq_w = np.full(len(w), 1.0 / len(w))
    lw = np.delete(w, top); lw /= lw.sum()
    idx_lv = np.delete(idx, top)
    sens = sets[:N_SENS]
    eq_obs_a = np.array([eng.nan_weighted_dot(res_sg["sec_resid"][a][idx], eq_w) for a in arrays])
    lv_obs_a = np.array([eng.nan_weighted_dot(res_sg["sec_resid"][a][idx_lv], lw) for a in arrays])
    eq_null_a = np.column_stack([eng.nan_weighted_mat(res_sg["sec_resid"][a][sens], eq_w) for a in arrays])
    lv_null_a = np.column_stack([eng.nan_weighted_mat(res_sg["sec_resid"][a][np.delete(sens, top, axis=1)], lw) for a in arrays])

    a2d = obs[["sample_id", "individual", "condition"]].drop_duplicates().set_index("sample_id")
    rows = []
    for j, a in enumerate(arrays):
        I = float(obs_arrays[j]); nl = null_per_array[:, j]; nl = nl[np.isfinite(nl)]
        rows.append({
            "dataset": dataset, "mode": mode, "array_id": a,
            "individual": str(a2d.loc[a, "individual"]), "condition": str(a2d.loc[a, "condition"]),
            "residual_I": I, "null_mean": float(nl.mean()), "null_sd": float(nl.std(ddof=1)),
            "null_q95": float(np.quantile(nl, 0.95)), "centered_d": I - float(nl.mean()),
            "matched_null_p": (1 + int(np.sum(nl >= I))) / (len(nl) + 1),
            "equal_centered": float(eq_obs_a[j] - np.nanmean(eq_null_a[:, j])),
            "leave_top_centered": float(lv_obs_a[j] - np.nanmean(lv_null_a[:, j])),
        })
    arr = pd.DataFrame(rows)
    info["observed_donor_collapsed_residual_I"] = float(np.dot(w, bar_res[idx]))
    info["cohort_null_mean"] = float(np.mean(bar_res[sets] @ w))
    log(f"{tag}: obs={info['observed_donor_collapsed_residual_I']:.5f} null_mean={info['cohort_null_mean']:.5f}")
    gene_tab = measured[["gene_symbol", "w", "weight"]].assign(dataset=dataset, mode=mode)
    return info, arr, (null_per_array, arrays, gaudit, gene_tab)


def donor_tests(arr: pd.DataFrame, null: np.ndarray, arrays: list, keep_individuals: set, label: str):
    sub = arr[arr["individual"].isin(keep_individuals)]
    cols = [arrays.index(a) for a in sub["array_id"]]
    donor_ids = sub["individual"].to_numpy()
    uniq = sorted(set(donor_ids))
    obs_I = np.array([np.nanmean(sub["residual_I"].to_numpy()[donor_ids == u]) for u in uniq], dtype=float)
    null_I = np.column_stack([
        np.nanmean(null[:, [cols[k] for k in range(len(cols)) if donor_ids[k] == u]], axis=1)
        for u in uniq
    ])
    d = obs_I - np.nanmean(null_I, axis=0)
    n = len(d)
    obs_sum = float(d.sum())

    if n <= 15:
        # Exact 2^n sign-flip enumeration
        n_total = 1 << n
        flips = np.array(list(itertools.product([-1.0, 1.0], repeat=n)), dtype=np.float64)
        flip_sums = flips @ d
        # Exact randomized/discrete count (non-strict >= obs_sum)
        p_one = float(np.sum(flip_sums >= obs_sum - 1e-12)) / n_total
        p_two = float(np.sum(np.abs(flip_sums) >= abs(obs_sum) - 1e-12)) / n_total
        flip_method = f"exact_2^{n}_{n_total}"
        n_flips_used = n_total
    else:
        flips = np.random.default_rng(SEED_FLIP).choice([-1.0, 1.0], size=(N_FLIP_MC, n))
        flip_sums = flips @ d
        p_one = (1 + int(np.sum(flip_sums >= obs_sum - 1e-12))) / (N_FLIP_MC + 1)
        p_two = (1 + int(np.sum(np.abs(flip_sums) >= abs(obs_sum) - 1e-12))) / (N_FLIP_MC + 1)
        flip_method = f"mc_{N_FLIP_MC}"
        n_flips_used = N_FLIP_MC

    boot = np.random.default_rng(SEED_BOOT).integers(0, n, size=(N_BOOT, n))
    boot_means = d[boot].mean(axis=1)
    cohort_obs = float(obs_I.mean())
    cohort_null = null_I.mean(axis=1)
    return {
        "analysis_set": label, "n_donors": n, "n_arrays": len(cols),
        "mean_centered_d": float(d.mean()), "median_centered_d": float(np.median(d)),
        "boot95_low": float(np.quantile(boot_means, 0.025)), "boot95_high": float(np.quantile(boot_means, 0.975)),
        "n_donors_d_positive": int(np.sum(d > 0)),
        "n_donors_above_own_null_q95": int(np.sum(obs_I > np.quantile(null_I, 0.95, axis=0))),
        "signflip_method": flip_method, "n_signflips": n_flips_used,
        "signflip_p_one_sided": p_one, "signflip_p_two_sided": p_two,
        "cohort_mean_I": cohort_obs, "cohort_null_mean": float(cohort_null.mean()),
        "cohort_null_sd": float(cohort_null.std(ddof=1)),
        "cohort_centered": cohort_obs - float(cohort_null.mean()),
        "cohort_matched_null_p": (1 + int(np.sum(cohort_null >= cohort_obs))) / (N_NULL + 1),
    }, pd.DataFrame({
        "analysis_set": label, "individual": uniq,
        "donor_residual_I": obs_I, "donor_null_mean": null_I.mean(axis=0), "donor_centered_d": d,
    })


def hmsma_control_vs_disease(arr: pd.DataFrame) -> pd.DataFrame:
    """Compare 3 inferred_control donors (27,28,29) vs 32 primary disease donors in HMSMA."""
    rows = []
    ctrl = arr[arr["condition"] == "inferred_control"]
    prim = arr[arr["condition"] == "primary"]
    n_c, n_p = len(ctrl), len(prim)
    all_idx = np.arange(n_c + n_p)
    comb_ctrl = list(itertools.combinations(all_idx, n_c))  # C(35,3) = 6,545 exact splits
    for metric in ("centered_d", "residual_I"):
        xc = ctrl[metric].to_numpy(float)
        xp = prim[metric].to_numpy(float)
        combined = np.concatenate([xc, xp])
        obs_diff = float(xp.mean() - xc.mean())
        total_sum = float(combined.sum())
        perm_diffs = []
        for idx_c in comb_ctrl:
            sum_c = float(combined[list(idx_c)].sum())
            sum_p = total_sum - sum_c
            perm_diffs.append(sum_p / n_p - sum_c / n_c)
        perm_diffs = np.asarray(perm_diffs, dtype=float)
        p_perm_two = float(np.mean(np.abs(perm_diffs) >= abs(obs_diff) - 1e-12))
        p_perm_one = float(np.mean(perm_diffs >= obs_diff - 1e-12))
        t_w, p_welch = stats.ttest_ind(xp, xc, equal_var=False)
        u_stat, p_mwu = stats.mannwhitneyu(xp, xc, alternative="two-sided")
        rows.append({
            "metric": metric,
            "n_primary": n_p, "n_inferred_control": n_c,
            "control_donors": ",".join(sorted(ctrl["individual"].astype(str))),
            "primary_mean": float(xp.mean()), "primary_median": float(np.median(xp)), "primary_sd": float(xp.std(ddof=1)),
            "control_mean": float(xc.mean()), "control_median": float(np.median(xc)), "control_sd": float(xc.std(ddof=1)),
            "diff_primary_minus_control": obs_diff,
            "exact_n_splits": len(comb_ctrl),
            "exact_perm_p_two_sided": p_perm_two,
            "exact_perm_p_one_sided_primary_gt_ctrl": p_perm_one,
            "welch_t": float(t_w), "welch_p_two_sided": float(p_welch),
            "mannwhitney_u": float(u_stat), "mannwhitney_p_two_sided": float(p_mwu),
        })
    return pd.DataFrame(rows)


def holm(p):
    p = np.asarray(p, float); order = np.argsort(p); m = len(p)
    adj = np.empty(m); running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * p[i])); adj[i] = running
    return adj


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    eng = load_engine()
    weights, union = load_program(eng)
    gm = pd.read_csv(eng.GENE_META_FILE, sep="\t", usecols=["gene_name", "gene_biotype", "chromosome"], compression="gzip")
    primary_chr = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}
    gm = gm[gm["chromosome"].isin(primary_chr)].drop_duplicates("gene_name")
    biotype = dict(zip(gm["gene_name"].astype(str), gm["gene_biotype"].astype(str)))
    donors = gsa_donor_map()
    donors.to_csv(OUT / "hmsma_gsa_donor_map.tsv", sep="\t", index=False)

    infos, arrays_all, audits, genes = [], [], [], []
    store = {}
    for dataset in ("GSE192741", "Vu_et_al_2025"):
        counts, obs, coords, var, hep = load_c2l(eng, dataset)
        for mode in ("c2l", "proxy"):
            info, arr, extra = run_dataset(eng, dataset, mode, counts, obs, coords, var, hep, weights, union, biotype)
            infos.append(info)
            if arr is not None:
                arrays_all.append(arr); audits.append(extra[2]); genes.append(extra[3]); store[(dataset, mode)] = (arr, extra)
        del counts

    counts, obs, coords, var, haudit = load_hmsma(donors)
    json.dump(haudit, open(OUT / "hmsma_load_audit.json", "w"), indent=2)
    info, arr, extra = run_dataset(eng, "HRA007511_HMSMA", "proxy", counts, obs, coords, var, None, weights, union, biotype)
    infos.append(info)
    pd.DataFrame(infos).to_csv(OUT / "run_info.tsv", sep="\t", index=False)
    arrays_all.append(arr); audits.append(extra[2]); genes.append(extra[3]); store[("HRA007511_HMSMA", "proxy")] = (arr, extra)

    pd.concat(arrays_all).to_csv(OUT / "per_array_results.tsv", sep="\t", index=False)
    pd.concat(audits).to_csv(OUT / "graph_audit.tsv", sep="\t", index=False)
    pd.concat(genes).to_csv(OUT / "measured_program_genes.tsv", sep="\t", index=False)

    # Control-vs-disease array stats in HMSMA
    ctrl_vs_dis = hmsma_control_vs_disease(arr)
    ctrl_vs_dis.to_csv(OUT / "hmsma_control_vs_disease_20261118.tsv", sep="\t", index=False)

    summaries, donor_rows = [], []
    for (dataset, mode), (arr_d, extra_d) in store.items():
        null, arrays = extra_d[0], extra_d[1]
        sets_map = {"all": set(arr_d["individual"])}
        if dataset == "HRA007511_HMSMA":
            prim = set(arr_d.loc[arr_d["condition"] == "primary", "individual"])
            sets_map = {"primary_32": prim, "all_35": set(arr_d["individual"])}
        for label, keep in sets_map.items():
            s, drows = donor_tests(arr_d, null, arrays, keep, label)
            s.update({"dataset": dataset, "mode": mode}); summaries.append(s)
            drows.insert(0, "mode", mode); drows.insert(0, "dataset", dataset); donor_rows.append(drows)

    prim_df = arr[arr["condition"] == "primary"]
    drop = set(prim_df.sort_values("centered_d", ascending=False)["individual"].head(4))
    s, drows = donor_tests(arr, extra[0], extra[1], set(prim_df["individual"]) - drop, "primary_drop_top4")
    s.update({"dataset": "HRA007511_HMSMA", "mode": "proxy", "dropped": ",".join(sorted(drop))}); summaries.append(s)
    summ = pd.DataFrame(summaries)
    fam = summ[((summ["dataset"] == "HRA007511_HMSMA") & (summ["analysis_set"] == "primary_32"))
               | ((summ["dataset"] == "Vu_et_al_2025") & (summ["mode"] == "c2l"))].copy()
    fam["holm_p"] = holm(fam["signflip_p_one_sided"].to_numpy())
    summ = summ.merge(fam[["dataset", "mode", "analysis_set", "holm_p"]], how="left", on=["dataset", "mode", "analysis_set"])
    summ["family_member"] = summ["holm_p"].notna()
    summ.to_csv(OUT / "summary_donor_tests.tsv", sep="\t", index=False)
    pd.concat(donor_rows).to_csv(OUT / "per_donor_results.tsv", sep="\t", index=False)
    print("=== V4 SUMMARY ===")
    print(summ[["dataset", "mode", "analysis_set", "n_donors", "mean_centered_d", "signflip_method", "signflip_p_one_sided", "holm_p"]].to_string(index=False))
    print("=== HMSMA CTRL VS DISEASE ===")
    print(ctrl_vs_dis.to_string(index=False))
    log("V4_20261118 done")


if __name__ == "__main__":
    main()
