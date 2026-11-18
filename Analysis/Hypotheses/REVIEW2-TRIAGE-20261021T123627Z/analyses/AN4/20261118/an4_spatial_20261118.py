#!/usr/bin/env python3
"""AN4 (20261118 candidate, N08): All 117 frozen v2 Hotspot programs under all-16
cell2location within-section nonlinear ridge residualization + weighted single-gene mean
Moran's I across 9,999 matched-gene draws.

Also runs focal lineage-matching sensitivity (matching Stromal ECM controls on
Fibroblasts and Ductular injury controls on Cholangiocytes).
"""

from __future__ import annotations

import hashlib
import importlib.util
import os
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


BASE = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
ENGINE_PATH = BASE / "Analysis/Multimodal_Program_Projection/scripts/04_spatial_projection_20261118.py"
OUT = BASE / "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN4/20261118"

FOCAL_ALT_LINEAGE = {
    "hotspot_hepatocytes_f05c535ae5bbc0b9": "fibroblasts",   # Stromal ECM
    "hotspot_hepatocytes_48f39dd4d817a10e": "cholangiocytes",  # Ductular injury
}


def load_engine():
    spec = importlib.util.spec_from_file_location("spatial_engine_20261118", ENGINE_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    eng = load_engine()

    gene_meta = pd.read_csv(eng.GENE_META_FILE, sep="\t", usecols=["gene_name", "gene_biotype", "chromosome"], compression="gzip")
    primary_chr = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}
    gene_meta_prim = gene_meta[gene_meta["chromosome"].isin(primary_chr)].drop_duplicates("gene_name")
    biotype = dict(zip(gene_meta_prim["gene_name"].astype(str), gene_meta_prim["gene_biotype"].astype(str)))

    reg = pd.read_csv(eng.V2_REGISTRY_FILE, sep="\t")
    reg["program_id"] = reg["program_uid"]
    reg["program_name"] = reg["module_name"]
    reg = reg.sort_values(["cell_type", "module"]).reset_index(drop=True)
    reg["display_order"] = np.arange(1, len(reg) + 1)

    mem = pd.read_csv(eng.V2_MEMBERSHIP_FILE, sep="\t")
    mem["program_id"] = mem["program_uid"]
    prim_counts = gene_meta_prim.groupby("gene_name").size().to_dict()
    sym = mem["mapped_symbol"].fillna("").astype(str).str.strip()
    canon = mem["canonical_gene"].fillna("").astype(str).str.strip()
    mem["gene_symbol"] = np.where(sym != "", sym, canon)
    in_prim = (
        (mem["mapped_symbol_status"] == "gencode_v49_unique_symbol_confirmed")
        | (mem["mapped_symbol_status"] == "gencode_v49_unambiguous_ensembl_to_symbol")
        | ((mem["mapped_symbol_status"] == "unmapped_or_ambiguous") & (mem["gene_symbol"].map(lambda g: prim_counts.get(g, 0) == 1)))
    )
    mem = mem[in_prim & (mem["gene_symbol"] != "")].copy()
    union = set(mem["gene_symbol"])

    all_rows, focal_rows = [], []
    for ds_idx, (dataset, path) in enumerate(eng.DATASETS.items()):
        t0 = time.time()
        adata = ad.read_h5ad(path)
        adata.var_names = adata.var_names.astype(str)
        obs = adata.obs.copy()
        obs["sample_id"] = obs["sample_id"].astype(str)
        obs["individual"] = obs["individual"].astype(str)
        counts = eng.dense_counts(adata)
        library = pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(float)
        norm = eng.log_normalize(counts, library)
        mean = np.asarray(norm.mean(axis=0)).ravel()
        detection = np.asarray((counts > 0).mean(axis=0)).ravel()
        base_pool = eng.matched_pool_table(adata.var_names, mean, detection, biotype)

        factor_names = [str(x).replace("means_per_cluster_mu_fg_", "") for x in adata.uns["mod"]["factor_names"]]
        q05_16 = np.asarray(adata.obsm["q05_cell_abundance_w_sf"], dtype=float)
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        graphs, _ = eng.section_graphs(obs, coords, k=6)

        res_16 = eng.compute_within_section_single_gene_moran(norm, obs, graphs, q05_16, extra_cov=None, compute_raw=True)
        bar_res = res_16["bar_I_resid"]
        bar_raw = res_16["bar_I_raw"]
        finite_m = np.isfinite(bar_res) & np.isfinite(bar_raw)
        g2i = {str(g): i for i, g in enumerate(adata.var_names)}

        lineage_pools = {}
        for lin, c2l_name in eng.LINEAGE_TO_C2L.items():
            p_df = base_pool.copy()
            p_df["lineage_corr"] = eng.abundance_correlation(norm, q05_16[:, factor_names.index(c2l_name)])
            p_df["lineage_corr_bin"] = pd.qcut(p_df["lineage_corr"].rank(method="first"), 10, labels=False, duplicates="drop")
            lineage_pools[lin] = p_df

        rng = np.random.default_rng(eng.SEED + 1000 * ds_idx)
        ds_rows = []
        for row in reg.itertuples(index=False):
            pid = str(row.program_id)
            lin = str(row.cell_type)
            m = mem[mem["program_id"] == pid].groupby("gene_symbol", as_index=False)["original_l1_weight"].sum()
            m["present"] = m["gene_symbol"].isin(g2i)
            m["detected"] = m["gene_symbol"].map(lambda g: (detection[g2i[g]] >= 0.01 and bool(finite_m[g2i[g]])) if g in g2i else False)
            measured = m[m["present"] & m["detected"]].copy()
            retained = float(measured["original_l1_weight"].sum())
            testable = len(measured) >= 8 and retained >= 0.20

            rec = {
                "dataset": dataset, "program_id": pid, "display_order": int(row.display_order),
                "cell_type": lin, "module": int(row.module), "program_name": str(row.program_name),
                "n_measured": len(measured), "retained_l1_weight": retained, "testable": testable,
                "n_null": eng.N_NULL,
            }
            if not testable:
                ds_rows.append(rec)
                continue

            measured["weight"] = measured["original_l1_weight"] / retained
            idx = np.array([g2i[g] for g in measured["gene_symbol"]], dtype=int)
            w = measured["weight"].to_numpy(np.float64)
            sets, er, lr = eng.build_matched_sets(measured, lineage_pools[lin], union, eng.N_NULL, rng, finite_mask=finite_m)
            stats = eng._score_weighted_mean_moran(bar_res, bar_raw, None, idx, w, sets, eng.N_SENSITIVITY_NULL)
            stats.pop("null_res")
            stats["excess"] = stats["residual_moran_i"] - stats["residual_null_mean"]
            stats["matched_set_sha256"] = hashlib.sha256(sets.tobytes()).hexdigest()
            ds_rows.append({**rec, **stats})

            if pid in FOCAL_ALT_LINEAGE:
                alt_lin = FOCAL_ALT_LINEAGE[pid]
                rng_alt = np.random.default_rng(eng.SEED + 5000 + 1000 * ds_idx + int(row.module))
                sets_alt, _, _ = eng.build_matched_sets(measured, lineage_pools[alt_lin], union, eng.N_NULL, rng_alt, finite_mask=finite_m)
                st_alt = eng._score_weighted_mean_moran(bar_res, bar_raw, None, idx, w, sets_alt, eng.N_SENSITIVITY_NULL)
                st_alt.pop("null_res")
                focal_rows.append({
                    "dataset": dataset, "program_id": pid, "program_name": str(row.program_name),
                    "match_lineage": lin, "residual_moran_i": stats["residual_moran_i"],
                    "residual_null_mean": stats["residual_null_mean"], "residual_moran_z": stats["residual_moran_z"],
                    "residual_pvalue": stats["residual_pvalue"], "sensitivity_sign_agree": stats["sensitivity_sign_agree"],
                })
                focal_rows.append({
                    "dataset": dataset, "program_id": pid, "program_name": str(row.program_name),
                    "match_lineage": alt_lin, "residual_moran_i": st_alt["residual_moran_i"],
                    "residual_null_mean": st_alt["residual_null_mean"], "residual_moran_z": st_alt["residual_moran_z"],
                    "residual_pvalue": st_alt["residual_pvalue"], "sensitivity_sign_agree": st_alt["sensitivity_sign_agree"],
                })

        ds_df = pd.DataFrame(ds_rows)
        ds_df["residual_padj"] = eng.bh(ds_df["residual_pvalue"])
        ds_df["robust"] = (
            ds_df["testable"].astype("boolean").fillna(False).to_numpy(bool)
            & (ds_df["residual_padj"] < 0.05).fillna(False).to_numpy(bool)
            & ds_df["sensitivity_sign_agree"].astype("boolean").fillna(False).to_numpy(bool)
        )
        all_rows.append(ds_df)
        print(f"[AN4_20261118] {dataset}: testable={int(ds_df['testable'].sum())} q<0.05={int((ds_df['residual_padj']<0.05).sum())} robust={int(ds_df['robust'].sum())} ({time.time()-t0:.1f}s)", flush=True)

    fam_df = pd.concat(all_rows, ignore_index=True)
    fam_df.to_csv(OUT / "family_all_programs_20261118.tsv", sep="\t", index=False)
    focal_df = pd.DataFrame(focal_rows)
    focal_df.to_csv(OUT / "focal_lineage_sensitivity_20261118.tsv", sep="\t", index=False)

    summ = fam_df.groupby(["dataset", "cell_type"], as_index=False).agg(
        n_programs=("program_id", "count"),
        n_testable=("testable", "sum"),
        n_q05=("residual_padj", lambda s: int((s < 0.05).sum())),
        n_robust=("robust", "sum"),
        n_positive_excess=("excess", lambda s: int((s > 0).sum())),
    )
    summ.to_csv(OUT / "family_summary_20261118.tsv", sep="\t", index=False)
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
