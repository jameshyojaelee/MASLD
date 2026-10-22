#!/usr/bin/env python3
"""AN4 (REVIEW2 triage, PLAN section 4 A4): spatial sensitivity for review item B-S5.

Uses the byte-pinned Visium engine (04_spatial_projection.py, loaded and hash
checked through spatial_context_v2/visium_rerun_lib.py) and the all-16
cell2location design built exactly as in
spatial_context_v2/31_run_full_composition_sensitivity.py (load_full_designs).

Modes
-----
validate_lineage : (a) rerun the two-program all-16 primary (hepatocyte-lineage
    matched sets, seeds 42/1042, exclusion = union of the two programs) and
    require the matched-set SHA256, observed residual Moran's I, null mean and
    P to equal the r2 sensitivity outputs; (b) rerun the two focal programs with
    control genes matched on correlation with their content lineage
    (fibroblasts for Stromal ECM, cholangiocytes for Ductular injury).
family : all 117 frozen v2 Hotspot programs in one source, engine semantics:
    one RNG per source (seed 42 GSE192741, 1042 Vu), programs in registry
    order, exclusion = union of all program genes, each program matched on its
    own discovery lineage (tcells -> 'T cells' factor). Matched sets are built
    for every program (so the RNG stream does not depend on slicing); the
    residual null is computed only for the programs in the requested slice.

Only the residual statistic (all-16 + log library + log detected genes) is
computed for the nulls, plus the equal-weight and leave-top-gene sensitivity
nulls (999 draws), which is what the primary support rule uses.
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

PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
V2DIR = PROJECT / "Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
sys.path.insert(0, str(V2DIR))
import visium_rerun_lib as lib  # noqa: E402

R2 = PROJECT / (
    "Analysis/Multimodal_Program_Projection/candidates/"
    "spatial-full-composition-sensitivity-candidate-2026-08-11-r2"
)
FOCAL = {
    "hotspot_hepatocytes_f05c535ae5bbc0b9": "fibroblasts",  # Stromal ECM (IGFBP7), hepatocytes::8
    "hotspot_hepatocytes_48f39dd4d817a10e": "cholangiocytes",  # Ductular injury (BICC1), hepatocytes::20
}
LINEAGE_TO_C2L = {
    "hepatocytes": "Hepatocytes",
    "fibroblasts": "Fibroblasts",
    "macrophages": "Macrophages",
    "cholangiocytes": "Cholangiocytes",
    "tcells": "T cells",
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_all16_module():
    spec = importlib.util.spec_from_file_location(
        "_all16", V2DIR / "31_run_full_composition_sensitivity.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def full_registry(paths) -> tuple[pd.DataFrame, pd.DataFrame]:
    """All 117 frozen v2 programs, ordered and filtered exactly as v2_universe /
    selected_v2_membership_rows, without the robust_display restriction."""
    lib.verify_hotspot_ready(paths)
    reg = pd.read_csv(paths.hotspot_root / "program_registry_v2.tsv", sep="\t")
    # v2_universe sorts on (cell_type, int(module), program_uid)
    reg["module"] = reg["module"].astype(int)
    reg = reg.sort_values(["cell_type", "module", "program_uid"]).reset_index(drop=True)
    registry = pd.DataFrame(
        {
            "program_id": reg["program_uid"],
            "display_order": np.arange(1, len(reg) + 1),
            "cell_type": reg["cell_type"],
            "module": reg["module"],
            "program_name": reg["module_name"],
            "robust_display": reg["robust_display"].astype(str).str.upper() == "TRUE",
            "n_source_genes": reg["n_source_genes"],
        }
    )
    mem = pd.read_csv(paths.hotspot_root / "program_membership_v2.tsv", sep="\t")
    mem = mem[mem["mapped_symbol_status"] == "gencode_v49_unique_symbol_confirmed"].copy()
    membership = pd.DataFrame(
        {
            "program_id": mem["program_uid"],
            "cell_type": mem["cell_type"],
            "module": mem["module"].astype(int),
            "source_gene": mem["source_gene"],
            "gene_symbol": mem["mapped_symbol"].astype(str).str.strip(),
            "mapped_symbol": True,
            "original_l1_weight": mem["original_l1_weight"].astype(float),
        }
    )
    return registry, membership


class Source:
    """Per-source state built with the pinned engine's functions."""

    def __init__(self, engine, dataset: str, path: Path, design: np.ndarray, design_obs_index, biotype):
        import anndata as ad

        log(f"loading {dataset}")
        adata = ad.read_h5ad(path)
        adata.var_names = adata.var_names.astype(str)
        obs = adata.obs.copy()
        if not obs.index.equals(design_obs_index):
            raise RuntimeError(f"{dataset}: obs order differs from all-16 design")
        obs["sample_id"] = obs["sample_id"].astype(str)
        obs["individual"] = obs["individual"].astype(str)
        counts = engine.dense_counts(adata)
        library = pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(float)
        self.norm = engine.log_normalize(counts, library)
        mean = np.asarray(self.norm.mean(axis=0)).ravel()
        self.detection = np.asarray((counts > 0).mean(axis=0)).ravel()
        self.pool = engine.matched_pool_table(adata.var_names, mean, self.detection, biotype)
        self.gene_to_index = {g: i for i, g in enumerate(adata.var_names)}
        factor_names = [
            str(x).replace("means_per_cluster_mu_fg_", "") for x in adata.uns["mod"]["factor_names"]
        ]
        abundance_mat = np.asarray(adata.obsm["q05_cell_abundance_w_sf"])
        self.lineage_pools = {}
        self.lineage_corr = {}
        for lineage, label in LINEAGE_TO_C2L.items():
            lp = self.pool.copy()
            lp["lineage_corr"] = engine.abundance_correlation(
                self.norm, abundance_mat[:, factor_names.index(label)]
            )
            lp["lineage_corr_bin"] = pd.qcut(
                lp["lineage_corr"].rank(method="first"), 10, labels=False, duplicates="drop"
            )
            self.lineage_pools[lineage] = lp
            self.lineage_corr[lineage] = lp["lineage_corr"].to_numpy(float)
        coords = np.asarray(adata.obsm["spatial"])
        self.graphs, self.graph_audit = engine.section_graphs(obs, coords, k=6)
        missing = sorted(set(obs["sample_id"]).difference(self.graphs))
        if missing:
            raise RuntimeError(f"{dataset}: arrays without graph {missing}")
        self.obs = obs
        self.design = design
        self.dataset = dataset
        del adata, counts

    def measured(self, membership: pd.DataFrame, pid: str) -> pd.DataFrame:
        m = membership[membership["program_id"] == pid].copy()
        m = m.groupby("gene_symbol", as_index=False)["original_l1_weight"].sum()
        m["present"] = m["gene_symbol"].isin(self.gene_to_index)
        m["detected"] = m["gene_symbol"].map(
            lambda g: self.detection[self.gene_to_index[g]] >= 0.01 if g in self.gene_to_index else False
        )
        return m[m["present"] & m["detected"]].copy()


def observed_stats(engine, src: Source, measured: pd.DataFrame):
    retained = float(measured["original_l1_weight"].sum())
    measured = measured.copy()
    measured["weight"] = measured["original_l1_weight"] / retained
    indices = np.array([src.gene_to_index[g] for g in measured["gene_symbol"]], dtype=int)
    weights = measured["weight"].to_numpy(float)
    z = engine.extract_z(src.norm, indices)
    score = z @ weights
    equal_score = np.mean(z, axis=1)
    top = int(np.argmax(weights))
    leave_weights = np.delete(weights, top)
    leave_weights = leave_weights / leave_weights.sum()
    leave_score = np.delete(z, top, axis=1) @ leave_weights
    res_by = engine.moran_columns(engine.residualize(score, src.design), src.graphs)
    eq_by = engine.moran_columns(engine.residualize(equal_score, src.design), src.graphs)
    lv_by = engine.moran_columns(engine.residualize(leave_score, src.design), src.graphs)
    out = {
        "residual_moran_i": float(engine.donor_collapse(res_by, src.obs, src.dataset)[0]),
        "equal_residual_moran_i": float(engine.donor_collapse(eq_by, src.obs, src.dataset)[0]),
        "leave_top_residual_moran_i": float(engine.donor_collapse(lv_by, src.obs, src.dataset)[0]),
    }
    unit = {s: float(v[0]) for s, v in res_by.items()}
    return out, weights, indices, top, leave_weights, unit


def lineage_weighted_corr(src: Source, idx: np.ndarray, weights: np.ndarray) -> dict:
    return {f"wcorr_{lin}": float(np.sum(src.lineage_corr[lin][idx] * weights)) for lin in LINEAGE_TO_C2L}


def run_program(engine, src: Source, row, measured, sets, n_sens: int):
    obs_stats, weights, indices, top, leave_weights, unit = observed_stats(engine, src, measured)
    null_res = engine.null_residual_moran(
        src.norm, sets, weights, src.design, src.graphs, src.obs, src.dataset
    )
    sens = sets[: min(n_sens, len(sets))]
    eq_null = engine.null_residual_moran(
        src.norm, sens, np.repeat(1.0 / len(weights), len(weights)), src.design, src.graphs, src.obs, src.dataset
    )
    lv_null = engine.null_residual_moran(
        src.norm, np.delete(sens, top, axis=1), leave_weights, src.design, src.graphs, src.obs, src.dataset
    )
    z, p, mu, sd = engine.standardized_against_null(obs_stats["residual_moran_i"], null_res)
    # per-unit excess over each unit's own matched-null mean (as Fig 5F unit panel)
    unit_null = {}
    rec = {
        **obs_stats,
        "residual_null_mean": mu,
        "residual_null_sd": sd,
        "residual_moran_z": z,
        "residual_pvalue": p,
        "n_null_ge_observed": int(np.sum(np.asarray(null_res) >= obs_stats["residual_moran_i"])),
        "excess": obs_stats["residual_moran_i"] - mu,
        "equal_residual_null_mean": float(np.nanmean(eq_null)),
        "leave_top_residual_null_mean": float(np.nanmean(lv_null)),
        "matched_set_sha256": hashlib.sha256(sets.tobytes()).hexdigest(),
    }
    rec["sensitivity_sign_agree"] = bool(
        np.isfinite(z)
        and np.sign(z) == np.sign(rec["equal_residual_moran_i"] - rec["equal_residual_null_mean"])
        and np.sign(z) == np.sign(rec["leave_top_residual_moran_i"] - rec["leave_top_residual_null_mean"])
    )
    prog = lineage_weighted_corr(src, indices, weights)
    rec.update({f"program_{k}": v for k, v in prog.items()})
    # null-set lineage correlations (weighted, same weights)
    for lin in LINEAGE_TO_C2L:
        vals = src.lineage_corr[lin][sets] @ weights
        rec[f"null_wcorr_{lin}_mean"] = float(np.mean(vals))
        rec[f"null_wcorr_{lin}_q025"] = float(np.quantile(vals, 0.025))
        rec[f"null_wcorr_{lin}_q975"] = float(np.quantile(vals, 0.975))
    return rec, np.asarray(null_res, dtype=float), unit


def per_unit_null(engine, src: Source, sets, weights):
    """Per-array residual Moran null means (no donor collapse) for unit-level excess."""
    sums = {}
    n = 0
    for start in range(0, sets.shape[0], 64):
        block = sets[start : start + 64]
        uniq, inverse = np.unique(block, return_inverse=True)
        z = engine.extract_z(src.norm, uniq)
        mapped = inverse.reshape(block.shape)
        scores = np.zeros((src.norm.shape[0], len(block)), dtype=np.float32)
        for j, w in enumerate(weights):
            scores += z[:, mapped[:, j]] * np.float32(w)
        by = engine.moran_columns(engine.residualize(scores, src.design), src.graphs)
        for s, v in by.items():
            sums[s] = sums.get(s, 0.0) + float(np.sum(v))
        n += len(block)
    return {s: v / n for s, v in sums.items()}


def balanced_slices(registry: pd.DataFrame, membership: pd.DataFrame, k: int) -> dict[str, int]:
    sizes = membership.groupby("program_id").size()
    order = registry.assign(size=registry["program_id"].map(sizes).fillna(0)).sort_values(
        ["size", "display_order"], ascending=[False, True]
    )
    loads = [0.0] * k
    assign = {}
    for r in order.itertuples(index=False):
        j = int(np.argmin(loads))
        assign[r.program_id] = j
        loads[j] += float(r.size)
    return assign


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["validate_lineage", "family"], required=True)
    ap.add_argument("--dataset", choices=["GSE192741", "Vu_et_al_2025", "both"], default="both")
    ap.add_argument("--slice", type=int, default=0)
    ap.add_argument("--n-slices", type=int, default=1)
    ap.add_argument("--n-null", type=int, default=lib.N_NULL)
    ap.add_argument("--n-sens", type=int, default=lib.N_SENSITIVITY_NULL)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)

    paths = lib.build_paths(PROJECT)
    engine = lib.load_legacy_engine(paths)  # verifies pinned engine SHA256
    all16 = load_all16_module()
    designs, factor_names = all16.load_full_designs(paths)
    gm = pd.read_csv(paths.gene_metadata, sep="\t", usecols=["gene_name", "gene_biotype"], compression="gzip")
    gm = gm.drop_duplicates("gene_name")
    biotype = dict(zip(gm["gene_name"].astype(str), gm["gene_biotype"].astype(str)))
    datasets = list(paths.datasets) if args.dataset == "both" else [args.dataset]
    meta = {
        "mode": args.mode,
        "datasets": datasets,
        "n_null": args.n_null,
        "n_sens": args.n_sens,
        "slice": args.slice,
        "n_slices": args.n_slices,
        "engine_sha256": lib.sha256_file(paths.v1_engine),
        "all16_factor_order": factor_names,
        "seed_rule": "np.random.default_rng(42 + (0 if GSE192741 else 1000)), sequential over programs in registry order",
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", "none"),
    }
    (args.out / "run_meta.json").write_text(json.dumps(meta, indent=2))

    if args.mode == "validate_lineage":
        reg2, mem2, _ = lib.prepare_engine_inputs(paths, "v2")
        r2 = pd.read_csv(R2 / "full_composition_program_results.tsv", sep="\t")
        rows, checks, units = [], [], []
        for ds in datasets:
            src = Source(engine, ds, paths.datasets[ds], designs[ds]["design"], designs[ds]["obs"].index, biotype)
            union2 = set(mem2["gene_symbol"].dropna().astype(str))
            for scheme in ("hepatocyte_primary_repro", "content_lineage"):
                rng = np.random.default_rng(42 + (0 if ds == "GSE192741" else 1000))
                for row in reg2.itertuples(index=False):
                    measured = src.measured(mem2, row.program_id)
                    lineage = row.cell_type if scheme == "hepatocyte_primary_repro" else FOCAL[row.program_id]
                    sets, er, lr = engine.build_matched_sets(
                        measured, src.lineage_pools[lineage], union2, args.n_null, rng
                    )
                    log(f"{ds} {scheme} {row.program_name}: n={len(measured)} lineage={lineage}")
                    rec, null, unit = run_program(engine, src, row, measured, sets, args.n_sens)
                    w = (measured["original_l1_weight"] / measured["original_l1_weight"].sum()).to_numpy(float)
                    unull = per_unit_null(engine, src, sets, w)
                    for s, v in unit.items():
                        units.append({"scheme": scheme, "dataset": ds, "program_id": row.program_id,
                                      "sample_id": s, "residual_moran_i": v, "unit_null_mean": unull[s],
                                      "unit_excess": v - unull[s]})
                    rec.update({
                        "scheme": scheme, "dataset": ds, "program_id": row.program_id,
                        "program_name": row.program_name, "match_lineage": lineage,
                        "n_measured": len(measured),
                        "retained_l1_weight": float(measured["original_l1_weight"].sum()),
                        "n_genes_expression_relaxed": int(np.sum(np.asarray(er) > 0)),
                        "n_genes_lineage_relaxed": int(np.sum(np.asarray(lr) > 0)),
                    })
                    rows.append(rec)
                    if scheme == "hepatocyte_primary_repro":
                        ref = r2[(r2.dataset == ds) & (r2.program_id == row.program_id)].iloc[0]
                        for col in ("matched_set_sha256", "residual_moran_i", "residual_null_mean", "residual_pvalue",
                                    "n_measured", "retained_l1_weight", "equal_residual_null_mean",
                                    "leave_top_residual_null_mean"):
                            a, b = rec[col], ref[col]
                            ok = (a == b) if isinstance(b, str) else bool(np.isclose(float(a), float(b), rtol=0, atol=1e-9))
                            checks.append({"dataset": ds, "program_id": row.program_id, "field": col,
                                           "this_run": a, "r2_file": b, "match": ok})
            del src
        res = pd.DataFrame(rows)
        for scheme in res["scheme"].unique():
            for ds in datasets:
                m = (res.scheme == scheme) & (res.dataset == ds)
                res.loc[m, "residual_padj_two_program"] = engine.bh(res.loc[m, "residual_pvalue"]).values
        res["robust"] = (res["residual_padj_two_program"] < 0.05) & res["sensitivity_sign_agree"]
        res.to_csv(args.out / "focal_lineage_results.tsv", sep="\t", index=False)
        pd.DataFrame(units).to_csv(args.out / "focal_unit_results.tsv", sep="\t", index=False)
        ck = pd.DataFrame(checks)
        ck.to_csv(args.out / "reproduction_check.tsv", sep="\t", index=False)
        log(f"reproduction checks passed {int(ck['match'].sum())}/{len(ck)}")
        return 0 if ck["match"].all() else 2

    # family mode
    registry, membership = full_registry(paths)
    union = set(membership["gene_symbol"].dropna().astype(str))
    assign = balanced_slices(registry, membership, args.n_slices)
    registry.to_csv(args.out / "family_registry.tsv", sep="\t", index=False)
    for ds in datasets:
        src = Source(engine, ds, paths.datasets[ds], designs[ds]["design"], designs[ds]["obs"].index, biotype)
        rng = np.random.default_rng(42 + (0 if ds == "GSE192741" else 1000))
        rows, nulls, set_hashes = [], [], []
        t0 = time.time()
        plan = []
        for row in registry.itertuples(index=False):
            measured = src.measured(membership, row.program_id)
            retained = float(measured["original_l1_weight"].sum())
            testable = len(measured) >= 8 and retained >= 0.20
            base = {"dataset": ds, "program_id": row.program_id, "display_order": row.display_order,
                    "cell_type": row.cell_type, "module": row.module,
                    "legacy_program_id": f"{row.cell_type}::{row.module}", "program_name": row.program_name,
                    "robust_display": row.robust_display, "n_source_genes": row.n_source_genes,
                    "n_confirmed_symbols": int((membership.program_id == row.program_id).sum()),
                    "n_measured": len(measured), "retained_l1_weight": retained, "testable": testable,
                    "slice": assign[row.program_id], "n_null": args.n_null}
            if not testable:
                if assign[row.program_id] == args.slice:
                    rows.append(base)
                continue
            # engine semantics: build matched sets for every testable program in order
            sets, er, lr = engine.build_matched_sets(
                measured, src.lineage_pools[row.cell_type], union, args.n_null, rng
            )
            set_hashes.append({"dataset": ds, "program_id": row.program_id,
                               "matched_set_sha256": hashlib.sha256(sets.tobytes()).hexdigest()})
            base["n_genes_expression_relaxed"] = int(np.sum(np.asarray(er) > 0))
            base["n_genes_lineage_relaxed"] = int(np.sum(np.asarray(lr) > 0))
            if assign[row.program_id] == args.slice:
                plan.append((row, measured, sets, base))
        log(f"{ds}: matched sets built for all testable programs in {time.time()-t0:.0f}s; slice has {len(plan)}")
        pd.DataFrame(set_hashes).to_csv(args.out / f"matched_set_hashes_{ds}.tsv", sep="\t", index=False)
        for row, measured, sets, base in plan:
            t1 = time.time()
            rec, null, _ = run_program(engine, src, row, measured, sets, args.n_sens)
            rows.append({**base, **rec})
            nulls.append(pd.DataFrame({"dataset": ds, "program_id": row.program_id,
                                       "draw_id": np.arange(1, len(null) + 1), "residual_moran_i": null}))
            log(f"{ds} {base['legacy_program_id']} {row.program_name}: n={len(measured)} excess={rec['excess']:.4f} "
                f"p={rec['residual_pvalue']:.4g} ({time.time()-t1:.0f}s)")
            pd.DataFrame(rows).to_csv(args.out / f"family_results_{ds}.partial.tsv", sep="\t", index=False)
        pd.DataFrame(rows).to_csv(args.out / f"family_results_{ds}.tsv", sep="\t", index=False)
        if nulls:
            pd.concat(nulls).to_csv(args.out / f"family_null_draws_{ds}.tsv.gz", sep="\t", index=False,
                                    compression={"method": "gzip", "mtime": 0})
        del src
    log("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
