#!/usr/bin/env python
"""12: spatial coherence of each module, against a matched-gene null.

WHAT THIS COLUMN IS AND IS NOT. It does not ask whether a module tracks disease
stage; four donors cannot answer that. It asks whether a module's score is
spatially organised in tissue beyond what genes of the same expression level
and detection rate produce on the same sections. Matching is on expression and
detection only. This is a bounded diagnostic and cannot rescue a disease
association or a regulatory claim.

v2 corrections. The neighbour graph is built within each physical section
(`sample_id`), never across sections of one donor, whose array coordinates
restart; section statistics are collapsed to donors afterwards. Matched draws
are without replacement inside a draw, and a candidate bin that is exhausted is
widened to the nearest bins by rank distance rather than falling back to the
whole pool.
"""
import argparse
import json
import pathlib
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CONTRACT = json.loads((ROOT / "scripts/analysis/cross_assay_modules/00_contract.json").read_text())
N_NULL = 1000
K_NEIGHBOURS = 6
MIN_SPOTS_PER_SECTION = 50


def spatial_graph(coords: np.ndarray) -> sparse.csr_matrix:
    """Row-standardised k-nearest-neighbour weights on one section's coordinates."""
    k = min(K_NEIGHBOURS, len(coords) - 1)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    rows = np.repeat(np.arange(len(coords)), k)
    cols = idx[:, 1:].reshape(-1)
    w = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(coords),) * 2)
    rs = np.asarray(w.sum(1)).ravel()
    rs[rs == 0] = 1
    return sparse.diags(1 / rs) @ w


def moran(scores: np.ndarray, w: sparse.csr_matrix) -> np.ndarray:
    x = scores - scores.mean(0, keepdims=True)
    denom = (x * x).sum(0)
    denom[denom == 0] = np.nan
    return np.einsum("ij,ij->j", x, w @ x) / denom


def residualise(y: np.ndarray, design: np.ndarray) -> np.ndarray:
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ coef


def draw_matched_set(members, pool_bins, bin_members, excluded, rng):
    """One matched set: per member, a distinct gene from the same (mean, detection)
    bin, widening by rank distance if the bin is exhausted. Never repeats a gene
    inside a draw and never uses a module member."""
    picked = []
    taken = set(excluded)
    for g in members:
        mb, db = pool_bins[g]
        for radius in range(0, 20):
            cand = [c for (m2, d2), cs in bin_members.items()
                    if abs(m2 - mb) + abs(d2 - db) == radius for c in cs if c not in taken]
            if cand:
                choice = int(rng.choice(cand))
                picked.append(choice)
                taken.add(choice)
                break
        else:
            raise RuntimeError("matched draw exhausted the pool")
    return np.array(picked)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", required=True, type=pathlib.Path)
    args = ap.parse_args()
    out = args.out_root / "assays"
    out.mkdir(parents=True, exist_ok=True)
    target = out / "spatial_results.tsv"
    if target.exists():
        print("[12] spatial results exist; refusing to overwrite", flush=True)
        return 0

    memb = pd.read_csv(args.out_root / "modules/module_membership.tsv", sep="\t")
    modules = memb.groupby("module_id")["gene_symbol"].apply(list).to_dict()
    module_ids = sorted(modules)

    rows = []
    unit_census = []
    for tag, key, unit in (("visium_gse192741", "visium_guilliams", "biological_donor"),
                           ("visium_vu", "visium_vu", "array")):
        a = ad.read_h5ad(ROOT / CONTRACT["inputs"][key])
        X = a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)
        var = pd.Index([str(v) for v in a.var_names])
        sections = a.obs["sample_id"].astype(str).values
        donors = (a.obs["individual"].astype(str).values if unit == "biological_donor"
                  else a.obs["sample_id"].astype(str).values)
        abundance = np.asarray(a.obsm["q05_cell_abundance_w_sf"])
        libsize = np.asarray(a.obs["log1p_total_counts"]).reshape(-1, 1)
        ngenes = np.asarray(a.obs["log1p_n_genes_by_counts"]).reshape(-1, 1)

        mean_expr = X.mean(0)
        detection = (X > 0).mean(0)
        mean_bin = pd.qcut(pd.Series(mean_expr).rank(method="first"), 10, labels=False, duplicates="drop").values
        det_bin = pd.qcut(pd.Series(detection).rank(method="first"), 10, labels=False, duplicates="drop").values
        gene_pos = {g: i for i, g in enumerate(var)}
        pool_bins = {g: (int(mean_bin[i]), int(det_bin[i])) for g, i in gene_pos.items()}
        bin_members: dict = {}
        for i in range(len(var)):
            bin_members.setdefault((int(mean_bin[i]), int(det_bin[i])), []).append(i)
        rng = np.random.default_rng(CONTRACT["seed"])

        # One graph and design per SECTION; donor of each section recorded.
        blocks = []
        for sec in pd.unique(sections):
            sel = np.where(sections == sec)[0]
            if len(sel) < MIN_SPOTS_PER_SECTION:
                continue
            donor = pd.unique(donors[sel])
            assert len(donor) == 1, f"section {sec} maps to several donors"
            coords = np.asarray(a.obsm["spatial"])[sel]
            design = np.hstack([np.ones((len(sel), 1)), abundance[sel], libsize[sel], ngenes[sel]])
            blocks.append((sec, str(donor[0]), sel, spatial_graph(coords), design))
        block_donor = np.array([b[1] for b in blocks])
        donor_ids = pd.unique(block_donor)
        unit_census.append(dict(assay=tag, n_sections=len(blocks), n_donors=len(donor_ids),
                                sections_per_donor=";".join(f"{d}:{(block_donor == d).sum()}" for d in donor_ids)))
        print(f"[12] {tag}: {len(blocks)} sections, {len(donor_ids)} donors, {X.shape[1]} genes", flush=True)

        Z = (X - X.mean(0)) / np.where(X.std(0) == 0, np.nan, X.std(0))

        def collapsed_moran(cols_per_set):
            """Residual Moran per section, mean within donor, then mean over donors."""
            per_sec = np.full((len(blocks), len(cols_per_set)), np.nan)
            for bi, (_sec, _d, sel, w, design) in enumerate(blocks):
                scores = np.column_stack([
                    np.nanmean(Z[np.ix_(sel, c)], axis=1) if len(c) else np.full(len(sel), np.nan)
                    for c in cols_per_set])
                scores = np.nan_to_num(scores, nan=0.0)
                per_sec[bi] = moran(residualise(scores, design), w)
            per_donor = np.vstack([np.nanmean(per_sec[block_donor == d], axis=0) for d in donor_ids])
            return np.nanmean(per_donor, axis=0)

        for mid in module_ids:
            members = [g for g in modules[mid] if g in gene_pos]
            n_meas = len(members)
            frac = n_meas / len(modules[mid])
            testable = (n_meas >= CONTRACT["scoring"]["testability"]["min_members_measured"]
                        and frac >= CONTRACT["scoring"]["testability"]["min_fraction_measured"])
            base = dict(assay=tag, module_id=mid, unit=unit, n_members=len(modules[mid]),
                        n_measured=n_meas, fraction_measured=frac,
                        n_sections=len(blocks), n_units=len(donor_ids))
            if not testable:
                rows.append(dict(base, testable=False, observed_residual_moran=np.nan,
                                 null_mean=np.nan, null_sd=np.nan, z_vs_null=np.nan,
                                 empirical_p=np.nan, n_null=0))
                continue
            cols = np.array([gene_pos[g] for g in members])
            obs = float(collapsed_moran([cols])[0])
            excluded = set(cols.tolist())
            null_sets = [draw_matched_set(members, pool_bins, bin_members, excluded, rng)
                         for _ in range(N_NULL)]
            assert all(len(set(s.tolist())) == len(s) for s in null_sets), "duplicate gene in a draw"
            nv = np.asarray(collapsed_moran(null_sets), dtype=float)
            nv = nv[np.isfinite(nv)]
            mu, sd = float(np.mean(nv)), float(np.std(nv, ddof=1))
            rows.append(dict(base, testable=True, observed_residual_moran=obs, null_mean=mu, null_sd=sd,
                             z_vs_null=(obs - mu) / sd if sd > 0 else np.nan,
                             empirical_p=(1.0 + float(np.sum(nv >= obs))) / (len(nv) + 1.0),
                             n_null=len(nv)))
        del a, X, Z

    df = pd.DataFrame(rows)
    banned = set(CONTRACT["prohibited_columns"]) & set(df.columns)
    if banned:
        raise SystemExit(f"prohibited columns: {banned}")
    df.to_csv(target, sep="\t", index=False)
    pd.DataFrame(unit_census).to_csv(out / "spatial_unit_census.tsv", sep="\t", index=False)
    (out / "spatial_summary.json").write_text(json.dumps({
        "n_null_draws": N_NULL, "k_neighbours": K_NEIGHBOURS,
        "assays": {t: int(df[(df.assay == t) & df.testable].shape[0]) for t in df.assay.unique()},
        "unit_census": unit_census,
        "statistic": "residual Moran's I per section, mean within donor, mean over donors, vs matched genes",
        "matching": "expression decile and detection decile; without replacement; nearest bins if exhausted",
        "residualised_on": "all cell-abundance factors, log library size, log detected genes",
        "graph": "k-nearest neighbours within each section only",
    }, indent=2))
    print("[12] done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
