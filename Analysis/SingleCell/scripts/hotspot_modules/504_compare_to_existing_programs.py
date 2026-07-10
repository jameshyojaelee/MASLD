"""Compute Jaccard + cosine + per-cell Pearson overlap of each Hotspot module
against 5 reference panels.

Panels:
  - cnmf_k16        (gene_spectra_score weights + per-cell usages for Pearson)
  - bulk_nmf_k6     (R NMF basis W; pre-extracted to a TSV by _extract_bulk_nmf_W.R)
  - hallmark        (50 gene sets)
  - scenic_hep      (TF -> targets)
  - curated_liver   (Govaere/Feng/Tzouanas)

Match rule (per plan section 6): a module is MATCHED if ANY of
  jaccard_top50 >= 0.20
  cosine        >= 0.35
  per_cell_pearson >= 0.40   (cNMF only — only panel with per-cell scores)
holds against ANY program in ANY panel. NOVEL otherwise.

Output:
  novelty_matches_full.tsv   one row per (cell_type, module, panel, program)
  novelty_matches.tsv        one row per (cell_type, module) with best match + is_novel
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from hotspot_io import (
    RESULTS, REF_DIR, RUN_ORDER, CNMF_GENE_SPECTRA_K16, CNMF_USAGES_K16,
)

TOPN = 50
JACCARD_THRESH = 0.20
COSINE_THRESH = 0.35
PEARSON_THRESH = 0.40


def jaccard(a: set, b: set) -> float:
    return len(a & b) / max(len(a | b), 1)


def cosine(u: pd.Series, v: pd.Series) -> float:
    """Cosine on the UNION of gene indices with zero-padding (audit fix
    2026-05-17). The earlier intersection-based implementation forced
    cosine=1.0 whenever a constant-1 membership reference vector
    (Hallmark / SCENIC+ / curated_liver panels) collapsed to a 1-gene
    intersection — a structural artifact rather than a similarity signal.
    Computing over the union with explicit zero-padding restores meaningful
    weighting and forces vectors with no real overlap toward 0.

    Additionally, gate sparse-vector comparisons with a min-overlap requirement
    (<5 shared genes after taking nonzero entries → return 0). Without this,
    weight vectors that happen to share 1-2 high-weight genes can still
    produce inflated scores."""
    idx = u.index.union(v.index)
    u_ = u.reindex(idx, fill_value=0.0).values
    v_ = v.reindex(idx, fill_value=0.0).values
    n_shared_nonzero = int(((u_ != 0) & (v_ != 0)).sum())
    if n_shared_nonzero < 5:
        return 0.0
    nu, nv = np.linalg.norm(u_), np.linalg.norm(v_)
    if nu == 0 or nv == 0:
        return 0.0
    return float(np.dot(u_, v_) / (nu * nv))


def load_hotspot_modules() -> dict[tuple[str, int], pd.Series]:
    """{(cell_type, module_int): Series gene->weight}"""
    out = {}
    for ct in RUN_ORDER:
        mg_path = RESULTS / ct / "module_genes.tsv"
        if not mg_path.exists():
            print(f"[WARN] {mg_path} missing; skipping {ct}")
            continue
        mg = pd.read_csv(mg_path, sep="\t")
        for mod, df in mg.groupby("module"):
            out[(ct, int(mod))] = df.set_index("gene")["weight"]
    return out


def load_hotspot_cell_scores() -> pd.DataFrame:
    """Long-form cell_id, cell_type, module, score across all 7 CTs.

    Used for the per-cell Pearson criterion only. Returns empty DF if 502 has
    not yet been run (cell_scores_all.parquet absent) or per-CT files missing.
    """
    consolidated = RESULTS / "cell_scores_all.parquet"
    if consolidated.exists():
        try:
            df = pd.read_parquet(consolidated)
            # 502 namespaces module as "<ct>__<int>"; canonicalize to (ct, int)
            if df["module"].dtype == object and df["module"].str.contains("__").all():
                split = df["module"].str.split("__", n=1, expand=True)
                df = df.assign(cell_type=split[0], module=split[1].astype(int))
            return df[["cell_id", "cell_type", "module", "score"]]
        except Exception as e:
            # Present-but-unreadable (e.g. corrupt/incompatible parquet) — fall
            # through to the per-CT files, which after the clean re-runs are the
            # current scores anyway.
            print(f"[WARN] could not read consolidated {consolidated} ({e}); "
                  f"falling back to per-CT cell_scores.parquet")
    # Fallback: read per-CT files directly
    frames = []
    for ct in RUN_ORDER:
        p = RESULTS / ct / "cell_scores.parquet"
        if not p.exists():
            continue
        try:
            df = pd.read_parquet(p)
        except Exception as e:
            print(f"[WARN] could not read {p} ({e}); skipping {ct} for Pearson criterion")
            continue
        df["cell_type"] = ct
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["cell_id", "cell_type", "module", "score"])
    return pd.concat(frames, ignore_index=True)


def load_cnmf_usages() -> pd.DataFrame | None:
    """cell_id, program (1..16) -> usage. None if file missing."""
    if not CNMF_USAGES_K16.exists():
        print(f"[WARN] {CNMF_USAGES_K16} missing; Pearson criterion skipped")
        return None
    # File is cell-rows x 16-program-cols (header is "index\t1\t2\t...\t16")
    usg = pd.read_csv(CNMF_USAGES_K16, sep="\t", index_col=0)
    usg.index.name = "cell_id"
    print(f"   cNMF usages: {usg.shape[0]:,} cells x {usg.shape[1]} programs")
    return usg


def load_reference_panels() -> dict[str, dict[str, pd.Series]]:
    """Each panel maps program_name -> Series gene -> weight (1.0 for set membership)."""
    panels: dict[str, dict[str, pd.Series]] = {}

    # cNMF k=16: file is laid out as program-rows × gene-cols (16 x ~37K).
    # Transpose so each program becomes a Series indexed by gene.
    cnmf = pd.read_csv(CNMF_GENE_SPECTRA_K16, sep="\t", index_col=0).T
    panels["cnmf_k16"] = {
        f"P{int(col)}": cnmf[col].dropna() for col in cnmf.columns
    }

    # bulk_nmf_k6 — extract via R helper into a TSV before running this script
    bulk_path = REF_DIR / "bulk_nmf_k6_W.tsv"
    if bulk_path.exists():
        bulk = pd.read_csv(bulk_path, sep="\t", index_col=0)
        panels["bulk_nmf_k6"] = {col: bulk[col].dropna() for col in bulk.columns}
    else:
        print(f"[WARN] {bulk_path} missing; skipping bulk_nmf_k6 panel")

    # Membership panels (1.0 weights)
    for nm in ["hallmark_v2025_1", "scenic_hep_regulons", "liver_curated"]:
        path = REF_DIR / f"{nm}.tsv"
        if not path.exists() or path.stat().st_size < 64:
            print(f"[WARN] {path} missing or empty; skipping")
            continue
        df = pd.read_csv(path, sep="\t")
        panel_name = df["panel"].iloc[0]
        if df.empty:
            print(f"[WARN] {path} has only a header; skipping")
            continue
        panels[panel_name] = {
            prog: pd.Series(1.0, index=g["gene"].unique())
            for prog, g in df.groupby("program")
        }
    return panels


def compare_one(hot_w: pd.Series, ref_w: pd.Series) -> tuple[float, float]:
    top_hot = set(hot_w.sort_values(ascending=False).head(TOPN).index)
    top_ref = set(ref_w.sort_values(ascending=False).head(TOPN).index)
    return jaccard(top_hot, top_ref), cosine(hot_w, ref_w)


def compute_cnmf_pearson(
    cell_scores: pd.DataFrame,
    cnmf_usages: pd.DataFrame,
) -> dict[tuple[str, int, str], float]:
    """For each Hotspot module x cNMF program pair, Pearson on shared cells.

    Returns {(cell_type, hotspot_module, cnmf_program_name): rho}.
    """
    out: dict[tuple[str, int, str], float] = {}
    if cnmf_usages is None or cell_scores.empty:
        return out
    # cNMF usages indexed by cell_id; columns are program ints
    usg = cnmf_usages.copy()
    usg.columns = [f"P{int(c)}" for c in usg.columns]

    for (ct, mod), df in cell_scores.groupby(["cell_type", "module"]):
        # Pivot Hotspot scores to a Series indexed by cell_id
        hs_scores = df.set_index("cell_id")["score"]
        common = hs_scores.index.intersection(usg.index)
        if len(common) < 50:
            for prog in usg.columns:
                out[(ct, int(mod), prog)] = float("nan")
            continue
        hs_v = hs_scores.loc[common].values
        # Vectorized Pearson against all 16 cNMF programs
        u = usg.loc[common].values
        # Center
        hs_c = hs_v - hs_v.mean()
        u_c = u - u.mean(axis=0, keepdims=True)
        num = u_c.T @ hs_c
        denom = np.sqrt((u_c ** 2).sum(axis=0) * (hs_c ** 2).sum())
        denom[denom == 0] = np.nan
        rhos = num / denom
        for prog, rho in zip(usg.columns, rhos):
            out[(ct, int(mod), prog)] = float(rho)
    return out


def main() -> None:
    print("[1/4] Loading Hotspot modules")
    hot = load_hotspot_modules()
    print(f"   {len(hot)} (cell_type, module) pairs")

    print("[2/4] Loading reference panels")
    panels = load_reference_panels()
    print(f"   {len(panels)} panels: {list(panels)}")

    print("[3/4] Loading per-cell scores + cNMF usages for Pearson criterion")
    cell_scores = load_hotspot_cell_scores()
    cnmf_usages = load_cnmf_usages()
    pearson_lookup = compute_cnmf_pearson(cell_scores, cnmf_usages)
    print(f"   Pearson pairs computed: {len(pearson_lookup):,}")

    print("[4/4] Computing overlaps")
    rows = []
    for (ct, mod), hw in hot.items():
        for panel_name, programs in panels.items():
            for prog, rw in programs.items():
                jac, cos = compare_one(hw, rw)
                rho = float("nan")
                if panel_name == "cnmf_k16":
                    rho = pearson_lookup.get((ct, mod, prog), float("nan"))
                rows.append({
                    "cell_type": ct, "module": mod,
                    "panel": panel_name, "program": prog,
                    "jaccard_top50": jac, "cosine": cos,
                    "per_cell_pearson": rho,
                })

    matches = pd.DataFrame(rows)
    matches.to_csv(RESULTS / "novelty_matches_full.tsv", sep="\t", index=False)

    matches["matched"] = (
        (matches["jaccard_top50"] >= JACCARD_THRESH)
        | (matches["cosine"] >= COSINE_THRESH)
        | (matches["per_cell_pearson"] >= PEARSON_THRESH)
    )

    # Pick best-match row per (cell_type, module): max jaccard wins
    summary = (
        matches.sort_values(["cell_type", "module", "jaccard_top50"], ascending=[True, True, False])
        .groupby(["cell_type", "module"], as_index=False)
        .first()
        .rename(columns={
            "panel": "best_match_panel",
            "program": "best_match_program",
            "jaccard_top50": "best_match_jaccard",
            "cosine": "best_match_cosine",
            "per_cell_pearson": "best_match_pearson",
        })
    )
    is_novel = (
        ~matches.groupby(["cell_type", "module"])["matched"].any()
    ).reset_index().rename(columns={"matched": "is_novel"})
    summary = summary.merge(is_novel, on=["cell_type", "module"], how="left")
    summary.to_csv(RESULTS / "novelty_matches.tsv", sep="\t", index=False)
    print(f"   wrote {len(summary)} module summaries; {summary['is_novel'].sum()} novel")


if __name__ == "__main__":
    main()
