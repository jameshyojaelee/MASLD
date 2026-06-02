"""Olink (secreted proteins) + Visium (spatial zonation) validation of Hotspot modules.

Olink: For each module, count genes encoding secreted proteins (from UniProt
       subcellular annotation cache; if missing, fall back to a hard-coded list).
       Test plasma protein level vs disease_stage via Spearman per protein,
       summarize concordance.

Visium: For each module, score spots and test zonation enrichment.
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from hotspot_io import PROJECT_ROOT, RESULTS, RUN_ORDER

OLINK_DIR = PROJECT_ROOT / "data/olink"
VISIUM_DIR = PROJECT_ROOT / "Analysis/Spatial/results"


def run_olink() -> pd.DataFrame:
    """Spearman of plasma protein level vs disease_stage per module-member protein."""
    olink_long = OLINK_DIR / "olink_long.tsv"        # adjust if your loader uses a different name
    if not olink_long.exists():
        print(f"[WARN] {olink_long} not found; skipping Olink arm")
        return pd.DataFrame()
    olink = pd.read_csv(olink_long, sep="\t")        # expected cols: sample, protein, npx, disease_stage_ordinal
    rows = []
    for ct in RUN_ORDER:
        mg_p = RESULTS / ct / "module_genes.tsv"
        if not mg_p.exists():
            print(f"[WARN] skipping {ct} (missing {mg_p})")
            continue
        mg = pd.read_csv(mg_p, sep="\t")
        for mod, g in mg.groupby("module"):
            top = g.sort_values("weight", ascending=False).head(200)["gene"]
            sub = olink[olink["protein"].isin(top)]
            if sub.empty:
                continue
            per_prot = (
                sub.groupby("protein")
                .apply(lambda d: spearmanr(d["npx"], d["disease_stage_ordinal"]).statistic)
                .reset_index(name="rho")
            )
            rows.append({
                "cell_type": ct, "module": int(mod),
                "n_proteins_tested": len(per_prot),
                "olink_pct_positive_rho": float((per_prot["rho"] > 0).mean()),
                "olink_mean_rho": float(per_prot["rho"].mean()),
            })
    return pd.DataFrame(rows)


def run_visium() -> pd.DataFrame:
    """Score hepatocyte/stellate modules per Visium spot; test zonation enrichment."""
    visium_spots = VISIUM_DIR / "gse192741_spot_metadata.tsv"   # adjust to your spatial output
    visium_expr = VISIUM_DIR / "gse192741_spot_expression.parquet"
    if not visium_spots.exists() or not visium_expr.exists():
        print("[WARN] Visium inputs not found; skipping spatial arm")
        return pd.DataFrame()
    spots = pd.read_csv(visium_spots, sep="\t").set_index("spot_id")
    expr = pd.read_parquet(visium_expr).set_index("spot_id")    # cells = spots, cols = genes

    rows = []
    for ct in ["hepatocytes", "fibroblasts", "macrophages"]:
        mg_p = RESULTS / ct / "module_genes.tsv"
        if not mg_p.exists():
            print(f"[WARN] skipping {ct} (missing {mg_p})")
            continue
        mg = pd.read_csv(mg_p, sep="\t")
        for mod, g in mg.groupby("module"):
            top = g.sort_values("weight", ascending=False).head(100)
            common = expr.columns.intersection(top["gene"])
            if len(common) < 5:
                continue
            w = top.set_index("gene").loc[common, "weight"].values
            score = (expr[common].values * w).sum(axis=1) / w.sum()
            score_s = pd.Series(score, index=expr.index)
            merged = spots.join(score_s.rename("module_score"), how="inner")

            zone = "spatial_domain"     # adjust if your spatial labels use another column
            if zone in merged.columns:
                periportal = merged[merged[zone].str.contains("portal", case=False, na=False)]["module_score"]
                pericentral = merged[merged[zone].str.contains("central", case=False, na=False)]["module_score"]
                rho = (periportal.mean() - pericentral.mean()) / merged["module_score"].std()
                rows.append({"cell_type": ct, "module": int(mod),
                             "visium_periportal_minus_pericentral_z": float(rho),
                             "n_spots": int(len(merged))})
    return pd.DataFrame(rows)


def main() -> None:
    olink_df = run_olink()
    visium_df = run_visium()
    # Ensure both have merge keys even when empty (graceful-skip arms).
    base_cols = ["cell_type", "module"]
    if olink_df.empty:
        olink_df = pd.DataFrame(columns=base_cols)
    if visium_df.empty:
        visium_df = pd.DataFrame(columns=base_cols)
    if olink_df.empty and visium_df.empty:
        out = pd.DataFrame(columns=base_cols)
    else:
        out = olink_df.merge(visium_df, on=base_cols, how="outer")
    out.to_csv(RESULTS / "olink_visium_validation.tsv", sep="\t", index=False)
    print(f"Wrote {len(out)} validation rows")


if __name__ == "__main__":
    main()
