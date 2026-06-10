#!/usr/bin/env python3
# ==========================================================================
# 269b: Aggregate Bayesian-network figure inputs
#
# Producer glue step. The Bayesian multiplex network pipeline (260-269)
# writes its per-layer posterior edge tables as Parquet files under
# results/network/posterior_edges/ (each named edges_<layer>_posterior.csv
# despite the .csv extension), and per-gene community membership under
# results/network/communities/community_assignments.csv. The supplementary
# figure scripts (figS_network_stats.R Panel B) instead read two *aggregated*
# CSVs that no pipeline step currently emits:
#
#   results/network/posterior_edges.csv               (gene_a, gene_b, layer, posterior)
#   results/network/communities/community_labels.csv  (gene, community, label)
#
# This script materialises those two CSVs from the genuine producer outputs.
# It fabricates nothing: posteriors come verbatim from 262_compute_posteriors,
# community membership from the Leiden assignment, and the per-community label
# is the top-OR Hallmark term per macro community from 268's enrichment table.
#
# It also rewrites bayesian_diagnostics.csv (a Parquet mislabelled .csv that
# crashes data.table::fread with "embedded nul") into a genuine CSV so that
# the figure's Panel F degrades gracefully instead of halting the whole script.
# The original Parquet is preserved as bayesian_diagnostics.parquet.
#
# Env: rnaseq (pandas + pyarrow). CPU only, light.
# ==========================================================================
import os
import sys
import shutil
from pathlib import Path

import pandas as pd

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
NET_DIR = BASE / "RNA-seq/results/network"
POST_DIR = NET_DIR / "posterior_edges"
COMM_DIR = NET_DIR / "communities"

LAYERS = ["ppi", "coexpr", "regulon", "lr", "genetic", "pathway", "cerna"]


def log(msg):
    print(f"[269b] {msg}", flush=True)


def aggregate_posterior_edges():
    """Concat per-layer posterior parquets into one long-form CSV."""
    out_path = NET_DIR / "posterior_edges.csv"
    frames = []
    for layer in LAYERS:
        p = POST_DIR / f"edges_{layer}_posterior.csv"
        if not p.exists():
            log(f"  WARNING: missing per-layer posterior file: {p}")
            continue
        # These files are Parquet despite the .csv name.
        try:
            df = pd.read_parquet(p)
        except Exception as e:
            log(f"  ERROR reading {p} as parquet: {e}")
            continue
        keep = [c for c in ["gene_a", "gene_b", "posterior"] if c in df.columns]
        if not {"gene_a", "gene_b", "posterior"}.issubset(df.columns):
            log(f"  WARNING: {layer} missing required cols; has {list(df.columns)}")
            continue
        sub = df[keep].copy()
        sub["layer"] = layer
        frames.append(sub[["gene_a", "gene_b", "layer", "posterior"]])
        log(f"  {layer}: {len(sub):,} edges")
    if not frames:
        log("  FATAL: no per-layer posterior frames assembled")
        sys.exit(1)
    allp = pd.concat(frames, ignore_index=True)
    allp.to_csv(out_path, index=False)
    log(f"  wrote {out_path}  ({len(allp):,} rows x {allp.shape[1]} cols)")
    log(f"  layers present: {sorted(allp['layer'].unique())}")
    return out_path


def build_community_labels():
    """Per-gene gene -> community(-> label) table for the figure.

    community = macro_id from community_assignments.csv.
    label     = top odds-ratio Hallmark term per macro community from
                community_enrichment.csv (falls back to 'Community <id>').
    """
    out_path = COMM_DIR / "community_labels.csv"
    assign_path = COMM_DIR / "community_assignments.csv"
    enrich_path = COMM_DIR / "community_enrichment.csv"

    assign = pd.read_csv(assign_path)
    if "macro_id" not in assign.columns:
        log(f"  FATAL: community_assignments.csv lacks macro_id; "
            f"has {list(assign.columns)}")
        sys.exit(1)
    out = assign[["gene", "macro_id"]].rename(columns={"macro_id": "community"})

    # Derive a human-readable label per macro community from the enrichment
    # table: top Hallmark term by odds_ratio among significant macro rows.
    label_map = {}
    if enrich_path.exists():
        enr = pd.read_csv(enrich_path)
        macro = enr[(enr.get("level") == "macro")]
        if "category" in macro.columns:
            hm = macro[macro["category"] == "hallmark"].copy()
        else:
            hm = macro.copy()
        if "odds_ratio" in hm.columns and len(hm):
            hm = hm.sort_values("odds_ratio", ascending=False)
            for cid, grp in hm.groupby("community_id"):
                term = str(grp.iloc[0]["term"])
                term = term.replace("HALLMARK_", "").replace("_", " ").title()
                label_map[int(cid)] = term
    out["label"] = out["community"].map(
        lambda c: label_map.get(int(c), f"Community {int(c)}"))

    out.to_csv(out_path, index=False)
    log(f"  wrote {out_path}  ({len(out):,} genes, "
        f"{out['community'].nunique()} communities, "
        f"{len(label_map)} labelled macro communities)")
    return out_path


def fix_bayesian_diagnostics_csv():
    """Rewrite the parquet-as-.csv diagnostics file into a genuine CSV.

    figS_network_stats Panel F freads bayesian_diagnostics.csv; the file is
    actually Parquet (magic PAR1) and crashes fread. Preserve the original as
    .parquet and emit a real CSV so Panel F degrades gracefully.
    """
    diag = NET_DIR / "bayesian_diagnostics.csv"
    if not diag.exists():
        log("  bayesian_diagnostics.csv absent; skipping")
        return None
    with open(diag, "rb") as fh:
        magic = fh.read(4)
    if magic != b"PAR1":
        log("  bayesian_diagnostics.csv already a non-parquet file; leaving as-is")
        return diag
    bak = NET_DIR / "bayesian_diagnostics.parquet"
    if not bak.exists():
        shutil.copy2(diag, bak)
        log(f"  preserved original parquet -> {bak}")
    df = pd.read_parquet(bak)
    df.to_csv(diag, index=False)
    log(f"  rewrote {diag} as CSV ({df.shape[0]} rows x {df.shape[1]} cols; "
        f"cols={list(df.columns)})")
    return diag


def main():
    log(f"NET_DIR = {NET_DIR}")
    aggregate_posterior_edges()
    build_community_labels()
    fix_bayesian_diagnostics_csv()
    log("done")


if __name__ == "__main__":
    main()
