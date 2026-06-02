#!/usr/bin/env python3
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --job-name=net_276_lls
#SBATCH --output=logs/net_276_lls_%j.out
#SBATCH --error=logs/net_276_lls_%j.err
# ===========================================================================
# Script 276: FunCoup-style LLS calibration per evidence layer
# ===========================================================================
# Purpose
# -------
# The current composite score (Script 263) combines raw confidence values
# across layers with a noisy-OR. That is a pragmatic heuristic. FunCoup 6
# (Persson et al. 2023), HumanNet v3, and STRING use a more principled
# Bayesian calibration:
#
#   For each layer m and each discretized score bin b:
#       LLS_m(b) = log( P(pair in gold+ | score in b) /
#                       P(pair in gold- | score in b) )
#
#   For a gene pair with scores (s_1, ..., s_M) across M layers:
#       integrated_LLS = sum_m LLS_m(bin(s_m))
#       P_integrated   = sigmoid(integrated_LLS)
#
# Higher integrated_LLS = more likely to be functionally linked. Missing
# layers contribute 0 (non-informative prior).
#
# Gold standard
# -------------
# We use CORUM protein complexes (or Reactome complex-like pathways if CORUM
# is unreachable). KEGG / MSigDB pathway co-membership is NOT usable as the
# pathway evidence layer was itself constructed from MSigDB Hallmark + KEGG
# + Reactome --- using KEGG as ground truth would be circular for that layer.
#
# CORUM:
#   https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt
#
# Gold negatives: random node-set gene pairs that do NOT share a CORUM
# complex, sampled to match the node-degree distribution of the positives
# (approximately, via matched random draws).
#
# Inputs
# ------
#   - RNA-seq/results/network/network_nodes.csv
#   - RNA-seq/results/network/edges_{layer}.csv      (raw edges, 7 layers)
#   - RNA-seq/results/network/composite_edges.csv    (for noisy-OR comparison)
#   - RNA-seq/results/network/id_normalization_summary.csv  (precondition)
#   - results/library/positive_control.csv            (clinical drug targets)
#   - data/corum/allComplexes.txt                     (fetched if missing)
#
# Outputs
# -------
#   results/network/lls_edges_per_layer/edges_{layer}_lls.csv
#       (gene_a, gene_b, raw_score, bin, lls)
#   results/network/lls_composite_edges.csv
#       (gene_a, gene_b, integrated_lls, p_integrated, k_contributing_layers,
#        lls_{layer} ...)
#   results/network/lls_calibration_summary.csv
#       (layer, bin, bin_low, bin_high, n_pos, n_neg, p_pos, p_neg, lls)
#   figures/misc/network_lls/lls_per_layer_calibration_curves.pdf
#   figures/misc/network_lls/lls_vs_noisyor_comparison.csv
#
# Environment
# -----------
# spatial env (numpy, pandas, scipy, matplotlib, urllib)
#
# Runtime
# -------
# ~30-60 min depending on gold-standard size and number of edges.
# ===========================================================================

from __future__ import annotations

import json
import os
import shutil
import sys
import time
import urllib.error
import urllib.request
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

NET_DIR            = BASE / "RNA-seq" / "results" / "network"
NODES_CSV          = NET_DIR / "network_nodes.csv"
PRECOND            = NET_DIR / "id_normalization_summary.csv"
COMPOSITE_CSV      = NET_DIR / "composite_edges.csv"

LLS_EDGE_DIR       = NET_DIR / "lls_edges_per_layer"
LLS_COMPOSITE_CSV  = NET_DIR / "lls_composite_edges.csv"
LLS_SUMMARY_CSV    = NET_DIR / "lls_calibration_summary.csv"

FIG_DIR            = BASE / "figures" / "misc" / "network_lls"
LLS_CURVES_PDF     = FIG_DIR / "lls_per_layer_calibration_curves.pdf"
LLS_CMP_CSV        = FIG_DIR / "lls_vs_noisyor_comparison.csv"

CORUM_DIR          = BASE / "data" / "corum"
CORUM_FILE         = CORUM_DIR / "allComplexes.txt"
CORUM_URL          = "https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt"

MSIGDB_GMT         = BASE / "downstream_analysis" / "pathway_analysis" / "data" / "msigdb.v2025.1.Hs.symbols.gmt"

POS_CONTROL_CSV    = BASE / "results" / "library" / "positive_control.csv"

LAYER_NAMES = [
    "ppi", "coexpr", "regulon", "lr", "genetic", "pathway", "cerna",
]

N_BINS           = 20            # number of quantile bins
N_POS_TARGET     = 50_000        # number of gold positive pairs
N_NEG_TARGET     = 500_000       # number of gold negative pairs
LAPLACE_PSEUDO   = 1.0           # Laplace smoothing pseudocount
RANDOM_SEED      = 42

# Govaere 25-gene prognostic panel (Govaere et al. 2020, Sci Transl Med).
# Mirror of the list used in 272b_robust_benchmark.R for reproducibility.
GOVAERE_25 = [
    "AKR1B10", "SOX4", "TTC39A", "VIL1", "CLIC6", "ABCB4", "LTBP2", "STMN2",
    "AGMO", "ANGPTL4", "CA3", "CCBE1", "CFB", "CSAD", "FABP4", "FADS1",
    "IGFBP7", "MME", "PANX2", "PLS1", "RP11-1399P15.1", "SDC1", "SRPX",
    "TMEM9", "TPBG",
]

RESMETIROM_PATHWAY = ["THRB", "PCSK9", "PPARA", "APOA1"]


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
t0_global = time.time()


def log(msg: str) -> None:
    elapsed = time.time() - t0_global
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


def canonical_pair(a: str, b: str) -> tuple[str, str]:
    """Return a sorted (min, max) tuple so (A, B) == (B, A)."""
    return (a, b) if a <= b else (b, a)


def pairs_to_frame(pairs: set[tuple[str, str]]) -> pd.DataFrame:
    """Convert a set of canonical pairs into a dataframe."""
    if not pairs:
        return pd.DataFrame(columns=["gene_a", "gene_b"])
    a, b = zip(*pairs)
    return pd.DataFrame({"gene_a": a, "gene_b": b})


def canonicalize_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure gene_a <= gene_b for every row and drop self-pairs / duplicates."""
    df = df.dropna(subset=["gene_a", "gene_b"]).copy()
    df["gene_a"] = df["gene_a"].astype(str)
    df["gene_b"] = df["gene_b"].astype(str)
    df = df[df["gene_a"] != df["gene_b"]]
    swap = df["gene_a"] > df["gene_b"]
    if swap.any():
        a = df["gene_a"].where(~swap, df["gene_b"])
        b = df["gene_b"].where(~swap, df["gene_a"])
        df["gene_a"] = a
        df["gene_b"] = b
    df = df.drop_duplicates(subset=["gene_a", "gene_b"])
    return df


# ---------------------------------------------------------------------------
# Preconditions
# ---------------------------------------------------------------------------
def check_preconditions() -> None:
    if not PRECOND.exists():
        print(f"ERROR: {PRECOND} not found — run 273_normalize_ids.py first.",
              flush=True)
        sys.exit(1)
    if not NODES_CSV.exists():
        print(f"ERROR: {NODES_CSV} not found.", flush=True)
        sys.exit(1)
    missing = [l for l in LAYER_NAMES
               if not (NET_DIR / f"edges_{l}.csv").exists()]
    if missing:
        print(f"ERROR: missing layer edge files: {missing}", flush=True)
        sys.exit(1)
    if not COMPOSITE_CSV.exists():
        log(f"WARNING: {COMPOSITE_CSV} not found — noisy-OR comparison will be skipped.")


# ---------------------------------------------------------------------------
# Gold standard: CORUM complexes (primary) with Reactome fallback
# ---------------------------------------------------------------------------
def try_download_corum() -> bool:
    """Attempt to download CORUM allComplexes.txt. Returns True on success."""
    CORUM_DIR.mkdir(parents=True, exist_ok=True)
    if CORUM_FILE.exists() and CORUM_FILE.stat().st_size > 1000:
        log(f"CORUM already cached at {CORUM_FILE}")
        return True
    log(f"Attempting CORUM download from {CORUM_URL} ...")
    try:
        req = urllib.request.Request(
            CORUM_URL,
            headers={"User-Agent": "Mozilla/5.0 (masld-network-pipeline)"},
        )
        with urllib.request.urlopen(req, timeout=60) as resp, \
                open(CORUM_FILE, "wb") as fh:
            shutil.copyfileobj(resp, fh)
        if CORUM_FILE.stat().st_size > 1000:
            log(f"CORUM downloaded ({CORUM_FILE.stat().st_size:,} bytes) -> {CORUM_FILE}")
            return True
        CORUM_FILE.unlink(missing_ok=True)
        return False
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError,
            OSError) as exc:
        log(f"CORUM download failed: {exc}")
        if CORUM_FILE.exists() and CORUM_FILE.stat().st_size <= 1000:
            CORUM_FILE.unlink(missing_ok=True)
        return False


def parse_corum(nodes: set[str]) -> dict[str, list[str]]:
    """Parse CORUM allComplexes.txt and return {complex_id: [gene_symbol,...]}.

    Filters complexes to genes in the node set. Keeps only human (organism=Human).
    """
    df = pd.read_csv(CORUM_FILE, sep="\t", dtype=str, low_memory=False)
    name_col = next((c for c in df.columns
                     if c.lower() in {"organism", "organism_name"}), None)
    if name_col is not None:
        df = df[df[name_col].astype(str).str.contains("Human", case=False, na=False)]

    subunit_col = None
    for cand in ["subunits(Gene name)", "subunits_gene_name",
                 "subunits(Gene name syn)", "Genes", "Gene name"]:
        if cand in df.columns:
            subunit_col = cand
            break
    if subunit_col is None:
        raise RuntimeError(f"Cannot find subunit column in CORUM columns: {df.columns.tolist()[:20]}")

    id_col = None
    for cand in ["ComplexID", "Complex ID", "complex_id"]:
        if cand in df.columns:
            id_col = cand
            break
    if id_col is None:
        df = df.reset_index(drop=True)
        df["ComplexID"] = df.index.astype(str)
        id_col = "ComplexID"

    complexes: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        cid = str(row[id_col])
        raw = row[subunit_col]
        if not isinstance(raw, str):
            continue
        # CORUM separates subunits with ';' (sometimes ',')
        sep = ";" if ";" in raw else ","
        members = [g.strip() for g in raw.split(sep) if g.strip()]
        members = [g for g in members if g in nodes]
        if len(members) >= 2:
            complexes[cid] = members
    return complexes


def reactome_complex_fallback(nodes: set[str]) -> dict[str, list[str]]:
    """Fallback gold standard: Reactome gene sets in MSigDB that look like
    complexes (small size, names containing COMPLEX / SUBUNIT / ASSEMBLY).

    Small pathways (<=50 members) are used to approximate CORUM-style
    complexes. This is less orthogonal than CORUM but still independent
    enough from the pathway layer's Jaccard-filtered MSigDB membership.
    """
    if not MSIGDB_GMT.exists():
        raise RuntimeError(
            f"No CORUM and no MSigDB GMT ({MSIGDB_GMT}). Cannot build gold standard."
        )
    complexes: dict[str, list[str]] = {}
    keep_keywords = ("COMPLEX", "SUBUNIT", "ASSEMBLY", "RIBOSOME", "SPLICEOSOME",
                     "PROTEASOME", "NUCLEOSOME")
    with open(MSIGDB_GMT) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            members = [g for g in parts[2:] if g in nodes]
            if not name.startswith("REACTOME_"):
                continue
            if not any(k in name for k in keep_keywords):
                continue
            if 2 <= len(members) <= 50:
                complexes[name] = members
    return complexes


def build_gold_standard(nodes: set[str], rng: np.random.Generator
                        ) -> tuple[set[tuple[str, str]], set[tuple[str, str]], str]:
    """Build gold positive and negative pair sets."""
    log("Building gold standard ...")
    source = "CORUM"
    complexes = None
    if try_download_corum():
        try:
            complexes = parse_corum(nodes)
            log(f"CORUM: {len(complexes):,} human complexes with >=2 node-set members")
        except Exception as exc:
            log(f"Failed to parse CORUM ({exc}) — falling back to Reactome.")
            complexes = None
    if not complexes:
        source = "Reactome_complex_fallback"
        complexes = reactome_complex_fallback(nodes)
        log(f"Reactome complex fallback: {len(complexes):,} gene sets")

    if not complexes:
        raise RuntimeError("No gold-standard complexes built.")

    # Positives: all within-complex pairs, canonicalised
    positives: set[tuple[str, str]] = set()
    for cid, members in complexes.items():
        members = list(dict.fromkeys(members))  # dedupe, keep order
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                positives.add(canonical_pair(members[i], members[j]))
    log(f"Gold positives (raw): {len(positives):,}")

    # Subsample positives to N_POS_TARGET if huge
    positives_list = sorted(positives)
    if len(positives_list) > N_POS_TARGET:
        idx = rng.choice(len(positives_list), size=N_POS_TARGET, replace=False)
        positives = {positives_list[i] for i in idx}
        log(f"Gold positives (subsampled): {len(positives):,}")
    # Negatives: random pairs from nodes NOT in positives
    nodes_list = np.array(sorted(nodes))
    n_nodes = len(nodes_list)
    log(f"Drawing {N_NEG_TARGET:,} gold negatives from {n_nodes:,} nodes ...")

    negatives: set[tuple[str, str]] = set()
    # Generous oversampling factor to account for collisions with positives
    # and with previously-drawn negatives.
    target = N_NEG_TARGET
    draws_done = 0
    max_rounds = 30
    for _ in range(max_rounds):
        if len(negatives) >= target:
            break
        need = max(target - len(negatives), 50_000)
        ia = rng.integers(0, n_nodes, size=need * 2)
        ib = rng.integers(0, n_nodes, size=need * 2)
        for k in range(len(ia)):
            a = nodes_list[ia[k]]
            b = nodes_list[ib[k]]
            if a == b:
                continue
            pair = canonical_pair(a, b)
            if pair in positives:
                continue
            negatives.add(pair)
            if len(negatives) >= target:
                break
        draws_done += 1
    log(f"Gold negatives: {len(negatives):,} (rounds={draws_done})")

    return positives, negatives, source


# ---------------------------------------------------------------------------
# Layer ingestion
# ---------------------------------------------------------------------------
def load_layer(layer: str) -> pd.DataFrame:
    path = NET_DIR / f"edges_{layer}.csv"
    df = pd.read_csv(path, usecols=["gene_a", "gene_b", "raw_score"])
    df = canonicalize_frame(df)
    df["raw_score"] = pd.to_numeric(df["raw_score"], errors="coerce")
    df = df.dropna(subset=["raw_score"]).reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# LLS calibration per layer (quantile binning + Laplace smoothing)
# ---------------------------------------------------------------------------
def _calibrate_layer(args):
    layer, edges, pos_set, neg_set = args
    n_edges = len(edges)
    if n_edges == 0:
        return layer, None, None

    pair_tuples = list(zip(edges["gene_a"].tolist(), edges["gene_b"].tolist()))
    in_pos = np.fromiter((p in pos_set for p in pair_tuples),
                         dtype=bool, count=n_edges)
    in_neg = np.fromiter((p in neg_set for p in pair_tuples),
                         dtype=bool, count=n_edges)

    scores = edges["raw_score"].to_numpy()

    # Quantile edges based on unique scores so bins aren't degenerate
    quantiles = np.linspace(0, 1, N_BINS + 1)
    edges_q = np.quantile(scores, quantiles)
    # Monotonise & make unique to avoid np.digitize collapses
    edges_q[0]  = -np.inf
    edges_q[-1] = np.inf
    for k in range(1, len(edges_q) - 1):
        if edges_q[k] <= edges_q[k - 1]:
            edges_q[k] = edges_q[k - 1] + 1e-12

    bin_ix = np.digitize(scores, edges_q[1:-1], right=False)
    bin_ix = np.clip(bin_ix, 0, N_BINS - 1)

    # Total pos/neg pairs observed in THIS layer's edges (denominator base).
    # FunCoup-style LLS uses P(score bin | gold+) vs P(score bin | gold-) — we
    # condition on the presence of an edge in the layer, so p_pos and p_neg
    # are computed relative to layer-level totals of annotated edges only.
    total_pos = int(in_pos.sum())
    total_neg = int(in_neg.sum())
    # If the layer lacks annotated gold pairs, fall back to equal weights
    # (LLS = 0) so the layer contributes nothing.
    if total_pos == 0 or total_neg == 0:
        lls_per_bin = np.zeros(N_BINS, dtype=float)
        summary = pd.DataFrame({
            "layer":    layer,
            "bin":      np.arange(N_BINS),
            "bin_low":  edges_q[:-1],
            "bin_high": edges_q[1:],
            "n_pos":    0,
            "n_neg":    0,
            "p_pos":    np.nan,
            "p_neg":    np.nan,
            "lls":      0.0,
            "n_edges":  np.bincount(bin_ix, minlength=N_BINS),
        })
        edges_out = edges.copy()
        edges_out["bin"] = bin_ix
        edges_out["lls"] = 0.0
        return layer, edges_out, summary

    n_pos_bin = np.bincount(bin_ix[in_pos], minlength=N_BINS)
    n_neg_bin = np.bincount(bin_ix[in_neg], minlength=N_BINS)
    n_edges_bin = np.bincount(bin_ix, minlength=N_BINS)

    # Laplace smoothing: +1 pseudocount per bin per class, normalised across bins
    p_pos = (n_pos_bin + LAPLACE_PSEUDO) / (total_pos + LAPLACE_PSEUDO * N_BINS)
    p_neg = (n_neg_bin + LAPLACE_PSEUDO) / (total_neg + LAPLACE_PSEUDO * N_BINS)
    lls_per_bin = np.log(p_pos / p_neg)

    summary = pd.DataFrame({
        "layer":    layer,
        "bin":      np.arange(N_BINS),
        "bin_low":  edges_q[:-1],
        "bin_high": edges_q[1:],
        "n_pos":    n_pos_bin,
        "n_neg":    n_neg_bin,
        "p_pos":    p_pos,
        "p_neg":    p_neg,
        "lls":      lls_per_bin,
        "n_edges":  n_edges_bin,
    })

    edges_out = edges.copy()
    edges_out["bin"] = bin_ix
    edges_out["lls"] = lls_per_bin[bin_ix]
    return layer, edges_out, summary


def calibrate_all_layers(layer_edges: dict, positives: set, negatives: set,
                         n_workers: int) -> tuple[dict, pd.DataFrame]:
    log(f"Calibrating {len(layer_edges)} layers with {n_workers} workers ...")
    args = [(layer, edges, positives, negatives)
            for layer, edges in layer_edges.items()]

    results = {}
    summaries = []
    if n_workers <= 1 or len(args) == 1:
        out = [_calibrate_layer(a) for a in args]
    else:
        with Pool(processes=min(n_workers, len(args))) as pool:
            out = pool.map(_calibrate_layer, args)

    for layer, edges_out, summary in out:
        if edges_out is None:
            log(f"  Layer {layer}: no edges — skipping.")
            continue
        results[layer] = edges_out
        summaries.append(summary)
        pos_frac = (summary["n_pos"].sum() / max(len(layer_edges[layer]), 1))
        log(f"  Layer {layer:8s}: n_edges={len(edges_out):,}  "
            f"gold+ in layer={int(summary['n_pos'].sum()):,}  "
            f"gold- in layer={int(summary['n_neg'].sum()):,}  "
            f"LLS range [{summary['lls'].min():.3f}, {summary['lls'].max():.3f}]")
    if not summaries:
        raise RuntimeError("No layers produced a calibration summary.")
    summary_df = pd.concat(summaries, ignore_index=True)
    return results, summary_df


# ---------------------------------------------------------------------------
# Write per-layer LLS edge files
# ---------------------------------------------------------------------------
def write_layer_lls(results: dict) -> None:
    LLS_EDGE_DIR.mkdir(parents=True, exist_ok=True)
    for layer, edges_out in results.items():
        out = edges_out[["gene_a", "gene_b", "raw_score", "bin", "lls"]]
        path = LLS_EDGE_DIR / f"edges_{layer}_lls.csv"
        out.to_csv(path, index=False)
        log(f"  Wrote {path} ({len(out):,} rows)")


# ---------------------------------------------------------------------------
# Composite integrated LLS
# ---------------------------------------------------------------------------
def build_composite(results: dict) -> pd.DataFrame:
    log("Integrating per-layer LLS into composite ...")
    merged = None
    for layer, edges_out in results.items():
        df = edges_out[["gene_a", "gene_b", "lls"]].rename(
            columns={"lls": f"lls_{layer}"}
        )
        merged = df if merged is None else merged.merge(
            df, on=["gene_a", "gene_b"], how="outer"
        )
        log(f"  After merging {layer}: {len(merged):,} unique pairs")

    # Fill missing layers with 0 (non-informative prior)
    lls_cols = [f"lls_{l}" for l in results.keys()]
    for c in lls_cols:
        merged[c] = merged[c].fillna(0.0)

    merged["integrated_lls"] = merged[lls_cols].sum(axis=1)
    merged["k_contributing_layers"] = (merged[lls_cols].abs() > 1e-9).sum(axis=1).astype(np.int32)
    # Sigmoid with overflow guard
    x = np.clip(merged["integrated_lls"].to_numpy(), -50.0, 50.0)
    merged["p_integrated"] = 1.0 / (1.0 + np.exp(-x))

    out_cols = ["gene_a", "gene_b", "integrated_lls", "p_integrated",
                "k_contributing_layers"] + lls_cols
    return merged[out_cols].sort_values("integrated_lls", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Calibration curves figure
# ---------------------------------------------------------------------------
def plot_calibration_curves(summary: pd.DataFrame) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    layers = sorted(summary["layer"].unique())
    n = len(layers)
    ncols = 3
    nrows = int(np.ceil(n / ncols))

    with PdfPages(LLS_CURVES_PDF) as pdf:
        # Page 1: per-layer LLS-vs-bin curves
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(ncols * 4.0, nrows * 3.0),
                                 squeeze=False)
        for ax, layer in zip(axes.flat, layers):
            s = summary[summary["layer"] == layer].sort_values("bin")
            ax.plot(s["bin"], s["lls"], marker="o", color="steelblue")
            ax.axhline(0, color="grey", linewidth=0.8, linestyle="--")
            ax.set_title(f"{layer} (n={int(s['n_edges'].sum()):,})", fontsize=10)
            ax.set_xlabel("Score bin (low -> high)")
            ax.set_ylabel("LLS = log P(+|bin)/P(-|bin)")
            ax.grid(alpha=0.25)
        for ax in axes.flat[n:]:
            ax.set_visible(False)
        fig.suptitle("Per-layer FunCoup-style LLS calibration",
                     fontsize=12, y=1.02)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # Page 2: combined LLS curves
        fig, ax = plt.subplots(figsize=(7, 5))
        for layer in layers:
            s = summary[summary["layer"] == layer].sort_values("bin")
            ax.plot(s["bin"], s["lls"], marker="o", label=layer)
        ax.axhline(0, color="grey", linewidth=0.8, linestyle="--")
        ax.set_xlabel("Score bin (low -> high)")
        ax.set_ylabel("LLS")
        ax.set_title("LLS curves across layers")
        ax.legend(fontsize=8, loc="best")
        ax.grid(alpha=0.25)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
    log(f"Saved calibration curves: {LLS_CURVES_PDF}")


# ---------------------------------------------------------------------------
# Benchmark: LLS composite vs noisy-OR composite
# ---------------------------------------------------------------------------
def _auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Mann-Whitney-based AUROC (no sklearn dependency)."""
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    all_scores = np.concatenate([pos, neg])
    ranks = pd.Series(all_scores).rank(method="average").to_numpy()
    rank_pos = ranks[:len(pos)].sum()
    return float((rank_pos - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def _fold_enrichment(top_pairs: set, gold_pairs: set, base_rate: float) -> float:
    if not top_pairs or base_rate <= 0:
        return float("nan")
    hits = len(top_pairs & gold_pairs)
    rate = hits / len(top_pairs)
    return rate / base_rate if base_rate > 0 else float("nan")


def build_gold_pairs_for_gene_set(gene_set: set[str],
                                  nodes: set[str]) -> set[tuple[str, str]]:
    """All within-set pairs restricted to the node set."""
    members = sorted(g for g in gene_set if g in nodes)
    pairs = set()
    for i in range(len(members)):
        for j in range(i + 1, len(members)):
            pairs.add(canonical_pair(members[i], members[j]))
    return pairs


def benchmark_lls_vs_noisyor(lls_composite: pd.DataFrame,
                             nodes: set[str]) -> pd.DataFrame:
    log("Benchmarking LLS vs noisy-OR composite ...")

    # Gold panels
    gold_panels: dict[str, set[tuple[str, str]]] = {}

    # (a) Clinical drug targets
    if POS_CONTROL_CSV.exists():
        pc = pd.read_csv(POS_CONTROL_CSV)
        clin_col = next((c for c in pc.columns if c.lower() == "clinical_drug"), None)
        sym_col  = "Gene symbol" if "Gene symbol" in pc.columns else pc.columns[0]
        if clin_col:
            clin = pc.loc[pc[clin_col].notna() & (pc[clin_col].astype(str) != ""), sym_col]
            gold_panels["clinical_drug"] = build_gold_pairs_for_gene_set(set(clin), nodes)
    # (b) Govaere-25
    gold_panels["govaere25"]          = build_gold_pairs_for_gene_set(set(GOVAERE_25), nodes)
    # (c) Resmetirom pathway
    gold_panels["resmetirom_pathway"] = build_gold_pairs_for_gene_set(set(RESMETIROM_PATHWAY), nodes)

    for name, pairs in gold_panels.items():
        log(f"  Gold panel '{name}': {len(pairs):,} pairs in node set")

    # Load noisy-OR composite for comparison (if present)
    noisyor = None
    if COMPOSITE_CSV.exists():
        noisyor = pd.read_csv(COMPOSITE_CSV, usecols=["gene_a", "gene_b", "p_composite"])
        noisyor = canonicalize_frame(noisyor).drop_duplicates(subset=["gene_a", "gene_b"])
    else:
        log("  WARNING: noisy-OR composite missing — LLS benchmark will report LLS only.")

    # Universe: union of pairs in either composite
    lls_pairs = set(zip(lls_composite["gene_a"], lls_composite["gene_b"]))
    universe = set(lls_pairs)
    if noisyor is not None:
        universe |= set(zip(noisyor["gene_a"], noisyor["gene_b"]))

    rows = []
    for name, gold in gold_panels.items():
        gold_in_univ = gold & universe
        if not gold_in_univ:
            log(f"  Panel '{name}' has 0 gold pairs in universe — skipping.")
            continue
        base_rate = len(gold_in_univ) / max(len(universe), 1)

        # ---- LLS composite
        lls_idx = lls_composite.index
        lls_labels = np.zeros(len(lls_composite), dtype=int)
        lls_pairs_iter = list(zip(lls_composite["gene_a"], lls_composite["gene_b"]))
        for k, p in enumerate(lls_pairs_iter):
            if p in gold:
                lls_labels[k] = 1
        lls_scores = lls_composite["integrated_lls"].to_numpy()
        auroc_lls = _auroc(lls_scores, lls_labels)
        for top_n in (100, 1_000, 10_000):
            top = lls_composite.nlargest(top_n, "integrated_lls")
            top_pairs = set(zip(top["gene_a"], top["gene_b"]))
            fe = _fold_enrichment(top_pairs, gold, base_rate)
            rows.append({
                "panel":        name,
                "method":       "LLS",
                "top_n":        top_n,
                "hits_top":     len(top_pairs & gold),
                "fold_enrich":  fe,
                "auroc":        auroc_lls,
                "n_gold_univ":  len(gold_in_univ),
                "n_universe":   len(universe),
            })

        # ---- Noisy-OR composite
        if noisyor is not None:
            nor_labels = np.zeros(len(noisyor), dtype=int)
            nor_pairs_iter = list(zip(noisyor["gene_a"], noisyor["gene_b"]))
            for k, p in enumerate(nor_pairs_iter):
                if p in gold:
                    nor_labels[k] = 1
            nor_scores = noisyor["p_composite"].to_numpy()
            auroc_nor = _auroc(nor_scores, nor_labels)
            for top_n in (100, 1_000, 10_000):
                top = noisyor.nlargest(top_n, "p_composite")
                top_pairs = set(zip(top["gene_a"], top["gene_b"]))
                fe = _fold_enrichment(top_pairs, gold, base_rate)
                rows.append({
                    "panel":        name,
                    "method":       "noisy_OR",
                    "top_n":        top_n,
                    "hits_top":     len(top_pairs & gold),
                    "fold_enrich":  fe,
                    "auroc":        auroc_nor,
                    "n_gold_univ":  len(gold_in_univ),
                    "n_universe":   len(universe),
                })
    cmp_df = pd.DataFrame(rows)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    cmp_df.to_csv(LLS_CMP_CSV, index=False)
    log(f"Saved benchmark comparison: {LLS_CMP_CSV}")
    log("\nBenchmark summary:")
    if not cmp_df.empty:
        log(cmp_df.to_string(index=False))
    return cmp_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    log("=== 276: FunCoup-style LLS Calibration ===")
    log(f"BASE: {BASE}")
    log(f"Network dir: {NET_DIR}")
    log(f"Figure dir:  {FIG_DIR}")

    check_preconditions()

    # Node set
    nodes_df = pd.read_csv(NODES_CSV, usecols=["human_symbol"])
    nodes_df = nodes_df.dropna()
    nodes = set(nodes_df["human_symbol"].astype(str).tolist())
    log(f"Loaded {len(nodes):,} network nodes")

    rng = np.random.default_rng(RANDOM_SEED)
    positives, negatives, gold_source = build_gold_standard(nodes, rng)
    log(f"Gold standard: positives={len(positives):,} negatives={len(negatives):,} "
        f"source={gold_source}")

    # Load layer edges
    log("Loading layer edge files ...")
    layer_edges: dict[str, pd.DataFrame] = {}
    for layer in LAYER_NAMES:
        df = load_layer(layer)
        if len(df) == 0:
            log(f"  Layer {layer}: 0 rows after canonicalisation — skipping.")
            continue
        layer_edges[layer] = df
        log(f"  Layer {layer:8s}: {len(df):,} edges")

    n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
    results, summary = calibrate_all_layers(layer_edges, positives, negatives, n_workers)

    # Persist per-layer summary
    NET_DIR.mkdir(parents=True, exist_ok=True)
    summary_out = summary.copy()
    summary_out["gold_source"] = gold_source
    summary_out["n_gold_pos"] = len(positives)
    summary_out["n_gold_neg"] = len(negatives)
    summary_out.to_csv(LLS_SUMMARY_CSV, index=False)
    log(f"Saved calibration summary: {LLS_SUMMARY_CSV}")

    # Persist per-layer LLS edges
    write_layer_lls(results)

    # Composite integrated LLS
    composite = build_composite(results)
    composite.to_csv(LLS_COMPOSITE_CSV, index=False)
    log(f"Saved integrated LLS composite: {LLS_COMPOSITE_CSV} "
        f"({len(composite):,} pairs)")

    # Figures
    plot_calibration_curves(summary)

    # Benchmark vs noisy-OR
    benchmark_lls_vs_noisyor(composite, nodes)

    elapsed = time.time() - t0_global
    log(f"\n=== 276 Complete in {elapsed:.1f}s ({elapsed/60:.1f} min) ===")
    log(f"Layers calibrated:       {len(results)}")
    log(f"Integrated-LLS pairs:    {len(composite):,}")
    log(f"Max integrated_LLS:      {composite['integrated_lls'].max():.3f}")
    log(f"Pairs with LLS>0:        {(composite['integrated_lls'] > 0).sum():,}")
    log(f"Pairs with k_layers>=2:  {(composite['k_contributing_layers'] >= 2).sum():,}")


if __name__ == "__main__":
    main()
