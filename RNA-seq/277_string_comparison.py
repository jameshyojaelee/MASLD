#!/usr/bin/env python
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_277_string
#SBATCH --output=logs/net_277_string_%j.out
#SBATCH --error=logs/net_277_string_%j.err
#
# 277_string_comparison.py
# ---------------------------------------------------------------------------
# Head-to-head comparison of our MASLD multi-evidence composite network vs.
# STRING v12 (generic human interactome) for MASLD gene-set recovery.
#
# Purpose: quantitative evidence for the "novelty vs. prior art" section of
# the Cell Metabolism Resource paper. The adversarial review challenged:
#   "FunCoup / HumanNet / STRING already do multi-modal integration -- how is
#    ours better for MASLD?"
# This script answers by benchmarking both networks on the SAME node set and
# SAME gold standards, computing:
#   - fraction of gold-gold pairs connected  (pairwise recall)
#   - mean degree of gold genes vs background (hub prominence + Wilcoxon U)
#   - ROC-AUC for gold-gold pair classification
#   - precision @ top-k (k = 1000, 5000, 10000) by edge weight
#   - fold-over-random with degree-preserving bootstrap 95% CIs + permutation p
#   - Jaccard overlap of edge sets (is composite just a STRING subset?)
#
# Environment: micromamba activate spatial
#
# Inputs:
#   data/string_ppi/9606.protein.links.v12.0.txt.gz
#   data/string_ppi/9606.protein.info.v12.0.txt.gz
#   RNA-seq/results/network/composite_edges.csv
#   RNA-seq/results/network/network_nodes.csv
#   results/library/positive_control.csv
#
# Outputs:
#   RNA-seq/results/network/string_v12_normalized.csv
#   figures/misc/network_vs_string/string_comparison_table.csv
#   figures/misc/network_vs_string/string_comparison_figure.pdf
#   figures/misc/network_vs_string/jaccard_overlap.csv
#   figures/misc/network_vs_string/roc_curves_data.csv
# ---------------------------------------------------------------------------

from __future__ import annotations

import argparse
import gzip
import os
import sys
import time
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
STRING_LINKS = BASE / "data/string_ppi/9606.protein.links.v12.0.txt.gz"
STRING_INFO = BASE / "data/string_ppi/9606.protein.info.v12.0.txt.gz"
NETDIR = BASE / "RNA-seq/results/network"
COMPOSITE_CSV = NETDIR / "composite_edges.csv"
NODES_CSV = NETDIR / "network_nodes.csv"
POS_CTRL_CSV = BASE / "results/library/positive_control.csv"

FIG_DIR = BASE / "figures/misc/network_vs_string"
FIG_DIR.mkdir(parents=True, exist_ok=True)

STRING_NORM_CSV = NETDIR / "string_v12_normalized.csv"
TABLE_CSV = FIG_DIR / "string_comparison_table.csv"
FIGURE_PDF = FIG_DIR / "string_comparison_figure.pdf"
JACCARD_CSV = FIG_DIR / "jaccard_overlap.csv"
ROC_CSV = FIG_DIR / "roc_curves_data.csv"

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
STRING_THRESHOLDS = {
    "ge400": 400,  # medium confidence
    "ge700": 700,  # high confidence
}
COMPOSITE_THRESHOLDS = {
    "p_ge_0.5": 0.5,
    "p_ge_0.7": 0.7,
    "p_ge_0.9": 0.9,
    "p_ge_0.95": 0.95,
}
TOP_K_VALUES = [1000, 5000, 10000]
N_PERMUTATIONS = 1000
N_BOOTSTRAP = 1000
SEED = 42
N_CORES = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
if N_CORES < 1:
    N_CORES = 8


# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------
def preflight_check() -> None:
    """Verify all required inputs exist. Exit loudly if any are missing."""
    missing = []
    for f in (
        STRING_LINKS,
        STRING_INFO,
        COMPOSITE_CSV,
        NODES_CSV,
        POS_CTRL_CSV,
    ):
        if not f.exists():
            missing.append(str(f))
    if missing:
        print(
            "ERROR: required input file(s) missing. Cannot continue.",
            file=sys.stderr,
        )
        for m in missing:
            print(f"  [missing] {m}", file=sys.stderr)
        if str(STRING_INFO) in missing or str(STRING_LINKS) in missing:
            print(
                "\nDownload STRING v12 human files to "
                f"{STRING_LINKS.parent} with:\n"
                "  wget https://stringdb-downloads.org/download/"
                "protein.links.v12.0/9606.protein.links.v12.0.txt.gz\n"
                "  wget https://stringdb-downloads.org/download/"
                "protein.info.v12.0/9606.protein.info.v12.0.txt.gz",
                file=sys.stderr,
            )
        sys.exit(1)


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
def load_nodes() -> pd.DataFrame:
    nodes = pd.read_csv(NODES_CSV, low_memory=False)
    if "human_symbol" not in nodes.columns:
        raise ValueError("network_nodes.csv missing 'human_symbol' column")
    return nodes


def load_string(node_set: set[str]) -> pd.DataFrame:
    """Load STRING v12, map ENSP->symbol, filter to genes in our node set.

    Returns dataframe with columns (gene_a, gene_b, string_score, string_ge400,
    string_ge700). Pairs are unordered: gene_a < gene_b by string sort.
    """
    print(f"[string] reading info: {STRING_INFO}", flush=True)
    info = pd.read_csv(
        STRING_INFO,
        sep="\t",
        compression="gzip",
        usecols=["#string_protein_id", "preferred_name"],
    )
    info = info.rename(
        columns={"#string_protein_id": "ensp", "preferred_name": "symbol"}
    )
    ensp2sym = dict(zip(info["ensp"], info["symbol"]))
    print(f"[string] {len(ensp2sym):,} ENSP -> symbol mappings", flush=True)

    print(f"[string] reading links: {STRING_LINKS}", flush=True)
    links = pd.read_csv(
        STRING_LINKS,
        sep=" ",
        compression="gzip",
    )
    links.columns = [c.strip() for c in links.columns]
    print(f"[string] raw edges: {len(links):,}", flush=True)

    links["gene_a_raw"] = links["protein1"].map(ensp2sym)
    links["gene_b_raw"] = links["protein2"].map(ensp2sym)
    links = links.dropna(subset=["gene_a_raw", "gene_b_raw"])
    print(f"[string] edges w/ symbols: {len(links):,}", flush=True)

    in_nodes = links["gene_a_raw"].isin(node_set) & links["gene_b_raw"].isin(
        node_set
    )
    links = links.loc[in_nodes].copy()
    print(
        f"[string] edges both-in-node-set: {len(links):,}",
        flush=True,
    )

    # Drop self-loops; normalize pair order so gene_a < gene_b
    links = links.loc[links["gene_a_raw"] != links["gene_b_raw"]].copy()
    pair = np.sort(links[["gene_a_raw", "gene_b_raw"]].to_numpy(), axis=1)
    links["gene_a"] = pair[:, 0]
    links["gene_b"] = pair[:, 1]

    # STRING links file is already pair-symmetric; collapse to one row per pair
    # by taking max combined_score (defensive; it should be identical).
    links = (
        links.groupby(["gene_a", "gene_b"], as_index=False)["combined_score"]
        .max()
        .rename(columns={"combined_score": "string_score"})
    )
    links["string_ge400"] = links["string_score"] >= STRING_THRESHOLDS["ge400"]
    links["string_ge700"] = links["string_score"] >= STRING_THRESHOLDS["ge700"]
    print(
        f"[string] unique pairs: {len(links):,} | "
        f"ge400: {links['string_ge400'].sum():,} | "
        f"ge700: {links['string_ge700'].sum():,}",
        flush=True,
    )
    links.to_csv(STRING_NORM_CSV, index=False)
    print(f"[string] wrote {STRING_NORM_CSV}", flush=True)
    return links


def load_composite(node_set: set[str]) -> pd.DataFrame:
    print(f"[composite] reading {COMPOSITE_CSV}", flush=True)
    comp = pd.read_csv(
        COMPOSITE_CSV, usecols=["gene_a", "gene_b", "p_composite"]
    )
    # Filter to node set and normalize pair order
    in_nodes = comp["gene_a"].isin(node_set) & comp["gene_b"].isin(node_set)
    comp = comp.loc[in_nodes].copy()
    comp = comp.loc[comp["gene_a"] != comp["gene_b"]].copy()
    pair = np.sort(comp[["gene_a", "gene_b"]].to_numpy(), axis=1)
    comp["gene_a"] = pair[:, 0]
    comp["gene_b"] = pair[:, 1]
    comp = comp.groupby(["gene_a", "gene_b"], as_index=False)[
        "p_composite"
    ].max()
    for name, thr in COMPOSITE_THRESHOLDS.items():
        comp[f"composite_{name}"] = comp["p_composite"] >= thr
    print(
        f"[composite] unique pairs: {len(comp):,} | "
        + " | ".join(
            f"{name}: {comp[f'composite_{name}'].sum():,}"
            for name in COMPOSITE_THRESHOLDS
        ),
        flush=True,
    )
    return comp


# ---------------------------------------------------------------------------
# Gold standards
# ---------------------------------------------------------------------------
GOVAERE25 = [
    "AKR1B10", "SOX4", "TTC39A", "VIL1", "CLIC6", "ABCB4", "LTBP2", "STMN2",
    "AGMO", "ANGPTL4", "CA3", "CCBE1", "CFB", "CSAD", "FABP4", "FADS1",
    "IGFBP7", "MME", "PANX2", "PLS1", "RP11-1399P15.1", "SDC1", "SRPX",
    "TMEM9", "TPBG",
]

RESMETIROM_PATHWAY = ["THRB", "PCSK9", "PPARA", "APOA1", "CYP7A1"]


def build_gold_sets(nodes: pd.DataFrame) -> Dict[str, List[str]]:
    all_genes = set(nodes["human_symbol"].unique())
    pc = pd.read_csv(POS_CTRL_CSV)

    # clinical_drug: any row with a non-empty Clinical_Drug
    clin_col = None
    for c in ("Clinical_Drug", "clinical_drug"):
        if c in pc.columns:
            clin_col = c
            break
    if clin_col:
        clin_mask = pc[clin_col].astype(str).str.strip().replace("nan", "") != ""
        clin_mask &= pc[clin_col].notna()
        gold_clin = pc.loc[clin_mask, "Gene symbol"].dropna().unique().tolist()
    else:
        gold_clin = []

    gold_clin = sorted([g for g in gold_clin if g in all_genes])
    gold_gov = sorted([g for g in GOVAERE25 if g in all_genes])
    gold_resm = sorted([g for g in RESMETIROM_PATHWAY if g in all_genes])

    # All positive controls (union across all columns of positive_control.csv)
    gold_pc = sorted(
        [
            g
            for g in pc["Gene symbol"].dropna().unique().tolist()
            if g in all_genes
        ]
    )

    # SuSiE-COLOC genes: PP.H4 > 0.5 from nodes$coloc_susie_best_pp4
    if "coloc_susie_best_pp4" in nodes.columns:
        pp4 = pd.to_numeric(nodes["coloc_susie_best_pp4"], errors="coerce")
        gold_susie = sorted(
            nodes.loc[pp4 > 0.5, "human_symbol"].dropna().unique().tolist()
        )
    else:
        gold_susie = []
        print(
            "[gold] WARNING: 'coloc_susie_best_pp4' not in nodes; "
            "skipping SuSiE-COLOC gold set",
            flush=True,
        )

    gold_sets = {
        "clinical_drug": gold_clin,
        "govaere25": gold_gov,
        "resmetirom_pathway": gold_resm,
        "positive_controls_all": gold_pc,
        "susie_coloc_pp4_gt_0.5": gold_susie,
    }
    for k, v in gold_sets.items():
        print(f"[gold] {k}: n={len(v)}", flush=True)
    return gold_sets


# ---------------------------------------------------------------------------
# Per-network edge accessors
# ---------------------------------------------------------------------------
@dataclass
class EdgeSet:
    """One (network, threshold) slice, with score vector for ranking."""

    network: str
    threshold_label: str
    edges: pd.DataFrame  # cols: gene_a, gene_b, score
    full_edges: pd.DataFrame  # full weighted edge list (all edges of network)

    @property
    def pair_set(self) -> set[tuple[str, str]]:
        return set(map(tuple, self.edges[["gene_a", "gene_b"]].to_numpy()))

    @property
    def n_edges(self) -> int:
        return len(self.edges)


def build_edge_sets(
    string_df: pd.DataFrame, composite_df: pd.DataFrame
) -> List[EdgeSet]:
    out: List[EdgeSet] = []
    for label, thr in STRING_THRESHOLDS.items():
        sub = string_df.loc[string_df["string_score"] >= thr, :].copy()
        sub = sub[["gene_a", "gene_b", "string_score"]].rename(
            columns={"string_score": "score"}
        )
        full = string_df[["gene_a", "gene_b", "string_score"]].rename(
            columns={"string_score": "score"}
        )
        out.append(
            EdgeSet(
                network="STRING_v12",
                threshold_label=label,
                edges=sub,
                full_edges=full,
            )
        )
    for label, thr in COMPOSITE_THRESHOLDS.items():
        sub = composite_df.loc[composite_df["p_composite"] >= thr, :].copy()
        sub = sub[["gene_a", "gene_b", "p_composite"]].rename(
            columns={"p_composite": "score"}
        )
        full = composite_df[["gene_a", "gene_b", "p_composite"]].rename(
            columns={"p_composite": "score"}
        )
        out.append(
            EdgeSet(
                network="our_composite",
                threshold_label=label,
                edges=sub,
                full_edges=full,
            )
        )
    return out


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def compute_degrees(edges: pd.DataFrame, all_genes: Sequence[str]) -> pd.Series:
    """Gene degree given an edge table with gene_a, gene_b."""
    deg = pd.concat(
        [edges["gene_a"], edges["gene_b"]], ignore_index=True
    ).value_counts()
    return deg.reindex(all_genes, fill_value=0).astype(int)


def gold_pair_fraction_connected(
    pair_set: set[tuple[str, str]], gold: Sequence[str]
) -> Tuple[int, int, float]:
    """Fraction of (unordered) gold-gold pairs that appear as edges."""
    gold_sorted = sorted(set(gold))
    n_pairs = 0
    hits = 0
    for i in range(len(gold_sorted)):
        for j in range(i + 1, len(gold_sorted)):
            a, b = gold_sorted[i], gold_sorted[j]
            n_pairs += 1
            if (a, b) in pair_set:
                hits += 1
    frac = hits / n_pairs if n_pairs else float("nan")
    return hits, n_pairs, frac


def _perm_worker(args):
    rng_seed, all_genes, gold_size, pair_set, n_iters = args
    rng = np.random.default_rng(rng_seed)
    genes_arr = np.asarray(all_genes)
    out = np.empty(n_iters, dtype=np.float64)
    for i in range(n_iters):
        sample = rng.choice(genes_arr, size=gold_size, replace=False)
        s = np.sort(sample)
        hits = 0
        n_pairs = 0
        for a_i in range(len(s)):
            a = s[a_i]
            for b_i in range(a_i + 1, len(s)):
                b = s[b_i]
                n_pairs += 1
                if (a, b) in pair_set:
                    hits += 1
        out[i] = hits / n_pairs if n_pairs else 0.0
    return out


def permutation_null_fraction(
    pair_set: set[tuple[str, str]],
    gold_size: int,
    all_genes: Sequence[str],
    n_perm: int = N_PERMUTATIONS,
    n_cores: int = N_CORES,
    seed: int = SEED,
) -> np.ndarray:
    """Null: draw random gene sets of same size; compute pair fraction.

    Returns array of null fractions.
    """
    if gold_size < 2:
        return np.zeros(n_perm)
    per_worker = max(1, n_perm // n_cores)
    jobs = []
    base_seed = np.random.SeedSequence(seed).spawn(n_cores)
    remaining = n_perm
    for i in range(n_cores):
        k = per_worker if i < n_cores - 1 else remaining
        jobs.append(
            (int(base_seed[i].generate_state(1)[0]), all_genes, gold_size,
             pair_set, k)
        )
        remaining -= k
        if remaining <= 0:
            break
    with Pool(min(n_cores, len(jobs))) as p:
        chunks = p.map(_perm_worker, jobs)
    return np.concatenate(chunks)


def degree_based_metrics(
    deg: pd.Series, gold: Sequence[str]
) -> Dict[str, float]:
    gold = [g for g in gold if g in deg.index]
    if len(gold) < 2:
        return {
            "mean_degree_gold": float("nan"),
            "mean_degree_bg": float("nan"),
            "wilcoxon_U": float("nan"),
            "wilcoxon_p": float("nan"),
        }
    gold_mask = deg.index.isin(gold)
    gold_deg = deg.loc[gold_mask].values
    bg_deg = deg.loc[~gold_mask].values
    U, p = stats.mannwhitneyu(gold_deg, bg_deg, alternative="greater")
    return {
        "mean_degree_gold": float(np.mean(gold_deg)),
        "mean_degree_bg": float(np.mean(bg_deg)),
        "wilcoxon_U": float(U),
        "wilcoxon_p": float(p),
    }


def roc_auc_gold_pair(
    edges: pd.DataFrame,
    gold: Sequence[str],
    all_genes: Sequence[str],
    n_neg_per_pos: int = 10,
    seed: int = SEED,
) -> Tuple[float, Optional[np.ndarray], Optional[np.ndarray]]:
    """ROC-AUC for predicting gold-gold vs non-gold pairs using edge score.

    Positives: gold-gold pairs that appear as edges.
    Negatives: non-gold pairs sampled from edges (so absence isn't trivial).
    Returns (auc, fpr, tpr).
    """
    gold_set = set(gold)
    if len(gold_set) < 2 or len(edges) == 0:
        return float("nan"), None, None
    g_a = edges["gene_a"].values
    g_b = edges["gene_b"].values
    is_gold_gold = np.isin(g_a, list(gold_set)) & np.isin(g_b, list(gold_set))
    scores = edges["score"].values
    n_pos = int(is_gold_gold.sum())
    if n_pos < 2:
        return float("nan"), None, None
    rng = np.random.default_rng(seed)
    neg_idx = np.where(~is_gold_gold)[0]
    n_neg = min(len(neg_idx), n_pos * n_neg_per_pos)
    sampled_neg = rng.choice(neg_idx, size=n_neg, replace=False)
    pos_idx = np.where(is_gold_gold)[0]
    keep_idx = np.concatenate([pos_idx, sampled_neg])
    y = np.zeros(len(keep_idx), dtype=int)
    y[: len(pos_idx)] = 1
    s = scores[keep_idx]
    try:
        auc = roc_auc_score(y, s)
        fpr, tpr, _ = roc_curve(y, s)
    except ValueError:
        return float("nan"), None, None
    return float(auc), fpr, tpr


def precision_at_top_k(
    full_edges: pd.DataFrame, gold: Sequence[str], ks: Sequence[int]
) -> Dict[int, float]:
    """Precision @ top-k: of the top-k edges by score, fraction that connect
    two gold genes."""
    gold_set = set(gold)
    ordered = full_edges.sort_values("score", ascending=False)
    ga = ordered["gene_a"].values
    gb = ordered["gene_b"].values
    is_gg = np.isin(ga, list(gold_set)) & np.isin(gb, list(gold_set))
    out = {}
    cum = np.cumsum(is_gg)
    for k in ks:
        k_eff = min(k, len(ordered))
        if k_eff == 0:
            out[k] = float("nan")
        else:
            out[k] = float(cum[k_eff - 1] / k_eff)
    return out


def bootstrap_fraction_ci(
    pair_set: set[tuple[str, str]],
    gold: Sequence[str],
    n_boot: int = N_BOOTSTRAP,
    seed: int = SEED,
) -> Tuple[float, float]:
    """95% CI for fraction_connected by resampling gold genes with replacement."""
    gold = list(gold)
    if len(gold) < 2:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for b in range(n_boot):
        samp = rng.choice(gold, size=len(gold), replace=True)
        s = np.unique(samp)
        if len(s) < 2:
            boots[b] = float("nan")
            continue
        n_pairs = 0
        hits = 0
        s = np.sort(s)
        for i in range(len(s)):
            for j in range(i + 1, len(s)):
                a, b_ = s[i], s[j]
                n_pairs += 1
                if (a, b_) in pair_set:
                    hits += 1
        boots[b] = hits / n_pairs if n_pairs else float("nan")
    boots = boots[~np.isnan(boots)]
    if len(boots) == 0:
        return (float("nan"), float("nan"))
    return (float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975)))


# ---------------------------------------------------------------------------
# Main benchmark
# ---------------------------------------------------------------------------
def run_benchmark(
    edge_sets: List[EdgeSet],
    gold_sets: Dict[str, List[str]],
    all_genes: List[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    roc_rows = []
    for es in edge_sets:
        t0 = time.time()
        print(
            f"\n[bench] {es.network} / {es.threshold_label} "
            f"(n_edges={es.n_edges:,})",
            flush=True,
        )
        pair_set = es.pair_set
        deg = compute_degrees(es.edges, all_genes)
        for gold_name, gold in gold_sets.items():
            gold_sz = len(gold)
            if gold_sz < 2:
                print(f"  [skip] {gold_name}: fewer than 2 gold genes")
                continue

            hits, n_pairs, frac = gold_pair_fraction_connected(pair_set, gold)
            null = permutation_null_fraction(
                pair_set, gold_sz, all_genes,
                n_perm=N_PERMUTATIONS, n_cores=N_CORES,
            )
            null_mean = float(np.mean(null)) if len(null) else float("nan")
            fold = (
                frac / null_mean if null_mean and null_mean > 0 else float("inf")
            )
            p_perm = (
                float((null >= frac).sum() + 1) / (len(null) + 1)
                if len(null)
                else float("nan")
            )
            ci_lo, ci_hi = bootstrap_fraction_ci(pair_set, gold)

            deg_m = degree_based_metrics(deg, gold)
            auc, fpr, tpr = roc_auc_gold_pair(
                es.edges, gold, all_genes, n_neg_per_pos=10
            )
            prec_at_k = precision_at_top_k(es.full_edges, gold, TOP_K_VALUES)

            row = {
                "network": es.network,
                "threshold": es.threshold_label,
                "gold_set": gold_name,
                "n_gold": gold_sz,
                "n_edges_in_network": es.n_edges,
                "hits": hits,
                "n_gold_pairs": n_pairs,
                "fraction_connected": frac,
                "fraction_ci_lo": ci_lo,
                "fraction_ci_hi": ci_hi,
                "null_mean_fraction": null_mean,
                "fold_over_random": fold,
                "p_perm": p_perm,
                **deg_m,
                "auc_roc": auc,
            }
            for k, v in prec_at_k.items():
                row[f"precision_at_{k}"] = v
            rows.append(row)

            # Store ROC data for plotting
            if fpr is not None and tpr is not None:
                sample_n = min(500, len(fpr))
                if sample_n < len(fpr):
                    idx = np.linspace(0, len(fpr) - 1, sample_n).astype(int)
                    fpr_s = fpr[idx]
                    tpr_s = tpr[idx]
                else:
                    fpr_s = fpr
                    tpr_s = tpr
                for f_, t_ in zip(fpr_s, tpr_s):
                    roc_rows.append(
                        {
                            "network": es.network,
                            "threshold": es.threshold_label,
                            "gold_set": gold_name,
                            "auc": auc,
                            "fpr": float(f_),
                            "tpr": float(t_),
                        }
                    )
        print(f"  [done] in {time.time() - t0:.1f}s", flush=True)

    return pd.DataFrame(rows), pd.DataFrame(roc_rows)


# ---------------------------------------------------------------------------
# Jaccard
# ---------------------------------------------------------------------------
def compute_jaccard(
    string_df: pd.DataFrame, composite_df: pd.DataFrame
) -> pd.DataFrame:
    jacc_rows = []
    comp_sets = {
        name: set(
            map(
                tuple,
                composite_df.loc[
                    composite_df["p_composite"] >= thr,
                    ["gene_a", "gene_b"],
                ].to_numpy(),
            )
        )
        for name, thr in COMPOSITE_THRESHOLDS.items()
    }
    string_sets = {
        name: set(
            map(
                tuple,
                string_df.loc[
                    string_df["string_score"] >= thr,
                    ["gene_a", "gene_b"],
                ].to_numpy(),
            )
        )
        for name, thr in STRING_THRESHOLDS.items()
    }
    for cname, cset in comp_sets.items():
        for sname, sset in string_sets.items():
            inter = len(cset & sset)
            uni = len(cset | sset)
            jacc = inter / uni if uni else float("nan")
            jacc_rows.append(
                {
                    "composite_threshold": cname,
                    "string_threshold": sname,
                    "n_composite": len(cset),
                    "n_string": len(sset),
                    "intersection": inter,
                    "union": uni,
                    "jaccard": jacc,
                    "frac_composite_in_string": (
                        inter / len(cset) if cset else float("nan")
                    ),
                    "frac_string_in_composite": (
                        inter / len(sset) if sset else float("nan")
                    ),
                }
            )
    return pd.DataFrame(jacc_rows)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def plot_figure(bench: pd.DataFrame, roc_df: pd.DataFrame) -> None:
    with PdfPages(FIGURE_PDF) as pdf:
        # Panel 1: Fold-over-random bar chart
        fig, ax = plt.subplots(figsize=(12, 6))
        gold_order = sorted(bench["gold_set"].unique())
        x = np.arange(len(gold_order))
        combos = (
            bench[["network", "threshold"]]
            .drop_duplicates()
            .sort_values(["network", "threshold"])
            .values
        )
        width = 0.8 / max(1, len(combos))
        for i, (net, thr) in enumerate(combos):
            sub = bench[(bench.network == net) & (bench.threshold == thr)]
            sub = sub.set_index("gold_set").reindex(gold_order)
            heights = sub["fold_over_random"].values.astype(float)
            heights = np.where(np.isinf(heights), np.nan, heights)
            color = (
                "#d95f02" if net == "STRING_v12" else "#1b9e77"
            )
            alpha = 0.5 + 0.5 * (i / max(1, len(combos) - 1))
            ax.bar(
                x + i * width - 0.4 + width / 2,
                heights,
                width=width,
                label=f"{net} / {thr}",
                color=color,
                alpha=alpha,
                edgecolor="black",
                linewidth=0.3,
            )
        ax.axhline(1, color="gray", linestyle="--", linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(gold_order, rotation=30, ha="right")
        ax.set_ylabel("Fold over random")
        ax.set_title("Fold-over-random: composite vs STRING v12")
        ax.legend(loc="upper right", fontsize=7, ncol=2)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Panel 2: Forest plot of fraction_connected with 95% bootstrap CI
        fig, ax = plt.subplots(figsize=(10, max(6, 0.35 * len(bench))))
        bench_sorted = bench.sort_values(["gold_set", "network", "threshold"])
        ypos = np.arange(len(bench_sorted))
        ax.errorbar(
            bench_sorted["fraction_connected"],
            ypos,
            xerr=np.vstack(
                [
                    bench_sorted["fraction_connected"].values
                    - bench_sorted["fraction_ci_lo"].values,
                    bench_sorted["fraction_ci_hi"].values
                    - bench_sorted["fraction_connected"].values,
                ]
            ),
            fmt="o",
            capsize=2,
            color="black",
            ecolor="gray",
        )
        ax.set_yticks(ypos)
        ax.set_yticklabels(
            [
                f"{r.gold_set} | {r.network} | {r.threshold}"
                for _, r in bench_sorted.iterrows()
            ],
            fontsize=7,
        )
        ax.set_xlabel("Fraction of gold-gold pairs connected (95% CI)")
        ax.set_title("Forest plot: gold-pair recall")
        ax.invert_yaxis()
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Panel 3: ROC curves per gold set
        for gs in sorted(roc_df["gold_set"].unique()):
            fig, ax = plt.subplots(figsize=(7, 7))
            sub = roc_df[roc_df["gold_set"] == gs]
            for (net, thr), g in sub.groupby(["network", "threshold"]):
                g = g.sort_values("fpr")
                auc = g["auc"].iloc[0]
                color = (
                    "#d95f02" if net == "STRING_v12" else "#1b9e77"
                )
                ax.plot(
                    g["fpr"],
                    g["tpr"],
                    label=f"{net}/{thr} (AUC={auc:.3f})",
                    color=color,
                    alpha=0.7,
                )
            ax.plot([0, 1], [0, 1], "--", color="gray", linewidth=0.7)
            ax.set_xlabel("FPR")
            ax.set_ylabel("TPR")
            ax.set_title(f"ROC: gold-gold pair classification — {gs}")
            ax.legend(fontsize=7, loc="lower right")
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        # Panel 4: precision @ top-k
        fig, axes = plt.subplots(1, len(TOP_K_VALUES), figsize=(5 * len(TOP_K_VALUES), 5))
        if len(TOP_K_VALUES) == 1:
            axes = [axes]
        for ax, k in zip(axes, TOP_K_VALUES):
            piv = bench.pivot_table(
                index="gold_set",
                columns=["network", "threshold"],
                values=f"precision_at_{k}",
            )
            piv.plot(kind="bar", ax=ax)
            ax.set_title(f"Precision @ top-{k}")
            ax.set_ylabel("Precision")
            ax.legend(fontsize=6, loc="upper right")
            ax.tick_params(axis="x", rotation=30)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    print(f"[fig] wrote {FIGURE_PDF}", flush=True)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main() -> int:
    global N_PERMUTATIONS, N_BOOTSTRAP, N_CORES
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_perm", type=int, default=N_PERMUTATIONS)
    ap.add_argument("--n_boot", type=int, default=N_BOOTSTRAP)
    ap.add_argument("--n_cores", type=int, default=N_CORES)
    args = ap.parse_args()
    N_PERMUTATIONS = args.n_perm
    N_BOOTSTRAP = args.n_boot
    N_CORES = args.n_cores

    print("=" * 72)
    print("277 STRING v12 vs our composite head-to-head benchmark")
    print(f"  N_PERMUTATIONS={N_PERMUTATIONS} | N_BOOTSTRAP={N_BOOTSTRAP} | "
          f"N_CORES={N_CORES}")
    print("=" * 72)
    preflight_check()

    nodes = load_nodes()
    all_genes = sorted(nodes["human_symbol"].dropna().unique().tolist())
    node_set = set(all_genes)
    print(f"[nodes] {len(all_genes):,} genes", flush=True)

    gold_sets = build_gold_sets(nodes)

    string_df = load_string(node_set)
    composite_df = load_composite(node_set)

    edge_sets = build_edge_sets(string_df, composite_df)

    bench, roc_df = run_benchmark(edge_sets, gold_sets, all_genes)
    bench.to_csv(TABLE_CSV, index=False)
    roc_df.to_csv(ROC_CSV, index=False)
    print(f"[out] wrote {TABLE_CSV}", flush=True)
    print(f"[out] wrote {ROC_CSV}", flush=True)

    jacc = compute_jaccard(string_df, composite_df)
    jacc.to_csv(JACCARD_CSV, index=False)
    print(f"[out] wrote {JACCARD_CSV}", flush=True)
    print("\n[jaccard summary]")
    print(jacc.to_string(index=False))

    plot_figure(bench, roc_df)

    print("\n[summary table]")
    show = bench[
        [
            "network", "threshold", "gold_set",
            "fraction_connected", "fold_over_random",
            "p_perm", "auc_roc", "precision_at_1000",
        ]
    ]
    print(show.to_string(index=False))
    print("\n[done]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
