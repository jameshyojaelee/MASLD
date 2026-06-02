#!/usr/bin/env python3
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_292_f2_leiden
#SBATCH --output=logs/net_292_f2_leiden_%j.out
#SBATCH --error=logs/net_292_f2_leiden_%j.err
# ===========================================================================
# Script 292: F2-switch weighted Leiden community detection
# ===========================================================================
# Purpose:
#   Run Leiden TWICE on the STRING backbone (raw_score >= 0.7), weighted by
#   per-stage DE magnitude for F0-F1 ("pre-switch") vs F3-F4 ("post-switch").
#   Identify dissolving (F0-F1 only) and emerging (F3-F4 only) communities
#   to confirm the metabolic-to-inflammatory switch at F2.
#
# DE source used (documented here):
#   Stage-specific F0-F1 vs F3-F4 mega-analysis DE is not available in a
#   single pre-computed file. The progression transition_fib_dream_results.csv
#   provides adjacent-stage contrasts (F0->F1, F1->F2, F2->F3, F3->F4), and
#   the atlas column `f2_inflection_logFC` contrasts F0-F2 vs F3-F4 directly.
#
#   We therefore construct two stage-weight vectors from `f2_inflection_logFC`
#   (atlas) which is the closest proxy for the F0-F1 vs F3-F4 axis:
#     - W_F0F1[gene] = max(-logFC, 0)    (magnitude of DOWN-regulation at F3-F4;
#                                          i.e., genes higher in F0-F1)
#     - W_F3F4[gene] = max(+logFC, 0)    (magnitude of UP-regulation at F3-F4)
#   For each edge (i,j):
#     weight_F0F1[i,j] = string[i,j] * exp(0.5 * (W_F0F1[i] + W_F0F1[j]))
#     weight_F3F4[i,j] = string[i,j] * exp(0.5 * (W_F3F4[i] + W_F3F4[j]))
#   Genes absent from the DE table receive weight multiplier exp(0) = 1.
#
# Inputs:
#   RNA-seq/results/network/edges_ppi.csv          (gene_a, gene_b, raw_score)
#   RNA-seq/results/network/network_nodes.csv
#   RNA-seq/results/multi_evidence/multi_evidence_atlas.csv  (f2_inflection_logFC)
#
# Outputs (RNA-seq/results/network/communities_f2/):
#   communities_F01.csv, communities_F34.csv
#   community_transition_table.csv
#   community_enrichment_F01.csv, community_enrichment_F34.csv
#
# Environment: spatial (igraph, leidenalg, pandas, numpy, scipy)
# ===========================================================================

import os
import sys
import time
import warnings
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import igraph as ig
import leidenalg
from scipy.stats import hypergeom

warnings.filterwarnings("ignore", category=FutureWarning)

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
NET_DIR = BASE / "RNA-seq" / "results" / "network"
ATLAS = BASE / "RNA-seq" / "results" / "multi_evidence" / "multi_evidence_atlas.csv"
OUT_DIR = NET_DIR / "communities_f2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STRING_THRESH = 0.7
LEIDEN_RES = 1.0
SEED = 42
MIN_COMM_SIZE = 10          # communities smaller than this are ignored for classification
DISSOLVE_THRESH = 0.5       # "no single target holds >50% of source members"
EMERGE_THRESH = 0.5         # "at least 50% members come from no-large-source"
ENRICH_P = 1e-5

t0_global = time.time()


def log(msg):
    elapsed = time.time() - t0_global
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


# ---------------------------------------------------------------------------
# MSigDB Hallmark gene lists (hardcoded subset for falsification criteria).
# Source: MSigDB v2023.2.Hs — HALLMARK_OXIDATIVE_PHOSPHORYLATION (200 genes),
# HALLMARK_FATTY_ACID_METABOLISM (158), HALLMARK_INFLAMMATORY_RESPONSE (200),
# HALLMARK_TNFA_SIGNALING_VIA_NFKB (200). Trimmed to canonical members.
# ---------------------------------------------------------------------------
HALLMARK = {
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION": {
        "ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E","ATP5MC1","ATP5MC2","ATP5MC3",
        "ATP5ME","ATP5MF","ATP5MG","ATP5PB","ATP5PD","ATP5PF","ATP5PO","ATP6AP1","ATP6V0B",
        "ATP6V0C","ATP6V0E1","ATP6V1C1","ATP6V1D","ATP6V1E1","ATP6V1F","ATP6V1G1","ATP6V1H",
        "COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX6C","COX7A2","COX7A2L","COX7B","COX7C",
        "COX8A","CYC1","CYCS","DLAT","DLD","ECH1","ETFA","ETFB","ETFDH","FH","FXN","GPX4",
        "GRPEL1","HADHA","HADHB","IDH1","IDH2","IDH3A","IDH3B","IDH3G","MDH1","MDH2","MFN2",
        "NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9",
        "NDUFAB1","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8",
        "NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS6","NDUFS7","NDUFS8",
        "NDUFV1","NDUFV2","OGDH","PDHA1","PDHB","PDHX","PDK4","POLR2F","PPA1","PRDX3","SDHA",
        "SDHB","SDHC","SDHD","SLC25A3","SLC25A4","SLC25A5","SLC25A6","SLC25A11","SLC25A12",
        "SUCLA2","SUCLG1","SUCLG2","TIMM10","TIMM13","TIMM17A","TIMM22","TIMM50","TIMM8B",
        "TIMM9","TOMM22","TOMM70","UQCR10","UQCR11","UQCRB","UQCRC1","UQCRC2","UQCRFS1",
        "UQCRH","UQCRQ","VDAC1","VDAC2","VDAC3","ACAA2","ACADM","ACADSB","ACADVL","ACAT1",
        "ACO2","AFG3L2","ALDH6A1","ATP1B1","BAX","CASP7","CPT1A","CS","DECR1","ECHS1","GLUD1",
        "GOT2","HCCS","HSD17B10","HSPA9","ISCA1","ISCU","LDHA","LDHB","LRPPRC","MAOB","MPC1",
        "MRPL11","MRPL15","MRPL34","MRPL35","MRPS11","MRPS12","MRPS15","MRPS22","MRPS30",
        "MRPS7","MTRR","MTX2","NNT","OAT","OPA1","PHB2","PHYH","POR","RHOT1","RHOT2","SLC25A1",
        "SLC25A20","SUPV3L1","SURF1","TCIRG1","UCP2","ACADS","ACAA1","ACAT2","HADH","IMMT",
        "MGST3","OXA1L","PMPCA","PMPCB","POLR2K","RETSAT","MAOA","HSD17B7","GLS","AIFM1",
        "MRPL32","MRPS33","MRPS34","MRPS2","DLST","ATP5IF1","ETFDH"
    },
    "HALLMARK_FATTY_ACID_METABOLISM": {
        "ACAA1","ACAA2","ACADL","ACADM","ACADS","ACADSB","ACADVL","ACAT1","ACAT2","ACOX1",
        "ACOX3","ACSL1","ACSL3","ACSL4","ACSL5","ACSM3","ADH1C","ADIPOR2","ADSL","AKR1C3",
        "ALAD","ALDH1A1","ALDH3A1","ALDH3A2","ALDH9A1","AOC3","APEX1","APOA1","APOA2",
        "AQP7","AUH","BCKDHB","BLVRA","BPHL","CA2","CA4","CBR1","CBR3","CD1D","CD36","CEL",
        "CIDEA","CPOX","CPT1A","CPT2","CRAT","CRYZ","CYP1A1","CYP4A11","CYP4A22","D2HGDH",
        "DECR1","DHCR24","DLD","DLST","ECH1","ECHS1","ECI1","ECI2","EHHADH","ELOVL5","ENO2",
        "ENO3","EPHX1","EPHX2","ERP29","ETFDH","FABP1","FASN","FH","G0S2","GAD2","GCDH",
        "GLUL","GPD1","GPX4","GRHPR","H2AFV","HADH","HADHA","HADHB","HCCS","HIBCH","HMGCL",
        "HMGCS1","HMGCS2","HPGD","HSD17B10","HSD17B11","HSD17B4","HSD17B7","HSDL2","HSPH1",
        "IDH1","IDH3A","IDH3G","IDI1","IL4I1","INMT","LDHA","LGALS1","LTC4S","MAOA","MAOB",
        "MCEE","MDH1","ME1","MECR","METAP1","MGLL","MIF","MLYCD","NBN","NCAPH2","NTHL1",
        "ODC1","OSTC","PCBD1","PCCA","PDHA1","PDHB","POR","PPARA","PRDX6","PSME1","PTGR1",
        "PTS","RAP1GDS1","RDH11","RDH16","REEP6","RETSAT","S100A10","SDHA","SDHB","SDHC",
        "SDHD","SERINC1","SMS","SUCLA2","SUCLG1","SUCLG2","TDO2","TP53INP2","UBE2L6","UGDH",
        "UROD","UROS","XIST","YWHAH"
    },
    "HALLMARK_INFLAMMATORY_RESPONSE": {
        "ABCA1","ABI1","ACVR1B","ACVR2A","ADGRE1","ADM","ADORA2B","ADRM1","AHR","APLNR",
        "AQP9","ATP2A2","ATP2B1","ATP2C1","AXL","BDKRB1","BEST1","BST2","BTG2","C3AR1",
        "C5AR1","CALCRL","CCL17","CCL2","CCL20","CCL22","CCL24","CCL5","CCL7","CCR7","CCRL2",
        "CD14","CD40","CD48","CD55","CD69","CD70","CD82","CDKN1A","CHST2","CLEC5A","CMKLR1",
        "CSF1","CSF3","CSF3R","CX3CL1","CXCL10","CXCL11","CXCL6","CXCL8","CXCL9","CXCR6",
        "CYBB","DCBLD2","EBI3","EDN1","EIF2AK2","EMP3","EREG","F3","FFAR2","FPR1","FZD5",
        "GABBR1","GCH1","GNA15","GNAI3","GP1BA","GPC3","GPR132","GPR183","HAS2","HBEGF",
        "HIF1A","HPN","HRH1","ICAM1","ICAM4","ICOSLG","IFITM1","IFNAR1","IFNGR2","IL10",
        "IL10RA","IL12B","IL15","IL15RA","IL18","IL18R1","IL18RAP","IL1A","IL1B","IL1R1",
        "IL2RB","IL4R","IL6","IL7R","INHBA","IRAK2","IRF1","IRF7","ITGA5","ITGB3","ITGB8",
        "KCNA3","KCNJ2","KCNMB2","KIF1B","KLF6","LAMP3","LCK","LCP2","LDLR","LIF","LPAR1",
        "LTA","LY6E","LYN","MARCO","MEFV","MEP1A","MET","MMP14","MSR1","MXD1","MYC","NAMPT",
        "NDP","NFKB1","NFKBIA","NLRP3","NMI","NMUR1","NOD2","NPFFR2","OLR1","OPRK1","OSM",
        "OSMR","P2RX4","P2RX7","P2RY2","PCDH7","PDE4B","PDPN","PIK3R5","PLAUR","PROK2",
        "PSEN1","PTAFR","PTGER2","PTGER4","PTGIR","PTPRE","PVR","RAF1","RASGRP1","RELA",
        "RGS1","RGS16","RHOG","RIPK2","RNF144B","ROS1","RTP4","S100A7","S100A8","S100A9",
        "SCARF1","SCN1B","SELE","SELENOS","SELL","SEMA4D","SERPINE1","SGMS2","SLAMF1",
        "SLC11A2","SLC1A2","SLC28A2","SLC31A1","SLC31A2","SLC4A4","SLC7A1","SLC7A2","SPHK1",
        "SRI","STAB1","TACR1","TACR3","TAPBP","TIMP1","TLR1","TLR2","TLR3","TNFAIP6","TNFRSF1B",
        "TNFRSF9","TNFSF10","TNFSF15","TPBG","VIP"
    },
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": {
        "ABCA1","ACKR3","AREG","ATF3","ATP2B1","B4GALT1","B4GALT5","BCL2A1","BCL3","BCL6",
        "BHLHE40","BIRC2","BIRC3","BMP2","BTG1","BTG2","BTG3","CCL2","CCL20","CCL4","CCL5",
        "CCN1","CCND1","CCNL1","CCRL2","CD44","CD69","CD80","CD83","CDKN1A","CEBPB","CEBPD",
        "CFLAR","CLCF1","CSF1","CSF2","CXCL1","CXCL10","CXCL11","CXCL2","CXCL3","CXCL6",
        "DDX58","DENND5A","DNAJB4","DRAM1","DUSP1","DUSP2","DUSP4","DUSP5","EDN1","EFNA1",
        "EGR1","EGR2","EGR3","EHD1","EIF1","ETS2","F2RL1","F3","FJX1","FOS","FOSB","FOSL1",
        "FOSL2","FUT4","G0S2","GADD45A","GADD45B","GCH1","GEM","GFPT2","HBEGF","HES1",
        "ICAM1","ICOSLG","ID2","IER2","IER3","IER5","IFIH1","IFIT2","IFNGR2","IL12B","IL15RA",
        "IL18","IL1A","IL1B","IL23A","IL6","IL6ST","IL7R","INHBA","IRF1","IRS2","JAG1","JUN",
        "JUNB","KDM6B","KLF10","KLF2","KLF4","KLF6","KLF9","KYNU","LAMB3","LDLR","LIF",
        "LITAF","MAFF","MAP2K3","MAP3K8","MARCKS","MCL1","MSC","MXD1","MYC","NAMPT","NFAT5",
        "NFE2L2","NFIL3","NFKB1","NFKB2","NFKBIA","NFKBIE","NINJ1","NR4A1","NR4A2","NR4A3",
        "OLR1","PANX1","PDE4B","PDLIM7","PER1","PFKFB3","PHLDA1","PHLDA2","PLAU","PLAUR",
        "PLEK","PLK2","PLPP3","PMEPA1","PNRC1","PPAP2B","PPP1R15A","PTGER4","PTGS2","PTPRE",
        "PTX3","RCAN1","REL","RELA","RELB","RHOB","RIPK2","RNF19B","SAT1","SDC4","SERPINB2",
        "SERPINB8","SERPINE1","SGK1","SIK1","SLC16A6","SLC2A3","SLC2A6","SMAD3","SNN","SOCS3",
        "SOD2","SPHK1","SPSB1","SQSTM1","STAT5A","TANK","TAP1","TGIF1","TIPARP","TLR2",
        "TNC","TNF","TNFAIP2","TNFAIP3","TNFAIP6","TNFRSF9","TNFSF9","TNIP1","TNIP2","TRAF1",
        "TRAF3","TRIB1","TRIP10","TSC22D1","TUBB2A","VEGFA","YRDC","ZBTB10","ZC3H12A","ZFP36"
    },
}

DISSOLVE_TARGET_SETS = ["HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_FATTY_ACID_METABOLISM"]
EMERGE_TARGET_SETS = ["HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]


# ---------------------------------------------------------------------------
# Build weighted graph
# ---------------------------------------------------------------------------
def load_de_weights():
    log("Loading atlas for f2_inflection_logFC...")
    df = pd.read_csv(ATLAS, usecols=lambda c: c in ("human_symbol", "f2_inflection_logFC"))
    df = df.dropna(subset=["human_symbol"])
    df["f2_inflection_logFC"] = pd.to_numeric(df["f2_inflection_logFC"], errors="coerce").fillna(0.0)
    # Per-gene stage-weight vectors
    w_f01 = {}
    w_f34 = {}
    for sym, lfc in zip(df["human_symbol"], df["f2_inflection_logFC"]):
        if lfc < 0:
            w_f01[sym] = -lfc
        elif lfc > 0:
            w_f34[sym] = lfc
    log(f"  Genes with f2_inflection_logFC: {(df['f2_inflection_logFC'] != 0).sum():,}")
    log(f"  F0-F1 weighted (downreg at F3-F4): {len(w_f01):,}")
    log(f"  F3-F4 weighted (upreg at F3-F4): {len(w_f34):,}")
    return w_f01, w_f34


def build_graph(edges_df, w_gene):
    """Build an igraph from STRING edges with per-edge weights.
    edge_weight = string_score * exp( 0.5 * (w_gene[a] + w_gene[b]) )
    """
    genes = pd.unique(edges_df[["gene_a", "gene_b"]].values.ravel("K"))
    g = ig.Graph()
    g.add_vertices(list(genes))
    name_to_idx = {n: i for i, n in enumerate(g.vs["name"])}

    wa = np.array([w_gene.get(a, 0.0) for a in edges_df["gene_a"]])
    wb = np.array([w_gene.get(b, 0.0) for b in edges_df["gene_b"]])
    mult = np.exp(0.5 * (wa + wb))
    weights = edges_df["raw_score"].to_numpy() * mult

    src = [name_to_idx[a] for a in edges_df["gene_a"]]
    dst = [name_to_idx[b] for b in edges_df["gene_b"]]
    g.add_edges(list(zip(src, dst)))
    g.es["weight"] = weights.tolist()
    g.simplify(combine_edges={"weight": "max"})
    return g


def run_leiden(g, label):
    log(f"Running Leiden ({label}) on {g.vcount():,} nodes, {g.ecount():,} edges...")
    t = time.time()
    part = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=LEIDEN_RES,
        seed=SEED,
        n_iterations=-1,
    )
    mod = g.modularity(part.membership, weights="weight")
    log(f"  {label}: {len(part)} communities, Q={mod:.4f}, time={time.time()-t:.1f}s")
    return part.membership, mod


def relabel(membership):
    counter = Counter(membership)
    order = sorted(counter.keys(), key=lambda x: (-counter[x], x))
    m = {old: new for new, old in enumerate(order)}
    return [m[x] for x in membership]


# ---------------------------------------------------------------------------
# Enrichment (hypergeometric)
# ---------------------------------------------------------------------------
def enrich_communities(assign_df, universe, label):
    rows = []
    N = len(universe)
    for cid, sub in assign_df.groupby("community_id"):
        comm_genes = set(sub["gene"]) & universe
        if len(comm_genes) < 5:
            continue
        for set_name, gene_set in HALLMARK.items():
            K = len(gene_set & universe)
            n = len(comm_genes)
            k = len(comm_genes & gene_set)
            if k == 0 or K == 0:
                continue
            # P(X >= k) where X ~ Hypergeom(N, K, n)
            p = hypergeom.sf(k - 1, N, K, n)
            rows.append({
                "weighting": label,
                "community_id": cid,
                "community_size": len(sub),
                "gene_set": set_name,
                "overlap": k,
                "set_size": K,
                "comm_size_in_universe": n,
                "universe": N,
                "pvalue": p,
            })
    df = pd.DataFrame(rows)
    if len(df) > 0:
        df = df.sort_values(["community_id", "pvalue"]).reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# Transition classification
# ---------------------------------------------------------------------------
def build_transition(assign_f01, assign_f34):
    """Compute transition table (F01_id x F34_id overlap & Jaccard)."""
    a = assign_f01.rename(columns={"community_id": "F01_id"})
    b = assign_f34.rename(columns={"community_id": "F34_id"})
    merged = a.merge(b, on="gene", how="inner")

    f01_sizes = merged.groupby("F01_id")["gene"].size()
    f34_sizes = merged.groupby("F34_id")["gene"].size()

    pair = merged.groupby(["F01_id", "F34_id"]).size().reset_index(name="n_shared_genes")
    pair["n_F01"] = pair["F01_id"].map(f01_sizes)
    pair["n_F34"] = pair["F34_id"].map(f34_sizes)
    pair["jaccard"] = pair["n_shared_genes"] / (
        pair["n_F01"] + pair["n_F34"] - pair["n_shared_genes"]
    )

    # Classify each F01 community
    f01_class = {}
    for f01, grp in pair.groupby("F01_id"):
        total = grp["n_shared_genes"].sum()
        if total == 0:
            f01_class[f01] = "missing"
            continue
        max_frac = grp["n_shared_genes"].max() / total
        if f01_sizes.get(f01, 0) < MIN_COMM_SIZE:
            f01_class[f01] = "small"
        elif max_frac < DISSOLVE_THRESH:
            f01_class[f01] = "dissolving"
        elif max_frac >= 0.9:
            f01_class[f01] = "invariant"
        else:
            f01_class[f01] = "fragmenting"

    # Classify each F34 community: emerging if most members came from
    # small/absent F01 communities (i.e., weren't in any large F01 module).
    small_f01 = {c for c, s in f01_sizes.items() if s < MIN_COMM_SIZE}
    f34_class = {}
    for f34, grp in pair.groupby("F34_id"):
        total = grp["n_shared_genes"].sum()
        if total == 0:
            f34_class[f34] = "missing"
            continue
        scattered = grp[grp["F01_id"].isin(small_f01)]["n_shared_genes"].sum()
        max_frac = grp["n_shared_genes"].max() / total
        if f34_sizes.get(f34, 0) < MIN_COMM_SIZE:
            f34_class[f34] = "small"
        elif scattered / total >= EMERGE_THRESH or max_frac < DISSOLVE_THRESH:
            f34_class[f34] = "emerging"
        elif max_frac >= 0.9:
            f34_class[f34] = "invariant"
        else:
            f34_class[f34] = "fragmenting"

    pair["F01_class"] = pair["F01_id"].map(f01_class)
    pair["F34_class"] = pair["F34_id"].map(f34_class)
    # Top-level edge classification
    def _edge_class(r):
        if r["F01_class"] == "dissolving" or r["F34_class"] == "emerging":
            return r["F34_class"] if r["F34_class"] == "emerging" else "dissolving"
        if r["F01_class"] == "invariant" and r["F34_class"] == "invariant":
            return "invariant"
        return "fragmenting"
    pair["classification"] = pair.apply(_edge_class, axis=1)
    return pair, f01_class, f34_class


# ---------------------------------------------------------------------------
# Criterion 4 check
# ---------------------------------------------------------------------------
def check_criterion4(enr_f01, enr_f34, f01_class, f34_class):
    passed = True
    log("\n--- Criterion 4 check ---")

    # 4a: at least one emerging F3-F4 community enriched for inflammation/NF-kB
    emerging_ids = {cid for cid, c in f34_class.items() if c == "emerging"}
    hit_4a = None
    if len(enr_f34) > 0:
        cand = enr_f34[
            enr_f34["community_id"].isin(emerging_ids)
            & enr_f34["gene_set"].isin(EMERGE_TARGET_SETS)
            & (enr_f34["pvalue"] < ENRICH_P)
        ]
        if len(cand) > 0:
            hit_4a = cand.iloc[0]
    if hit_4a is not None:
        log(f"  4a PASS: F34 comm {int(hit_4a['community_id'])} enriched for "
            f"{hit_4a['gene_set']} (p={hit_4a['pvalue']:.2e})")
    else:
        log("  4a FAIL: no emerging F3-F4 community enriched for "
            "INFLAMMATORY_RESPONSE / TNFA_SIGNALING_VIA_NFKB at p<1e-5")
        passed = False

    # 4b: at least one dissolving F0-F1 community enriched for OXPHOS/FAO
    dissolving_ids = {cid for cid, c in f01_class.items() if c == "dissolving"}
    hit_4b = None
    if len(enr_f01) > 0:
        cand = enr_f01[
            enr_f01["community_id"].isin(dissolving_ids)
            & enr_f01["gene_set"].isin(DISSOLVE_TARGET_SETS)
            & (enr_f01["pvalue"] < ENRICH_P)
        ]
        if len(cand) > 0:
            hit_4b = cand.iloc[0]
    if hit_4b is not None:
        log(f"  4b PASS: F01 comm {int(hit_4b['community_id'])} enriched for "
            f"{hit_4b['gene_set']} (p={hit_4b['pvalue']:.2e})")
    else:
        log("  4b FAIL: no dissolving F0-F1 community enriched for "
            "OXIDATIVE_PHOSPHORYLATION / FATTY_ACID_METABOLISM at p<1e-5")
        passed = False

    log(f"\nCRITERION 4: {'PASS' if passed else 'FAIL'}")
    return passed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log("=== 292: F2-switch weighted Leiden ===")

    edges = pd.read_csv(NET_DIR / "edges_ppi.csv", usecols=["gene_a", "gene_b", "raw_score"])
    n0 = len(edges)
    edges = edges[edges["raw_score"] >= STRING_THRESH].reset_index(drop=True)
    log(f"STRING edges: {n0:,} total; {len(edges):,} at raw_score >= {STRING_THRESH}")
    if len(edges) == 0:
        log("ERROR: no edges after STRING threshold filter.")
        sys.exit(1)

    w_f01, w_f34 = load_de_weights()

    g_f01 = build_graph(edges, w_f01)
    mem_f01, mod_f01 = run_leiden(g_f01, "F0-F1")
    mem_f01 = relabel(mem_f01)

    g_f34 = build_graph(edges, w_f34)
    mem_f34, mod_f34 = run_leiden(g_f34, "F3-F4")
    mem_f34 = relabel(mem_f34)

    # Use F0-F1 graph gene list as canonical vertex order (identical to F3-F4)
    genes_f01 = g_f01.vs["name"]
    genes_f34 = g_f34.vs["name"]

    assign_f01 = pd.DataFrame({
        "gene": genes_f01,
        "community_id": mem_f01,
        "weighting": "F0-F1",
        "modularity": mod_f01,
    })
    assign_f34 = pd.DataFrame({
        "gene": genes_f34,
        "community_id": mem_f34,
        "weighting": "F3-F4",
        "modularity": mod_f34,
    })
    # Hierarchical columns (macro=meso=micro for single-resolution; filled for
    # schema consistency with downstream figure scripts).
    for df in (assign_f01, assign_f34):
        df["macro_id"] = df["community_id"]
        df["meso_id"] = df["community_id"]
        df["micro_id"] = df["community_id"]

    out_f01 = OUT_DIR / "communities_F01.csv"
    out_f34 = OUT_DIR / "communities_F34.csv"
    assign_f01.to_csv(out_f01, index=False)
    assign_f34.to_csv(out_f34, index=False)
    log(f"Wrote {out_f01} ({len(assign_f01):,} rows)")
    log(f"Wrote {out_f34} ({len(assign_f34):,} rows)")

    # Transition table
    trans, f01_class, f34_class = build_transition(assign_f01, assign_f34)
    trans_out = OUT_DIR / "community_transition_table.csv"
    trans.to_csv(trans_out, index=False)
    log(f"Wrote {trans_out} ({len(trans):,} rows)")

    # Enrichment
    universe = set(assign_f01["gene"]) | set(assign_f34["gene"])
    enr_f01 = enrich_communities(assign_f01, universe, "F0-F1")
    enr_f34 = enrich_communities(assign_f34, universe, "F3-F4")
    enr_f01_out = OUT_DIR / "community_enrichment_F01.csv"
    enr_f34_out = OUT_DIR / "community_enrichment_F34.csv"
    enr_f01.to_csv(enr_f01_out, index=False)
    enr_f34.to_csv(enr_f34_out, index=False)
    log(f"Wrote {enr_f01_out} ({len(enr_f01)} rows)")
    log(f"Wrote {enr_f34_out} ({len(enr_f34)} rows)")

    # Criterion 4
    check_criterion4(enr_f01, enr_f34, f01_class, f34_class)

    log(f"\n=== 292 Complete === total time {(time.time()-t0_global)/60:.1f} min")


if __name__ == "__main__":
    main()
