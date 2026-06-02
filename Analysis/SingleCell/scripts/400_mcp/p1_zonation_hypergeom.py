#!/usr/bin/env python
"""
p1_zonation_hypergeom.py — Charlie Task 14 (Action A4, blocker B4).

Test whether cNMF P1 (k=16) top-100 genes are significantly enriched for
Halpern 2017 (Nature 542:352) periportal / pericentral marker panels, with
universe = cNMF global HVG list (N=2499).

If periportal OR >= 3 AND p < 0.001 -> relabel P1 from "hepatocyte metabolic
loss" to "mature-periportal hepatocyte metabolic program" in RESULTS_FINAL_k16.md.

Outputs: reviewer_defense/p1_zonation_enrichment.tsv
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
OUT = MCP / "reviewer_defense"
OUT.mkdir(exist_ok=True, parents=True)

K_HEAD = 16

# Halpern 2017 Nature 542:352 Supplementary Table S3 canonical zonated markers
# Supplemented with Ben-Moshe 2019 Nat Rev Gastro Hepatol review panels.
# Source: PMID 28166538 + PMID 31597341.
PERIPORTAL = [
    "GHR", "ADH4", "FMO3", "ACSM2A", "ACSM2B", "AGMO", "BHMT", "GLYAT", "AFM",
    "A1CF", "ALB", "CYP3A4", "ARG1", "HAL", "ASS1", "SDS", "SERPINA1",
    "PCK1", "G6PC", "TAT", "ASL", "CPS1", "OTC", "HNF4A", "NNMT", "FBP1",
    "SLC27A5", "APOF", "CYP2A6", "CYP2A7", "HSD11B1", "APOM", "TDO2",
    "CYP2C8", "CYP2D6", "HAO1", "CYP4A11", "AADAT", "CYP1A2",
    "IGFBP3", "MGLL", "UPB1", "SULT2A1", "MPST", "KHK", "GNMT",
    "LDHD", "BAAT", "PON3", "ABCA1",
]

PERICENTRAL = [
    "CYP2E1", "GLUL", "CYP7A1", "CYP1A2", "OAT", "BCHE", "AXIN2", "LGR5",
    "RSPO3", "NTRK2", "SLCO1B3", "ADH1A", "ADH1B", "ADH1C", "CYP3A5",
    "CYP2C9", "CYP2C19", "CYP2B6", "CYP3A4", "AKR1C4", "AKR1D1",
    "GSTA1", "GSTA2", "GSTA3", "WNT2", "NOTUM", "TBX3", "SLC1A2",
    "SPP1", "CCL2", "GPR65", "SLC7A2", "HSD17B13", "GDF15",
    "GULP1", "CES1", "CYP2C18", "CYP2J2", "RHBG",
    "ANGPTL6", "CFI", "SULF2", "TAF4B", "ARG2", "MT1X",
]


def hypergeom_test(observed: int, n_sample: int, n_success: int, N: int):
    """Fisher-style hypergeometric upper-tail p-value + OR + expected.
    observed: overlap of sample with panel
    n_sample: |cNMF top-100 in universe|
    n_success: |panel in universe|
    N: |universe|
    """
    if N == 0 or n_sample == 0 or n_success == 0:
        return {"observed": observed, "expected": 0.0, "OR": np.nan, "p": 1.0}
    expected = n_sample * n_success / N
    # Upper tail: P(X >= observed)
    p = 1.0 - stats.hypergeom.cdf(observed - 1, N, n_success, n_sample)
    # Odds ratio: [a / (n_sample-a)] / [(n_success-a) / (N-n_success-(n_sample-a))]
    a = observed; b = n_sample - a; c = n_success - a; d = N - n_success - b
    if b == 0 or c == 0 or d == 0:
        # Avoid div-by-zero: add 0.5 continuity correction
        a2, b2, c2, d2 = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    else:
        a2, b2, c2, d2 = a, b, c, d
    OR = (a2 * d2) / (b2 * c2) if (b2 * c2) > 0 else np.nan
    return {"observed": int(observed), "expected": float(expected),
            "OR": float(OR), "p": float(p)}


def bh_q(ps):
    ps = np.asarray(ps, float)
    n = ps.size
    order = np.argsort(ps)
    ranked = np.empty(n, float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        j = order[i]
        prev = min(prev, ps[j] * n / (i + 1))
        ranked[j] = prev
    return ranked


def main() -> None:
    # Universe: cNMF global HVG
    hvg = set(pd.read_csv(MCP / "cnmf_runs/global/global.overdispersed_genes.txt",
                          header=None)[0].astype(str).tolist())
    N = len(hvg)
    print(f"[P1-zon] universe (HVG) = {N}")

    # P1 top-100
    topg = pd.read_csv(MCP / f"cnmf_annot/global/program_topgenes.k{K_HEAD}.tsv", sep="\t")
    p1 = topg.loc[topg["program"] == 1].nlargest(100, "spectra_score")["gene_name"].astype(str).tolist()
    p1_in_univ = [g for g in p1 if g in hvg]
    n_sample = len(p1_in_univ)
    print(f"[P1-zon] P1 top-100 in universe: {n_sample}/100")

    rows = []
    for panel_name, panel in [("periportal_Halpern2017", PERIPORTAL),
                              ("pericentral_Halpern2017", PERICENTRAL)]:
        panel_in_univ = set(panel) & hvg
        panel_in_univ_n = len(panel_in_univ)
        overlap = len(set(p1_in_univ) & panel_in_univ)
        overlap_genes = sorted(set(p1_in_univ) & panel_in_univ)
        stats_d = hypergeom_test(overlap, n_sample, panel_in_univ_n, N)
        rows.append({
            "panel": panel_name,
            "panel_size": len(panel),
            "panel_in_universe": panel_in_univ_n,
            "p1_top100_in_universe": n_sample,
            "observed_overlap": stats_d["observed"],
            "expected_overlap": round(stats_d["expected"], 3),
            "OR": round(stats_d["OR"], 3) if stats_d["OR"] == stats_d["OR"] else np.nan,
            "p": stats_d["p"],
            "overlap_genes": ",".join(overlap_genes),
        })
    df = pd.DataFrame(rows)
    df["q"] = bh_q(df["p"].values)
    df.to_csv(OUT / "p1_zonation_enrichment.tsv", sep="\t", index=False)
    print(df.to_string(index=False))

    # Decision flag
    pp = df.loc[df["panel"] == "periportal_Halpern2017"].iloc[0]
    pc = df.loc[df["panel"] == "pericentral_Halpern2017"].iloc[0]
    relabel = (pp["OR"] >= 3.0) and (pp["p"] < 0.001)
    print(f"\n[P1-zon] periportal OR={pp['OR']:.2f}, p={pp['p']:.2e} → relabel={relabel}")
    print(f"[P1-zon] pericentral OR={pc['OR']:.2f}, p={pc['p']:.2e}")
    # Write a flag file for the downstream text-edit step
    with open(OUT / "p1_relabel_flag.txt", "w") as fh:
        fh.write("RELABEL\n" if relabel else "KEEP\n")
        fh.write(f"periportal_OR\t{pp['OR']}\nperipportal_p\t{pp['p']}\n")
        fh.write(f"pericentral_OR\t{pc['OR']}\npericentral_p\t{pc['p']}\n")


if __name__ == "__main__":
    main()
