#!/usr/bin/env python
"""506b_liger_cnmf_hungarian.py

Hungarian (linear-sum-assignment) matching of LIGER k=16 factors to cNMF k=16
programs based on top-100 gene Jaccard. Closes adversarial-#4 Item 3 (LIGER vs
cNMF batch-program correspondence).

Inputs
------
- benchmarks/liger_factors_k16.tsv          (gene loadings, gene rows x 16 factor cols)
- benchmarks/liger_dataset_specificity_k16.tsv  (per-factor max dataset_specificity)
- benchmarks/program_jaccard_liger_vs_cnmf.tsv  (LIGER per-factor dataset_loading_flag)
- cnmf_annot/global/program_topgenes.k16.tsv   (cNMF top-N per program)
- reviewer_defense/dataset_confound.tsv         (cNMF per-program dataset_frac_var)

Output
------
- benchmarks/liger_cnmf_hungarian_match_k16.tsv  (16 matched pairs + Fisher 2x2)
"""
from __future__ import annotations

import os
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import fisher_exact

ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
TOP_N = 100
CNMF_BATCH_THRESH = 0.50  # dataset_frac_var threshold to flag cNMF program as batch-driven


def load_liger_top_genes(path: Path, top_n: int = TOP_N) -> dict[int, set[str]]:
    df = pd.read_csv(path, sep="\t", index_col=0)
    out: dict[int, set[str]] = {}
    for col in df.columns:
        # Higher loading = more characteristic; take top-N by descending value
        s = df[col].astype(float)
        top = s.sort_values(ascending=False).head(top_n).index.tolist()
        # Strip Ensembl base if present, keep symbol for Jaccard matching with cNMF
        out[int(col.replace("factor_", ""))] = set(top)
    return out


def load_cnmf_top_genes(path: Path, top_n: int = TOP_N) -> dict[int, set[str]]:
    df = pd.read_csv(path, sep="\t")
    df = df.sort_values(["program", "rank"])
    out: dict[int, set[str]] = {}
    for prog, grp in df.groupby("program"):
        out[int(prog)] = set(grp.head(top_n)["gene_name"].astype(str).tolist())
    return out


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def main() -> None:
    liger_factors = load_liger_top_genes(MCP / "benchmarks/liger_factors_k16.tsv")
    cnmf_progs = load_cnmf_top_genes(MCP / "cnmf_annot/global/program_topgenes.k16.tsv")

    liger_keys = sorted(liger_factors)
    cnmf_keys = sorted(cnmf_progs)
    assert len(liger_keys) == 16 and len(cnmf_keys) == 16, (
        f"expected 16x16, got {len(liger_keys)}x{len(cnmf_keys)}"
    )

    # 16x16 Jaccard matrix
    J = np.zeros((16, 16))
    for i, lf in enumerate(liger_keys):
        for j, cp in enumerate(cnmf_keys):
            J[i, j] = jaccard(liger_factors[lf], cnmf_progs[cp])

    # Hungarian on -J to MAXIMIZE total Jaccard
    row_ind, col_ind = linear_sum_assignment(-J)

    # LIGER per-factor metadata
    spec_df = pd.read_csv(MCP / "benchmarks/liger_dataset_specificity_k16.tsv", sep="\t")
    spec_df = spec_df[spec_df["k"] == 16].copy()
    # max specificity per factor across datasets
    spec_max = spec_df.groupby("factor")["dataset_specificity"].max()

    flag_df = pd.read_csv(MCP / "benchmarks/program_jaccard_liger_vs_cnmf.tsv", sep="\t")
    flag_df = flag_df[(flag_df["method"] == "LIGER(pyliger)") & (flag_df["k"] == 16)].copy()
    flag_lookup = dict(zip(flag_df["factor"], flag_df["dataset_loading_flag"].astype(str).str.lower() == "true"))

    # cNMF batch-artifact flag (dataset_frac_var > threshold)
    conf_df = pd.read_csv(MCP / "reviewer_defense/dataset_confound.tsv", sep="\t")
    # Map cnmf_global_k16_P{n} -> n
    conf_df["program"] = conf_df["program"].str.extract(r"P(\d+)$").astype(int)
    conf_lookup = dict(zip(conf_df["program"], conf_df["dataset_frac_var"]))

    rows = []
    for i, j in zip(row_ind, col_ind):
        lf = liger_keys[i]
        cp = cnmf_keys[j]
        liger_factor_str = f"factor_{lf}"
        liger_flag = bool(flag_lookup.get(liger_factor_str, False))
        cnmf_var = float(conf_lookup.get(cp, np.nan))
        cnmf_flag = bool(cnmf_var > CNMF_BATCH_THRESH)
        rows.append({
            "liger_factor": liger_factor_str,
            "cnmf_program": f"P{cp}",
            "jaccard_top100": round(float(J[i, j]), 4),
            "liger_dataset_specificity": round(float(spec_max.get(liger_factor_str, np.nan)), 4),
            "liger_dataset_loading_flag": liger_flag,
            "cnmf_dataset_variance_pct": round(cnmf_var * 100, 2),
            "cnmf_batch_artifact_flag": cnmf_flag,
        })

    out_df = pd.DataFrame(rows).sort_values("liger_factor", key=lambda s: s.str.replace("factor_", "").astype(int))

    # 2x2 contingency: liger-flagged x cnmf-flagged
    a = int(((out_df["liger_dataset_loading_flag"]) & (out_df["cnmf_batch_artifact_flag"])).sum())
    b = int(((out_df["liger_dataset_loading_flag"]) & (~out_df["cnmf_batch_artifact_flag"])).sum())
    c = int(((~out_df["liger_dataset_loading_flag"]) & (out_df["cnmf_batch_artifact_flag"])).sum())
    d = int(((~out_df["liger_dataset_loading_flag"]) & (~out_df["cnmf_batch_artifact_flag"])).sum())

    contingency = np.array([[a, b], [c, d]])
    odds_ratio, p_val = fisher_exact(contingency, alternative="greater")

    n_liger_flagged = a + b
    n_cnmf_flagged = a + c

    summary_lines = [
        f"# 506b_liger_cnmf_hungarian.py — k=16 matched pairs (Hungarian on top-{TOP_N} Jaccard)",
        f"# total Jaccard sum = {J[row_ind, col_ind].sum():.4f} (mean per pair = {J[row_ind, col_ind].mean():.4f})",
        f"# LIGER batch flag uses report-level `dataset_loading_flag` (max_dataset_specificity > 0.4 in upstream).",
        f"# cNMF batch flag uses dataset_frac_var > {CNMF_BATCH_THRESH:.2f} from reviewer_defense/dataset_confound.tsv",
        f"#   => cNMF-flagged programs at this threshold: " + ", ".join(
            sorted(out_df.loc[out_df['cnmf_batch_artifact_flag'], 'cnmf_program'].tolist())
        ),
        f"#   => LIGER-flagged factors: {n_liger_flagged}/16; cNMF-flagged programs: {n_cnmf_flagged}/16",
        f"# Contingency 2x2 (LIGER-flag x cNMF-flag): a={a} b={b} c={c} d={d}",
        f"# Fisher's exact (one-sided greater): OR={odds_ratio:.3g}, p={p_val:.3g}",
    ]
    if (odds_ratio > 5) and (p_val < 0.05):
        verdict = "STRENGTHEN: LIGER and cNMF jointly identify the same batch-driven programs."
    else:
        verdict = "STAND-AS-IS: count-level corroboration only; LIGER and cNMF disagree on WHICH programs are batch-driven."
    summary_lines.append(f"# verdict: {verdict}")

    out_path = MCP / "benchmarks/liger_cnmf_hungarian_match_k16.tsv"
    with open(out_path, "w") as fh:
        fh.write("\n".join(summary_lines) + "\n")
        out_df.to_csv(fh, sep="\t", index=False)

    print("\n".join(summary_lines))
    print(out_df.to_string(index=False))
    print(f"\nWrote: {out_path}")

    # Append summary to alpha REPORT.md
    alpha_report = MCP / "reviewer_audit/remediation/alpha/REPORT.md"
    if alpha_report.parent.exists():
        with open(alpha_report, "a") as fh:
            fh.write(textwrap.dedent(f"""

                ## 2026-04-25 Hungarian LIGER↔cNMF matching (adversarial-#4 Item 3)
                - 16x16 top-{TOP_N} Jaccard, total = {J[row_ind, col_ind].sum():.3f}, mean/pair = {J[row_ind, col_ind].mean():.3f}.
                - Contingency (LIGER-flag x cNMF-flag, dataset_frac_var>{CNMF_BATCH_THRESH:.2f}): a={a} b={b} c={c} d={d}; Fisher OR={odds_ratio:.3g}, p={p_val:.3g}.
                - {verdict}
                """).strip() + "\n")


if __name__ == "__main__":
    main()
