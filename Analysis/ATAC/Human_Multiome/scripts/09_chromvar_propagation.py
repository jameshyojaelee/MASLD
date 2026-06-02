#!/usr/bin/env python3
"""
09_chromvar_propagation.py — Propagate chromVAR per-donor TF activity into
downstream tests (B5 of the ATAC improvement plan).

Inputs
------
  - Analysis/ATAC/Human_Multiome/results/chromvar_v2/chromvar_per_donor.tsv.gz
        (output of 03b_chromvar_per_donor_export.py)
  - Analysis/ATAC/Human_Multiome/metadata/donor_metadata_curated.tsv
        (B1 bridge MM_ <-> JB_ <-> D## donor identifiers + condition)
  - Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv
        (24 SCENIC+ disease regulons; regulon_activity_diff sign for concordance)
  - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
        (progression_gene_class + sex_class for fgsea sets)
  - RNA-seq/results/drug_repurposing/lincs_reversal_compounds.csv
        (drug target gene symbols for the LINCS-target fgsea)

Outputs (Analysis/ATAC/Human_Multiome/results/chromvar_v2/)
-------
  - chromvar_disease_wilcoxon_per_donor.csv  (TF x cell_type, disease vs healthy
        Mann-Whitney across donors, BH-adjusted within cell type)
  - chromvar_scenic_concordance.csv          (chromVAR direction vs SCENIC+
        regulon_activity_diff sign, for the 24 disease regulons in each CT)
  - chromvar_lincs_fgsea.csv                 (TF symbols ranked by disease
        Wilcoxon Z; pre-ranked GSEA vs LINCS reversal drug target genes;
        skipped with a marker file if LINCS set is empty)
  - chromvar_progression_fgsea.csv           (vs atlas progression_gene_class
        gene sets: onset_only / fibrosis_specific / inflammation_specific /
        ubiquitous / progression_only / both_onset_and_progression / extreme_only)
  - chromvar_sex_fgsea.csv                   (vs atlas sex_class:
        Female_biased / Male_biased / Divergent / Concordant)

Notes
-----
* TF -> gene-symbol mapping. JASPAR motif IDs include monomers (HNF4A),
  rodent paralogs lowercased (Hnf4a -> mouse), and dimers (RXRA::VDR,
  ARNT::HIF1A). We split on '::' and upper-case to get human symbols
  for cross-method joins. For fgsea ranking we deduplicate per gene by
  taking the strongest-effect TF row.
* Disease grouping for the per-donor Wilcoxon: disease := MASL + MASH
  (n=12), healthy := NORMAL (n=6). At the lower bound this is n=2 vs
  n>=2 — for that we still report a p-value but the result is informational.
* The pre-ranked GSEA is implemented in-house (running-sum statistic +
  permutation null) because gseapy is not in the snapatac2 env. Standard
  Subramanian 2005 formula, n_perm=1000 by default.

Environment: snapatac2 (with PYTHONNOUSERSITE=1).
"""
from __future__ import annotations

import logging
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATAC_DIR = PROJECT_ROOT / "Analysis" / "ATAC" / "Human_Multiome"
OUT_DIR = ATAC_DIR / "results" / "chromvar_v2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PER_DONOR_TSV = OUT_DIR / "chromvar_per_donor.tsv.gz"
DONOR_META_TSV = ATAC_DIR / "metadata" / "donor_metadata_curated.tsv"
SCENIC_CSV = ATAC_DIR / "scenic_plus" / "disease_regulons.csv"

ATLAS_CSV = PROJECT_ROOT / "RNA-seq" / "results" / "multi_evidence" / "multi_evidence_atlas.csv"
LINCS_CSV = PROJECT_ROOT / "RNA-seq" / "results" / "drug_repurposing" / "lincs_reversal_compounds.csv"

OUT_WILCOX = OUT_DIR / "chromvar_disease_wilcoxon_per_donor.csv"
OUT_CONCORD = OUT_DIR / "chromvar_scenic_concordance.csv"
OUT_FGSEA_LINCS = OUT_DIR / "chromvar_lincs_fgsea.csv"
OUT_FGSEA_PROG = OUT_DIR / "chromvar_progression_fgsea.csv"
OUT_FGSEA_SEX = OUT_DIR / "chromvar_sex_fgsea.csv"

MIN_DONORS_PER_GROUP = 2  # disease/healthy floor for Wilcoxon
MIN_DONORS_TOTAL = 5      # require >=5 donors with non-NaN values for the test
RNG_SEED = 42
N_PERM_FGSEA = 1000
MIN_SET_SIZE = 5
MAX_SET_SIZE = 500


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def split_dimer(tf: str) -> List[str]:
    """JASPAR motif id -> human gene symbols.

    Splits dimer IDs (e.g., 'RXRA::VDR'), uppercases to match human symbols,
    and strips suffix tags added during deduplication ('ATF3.1' -> 'ATF3').
    """
    parts: List[str] = []
    for piece in str(tf).split("::"):
        sym = piece.split(".")[0].upper()
        if sym:
            parts.append(sym)
    return parts


def collapse_to_gene_rank(df: pd.DataFrame, score_col: str, gene_col: str = "TF_symbol_h") -> pd.Series:
    """One score per gene by abs-max across all motif rows mapping to that gene.

    Returns a Series indexed by gene symbol, sorted descending by signed score.
    """
    # signed-value-of-max-abs across motif rows
    by_gene = df.groupby(gene_col)[score_col].apply(
        lambda s: s.iloc[s.abs().argmax()] if s.notna().any() else np.nan
    )
    by_gene = by_gene.dropna().sort_values(ascending=False)
    return by_gene


def fgsea_preranked(
    ranks: pd.Series,
    gene_sets: Dict[str, set],
    n_perm: int = N_PERM_FGSEA,
    seed: int = RNG_SEED,
    min_size: int = MIN_SET_SIZE,
    max_size: int = MAX_SET_SIZE,
) -> pd.DataFrame:
    """Pre-ranked GSEA per Subramanian 2005 with permutation null.

    Args:
        ranks: pd.Series with gene symbol index, signed scores, sorted descending.
        gene_sets: dict[set_name] -> set of gene symbols.
        n_perm: permutation iterations for null ES distribution.

    Returns:
        DataFrame: gene_set | size | ES | NES | pvalue | padj | leading_edge_size
    """
    if ranks.empty:
        return pd.DataFrame()

    rng = np.random.default_rng(seed)
    n_genes = len(ranks)
    gene_index = {g: i for i, g in enumerate(ranks.index)}
    scores = ranks.values.astype(np.float64)
    abs_scores = np.abs(scores)
    abs_sum_total = abs_scores.sum()
    if abs_sum_total == 0:
        return pd.DataFrame()

    def _es(member_mask: np.ndarray) -> Tuple[float, int]:
        """Running-sum enrichment + leading-edge size (members at/before peak)."""
        hits = member_mask
        misses = ~member_mask
        n_hits = int(hits.sum())
        n_miss = int(misses.sum())
        if n_hits == 0 or n_miss == 0:
            return 0.0, 0
        hit_w = np.where(hits, abs_scores, 0.0)
        norm = hit_w.sum()
        if norm == 0:
            return 0.0, 0
        p_hit = hit_w / norm
        p_miss = misses.astype(np.float64) / n_miss
        running = np.cumsum(p_hit - p_miss)
        max_pos = float(running.max())
        max_neg = float(running.min())
        es = max_pos if abs(max_pos) >= abs(max_neg) else max_neg
        # leading edge: members up to (and including) the index of |running| max
        peak_idx = int(np.argmax(np.abs(running)))
        if es >= 0:
            leading = int(np.sum(hits[: peak_idx + 1]))
        else:
            leading = int(np.sum(hits[peak_idx:]))
        return es, leading

    rows = []
    for set_name, members in gene_sets.items():
        # observed
        in_universe = [g for g in members if g in gene_index]
        m = len(in_universe)
        if m < min_size or m > max_size:
            continue
        mask = np.zeros(n_genes, dtype=bool)
        mask[[gene_index[g] for g in in_universe]] = True
        es_obs, le_size = _es(mask)

        # permutation null (gene-set permutation: random m positions)
        null_es = np.empty(n_perm, dtype=np.float64)
        for i in range(n_perm):
            perm = rng.choice(n_genes, size=m, replace=False)
            mask_p = np.zeros(n_genes, dtype=bool)
            mask_p[perm] = True
            null_es[i], _ = _es(mask_p)

        # p-value: directional one-sided null counts
        if es_obs >= 0:
            n_ge = int((null_es >= es_obs).sum()) + 1
            p = n_ge / (n_perm + 1)
            pos = null_es[null_es > 0]
            mean_pos = pos.mean() if pos.size else 1e-9
            nes = es_obs / max(abs(mean_pos), 1e-9)
        else:
            n_le = int((null_es <= es_obs).sum()) + 1
            p = n_le / (n_perm + 1)
            neg = null_es[null_es < 0]
            mean_neg = neg.mean() if neg.size else -1e-9
            nes = es_obs / max(abs(mean_neg), 1e-9)
        rows.append((set_name, m, es_obs, nes, p, le_size))

    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows, columns=["gene_set", "size", "ES", "NES", "pvalue", "leading_edge_size"])
    _, padj, _, _ = multipletests(out["pvalue"], method="fdr_bh")
    out["padj"] = padj
    return out.sort_values("padj").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------
def main() -> None:
    t0 = time.time()
    log.info("=" * 72)
    log.info("chromVAR propagation (B5)")
    log.info("=" * 72)

    # -- Load per-donor TF deviations ----------------------------------------
    log.info(f"Reading {PER_DONOR_TSV}")
    per_donor = pd.read_csv(PER_DONOR_TSV, sep="\t")
    log.info(f"  rows: {len(per_donor):,}; donors: {per_donor['donor_id_atac'].nunique()}; "
             f"cell types: {per_donor['cell_type'].nunique()}; TFs: {per_donor['TF'].nunique()}")

    # -- Donor metadata -------------------------------------------------------
    log.info(f"Reading {DONOR_META_TSV}")
    meta = pd.read_csv(DONOR_META_TSV, sep="\t")
    log.info(f"  donors: {meta['donor_id'].nunique()}; conditions: {meta['condition'].value_counts().to_dict()}")
    # The per-donor TSV uses donor_id_atac = D## (matches h5ad obs.donor_id);
    # donor_metadata_curated has both donor_id (D##) and donor_id_atac (MM_##).
    # Join on D##.
    meta_keep = meta[["donor_id", "donor_id_atac", "donor_id_rna", "condition", "disease_stage_coarse"]].rename(
        columns={"donor_id": "donor_id_atac_local", "donor_id_atac": "donor_id_mm"}
    )
    joined = per_donor.merge(
        meta_keep, how="left", left_on="donor_id_atac", right_on="donor_id_atac_local"
    ).drop(columns=["donor_id_atac_local"])
    n_miss = joined["condition"].isna().sum()
    if n_miss:
        log.warning(f"  {n_miss:,} per-donor rows have no condition match")
    else:
        log.info("  all per-donor rows matched to a donor in the curated bridge")

    # Boolean disease label
    joined["is_disease"] = joined["condition"].isin(["MASL", "MASH"])

    # ------------------------------------------------------------------------
    # Step 1: Wilcoxon disease vs healthy across donors, per (TF, cell type)
    # ------------------------------------------------------------------------
    log.info("-" * 60)
    log.info("Step 1: per-donor Mann-Whitney disease vs healthy")
    wilcox_records = []
    for (tf, ct), grp in joined.groupby(["TF", "cell_type"], observed=True):
        # Drop rows missing mean_deviation (rare; cell-type absent in a donor)
        valid = grp.dropna(subset=["mean_deviation", "is_disease"])
        if valid["donor_id_atac"].nunique() < MIN_DONORS_TOTAL:
            continue
        d_vals = valid.loc[valid["is_disease"], "mean_deviation"].values
        h_vals = valid.loc[~valid["is_disease"], "mean_deviation"].values
        if len(d_vals) < MIN_DONORS_PER_GROUP or len(h_vals) < MIN_DONORS_PER_GROUP:
            continue
        try:
            stat, pval = mannwhitneyu(d_vals, h_vals, alternative="two-sided")
        except ValueError:
            pval = 1.0
            stat = np.nan
        mean_d = float(np.mean(d_vals))
        mean_h = float(np.mean(h_vals))
        logFC = mean_d - mean_h  # already z-scored deviations
        wilcox_records.append(
            {
                "TF": tf,
                "TF_symbol": grp["TF_symbol"].iloc[0],
                "cell_type": ct,
                "n_disease": len(d_vals),
                "n_healthy": len(h_vals),
                "mean_disease": round(mean_d, 6),
                "mean_healthy": round(mean_h, 6),
                "logFC": round(logFC, 6),
                "U_stat": float(stat) if not np.isnan(stat) else np.nan,
                "pvalue": pval,
            }
        )
    wilcox_df = pd.DataFrame(wilcox_records)

    # BH within cell type
    if not wilcox_df.empty:
        wilcox_df["padj"] = np.nan
        for ct in wilcox_df["cell_type"].unique():
            mask = wilcox_df["cell_type"] == ct
            if mask.sum() > 0:
                _, padj, _, _ = multipletests(wilcox_df.loc[mask, "pvalue"], method="fdr_bh")
                wilcox_df.loc[mask, "padj"] = padj
        wilcox_df = wilcox_df.sort_values(["cell_type", "padj"]).reset_index(drop=True)
        log.info(f"  Wilcoxon table: {len(wilcox_df):,} rows; "
                 f"padj<0.05: {(wilcox_df['padj'] < 0.05).sum():,}")
        wilcox_df.to_csv(OUT_WILCOX, index=False)
        log.info(f"  Wrote {OUT_WILCOX}")
    else:
        log.warning("  Wilcoxon table empty; skipping output")

    # ------------------------------------------------------------------------
    # Step 2: chromVAR vs SCENIC+ regulon direction concordance
    # ------------------------------------------------------------------------
    log.info("-" * 60)
    log.info("Step 2: chromVAR vs SCENIC+ regulon direction concordance")
    scenic = pd.read_csv(SCENIC_CSV)
    log.info(f"  SCENIC+ disease regulons: {len(scenic):,}")
    # scenic.regulon_activity_diff signed: positive = up in MASLD

    # For each scenic TF, take the TF symbol (e.g., HNF4A) and find chromVAR
    # rows whose TF_symbol matches (case-insensitive) — keep all dimers that
    # include the TF as a partner so that the concordance test reflects
    # all dimer contexts.
    wilcox_df["TF_partners"] = wilcox_df["TF"].apply(split_dimer)
    concord_records = []
    if not wilcox_df.empty:
        for _, srow in scenic.iterrows():
            tf_q = str(srow["tf_name"]).upper()
            scenic_dir = float(srow["regulon_activity_diff"])
            for ct in wilcox_df["cell_type"].unique():
                ct_rows = wilcox_df[
                    (wilcox_df["cell_type"] == ct)
                    & wilcox_df["TF_partners"].apply(lambda x: tf_q in x)
                ]
                if ct_rows.empty:
                    continue
                # use the largest |logFC| match in that cell type
                best = ct_rows.iloc[ct_rows["logFC"].abs().argmax()]
                chrom_dir = float(best["logFC"])
                agree = int(np.sign(chrom_dir) == np.sign(scenic_dir)) if (chrom_dir != 0 and scenic_dir != 0) else 0
                concord_records.append(
                    {
                        "regulon_id": srow["regulon_id"],
                        "tf_name": srow["tf_name"],
                        "cell_type": ct,
                        "chromvar_TF_motif": best["TF"],
                        "chromvar_logFC": round(chrom_dir, 6),
                        "chromvar_padj": float(best["padj"]) if pd.notna(best["padj"]) else np.nan,
                        "scenic_regulon_diff": round(scenic_dir, 6),
                        "scenic_padj": float(srow["activity_padj"]),
                        "direction_agree": agree,
                    }
                )
        concord_df = pd.DataFrame(concord_records)
        if not concord_df.empty:
            concord_df = concord_df.sort_values(["regulon_id", "cell_type"]).reset_index(drop=True)
            concord_df.to_csv(OUT_CONCORD, index=False)
            log.info(f"  Concordance rows: {len(concord_df):,}")
            agree_rate = concord_df["direction_agree"].mean()
            log.info(f"  Direction agreement (all 24 regulons x all cell types where chromVAR present): {agree_rate * 100:.1f}%")
            # Per-regulon summary (best across cell types — i.e., is there *any* CT
            # where chromVAR agrees with SCENIC+?)
            per_reg = concord_df.groupby(["regulon_id", "tf_name"])["direction_agree"].agg(["max", "mean", "count"]).reset_index()
            n_any_agree = (per_reg["max"] == 1).sum()
            log.info(f"  {n_any_agree}/{len(per_reg)} regulons agree in at least one CT")
            log.info(f"  Wrote {OUT_CONCORD}")
        else:
            log.warning("  No concordance rows (no Wilcoxon overlap with SCENIC+ TFs)")
    else:
        log.warning("  Skipping concordance — no Wilcoxon rows")

    # ------------------------------------------------------------------------
    # Step 3: pre-ranked GSEA across cell types
    # ------------------------------------------------------------------------
    log.info("-" * 60)
    log.info("Step 3: pre-ranked GSEA (chromVAR Wilcoxon Z by cell type)")

    # Build per-cell-type signed-z ranking of TFs (use Z = sign(logFC) * -log10(pvalue)
    # to combine magnitude + direction; avoids the all-near-zero problem of raw logFC)
    def signed_z(row) -> float:
        p = row["pvalue"]
        s = row["logFC"]
        if pd.isna(p) or p <= 0:
            return 0.0
        return np.sign(s) * -np.log10(max(p, 1e-300))

    if not wilcox_df.empty:
        wilcox_df["signed_z"] = wilcox_df.apply(signed_z, axis=1)
        # TF -> primary gene symbol for ranking (use first partner only;
        # most TFs in chromVAR are monomers anyway)
        wilcox_df["TF_symbol_h"] = wilcox_df["TF"].apply(
            lambda x: split_dimer(x)[0] if split_dimer(x) else x.upper()
        )

    # Load atlas gene set definitions
    log.info("  Loading atlas progression / sex class columns")
    atlas = pd.read_csv(
        ATLAS_CSV,
        usecols=["human_symbol", "sex_class", "progression_gene_class"],
        low_memory=False,
    )
    atlas["human_symbol"] = atlas["human_symbol"].astype(str).str.upper()

    progression_sets: Dict[str, set] = {}
    for cls, grp in atlas.dropna(subset=["progression_gene_class"]).groupby("progression_gene_class"):
        if cls in {"not_significant"}:
            continue
        progression_sets[f"prog::{cls}"] = set(grp["human_symbol"])
    for k, v in progression_sets.items():
        log.info(f"    {k}: {len(v):,} genes")

    sex_sets: Dict[str, set] = {}
    for cls, grp in atlas.dropna(subset=["sex_class"]).groupby("sex_class"):
        if cls == "Uncertain" or cls == "Concordant":
            # Concordant is the null background; including it is uninformative
            # (would be the whole universe). Uncertain similarly background.
            continue
        sex_sets[f"sex::{cls}"] = set(grp["human_symbol"])
    for k, v in sex_sets.items():
        log.info(f"    {k}: {len(v):,} genes")

    # LINCS drug target set
    lincs_sets: Dict[str, set] = {}
    if LINCS_CSV.exists():
        lincs = pd.read_csv(LINCS_CSV, low_memory=False)
        # split multi-target strings
        def explode(s):
            if pd.isna(s):
                return []
            return [x.strip().upper() for x in str(s).replace("|", ";").replace(",", ";").split(";") if x.strip()]
        all_targets = lincs["target"].apply(explode).explode().dropna()
        all_targets = all_targets[all_targets.str.match(r"^[A-Z0-9_-]+$")]
        if all_targets.size:
            lincs_sets["lincs::reversal_drug_targets"] = set(all_targets.unique())
            # Also a stratified set by MOA when MOA non-empty
            if "moa" in lincs.columns:
                moa_grp = lincs.dropna(subset=["target", "moa"])
                moa_grp = moa_grp[moa_grp["moa"].astype(str).str.strip().str.strip('"').str.len() > 0]
                for moa, g in moa_grp.groupby("moa"):
                    if str(moa).strip() in ("", '"', '""', '""""'):
                        continue
                    targets = g["target"].apply(explode).explode().dropna()
                    targets = targets[targets.str.match(r"^[A-Z0-9_-]+$")]
                    if targets.nunique() >= MIN_SET_SIZE:
                        lincs_sets[f"lincs_moa::{moa}"] = set(targets.unique())
            log.info(f"    LINCS reversal-target sets prepared: {len(lincs_sets)}")
            log.info(f"      universe: {len(lincs_sets.get('lincs::reversal_drug_targets', set())):,} unique targets")
        else:
            log.warning("    LINCS reversal compound 'target' column empty after parsing; skipping LINCS fgsea")
    else:
        log.warning(f"    {LINCS_CSV} missing; skipping LINCS fgsea")

    def run_fgsea_per_celltype(sets: Dict[str, set], label: str) -> pd.DataFrame:
        if not sets:
            return pd.DataFrame()
        cts = wilcox_df["cell_type"].unique()
        rows: List[pd.DataFrame] = []
        for ct in cts:
            sub = wilcox_df[wilcox_df["cell_type"] == ct].copy()
            if sub.empty:
                continue
            ranks = collapse_to_gene_rank(sub, score_col="signed_z", gene_col="TF_symbol_h")
            if ranks.size < MIN_SET_SIZE:
                continue
            r = fgsea_preranked(ranks, sets, n_perm=N_PERM_FGSEA, seed=RNG_SEED)
            if not r.empty:
                r.insert(0, "cell_type", ct)
                rows.append(r)
        if not rows:
            return pd.DataFrame()
        out = pd.concat(rows, ignore_index=True)
        log.info(f"  {label}: {len(out):,} (cell_type, gene_set) tests; significant padj<0.1: {(out['padj'] < 0.1).sum()}")
        return out

    prog_fgsea = run_fgsea_per_celltype(progression_sets, "progression fgsea")
    if not prog_fgsea.empty:
        prog_fgsea.to_csv(OUT_FGSEA_PROG, index=False)
        log.info(f"  Wrote {OUT_FGSEA_PROG}")
    sex_fgsea = run_fgsea_per_celltype(sex_sets, "sex fgsea")
    if not sex_fgsea.empty:
        sex_fgsea.to_csv(OUT_FGSEA_SEX, index=False)
        log.info(f"  Wrote {OUT_FGSEA_SEX}")
    lincs_fgsea = run_fgsea_per_celltype(lincs_sets, "LINCS fgsea")
    if not lincs_fgsea.empty:
        lincs_fgsea.to_csv(OUT_FGSEA_LINCS, index=False)
        log.info(f"  Wrote {OUT_FGSEA_LINCS}")
    elif not lincs_sets:
        # write marker so downstream knows it was skipped
        marker = OUT_FGSEA_LINCS.with_suffix(".skipped")
        marker.write_text("LINCS gene set was empty or LINCS file missing — see log\n")
        log.info(f"  Wrote skip marker {marker}")

    # ------------------------------------------------------------------------
    # Summaries
    # ------------------------------------------------------------------------
    log.info("=" * 72)
    log.info("Summaries")
    log.info("=" * 72)
    if not wilcox_df.empty:
        log.info("Top 5 TFs by disease Wilcoxon padj per cell type:")
        for ct, sub in wilcox_df.groupby("cell_type"):
            top = sub.sort_values("padj").head(5)
            for _, r in top.iterrows():
                log.info(
                    f"  {ct:>22s}  {r['TF']:>22s}  logFC={r['logFC']:+.3f}  padj={r['padj']:.2e}  "
                    f"n_dis/heal={int(r['n_disease'])}/{int(r['n_healthy'])}"
                )

    log.info(f"Total runtime: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
