#!/usr/bin/env python3
"""
Proteo-transcriptomic integration: population-level fold-change concordance
between protein and transcript abundances.

Concordance is computed at the population (per-gene fold-change) level, not by
positional per-sample pairing: the plasma proteomics (PXD052937) and the scRNA
hepatocyte pseudobulk come from independent cohorts with no shared subjects, so
correlating samples by array position is meaningless. We instead correlate the
disease-vs-control protein logFC against the dream mega-analysis transcript
logFC across the shared gene set (the pattern used by differential_proteomics.R
and mrna_protein_concordance.R).

Inputs:
  - Hepatocyte pseudobulk from scVI
    (Analysis/SingleCell/results_gpu_v2/pseudobulk/Hepatocytes_pseudobulk.csv;
    HGNC row index, SRR-accession columns)
  - Protein matrix (Analysis/Proteomics/results/pxd052937_protein_matrix.csv;
    UniProt-accession row index, remapped to HGNC via Candidates.tsv)
  - PXD052937 disease metadata (group labels for protein logFC)
  - dream mega-analysis transcript logFC (transcript side of the concordance)
  - Multi-evidence scores (RNA-seq/results/multi_evidence/multi_evidence_scored_genes.csv)
  - TWAS/COLOC-convergent drug targets
    (RNA-seq/results/drug_repurposing/convergent_drug_targets.csv)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def load_uniprot_to_hgnc(candidates_file: Path) -> dict:
    """Build a UniProt-accession -> HGNC-symbol map from the Spectronaut
    Candidates.tsv (same pattern as differential_proteomics.R ~lines 213-225):
    protein_id = first ';'-token of ProteinGroups -> gene = first token of Genes,
    deduped by protein_id.
    """
    if not candidates_file.exists():
        raise FileNotFoundError(
            f"Candidates.tsv not found (needed for UniProt->HGNC remap): {candidates_file}")

    cand = pd.read_csv(candidates_file, sep="\t", usecols=["ProteinGroups", "Genes"])
    cand = cand.dropna(subset=["ProteinGroups"]).copy()
    cand["protein_id"] = cand["ProteinGroups"].astype(str).str.replace(r";.*", "", regex=True)
    cand["gene"] = cand["Genes"].astype(str).str.replace(r";.*", "", regex=True)
    cand = cand[(cand["gene"] != "") & (cand["gene"].str.lower() != "nan")]
    cand = cand.drop_duplicates(subset="protein_id")
    return dict(zip(cand["protein_id"], cand["gene"]))


def remap_protein_index_to_hgnc(protein_df: pd.DataFrame, uniprot_map: dict) -> pd.DataFrame:
    """Remap a UniProt-accession-indexed protein matrix to HGNC symbols BEFORE
    intersecting with HGNC-keyed transcript data. Rows without a mapping are
    dropped; collisions (multiple accessions -> same gene) are collapsed by mean.
    """
    mapped = protein_df.copy()
    mapped.index = mapped.index.map(lambda pid: uniprot_map.get(pid, pid))
    mapped = mapped.loc[mapped.index.isin(set(uniprot_map.values()))]
    if mapped.index.duplicated().any():
        mapped = mapped.groupby(level=0).mean()
    return mapped


def _group_logfc(df: pd.DataFrame, disease_cols, control_cols) -> pd.Series:
    """Per-gene log fold-change = mean(disease) - mean(control) on the
    log-abundance scale (protein matrix is already log-transformed; pseudobulk
    counts are log1p-CPM transformed before calling this)."""
    return df[disease_cols].mean(axis=1) - df[control_cols].mean(axis=1)


def concordance_protein_transcript(protein_logfc: pd.Series, transcript_logfc: pd.Series,
                                   dataset_name: str) -> pd.DataFrame:
    """Population-level fold-change concordance: correlate per-gene protein logFC
    against per-gene transcript logFC across the shared gene set.

    This replaces the previous positional per-sample Spearman, which paired
    unrelated samples by array position across two independent cohorts and was
    therefore meaningless. No paired-sample design exists between the PXD052937
    plasma proteomics and the scRNA hepatocyte pseudobulk, so concordance is
    reported strictly at the gene (population fold-change) level.
    """
    common_genes = sorted(set(protein_logfc.dropna().index) &
                          set(transcript_logfc.dropna().index))
    print(f"  {dataset_name}: {len(common_genes)} genes with both protein + transcript logFC")

    if len(common_genes) < 5:
        return pd.DataFrame()

    out = pd.DataFrame({
        "gene": common_genes,
        "protein_logFC": protein_logfc.loc[common_genes].values,
        "transcript_logFC": transcript_logfc.loc[common_genes].values,
        "dataset": dataset_name,
    })
    out["direction_match"] = np.sign(out["protein_logFC"]) == np.sign(out["transcript_logFC"])
    return out


def validate_drug_targets(protein_df: pd.DataFrame, targets_file: Path,
                           dataset_name: str) -> pd.DataFrame:
    """Check if TWAS/COLOC-convergent drug targets are detectable at protein level."""
    if not targets_file.exists():
        print(f"  Drug targets file not found: {targets_file}")
        return pd.DataFrame()

    targets = pd.read_csv(targets_file)
    gene_col = None
    for col in ["gene", "Gene", "symbol", "gene_symbol", "HGNC_symbol", "target"]:
        if col in targets.columns:
            gene_col = col
            break
    if gene_col is None:
        print(f"  Could not identify gene column in {targets_file}")
        return pd.DataFrame()

    target_genes = targets[gene_col].dropna().unique()
    found = [g for g in target_genes if g in protein_df.index]
    missing = [g for g in target_genes if g not in protein_df.index]

    print(f"  {dataset_name}: {len(found)}/{len(target_genes)} drug targets detected at protein level")
    if missing:
        print(f"    Missing: {', '.join(missing[:10])}{'...' if len(missing) > 10 else ''}")

    results = []
    for gene in found:
        vals = protein_df.loc[gene].dropna()
        results.append({
            "gene": gene,
            "detected": True,
            "mean_abundance": vals.mean(),
            "std_abundance": vals.std(),
            "n_samples": len(vals),
            "dataset": dataset_name,
        })
    for gene in missing:
        results.append({
            "gene": gene,
            "detected": False,
            "mean_abundance": np.nan,
            "std_abundance": np.nan,
            "n_samples": 0,
            "dataset": dataset_name,
        })

    return pd.DataFrame(results)


def compute_l8_scores(corr_df: pd.DataFrame) -> pd.DataFrame:
    """Compute optional L8_protein evidence scores from gene-level fold-change
    concordance (protein logFC vs transcript logFC per gene)."""
    if corr_df.empty:
        return pd.DataFrame()

    gene_corr = corr_df.groupby("gene").agg(
        mean_protein_logFC=("protein_logFC", "mean"),
        mean_transcript_logFC=("transcript_logFC", "mean"),
        direction_match=("direction_match", "any"),
        n_datasets=("dataset", "nunique"),
    ).reset_index()

    # Concordance magnitude: smaller |protein logFC - transcript logFC| = stronger
    gene_corr["logfc_gap"] = (
        gene_corr["mean_protein_logFC"] - gene_corr["mean_transcript_logFC"]).abs()

    gene_corr["l8_score"] = 0.0
    gene_corr.loc[gene_corr["direction_match"], "l8_score"] += 3
    gene_corr.loc[gene_corr["logfc_gap"] < 0.5, "l8_score"] += 3
    gene_corr.loc[gene_corr["logfc_gap"] < 0.25, "l8_score"] += 2
    gene_corr.loc[gene_corr["n_datasets"] > 1, "l8_score"] += 1

    return gene_corr


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root",
                        default="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
    args = parser.parse_args()

    root = Path(args.project_root)
    prot_dir = root / "Analysis" / "Proteomics" / "results"
    sc_dir = root / "Analysis" / "SingleCell" / "results_gpu_v2"
    me_dir = root / "RNA-seq" / "results" / "multi_evidence"
    drug_dir = root / "RNA-seq" / "results" / "drug_repurposing"

    # Load protein matrix (UniProt-accession indexed) + UniProt->HGNC map
    print("Loading protein matrix (PXD052937)...")
    pxd = pd.read_csv(prot_dir / "pxd052937_protein_matrix.csv", index_col=0)
    print(f"  PXD052937 (UniProt-indexed): {pxd.shape}")

    candidates_f = (root / "data" / "PXD052937" /
                    "20220805_163627_Plasma_liver2" / "Candidates.tsv")
    uniprot_map = load_uniprot_to_hgnc(candidates_f)
    print(f"  Loaded UniProt->HGNC map: {len(uniprot_map)} entries")
    pxd = remap_protein_index_to_hgnc(pxd, uniprot_map)
    print(f"  PXD052937 (HGNC-remapped): {pxd.shape}")

    # Disease/control group labels for protein logFC
    meta_f = prot_dir / "pxd052937_disease_metadata.csv"
    if not meta_f.exists():
        raise FileNotFoundError(f"PXD052937 disease metadata not found: {meta_f}")
    meta = pd.read_csv(meta_f)
    meta["is_masld"] = meta["is_masld"].astype(str).str.upper().isin(["TRUE", "1"])
    disease_cols = [c for c in pxd.columns
                    if c in set(meta.loc[meta["is_masld"], "file_name"])]
    control_cols = [c for c in pxd.columns
                    if c in set(meta.loc[~meta["is_masld"], "file_name"])]
    print(f"  Protein groups: {len(disease_cols)} MASLD, {len(control_cols)} control")
    if not disease_cols or not control_cols:
        raise ValueError("Could not assign MASLD/control protein columns from metadata")
    protein_logfc = _group_logfc(pxd, disease_cols, control_cols)

    # Load hepatocyte pseudobulk (HGNC row index, SRR-accession columns).
    # Explicit filename, RAISE if missing — no silent glob[0] fallback to a
    # non-hepatocyte matrix.
    print("\nLoading hepatocyte pseudobulk transcriptomic data...")
    hep_file = sc_dir / "pseudobulk" / "Hepatocytes_pseudobulk.csv"
    if not hep_file.exists():
        raise FileNotFoundError(f"Hepatocyte pseudobulk not found: {hep_file}")
    hep_pb = pd.read_csv(hep_file, index_col=0)
    print(f"  Hepatocyte pseudobulk: {hep_pb.shape}")

    # Transcript logFC = canonical bulk disease-vs-control logFC (population
    # level). The pseudobulk above has no shared subjects with the plasma
    # proteomics, so we use the canonical bulk logFC for the transcript side
    # rather than positionally pairing pseudobulk columns to protein samples.
    # C2 swap 2026-06-08: dream_results.csv -> canonical_deg_results.csv (LVQW C2).
    dream_f = (root / "RNA-seq" / "Human" / "Patient_Cohorts" / "analysis" /
               "integration" / "results" / "integration" / "canonical_deg_results.csv")
    if not dream_f.exists():
        raise FileNotFoundError(f"canonical bulk DEG results not found: {dream_f}")
    dream = pd.read_csv(dream_f)
    gene_col = "symbol" if "symbol" in dream.columns else "gene"
    transcript_logfc = (dream.dropna(subset=[gene_col])
                             .drop_duplicates(subset=gene_col)
                             .set_index(gene_col)["logFC"])
    # Restrict transcript side to genes detected in the hepatocyte pseudobulk
    transcript_logfc = transcript_logfc[transcript_logfc.index.isin(hep_pb.index)]
    print(f"  Transcript logFC genes (dream ∩ hepatocyte pseudobulk): {len(transcript_logfc)}")

    # Population-level fold-change concordance (NOT positional sample pairing)
    print("\nComputing population-level protein-transcript fold-change concordance...")
    all_corr = concordance_protein_transcript(protein_logfc, transcript_logfc, "PXD052937")

    if not all_corr.empty:
        rho, pval = stats.spearmanr(all_corr["protein_logFC"], all_corr["transcript_logFC"])
        dir_pct = 100 * all_corr["direction_match"].mean()
        print(f"\n  Genes compared: {len(all_corr)}")
        print(f"  Gene-level Spearman rho (protein logFC vs transcript logFC): {rho:.3f} (p={pval:.2e})")
        print(f"  Direction concordance: {dir_pct:.1f}%")
        all_corr.to_csv(prot_dir / "protein_transcript_correlation.csv", index=False)

    # Drug target validation (TWAS/COLOC-convergent targets; MR retired 2026-04-22)
    print("\nValidating drug targets at protein level...")
    targets_file = drug_dir / "convergent_drug_targets.csv"
    all_val = validate_drug_targets(pxd, targets_file, "PXD052937")
    if not all_val.empty:
        all_val.to_csv(prot_dir / "protein_validated_targets.csv", index=False)

    # L8 evaluation
    print("\nEvaluating optional L8_protein evidence layer...")
    l8_scores = compute_l8_scores(all_corr)
    if not l8_scores.empty:
        strong = l8_scores[l8_scores["direction_match"] & (l8_scores["logfc_gap"] < 0.5)]
        pct_strong = 100 * len(strong) / len(l8_scores) if len(l8_scores) > 0 else 0
        print(f"  Genes concordant in direction with logFC gap < 0.5: {len(strong)} ({pct_strong:.1f}%)")

        if pct_strong > 20:
            print("  RECOMMENDATION: Add L8_protein to multi-evidence framework")
            me_file = me_dir / "multi_evidence_scored_genes.csv"
            if me_file.exists():
                me = pd.read_csv(me_file)
                me = me.merge(l8_scores[["gene", "l8_score"]], on="gene", how="left")
                me["l8_score"] = me["l8_score"].fillna(0)
                me.to_csv(me_dir / "multi_evidence_scored_genes_v2.csv", index=False)
                print(f"  Updated scores saved to multi_evidence_scored_genes_v2.csv")
        else:
            print("  RECOMMENDATION: Do NOT add L8 — concordance too weak for formal scoring")

        l8_scores.to_csv(prot_dir / "l8_protein_scores.csv", index=False)

    print("\n=== PROTEO-TRANSCRIPTOMIC INTEGRATION COMPLETE ===")


if __name__ == "__main__":
    main()
