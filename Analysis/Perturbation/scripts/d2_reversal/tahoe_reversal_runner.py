#!/usr/bin/env python
"""Tahoe-100M runner for D2 reversal — drug-pseudobulk vs ref_LFC cosine.

D2 task: rank atlas genes by predicted reversal of the Diseased→Healthy
hepatocyte signature. Tahoe-100M is a *drug-perturbation* atlas (no free-form
gene KO), so the proxy is:

  1. Map each hit gene → drug(s) that target it
     (RNA-seq/results/drug_repurposing/convergent_drug_targets.csv, DGIdb).
  2. Intersect drug list with the Tahoe drug vocabulary (~380 compounds).
  3. For each (gene, drug-hitting-it) tuple, retrieve the drug's pseudobulk
     DE profile across the genome from Tahoe's pseudobulk DE shards. We
     **prefer HepG2/C3A** rows (Tahoe's only liver-derived line; lives in
     shards 432-451) for tissue-relevance; if a drug has no HepG2/C3A row we
     fall back to mean across all available cell lines for that drug.
  4. The drug-induced DE profile is (gene_name -> log2FoldChange). Compute
     cosine similarity with the reference LFC vector on the shared gene
     intersection. Reversal score = - cosine. High positive = drug pushes
     expression *opposite* to the disease axis → reverses signature.

Asymmetric consensus: genes with NO Tahoe-vocab drug → return []. These
absences do NOT count as disagreement during D2 consensus (per plan §5.1).

Set TAHOE_HEPG2_ONLY=1 to restrict to HepG2/C3A shards only (default behaviour;
gives ~80 drugs with coverage). Set TAHOE_HEPG2_ONLY=0 to scan all 1,026 shards
across all cell lines (slow ~hours; broader drug coverage but less tissue-relevant).

CLI:
    python tahoe_reversal_runner.py --modality zero_shot --context all --reference ref_a
"""
from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from _runner_template import D2Adapter, main_for_model

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
TAHOE_DATASET_DIR = PROJECT_ROOT / "data/perturbation/datasets/tahoe-100m"
PSEUDOBULK_DE_DIR = (
    TAHOE_DATASET_DIR / "metadata" / "pseudobulk_differential_expression"
)
DGIDB_PATH = PROJECT_ROOT / "RNA-seq/results/drug_repurposing/convergent_drug_targets.csv"

# HepG2/C3A is on shards 432-451 (20 contiguous shards, see scan in task notes).
# Sticking to liver lines keeps the proxy tissue-relevant; other lines (A549,
# HEC-1-A, etc) are kept off by default.
HEPG2_SHARD_RANGE = (432, 452)


class TahoeReversalAdapter(D2Adapter):
    model_name = "tahoe_reversal"
    default_checkpoint = str(TAHOE_DATASET_DIR)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Initialise trackers + lookup tables up-front so dry-run (which skips
        # load_checkpoint) doesn't AttributeError on first predict_for_hit.
        self._n_attempted = 0
        self._n_skipped_no_drug = 0
        self._n_skipped_no_de = 0
        self._n_skipped_empty_overlap = 0
        self._n_shards_scanned = 0
        self._n_de_rows = 0
        self._gene_to_drugs: dict[str, list[str]] = {}
        self._tahoe_drug_index: set[str] = set()
        self._needed_drugs: set[str] = set()
        self._drug_de: dict[str, pd.Series] = {}
        self._ref_lfc: dict[str, float] = {}
        self._ref_series: pd.Series = pd.Series(dtype=float)

    # ------------------------------------------------------------------
    # Load: gene→drug map + Tahoe vocab + pre-scan pseudobulk DE shards
    # ------------------------------------------------------------------
    def load_checkpoint(self) -> None:
        if not TAHOE_DATASET_DIR.exists():
            raise FileNotFoundError(f"Tahoe-100M dataset not at {TAHOE_DATASET_DIR}")
        if not PSEUDOBULK_DE_DIR.exists():
            raise FileNotFoundError(f"Pseudobulk DE dir not at {PSEUDOBULK_DE_DIR}")

        # ----- Reference vector (gene -> LFC) ----------------------------
        # The D2Adapter base loaded self._ref_df. Build a per-gene LFC dict
        # restricted to gene-symbol rows (drop ENSG… leftovers in ref_b).
        ref = self._ref_df.dropna(subset=["LFC", "gene"]).copy()
        ref["gene"] = ref["gene"].astype(str)
        ref = ref[~ref["gene"].str.startswith("ENSG")]
        self._ref_lfc: dict[str, float] = dict(zip(ref["gene"], ref["LFC"].astype(float)))
        print(
            f"[tahoe_reversal] ref={self.reference_signature} usable LFC entries: "
            f"{len(self._ref_lfc):,}"
        )

        # ----- gene -> drugs via DGIdb / convergent_drug_targets ---------
        self._gene_to_drugs: dict[str, list[str]] = {}
        if DGIDB_PATH.exists():
            try:
                dgi = pd.read_csv(DGIDB_PATH, low_memory=False)
                gene_col = next(
                    (c for c in dgi.columns
                     if c.lower() in ("symbol", "gene", "gene_symbol", "human_symbol", "target")),
                    None,
                )
                drug_col = next(
                    (c for c in dgi.columns
                     if c.lower() in ("dgidb_drugs", "drug", "drug_name",
                                      "compound", "drug_id", "compound_name")),
                    None,
                )
                if gene_col and drug_col:
                    for _, row in dgi[[gene_col, drug_col]].dropna().iterrows():
                        g = str(row[gene_col]).strip()
                        drugs_raw = str(row[drug_col])
                        drugs = [
                            d.strip().lower()
                            for d in drugs_raw.replace(",", ";").split(";")
                            if d.strip()
                        ]
                        if g and drugs:
                            self._gene_to_drugs[g] = drugs
                    print(
                        f"[tahoe_reversal] gene→drug map: "
                        f"{len(self._gene_to_drugs)} genes from DGIdb"
                    )
                else:
                    warnings.warn(f"DGIdb gene/drug cols not found: {list(dgi.columns)}")
            except Exception as e:
                warnings.warn(f"Failed to parse DGIdb: {e}")

        # ----- Tahoe drug vocab ------------------------------------------
        self._tahoe_drug_index: set[str] = self._load_tahoe_drug_vocab()
        print(
            f"[tahoe_reversal] Tahoe drug vocabulary: "
            f"{len(self._tahoe_drug_index)} compounds"
        )

        # ----- Identify needed drugs (gene→drug ∩ Tahoe vocab) -----------
        self._needed_drugs: set[str] = set()
        for gene, drugs in self._gene_to_drugs.items():
            for d in drugs:
                if d in self._tahoe_drug_index:
                    self._needed_drugs.add(d)
        print(
            f"[tahoe_reversal] gene→drug ∩ Tahoe vocab: "
            f"{len(self._needed_drugs)} candidate drugs"
        )

        # ----- Pre-scan pseudobulk DE shards once ------------------------
        # Default: HepG2/C3A only (shards 432-451). TAHOE_HEPG2_ONLY=0 → all shards.
        hepg2_only = os.environ.get("TAHOE_HEPG2_ONLY", "1") not in ("0", "false", "False")
        all_shards = sorted(PSEUDOBULK_DE_DIR.glob("train-*.parquet"))
        if hepg2_only:
            shards = all_shards[HEPG2_SHARD_RANGE[0]:HEPG2_SHARD_RANGE[1]]
            print(
                f"[tahoe_reversal] HepG2-only mode: {len(shards)} shards "
                f"(indices {HEPG2_SHARD_RANGE[0]}-{HEPG2_SHARD_RANGE[1]-1})"
            )
        else:
            shards = all_shards
            print(
                f"[tahoe_reversal] all-cell-line mode: {len(shards)} shards "
                "(slow; broader drug coverage)"
            )

        self._n_shards_scanned = 0
        self._n_de_rows = 0
        # _drug_de[drug] = pd.Series(index=gene_name, values=mean log2FC across rows)
        self._drug_de: dict[str, pd.Series] = {}
        if self._needed_drugs:
            self._scan_pseudobulk(shards, hepg2_only=hepg2_only)
        print(
            f"[tahoe_reversal] DE cache: {len(self._drug_de)} / "
            f"{len(self._needed_drugs)} drugs covered ({self._n_de_rows:,} DE rows)"
        )

        # ----- Pre-compute ref vector aligned with pooled Tahoe gene index
        # We'll align per-hit (since different drugs may cover different genes),
        # but cache the full ref-LFC pandas Series for fast intersection.
        self._ref_series = pd.Series(self._ref_lfc, name="ref_LFC")

    def _load_tahoe_drug_vocab(self) -> set[str]:
        candidates = [
            TAHOE_DATASET_DIR / "metadata" / "drug_metadata.parquet",
            TAHOE_DATASET_DIR / "metadata" / "drugs.csv",
        ]
        for path in candidates:
            if not path.exists():
                continue
            try:
                df = (
                    pd.read_parquet(path)
                    if path.suffix == ".parquet"
                    else pd.read_csv(path)
                )
                drug_col = next(
                    (c for c in df.columns
                     if "drug" in c.lower() or "compound" in c.lower() or "chemical" in c.lower()),
                    None,
                )
                if drug_col:
                    return set(df[drug_col].dropna().astype(str).str.lower())
            except Exception as e:
                warnings.warn(f"failed reading {path}: {e}")
        warnings.warn(f"No drug-metadata file under {TAHOE_DATASET_DIR}")
        return set()

    def _scan_pseudobulk(self, shards: list[Path], hepg2_only: bool) -> None:
        """One-pass scan: build {drug -> {gene -> mean log2FoldChange}}.

        Aggregation is mean log2FC across whatever (cell_line × concentration ×
        plate) combinations land in the scanned shards. In HepG2-only mode this
        gives a clean HepG2-specific DE profile per drug.
        """
        per_drug_buckets: dict[str, dict[str, list[float]]] = {
            d: {} for d in self._needed_drugs
        }
        needed = self._needed_drugs

        cols = ["gene_name", "log2FoldChange", "drug", "Cell_Name_Vevo"]
        for i, shard in enumerate(shards):
            try:
                df = pd.read_parquet(shard, columns=cols)
            except Exception as e:
                warnings.warn(f"failed reading shard {shard.name}: {e}")
                continue

            if hepg2_only:
                df = df[df["Cell_Name_Vevo"].astype(str) == "HepG2/C3A"]
                if df.empty:
                    self._n_shards_scanned += 1
                    continue

            df["drug_lc"] = df["drug"].astype(str).str.lower()
            sub = df[df["drug_lc"].isin(needed)]
            sub = sub.dropna(subset=["log2FoldChange", "gene_name"])

            for drug, ddf in sub.groupby("drug_lc"):
                buckets = per_drug_buckets.setdefault(drug, {})
                for g, lfc in zip(
                    ddf["gene_name"].astype(str),
                    ddf["log2FoldChange"].astype(float),
                ):
                    buckets.setdefault(g, []).append(lfc)
                self._n_de_rows += len(ddf)
            self._n_shards_scanned += 1

            if (i + 1) % 10 == 0:
                covered = sum(1 for b in per_drug_buckets.values() if b)
                print(
                    f"[tahoe_reversal] scanned {i+1}/{len(shards)} shards; "
                    f"{covered} drugs with ≥1 row"
                )

        for drug, buckets in per_drug_buckets.items():
            if not buckets:
                continue
            series = pd.Series(
                {g: float(np.mean(v)) for g, v in buckets.items()},
                name="log2FC_mean",
            )
            self._drug_de[drug] = series

    @staticmethod
    def _cosine(a: np.ndarray, b: np.ndarray) -> float:
        na = float(np.linalg.norm(a)) + 1e-12
        nb = float(np.linalg.norm(b)) + 1e-12
        return float(np.dot(a, b) / (na * nb))

    # ------------------------------------------------------------------
    # Prediction
    # ------------------------------------------------------------------
    def predict_for_hit(self, *, hit_row, context, substrate):
        gene = hit_row.get("gene")
        if not isinstance(gene, str) or not gene:
            return []

        self._n_attempted += 1
        if self._n_attempted % 100 == 0:
            print(
                f"[tahoe_reversal] processed {self._n_attempted} "
                f"(no_drug={self._n_skipped_no_drug}, "
                f"no_de={self._n_skipped_no_de}, "
                f"no_overlap={self._n_skipped_empty_overlap})"
            )

        drugs = self._gene_to_drugs.get(gene, [])
        if not drugs:
            self._n_skipped_no_drug += 1
            return []
        drugs_tahoe = [d for d in drugs if d in self._tahoe_drug_index]
        if not drugs_tahoe:
            self._n_skipped_no_drug += 1
            return []
        drugs_with_de = [d for d in drugs_tahoe if d in self._drug_de]
        if not drugs_with_de:
            self._n_skipped_no_de += 1
            return []

        # Aggregate DE across all matching drugs (per-gene mean of means)
        frames = [self._drug_de[d] for d in drugs_with_de]
        combo = pd.concat(frames, axis=1).mean(axis=1)
        # Intersect with reference vector
        common = combo.index.intersection(self._ref_series.index)
        if len(common) < 50:
            # Too few overlapping genes → unstable cosine; skip.
            self._n_skipped_empty_overlap += 1
            return []
        drug_vec = combo.loc[common].astype(np.float32).values
        ref_vec = self._ref_series.loc[common].astype(np.float32).values

        cos = self._cosine(drug_vec, ref_vec)
        # Reversal score: -cosine. Positive = drug pushes expression *opposite*
        # to disease direction → reversal candidate.
        reversal_score = float(-cos)

        return [
            dict(
                gene=str(gene),
                reversal_score=reversal_score,
                reversal_rank=0,  # rank assigned by template post-processor
                stage_specific="—",
                cell_type="hepatocyte_progressor",
                reference_signature=self.reference_signature,
            )
        ]


if __name__ == "__main__":
    main_for_model(TahoeReversalAdapter)
