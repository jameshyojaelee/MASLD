#!/usr/bin/env python
"""Tahoe-100M runner for D1 mechanism (drug-induced expression delta).

D1 task: for each hit gene, predict top-K downstream genes affected by
perturbing it. Tahoe-100M is a drug-perturbation atlas (not a free-form gene-KO
predictor), so the proxy is:

  1. Map each hit gene → drug(s) targeting it via
     `RNA-seq/results/drug_repurposing/convergent_drug_targets.csv` (DGIdb).
  2. Intersect drug list with Tahoe drug vocabulary (~1,100 compounds).
  3. If at least one drug is in Tahoe vocab, scan the pseudobulk DE shards
     at `data/perturbation/datasets/tahoe-100m/metadata/pseudobulk_differential_expression/`
     for rows matching those drugs.
  4. Aggregate log2FC across (cell-line × concentration × plate) by absolute
     mean, take top-K by |aggregated log2FC|, emit as downstream predictions.
  5. If a gene has NO drug in Tahoe vocab → return [] (asymmetric consensus
     coverage; this is handled correctly by consensus_aggregator.py).

For Phase 1 we scan a configurable subset of shards (default 50 of 1026 ~5%)
to keep runtime tractable. Phase 2 (state.tx unlock) will replace this with
real Tahoe model inference and full-shard coverage.

Set TAHOE_N_SHARDS=ALL to scan every shard (slow, ~hours).
"""
from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from _runner_template import D1Adapter, TOP_N_DOWNSTREAM, main_for_model

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
TAHOE_DATASET_DIR = PROJECT_ROOT / "data/perturbation/datasets/tahoe-100m"
PSEUDOBULK_DE_DIR = (
    TAHOE_DATASET_DIR / "metadata" / "pseudobulk_differential_expression"
)
DGIDB_PATH = PROJECT_ROOT / "RNA-seq/results/drug_repurposing/convergent_drug_targets.csv"


class TahoeAdapter(D1Adapter):
    model_name = "tahoe"
    default_checkpoint = str(TAHOE_DATASET_DIR)

    # ------------------------------------------------------------------
    # Load: build gene→drug map + Tahoe vocab + per-drug DE cache
    # ------------------------------------------------------------------
    def load_checkpoint(self) -> None:
        if not TAHOE_DATASET_DIR.exists():
            raise FileNotFoundError(f"Tahoe-100M dataset not at {TAHOE_DATASET_DIR}")
        if not PSEUDOBULK_DE_DIR.exists():
            raise FileNotFoundError(f"Pseudobulk DE dir not at {PSEUDOBULK_DE_DIR}")

        # ----- gene -> drugs from DGIdb / convergent_drug_targets ----------
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
                     if c.lower() in ("dgidb_drugs", "drug", "drug_name", "compound", "drug_id", "compound_name")),
                    None,
                )
                if gene_col and drug_col:
                    for _, row in dgi[[gene_col, drug_col]].dropna().iterrows():
                        g = str(row[gene_col]).strip()
                        drugs_raw = str(row[drug_col])
                        # convergent_drug_targets stores '; '-delimited list
                        drugs = [
                            d.strip().lower() for d in drugs_raw.replace(",", ";").split(";")
                            if d.strip()
                        ]
                        if g and drugs:
                            self._gene_to_drugs[g] = drugs
                    print(f"[tahoe] gene→drug map: {len(self._gene_to_drugs)} genes from DGIdb")
                else:
                    warnings.warn(f"DGIdb gene/drug cols not found: {list(dgi.columns)}")
            except Exception as e:
                warnings.warn(f"Failed to parse DGIdb: {e}")

        # ----- Tahoe drug vocab from drug_metadata.parquet ----------------
        self._tahoe_drug_index: set[str] = self._load_tahoe_drug_vocab()
        print(f"[tahoe] Tahoe drug vocabulary: {len(self._tahoe_drug_index)} compounds")

        # ----- Identify hits with drug coverage in Tahoe -----------------
        # We only need pseudobulk DE for these drugs. Build the union of all
        # relevant drugs across our hits (will be filled lazily on first use).
        self._drug_de_cache: dict[str, pd.DataFrame] = {}
        self._needed_drugs: set[str] = set()
        for gene, drugs in self._gene_to_drugs.items():
            for d in drugs:
                if d in self._tahoe_drug_index:
                    self._needed_drugs.add(d)
        print(
            f"[tahoe] gene→drug intersect Tahoe vocab: {len(self._needed_drugs)} "
            "candidate drugs across gene→drug map"
        )

        # ----- Pre-scan pseudobulk DE shards once -------------------------
        # Single pass to materialize per-drug DE tables for `_needed_drugs`.
        # Cap shards for tractability; default 50/1026 (~5%).
        n_shards_env = os.environ.get("TAHOE_N_SHARDS", "50")
        all_shards = sorted(PSEUDOBULK_DE_DIR.glob("train-*.parquet"))
        if n_shards_env.upper() == "ALL":
            shards = all_shards
        else:
            try:
                n_shards = int(n_shards_env)
            except ValueError:
                n_shards = 50
            shards = all_shards[:max(1, n_shards)]
        print(
            f"[tahoe] scanning {len(shards)} / {len(all_shards)} pseudobulk DE shards "
            f"(TAHOE_N_SHARDS={n_shards_env})"
        )

        self._n_shards_scanned = 0
        self._n_de_rows = 0
        if self._needed_drugs:
            self._scan_pseudobulk_for_drugs(shards)
        print(
            f"[tahoe] DE cache: {len(self._drug_de_cache)} / {len(self._needed_drugs)} "
            f"drugs covered ({self._n_de_rows:,} DE rows total)"
        )

        # Stats trackers
        self._n_attempted = 0
        self._n_skipped_no_drug = 0
        self._n_skipped_no_de = 0

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
        warnings.warn(f"No drug-metadata file found under {TAHOE_DATASET_DIR}")
        return set()

    def _scan_pseudobulk_for_drugs(self, shards: list[Path]) -> None:
        """Single-pass scan: per-shard groupby on (drug, gene), accumulate."""
        # Per-drug running sums: dict[drug] -> dict[gene] -> [sum_lfc, sum_abs, n]
        per_drug: dict[str, dict[str, list[float]]] = {}
        needed_lc = self._needed_drugs

        for i, shard in enumerate(shards):
            try:
                df = pd.read_parquet(
                    shard,
                    columns=["gene_name", "log2FoldChange", "drug"],
                )
            except Exception as e:
                warnings.warn(f"failed reading shard {shard.name}: {e}")
                continue

            df["drug_lc"] = df["drug"].astype(str).str.lower()
            sub = df[df["drug_lc"].isin(needed_lc)]
            sub = sub.dropna(subset=["log2FoldChange", "gene_name"])
            if sub.empty:
                self._n_shards_scanned += 1
                continue

            # Vectorized per-shard aggregation: groupby (drug, gene_name)
            agg = (
                sub.assign(abs_lfc=sub["log2FoldChange"].abs())
                .groupby(["drug_lc", "gene_name"], sort=False)
                .agg(
                    sum_lfc=("log2FoldChange", "sum"),
                    sum_abs=("abs_lfc", "sum"),
                    n=("log2FoldChange", "count"),
                )
                .reset_index()
            )
            self._n_de_rows += len(sub)

            # Accumulate into per_drug running sums
            for drug, sub_d in agg.groupby("drug_lc", sort=False):
                buckets = per_drug.setdefault(drug, {})
                for g, sl, sa, n in zip(
                    sub_d["gene_name"].astype(str),
                    sub_d["sum_lfc"].astype(float),
                    sub_d["sum_abs"].astype(float),
                    sub_d["n"].astype(int),
                ):
                    if g in buckets:
                        b = buckets[g]
                        b[0] += sl
                        b[1] += sa
                        b[2] += n
                    else:
                        buckets[g] = [sl, sa, n]
            self._n_shards_scanned += 1

            if (i + 1) % 10 == 0:
                covered = sum(1 for d, b in per_drug.items() if b)
                print(
                    f"[tahoe] scanned {i+1}/{len(shards)} shards; "
                    f"{covered} drugs with ≥1 row; {self._n_de_rows:,} rows so far"
                )

        # Materialize per-drug DE tables
        for drug, buckets in per_drug.items():
            if not buckets:
                continue
            genes, sl_arr, sa_arr, n_arr = zip(*[
                (g, b[0], b[1], b[2]) for g, b in buckets.items()
            ])
            sl = np.asarray(sl_arr, dtype=np.float64)
            sa = np.asarray(sa_arr, dtype=np.float64)
            n = np.asarray(n_arr, dtype=np.int64)
            self._drug_de_cache[drug] = pd.DataFrame(
                dict(
                    gene_name=list(genes),
                    log2FC_mean=sl / np.maximum(n, 1),
                    abs_lfc_mean=sa / np.maximum(n, 1),
                )
            )

    # ------------------------------------------------------------------
    # Prediction
    # ------------------------------------------------------------------
    def predict_for_hit(self, *, hit_row, context, substrate):
        target = hit_row.get("gene")
        if not isinstance(target, str) or not target:
            return []

        self._n_attempted += 1
        if self._n_attempted % 50 == 0:
            print(
                f"[tahoe] processed {self._n_attempted} "
                f"(no_drug={self._n_skipped_no_drug}, no_de={self._n_skipped_no_de})"
            )

        drugs = self._gene_to_drugs.get(target, [])
        if not drugs:
            self._n_skipped_no_drug += 1
            return []

        # Restrict to drugs in Tahoe vocab
        drugs_tahoe = [d for d in drugs if d in self._tahoe_drug_index]
        if not drugs_tahoe:
            self._n_skipped_no_drug += 1
            return []

        # Restrict to drugs we actually have DE rows for
        drugs_with_de = [d for d in drugs_tahoe if d in self._drug_de_cache]
        if not drugs_with_de:
            self._n_skipped_no_de += 1
            return []

        # Aggregate DE across all matching drugs (mean of means, weighted equally)
        frames = [self._drug_de_cache[d] for d in drugs_with_de]
        all_de = pd.concat(frames, ignore_index=True)
        agg = (
            all_de.groupby("gene_name")
            .agg(log2FC_mean=("log2FC_mean", "mean"), abs_lfc_mean=("abs_lfc_mean", "mean"))
            .reset_index()
        )

        # Drop self (drug target's own gene), drop NaN, sort by |LFC|
        agg = agg[agg["gene_name"].astype(str) != target]
        agg = agg.dropna(subset=["log2FC_mean"])
        agg["abs_lfc"] = agg["log2FC_mean"].abs()
        top = agg.sort_values("abs_lfc", ascending=False).head(TOP_N_DOWNSTREAM)

        if top.empty:
            self._n_skipped_no_de += 1
            return []

        # Confidence: rescale abs_lfc to [0,1] within this gene's top-K.
        # Use tanh squashing so |LFC|~1 -> ~0.76 confidence; |LFC|~3 -> ~0.99.
        rows: list[dict] = []
        max_abs = float(top["abs_lfc"].max()) or 1.0
        for rank, (_, r) in enumerate(top.iterrows(), start=1):
            lfc = float(r["log2FC_mean"])
            conf = float(np.tanh(abs(lfc) / max(0.5, max_abs * 0.5)))
            rows.append(
                dict(
                    target_gene=str(target),
                    downstream_gene=str(r["gene_name"]),
                    logFC_predicted=lfc,
                    abs_rank=int(rank),
                    direction="up" if lfc > 0 else "down",
                    confidence=max(0.0, min(1.0, conf)),
                )
            )
        return rows


if __name__ == "__main__":
    main_for_model(TahoeAdapter)
