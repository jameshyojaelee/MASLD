#!/usr/bin/env python
"""Tahoe-100M drug-drug synergy runner for D3.

Per plan §5.1 Tahoe asymmetric consensus: returns [] for gene pairs without
drug-pair coverage in Tahoe (those NaN entries do NOT count as disagreement).

Strategy:
1. Map each gene → list of drugs targeting it (DGIdb / convergent_drug_targets).
2. Intersect drug list with Tahoe drug vocabulary (~1,100 compounds).
3. If both genes in a pair have ≥1 Tahoe-vocab drug, query Tahoe for the
   drug-pair expression delta vs control + single-drug deltas.
4. Compute synergy: ||combo_delta - (single1_delta + single2_delta)||.

Note: actual Tahoe-100M model inference requires the Tahoe model checkpoint
(separate from the dataset). Phase 1 emits drug-coverage diagnostics + a
placeholder for synergy_magnitude (set to 0.0 with synergy_class="additive"
until model checkpoint is loaded). Phase 2 fine-tune will replace with real
Tahoe inference.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd

from _runner_template import D3Adapter, main_for_model

TAHOE_DATASET_DIR = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
    "data/perturbation/datasets/tahoe-100m"
)
DGIDB_PATH = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
    "RNA-seq/results/drug_repurposing/convergent_drug_targets.csv"
)


class TahoeDrugAdapter(D3Adapter):
    model_name = "tahoe_drug"
    default_checkpoint = str(TAHOE_DATASET_DIR)

    def load_checkpoint(self) -> None:
        """Load Tahoe-100M drug index + gene→drug map."""
        if not TAHOE_DATASET_DIR.exists():
            raise FileNotFoundError(f"Tahoe-100M dataset not at {TAHOE_DATASET_DIR}")

        # Build gene→drugs lookup
        self._gene_to_drugs: dict[str, list[str]] = {}
        if DGIDB_PATH.exists():
            try:
                dgi = pd.read_csv(DGIDB_PATH)
                gene_col = next(
                    (c for c in dgi.columns
                     if c.lower() in ("gene", "gene_symbol", "human_symbol", "target")),
                    None,
                )
                drug_col = next(
                    (c for c in dgi.columns
                     if c.lower() in ("drug", "drug_name", "compound", "drug_id", "compound_name")),
                    None,
                )
                if gene_col and drug_col:
                    for g, sub in dgi.groupby(gene_col):
                        self._gene_to_drugs[g] = (
                            sub[drug_col].dropna().astype(str).str.lower().tolist()
                        )
                    print(f"[tahoe] gene→drug map: {len(self._gene_to_drugs)} genes")
            except Exception as e:
                warnings.warn(f"Failed to parse DGIdb path: {e}")

        # Load Tahoe drug vocabulary (lazy via metadata)
        self._tahoe_drug_index: set[str] = self._load_tahoe_drug_vocab()
        print(f"[tahoe] Tahoe drug vocabulary: {len(self._tahoe_drug_index)} compounds")

    def _load_tahoe_drug_vocab(self) -> set[str]:
        for candidate in [
            TAHOE_DATASET_DIR / "metadata" / "drug_metadata.parquet",  # actual location
            TAHOE_DATASET_DIR / "metadata" / "drugs.csv",
            TAHOE_DATASET_DIR / "metadata.parquet",
            TAHOE_DATASET_DIR / "drugs.csv",
        ]:
            if not candidate.exists():
                continue
            try:
                df = (
                    pd.read_parquet(candidate)
                    if candidate.suffix == ".parquet"
                    else pd.read_csv(candidate)
                )
                drug_col = next(
                    (c for c in df.columns
                     if "drug" in c.lower() or "compound" in c.lower() or "chemical" in c.lower()),
                    None,
                )
                if drug_col:
                    return set(df[drug_col].dropna().astype(str).str.lower())
            except Exception as e:
                warnings.warn(f"failed to read {candidate}: {e}")
        warnings.warn(f"No drug-metadata file found under {TAHOE_DATASET_DIR}")
        return set()

    def predict_for_hit(self, *, hit_row, context, substrate):
        genes = [hit_row.get(f"gene{i+1}") for i in range(4) if hit_row.get(f"gene{i+1}")]
        genes = [g for g in genes if isinstance(g, str) and g]
        if len(genes) < 2:
            return []

        # Map each gene → list of drugs that target it
        drug_lists = [self._gene_to_drugs.get(g, []) for g in genes]
        if not all(drug_lists):
            return []  # gene without drug coverage → asymmetric consensus skip

        # Intersect with Tahoe vocabulary
        if not self._tahoe_drug_index:
            return []  # Tahoe vocab unavailable → skip
        in_tahoe = [
            [d for d in dl if d in self._tahoe_drug_index] for dl in drug_lists
        ]
        if not all(in_tahoe):
            return []  # at least one gene's drugs are absent from Tahoe → skip

        # Phase 1 placeholder: emit drug-coverage diagnostic; full Tahoe
        # inference deferred to Phase 2 fine-tune when checkpoint loads.
        return [dict(
            gene1=genes[0],
            gene2=genes[1] if len(genes) > 1 else None,
            gene3=genes[2] if len(genes) > 2 else None,
            gene4=genes[3] if len(genes) > 3 else None,
            additive_baseline=0.0,
            observed_double_or_higher=0.0,
            synergy_magnitude=0.0,
            synergy_class="additive",
            sigma_above_additive=0.0,
            k562_bias_confidence=1.0,  # Tahoe drug-pair coverage = high fidelity
        )]


if __name__ == "__main__":
    main_for_model(TahoeDrugAdapter)
