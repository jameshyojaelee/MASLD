#!/usr/bin/env python3
"""Regenerate label transfer comparison plots with harmonized major-type labels."""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import confusion_matrix

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# Label harmonization (same as 05_annotation_sensitivity.py)
LABEL_MAP = {
    "Hepatocyte": "Hepatocyte", "Hepatocytes": "Hepatocyte",
    "Cholangiocyte": "Cholangiocyte", "Cholangiocytes": "Cholangiocyte",
    "Stellate_Cell": "Mesenchymal", "Fibroblasts": "Mesenchymal",
    "Endothelial": "Endothelial", "LSEC": "Endothelial",
    "Endothelial cells": "Endothelial", "Endothelial_cells": "Endothelial",
    "Kupffer_Cell": "Myeloid", "Kupffer cells": "Myeloid",
    "Macrophage": "Myeloid", "Macrophages": "Myeloid",
    "Mono+mono_derived_cells": "Myeloid",
    "NK_T_Cell": "Lymphoid", "B_Cell": "Lymphoid",
    "T cells": "Lymphoid", "T_cells": "Lymphoid",
    "Circulating NK/NKT": "Lymphoid", "Circulating_NK_NKT": "Lymphoid",
    "Resident NK": "Lymphoid", "Resident_NK": "Lymphoid",
    "B cells": "Lymphoid", "Plasma_Cell": "Plasma_cell",
    "Plasma cells": "Plasma_cell", "Plasma_cells": "Plasma_cell",
    "cDC1s": "Myeloid", "cDC2s": "Myeloid", "pDCs": "Myeloid",
    "Neutrophils": "Myeloid", "Basophils": "Myeloid",
    "Mono+mono derived cells": "Myeloid",
}

out_dir = Path("results/label_transfer")
anno = pd.read_csv(out_dir / "cell_annotations_comparison.csv", index_col=0)

# Filter valid cells
valid = ~anno["cell_type_transferred"].isin(["Low_confidence", "Unassigned"])
anno_v = anno[valid].copy()

old_raw = anno_v["cell_type_gene_activity"].astype(str)
new_raw = anno_v["cell_type_transferred"].astype(str)
old_harm = old_raw.map(lambda x: LABEL_MAP.get(x, x))
new_harm = new_raw.map(lambda x: LABEL_MAP.get(x, x))

# --- 1. Harmonized confusion matrix heatmap ---
major_types = sorted(set(old_harm.unique()) | set(new_harm.unique()))
cm = confusion_matrix(old_harm, new_harm, labels=major_types)
cm_df = pd.DataFrame(cm, index=major_types, columns=major_types)
cm_norm = cm_df.div(cm_df.sum(axis=1), axis=0).fillna(0)

fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(cm_norm, annot=True, fmt=".2f", cmap="Blues", ax=ax,
            xticklabels=True, yticklabels=True, linewidths=0.5,
            cbar_kws={"label": "Fraction of old type"})
# Add cell counts as secondary annotation
for i, row_label in enumerate(major_types):
    for j, col_label in enumerate(major_types):
        count = cm_df.iloc[i, j]
        if count > 0:
            ax.text(j + 0.5, i + 0.75, f"n={count:,}", ha="center", va="center",
                    fontsize=7, color="gray")
ax.set_xlabel("Label Transfer (harmonized major type)", fontsize=12)
ax.set_ylabel("Gene Activity (harmonized major type)", fontsize=12)
ax.set_title("Cell Type Reclassification: Gene Activity → Label Transfer\n"
             f"(harmonized concordance = {(old_harm == new_harm).mean():.1%})", fontsize=13)
plt.tight_layout()
fig.savefig(out_dir / "confusion_matrix_harmonized.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out_dir / 'confusion_matrix_harmonized.pdf'}")

# --- 2. Raw cross-mapping heatmap (old → new, informative) ---
raw_labels_old = sorted(old_raw.unique())
raw_labels_new = sorted(new_raw.unique())
cm_raw = confusion_matrix(old_raw, new_raw, labels=raw_labels_old + [l for l in raw_labels_new if l not in raw_labels_old])
# Subset to old-rows × new-cols only (skip zero rows/cols)
old_with_cells = [l for l in raw_labels_old if (old_raw == l).sum() > 0]
new_with_cells = [l for l in raw_labels_new if (new_raw == l).sum() > 0]
cm_cross = pd.DataFrame(0, index=old_with_cells, columns=new_with_cells)
for old_l in old_with_cells:
    for new_l in new_with_cells:
        cm_cross.loc[old_l, new_l] = ((old_raw == old_l) & (new_raw == new_l)).sum()
cm_cross_norm = cm_cross.div(cm_cross.sum(axis=1), axis=0).fillna(0)

fig, ax = plt.subplots(figsize=(12, 8))
sns.heatmap(cm_cross_norm, annot=True, fmt=".2f", cmap="YlOrRd", ax=ax,
            xticklabels=True, yticklabels=True, linewidths=0.5)
ax.set_xlabel("Label Transfer cell type", fontsize=12)
ax.set_ylabel("Gene Activity cell type", fontsize=12)
ax.set_title("Cross-Mapping: Gene Activity → Label Transfer (row-normalized)", fontsize=13)
plt.tight_layout()
fig.savefig(out_dir / "confusion_matrix_crossmap.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out_dir / 'confusion_matrix_crossmap.pdf'}")

# --- 3. Harmonized Jaccard bar plot ---
jaccard = {}
for mt in major_types:
    old_set = set(anno_v.index[old_harm == mt])
    new_set = set(anno_v.index[new_harm == mt])
    union = len(old_set | new_set)
    jaccard[mt] = len(old_set & new_set) / union if union > 0 else 0.0

fig, ax = plt.subplots(figsize=(8, 5))
jac_s = pd.Series(jaccard).sort_values()
colors = ["#2166ac" if v > 0.3 else "#d6604d" for v in jac_s]
jac_s.plot.barh(ax=ax, color=colors)
ax.set_xlabel("Jaccard Index (harmonized major types)")
ax.set_title(f"Per-Major-Type Concordance\n(overall harmonized = {(old_harm == new_harm).mean():.1%})")
ax.axvline(0.5, color="black", linestyle="--", alpha=0.3, label="Jaccard=0.5")
for i, (mt, val) in enumerate(jac_s.items()):
    ax.text(val + 0.01, i, f"{val:.3f}", va="center", fontsize=9)
ax.legend()
plt.tight_layout()
fig.savefig(out_dir / "jaccard_harmonized.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out_dir / 'jaccard_harmonized.pdf'}")

# --- 4. Save harmonized confusion matrix CSV ---
cm_df.to_csv(out_dir / "confusion_matrix_harmonized.csv")
cm_cross.to_csv(out_dir / "confusion_matrix_crossmap.csv")
print(f"Saved CSVs")
print("Done!")
