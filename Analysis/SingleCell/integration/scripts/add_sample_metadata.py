#!/usr/bin/env python3
"""Add preparation_method and condition_harmonized to integrated atlas h5ad.

Fixes:
  1. Liver_Atlas condition: "Mixed" → "Healthy" (all 48 samples are healthy donors)
  2. Adds obs["preparation_method"]: nuclei, cd45_negative, cd45_positive,
     dc_sorted, monos_macs_sorted, unsorted
  3. Adds obs["condition_harmonized"]: Healthy vs MASLD (collapses NAFLD/NASH)

Uses h5py for direct HDF5 modification — avoids loading the full expression
matrix into memory. Safe for 25GB+ h5ad files on CPU nodes.

Also regenerates cell_type_proportions.csv with new metadata columns.
"""

import os
import sys
import csv
import numpy as np
import pandas as pd
import h5py

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)

SCVI_DIR = os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2")
H5AD_FILE = os.path.join(SCVI_DIR, "integrated_atlas.h5ad")
PROPS_FILE = os.path.join(SCVI_DIR, "cell_type_proportions.csv")

# ---------------------------------------------------------------------------
# SRR → preparation_method lookup (from SraRunTable.csv × GSE192740_sampleInfo)
# ---------------------------------------------------------------------------

NUCLEI_SRRS = {
    "SRR17375011",  # ABU8
    "SRR17375051",  # CS161
    "SRR17375052",  # CS162
    "SRR17375053",  # CS164
    "SRR17375054",  # CS166
    "SRR17375055",  # CS167
    "SRR17375056",  # CS169
    "SRR17375057",  # CS170
    "SRR17375058",  # CS171
}

CD45_NEG_SRRS = {
    "SRR17375020",  # CS31
    "SRR17375047",  # CS111
    "SRR17375050",  # CS127
}

DC_SORTED_SRRS = {
    "SRR17375021",  # CS32
    "SRR17375022",  # CS33
}

MONOS_MACS_SORTED_SRRS = {
    "SRR17375023",  # CS34
}

# Everything else in Liver_Atlas human is CD45+ sorted
# SRR17375012-019, SRR17375024-046 (excl 020/021/022/023/047/050),
# SRR17375048-049 = 33 SRRs

CONDITION_HARMONIZE = {
    "Healthy": "Healthy",
    "Mixed": "Healthy",  # Liver_Atlas healthy donors mislabeled as Mixed
    "MASLD": "MASLD",
    "NAFLD": "MASLD",
    "NASH": "MASLD",
    "Cirrhotic": "MASLD",  # GSE136103 (Ramachandran 2019) cirrhosis → MASLD
}


def read_obs_column(f, col_name):
    """Read an obs column from h5ad, handling both categorical and plain arrays."""
    group_path = f"/obs/{col_name}"
    obj = f[group_path]
    # Check if it's a categorical group (has codes + categories)
    if isinstance(obj, h5py.Group) and "codes" in obj and "categories" in obj:
        codes = obj["codes"][:]
        # Use .asstr() for string categories (h5py 3.x compatible)
        try:
            cats = obj["categories"].asstr()[:]
        except Exception:
            cats = obj["categories"][:]
            if cats.dtype.kind == "S":
                cats = np.array([c.decode("utf-8") for c in cats])
            elif cats.dtype.kind == "O":
                cats = np.array([str(c) for c in cats])
        return cats[codes]
    # Plain dataset
    try:
        vals = obj.asstr()[:]
    except Exception:
        vals = obj[:]
        if vals.dtype.kind in ("S", "O"):
            return np.array([
                v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in vals
            ])
    return vals


def write_categorical_column(f, col_name, values):
    """Write a categorical column to h5ad obs."""
    group_path = f"obs/{col_name}"
    # Remove existing
    if col_name in f["obs"]:
        del f["obs"][col_name]

    cats = np.unique(values)
    cat_to_code = {c: i for i, c in enumerate(cats)}
    codes = np.array([cat_to_code[v] for v in values], dtype=np.int8)

    grp = f["obs"].create_group(col_name)
    grp.create_dataset("codes", data=codes, compression="gzip")
    # Store categories as variable-length UTF-8 strings
    dt = h5py.string_dtype()
    grp.create_dataset("categories", data=cats.astype(object), dtype=dt)
    grp.attrs["encoding-type"] = "categorical"
    grp.attrs["encoding-version"] = "0.2.0"


## GSE136103 (Ramachandran 2019): GSM → preparation_method
## CD45+ = immune-enriched (sorted), CD45- = non-immune (sorted)
GSE136103_CD45_NEG = {
    "GSM4041151", "GSM4041152",  # Healthy1 CD45-A, CD45-B
    "GSM4041154",                # Healthy2 CD45-
    "GSM4041156", "GSM4041157",  # Healthy3 CD45-A, CD45-B
    "GSM4041159",                # Healthy4 CD45-
    "GSM4041162", "GSM4041163",  # Cirrhotic1 CD45-A, CD45-B
    "GSM4041165",                # Cirrhotic2 CD45-
    "GSM4041167",                # Cirrhotic3 CD45-
}
GSE136103_CD45_POS = {
    "GSM4041150",  # Healthy1 CD45+
    "GSM4041153",  # Healthy2 CD45+
    "GSM4041155",  # Healthy3 CD45+
    "GSM4041158",  # Healthy4 CD45+
    "GSM4041160",  # Healthy5 CD45+
    "GSM4041161",  # Cirrhotic1 CD45+
    "GSM4041164",  # Cirrhotic2 CD45+
    "GSM4041166",  # Cirrhotic3 CD45+
    "GSM4041168",  # Cirrhotic4 CD45+
    "GSM4041169",  # Cirrhotic5 CD45+
}


def classify_preparation(sample_id, dataset):
    """Classify a sample's preparation method."""
    # GSE136103: enzymatic dissociation + FACS (CD45+/CD45-)
    if dataset == "GSE136103":
        if sample_id in GSE136103_CD45_NEG:
            return "cd45_negative"
        if sample_id in GSE136103_CD45_POS:
            return "cd45_positive"
        return "unsorted"  # fallback
    if dataset != "Liver_Atlas":
        return "unsorted"
    if sample_id in NUCLEI_SRRS:
        return "nuclei"
    if sample_id in CD45_NEG_SRRS:
        return "cd45_negative"
    if sample_id in DC_SORTED_SRRS:
        return "dc_sorted"
    if sample_id in MONOS_MACS_SORTED_SRRS:
        return "monos_macs_sorted"
    return "cd45_positive"


def main():
    print(f"Reading h5ad: {H5AD_FILE}")
    assert os.path.exists(H5AD_FILE), f"h5ad not found: {H5AD_FILE}"

    # Resolve symlink for the actual file path
    real_path = os.path.realpath(H5AD_FILE)
    print(f"  Resolved: {real_path}")
    print(f"  Size: {os.path.getsize(real_path) / 1e9:.1f} GB")

    # --- Read existing obs columns ---
    with h5py.File(real_path, "r") as f:
        samples = read_obs_column(f, "sample")
        datasets = read_obs_column(f, "dataset")
        conditions = read_obs_column(f, "condition")

    n_cells = len(samples)
    print(f"  Total cells: {n_cells:,}")

    # --- Compute preparation_method ---
    print("Computing preparation_method...")
    prep_methods = np.array([
        classify_preparation(s, d) for s, d in zip(samples, datasets)
    ])
    from collections import Counter
    prep_counts = Counter(prep_methods)
    for k, v in sorted(prep_counts.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v:,} cells")

    # --- Fix condition: Liver_Atlas "Mixed" → "Healthy" ---
    print("Fixing condition labels...")
    new_conditions = conditions.copy()
    mixed_mask = (datasets == "Liver_Atlas") & (conditions == "Mixed")
    n_fixed = mixed_mask.sum()
    new_conditions[mixed_mask] = "Healthy"
    print(f"  Fixed {n_fixed:,} cells: Liver_Atlas Mixed → Healthy")

    cond_counts = Counter(new_conditions)
    for k, v in sorted(cond_counts.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v:,} cells")

    # --- Compute condition_harmonized ---
    print("Computing condition_harmonized...")
    cond_harmonized = np.array([
        CONDITION_HARMONIZE.get(c, "Unknown") for c in new_conditions
    ])
    harm_counts = Counter(cond_harmonized)
    for k, v in sorted(harm_counts.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v:,} cells")

    # --- Write new/updated columns to h5ad ---
    print(f"Writing columns to {real_path}...")
    with h5py.File(real_path, "r+") as f:
        write_categorical_column(f, "preparation_method", prep_methods)
        write_categorical_column(f, "condition", new_conditions)
        write_categorical_column(f, "condition_harmonized", cond_harmonized)
    print("  Done writing h5ad columns.")

    # --- Regenerate cell_type_proportions.csv ---
    print("Regenerating cell_type_proportions.csv...")
    with h5py.File(real_path, "r") as f:
        cell_types = read_obs_column(f, "cell_type")
        # Re-read updated columns
        samples_upd = read_obs_column(f, "sample")
        datasets_upd = read_obs_column(f, "dataset")
        conditions_upd = read_obs_column(f, "condition")
        preps_upd = read_obs_column(f, "preparation_method")
        harms_upd = read_obs_column(f, "condition_harmonized")

    # Build per-cell dataframe (just metadata, no expression)
    cell_df = pd.DataFrame({
        "sample": samples_upd,
        "dataset": datasets_upd,
        "condition": conditions_upd,
        "preparation_method": preps_upd,
        "condition_harmonized": harms_upd,
        "cell_type": cell_types,
    })

    # Compute proportions per sample
    ct_counts = cell_df.groupby(["sample", "cell_type"]).size().unstack(fill_value=0)
    ct_props = ct_counts.div(ct_counts.sum(axis=1), axis=0)

    # Get sample-level metadata (one row per sample, majority vote)
    sample_meta = cell_df.groupby("sample").agg({
        "dataset": "first",
        "condition": lambda x: x.mode().iloc[0],
        "preparation_method": "first",
        "condition_harmonized": lambda x: x.mode().iloc[0],
    }).reset_index()

    props_out = sample_meta.merge(ct_props, left_on="sample", right_index=True)
    props_out.to_csv(PROPS_FILE, index=False)
    print(f"  Saved: {PROPS_FILE} ({len(props_out)} samples × {ct_props.shape[1]} cell types)")

    # --- Summary ---
    print("\n=== SUMMARY ===")
    print(f"Total cells: {n_cells:,}")
    print(f"Total samples: {len(props_out)}")
    print(f"\nPreparation method distribution (samples):")
    for k, v in sample_meta["preparation_method"].value_counts().items():
        print(f"  {k}: {v}")
    print(f"\nCondition harmonized distribution (samples):")
    for k, v in sample_meta["condition_harmonized"].value_counts().items():
        print(f"  {k}: {v}")

    # Verify: no remaining "Mixed" in condition
    remaining_mixed = (cell_df["condition"] == "Mixed").sum()
    print(f"\nRemaining 'Mixed' in condition: {remaining_mixed}")
    if remaining_mixed > 0:
        print("  WARNING: Some 'Mixed' cells remain!")
    else:
        print("  OK: All 'Mixed' labels fixed.")

    # Verify: unsorted samples for DE
    de_eligible = sample_meta[
        sample_meta["preparation_method"].isin(["nuclei", "cd45_negative", "unsorted"])
    ]
    print(f"\nDE-eligible samples (unsorted/nuclei/CD45-): {len(de_eligible)}")
    for k, v in de_eligible["condition_harmonized"].value_counts().items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
