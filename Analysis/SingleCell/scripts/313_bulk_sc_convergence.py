#!/usr/bin/env python3
"""
313_bulk_sc_convergence.py

Compute scRNA hepatocyte meta-subtype markers and project bulk NMF
signatures (S1 Metabolic / S2 Fibrogenic) onto every hepatocyte via
AUCell. Produces the data tables for figS_bulk_sc_convergence.R.

Outputs (in crossmodal/bulk_sc_convergence/):
    meta_subtype_markers.csv                         VIZ_ONLY (do not publish)
                                                     rank_genes_groups wilcoxon per meta-subtype
    aucell_hepatocyte_scores.csv                     cell_id, S1_AUC, S2_AUC, meta_subtype, UMAP coords
    run_summary.txt                                  mapping stats + validation
    pseudobulk_meta_subtype/sample_metadata.csv
    pseudobulk_meta_subtype/cells_per_donor_meta_subtype.csv
    pseudobulk_meta_subtype/subtype_{META}_pseudobulk.csv   raw-count input for limma-voom
    pseudobulk_meta_subtype/subtype_{META}_de.csv           PUBLISHED markers (written by
                                                             313b_pseudobulk_meta_subtype_markers.R)

Pachter P0 retrofit (Squair et al. 2021, Nat Commun): the rank_genes_groups
wilcoxon call at L122 is n=cells (5–10× FPR vs pseudobulk). It is retained
for VIZ_ONLY label/UMAP overlays but its output (meta_subtype_markers.csv)
MUST NOT be used as a publishable marker table. Downstream consumers
(figures, atlas integration) read pseudobulk_meta_subtype/subtype_{META}_de.csv.

Environment: spatial (scanpy, anndata, scipy, numpy, pandas, decoupler)
"""

import os
import sys
import re
import logging
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SUB_DIR = f"{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes"
# Use full-gene 100K subset (all ~37K genes) — the HVG-restricted atlas only has
# 2,931 genes and fails to intersect enough bulk markers for AUCell.
H5AD_PATH = f"{PROJECT}/Analysis/SingleCell/results_gpu_v2/pseudotime/Hepatocytes_subset.h5ad"
META_MAP_PATH = f"{SUB_DIR}/meta_subtype_mapping.csv"
# Fallback: subtype metadata (needed if subset lacks hepatocyte_subtype)
SUBTYPE_META_PATH = f"{SUB_DIR}/hepatocyte_subtype_metadata.csv"

BULK_MARKERS_PATH = f"{PROJECT}/RNA-seq/results/subtypes/subtype_markers.csv"
ATLAS_PATH = f"{PROJECT}/RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"

OUT_DIR = f"{SUB_DIR}/crossmodal/bulk_sc_convergence"
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

summary_lines: list[str] = []


def add_summary(line: str) -> None:
    log.info(line)
    summary_lines.append(line)


def write_pseudobulk_inputs(adata, group_col, groups, out_dir):
    """Aggregate raw counts per (sample, group_col) and write a
    subtype_{GROUP}_pseudobulk.csv (genes × samples) + sample_metadata.csv
    for the downstream R limma-voom step (313b_pseudobulk_meta_subtype_markers.R).

    Squair et al. 2021 (Nat Commun) → marker discovery on n=cells is
    anti-conservative (5–10× FPR). This builds the n=donors design.

    Raw counts are pulled from `adata.layers["counts"]` if present, else
    `adata.raw.X`. adata.X may be log-normalized so it MUST NOT be used.
    """
    os.makedirs(out_dir, exist_ok=True)

    if "counts" in adata.layers:
        X = adata.layers["counts"]
        var_names = adata.var_names
        log.info("Pseudobulk source: adata.layers['counts'] (%s)", type(X).__name__)
    elif adata.raw is not None:
        # Subset raw to current var_names so gene order/intersection matches adata.X
        common = adata.var_names.intersection(adata.raw.var_names)
        X = adata.raw[:, common].X
        var_names = common
        log.info("Pseudobulk source: adata.raw.X restricted to %d shared genes", len(common))
    else:
        log.warning("No .layers['counts'] or .raw — falling back to adata.X. "
                    "Pseudobulk counts WILL BE INCORRECT if X is normalized.")
        X = adata.X
        var_names = adata.var_names

    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    else:
        X = X.tocsr()

    samples = adata.obs["sample"].astype(str).values
    group_vec = adata.obs[group_col].astype(str).values
    datasets = (
        adata.obs["dataset"].astype(str).values
        if "dataset" in adata.obs.columns
        else np.array(["unknown"] * len(samples))
    )
    if "condition_binary" in adata.obs.columns:
        conditions = adata.obs["condition_binary"].astype(str).values
    elif "condition" in adata.obs.columns:
        conditions = adata.obs["condition"].astype(str).values
    else:
        conditions = np.array(["unknown"] * len(samples))

    # Sample-level metadata (one row per donor/sample).
    sample_meta_rows = {}
    for s, ds, cond in zip(samples, datasets, conditions):
        if s not in sample_meta_rows:
            sample_meta_rows[s] = (ds, cond)
    sample_meta_df = pd.DataFrame(
        [(s, ds, cond) for s, (ds, cond) in sample_meta_rows.items()],
        columns=["sample", "dataset", "condition"],
    )
    sample_meta_df.to_csv(os.path.join(out_dir, "sample_metadata.csv"), index=False)
    log.info("  Wrote sample_metadata.csv (%d samples)", len(sample_meta_df))

    # Cell-count contribution table (donor × group) for QC.
    cell_counts = (
        pd.DataFrame({"sample": samples, group_col: group_vec})
        .groupby(["sample", group_col]).size()
        .unstack(fill_value=0)
    )
    cell_counts.to_csv(
        os.path.join(out_dir, f"cells_per_donor_{group_col}.csv")
    )
    log.info("  Wrote cells_per_donor_%s.csv (%s)", group_col, cell_counts.shape)

    # Per group: aggregate raw counts across cells per donor → genes × donors.
    for g in groups:
        g_mask = group_vec == g
        if g_mask.sum() < 100:
            log.warning("  Group %s: only %d cells — skipping pseudobulk", g, g_mask.sum())
            continue

        g_samples = samples[g_mask]
        unique_donors = sorted(set(g_samples))

        donor_to_idx = {d: i for i, d in enumerate(unique_donors)}
        row_ind = np.array([donor_to_idx[s] for s in g_samples])
        col_ind = np.where(g_mask)[0]
        data = np.ones(g_mask.sum(), dtype=np.float32)
        M = sparse.csr_matrix(
            (data, (row_ind, np.arange(len(col_ind)))),
            shape=(len(unique_donors), len(col_ind)),
        )
        donor_gene = M @ X[col_ind, :]
        donor_gene = np.asarray(donor_gene.todense())
        donor_gene = np.rint(donor_gene).astype(np.int64)

        # Sanitize group label for filename (e.g., "Disease-Associated" → "Disease-Associated")
        safe_g = re.sub(r"[^A-Za-z0-9_\-]+", "_", str(g))
        out_df = pd.DataFrame(
            donor_gene.T,
            index=list(var_names),
            columns=unique_donors,
        )
        out_df.index.name = "gene"
        out_file = os.path.join(out_dir, f"subtype_{safe_g}_pseudobulk.csv")
        out_df.to_csv(out_file)
        log.info("  Group %s: %d cells → %d donors × %d genes → %s",
                 g, g_mask.sum(), len(unique_donors), len(var_names), out_file)


# ---------------------------------------------------------------------------
# 1. Load h5ad + attach meta_subtype
# ---------------------------------------------------------------------------
add_summary(f"Loading {H5AD_PATH}")
adata = sc.read_h5ad(H5AD_PATH)
add_summary(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

meta_map = pd.read_csv(META_MAP_PATH)
subtype_to_meta = dict(
    zip(meta_map["subtype"].astype(int), meta_map["meta_subtype"])
)

# If h5ad lacks hepatocyte_subtype, join from subtype metadata CSV
if "hepatocyte_subtype" not in adata.obs.columns:
    add_summary("  h5ad lacks hepatocyte_subtype — joining from subtype_metadata.csv")
    sm = pd.read_csv(SUBTYPE_META_PATH, index_col=0)
    # Keep only hepatocyte_subtype column
    sm = sm[["hepatocyte_subtype"]]
    adata.obs = adata.obs.join(sm, how="left")

adata.obs["meta_subtype"] = (
    pd.to_numeric(adata.obs["hepatocyte_subtype"], errors="coerce")
    .astype("Int64")
    .map(subtype_to_meta)
)

n_labelled = adata.obs["meta_subtype"].notna().sum()
add_summary(f"  labelled with meta_subtype: {n_labelled:,} / {adata.n_obs:,}")
add_summary("Meta-subtype distribution:")
for name, count in (
    adata.obs["meta_subtype"].value_counts().items()
):
    add_summary(f"    {name}: {count:,}")


# ---------------------------------------------------------------------------
# 2. Meta-subtype markers via rank_genes_groups (wilcoxon) — VIZ_ONLY
# ---------------------------------------------------------------------------
# VIZ_ONLY: rank_genes_groups is n=cells (Squair et al. 2021, Nat Commun:
# 5–10× FPR vs pseudobulk). It is retained here for provisional label
# overlays / UMAP coloring / dotplots ONLY. The PUBLISHED meta-subtype
# markers come from the pseudobulk + limma-voom retrofit at Section 2b
# below (pseudobulk_meta_subtype/subtype_{META}_de.csv). Downstream
# consumers (figS_bulk_sc_convergence, atlas integration) MUST read
# pseudobulk_meta_subtype/, NOT meta_subtype_markers.csv.
add_summary("")
add_summary("Computing meta-subtype markers (wilcoxon) — VIZ_ONLY for label assignment...")

# Subset to labelled cells (drop NaN) and work on log-normalized layer
adata_lab = adata[adata.obs["meta_subtype"].notna()].copy()
adata_lab.obs["meta_subtype"] = adata_lab.obs["meta_subtype"].astype(str)

# Check normalization state: use adata.X if it looks log-normalized,
# else normalize a working copy.
X_sample = adata_lab.X[:1000].toarray() if sparse.issparse(adata_lab.X) else adata_lab.X[:1000]
max_val = float(np.max(X_sample))
if max_val > 50:
    add_summary(f"  X looks raw (max={max_val:.1f}); normalizing...")
    sc.pp.normalize_total(adata_lab, target_sum=1e4)
    sc.pp.log1p(adata_lab)
else:
    add_summary(f"  X looks log-normalized (max={max_val:.2f}); using as-is")

# VIZ_ONLY rank_genes_groups — see note above; NOT a publishable marker table.
sc.tl.rank_genes_groups(
    adata_lab, groupby="meta_subtype", method="wilcoxon", pts=True
)

meta_subtypes = sorted(adata_lab.obs["meta_subtype"].unique())
marker_rows = []
for ms in meta_subtypes:
    df = sc.get.rank_genes_groups_df(adata_lab, group=ms)
    df["meta_subtype"] = ms
    marker_rows.append(df)

markers_df = pd.concat(marker_rows, ignore_index=True)
markers_out = f"{OUT_DIR}/meta_subtype_markers.csv"
markers_df.to_csv(markers_out, index=False)
add_summary(f"  wrote {markers_out} ({len(markers_df):,} rows) — VIZ_ONLY; DO NOT publish (see pseudobulk_meta_subtype/)")

# Report top markers per meta-subtype (padj < 0.05)
sig = markers_df[markers_df["pvals_adj"] < 0.05]
add_summary("Significant marker counts per meta-subtype (pvals_adj<0.05) [VIZ_ONLY]:")
for ms in meta_subtypes:
    add_summary(f"    {ms}: {(sig['meta_subtype']==ms).sum():,}")

# Free memory before AUCell
del adata_lab


# ---------------------------------------------------------------------------
# 2b. Pseudobulk meta-subtype markers (PUBLISHED; Pachter P0 Squair-2021 fix)
# ---------------------------------------------------------------------------
# Donor × meta_subtype raw-count pseudobulk + limma-voom (run separately via
# 313b_pseudobulk_meta_subtype_markers.R). Subtype-vs-rest contrast per
# meta-subtype. Filter padj<0.05, |logFC|>0.5. n=donors, not n=cells.
add_summary("")
add_summary("Writing pseudobulk inputs for limma-voom (Pachter P0 retrofit)...")
pseudobulk_dir = f"{OUT_DIR}/pseudobulk_meta_subtype"
# Subset to labelled cells and use raw counts (adata.raw)
adata_pb = adata[adata.obs["meta_subtype"].notna()].copy()
adata_pb.obs["meta_subtype"] = adata_pb.obs["meta_subtype"].astype(str)
write_pseudobulk_inputs(
    adata_pb,
    group_col="meta_subtype",
    groups=sorted(adata_pb.obs["meta_subtype"].unique()),
    out_dir=pseudobulk_dir,
)
add_summary(
    f"  Pseudobulk inputs written to {pseudobulk_dir}. "
    "Run 313b_pseudobulk_meta_subtype_markers.R to produce PUBLISHED limma-voom DE tables."
)
del adata_pb


# ---------------------------------------------------------------------------
# 3. Build bulk S1/S2 gene sets (ENSG → symbol via multi_evidence atlas)
# ---------------------------------------------------------------------------
add_summary("")
add_summary("Loading bulk NMF markers + gene ID map...")

bulk_markers = pd.read_csv(BULK_MARKERS_PATH)
add_summary(f"  bulk markers: {len(bulk_markers):,} rows")

# Strip ENSG version
bulk_markers["ensg_clean"] = bulk_markers["gene"].str.replace(
    r"\.\d+$", "", regex=True
)

atlas = pd.read_csv(
    ATLAS_PATH, usecols=["human_symbol", "ensembl_id"]
).dropna().drop_duplicates("ensembl_id")
atlas["ensg_clean"] = atlas["ensembl_id"].str.replace(r"\.\d+$", "", regex=True)
ensg_to_sym = dict(zip(atlas["ensg_clean"], atlas["human_symbol"]))

bulk_markers["symbol"] = bulk_markers["ensg_clean"].map(ensg_to_sym)
mapping_rate = bulk_markers["symbol"].notna().mean()
add_summary(f"  ENSG → symbol mapping rate: {mapping_rate:.1%}")

sc_genes = set(adata.var_names)
bulk_markers["in_sc"] = bulk_markers["symbol"].isin(sc_genes)
add_summary(
    f"  bulk markers present in sc var_names: {bulk_markers['in_sc'].mean():.1%}"
)

# Top 50 up-markers per S1, S2 by limma t-stat
top_n = 50
gene_sets: dict[str, list[str]] = {}
for subtype in ["S1", "S2"]:
    sub = bulk_markers[
        (bulk_markers["subtype"] == subtype)
        & (bulk_markers["direction"] == "up")
        & (bulk_markers["in_sc"])
    ].copy()
    sub = sub.sort_values("t", ascending=False).head(top_n)
    gene_sets[subtype] = sub["symbol"].dropna().unique().tolist()
    add_summary(
        f"  {subtype} gene set: {len(gene_sets[subtype])} symbols"
    )


# ---------------------------------------------------------------------------
# 4. AUCell scoring — use decoupler if available, else sc.tl.score_genes
# ---------------------------------------------------------------------------
add_summary("")
add_summary("Scoring bulk S1/S2 signatures on every hepatocyte...")

aucell_method = None
try:
    import decoupler as dc

    # decoupler's AUCell wants long-form gene set dataframe
    gs_df = pd.concat(
        [
            pd.DataFrame({"source": "S1_bulk", "target": gene_sets["S1"]}),
            pd.DataFrame({"source": "S2_bulk", "target": gene_sets["S2"]}),
        ]
    )
    try:
        dc.run_aucell(
            mat=adata,
            net=gs_df,
            source="source",
            target="target",
            min_n=5,
            use_raw=False,
            verbose=True,
        )
        est = adata.obsm["aucell_estimate"]
        s1_score = est["S1_bulk"].values
        s2_score = est["S2_bulk"].values
        aucell_method = "decoupler.run_aucell"
    except AttributeError:
        # Newer decoupler API: dc.mt.aucell
        result = dc.mt.aucell(
            data=adata, net=gs_df, tmin=5, verbose=True
        )
        est = adata.obsm["score_aucell"]
        s1_score = est["S1_bulk"].values
        s2_score = est["S2_bulk"].values
        aucell_method = "decoupler.mt.aucell"
except Exception as e:
    log.warning(f"  decoupler AUCell unavailable ({e}); falling back to score_genes")
    sc.tl.score_genes(
        adata, gene_list=gene_sets["S1"], score_name="_S1_score", use_raw=False
    )
    sc.tl.score_genes(
        adata, gene_list=gene_sets["S2"], score_name="_S2_score", use_raw=False
    )
    s1_score = adata.obs["_S1_score"].values
    s2_score = adata.obs["_S2_score"].values
    aucell_method = "scanpy.score_genes (fallback)"

add_summary(f"  method: {aucell_method}")
add_summary(f"  S1 score: mean={float(np.nanmean(s1_score)):.3f}, std={float(np.nanstd(s1_score)):.3f}")
add_summary(f"  S2 score: mean={float(np.nanmean(s2_score)):.3f}, std={float(np.nanstd(s2_score)):.3f}")


# ---------------------------------------------------------------------------
# 5. Assemble AUCell output CSV
# ---------------------------------------------------------------------------
umap_arr = adata.obsm.get("X_umap")
out = pd.DataFrame(
    {
        "cell_id": adata.obs_names,
        "meta_subtype": adata.obs["meta_subtype"].astype(str).values,
        "hepatocyte_subtype": adata.obs["hepatocyte_subtype"].astype(str).values,
        "condition": adata.obs["condition"].astype(str).values,
        "S1_AUC": s1_score,
        "S2_AUC": s2_score,
    }
)
if umap_arr is not None:
    out["UMAP_1"] = umap_arr[:, 0]
    out["UMAP_2"] = umap_arr[:, 1]

aucell_out = f"{OUT_DIR}/aucell_hepatocyte_scores.csv"
out.to_csv(aucell_out, index=False)
add_summary(f"  wrote {aucell_out} ({len(out):,} rows)")

# Sanity print: mean scores per meta-subtype
add_summary("")
add_summary("Mean AUC per meta-subtype (validation):")
for ms in sorted(out["meta_subtype"].dropna().unique()):
    sub = out[out["meta_subtype"] == ms]
    add_summary(
        f"    {ms:>22s}  S1={sub['S1_AUC'].mean():.4f}  S2={sub['S2_AUC'].mean():.4f}  n={len(sub):,}"
    )


# ---------------------------------------------------------------------------
# 6. Write run summary
# ---------------------------------------------------------------------------
summary_path = f"{OUT_DIR}/run_summary.txt"
with open(summary_path, "w") as f:
    f.write("\n".join(summary_lines) + "\n")
log.info(f"Wrote run summary to {summary_path}")
log.info("Done.")
