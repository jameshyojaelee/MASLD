#!/usr/bin/env python
"""
Audit: does Hotspot global__19 (ductular reaction) reflect trans-differentiation,
or is it a composition artifact of cholangiocyte expansion in MASLD?

Step 1: CT x stage mean module-19 score (with cell counts).
Step 2: lmer(donor_mean_19 ~ stage + (1|dataset)) before vs after frac_Cholangiocytes.
Step 3: Hepatocyte-only donor_mean_19 ~ stage + (1|dataset).
"""

import os, sys, json
import numpy as np
import pandas as pd
import anndata as ad

PROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUTDIR = f"{PROOT}/Analysis/SingleCell/scripts/hotspot_modules/audits"
os.makedirs(OUTDIR, exist_ok=True)

CELL_SCORES = f"{PROOT}/Analysis/SingleCell/results_gpu_v2/hotspot_modules/global/cell_scores.parquet"
DONOR_META  = f"{PROOT}/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
ATLAS_H5AD  = f"{PROOT}/Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"

MODULE_ID = 19

print("[1/5] Loading module-19 cell scores ...")
cs = pd.read_parquet(CELL_SCORES)
cs = cs[cs['module'] == MODULE_ID][['cell_id', 'score']].rename(columns={'score': 'module19'})
print(f"  cells with module-19 score: {len(cs):,}")

print("[2/5] Loading atlas obs (cell_type, sample) ...")
a = ad.read_h5ad(ATLAS_H5AD, backed='r')
obs = a.obs[['sample', 'cell_type']].copy()
obs.index.name = 'cell_id'
obs = obs.reset_index()
print(f"  atlas cells: {len(obs):,}; cell_type levels: {obs['cell_type'].nunique()}")
print("  top cell_types:", obs['cell_type'].value_counts().head(10).to_dict())

print("[3/5] Merge scores + obs + donor metadata ...")
df = cs.merge(obs, on='cell_id', how='inner')
print(f"  merged cells: {len(df):,}")

dm = pd.read_csv(DONOR_META, sep='\t')
keep_cols = ['sample','dataset','disease_stage_coarse','disease_stage_numeric',
             'frac_Cholangiocytes','frac_Hepatocytes','frac_Macrophages',
             'frac_Fibroblasts','frac_Endothelial cells']
dm = dm[[c for c in keep_cols if c in dm.columns]]
df = df.merge(dm, on='sample', how='inner')
print(f"  after donor-meta merge: {len(df):,} cells; donors: {df['sample'].nunique()}")

# ------------------------------------------------------------------
# Step 1: CT x stage table
# ------------------------------------------------------------------
print("\n[4/5] STEP 1: CT x disease_stage mean module-19 score")
stage_order = ['Healthy', 'Steatosis', 'Steatohepatitis', 'Cirrhosis']
df['disease_stage_coarse'] = pd.Categorical(df['disease_stage_coarse'],
                                            categories=stage_order, ordered=True)

mean_tbl = df.groupby(['cell_type', 'disease_stage_coarse'], observed=True)['module19'].agg(['mean','count']).reset_index()
mean_pvt = mean_tbl.pivot(index='cell_type', columns='disease_stage_coarse', values='mean')
count_pvt = mean_tbl.pivot(index='cell_type', columns='disease_stage_coarse', values='count')
# Order by mean across all
mean_pvt['overall_mean'] = df.groupby('cell_type', observed=True)['module19'].mean()
mean_pvt = mean_pvt.sort_values('overall_mean', ascending=False)
count_pvt = count_pvt.loc[mean_pvt.index]

print("\n--- Mean module-19 score by cell_type x stage ---")
print(mean_pvt.round(4).to_string())
print("\n--- Cell counts by cell_type x stage ---")
print(count_pvt.fillna(0).astype('Int64').to_string())

mean_pvt.to_csv(f"{OUTDIR}/step1_mean_by_ct_stage.tsv", sep='\t')
count_pvt.to_csv(f"{OUTDIR}/step1_count_by_ct_stage.tsv", sep='\t')

# Hepatocyte-specific monotonic trend check
print("\nHepatocyte trend across stages:")
if 'Hepatocytes' in mean_pvt.index:
    print(mean_pvt.loc['Hepatocytes', stage_order].round(4).to_dict())
if 'Cholangiocytes' in mean_pvt.index:
    print("Cholangiocyte trend across stages:")
    print(mean_pvt.loc['Cholangiocytes', stage_order].round(4).to_dict())

# ------------------------------------------------------------------
# Save per-donor aggregates for R lmer step
# ------------------------------------------------------------------
print("\n[5/5] Build donor-level aggregates for Step 2 / Step 3 ...")

# Step 2: donor-level overall mean of module-19 (across all cells)
donor_overall = (df.groupby(['sample','dataset','disease_stage_numeric','disease_stage_coarse',
                             'frac_Cholangiocytes','frac_Hepatocytes'], observed=True)['module19']
                   .mean().reset_index().rename(columns={'module19':'donor_mean_19'}))

# Step 3: hepatocyte-only donor mean
hep = df[df['cell_type'] == 'Hepatocytes']
donor_hep = (hep.groupby(['sample','dataset','disease_stage_numeric','disease_stage_coarse'], observed=True)['module19']
                .agg(['mean','count']).reset_index()
                .rename(columns={'mean':'donor_mean_19_hep','count':'n_hep_cells'}))

# Also save chol-only as comparison
chol = df[df['cell_type'] == 'Cholangiocytes']
donor_chol = (chol.groupby(['sample','dataset','disease_stage_numeric','disease_stage_coarse'], observed=True)['module19']
                .agg(['mean','count']).reset_index()
                .rename(columns={'mean':'donor_mean_19_chol','count':'n_chol_cells'}))

donor_overall.to_csv(f"{OUTDIR}/donor_overall_module19.tsv", sep='\t', index=False)
donor_hep.to_csv(f"{OUTDIR}/donor_hep_module19.tsv", sep='\t', index=False)
donor_chol.to_csv(f"{OUTDIR}/donor_chol_module19.tsv", sep='\t', index=False)

print(f"\nDonor-overall rows: {len(donor_overall)}; hep rows: {len(donor_hep)}; chol rows: {len(donor_chol)}")
print("Saved donor aggregates. Run audit_global19_lmer.R next.")
print("\n[DONE Python phase]")
