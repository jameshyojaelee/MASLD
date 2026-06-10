#!/usr/bin/env python
# ─────────────────────────────────────────────────────────────────────────────
# Tier 1C — single-cell mechanistic anchor (integration of EXISTING assets;
# no GPU retraining). Three deliverables:
#   A. Fibroblast meta-state taxonomy: collapse Hotspot fibroblast modules into
#      coherent meta-states (co-variation across donors) + stage trajectory + markers.
#   B. Cross-cell-type disease master-regulator ranking (decoupleR TF activity):
#      rank TFs by breadth (n cell types w/ significant, same-direction disease
#      activity shift) — the anchor to rival Tzouanas SOX4/RELB & Li MITF.
#   C. 1E cross-check: which cell types carry the regression "OFF/ON" program.
#
# Inputs (Analysis/SingleCell/results_gpu_v2/):
#   hotspot_modules/fibroblasts/{donor_scores.tsv,module_genes.tsv}
#   hotspot_modules/crn_transition_scores.tsv
#   fig2_data/tf_activity_per_celltype_condition.csv
#   ../../RNA-seq/results/reversal/de_regression_specific.csv
# Outputs: results_gpu_v2/master_anchor/
# ─────────────────────────────────────────────────────────────────────────────
import os, numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SC   = f"{ROOT}/Analysis/SingleCell/results_gpu_v2"
OUT  = f"{SC}/master_anchor"; os.makedirs(OUT, exist_ok=True)
FIB  = f"{SC}/hotspot_modules/fibroblasts"

# ============================ A. fibroblast meta-states =======================
ds = pd.read_csv(f"{FIB}/donor_scores.tsv", sep="\t")           # sample, module, score
mat = ds.pivot_table(index="module", columns="sample", values="score")
corr = mat.T.corr(method="spearman")                            # module x module
nmod = corr.shape[0]
k = min(5, max(3, nmod // 6))
dist = 1 - corr.values; np.fill_diagonal(dist, 0); dist = (dist + dist.T) / 2
Z = linkage(squareform(dist, checks=False), method="average")
lab = fcluster(Z, t=k, criterion="maxclust")
meta = pd.DataFrame({"module": corr.index, "metastate": lab})

# stage trajectory per module (fibroblasts only)
crn = pd.read_csv(f"{SC}/hotspot_modules/crn_transition_scores.tsv", sep="\t")
crn = crn[crn.cell_type == "fibroblasts"]
traj = crn.pivot_table(index="module", columns="transition", values="delta")
torder = [c for c in ["F0_to_F1","F1_to_F2","F2_to_F3","F3_to_F4"] if c in traj.columns]
traj = traj[torder]
meta = meta.merge(traj.reset_index(), on="module", how="left")

# top genes per module -> per meta-state
mg = pd.read_csv(f"{FIB}/module_genes.tsv", sep="\t")
top_per_mod = (mg.sort_values("weight", ascending=False)
                 .groupby("module").head(12).groupby("module")["gene"].apply(list))
meta["top_genes"] = meta["module"].map(lambda m: ",".join(top_per_mod.get(m, [])[:12]))

# meta-state summary: mean trajectory + pooled markers + monotonic stage rho
def stage_rho(row):
    v = row[torder].values.astype(float)
    if np.all(np.isnan(v)): return np.nan
    x = np.arange(len(v)); m = ~np.isnan(v)
    return np.corrcoef(x[m], v[m])[0,1] if m.sum() > 2 else np.nan
meta["traj_rho"] = meta.apply(stage_rho, axis=1)
summ = meta.groupby("metastate").agg(
    n_modules=("module","size"),
    mean_traj_rho=("traj_rho","mean"),
    **{t:(t,"mean") for t in torder})
ms_markers = (meta.groupby("metastate")["top_genes"]
                  .apply(lambda s: ",".join(pd.Series(",".join(s).split(","))
                  .value_counts().head(15).index)))
summ["top_markers"] = ms_markers
meta.to_csv(f"{OUT}/fibroblast_module_metastates.csv", index=False)
summ.reset_index().to_csv(f"{OUT}/fibroblast_metastate_summary.csv", index=False)
print(f"[A] fibroblast modules={nmod} -> {k} meta-states")
print(summ[["n_modules","mean_traj_rho","top_markers"]].to_string())

# ============================ B. master-regulator ranking =====================
tf = pd.read_csv(f"{SC}/fig2_data/tf_activity_per_celltype_condition.csv")
tf["sig"] = tf.padj < 0.05
tf["dir"] = np.sign(tf.activity_diff)
g = tf.groupby("TF")
rank = pd.DataFrame({
    "n_celltypes_tested": g.size(),
    "n_sig": g.apply(lambda d: int(d.sig.sum())),
    "n_sig_up": g.apply(lambda d: int(((d.sig) & (d.dir > 0)).sum())),
    "n_sig_dn": g.apply(lambda d: int(((d.sig) & (d.dir < 0)).sum())),
    "mean_activity_diff": g.activity_diff.mean(),
    "max_abs_activity_diff": g.activity_diff.apply(lambda x: x.abs().max()),
    "celltypes_sig": g.apply(lambda d: ";".join(sorted(d.loc[d.sig,"cell_type"]))),
}).reset_index()
# directional consistency (breadth in dominant direction)
rank["breadth"] = rank[["n_sig_up","n_sig_dn"]].max(axis=1)
rank["consistent_dir"] = np.where(rank.n_sig_up >= rank.n_sig_dn, "up", "down")
# cross-ref bulk DEG (is the TF itself disease-dysregulated?)
try:
    deg = pd.read_csv(f"{ROOT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
    degsig = deg[(deg.padj < 0.05) & (deg.logFC.abs() > 0.5)]
    rank["bulk_DEG_logFC"] = rank.TF.map(dict(zip(degsig.symbol, degsig.logFC)))
except Exception as e:
    rank["bulk_DEG_logFC"] = np.nan
# known liver / fibrosis / MASLD TFs (literature-curated) for biological grounding
KNOWN = set("""HNF4A HNF1A HNF1B RORA RORC THRB PPARA PPARG PPARD NR1H4 NR1H3 NR1I2 NR1I3
SREBF1 SREBF2 MLXIPL CEBPA CEBPB CEBPD NFKB1 RELA RELB REL STAT1 STAT3 STAT6 SMAD2 SMAD3 SMAD4
JUN JUNB FOS FOSB EGR1 ATF3 ATF4 ATF6 XBP1 DDIT3 FOXO1 FOXO3 FOXA1 FOXA2 ESRRA ESRRG NR5A2
ONECUT1 ONECUT2 GLI1 GLI2 GLI3 TCF7L2 RUNX1 RUNX2 CREB3L1 PRRX1 PRRX2 TWIST1 SNAI1 SNAI2 ZEB1
TEAD1 TEAD4 YAP1 WWTR1 SOX4 SOX9 SPI1 IRF1 IRF7 IRF8 MAF MAFB BHLHE40 KLF6 KLF15 NR3C1 AHR
GATA4 GATA6 MEF2A MEF2C BCL6 ETS1 ELK3 NFE2L2 HIF1A EPAS1""".split())
rank["known_liver_TF"] = rank.TF.isin(KNOWN)
# ubiquity artifact flag: significant in (almost) every cell type => likely global/technical
rank["likely_global_artifact"] = rank.breadth >= 10
rank = rank.sort_values(["breadth","max_abs_activity_diff"], ascending=False)
rank.to_csv(f"{OUT}/master_regulator_ranking.csv", index=False)

# credible master regulators: grounded in biology, NOT ubiquitous
cred = rank[((rank.bulk_DEG_logFC.notna()) | (rank.known_liver_TF)) &
            (rank.breadth.between(3, 9))].copy()
cred = cred.sort_values(["breadth", "max_abs_activity_diff"], ascending=False)
cred.to_csv(f"{OUT}/master_regulator_credible.csv", index=False)
print(f"\n[B] {int(rank.likely_global_artifact.sum())} TFs flagged ubiquitous (>=10/11 cell types) = likely global artifact (excluded).")
print("[B] Top credible cross-cell-type disease master-regulators (bulk-DEG or known-liver-TF, breadth 3-9):")
print(cred.head(20)[["TF","breadth","consistent_dir","n_sig","mean_activity_diff","bulk_DEG_logFC","known_liver_TF","celltypes_sig"]].to_string(index=False))

# ============================ C. 1E reversal cross-check ======================
# Which cell types' disease activity does the regression program reverse?
# Use the master-regulator table: TFs DOWN in regression_specific that are UP in disease.
try:
    rs = pd.read_csv(f"{ROOT}/RNA-seq/results/reversal/de_regression_specific.csv")
    rs_dn = set(rs.loc[(rs["P.Value"] < 0.05) & (rs.logFC < 0), "symbol"].dropna())
    rs_up = set(rs.loc[(rs["P.Value"] < 0.05) & (rs.logFC > 0), "symbol"].dropna())
    # disease-UP TFs that are regression-DOWN (resolution), by cell type
    res = []
    for ct, d in tf.groupby("cell_type"):
        up_dis = set(d.loc[(d.sig) & (d.activity_diff > 0), "TF"])
        dn_dis = set(d.loc[(d.sig) & (d.activity_diff < 0), "TF"])
        res.append(dict(cell_type=ct,
                        n_disease_up_TF=len(up_dis),
                        resolved_by_regression=len(up_dis & rs_dn),
                        regression_restores_down=len(dn_dis & rs_up)))
    cc = pd.DataFrame(res).sort_values("resolved_by_regression", ascending=False)
    cc.to_csv(f"{OUT}/reversal_celltype_crosscheck.csv", index=False)
    print("\n[C] Cell types whose disease-up regulators are reversed by regression:")
    print(cc.to_string(index=False))
except Exception as e:
    print(f"\n[C] cross-check skipped: {e}")
print(f"\nWrote outputs to {OUT}")
