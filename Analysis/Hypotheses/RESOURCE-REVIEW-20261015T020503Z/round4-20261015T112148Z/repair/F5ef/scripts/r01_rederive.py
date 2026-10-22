#!/usr/bin/env python3
"""F5ef repair: re-derive every value the blind checker disputed or flagged, from release files.

Light pandas on saved tables only (login node allowed by the standing rules). Reads:
  * the 2026-08-07 v2 program registry (the stage fit that fixed the Fig 5F pair) and the
    2026-10-14 source-diagnosis registry;
  * the adopted spatial full-composition release (Fig 5F values) and review job 1731 unit tables;
  * the r8 Catalog program-context rows;
  * the COLOC r2 lineage-context primary pairs (RORA in T/NK);
  * the first-round protein-histology table (seed index, not-scored rows, robust names).
Nothing is written except stdout.
"""
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
REV = P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z"
pd.set_option("display.width", 250)

print("== 1. Stage selection of the Fig 5F pair (R06-02) ==")
REG0807 = P / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv"
REG1014 = P / "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv"
cols = ["cell_type", "module", "module_name", "primary_qvalue", "primary_hc3_qvalue", "stability_median",
        "primary_selected", "robust_display"]
for tag, path in (("2026-08-07", REG0807), ("2026-10-14", REG1014)):
    reg = pd.read_csv(path, sep="\t")
    print(tag, path.name, "rows", len(reg), "primary_selected", int(reg.primary_selected.sum()),
          "robust_display", int(reg.robust_display.sum()))
    print(reg[(reg.cell_type == "hepatocytes") & reg.module.isin([8, 20])][cols].to_string(index=False))
h59 = (P / "Analysis/Hypotheses/H59-20261015T112432Z/RESULT_H59.md").read_text()
print("H59:", [l for l in h59.splitlines() if l.startswith("**Verdict")][0])
h90 = (P / "Analysis/Hypotheses/H90-20261015T120337Z/RESULT_H90.md").read_text()
print("H90:", [l for l in h90.splitlines() if l.startswith("**Primary verdict")][0])

print("\n== 2. Fig 5F release values (all 16 abundances) ==")
REL = P / "Analysis/Multimodal_Program_Projection/candidates/spatial-full-composition-sensitivity-candidate-2026-08-11-r2"
rp = pd.read_csv(REL / "full_composition_program_results.tsv", sep="\t")
for r in rp.itertuples(index=False):
    print(f"{r.program_name:24s} {r.dataset:14s} I {r.residual_moran_i:.7f} null_mean {r.residual_null_mean:.7f} "
          f"null_sd {r.residual_null_sd:.7f} (4 sf {r.residual_null_sd:.4g}) excess {r.residual_moran_i - r.residual_null_mean:.7f} "
          f"P {r.residual_pvalue} q {r.residual_padj} state {r.evidence_state}")
nd = pd.read_csv(REL / "full_composition_null_draws.tsv.gz", sep="\t")
for (ds, pid), g in nd.groupby(["dataset", "program_id"]):
    obs = rp[(rp.dataset == ds) & (rp.program_id == pid)].residual_moran_i.iloc[0]
    k = int((g.residual_moran_i >= obs).sum())
    print(f"recount {ds} {pid}: n {len(g)} k {k} P {(k + 1) / (len(g) + 1):.4f} sd(ddof=1) {g.residual_moran_i.std(ddof=1):.7f}")

print("\n== 3. Unit level (release per-section I; job 1731 per-donor excess) ==")
ru = pd.read_csv(REL / "full_composition_unit_results.tsv", sep="\t")
name = dict(zip(rp.program_id, rp.program_name.str.split(" ").str[0]))
ru["prog"] = ru.program_id.map(name)
don = ru.groupby(["dataset", "prog", "individual"]).residual_moran_i.mean().unstack("prog")
print(don.round(4).to_string())
print("units with residual I > 0:\n", (don > 0).groupby(level=0).sum().to_string())
u = pd.read_csv(REV / "jobs/1731/out/unit_level.tsv", sep="\t")


def signflip(x):
    x = np.asarray(x, float)
    obs = x.sum()
    hits = sum(1 for s in itertools.product([1, -1], repeat=len(x)) if (np.array(s) * np.abs(x)).sum() >= obs - 1e-15)
    return hits / 2 ** len(x)


d = u[(u.design == "all16") & (u.unit_type == "donor_or_array")]
w = d.pivot_table(index=["dataset", "unit_id"], columns="program_name", values="excess")
w.columns = [c.split(" ")[0] for c in w.columns]
w["diff"] = w["Stromal"] - w["Ductular"]
print(w.round(4).to_string())
for ds, g in w.groupby(level=0):
    print(f"{ds}: Stromal ECM excess > 0 in {(g.Stromal > 0).sum()}/{len(g)} (exact one-sided sign-flip P {signflip(g.Stromal):.4f}); "
          f"difference > 0 in {(g['diff'] > 0).sum()}/{len(g)} (P {signflip(g['diff']):.4f}); floor {1 / 2 ** len(g):.6f}")
s = u[(u.design == "all16") & (u.unit_type == "section") & (u.dataset == "GSE192741")]
print(s[["program_name", "unit_id", "individual", "residual_moran_i", "null_mean", "excess", "k_null_ge_obs"]].round(4).to_string(index=False))
lin = u[(u.design == "lineage_only") & (u.unit_type == "section") & (u.dataset == "GSE192741") & (u.individual == "H35")]
print("H35 sections, hepatocyte-only design:", lin[["program_name", "unit_id", "residual_moran_i"]].round(4).values.tolist())
src = (REV / "jobs/1731/rerun_spatial.py").read_text().splitlines()
print("lineage-only design (rerun_spatial.py):", [l.strip() for l in src if "lineage_only :" in l or 'factor_names.index("Hepatocytes")' in l])

print("\n== 4. r8 Catalog program-context rows for the spatial results ==")
pc = pd.read_parquet(P / "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_program_context.parquet")
print(pc[pc.dataset_id.isin(["GSE192741", "Vu_et_al_2025"])][["program_label", "dataset_id", "estimate", "q_value", "call_state", "biological_unit"]].to_string(index=False))

print("\n== 5. RORA lineage context in T/NK (5E legend shading) ==")
pairs = pd.read_csv(P / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv", sep="\t")
rr = pairs[(pairs.gene == "RORA") & (pairs.lineage == "t_nk")]
print(rr[["gwas_name", "trait_class", "any_accessible_mass", "evidence_state"]].round(3).to_string(index=False))
WS = REV / "round4-20261015T112148Z/ws/F5ef/out"
g117 = pd.read_csv(WS / "genetic_context_117.tsv", sep="\t")
print(g117[g117.open_genes.fillna("").str.contains("RORA")][["cell_type", "module", "program_name", "n_open_genes", "n_open_direct_genes", "n_open_enzyme_genes"]].to_string(index=False))
print("programs with >=1 open gene", int((g117.n_open_genes > 0).sum()), "; with a direct open gene", int((g117.n_open_direct_genes > 0).sum()))

print("\n== 6. Protein-histology (seed index, not scored, robust names) ==")
h = pd.read_csv(WS / "protein_histology_117.tsv", sep="\t")
t = h[h.testable].reset_index(drop=True)
print("testable", len(t), "; perm seed - 20260923 equals 1-based position among testable programs in registry order:",
      bool(((t.perm_seed - 20260923) == (t.index + 1)).all()),
      "; boot seed - 30260923 equals the same position:", bool(((t.boot_seed - 30260923) == (t.index + 1)).all()))
c = pd.read_csv(P / "figures/candidates/review-r3-fixes-20261015T112148Z/F5ef/fig5e_program_context_117_candidate.tsv", sep="\t")
dd = c[c.display]
print("displayed", len(dd), "; not scored among displayed:")
print(dd[~dd.testable][["cell_type", "module", "program_name", "n_measured", "retained_l1_weight"]].to_string(index=False))
print("robust", int(h.robust.sum()), "negative", int((h.robust & (h.rho < 0)).sum()), "positive", int((h.robust & (h.rho > 0)).sum()))
print(h[h.module_name.str.contains("Amino")][["cell_type", "module", "module_name", "rho", "padj", "robust"]].to_string(index=False))
print(h[h.robust & (h.rho > 0)][["cell_type", "module", "module_name", "rho", "padj"]].sort_values("rho").to_string(index=False))
both = int(((c.n_open_genes > 0) & c.robust).sum())
print("display both genetic and robust histology:", both)
