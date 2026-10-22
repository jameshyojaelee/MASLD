"""Independent recomputation of Figure 5E rows, bars and intervals.

Genetic column, own code from:
  ATAC v3 COLOC-r2 replay pairs: Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/
      genetics/context/genetic_lineage_context_primary_pairs.tsv
  117-program registry: figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv
  membership: Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv
  rule: member gene with >=1 primary pair whose normalized posterior has >=0.5 mass in peaks of the program lineage in >=1 source.
Protein column, from the analysis output protein_histology_117.tsv (job 22253064; 9,999 Freedman-Lane permutations, seeds recorded):
  BH recomputed over testable programs; robust = q<0.05 & sensitivity_sign_agree & batch_sign_agree.
"""
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = P / "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5/recompute"
pairs = pd.read_csv(P / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv", sep="\t")
reg = pd.read_csv(P / "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv", sep="\t",
                  usecols=["cell_type", "module", "program_uid", "module_name"])
mem = pd.read_csv(P / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv", sep="\t")
ph = pd.read_csv(P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef/out/protein_histology_117.tsv", sep="\t")
LIN = {"hepatocytes": "hepatocyte", "cholangiocytes": "cholangiocyte", "fibroblasts": "stellate", "macrophages": "macrophage", "tcells": "t_nk"}

pairs = pairs[pairs.primary_signal_pair.astype(str).str.upper() == "TRUE"].copy()
pairs["maxsrc"] = np.maximum(pairs.gse244832_accessible_mass, pairs.gse281367_accessible_mass) / pairs.mapped_posterior_mass
pairs["open"] = pairs.maxsrc >= 0.5
pairs["direct"] = pairs.trait_class.str.contains("direct")
op = pairs[pairs.open].groupby(["lineage", "gene"]).agg(direct=("direct", "max"), enzyme=("direct", lambda s: (~s).any())).reset_index()
coloc_genes = set(pairs.gene)

mem["sym"] = mem.mapped_symbol.where(mem.mapped_symbol.notna() & (mem.mapped_symbol.astype(str) != "False") & (mem.mapped_symbol.astype(str) != "True"), mem.canonical_gene)
# mapped_symbol may be a symbol string; fall back to canonical_gene
mem["sym"] = np.where(mem.mapped_symbol.astype(str).isin(["True", "False", "nan"]), mem.canonical_gene, mem.mapped_symbol)
STRICT = __import__("sys").argv[1:] == ["strict"]
if STRICT:  # panel rule as inferred: drop members whose symbol is flagged unmapped_or_ambiguous
    mem = mem[mem.mapped_symbol_status != "unmapped_or_ambiguous"]
rows = []
for _, r in reg.iterrows():
    lin = LIN[r.cell_type]
    g = set(mem.loc[mem.program_uid == r.program_uid, "sym"].astype(str))
    o = op[(op.lineage == lin) & op.gene.isin(g)]
    rows.append({"program_uid": r.program_uid, "cell_type": r.cell_type, "module": r.module, "name": r.module_name,
                 "n_members": len(g), "n_coloc_genes": len(g & coloc_genes), "n_open": len(o),
                 "n_open_direct": int(o.direct.sum()), "n_open_enzyme_only": int((~o.direct).sum()),
                 "open_genes": ";".join(sorted(o.gene))})
gen = pd.DataFrame(rows)
print("programs", len(gen), "with >=1 open gene", (gen.n_open > 0).sum(), "with a direct-open gene", (gen.n_open_direct > 0).sum())
print(gen[gen.n_open_direct > 0][["cell_type", "module", "name", "open_genes"]].to_string())

# protein column
ph = ph.copy()
t = ph.testable.astype(str).str.upper() == "TRUE"
q = pd.Series(np.nan, index=ph.index)
q[t] = stats.false_discovery_control(ph.loc[t, "pvalue"], method="bh")
ph["q_re"] = q
ph["robust_re"] = (ph.q_re < 0.05) & ph.sensitivity_sign_agree.astype(str).str.upper().eq("TRUE") & ph.batch_sign_agree.astype(str).str.upper().eq("TRUE")
print("testable", int(t.sum()), "max |q - padj|", float((ph.q_re - ph.padj).abs().max()), "robust", int(ph.robust_re.sum()),
      "robust neg/pos", int((ph.robust_re & (ph.rho < 0)).sum()), int((ph.robust_re & (ph.rho > 0)).sum()),
      "q<0.05 regardless of sign checks", int((ph.q_re < 0.05).sum()))
print("n per testable row", ph.loc[t, "n"].value_counts().to_dict(), "perm", ph.n_permutations.unique(), "seeds", ph.perm_seed.unique()[:3], ph.boot_seed.unique()[:3])
m = gen.merge(ph[["program_uid", "n_measured", "retained_l1_weight", "testable", "rho", "ci_low", "ci_high", "pvalue", "q_re", "robust_re",
                  "sensitivity_sign_agree", "batch_sign_agree"]], on="program_uid", how="left")
m["display"] = (m.n_open > 0) | m.robust_re
print("displayed", int(m.display.sum()), "open", int((m.n_open > 0).sum()), "robust", int(m.robust_re.sum()), "both", int(((m.n_open > 0) & m.robust_re).sum()))
print("displayed by cell type", m[m.display].cell_type.value_counts().to_dict())
print("displayed but not scored", m[m.display & ~(m.testable.astype(str).str.upper() == "TRUE")][["cell_type", "module", "name", "n_measured", "retained_l1_weight"]].to_string())
# CI excludes 0 but not robust (drawn grey)
m["ci_excl0"] = (m.ci_low > 0) | (m.ci_high < 0)
print("displayed, gray, CI excludes 0:", m[m.display & ~m.robust_re & m.ci_excl0][["cell_type", "module", "name", "rho", "ci_low", "ci_high", "q_re"]].round(3).to_string())
print("all testable with CI excluding 0:", int((m.ci_excl0 & (m.testable.astype(str).str.upper() == "TRUE")).sum()), "nominal p<0.05:", int((m.pvalue < 0.05).sum()))
# hidden testable programs
hid = m[~m.display & (m.testable.astype(str).str.upper() == "TRUE")]
print("hidden testable programs", len(hid), "of which q<0.05 but not robust", int((hid.q_re < 0.05).sum()))
m.to_csv(OUT / ("5e_recomputed_117_strict.tsv" if STRICT else "5e_recomputed_117.tsv"), sep="\t", index=False)
