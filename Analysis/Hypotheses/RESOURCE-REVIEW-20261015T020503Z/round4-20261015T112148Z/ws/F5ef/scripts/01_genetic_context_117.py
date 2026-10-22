#!/usr/bin/env python3
"""F5ef: genetic-chromatin context for all 117 frozen programs on the adopted COLOC r2 release.

Replaces the Fig 5E "Open-promoter n" column, which counted program genes whose promoter held a
July fine-mapping credible-set variant (GWAS/finemapping/results/credible_sets.csv, no MVP EUR
studies). The replacement uses the same object and release as Fig 5D: the 767 adopted r2
colocalized gene-study primary signal pairs (28 direct MASLD/liver fat, 739 liver enzyme) and
their normalized SNP.PP.H4 posterior mass in static open chromatin of five liver lineages
(GSE244832 and GSE281367 snATAC).

Per program (v2 membership, gencode_v49_unique_symbol_confirmed symbols) and its own lineage
(hepatocytes->hepatocyte, cholangiocytes->cholangiocyte, macrophages->macrophage,
fibroblasts->stellate, tcells->t_nk):
  n_coloc_genes       member genes with >=1 adopted r2 primary pair (any lineage state)
  n_open_genes        member genes with >=1 pair in state replicated_accessible or source_dependent
                      (>=0.5 posterior mass in peaks of at least one source) in the program lineage
  n_open_direct / n_open_enzyme  the same split by trait class
Descriptive counts; no test and no multiplicity family.
Also recomputes the same counts for the 22 superseded programs (v1 membership) for comparison
with the printed Fig 5E open-promoter column.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
WS = P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef"
OUT = WS / "out"
PAIRS = P / ("Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/"
             "genetics/context/genetic_lineage_context_primary_pairs.tsv")
REG = P / "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv"
MEM = P / ("Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/"
           "hotspot/program_membership_v2.tsv")
V1_MEM = P / "Analysis/Multimodal_Program_Projection/results/frozen_program_membership.tsv"
V1_WIDE = P / "Analysis/Multimodal_Program_Projection/results/program_context_wide.tsv"
LINEAGE = {"hepatocytes": "hepatocyte", "cholangiocytes": "cholangiocyte", "macrophages": "macrophage",
           "fibroblasts": "stellate", "tcells": "t_nk"}
OPEN = {"replicated_accessible", "source_dependent"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


pairs = pd.read_csv(PAIRS, sep="\t")
assert len(pairs) == 3835 and pairs["primary_signal_pair"].all()
pairs["pair_id"] = pairs["gwas_name"] + "|" + pairs["gene"]
assert pairs["pair_id"].nunique() == 767
assert pairs.drop_duplicates("pair_id")["trait_class"].value_counts().to_dict() == {"liver_enzyme": 739, "direct_MASLD": 28}

reg = pd.read_csv(REG, sep="\t", usecols=["cell_type", "module", "program_uid", "module_name"])
assert len(reg) == 117
mem = pd.read_csv(MEM, sep="\t")
mem = mem[mem["mapped_symbol_status"] == "gencode_v49_unique_symbol_confirmed"]


def count(genes: set[str], lineage: str) -> dict:
    sub = pairs[pairs["gene"].isin(genes)]
    lin = sub[sub["lineage"] == lineage]
    open_ = lin[lin["evidence_state"].isin(OPEN)]
    return {
        "n_member_genes": len(genes),
        "n_coloc_genes": sub["gene"].nunique(),
        "n_coloc_pairs": sub["pair_id"].nunique(),
        "n_open_genes": open_["gene"].nunique(),
        "n_open_pairs": open_["pair_id"].nunique(),
        "n_open_direct_genes": open_.loc[open_["trait_class"] == "direct_MASLD", "gene"].nunique(),
        "n_open_enzyme_genes": open_.loc[open_["trait_class"] == "liver_enzyme", "gene"].nunique(),
        "open_genes": ";".join(sorted(open_["gene"].unique())),
        "open_pairs": ";".join(sorted(open_["pair_id"].unique())),
    }


rows = []
for r in reg.itertuples(index=False):
    genes = set(mem.loc[mem["program_uid"] == r.program_uid, "mapped_symbol"].astype(str).str.strip())
    d = {"program_uid": r.program_uid, "cell_type": r.cell_type, "module": r.module,
         "program_name": r.module_name, "lineage": LINEAGE[r.cell_type]}
    d.update(count(genes, LINEAGE[r.cell_type]))
    rows.append(d)
g117 = pd.DataFrame(rows)
g117["release"] = "COLOC r2 adopted 2026-10-14; 767 primary pairs; atac-context-v3-coloc-r2-20260924T053244Z"
g117["definition"] = ("member genes with >=1 adopted r2 colocalized primary pair whose normalized SNP.PP.H4 "
                      "posterior has >=0.5 mass in static open chromatin of the program lineage in >=1 snATAC source")
g117.to_csv(OUT / "genetic_context_117.tsv", sep="\t", index=False)

# Comparison on the 22 superseded programs (v1 membership), against the printed open-promoter column
v1 = pd.read_csv(V1_MEM, sep="\t")
v1 = v1[v1["mapped_symbol"].fillna(False).astype(bool) & v1["gene_symbol"].notna()]
wide = pd.read_csv(V1_WIDE, sep="\t")
rows = []
for r in wide.itertuples(index=False):
    genes = set(v1.loc[v1["program_id"] == r.program_id, "gene_symbol"].astype(str).str.strip())
    d = {"program_id": r.program_id, "program_name": r.program_name, "lineage": LINEAGE[r.cell_type],
         "printed_open_promoter_enzyme": r.gwas_atac_enzyme_genes,
         "printed_open_promoter_direct": r.gwas_atac_direct_genes}
    d.update(count(genes, LINEAGE[r.cell_type]))
    rows.append(d)
c22 = pd.DataFrame(rows)
c22.to_csv(OUT / "genetic_context_22_comparison.tsv", sep="\t", index=False)

with open(OUT / "genetic_context_inputs.sha256", "w") as fh:
    for f in (PAIRS, REG, MEM, V1_MEM, V1_WIDE):
        fh.write(f"{sha(f)}  {f}\n")

print("117 programs: with >=1 open coloc gene", int((g117["n_open_genes"] > 0).sum()),
      "| direct-trait open genes in any program", int((g117["n_open_direct_genes"] > 0).sum()))
print(g117.loc[g117["n_open_genes"] > 0, ["program_uid", "cell_type", "module", "program_name", "n_open_genes",
                                          "n_open_direct_genes", "n_open_enzyme_genes", "open_genes"]].to_string())
print()
print(c22[["program_id", "program_name", "printed_open_promoter_enzyme", "printed_open_promoter_direct",
           "n_coloc_genes", "n_open_genes", "n_open_direct_genes", "open_genes"]].to_string())
