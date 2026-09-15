#!/usr/bin/env python3
"""Step 12 (Package G): place each colocalized signal in the existing biological context, link by link.

Link ladder (each link keeps its own evidence state; missing links stay visible):
  L1 genetic signal → colocalized gene (existing COLOC pp_h4; 'genetically anchored candidate')
  L2 lineage context: variant-consistent accessibility per lineage (atac-context-v3 evidence_state);
     lineage EXPRESSION is not derived here (no adopted per-gene lineage expression table)
  L3 membership in a read-only Hotspot program (program_membership_v2)
  L4 spatial coverage of that program (spatial-resource candidate coverage_status)
  L5 protein: translation-build gene_protein_native (log2FC, BH q) and catalog protein_state
Creative comparison: lineages where inherited regulation is chromatin-supported vs the cell types
of the programs the gene belongs to (compartment agreement / mismatch / undetermined).

Output (tables/): context_links.tsv (long, one row per signal × link), context_summary.tsv (one row per signal)
"""

from __future__ import annotations

from collections import defaultdict

import pandas as pd

import lib_atlas as la

TABLES = la.out_root() / "tables"
CTX = la.PROJECT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
PROG = la.PROJECT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
SPATIAL = la.PROJECT / "Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11/coverage/spatial_program_coverage.tsv"
TR = la.PROJECT / "RNA-seq/results/histology_anchored_continuum/translation/hac-translation-20260906T230500Z"
CELLTYPE_TO_LINEAGE = {"hepatocytes": "hepatocyte", "fibroblasts": "stellate", "macrophages": "macrophage", "cholangiocytes": "cholangiocyte", "tcells": "t_nk"}


def main() -> None:
    signals = [s for s in la.read_tsv(TABLES / "eligible_signals.tsv") if s["universe"] != "B_direct"]
    pairs = la.read_tsv(CTX / "genetics/context/genetic_lineage_context_primary_pairs.tsv")
    lineage_state = {(r["gwas_name"], r["ensembl"], r["signal_pair_index"], r["lineage"]): r for r in pairs}
    prog = pd.read_csv(PROG, sep="\t", keep_default_na=False)
    prog_by_gene = defaultdict(list)
    for _, r in prog.iterrows():
        prog_by_gene[r["canonical_gene"]].append((r["program_uid"], r["cell_type"], float(r["original_l1_weight"])))
    spatial = pd.read_csv(SPATIAL, sep="\t", keep_default_na=False)
    spatial_by_prog = defaultdict(list)
    for _, r in spatial.iterrows():
        spatial_by_prog[r["program_uid"]].append(r)
    protein = pd.read_csv(TR / "gene_protein_native.tsv", sep="\t", keep_default_na=False).drop_duplicates("gene_name").set_index("gene_name")
    cat = pd.read_csv(TR / "gene_catalog_translation.tsv.gz", sep="\t", keep_default_na=False).drop_duplicates("ensembl").set_index("ensembl")
    chrom_census = pd.read_csv(TABLES / "chromatin_census.tsv", sep="\t").set_index("signal_uid")

    links, summary = [], []
    for s in signals:
        sig, gene, ens = s["signal_uid"], s["gene"], s["ensembl"]
        base = {"signal_uid": sig, "universe": s["universe"], "gwas_name": s["gwas_name"], "gene": gene, "ensembl": ens}
        links.append({**base, "link": "L1_signal_to_gene", "target": ens, "evidence_state": "supported", "value": s["pp_h4"],
                      "unit": "COLOC PP.H4 (SuSiE signal pair)", "note": "genetically anchored candidate; single colocalized signal; combinations untested"})
        supported_lineages = []
        for lin in ("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"):
            r = lineage_state.get((s["gwas_name"], ens, s["signal_pair_index"], lin))
            st = r["evidence_state"] if r else "missing"
            if st == "replicated_accessible":
                supported_lineages.append(lin)
            links.append({**base, "link": "L2_lineage_chromatin_context", "target": lin, "evidence_state": st,
                          "value": r["joint_any_accessible_mass"] if r else "", "unit": "posterior mass in jointly testable accessible peaks",
                          "note": (r["testability_reason"] if r else "no row") + "; assay-native chromatin context, not expression"})
        links.append({**base, "link": "L2_lineage_expression", "target": "", "evidence_state": "untestable",
                      "value": "", "unit": "", "note": "no adopted per-gene lineage expression table; not derived in this extension"})
        programs = prog_by_gene.get(gene, [])
        prog_cell_types = sorted({CELLTYPE_TO_LINEAGE.get(ct, ct) for _, ct, _ in programs})
        if programs:
            for uid, ct, w in programs:
                links.append({**base, "link": "L3_program_membership", "target": uid, "evidence_state": "supported", "value": w,
                              "unit": "original L1 weight in read-only program", "note": f"cell_type={ct}; membership does not create a variant-to-gene link"})
                cov = spatial_by_prog.get(uid, [])
                if cov:
                    for c in cov:
                        links.append({**base, "link": "L4_program_spatial_coverage", "target": f"{uid}|{c['dataset_id']}", "evidence_state": c["coverage_status"],
                                      "value": c["detection_fraction"], "unit": "detected fraction of program genes",
                                      "note": f"{c['biological_unit']} ({c['biological_unit_resolution']}); {c['testability_reason']}"})
                else:
                    links.append({**base, "link": "L4_program_spatial_coverage", "target": uid, "evidence_state": "missing", "value": "", "unit": "", "note": "no spatial coverage row"})
        else:
            links.append({**base, "link": "L3_program_membership", "target": "", "evidence_state": "not_applicable", "value": "", "unit": "", "note": "gene in no read-only program"})
            links.append({**base, "link": "L4_program_spatial_coverage", "target": "", "evidence_state": "not_applicable", "value": "", "unit": "", "note": "no program"})
        if gene in protein.index:
            p = protein.loc[gene]
            links.append({**base, "link": "L5_protein", "target": gene, "evidence_state": cat.at[ens, "protein_state"] if ens in cat.index else "reported",
                          "value": p["protein_log2FC"], "unit": f"{p['protein_effect_unit']} ({p['protein_assay']}); BH q={p['protein_BH_q']}",
                          "note": f"n={p['protein_n_samples']} (control {p['protein_n_control']}, MASLD {p['protein_n_masld']}); liver tissue proteome; inapplicable to lncRNA products"})
        else:
            links.append({**base, "link": "L5_protein", "target": gene, "evidence_state": "untestable", "value": "", "unit": "", "note": "gene not quantified in the liver proteome table"})
        cc = chrom_census.loc[sig] if sig in chrom_census.index else None
        if not supported_lineages or not prog_cell_types:
            compartment = "undetermined"
        elif set(supported_lineages) & set(prog_cell_types):
            compartment = "agreement"
        else:
            compartment = "mismatch"
        summary.append({**base, "analysis_block": s["analysis_block"], "pp_h4": s["pp_h4"],
                        "chromatin_supported_lineages": ";".join(supported_lineages), "program_cell_types": ";".join(prog_cell_types),
                        "n_programs": len(programs), "compartment_relation": compartment,
                        "remodeling_state": cat.at[ens, "remodeling_state"] if ens in cat.index else "", "protein_state": cat.at[ens, "protein_state"] if ens in cat.index else "",
                        "target_readout_identity": cat.at[ens, "target_readout_identity"] if ens in cat.index else "",
                        "hepatocyte_any_mass": cc["hepatocyte_any_mass"] if cc is not None else "", "stellate_any_mass": cc["stellate_any_mass"] if cc is not None else "",
                        "n_missing_links": sum(1 for l in links if l["signal_uid"] == sig and l["evidence_state"] in ("missing", "untestable", "not_applicable"))})
    la.write_tsv_once(TABLES / "context_links.tsv", links, ["signal_uid", "universe", "gwas_name", "gene", "ensembl", "link", "target", "evidence_state", "value", "unit", "note"])
    la.write_tsv_once(TABLES / "context_summary.tsv", summary, list(summary[0].keys()))
    la.log(f"Package G: {len(summary)} signals, {len(links)} links")


if __name__ == "__main__":
    main()
