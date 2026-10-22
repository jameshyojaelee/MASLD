"""Build work/index_edits.json (I- edits). Each `old` is counted in the live file;
value_reread is joined from work/index_reread.tsv (produced by reread_index.py)."""
import json
import os

import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
WORK = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-QUALITY-20261017T020832Z/text_fixes-20261017T025625Z/work")
H59 = "Analysis/Hypotheses/H59-20261015T112432Z/out/outcome/results_primary_family.tsv"
REP4 = "figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/repair"
R9D = ("Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/"
       "catalog_r9d-20261016T105550Z/build_out")

edits = [
    {
        "id": "I-01",
        "file": "figures/main/INDEX.md",
        "old": "Figure 4C now reflects the corrected source diagnoses: BICC1 alone passes\n"
               "primary BH and no program passes HC3.",
        "new": "Figure 4C now shows the H59 primary fit (60 ambient-correctable donors, a\n"
               "separate residual variance per dataset): 0 of 117 programs pass the\n"
               "within-dataset permutation BH/stability rule (smallest q=0.323). In the\n"
               "supplementary 62-donor equal-variance sensitivity, only Ductular injury\n"
               "(BICC1 member) has q<0.05 (q=0.038), and no program passes HC3.",
        "reason": "The installed 4C is the H59 repair render (sha256 match); its primary result is 0 of 117 "
                  "programs, not a BICC1 pass. The BICC1-only call belongs to the 62-donor sensitivity.",
        "source": f"{H59} (selected_H59, perm_q, ols60_n, ols62_n, ols62_q, ols62_q_hc3; row hepatocytes_20)",
        "numbers_in_new": "60; 0; 117; 0.323; 62; 0.038",
    },
    {
        "id": "I-02",
        "file": "figures/main/fig4_singlecell_programs/CANDIDATE_PANEL_INDEX.tsv",
        "old": "panels/fig4c_hotspot_stage_heatmap.pdf\tscripts/analysis/bulk_refit_followup/render_verify_qualified_programs.py\t"
               "source_diagnosis_corrected_20261014;source_tables/source_diagnosis_correction_20261014/\t"
               "62_max_donors_fixed_three_stage_117_BH_BICC1_only_primary_selected",
        "new": f"panels/fig4c_hotspot_stage_heatmap.pdf\t{REP4}/producers/02r_render_fig4c_h59.R\t"
               f"H59_repair_20261015;{REP4}/source_tables/\t"
               "H59_60_donors_per_dataset_variance_117_permutation_BH_0_selected",
        "reason": "Row 4C gate says BICC1 is primary-selected; the installed panel (sha256 = H59 repair render) "
                  "shows 0 of 117 selected. Generator and state updated to the producer the hash identifies.",
        "source": f"{H59} (selected_H59, ols60_n); sha256 of installed panel vs {REP4}/fig4c_hotspot_stage_heatmap_H59_candidate_repair.pdf",
        "numbers_in_new": "60; 117; 0",
    },
    {
        "id": "I-03",
        "file": "figures/main/fig4_singlecell_programs/CANDIDATE_PANEL_INDEX.tsv",
        "old": "panels/supplementary/figs4c_hotspot_stage_heatmap_full.pdf\tscripts/analysis/bulk_refit_followup/render_verify_qualified_programs.py\t"
               "source_diagnosis_corrected_20261014\t62_max_donors_fixed_three_stage_117_BH_BICC1_only_primary_selected",
        "new": f"panels/supplementary/figs4c_hotspot_stage_heatmap_full.pdf\t{REP4}/producers/02r_render_fig4c_h59.R\t"
               "H59_repair_20261015\tH59_60_donors_per_dataset_variance_117_permutation_BH_0_selected",
        "reason": "Row S4C gate says BICC1 is primary-selected; the installed panel (sha256 = H59 repair render) "
                  "shows 0 of 117 selected.",
        "source": f"{H59} (selected_H59, ols60_n); sha256 of installed panel vs {REP4}/figs4c_hotspot_stage_heatmap_full_H59_candidate_repair.pdf",
        "numbers_in_new": "60; 117; 0",
    },
    {
        "id": "I-04",
        "file": "docs/manuscript/05_figure_legends.md",
        "old": "**(B)** Genetics-call by RNA-call matrix over 28,233 selected r8 gene identities.",
        "new": "**(B)** Genetics-call by RNA-call matrix over 28,233 gene identities with r8\n"
               "classes, unchanged in Catalog r9d (adopted 2026-10-16).",
        "reason": "Release label: the adopted Catalog is r9d, which keeps the 28,233 r8 identities and classes.",
        "source": f"{R9D}/catalog_r9d_gene_table.parquet (row count); {R9D}/invariants.tsv (r8 class counts); docs/RESULTS.md 'Catalog r9d was adopted on 2026-10-16'",
        "numbers_in_new": "28,233; 2026-10-16",
    },
    {
        "id": "I-05",
        "file": "docs/manuscript/05_figure_legends.md",
        "old": "The local Catalog uses the\nsame r8 gene calls as the panels.",
        "new": "The local Catalog uses the\nsame gene calls as the panels: r8 classes, unchanged in Catalog r9d\n(adopted 2026-10-16).",
        "reason": "Release label: the adopted Catalog is r9d; its classes are the r8 classes.",
        "source": f"{R9D}/invariants.tsv (r8 class counts passed); docs/RESULTS.md 'Catalog r9d was adopted on 2026-10-16'",
        "numbers_in_new": "2026-10-16",
    },
    {
        "id": "I-06",
        "file": "docs/manuscript/05_figure_legends.md",
        "old": "All five individual panels are installed on October 14;",
        "new": "Panels A and E were installed on October 14 and panels B, C and D are the\n"
               "October 15 repair renders;",
        "reason": "Installed fig6b/c/d.pdf are byte-identical to the October 15 repair renders; fig6a/e match the "
                  "October 14 r8 panel build.",
        "source": "work/fig6_sha256.tsv: figures/main/fig6_gene_catalog/panels/fig6{a..e}*.pdf vs "
                  "figures/candidates/review-r3-fixes-20261015T112148Z/F6/repair/panels_20261015T173321Z/ (B, C, D) and "
                  "figures/candidates/catalog-r8-panels-20261014T143412Z_22236670/panels/ (A, E)",
        "numbers_in_new": "October 14; October 15",
    },
]

rr = pd.read_csv(os.path.join(WORK, "index_reread.tsv"), sep="\t", dtype=str)
for e in edits:
    with open(os.path.join(ROOT, e["file"]), encoding="utf-8") as fh:
        live = fh.read()
    e["old_count"] = live.count(e["old"])
    sub = rr[rr["edit_id"] == e["id"]]
    e["value_reread"] = "; ".join(f"{c}={v}" for c, v in zip(sub["column_or_row"], sub["value_reread"]))
    assert e["old_count"] == 1, (e["id"], e["old_count"])
    assert "—" not in e["new"], e["id"]

with open(os.path.join(WORK, "index_edits.json"), "w", encoding="utf-8") as fh:
    json.dump(edits, fh, indent=1, ensure_ascii=False)
for e in edits:
    print(e["id"], e["file"], "old_count", e["old_count"])
