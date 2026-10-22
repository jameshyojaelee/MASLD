"""Re-read every number used in the I- edits (INDEX.md, fig4 CANDIDATE_PANEL_INDEX.tsv,
Fig 6 legend) and the Fig 5D/5E source-table facts. Read-only on all inputs.
Writes TSV to stdout: edit_id, source_path, column_or_row, value_reread.
Also writes work/fig6_sha256.tsv and work/fig5_sourcetables.tsv."""
import hashlib
import os
import sys

import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
WORK = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-QUALITY-20261017T020832Z/text_fixes-20261017T025625Z/work")
rows = []


def add(eid, path, col, val):
    rows.append((eid, os.path.relpath(path, ROOT) if path.startswith(ROOT) else path, col, str(val)))


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---- H59 primary family (I-01, I-02, I-03)
h59 = os.path.join(ROOT, "Analysis/Hypotheses/H59-20261015T112432Z/out/outcome/results_primary_family.tsv")
d = pd.read_csv(h59, sep="\t")
duct = d.loc[d["program"] == "hepatocytes_20"].iloc[0]
for eid in ("I-01", "I-02", "I-03"):
    add(eid, h59, "n rows (programs)", len(d))
    add(eid, h59, "sum(selected_H59)", int(d["selected_H59"].sum()))
    add(eid, h59, "min(perm_q)", f"{d['perm_q'].min():.6f}")
    add(eid, h59, "max(ols60_n) = donor universe of the primary fit", int(d["ols60_n"].max()))
add("I-01", h59, "max(ols62_n) = equal-variance sensitivity donors", int(d["ols62_n"].max()))
add("I-01", h59, "count(ols62_q<0.05)", int((d["ols62_q"] < 0.05).sum()))
add("I-01", h59, "hepatocytes_20 module_name; ols62_q", f"{duct['module_name']}; {duct['ols62_q']:.6f}")
add("I-01", h59, "count(ols62_q_hc3<0.05)", int((d["ols62_q_hc3"] < 0.05).sum()))

# installed 4C / S4C are the H59 repair renders
pairs = {
    "I-02": ("figures/main/fig4_singlecell_programs/panels/fig4c_hotspot_stage_heatmap.pdf",
             "figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/repair/fig4c_hotspot_stage_heatmap_H59_candidate_repair.pdf"),
    "I-03": ("figures/main/fig4_singlecell_programs/panels/supplementary/figs4c_hotspot_stage_heatmap_full.pdf",
             "figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/repair/figs4c_hotspot_stage_heatmap_full_H59_candidate_repair.pdf"),
}
for eid, (inst, cand) in pairs.items():
    a, b = sha(os.path.join(ROOT, inst)), sha(os.path.join(ROOT, cand))
    add(eid, os.path.join(ROOT, inst), "sha256 installed == repair render", f"{a == b} ({a[:12]})")
log = os.path.join(ROOT, "figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/repair/producers/02r_render.log")
with open(log) as fh:
    first = fh.readline().strip()
add("I-02", log, "render check line", first)

# ---- Catalog r9d (I-04, I-05)
r9d = os.path.join(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/"
                   "catalog_r9d-20261016T105550Z/build_out")
gt = pd.read_parquet(os.path.join(r9d, "catalog_r9d_gene_table.parquet"), columns=["ensembl_id"])
inv = pd.read_csv(os.path.join(r9d, "invariants.tsv"), sep="\t", dtype=str)
cls = inv.loc[inv["check"] == "r8 class counts"].iloc[0]
with open(os.path.join(ROOT, "docs/RESULTS.md")) as fh:
    adopted = [ln.strip() for ln in fh if ln.startswith("Catalog r9d was adopted on")]
for eid in ("I-04", "I-05"):
    add(eid, os.path.join(r9d, "catalog_r9d_gene_table.parquet"), "n rows", len(gt))
    add(eid, os.path.join(r9d, "invariants.tsv"), "r8 class counts passed", f"{cls['passed']} {cls['observed']}")
    add(eid, os.path.join(ROOT, "docs/RESULTS.md"), "line 'Catalog r9d was adopted on'", adopted[0][:60])

# ---- Fig 6 panel provenance by sha256 (I-06)
inst_dir = os.path.join(ROOT, "figures/main/fig6_gene_catalog/panels")
cands = {
    "fig6a_catalog_structure.pdf": "figures/candidates/catalog-r8-panels-20261014T143412Z_22236670/panels/fig6a_catalog_structure.pdf",
    "fig6b_coverage_observability.pdf": "figures/candidates/review-r3-fixes-20261015T112148Z/F6/repair/panels_20261015T173321Z/fig6b_coverage_observability_candidate.pdf",
    "fig6c_worked_entries.pdf": "figures/candidates/review-r3-fixes-20261015T112148Z/F6/repair/panels_20261015T173321Z/fig6c_worked_entries_candidate.pdf",
    "fig6d_next_experiment.pdf": "figures/candidates/review-r3-fixes-20261015T112148Z/F6/repair/panels_20261015T173321Z/fig6d_next_experiment_candidate.pdf",
    "fig6e_portal_source_access.pdf": "figures/candidates/catalog-r8-panels-20261014T143412Z_22236670/panels/fig6e_portal_source_access.pdf",
}
sha_rows = []
for name, cand in cands.items():
    ip = os.path.join(inst_dir, name)
    cp = os.path.join(ROOT, cand)
    a, b = sha(ip), sha(cp)
    sha_rows.append({"installed_panel": os.path.relpath(ip, ROOT), "installed_sha256": a,
                     "matched_candidate": cand, "candidate_sha256": b, "match": a == b,
                     "candidate_dir_timestamp": cand.split("/")[2] if "catalog-r8" in cand else cand.split("/")[5]})
    add("I-06", ip, f"sha256 == {sha_rows[-1]['candidate_dir_timestamp']}", f"{a == b} ({a[:12]})")
pd.DataFrame(sha_rows).to_csv(os.path.join(WORK, "fig6_sha256.tsv"), sep="\t", index=False)

# ---- Fig 5D / 5E source tables (needs_user, no edit)
tables = [
    ("installed", "figures/main/fig5_molecular_context/panels/data/fig5d_variant_consistent_accessibility.tsv"),
    ("installed", "figures/main/fig5_molecular_context/panels/data/fig5e_multimodal_program_summary_selection.tsv"),
    ("candidate", "figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/repair/fig5d_gene_level_display_counts.tsv"),
    ("candidate", "figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/repair/fig5d_locus_background_plotted_points.tsv"),
    ("candidate", "figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/fig5d_locus_background_summary.tsv"),
    ("candidate", "figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/fig5d_pair_rows_with_background.tsv"),
    ("candidate", "figures/candidates/review-r3-fixes-20261015T112148Z/F5ef/fig5e_program_context_117_candidate.tsv"),
]
trows = []
for kind, rel in tables:
    p = os.path.join(ROOT, rel)
    t = pd.read_csv(p, sep="\t")
    extra = ""
    if rel.endswith("fig5e_program_context_117_candidate.tsv"):
        extra = (f"display TRUE {int(t['display'].sum())}; n_open_genes>0 {int((t['n_open_genes'] > 0).sum())}; "
                 f"robust {int(t['robust'].fillna(False).astype(bool).sum())}")
    if rel.endswith("fig5e_multimodal_program_summary_selection.tsv"):
        extra = f"selection_rule {t['selection_rule'].iloc[0]}"
    if rel.endswith("fig5d_variant_consistent_accessibility.tsv"):
        extra = f"denominators {sorted(t['denominator'].unique().tolist())}"
    if rel.endswith("fig5d_gene_level_display_counts.tsv"):
        extra = f"denominators {sorted(t['denominator'].unique().tolist())}; partial_floor {t['partial_floor'].unique().tolist()}"
    if rel.endswith("fig5d_locus_background_summary.tsv"):
        extra = f"bootstrap {t['bootstrap'].unique().tolist()}"
    trows.append({"kind": kind, "path": rel, "sha256": sha(p), "n_rows": len(t), "note": extra})
    add("NU-5DE", p, "sha256; n_rows; note", f"{trows[-1]['sha256'][:12]}; {len(t)}; {extra}")
pd.DataFrame(trows).to_csv(os.path.join(WORK, "fig5_sourcetables.tsv"), sep="\t", index=False)
for inst, cand in [("figures/main/fig5_molecular_context/panels/fig5d_snatac_accessibility.pdf",
                    "figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/repair/fig5d_candidate_gene_level_mass_floor.pdf"),
                   ("figures/main/fig5_molecular_context/panels/fig5d_companion_locus_matched_background.pdf",
                    "figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/repair/fig5d_candidate_locus_matched_background.pdf"),
                   ("figures/main/fig5_molecular_context/panels/fig5e_multimodal_program_summary.pdf",
                    "figures/candidates/review-r3-fixes-20261015T112148Z/F5ef/fig5e_program_context_117_candidate.pdf")]:
    a, b = sha(os.path.join(ROOT, inst)), sha(os.path.join(ROOT, cand))
    add("NU-5DE", os.path.join(ROOT, inst), f"sha256 == {cand.split('/')[-1]}", f"{a == b} ({a[:12]})")

w = sys.stdout
w.write("edit_id\tsource_path\tcolumn_or_row\tvalue_reread\n")
for r in rows:
    w.write("\t".join(r) + "\n")
