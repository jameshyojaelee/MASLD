#!/usr/bin/env python3
# KEY MESSAGE: A checksum-linked candidate package exposes every panel, source table, proof, and unresolved promotion gate without touching canonical figures.

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
CANDIDATE = Path(os.environ["FIGURE_CANDIDATE_ROOT"])
GENETICS = Path(os.environ["GENETICS_RELEASE_ROOT"]) if os.environ.get("GENETICS_RELEASE_ROOT") else None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_new(source: Path, target: Path) -> None:
    if target.exists():
        raise RuntimeError(f"Refusing overwrite: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    if path.exists():
        raise RuntimeError(f"Refusing overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# Copy the current Figure 5 panels into the isolated PI-review candidate.
figure5 = CANDIDATE / "figure5" / "panels"
figure6 = CANDIDATE / "figure6" / "panels"
for name in [
    "fig5c_mrna_protein_composite.pdf",
    "fig5d_snatac_accessibility.pdf",
    "fig5e_multimodal_program_summary.pdf",
    "fig5f_spatial_program_calibration.pdf",
]:
    copy_new(ROOT / "figures/main/fig5_molecular_context/panels" / name, figure5 / name)
# Canonical panels were renamed fig5*->fig6* on 2026-08-23, so the historical
# old->new translation is now an identity copy; names are already fig6-correct.
for name in [
    "fig6b_coverage_observability.pdf",
    "fig6d_next_experiment.pdf",
]:
    copy_new(ROOT / "figures/main/fig6_gene_catalog/panels" / name, figure6 / name)

source_dir = CANDIDATE / "source_tables"
for source in sorted((ROOT / "figures/main/fig5_molecular_context/data").glob("composite_*.csv")):
    copy_new(source, source_dir / f"fig5_{source.name}")
for source in [
    ROOT / "figures/main/fig5_molecular_context/data/fig5d_snatac_accessibility_source.tsv",
    ROOT / "figures/main/fig5_molecular_context/panels/data/fig5e_multimodal_program_summary_selection.tsv",
    ROOT / "figures/main/fig5_molecular_context/panels/data/fig5f_spatial_program_calibration.tsv",
]:
    copy_new(source, source_dir / source.name)
for source in sorted((CANDIDATE / "figure5/panels/data").glob("fig5*.tsv")):
    copy_new(source, source_dir / source.name)
for name in [
    "fig5_passport_call_matrix.csv", "fig5_passport_domain_margin.csv",
    "fig5_passport_hero_states.csv", "fig5_passport_next_experiment.csv",
    "fig5_passport_provenance.csv",
]:
    copy_new(ROOT / "figures/main/fig6_gene_catalog/data" / name,
             source_dir / name.replace("fig5_", "fig6_", 1))

# The apparent 113/27 versus 114/26 mismatch is an estimand difference, not a
# drifting registry. Weighted projection requires at least 0.80 retained weight
# and tests 113 programs; the unweighted vocabulary check uses at least 0.80
# gene coverage and additionally admits T-cell module 9 (0.95 gene coverage but
# 0.783 retained weight). Its 26 calls use unweighted effects and a 228-test
# fibrosis/NAS family, whereas 27 uses weighted effects in the 226-test family.
reconciliation = [
    {
        "view": "weighted_axis_map",
        "n_registry": 117,
        "n_testable": 113,
        "n_fibrosis_supported": 27,
        "family_size": 226,
        "eligibility": "retained source-weight coverage >=0.80",
        "extra_or_excluded_program": "four programs below retained-weight threshold",
    },
    {
        "view": "unweighted_vocabulary_check",
        "n_registry": 117,
        "n_testable": 114,
        "n_fibrosis_supported": 26,
        "family_size": 228,
        "eligibility": "mapped gene coverage >=0.80",
        "extra_or_excluded_program": "hotspot_tcells_eaa9c262b58d9323 admitted at 0.95 gene coverage; retained-weight coverage=0.783",
    },
]
write_tsv(source_dir / "figs4_program_count_reconciliation.tsv", reconciliation,
          ["view", "n_registry", "n_testable", "n_fibrosis_supported", "family_size", "eligibility", "extra_or_excluded_program"])
axis_map_path = ROOT / (
    "RNA-seq/results/manuscript_release/candidates/"
    "resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/"
    "BULK-PROGRAM-MAP-v9/discovery/program_axis_map.tsv"
)
coverage_path = ROOT / (
    "RNA-seq/results/manuscript_release/candidates/"
    "resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/"
    "BULK-PROGRAM-MAP-v9/systems/T3_vocabulary/vocabulary_coverage.tsv"
)
axis_rows = read_tsv(axis_map_path)
coverage_by_id = {row["feature_id"]: row for row in read_tsv(coverage_path) if row["vocabulary"] == "hotspot_117"}
exceptions = []
for row in axis_rows:
    if row["testable"] != "FALSE":
        continue
    coverage = coverage_by_id[row["feature_id"]]
    exceptions.append({
        "program_uid": row["feature_id"],
        "cell_type": row["cell_type"],
        "module": row["module"],
        "module_name": row["module_name"],
        "retained_weight_coverage": row["weight_coverage"],
        "weighted_axis_state": "untestable_below_0.80_retained_weight",
        "unweighted_gene_coverage": coverage["gene_coverage"],
        "unweighted_vocabulary_state": "testable" if coverage["testable"] == "TRUE" else "untestable",
    })
if len(exceptions) != 4:
    raise RuntimeError("Weighted program-testability exception count drift")
write_tsv(source_dir / "figs4_program_testability_exceptions.tsv", exceptions,
          ["program_uid", "cell_type", "module", "module_name", "retained_weight_coverage",
           "weighted_axis_state", "unweighted_gene_coverage", "unweighted_vocabulary_state"])
systems_source = ROOT / "figures/candidates/program-systems-2026-08-12-v3"
systems_target = CANDIDATE / "supplementary" / "figureS4" / "panels"
singlecell_supplement_source = ROOT / "figures/supplementary/figS03_deconvolution/figS_singlecell.pdf"
copy_new(singlecell_supplement_source, systems_target / "figs4_singlecell_detail.pdf")
for old, new in [
    ("fig_system_t3_vocabulary.pdf", "figs4a_vocabulary_robustness.pdf"),
    ("fig_system_t4_coordination.pdf", "figs4b_residual_program_coordination.pdf"),
    ("fig_system_t5_composition.pdf", "figs4c_composition_coupling.pdf"),
    ("fig_system_t6_modifiers.pdf", "figs4d_sex_effect_modification.pdf"),
]:
    copy_new(systems_source / old, systems_target / new)
for source in sorted(systems_source.glob("panel_*_source.tsv")):
    copy_new(source, source_dir / f"figs4_{source.name}")

genetics_promoted = False
if GENETICS:
    marker = GENETICS / "PROMOTED.json"
    if marker.exists():
        payload = json.loads(marker.read_text())
        genetics_promoted = (
            payload.get("status") == "promoted"
            and payload.get("terminal") is True
            and payload.get("canonical_promotion_authorized") is True
        )

blockers = [
    {
        "gate": "corrected_coloc",
        "status": "pass" if genetics_promoted else "blocked",
        "affected": "3F; 4E genetic anchors",
        "reason": "requires promoted terminal corrected-COLOC manifest, checksums, and exact signal-pair export",
    },
    {
        "gate": "program_system_counts",
        "status": "pass",
        "affected": "supplementary vocabulary, coordination, composition, and sex-modification panels",
        "reason": "reconciled as weighted retained-mass eligibility versus unweighted gene-coverage eligibility; exact extra program documented",
    },
    {
        "gate": "promotion",
        "status": "blocked",
        "affected": "all candidate panels",
        "reason": "explicit release-promotion approval has not been given",
    },
]
write_tsv(CANDIDATE / "BLOCKERS.tsv", blockers, ["gate", "status", "affected", "reason"])

supplement_renumbering = [
    {
        "target": "Figure S3",
        "source": "candidate supplementary/figureS3",
        "status": "ready_non_genetic",
        "content": "full stage-remodeling panel with donor-aware composition; genetics alternatives remain gated",
    },
    {
        "target": "Figure S4",
        "source": "former single-cell supplement plus program-systems-2026-08-12-v3",
        "status": "candidate",
        "content": "single-cell detail; vocabulary robustness; residual coordination; composition coupling; sex effect modification",
    },
    {
        "target": "Figure S5",
        "source": "figures/supplementary/figS04_coloc",
        "status": "renumber_on_promotion",
        "content": "former Figure S4 COLOC supplement",
    },
    {
        "target": "Figure S6",
        "source": "figures/supplementary/figS05_epigenomic_spatial",
        "status": "renumber_on_promotion",
        "content": "former Figure S5 epigenomic and spatial supplement",
    },
]
write_tsv(CANDIDATE / "supplement_renumbering.tsv", supplement_renumbering,
          ["target", "source", "status", "content"])

callouts: list[dict[str, object]] = []
for figure, names in {
    "3": [
        "fig3a_cohort_metadata_matrix.pdf", "fig3b_nas_fib_grid.pdf",
        "fig3c_cohort_alluvial.pdf", "fig3d_pca_fibrosis_gradient.pdf",
        "fig3e_stage_remodeling.pdf", "fig3f_stage_deg_genetics_matrix.pdf",
    ],
    "4": [
        "fig4a_singlecell_atlas_schematic.pdf", "fig4b_scrna_umap_embeddable.pdf",
        "fig4c_program_landscape.pdf", "fig4d_communication.pdf",
        "fig4e_bulk_projection.pdf", "fig4f_tf_activity.pdf",
    ],
    "5": [
        "fig5a_input_firewall.pdf", "fig5b_protein_triage.pdf",
        "fig5c_mrna_protein_composite.pdf", "fig5d_snatac_accessibility.pdf",
        "fig5e_multimodal_program_summary.pdf", "fig5f_spatial_program_calibration.pdf",
    ],
    "6": [
        "fig6a_catalog_structure.pdf", "fig6b_coverage_observability.pdf",
        "fig6c_worked_entries.pdf", "fig6d_next_experiment.pdf",
        "fig6e_portal_source_graph.pdf",
    ],
}.items():
    for index, name in enumerate(names):
        path = CANDIDATE / f"figure{figure}" / "panels" / name
        if figure == "4" and name == "fig4e_bulk_projection.pdf" and not path.exists():
            partial = path.with_name("fig4e_bulk_projection_pre_genetics.pdf")
            if partial.exists():
                path = partial
                status = "partial"
            else:
                status = "blocked"
        else:
            status = "ready" if path.exists() else "blocked"
        callouts.append({
            "callout": f"{figure}{chr(65 + index)}",
            "relative_path": str(path.relative_to(CANDIDATE)),
            "status": status,
            "blocker": "" if status == "ready" else (
                "corrected_coloc" if figure in {"3", "4"} else "missing_candidate_contract"
            ),
        })
write_tsv(CANDIDATE / "callout_map.tsv", callouts, ["callout", "relative_path", "status", "blocker"])

(CANDIDATE / "captions").mkdir(parents=True, exist_ok=True)
(CANDIDATE / "captions" / "figure3.md").write_text(
    "Figure 3. Bulk transcriptomic remodeling and genetic overlap. A-C, synchronized five-cohort overview, stage/NAS coverage, and per-cohort DEG replication. D, covariate-adjusted PCA; control is gray and disease uses the fibrosis blue scale. E, adjacent-stage DEG counts and continuous prespecified NMF activity. F, adjacent-stage effects with separate exact-pair COLOC PP.H4 and lead-variant PIP lanes; pending adoption of the corrected COLOC release.\n",
    encoding="utf-8",
)
(CANDIDATE / "captions" / "figure4.md").write_text(
    "Figure 4. Single-cell programs and regulatory activity. A, seven datasets comprising 102 analyzed biological donors and 1.23 million analyzed cells define an annotated atlas and 117 programs in five tested lineages. B, directly labeled atlas UMAP with disease-density scale. C, all 117 donor-level stage effects and four selected three-stage profiles; outlines mark the two HC3-supported programs. D, fixed 13-pair donor-collapsed communication roster; zero of 13 pass donor-level family FDR, and failed or untestable rows remain visible as muted cells. E, two HC3-supported signatures projected into fragment-native F-versus-F0 models; genetic anchors await adoption of the corrected COLOC release. F, DoRothEA A/B/C TF activity from decoupleR weighted-mean inference; filled points meet within-contrast BH q<0.05.\n",
    encoding="utf-8",
)
(CANDIDATE / "captions" / "figure5.md").write_text(
    "Figure 5. Assay-native molecular and physical context. Panel B retains every measured prioritized gene in the scatter; its adjacent bar is restricted to protein-significant genes and reports concordant and discordant percentages within that denominator.\n",
    encoding="utf-8",
)
(CANDIDATE / "captions" / "figure6.md").write_text(
    "Figure 6. The MASLD Gene Catalog. Candidate components currently cover B, evidence coverage and observability states, and D, deterministic next-experiment routing. Panels A, C, and E remain explicit candidate gaps.\n",
    encoding="utf-8",
)
(CANDIDATE / "RETIREMENT_NOTES.md").write_text(
    "The prior five-figure callout maps remain in place. This candidate supersedes their numbering only after it is explicitly adopted. Single-cell supplementary callouts move to S4; prior S4 and S5 callouts move to S5 and S6. No historical map was deleted.\n",
    encoding="utf-8",
)

# Fail closed on scientific cardinality and checksum contracts before PDF
# validation. These checks use the exact candidate tables, not plot objects.
stage_dir = CANDIDATE / "analysis/stage_extensions"
stage_manifest = read_tsv(stage_dir / "output_manifest.tsv")
for row in stage_manifest:
    artifact = stage_dir / row["relative_path"]
    if not artifact.exists() or artifact.stat().st_size != int(row["size_bytes"]) or sha256(artifact) != row["sha256"]:
        raise RuntimeError(f"Stage output manifest mismatch: {artifact}")

stage_rows = read_tsv(stage_dir / "stage_extension_all_gene_results.tsv")
stage_groups: dict[tuple[str, str], int] = {}
for row in stage_rows:
    key = (row["axis"], row["contrast"])
    stage_groups[key] = stage_groups.get(key, 0) + 1
if len(stage_groups) != 14 or set(stage_groups.values()) != {23370}:
    raise RuntimeError(f"Stage family cardinality drift: {stage_groups}")

program_rows = read_tsv(source_dir / "fig4c_all_117_programs.tsv")
if len(program_rows) != 117 or len({row["program_uid"] for row in program_rows}) != 117:
    raise RuntimeError("Figure 4C does not preserve all 117 programs")
if sum(row["primary_selected"] == "TRUE" for row in program_rows) != 4 or sum(
    row["primary_selected"] == "TRUE" and row["hc3_supported"] == "TRUE" for row in program_rows
) != 2:
    raise RuntimeError("Figure 4C selected or HC3 state drift")

communication_rows = read_tsv(source_dir / "fig4d_fixed_13_pair_communication.tsv")
if len(communication_rows) != 39 or len({row["headline_label"] for row in communication_rows}) != 13:
    raise RuntimeError("Figure 4D fixed communication roster drift")

tf_rows = read_tsv(source_dir / "fig4f_all_tf_activity.tsv")
tf_contrasts = {row["contrast"] for row in tf_rows}
if len(tf_contrasts) != 4:
    raise RuntimeError("Figure 4F contrast cardinality drift")
for contrast in tf_contrasts:
    rows = [row for row in tf_rows if row["contrast"] == contrast]
    if any(int(row["bh_family_size"]) != len(rows) for row in rows):
        raise RuntimeError(f"Figure 4F BH family drift: {contrast}")

protein_rows = read_tsv(source_dir / "fig5b_protein_triage_summary.tsv")
if len(protein_rows) != 1:
    raise RuntimeError("Figure 5B protein summary cardinality drift")
protein = protein_rows[0]
if int(float(protein["n_protein_significant"])) != int(float(protein["n_significant_concordant"])) + int(float(protein["n_significant_discordant"])):
    raise RuntimeError("Figure 5B protein-significant denominator drift")
if int(float(protein["n_prioritized_measured"])) != int(float(protein["n_protein_significant"])) + int(float(protein["n_not_significant"])):
    raise RuntimeError("Figure 5B full-measured denominator drift")


def pdf_page_count(path: Path) -> int:
    result = subprocess.run(["pdfinfo", str(path)], check=True, text=True, capture_output=True)
    for line in result.stdout.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise RuntimeError(f"No page count: {path}")


validation_rows = []
panel_pdfs = sorted(CANDIDATE.glob("**/panels/*.pdf"))
for path in panel_pdfs:
    gs = subprocess.run(
        ["gs", "-q", "-dNOPAUSE", "-dBATCH", "-dPDFDEBUG", "-sDEVICE=nullpage", str(path)],
        capture_output=True,
    )
    raw = path.read_bytes()
    debug = (gs.stdout + gs.stderr).decode("latin-1")
    has_type3 = (
        b"/Subtype /Type3" in raw
        or b"/FontType 3" in raw
        or "/Subtype /Type3" in debug
        or "/FontType 3" in debug
    )
    validation_rows.append({
        "relative_path": str(path.relative_to(CANDIDATE)),
        "pages": pdf_page_count(path),
        "ghostscript": "pass" if gs.returncode == 0 else "fail",
        "type3_check": "fail" if has_type3 else "pass",
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    })
write_tsv(CANDIDATE / "pdf_validation.tsv", validation_rows,
          ["relative_path", "pages", "ghostscript", "type3_check", "size_bytes", "sha256"])
if any(row["pages"] != 1 or row["ghostscript"] != "pass" or row["type3_check"] != "pass"
       for row in validation_rows):
    raise RuntimeError("One or more individual panel PDFs failed structural validation")


def render_pdf(path: Path, temp: Path) -> Image.Image:
    prefix = temp / path.stem
    subprocess.run(["pdftoppm", "-f", "1", "-singlefile", "-r", "120", "-png", str(path), str(prefix)], check=True)
    return Image.open(prefix.with_suffix(".png")).convert("RGB")


def composite(name: str, entries: list[tuple[str, Path]], layout: list[tuple[int, int, int, int]]) -> None:
    canvas = Image.new("RGB", (1420, 1800), "white")
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 28)
    with tempfile.TemporaryDirectory(prefix="pi-proof-") as tmp:
        temp = Path(tmp)
        for (label, path), (x, y, w, h) in zip(entries, layout):
            if not path.exists():
                draw.rectangle((x, y, x + w, y + h), outline="#9E9E9E", width=2)
                draw.text((x + 12, y + 12), f"{label} blocked", fill="#666666", font=font)
                continue
            image = render_pdf(path, temp)
            image.thumbnail((w - 24, h - 24), Image.Resampling.LANCZOS)
            px = x + (w - image.width) // 2
            py = y + (h - image.height) // 2
            canvas.paste(image, (px, py))
            draw.text((x + 4, y + 4), label, fill="black", font=font)
    proof = CANDIDATE / "proofs" / f"{name}.pdf"
    proof.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(proof, "PDF", resolution=200.0)


def panel(figure: int, filename: str) -> Path:
    return CANDIDATE / f"figure{figure}" / "panels" / filename


composite("figure3_composite_proof", [(x[0], panel(3, x[1])) for x in [
    ("A", "fig3a_cohort_metadata_matrix.pdf"), ("B", "fig3b_nas_fib_grid.pdf"),
    ("C", "fig3c_cohort_alluvial.pdf"), ("D", "fig3d_pca_fibrosis_gradient.pdf"),
    ("E", "fig3e_stage_remodeling.pdf"), ("F", "fig3f_stage_deg_genetics_matrix.pdf")]],
    [(0, 0, 474, 720), (474, 0, 473, 720), (947, 0, 473, 720),
     (0, 720, 420, 1080), (420, 720, 650, 1080), (1070, 720, 350, 1080)])

e4 = "fig4e_bulk_projection.pdf" if panel(4, "fig4e_bulk_projection.pdf").exists() else "fig4e_bulk_projection_pre_genetics.pdf"
composite("figure4_composite_proof", [(x[0], panel(4, x[1])) for x in [
    ("A", "fig4a_singlecell_atlas_schematic.pdf"), ("B", "fig4b_scrna_umap_embeddable.pdf"),
    ("C", "fig4c_program_landscape.pdf"), ("D", "fig4d_communication.pdf"),
    ("E", e4), ("F", "fig4f_tf_activity.pdf")]],
    [(0, 0, 596, 500), (596, 0, 824, 500), (0, 500, 1420, 650),
     (0, 1150, 474, 650), (474, 1150, 473, 650), (947, 1150, 473, 650)])

grid6 = [(0, 0, 474, 900), (474, 0, 473, 900), (947, 0, 473, 900),
         (0, 900, 474, 900), (474, 900, 473, 900), (947, 900, 473, 900)]
composite("figure5_composite_proof", [(chr(65+i), p) for i, p in enumerate(sorted(figure5.glob("fig5*.pdf")))], grid6)
composite("figure6_composite_proof", [("B", panel(6, "fig6b_coverage_observability.pdf")),
                                      ("D", panel(6, "fig6d_next_experiment.pdf"))],
          [(0, 200, 710, 1300), (710, 200, 710, 1300)])

# Candidate-level manifests are written last so they cover every deliverable.
inputs = [
    ROOT / "scripts/figures/build_pi_stage_extensions.R",
    ROOT / "scripts/figures/nas_fib_grid.R",
    ROOT / "scripts/figures/render_pi_bulk_opening_panels.R",
    ROOT / "scripts/figures/figS_pca_definitive.R",
    ROOT / "scripts/figures/progression_cascade.R",
    ROOT / "scripts/figures/gen_scrna_umap_embeddable.R",
    ROOT / "scripts/figures/singlecell_module_heatmap.R",
    ROOT / "scripts/figures/render_pi_singlecell_panels.R",
    ROOT / "scripts/figures/stage_deg_genetics_matrix.R",
    ROOT / "scripts/figures/tf_convergence_vsF0.R",
    ROOT / "scripts/figures/fig5a_input_firewall.R",
    ROOT / "scripts/figures/fig5b_protein_triage.R",
    ROOT / "scripts/figures/finalize_pi_figure_candidate.py",
    ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv",
    ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/donor_program_scores_primary.tsv",
    ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv",
    ROOT / "figures/main/fig4_singlecell_programs/source_tables/legacy_from_fig3/ccc_trajectories_data_dc.csv",
    ROOT / "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2/ccc_dc_before_after_report.tsv",
    axis_map_path,
    coverage_path,
    systems_source / "output_manifest.tsv",
    singlecell_supplement_source,
]
inputs.extend(sorted((ROOT / "figures/main/fig5_molecular_context/panels").glob("fig5[cd-f]_*.pdf")))
inputs.extend(sorted((ROOT / "figures/main/fig6_gene_catalog/panels").glob("fig6[bd]_*.pdf")))
input_rows = [{"path": str(p), "size_bytes": p.stat().st_size, "sha256": sha256(p)} for p in inputs]
write_tsv(CANDIDATE / "input_manifest.tsv", input_rows, ["path", "size_bytes", "sha256"])

outputs = sorted(p for p in CANDIDATE.rglob("*") if p.is_file() and p.name != "output_manifest.tsv")
output_rows = [{"relative_path": str(p.relative_to(CANDIDATE)), "size_bytes": p.stat().st_size, "sha256": sha256(p)} for p in outputs]
write_tsv(CANDIDATE / "output_manifest.tsv", output_rows, ["relative_path", "size_bytes", "sha256"])
print(f"CANDIDATE_FINALIZED: {CANDIDATE}")
