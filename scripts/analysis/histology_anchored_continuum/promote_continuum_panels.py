#!/usr/bin/env python3
"""Promote the newest validated continuum panels into Figure 4/S3/S4 authorities."""

from __future__ import annotations

import csv
import hashlib
import os
import shutil
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", ""))
if not ROOT.is_dir():
    raise SystemExit("MASLD_PROJECT_ROOT is unset or invalid")

REVIEW = ROOT / "figures/candidates/histology-continuum-figure-review-20260818T151937Z/panels"
MOLECULAR = ROOT / (
    "RNA-seq/results/histology_anchored_continuum/molecular_layers/"
    "hac-molecular-layers-20260818T173348Z/figures/panels"
)
COMPARISON = ROOT / "figures/candidates/histology-vs-continuum-comparison-20260818T192630Z/panels"
WINDOWS = ROOT / "figures/candidates/histology-continuum-hotspot-windows-20260824T160938Z/panels"
FIG4 = ROOT / "figures/main/fig4_singlecell_programs"
FIGS3 = ROOT / "figures/supplementary/figS03_bulk_transcriptomics"


@dataclass(frozen=True)
class Panel:
    authority: str
    callout: str
    source: Path
    destination: Path
    expected_sha256: str
    role: str
    generator: str
    selection_reason: str


def panel(
    authority: str,
    callout: str,
    source: Path,
    destination: Path,
    expected_sha256: str,
    role: str,
    generator: str,
    selection_reason: str,
) -> Panel:
    return Panel(
        authority, callout, source, destination, expected_sha256, role,
        generator, selection_reason,
    )


PANELS = [
    panel(
        "Figure4", "4F", MOLECULAR / "fig4f_continuum_program_trajectories.pdf",
        FIG4 / "panels/fig4f_continuum_program_trajectories.pdf",
        "eefddbd0c72ae1591a5eaf169c70a546d4e52e2c8846522e6381ba7595bae73b",
        "stage/sex-adjusted focal-program continuum trajectories",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/60_render_figures.R",
        "latest validated complete molecular-layer candidate",
    ),
    panel(
        "FigureS3", "S3H", ROOT / "figures/main/fig3_bulk_transcriptomics/panels/sample_spearman_heatmap.pdf",
        FIGS3 / "panels/figs3h_sample_spearman_heatmap.pdf",
        "2b517c026972ef81a248c9060431a0d061a0bdfb92fb2b70adef189ac8392279",
        "sample-level Spearman clustering sensitivity",
        "scripts/figures/figS_sample_spearman_heatmap.R",
        "current byte-identical canonical sample-clustering render",
    ),
    panel(
        "FigureS3", "S3I", REVIEW / "fig3_candidate_continuum_by_stage.pdf",
        FIGS3 / "panels/figs3i_continuum_calibration_by_stage.pdf",
        "a71bd46f29f5dc523189b000fef7c4762e5c985ad705eb27475e7a0436f0835a",
        "continuum calibration against recorded fibrosis stage",
        "scripts/analysis/histology_anchored_continuum/10_render_figure3_review_suite.R",
        "latest validated core-scorer review render; routed to supplement",
    ),
    panel(
        "FigureS3", "S3J", REVIEW / "figs3_candidate_score_fibrosis_anchor_forest.pdf",
        FIGS3 / "panels/figs3j_continuum_score_fibrosis_anchor_forest.pdf",
        "65e684a2762097f88d030efec30521386d8bcea1a770a4e316f0723ff2937013",
        "continuum-score fibrosis anchoring control",
        "scripts/analysis/histology_anchored_continuum/10_render_figure3_review_suite.R",
        "latest validated core-scorer review render",
    ),
    panel(
        "FigureS3", "S3K", REVIEW / "figs3_candidate_all_scores_by_stage.pdf",
        FIGS3 / "panels/figs3k_all_continuum_scores_by_stage.pdf",
        "7e62e9535ea8d21b8149e500624f4238311bf6a401603ff66aa91efd281364b2",
        "all continuum scorers across fibrosis stage",
        "scripts/analysis/histology_anchored_continuum/10_render_figure3_review_suite.R",
        "latest validated core-scorer review render",
    ),
    panel(
        "FigureS3", "S3L", REVIEW / "figs3_reference_all_scorer_agreement.pdf",
        FIGS3 / "panels/figs3l_continuum_scorer_agreement.pdf",
        "3c336bc80042d489e2fc3110406b85dfbf9460a9f191acdd320fa9f88a0174d7",
        "continuum scorer agreement",
        "scripts/analysis/histology_anchored_continuum/10_render_figure3_review_suite.R",
        "validated byte-copy from the frozen scorer candidate",
    ),
    panel(
        "FigureS3", "S3M", COMPARISON / "s3_transcript_stage_vs_continuum_effects.pdf",
        FIGS3 / "panels/figs3m_transcript_stage_vs_continuum_effects.pdf",
        "220e320877d985d5565dc8850936478ef3db1ae3eaac1a3c2ad6f381c17ca651",
        "all-transcript stage-versus-continuum effect comparison",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/62_render_stage_continuum_comparison.R",
        "latest validated matched comparison render",
    ),
    panel(
        "FigureS3", "S3N", COMPARISON / "s3_adjacent_stage_degs_vs_continuum.pdf",
        FIGS3 / "panels/figs3n_adjacent_stage_degs_vs_continuum.pdf",
        "855df85348cf36d4d40ba52c1014c9fae3b582626443ec7401d029ab15c2f9d5",
        "adjacent-stage DEG versus continuum comparison",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/62_render_stage_continuum_comparison.R",
        "latest validated matched comparison render",
    ),
    panel(
        "FigureS3", "S3O", MOLECULAR / "s3_fixed_gene_adjusted_forest.pdf",
        FIGS3 / "panels/figs3o_fixed_gene_continuum_forest.pdf",
        "007cf7e7e7050bc03574fe55f6e12cb5eb9cb3ba8e781d5deec17d9b55814401",
        "fixed five-gene continuum forest",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/60_render_figures.R",
        "latest validated complete molecular-layer candidate",
    ),
    panel(
        "FigureS3", "S3P", MOLECULAR / "s3_fixed_gene_adjusted_trajectories.pdf",
        FIGS3 / "panels/figs3p_fixed_gene_continuum_trajectories.pdf",
        "14b477605588ac8991bd07b78fbb72292d91eef0f320a46dd97208bf7191c909",
        "fixed five-gene adjusted trajectories",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/60_render_figures.R",
        "latest validated complete molecular-layer candidate",
    ),
    panel(
        "FigureS3", "S3Q", MOLECULAR / "s3_all_gene_disease_continuum_concordance.pdf",
        FIGS3 / "panels/figs3q_all_gene_disease_continuum_concordance.pdf",
        "6fe763ab60609902661d7aad833e9e6b0359b18e5d88571cffe8f025f8ea38d8",
        "all-gene disease-control versus continuum concordance",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/60_render_figures.R",
        "latest validated complete molecular-layer candidate",
    ),
    panel(
        "FigureS3", "S3R", MOLECULAR / "s3_heldout_disease_signature_validation.pdf",
        FIGS3 / "panels/figs3r_heldout_disease_signature_validation.pdf",
        "0e831694c03e45c5203854d16a431baccb5f590ca42f0d92624f14cc8a0fd2ee",
        "held-out disease-signature validation",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/60_render_figures.R",
        "latest validated complete molecular-layer candidate",
    ),
    panel(
        "FigureS3", "S3S", MOLECULAR / "s3_all50_hallmark_heatmap.pdf",
        FIGS3 / "panels/figs3s_all50_hallmark_continuum_effects.pdf",
        "fc83ee2752ee6905241205c64d3e948f51b3e76745d907320806625c661f1fce",
        "complete 50-Hallmark continuum family",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/60_render_figures.R",
        "latest validated complete molecular-layer candidate",
    ),
    panel(
        "FigureS3", "S3T", COMPARISON / "s3_six_hallmarks_by_continuum_windows.pdf",
        FIGS3 / "panels/figs3t_six_hallmarks_by_continuum_windows.pdf",
        "c25d7577fefb196391cb8b547db1821d70e5f419c48c4b94052ee52d3ad288ff",
        "six prespecified Hallmarks across continuum windows",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/62_render_stage_continuum_comparison.R",
        "latest validated matched comparison render",
    ),
    panel(
        "FigureS3", "S3U", COMPARISON / "s3_six_hallmarks_by_fibrosis_stage.pdf",
        FIGS3 / "panels/figs3u_six_hallmarks_by_fibrosis_stage.pdf",
        "6ec47c1d6fc1b0206248818460d3865581db2606d4e58880dd9371f9ac1e7890",
        "six prespecified Hallmarks across fibrosis stage",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/62_render_stage_continuum_comparison.R",
        "latest validated matched comparison render",
    ),
    panel(
        "FigureS3", "S3V", COMPARISON / "s3_pathways_stage_vs_continuum_effects.pdf",
        FIGS3 / "panels/figs3v_pathways_stage_vs_continuum_effects.pdf",
        "235a8bdc4b97eaf79c5868e2190765134267811f72c60152418ab74a72f3f6cd",
        "complete pathway-family stage-versus-continuum comparison",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/62_render_stage_continuum_comparison.R",
        "latest validated matched comparison render",
    ),
    panel(
        "FigureS3", "S3W", MOLECULAR / "s3_all_pathway_family_skylines.pdf",
        FIGS3 / "panels/figs3w_all_pathway_family_skylines.pdf",
        "9a75b5dc67f722b339fb2f55e4e3a5cd9cc3994d392ecd28288f74f1029432ce",
        "complete pathway-family continuum skylines",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/60_render_figures.R",
        "latest validated complete molecular-layer candidate",
    ),
    panel(
        "FigureS4", "S4F", COMPARISON / "s4_focal_hotspot_by_continuum_windows.pdf",
        FIG4 / "panels/supplementary/figs4f_focal_hotspot_by_continuum_windows.pdf",
        "2a64e1bb3aec98369ab20a54ec701a0606cf7299bc7c12dcae7bca488bd5db3b",
        "focal Hotspot programs across continuum windows",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/62_render_stage_continuum_comparison.R",
        "latest validated matched comparison render; supersedes older consensus-window panel",
    ),
    panel(
        "FigureS4", "S4G", COMPARISON / "s4_focal_hotspot_by_fibrosis_stage.pdf",
        FIG4 / "panels/supplementary/figs4g_focal_hotspot_by_fibrosis_stage.pdf",
        "ca1b6682081600c2747b9d76d49c273a392ea87da388e58861ab8476b2c13874",
        "focal Hotspot programs across fibrosis stage",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/62_render_stage_continuum_comparison.R",
        "latest validated matched comparison render",
    ),
    panel(
        "FigureS4", "S4H", WINDOWS / "s4_all117_hotspot_continuum_windows_labeled.pdf",
        FIG4 / "panels/supplementary/figs4h_all117_hotspot_continuum_windows_labeled.pdf",
        "b6726b02a768ddeb9b9db0e0f83712c1bcb47630f1674872aaafbe4b7af8594d",
        "all 117 labeled Hotspot programs across continuum windows",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/63_render_all_hotspot_windows.R",
        "newest validated labeled complete-family render; replaces unlabeled effects heatmap",
    ),
    panel(
        "FigureS4", "S4I", WINDOWS / "s4_continuum_associated_hotspot_windows_labeled.pdf",
        FIG4 / "panels/supplementary/figs4i_supported_hotspot_continuum_windows_labeled.pdf",
        "2e325899643cbe782577b5669c5cbd7964aca22aa9b4505652d672e6531ef77a",
        "continuum-associated Hotspot programs across continuum windows",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/63_render_all_hotspot_windows.R",
        "newest validated readable supported-program render",
    ),
    panel(
        "FigureS4", "S4J", COMPARISON / "s4_all117_hotspot_stage_vs_continuum_effects.pdf",
        FIG4 / "panels/supplementary/figs4j_all117_hotspot_stage_vs_continuum_effects.pdf",
        "d1335349522051dbdbfce2eeb8a53ec368751a4118e62057b3a9a5e553b8f170",
        "all-117 Hotspot stage-versus-continuum comparison",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/62_render_stage_continuum_comparison.R",
        "latest validated matched comparison render",
    ),
    panel(
        "FigureS4", "S4K", MOLECULAR / "s4_all10_nmf_forest.pdf",
        FIG4 / "panels/supplementary/figs4k_all10_nmf_continuum_forest.pdf",
        "6dc5eb036007615906bd6cad0f2fe05523dcd0e40dbc2701dc91ec86fdea2887",
        "all ten NMF co-primary continuum effects",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/60_render_figures.R",
        "latest validated complete molecular-layer candidate",
    ),
    panel(
        "FigureS4", "S4L", COMPARISON / "s4_all10_nmf_stage_vs_continuum_effects.pdf",
        FIG4 / "panels/supplementary/figs4l_all10_nmf_stage_vs_continuum_effects.pdf",
        "d802c4da98b1dcbfc550264b41db6fd1b6de5bb06a0d00c088fbe521d98865b9",
        "all-ten NMF stage-versus-continuum comparison",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/62_render_stage_continuum_comparison.R",
        "latest validated matched comparison render",
    ),
    panel(
        "FigureS4", "S4M", COMPARISON / "s4_all10_nmf_by_continuum_windows.pdf",
        FIG4 / "panels/supplementary/figs4m_all10_nmf_by_continuum_windows.pdf",
        "03505dd112ac95f29c1728e0f2bf647ac6ab123451bec0514cf7aa31d13e4272",
        "all ten NMF axes across continuum windows",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/62_render_stage_continuum_comparison.R",
        "latest validated matched comparison render",
    ),
    panel(
        "FigureS4", "S4N", COMPARISON / "s4_all10_nmf_by_fibrosis_stage.pdf",
        FIG4 / "panels/supplementary/figs4n_all10_nmf_by_fibrosis_stage.pdf",
        "599301edfa21e5e7ea11c43dfc4e87c82da8d105766235fa7067dc8da68bd623",
        "all ten NMF axes across fibrosis stage",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/62_render_stage_continuum_comparison.R",
        "latest validated matched comparison render",
    ),
    panel(
        "FigureS4", "S4O", MOLECULAR / "s4_four_tf_regulon_trajectories.pdf",
        FIG4 / "panels/supplementary/figs4o_four_tf_regulon_continuum_windows.pdf",
        "d312d21961e6780719e163c655fa2b6f7a3240cd5af65d25067bc054025de6fd",
        "four fixed regulon trajectories across continuum windows",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/60_render_figures.R",
        "latest validated complete molecular-layer candidate",
    ),
    panel(
        "FigureS4", "S4P", MOLECULAR / "s4_paired_focal_programs.pdf",
        FIG4 / "panels/supplementary/figs4p_paired_focal_continuum_programs.pdf",
        "35731f1850045a2cb0ca55bf485774c662cce5ded58404e45379f1bbb2a62ee0",
        "paired GSE193066 focal-program continuum sensitivity",
        "scripts/analysis/histology_anchored_continuum/molecular_layers/60_render_figures.R",
        "latest validated complete molecular-layer candidate",
    ),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_verified(item: Panel) -> dict[str, str]:
    if not item.source.is_file():
        raise RuntimeError(f"Missing source: {item.source}")
    observed = sha256(item.source)
    if observed != item.expected_sha256:
        raise RuntimeError(
            f"Source hash drift for {item.callout}: {observed} != {item.expected_sha256}"
        )
    item.destination.parent.mkdir(parents=True, exist_ok=True)
    copy_state = "copied"
    if item.destination.exists():
        if sha256(item.destination) != observed:
            raise RuntimeError(f"Refusing to overwrite differing destination: {item.destination}")
        copy_state = "already_identical"
    else:
        temporary = item.destination.with_name(f".{item.destination.name}.tmp-{os.getpid()}")
        if temporary.exists():
            raise RuntimeError(f"Temporary path already exists: {temporary}")
        shutil.copy2(item.source, temporary)
        if sha256(temporary) != observed:
            raise RuntimeError(f"Copy verification failed: {temporary}")
        temporary.replace(item.destination)
    return {
        "authority": item.authority,
        "callout": item.callout,
        "destination": str(item.destination.relative_to(ROOT)),
        "source": str(item.source.relative_to(ROOT)),
        "source_sha256": observed,
        "destination_sha256": sha256(item.destination),
        "role": item.role,
        "generator": item.generator,
        "selection_reason": item.selection_reason,
        "copy_state": copy_state,
        "release_state": "promoted_working_panel_synchronized_release_pending",
    }


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite manifest: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    rows = [copy_verified(item) for item in PANELS]
    s3_rows = [row for row in rows if row["authority"] == "FigureS3"]
    fig4_rows = [row for row in rows if row["authority"] != "FigureS3"]
    write_manifest(
        FIGS3 / "manifests/continuum_panel_promotion_20260824.tsv", s3_rows
    )
    write_manifest(
        FIG4 / "manifests/continuum_panel_promotion_20260824.tsv", fig4_rows
    )
    print(f"PROMOTED\t{len(rows)}")
    print(f"FIGURE4_MAIN\t{sum(row['authority'] == 'Figure4' for row in rows)}")
    print(f"FIGURES3\t{len(s3_rows)}")
    print(f"FIGURES4\t{sum(row['authority'] == 'FigureS4' for row in rows)}")


if __name__ == "__main__":
    main()
