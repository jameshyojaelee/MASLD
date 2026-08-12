#!/usr/bin/env python3
"""Static portal regression checks for categorical spatial evidence semantics."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    replay = (ROOT / "masld-atlas-v2/src/components/evidence-replay.tsx").read_text(encoding="utf-8")
    colors = (ROOT / "masld-atlas-v2/src/lib/colors.ts").read_text(encoding="utf-8")
    palette = (ROOT / "masld-atlas-v2/src/lib/palette.ts").read_text(encoding="utf-8")
    gene_data = (ROOT / "masld-atlas-v2/src/app/gene/use-gene-data.ts").read_text(encoding="utf-8")
    preprocessing = (ROOT / "masld-atlas-v2/scripts/preprocess_atlas_data.py").read_text(encoding="utf-8")
    downloads = (ROOT / "masld-atlas-v2/src/app/downloads/page.tsx").read_text(encoding="utf-8")
    spatial_page = (ROOT / "masld-atlas-v2/src/app/single-cell/page.tsx").read_text(encoding="utf-8")
    gene_sections = (ROOT / "masld-atlas-v2/src/app/gene/gene-sections.tsx").read_text(encoding="utf-8")
    programs_page = (ROOT / "masld-atlas-v2/src/app/programs/programs-client.tsx").read_text(encoding="utf-8")
    drugs_page = (ROOT / "masld-atlas-v2/src/app/drugs/page.tsx").read_text(encoding="utf-8")
    modality_matrix = (ROOT / "masld-atlas-v2/src/components/convergence/modality-matrix.tsx").read_text(encoding="utf-8")
    gene_view = (ROOT / "masld-atlas-v2/src/app/gene/gene-view.tsx").read_text(encoding="utf-8")
    explore_page = (ROOT / "masld-atlas-v2/src/app/explore/page.tsx").read_text(encoding="utf-8")

    assert "Not spatially variable" not in replay
    assert "Spatial score:" not in replay
    assert 'key: "s5_spatial"' not in colors
    numeric_modalities = palette.split("export const MODALITIES", 1)[1].split("];", 1)[0]
    assert 'key: "s5_spatial"' not in numeric_modalities
    assert "s.s5_spatial =" not in gene_data
    assert 'evidence["s5_spatial"] = 0.0' in preprocessing
    assert '"s5_spatial", "s6_singlecell"' not in preprocessing.split('atlas["layers_active"] =', 1)[1].split(".sum(axis=1)", 1)[0]
    assert "150 SVGs" not in downloads
    assert "spatial_gene_context.parquet" in downloads
    assert "spatial_program_coverage.parquet" in downloads
    assert "spatial_program_effects.parquet" in downloads
    assert "spatial_dataset_registry.json" in spatial_page
    assert "biological n" in spatial_page and "technical n" in spatial_page
    assert "Missing coverage renders untestable or" in spatial_page
    assert "spatial_gene_context.parquet" in gene_sections
    assert "Missing coverage is untestable or indeterminate" in gene_sections
    assert "spatial_program_coverage.parquet" in programs_page
    assert "spatial_program_effects.parquet" in programs_page
    assert "Coverage only. This program was not part" in programs_page
    assert 'col: "Spatial"' not in drugs_page
    assert "ORDER BY layers_active" not in drugs_page
    assert 'key: "spatial"' not in modality_matrix
    assert "displayedCount" in modality_matrix
    assert 'label="Convergence rank"' not in gene_view
    assert 'label="Non-spatial layers"' in gene_view
    assert 'key: "layers_active"' not in explore_page
    print("PASS portal spatial semantics")


if __name__ == "__main__":
    main()
