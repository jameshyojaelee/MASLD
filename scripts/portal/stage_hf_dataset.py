#!/usr/bin/env python3
"""
Phase-7 packaging for the MASLD atlas web refresh (data contract §7).

LEGACY PRE-RESOURCE PRODUCER. It is blocked by default because the staged
rank-first schema is not the MASLD Gene Catalog contract. Set
ALLOW_LEGACY_PORTAL_REBUILD=true only for provenance-only regeneration; never
use that opt-in as publication or deployment authorization.

Assembles a Hugging Face **Dataset** staging folder from the built web-data dir
(`masld-atlas-v2/public/data`) using a strict ALLOW-LIST: loose parquet (→ Git
LFS) + compact JSON + the `network/portal_export_v2/` tree (minus its
`gene_graphs/`). Because staging is allow-list-driven, the heavy retired
artifacts are never copied:
  - `genes/`                             (33,943 per-gene JSONs — retired)
  - `network/portal_export_v2/gene_graphs/`  (16,601 per-gene graph JSONs — retired)
  - any `data.tar.gz`                    (old monolith)
  - anything not explicitly on the keep-list (e.g. `bayesian_ranking.json`)

Also writes the HF dataset card (`README.md`) and `.gitattributes`, then PRINTS
the exact `hf upload` command for the USER to run. This script does NO network
I/O and NO push (PUSH is out of scope; HF auth + egress are user-run).

Env-agnostic (stdlib only). Honors MASLD_PROJECT_ROOT.

Usage:
  python scripts/portal/stage_hf_dataset.py \
      [--web-dir masld-atlas-v2/public/data] \
      [--stage-dir masld-atlas-v2/hf_dataset_stage]

Keep-list files that don't exist yet are WARNED and skipped (some are produced
by sibling generators; this step runs last in the orchestrator).
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from datetime import date
from pathlib import Path

PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)

HF_REPO = "jameshyojaelee/masld-atlas-data"

# ── Keep-list (contract §7) ────────────────────────────────────────────────
# Loose parquet → tracked by Git LFS (see .gitattributes below).
KEEP_PARQUET = [
    "atlas.parquet",
    "atlas_core.parquet",
    "drug_targets.parquet",
    "coloc_by_gwas.parquet",
    "sc_pseudobulk_de.parquet",
    "sc_hep_markers.parquet",
    "spatial_zonation.parquet",
    "program_gene_membership.parquet",
    "sex_gene_classification.parquet",
    "gene_per_cohort_de.parquet",
    "gene_trajectories.parquet",
    "gene_pseudobulk_de.parquet",
    "gene_lincs.parquet",
    "gene_drugs.parquet",
    "network_edges.parquet",
    "network_nodes.parquet",
    "sc_umap_downsampled.parquet",  # optional; warn-and-skip if absent
]

# Small JSON → plain git (no LFS).
KEEP_JSON = [
    "gene_symbols.json",
    "gene_index.json",
    "atlas_summary.json",
    "featured_genes.json",
    "convergence_ranking.json",
    "convergence_matrix.json",
    "bayesian_ranking.json",
    "proteomics_summary.json",
    "progression_journey.json",
    "drug_pipeline.json",
    "gwas_atac_browser.json",
    "knowledge_graph.json",
    "pathway_genesets.json",
    "singlecell_summary.json",
    "spatial_summary.json",
    "coloc_ancestry_summary.json",
    "programs_summary.json",
    "sex_summary.json",
    "cross_species.json",
    "volcano_preview.json",
    "umap_thumbnail.json",
    "landing_ticker.json",
]

# Directory tree copied wholesale EXCEPT the named subdir(s).
NETWORK_TREE_REL = "network/portal_export_v2"
NETWORK_TREE_EXCLUDE = {"gene_graphs"}

GITATTRIBUTES = (
    "*.parquet filter=lfs diff=lfs merge=lfs -text\n"
    "*.tar.gz filter=lfs diff=lfs merge=lfs -text\n"
)


def _human_size(n: int) -> str:
    x = float(n)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if x < 1024 or unit == "TB":
            return f"{x:.1f} {unit}"
        x /= 1024
    return f"{x:.1f} TB"


def _copy_flat(web_dir: Path, stage_dir: Path, names: list[str]) -> tuple[list[str], list[str]]:
    """Copy allow-listed loose files; return (copied, missing)."""
    copied, missing = [], []
    for name in names:
        src = web_dir / name
        if src.is_file():
            shutil.copy2(src, stage_dir / name)
            copied.append(name)
        else:
            missing.append(name)
    return copied, missing


def _copy_network_tree(web_dir: Path, stage_dir: Path) -> tuple[int, bool]:
    """Copy network/portal_export_v2 minus gene_graphs/. Return (n_files, present)."""
    src = web_dir / NETWORK_TREE_REL
    if not src.is_dir():
        return 0, False
    dst = stage_dir / NETWORK_TREE_REL
    shutil.copytree(
        src,
        dst,
        ignore=shutil.ignore_patterns(*NETWORK_TREE_EXCLUDE),
        dirs_exist_ok=True,
    )
    n = sum(1 for p in dst.rglob("*") if p.is_file())
    return n, True


def _readme_card() -> str:
    """HF dataset card: YAML front-matter (viewer → atlas_core.parquet) + body."""
    gen = date.today().isoformat()
    front = f"""---
license: cc-by-4.0
pretty_name: MASLD Gene Catalog
tags:
  - genomics
  - MASLD
  - liver
  - transcriptomics
configs:
  - config_name: atlas_core
    data_files:
      - split: train
        path: atlas_core.parquet
---
"""
    body = f"""
# MASLD Gene Catalog web data

Backing data for the MASLD atlas web portal (Hugging Face Static Space). A
**human-only** metabolic-dysfunction-associated steatotic liver disease (MASLD)
transcriptomic and multi-omic MASLD Gene Catalog. Disease effects are framed in human
tissue only; **cross-ancestry** genetics (EUR / AFR / EAS / AMR / SAS) is
first-class. Mouse content is retained solely for the Cas13 perturbation-library
design and is never used to make disease claims.

Generated: {gen}. Column names in the parquet files are data (kept verbatim);
UI labels never reproduce internal statistical shorthand.

## Provenance (canonical sources + dates)

| Layer | Source (project-relative) | Canonical date | Headline |
|---|---|---|---|
| Bulk RNA-seq DEG | `canonical_deg_results.csv` (pooled limma-voom-qw, C2 design) | 2026-06-29 | 1,918 DEGs via the effect-size-aware interval-null FDR gate (`fdr < 0.05` at `lfc = 0.25`, McCarthy & Smyth 2009) |
| MASLD Gene Catalog | `multi_evidence_atlas.csv` | 2026-07-06 | 27,187 genes × 443 evidence columns |
| Convergence ranking | `convergence_evidence.csv` (46d) | 2026-07-06 | Tier-1 = 677 (evidence-weighted rank, not a posterior) |
| Genetics / COLOC | `susie_coloc/` (SuSiE + ABF, PolyFun EUR LD) | 2026-07-06 | 50-GWAS / 5-ancestry portfolio; 473 SuSiE / 1,031 union effector genes (main Tier-1/2) |
| Single-cell | `Analysis/SingleCell/results_gpu_v2/` | 2026 | ~1.23M cells; 16 cell types; pseudobulk + hepatocyte subtypes |
| Spatial | `Analysis/Spatial/results/` | 2026 | Visium + GeoMx + CosMx zonation, SVGs, zonation disruption |
| Molecular programs | `subtypes/` (NMF k=6, P1–P6) + Hotspot modules | 2026 | 6 continuous program axes (not discrete subtypes); autocorrelation modules |
| Sex | `sex_v3/` (LVQW-fixed, C2) | 2026-06-08 | 8 sex-dimorphic genes (7 male-biased / 1 female-biased) |
| Proteomics | DIA-MS + Olink | 2026 | plasma/tissue mRNA↔protein concordance |
| Drug repurposing | `drug_repurposing/` (LINCS/CGP + clinical) | 2026 | 2 approved / 208 clinical / 1,504 preclinical target-linked |
| Network | `network/edge_annotation_atlas.parquet` | 2026 | filtered gene–gene multi-channel edges + nodes |

## Primary table — `atlas_core.parquet`

The client's main table (~42 columns), keyed by `human_symbol`. The dataset
viewer above is configured to preview it. Column families:

- **Identity**: `human_symbol`, `ensembl_id`, `gene_biotype`, `mouse_ortholog`.
- **Bulk DE (pooled C2 limma-voom-qw)**: `bulk_logFC`, `bulk_padj`, `bulk_tstat`,
  `bulk_shrunk_logFC`, `bulk_lfsr`, `bulk_sig`, `is_deg` (canonical 2026-08-12 gate:
  padj<0.05 and |log2FC|>0.50; `bulk_sig` in the stored atlas is still TREAT-era),
  plus interval-null estimate columns and sex-stratified `bulk_logFC_M/F`,
  `sex_class`, `sex_interaction_padj`.
- **Genetics / COLOC**: `coloc_best_susie_pp4`, `coloc_best_susie_gwas`,
  `coloc_abf_best_pp4`, `coloc_abf_best_gwas`, `n_coloc_sources`,
  `n_ancestry_gwas`, `coloc_cross_ancestry_replicated`, `coloc_susie_conf_tier`,
  `twas_z`, `twas_pval`.
- **Convergence**: `convergence_rank`, `convergence_score`, `convergence_tier`,
  `concordance_state`.
- **Essentiality / spatial / zonation / ferroptosis**: `essentiality_chronos`,
  `is_essential`, `spatial_is_svg`, `spatial_morans_i`,
  `spatial_consensus_direction`, `zonation_class`, `ferroptosis_class`.
- **Druggability / programs**: `dgidb_druggable`, `opentargets_drug`,
  `max_phase_masld`, `drug_dev_status`, `pharos_tdl`,
  `dominant_program_for_gene`, `dominant_program_logFC`, `layers_active`,
  `is_conserved`.

The full 443-column `atlas.parquet` is retained as a bulk download artifact.
Long-format `gene_*` / `coloc_by_gwas` / `sc_*` / `spatial_*` /
`program_gene_membership` / `sex_gene_classification` parquets are sorted by
`symbol` for `WHERE symbol = ?` gene-detail lookups.

## How the web app consumes this

The portal is a static Next.js export on a Hugging Face **Static Space**. It
fetches every file from this Dataset's `resolve/main/` CDN via a whole-file
`GET` (`DATA_BASE = https://huggingface.co/datasets/{HF_REPO}/resolve/main`),
registers the bytes into **DuckDB-WASM**, and queries client-side. HF's
`resolve` endpoint serves `Access-Control-Allow-Origin: *`, so cross-origin
whole-file reads work without a proxy. Parquet is tracked by Git LFS; small
JSON is plain git.

## Files

- Loose parquet (LFS): atlas, atlas_core, coloc_by_gwas, gene_* detail tables,
  sc_*, spatial_zonation, program_gene_membership, sex_gene_classification,
  network_edges, network_nodes, (optional) sc_umap_downsampled.
- Compact JSON: per-layer summaries + rankings + gene_symbols (fuse.js index).
- `network/portal_export_v2/`: community + layer metadata + filter indexes
  (per-gene graph JSONs are intentionally excluded — the client reads
  `network_edges.parquet` instead).
"""
    return front + body


def main() -> int:
    if os.environ.get("ALLOW_LEGACY_PORTAL_REBUILD", "false").lower() != "true":
        raise SystemExit(
            "REFUSED: legacy HF staging is outside the standalone Resource "
            "contract. Use Plans 50/60."
        )
    ap = argparse.ArgumentParser()
    ap.add_argument("--web-dir", default="masld-atlas-v2/public/data")
    ap.add_argument("--stage-dir", default="masld-atlas-v2/hf_dataset_stage")
    args = ap.parse_args()

    web_dir = Path(args.web_dir)
    stage_dir = Path(args.stage_dir)
    if not web_dir.is_absolute():
        web_dir = PROJECT_ROOT / web_dir
    if not stage_dir.is_absolute():
        stage_dir = PROJECT_ROOT / stage_dir

    if not web_dir.is_dir():
        print(f"ERROR: web-dir not found: {web_dir}", file=sys.stderr)
        return 1

    # Safety: refuse to wipe a real git clone; otherwise recreate for idempotency.
    if stage_dir.exists():
        if (stage_dir / ".git").exists():
            print(
                f"ERROR: stage-dir looks like a git clone (has .git): {stage_dir}\n"
                "Refusing to wipe. Point --stage-dir at a scratch folder.",
                file=sys.stderr,
            )
            return 1
        shutil.rmtree(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)

    print(f"web-dir  : {web_dir}")
    print(f"stage-dir: {stage_dir}\n")

    pq_copied, pq_missing = _copy_flat(web_dir, stage_dir, KEEP_PARQUET)
    js_copied, js_missing = _copy_flat(web_dir, stage_dir, KEEP_JSON)
    net_files, net_present = _copy_network_tree(web_dir, stage_dir)

    (stage_dir / ".gitattributes").write_text(GITATTRIBUTES)
    (stage_dir / "README.md").write_text(_readme_card())

    # ── Report ────────────────────────────────────────────────────────────
    print(f"Parquet copied ({len(pq_copied)}/{len(KEEP_PARQUET)}): {pq_copied}")
    print(f"JSON copied    ({len(js_copied)}/{len(KEEP_JSON)}): {js_copied}")
    print(f"network/portal_export_v2 tree: "
          f"{'copied ' + str(net_files) + ' files (gene_graphs/ excluded)' if net_present else 'NOT PRESENT'}")

    warned = pq_missing + js_missing + ([] if net_present else [NETWORK_TREE_REL])
    if warned:
        print("\nWARN — keep-list entries not present (skipped; produced by sibling generators):")
        for w in warned:
            print(f"  - {w}")

    total_bytes = sum(p.stat().st_size for p in stage_dir.rglob("*") if p.is_file())
    total_files = sum(1 for p in stage_dir.rglob("*") if p.is_file())
    print(f"\nStaged: {total_files} files, {_human_size(total_bytes)} total")

    commit_msg = f"Refresh MASLD atlas web data {date.today().isoformat()}"
    print("\n" + "=" * 72)
    print("USER-RUN upload command (no push performed by this script):")
    print("=" * 72)
    print(
        f'hf upload {HF_REPO} "{stage_dir}" . '
        f'--repo-type=dataset --commit-message="{commit_msg}"'
    )
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
