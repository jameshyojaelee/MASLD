#!/usr/bin/env python3
"""Convert Cell Ranger h5 outputs to per-sample h5ad files for ScaleSC.

Creates one h5ad per sample with metadata columns:
  - sample: SRR accession
  - dataset: GSE ID
  - species: human/mouse
  - condition: parsed from sample titles where available

Output: {output_dir}/human/ and {output_dir}/mouse/ directories of h5ad files.
"""

import csv
import logging
import os
import sys
import time
from pathlib import Path

import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTPUT_DIR = BASE / "Analysis/SingleCell/integration/input_h5ad"


def parse_condition_gse174748():
    """GSE174748: Healthy vs NAFLD from sample titles."""
    mapping = {}
    runinfo = BASE / "data/GSE174748/metadata/ena_runinfo.tsv"
    with runinfo.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            title = row.get("sample_title", "")
            srr = row["run_accession"]
            if "Healthy" in title:
                mapping[srr] = "Healthy"
            elif "NAFLD" in title:
                mapping[srr] = "NAFLD"
            else:
                mapping[srr] = "Unknown"
    return mapping


def parse_condition_gse189600():
    """GSE189600: Parse condition from sample titles."""
    mapping = {}
    runinfo = BASE / "data/GSE189600/metadata/ena_runinfo.tsv"
    with runinfo.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            title = row.get("sample_title", "")
            srr = row["run_accession"]
            if "healthy" in title.lower():
                mapping[srr] = "Healthy"
            elif "NASH" in title:
                mapping[srr] = "NASH"
            elif "NC_" in title:
                mapping[srr] = "Normal_Chow"
            elif "ALIOS" in title:
                mapping[srr] = "ALIOS"
            elif "GFP" in title or "Ephb2" in title:
                mapping[srr] = "Bulk_RNA"  # skip
            else:
                mapping[srr] = "Unknown"
    return mapping


def parse_condition_gse244832():
    """GSE244832: Sample titles are just IDs (JB_288, MM_456). Need GEO supp."""
    return {}  # Will default to "Unknown" - need GEO metadata


def parse_condition_gse185477():
    """GSE185477: Single cell/nuclei donors - healthy liver."""
    return {}  # All healthy liver donors


def collect_samples():
    """Collect all Cell Ranger h5 files with metadata."""
    samples = []

    # --- Liver Atlas (human) ---
    la_human = BASE / "Liver_Atlas/cellranger/human"
    for h5 in sorted(la_human.glob("*/outs/filtered_feature_bc_matrix.h5")):
        srr = h5.parent.parent.name
        samples.append({
            "srr": srr, "dataset": "Liver_Atlas", "species": "human",
            "condition": "Mixed", "h5_path": str(h5),
        })

    # --- Liver Atlas (mouse) ---
    la_mouse = BASE / "Liver_Atlas/cellranger/mouse"
    for h5 in sorted(la_mouse.glob("*/outs/filtered_feature_bc_matrix.h5")):
        srr = h5.parent.parent.name
        samples.append({
            "srr": srr, "dataset": "Liver_Atlas", "species": "mouse",
            "condition": "Mixed", "h5_path": str(h5),
        })

    # --- GSE202379 (human, MASLD cohort) ---
    cr_dir = BASE / "data/GSE202379/cellranger"
    for h5 in sorted(cr_dir.glob("*/outs/filtered_feature_bc_matrix.h5")):
        srr = h5.parent.parent.name
        samples.append({
            "srr": srr, "dataset": "GSE202379", "species": "human",
            "condition": "MASLD", "h5_path": str(h5),
        })

    # --- GSE244832 (human, MASLD spectrum + HSC focus) ---
    cond_map = parse_condition_gse244832()
    cr_dir = BASE / "data/GSE244832/cellranger"
    for h5 in sorted(cr_dir.glob("*/outs/filtered_feature_bc_matrix.h5")):
        srr = h5.parent.parent.name
        samples.append({
            "srr": srr, "dataset": "GSE244832", "species": "human",
            "condition": cond_map.get(srr, "MASLD"), "h5_path": str(h5),
        })

    # --- GSE174748 (human, Healthy + NAFLD) ---
    cond_map = parse_condition_gse174748()
    cr_dir = BASE / "data/GSE174748/cellranger"
    for h5 in sorted(cr_dir.glob("*/outs/filtered_feature_bc_matrix.h5")):
        srr = h5.parent.parent.name
        samples.append({
            "srr": srr, "dataset": "GSE174748", "species": "human",
            "condition": cond_map.get(srr, "Unknown"), "h5_path": str(h5),
        })

    # --- GSE185477 (human, healthy liver donors) ---
    cr_dir = BASE / "data/GSE185477/cellranger"
    for h5 in sorted(cr_dir.glob("*/outs/filtered_feature_bc_matrix.h5")):
        srr = h5.parent.parent.name
        samples.append({
            "srr": srr, "dataset": "GSE185477", "species": "human",
            "condition": "Healthy", "h5_path": str(h5),
        })

    # --- GSE189600 (mixed human + mouse) ---
    cond_map = parse_condition_gse189600()
    # Human snRNA
    cr_human = BASE / "data/GSE189600/cellranger/human"
    for h5 in sorted(cr_human.glob("*/outs/filtered_feature_bc_matrix.h5")):
        srr = h5.parent.parent.name
        samples.append({
            "srr": srr, "dataset": "GSE189600", "species": "human",
            "condition": cond_map.get(srr, "Unknown"), "h5_path": str(h5),
        })
    # Mouse snRNA
    cr_mouse = BASE / "data/GSE189600/cellranger/mouse"
    for h5 in sorted(cr_mouse.glob("*/outs/filtered_feature_bc_matrix.h5")):
        srr = h5.parent.parent.name
        samples.append({
            "srr": srr, "dataset": "GSE189600", "species": "mouse",
            "condition": cond_map.get(srr, "Unknown"), "h5_path": str(h5),
        })

    return samples


def convert_h5_to_h5ad(sample: dict, output_dir: Path) -> str:
    """Convert a Cell Ranger h5 to an annotated h5ad file."""
    species_dir = output_dir / sample["species"]
    species_dir.mkdir(parents=True, exist_ok=True)
    out_path = species_dir / f"{sample['srr']}.h5ad"

    if out_path.exists():
        return f"  {sample['srr']}: exists, skipping"

    t0 = time.time()
    adata = sc.read_10x_h5(sample["h5_path"])
    adata.var_names_make_unique()

    # Add metadata
    adata.obs["sample"] = sample["srr"]
    adata.obs["dataset"] = sample["dataset"]
    adata.obs["species"] = sample["species"]
    adata.obs["condition"] = sample["condition"]

    # Store raw counts
    adata.write_h5ad(out_path)
    dt = time.time() - t0
    return f"  {sample['srr']}: {adata.n_obs:,} cells ({dt:.1f}s)"


def main():
    t_start = time.time()
    samples = collect_samples()

    human = [s for s in samples if s["species"] == "human"]
    mouse = [s for s in samples if s["species"] == "mouse"]
    log.info(f"Found {len(samples)} samples: {len(human)} human, {len(mouse)} mouse")

    # Summary by dataset
    from collections import Counter
    ds_counts = Counter(s["dataset"] for s in samples)
    for ds, n in sorted(ds_counts.items()):
        sp = Counter(s["species"] for s in samples if s["dataset"] == ds)
        log.info(f"  {ds}: {n} samples ({dict(sp)})")

    # Convert
    log.info(f"\nConverting to h5ad → {OUTPUT_DIR}")
    for i, sample in enumerate(samples, 1):
        result = convert_h5_to_h5ad(sample, OUTPUT_DIR)
        if i % 20 == 0 or i == len(samples):
            log.info(f"  Progress: {i}/{len(samples)}")
        log.info(result)

    # Write sample manifest
    manifest_path = OUTPUT_DIR / "sample_manifest.csv"
    with manifest_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["srr", "dataset", "species", "condition", "h5_path"])
        writer.writeheader()
        writer.writerows(samples)
    log.info(f"Manifest: {manifest_path}")

    dt = time.time() - t_start
    n_human = len(list((OUTPUT_DIR / "human").glob("*.h5ad"))) if (OUTPUT_DIR / "human").exists() else 0
    n_mouse = len(list((OUTPUT_DIR / "mouse").glob("*.h5ad"))) if (OUTPUT_DIR / "mouse").exists() else 0
    log.info(f"\nDone in {dt:.0f}s: {n_human} human + {n_mouse} mouse h5ad files")


if __name__ == "__main__":
    main()
