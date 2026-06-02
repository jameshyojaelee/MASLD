#!/usr/bin/env python3
"""Organize HMSMA processed Visium data into scanpy.read_visium()-compatible structure.

HMSMA naming: {CONDITION}-{PATIENT_ID}_{filetype}.{ext}
  e.g., CTRL-161_filtered_feature_bc_matrix.h5
        CTRL-161_tissue_positions.csv

Target structure:
  results/spaceranger/HRA007511_HMSMA/{SAMPLE}/outs/
    ├── filtered_feature_bc_matrix.h5
    └── spatial/
        ├── tissue_positions_list.csv   (renamed from tissue_positions.csv)
        ├── scalefactors_json.json
        ├── tissue_hires_image.png
        └── tissue_lowres_image.png

Key differences from GSE192741:
  - HMSMA uses tissue_positions.csv (SpaceRanger v2 format)
  - GSE192741 uses tissue_positions_list.csv (SpaceRanger v1 format)
  - scanpy ≥1.9 handles both, but we symlink for consistency
  - HMSMA files are not compressed (.gz)
"""

import pathlib
import shutil
import sys
import csv
import json
import argparse

PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HMSMA_DATA = PROJECT_ROOT / "data" / "HRA007511" / "hmsma_processed" / "stomics"
HE_DATA = PROJECT_ROOT / "data" / "HRA007511" / "hmsma_processed" / "he_images"
SR_OUTPUT = PROJECT_ROOT / "Analysis" / "Spatial" / "results" / "spaceranger" / "HRA007511_HMSMA"
METADATA_DIR = PROJECT_ROOT / "Analysis" / "Spatial" / "metadata"
HISTOLOGY_DIR = METADATA_DIR / "histology" / "HRA007511_HMSMA"

# ─── Sample Registry ────────────────────────────────────────────────────────
# Disease condition mapping (extracted from HMSMA naming convention)
CONDITION_MAP = {
    "CTRL": "Healthy",
    "MASLD": "Steatotic",  # MASL in H&E naming = MASLD in Stomics naming
    "MASH": "NASH",
}

# All 35 Visium samples
SAMPLES = [
    # CTRL (7)
    "CTRL-161", "CTRL-179", "CTRL-180", "CTRL-5113", "CTRL-5759", "CTRL-8715", "CTRL-8395",
    # MASLD (10)
    "MASLD-0966", "MASLD-1479", "MASLD-1492", "MASLD-1493", "MASLD-1495",
    "MASLD-1497", "MASLD-1498", "MASLD-2768", "MASLD-4973", "MASLD-9993",
    # MASH (18)
    "MASH-0413", "MASH-0422", "MASH-0835", "MASH-1086", "MASH-1475",
    "MASH-1478", "MASH-1480", "MASH-1481", "MASH-1494", "MASH-1501",
    "MASH-2534", "MASH-3096", "MASH-3344", "MASH-4426", "MASH-7866",
    "MASH-8684", "MASH-9136", "MASH-9440",
]

# File type → (target subdir relative to outs/, target filename)
FILE_MAP = {
    "filtered_feature_bc_matrix.h5": (".", "filtered_feature_bc_matrix.h5"),
    "scalefactors_json.json": ("spatial", "scalefactors_json.json"),
    "tissue_lowres_image.png": ("spatial", "tissue_lowres_image.png"),
    "tissue_hires_image.png": ("spatial", "tissue_hires_image.png"),
    # tissue_positions.csv → tissue_positions_list.csv for SpaceRanger v1 compat
    "tissue_positions.csv": ("spatial", "tissue_positions_list.csv"),
}


def get_condition(sample_id: str) -> str:
    """Extract disease condition from sample ID prefix."""
    prefix = sample_id.split("-")[0]
    return CONDITION_MAP.get(prefix, "unknown")


def get_individual(sample_id: str) -> str:
    """Extract patient ID from sample ID."""
    return sample_id.split("-")[1]


def fix_tissue_positions(src: pathlib.Path, dst: pathlib.Path):
    """Convert SpaceRanger v2 tissue_positions.csv to v1 tissue_positions_list.csv.

    v2 format: barcode,in_tissue,array_row,array_col,pxl_col_in_fullres,pxl_row_in_fullres (with header)
    v1 format: barcode,in_tissue,array_row,array_col,pxl_col_in_fullres,pxl_row_in_fullres (no header)

    scanpy.read_visium() auto-detects format, but we ensure consistency.
    """
    with open(src) as f:
        reader = csv.reader(f)
        rows = list(reader)

    # Check if first row is a header (non-numeric second column)
    has_header = False
    if rows and not rows[0][1].isdigit():
        has_header = True

    with open(dst, "w", newline="") as f:
        writer = csv.writer(f)
        start = 1 if has_header else 0
        for row in rows[start:]:
            writer.writerow(row)


def organize_sample(sample_id: str, dry_run: bool = False) -> dict:
    """Organize one sample's files into SpaceRanger directory structure."""
    outs_dir = SR_OUTPUT / sample_id / "outs"
    spatial_dir = outs_dir / "spatial"

    status = {"sample": sample_id, "condition": get_condition(sample_id)}
    missing = []
    copied = []

    if not dry_run:
        outs_dir.mkdir(parents=True, exist_ok=True)
        spatial_dir.mkdir(parents=True, exist_ok=True)

    for src_suffix, (subdir, target_name) in FILE_MAP.items():
        src = HMSMA_DATA / f"{sample_id}_{src_suffix}"
        target_dir = outs_dir / subdir if subdir != "." else outs_dir
        target = target_dir / target_name

        if not src.exists():
            missing.append(src_suffix)
            continue

        if target.exists():
            copied.append(f"{target_name} (exists)")
            continue

        if dry_run:
            copied.append(f"{target_name} (would copy)")
            continue

        # Special handling for tissue_positions
        if "tissue_positions" in src_suffix:
            fix_tissue_positions(src, target)
            copied.append(f"{target_name} (converted)")
        else:
            shutil.copy2(src, target)
            copied.append(f"{target_name} (copied)")

    status["missing"] = missing
    status["files"] = copied
    status["complete"] = len(missing) == 0
    return status


def organize_he_images(dry_run: bool = False):
    """Copy H&E images to histology directory for figures."""
    if not dry_run:
        HISTOLOGY_DIR.mkdir(parents=True, exist_ok=True)

    n_ok = 0
    for sample_id in SAMPLES:
        # MASLD samples use MASL prefix in H&E naming
        he_name = sample_id
        if sample_id.startswith("MASLD-"):
            he_name = "MASL-" + sample_id[len("MASLD-"):]

        # Try both .tif and .jpg
        src = None
        for ext in [".tif", ".jpg"]:
            candidate = HE_DATA / f"{he_name}{ext}"
            if candidate.exists():
                src = candidate
                break

        if src is None:
            continue

        dst = HISTOLOGY_DIR / f"{sample_id}{src.suffix}"
        if dst.exists() or dry_run:
            n_ok += 1
            continue

        shutil.copy2(src, dst)
        n_ok += 1

    return n_ok


def write_sample_list():
    """Write sample list for pipeline use."""
    samples_file = METADATA_DIR / "HRA007511_HMSMA_samples.txt"
    samples_file.write_text("\n".join(sorted(SAMPLES)) + "\n")
    print(f"\nSample list: {samples_file} ({len(SAMPLES)} samples)")


def main():
    parser = argparse.ArgumentParser(
        description="Organize HMSMA data into SpaceRanger-compatible structure"
    )
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be done without copying")
    parser.add_argument("--validate-only", action="store_true",
                        help="Only check if target structure is complete")
    args = parser.parse_args()

    print("=" * 70)
    print("00h: Organize HMSMA Processed Visium Data")
    print(f"  Source: {HMSMA_DATA}")
    print(f"  Target: {SR_OUTPUT}")
    print(f"  Samples: {len(SAMPLES)} (7 CTRL, 10 MASLD, 18 MASH)")
    if args.dry_run:
        print("  MODE: DRY RUN")
    print("=" * 70)

    # Check source directory
    if not HMSMA_DATA.exists():
        print(f"\nERROR: Source directory not found: {HMSMA_DATA}")
        print("Please download HMSMA data first (see 00g_download_hmsma_processed.sh)")
        sys.exit(1)

    src_files = list(HMSMA_DATA.glob("*_filtered_feature_bc_matrix.h5"))
    print(f"\nFound {len(src_files)} h5 files in source directory")

    if not src_files and not args.validate_only:
        print("WARNING: No h5 files found. Is the data downloaded?")

    # Organize each sample
    results = []
    n_complete = 0
    for i, sample_id in enumerate(sorted(SAMPLES)):
        status = organize_sample(sample_id, dry_run=args.dry_run or args.validate_only)
        results.append(status)
        if status["complete"]:
            n_complete += 1
            h5 = SR_OUTPUT / sample_id / "outs" / "filtered_feature_bc_matrix.h5"
            size = f" ({h5.stat().st_size / 1e6:.1f} MB)" if h5.exists() else ""
            print(f"  [{i+1:2d}/{len(SAMPLES)}] {sample_id} ({status['condition']}): OK{size}")
        else:
            print(f"  [{i+1:2d}/{len(SAMPLES)}] {sample_id} ({status['condition']}): "
                  f"MISSING {status['missing']}")

    # Organize H&E images
    print("\n--- H&E Images ---")
    n_he = organize_he_images(dry_run=args.dry_run or args.validate_only)
    print(f"  H&E images: {n_he}/35")

    # Write sample list
    if not args.dry_run:
        write_sample_list()

    # Summary
    print("\n=== Summary ===")
    print(f"  Complete: {n_complete}/{len(SAMPLES)} samples")
    by_condition = {}
    for r in results:
        c = r["condition"]
        by_condition.setdefault(c, {"ok": 0, "total": 0})
        by_condition[c]["total"] += 1
        if r["complete"]:
            by_condition[c]["ok"] += 1
    for c, counts in sorted(by_condition.items()):
        print(f"    {c}: {counts['ok']}/{counts['total']}")

    if n_complete == len(SAMPLES):
        print("\n  ALL SAMPLES ORGANIZED SUCCESSFULLY!")
        print(f"  Ready for pipeline: 02_build_anndata.py --dataset HRA007511_HMSMA")
    elif n_complete > 0:
        print(f"\n  Partial: {n_complete} samples ready, {len(SAMPLES) - n_complete} incomplete")
    else:
        print("\n  No samples organized yet. Download data first.")

    return 0 if n_complete == len(SAMPLES) else 1


if __name__ == "__main__":
    sys.exit(main())
