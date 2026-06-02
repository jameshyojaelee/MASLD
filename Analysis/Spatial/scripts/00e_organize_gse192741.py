#!/usr/bin/env python3
"""Organize downloaded GSE192741 files into scanpy.read_visium()-compatible structure.

GEO naming: GSM{ID}_{filetype}_{SAMPLE}.{ext}
Target:     results/spaceranger/GSE192741/{SAMPLE}/outs/
              ├── filtered_feature_bc_matrix.h5
              └── spatial/
                  ├── tissue_positions_list.csv
                  ├── scalefactors_json.json
                  ├── tissue_hires_image.png
                  └── tissue_lowres_image.png
"""

import pathlib, re, gzip, shutil
from collections import defaultdict

DOWNLOAD_DIR = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/data/GSE192741")
SR_OUTPUT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/results/spaceranger/GSE192741")
METADATA_DIR = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/metadata")

# File type patterns → (target directory relative to outs/, target filename)
FILE_PATTERNS = {
    "filtered_feature_bc_matrix": (".", "filtered_feature_bc_matrix.h5"),
    "tissue_positions_list": ("spatial", "tissue_positions_list.csv"),
    "scalefactors_json": ("spatial", "scalefactors_json.json"),
    "tissue_hires_image": ("spatial", "tissue_hires_image.png"),
    "tissue_lowres_image": ("spatial", "tissue_lowres_image.png"),
}

# Group files by GSM ID
gsm_files = defaultdict(list)
for f in sorted(DOWNLOAD_DIR.iterdir()):
    if not f.is_file() or f.name.startswith("GSE") or f.name.startswith("."):
        continue
    match = re.match(r"(GSM\d+)_(.+)", f.name)
    if match:
        gsm_files[match.group(1)].append(f)

print(f"Found {len(gsm_files)} GSM samples")

# Extract sample names and organize
samples = []
for gsm_id in sorted(gsm_files.keys()):
    files = gsm_files[gsm_id]

    # Extract sample name from any file (e.g., _JBO001 from filtered_feature_bc_matrix_JBO001.h5)
    sample_name = None
    for f in files:
        for pattern in FILE_PATTERNS:
            if pattern in f.name:
                # Get everything after {GSM_ID}_{pattern}_
                prefix = f"{gsm_id}_{pattern}_"
                if f.name.startswith(prefix):
                    rest = f.name[len(prefix):]
                    # Remove extensions (.h5, .csv.gz, .json.gz, .png.gz)
                    sample_name = re.sub(r'\.(h5|csv|json|png)(\.gz)?$', '', rest)
                    break
        if sample_name:
            break

    if not sample_name:
        print(f"  WARNING: Could not extract sample name for {gsm_id}")
        continue

    print(f"\n  {gsm_id} → {sample_name}")
    samples.append(sample_name)

    # Create directory structure
    outs_dir = SR_OUTPUT / sample_name / "outs"
    spatial_dir = outs_dir / "spatial"
    outs_dir.mkdir(parents=True, exist_ok=True)
    spatial_dir.mkdir(parents=True, exist_ok=True)

    # Copy/decompress each file to its target
    for f in files:
        for pattern, (subdir, target_name) in FILE_PATTERNS.items():
            if pattern in f.name:
                target_dir = outs_dir / subdir if subdir != "." else outs_dir
                target = target_dir / target_name

                if target.exists():
                    print(f"    EXISTS: {target_name}")
                    break

                if f.name.endswith(".gz"):
                    print(f"    {f.name} → {subdir}/{target_name} (decompress)")
                    with gzip.open(f, 'rb') as fin, open(target, 'wb') as fout:
                        shutil.copyfileobj(fin, fout)
                else:
                    print(f"    {f.name} → {subdir}/{target_name} (copy)")
                    shutil.copy2(f, target)
                break

# Verify all samples have required files
print(f"\n=== Verification ===")
all_ok = True
for sample in sorted(samples):
    outs = SR_OUTPUT / sample / "outs"
    h5 = outs / "filtered_feature_bc_matrix.h5"
    tp = outs / "spatial" / "tissue_positions_list.csv"
    sf = outs / "spatial" / "scalefactors_json.json"
    hi = outs / "spatial" / "tissue_hires_image.png"
    lo = outs / "spatial" / "tissue_lowres_image.png"

    status = {
        "h5": h5.exists(),
        "positions": tp.exists(),
        "scalefactors": sf.exists(),
        "hires": hi.exists(),
        "lowres": lo.exists(),
    }
    ok = all(status.values())
    all_ok = all_ok and ok
    missing = [k for k, v in status.items() if not v]
    if ok:
        h5_size = h5.stat().st_size / 1e6
        print(f"  {sample}: OK ({h5_size:.1f} MB)")
    else:
        print(f"  {sample}: MISSING {missing}")

# Write sample list
samples_file = METADATA_DIR / "GSE192741_samples.txt"
samples_file.write_text("\n".join(sorted(samples)) + "\n")
print(f"\nSample list: {samples_file} ({len(samples)} samples)")
print(f"\n{'ALL SAMPLES OK!' if all_ok else 'Some samples have issues'}")
