#!/usr/bin/env bash
#SBATCH --job-name=dl_visium_proc
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/dl_visium_proc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/dl_visium_proc_%j.err
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SPATIAL_ROOT="Analysis/Spatial"
DOWNLOAD_DIR="${SPATIAL_ROOT}/data/GSE192741"
SR_OUTPUT="${SPATIAL_ROOT}/results/spaceranger/GSE192741"

mkdir -p "${DOWNLOAD_DIR}" "${SR_OUTPUT}"

###############################################################################
# 1. Download GSE192741_RAW.tar (SpaceRanger per-sample outputs: H5, CSV, JSON, PNG)
###############################################################################
echo "=== Downloading GSE192741_RAW.tar ==="
RAW_TAR="${DOWNLOAD_DIR}/GSE192741_RAW.tar"

if [ ! -f "${RAW_TAR}" ]; then
    wget -q --show-progress -O "${RAW_TAR}" \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE192nnn/GSE192741/suppl/GSE192741_RAW.tar" \
        2>&1 || {
        echo "wget failed, trying curl..."
        curl -L -o "${RAW_TAR}" \
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE192nnn/GSE192741/suppl/GSE192741_RAW.tar"
    }
    echo "Downloaded: $(ls -lh "${RAW_TAR}" | awk '{print $5}')"
else
    echo "Already downloaded: ${RAW_TAR}"
fi

echo ""
echo "=== Extracting RAW.tar ==="
cd "${DOWNLOAD_DIR}"
tar xvf GSE192741_RAW.tar 2>&1

echo ""
echo "=== Listing extracted files ==="
find . -type f | sort | head -100

echo ""
echo "=== Organizing into SpaceRanger-compatible directory structure ==="
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Examine what we got to figure out sample naming
python3 << 'PYEOF'
import os, pathlib, re, shutil, gzip

download_dir = pathlib.Path("Analysis/Spatial/data/GSE192741")
sr_output = pathlib.Path("Analysis/Spatial/results/spaceranger/GSE192741")

# List all files
files = sorted(download_dir.glob("*"))
print(f"Found {len(files)} files in download dir:")
for f in files:
    if f.is_file():
        print(f"  {f.name} ({f.stat().st_size / 1e6:.1f} MB)")

# Group files by sample (GSM accession prefix)
# Typical naming: GSM{ID}_{samplename}_{filetype}.{ext}
# SpaceRanger outputs: filtered_feature_bc_matrix.h5, tissue_positions_list.csv,
#                      scalefactors_json.json, tissue_hires_image.png, tissue_lowres_image.png
from collections import defaultdict
sample_files = defaultdict(list)

for f in download_dir.iterdir():
    if not f.is_file() or f.name.startswith(".") or f.name == "GSE192741_RAW.tar":
        continue
    # Extract GSM ID
    match = re.match(r"(GSM\d+)", f.name)
    if match:
        gsm_id = match.group(1)
        sample_files[gsm_id].append(f)
    else:
        print(f"  WARNING: Can't parse GSM from: {f.name}")

print(f"\nFound {len(sample_files)} GSM samples")
for gsm_id in sorted(sample_files.keys()):
    print(f"  {gsm_id}: {len(sample_files[gsm_id])} files")
    for f in sorted(sample_files[gsm_id]):
        print(f"    {f.name}")

# For each sample, create SpaceRanger-compatible outs/ directory
# scanpy.read_visium() expects:
#   {sample}/outs/filtered_feature_bc_matrix.h5
#   {sample}/outs/spatial/tissue_positions_list.csv (or tissue_positions.csv)
#   {sample}/outs/spatial/scalefactors_json.json
#   {sample}/outs/spatial/tissue_hires_image.png
#   {sample}/outs/spatial/tissue_lowres_image.png

for gsm_id, flist in sorted(sample_files.items()):
    # Determine sample name from files
    # Try to extract the part after GSM{ID}_
    sample_name = gsm_id  # fallback
    for f in flist:
        # e.g., GSM5765399_SampleName_filtered_feature_bc_matrix.h5.gz
        rest = f.name[len(gsm_id) + 1:] if f.name.startswith(gsm_id + "_") else ""
        if "filtered_feature_bc_matrix" in rest:
            # Extract sample name (everything before _filtered_feature_bc_matrix)
            parts = rest.split("_filtered_feature_bc_matrix")
            if parts[0]:
                sample_name = parts[0]
                break
        elif "barcodes" in rest:
            parts = rest.split("_barcodes")
            if parts[0]:
                sample_name = parts[0]
                break

    # Create directory structure
    outs_dir = sr_output / sample_name / "outs"
    spatial_dir = outs_dir / "spatial"
    outs_dir.mkdir(parents=True, exist_ok=True)
    spatial_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n  Organizing {gsm_id} → {sample_name}")

    for f in flist:
        fname = f.name
        # Determine target path based on file type
        # Handle .gz compression
        target = None
        if "filtered_feature_bc_matrix.h5" in fname:
            target = outs_dir / "filtered_feature_bc_matrix.h5"
        elif "tissue_positions" in fname:
            target = spatial_dir / "tissue_positions_list.csv"
        elif "scalefactors" in fname:
            target = spatial_dir / "scalefactors_json.json"
        elif "hires" in fname and (".png" in fname or ".jpg" in fname):
            target = spatial_dir / "tissue_hires_image.png"
        elif "lowres" in fname and (".png" in fname or ".jpg" in fname):
            target = spatial_dir / "tissue_lowres_image.png"
        elif "barcodes" in fname:
            target = outs_dir / "barcodes.tsv"
        elif "features" in fname or "genes" in fname:
            target = outs_dir / "features.tsv"
        elif "matrix" in fname and ".mtx" in fname:
            target = outs_dir / "matrix.mtx"

        if target is None:
            print(f"    SKIP: {fname} (unrecognized)")
            continue

        # Copy/decompress
        if fname.endswith(".gz"):
            print(f"    {fname} → {target.name} (decompress)")
            with gzip.open(f, 'rb') as fin, open(target, 'wb') as fout:
                shutil.copyfileobj(fin, fout)
        else:
            print(f"    {fname} → {target.name} (copy)")
            shutil.copy2(f, target)

# Generate sample list
samples = sorted([d.name for d in sr_output.iterdir() if d.is_dir() and (d / "outs").exists()])
samples_file = pathlib.Path("Analysis/Spatial/metadata/GSE192741_samples.txt")
samples_file.write_text("\n".join(samples) + "\n")
print(f"\nSample list updated: {samples_file} ({len(samples)} samples)")
for s in samples:
    outs = sr_output / s / "outs"
    h5 = outs / "filtered_feature_bc_matrix.h5"
    spatial = outs / "spatial" / "tissue_positions_list.csv"
    print(f"  {s}: h5={'OK' if h5.exists() else 'MISSING'}, spatial={'OK' if spatial.exists() else 'MISSING'}")
PYEOF

echo ""
echo "=== Also download Seurat object as backup ==="
SEURAT_RDS="${DOWNLOAD_DIR}/GSE192741_seuratObj_humanVisium.rds.gz"
if [ ! -f "${SEURAT_RDS}" ]; then
    wget -q --show-progress -O "${SEURAT_RDS}" \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE192nnn/GSE192741/suppl/GSE192741_seuratObj_humanVisium.rds.gz" \
        2>&1 || echo "Seurat download failed (non-critical)"
    echo "Downloaded Seurat: $(ls -lh "${SEURAT_RDS}" 2>/dev/null | awk '{print $5}')"
else
    echo "Already downloaded: ${SEURAT_RDS}"
fi

echo ""
echo "=== Also download sample info ==="
SAMPLE_INFO="${DOWNLOAD_DIR}/GSE192741_sampleInfo_spatialRNAseq.tsv.gz"
if [ ! -f "${SAMPLE_INFO}" ]; then
    wget -q --show-progress -O "${SAMPLE_INFO}" \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE192nnn/GSE192741/suppl/GSE192741_sampleInfo_spatialRNAseq.tsv.gz" \
        2>&1 || echo "Sample info download failed (non-critical)"
fi

echo ""
echo "=== Complete at $(date) ==="
