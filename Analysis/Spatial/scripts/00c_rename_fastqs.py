#!/usr/bin/env python3
"""
00c_rename_fastqs.py — Create SpaceRanger-compatible FASTQ symlinks.

SpaceRanger expects: {SampleName}_S{N}_L{Lane}_R{Read}_001.fastq.gz
Current state:
  GSE192741: SRR17375074_1.fastq.gz / SRR17375074_2.fastq.gz
  HRA007511: HRR1782781_f1.fastq.gz / HRR1782781_r2.fastq.gz

Creates symlinks (not copies) to save disk space.
Also generates sample list files for SLURM array jobs.
"""

import pathlib
import sys
import yaml

# Add scripts dir to path
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import PROJECT_ROOT, SPATIAL_ROOT, load_dataset_config, print_header, print_step


def create_spaceranger_symlinks(
    fastq_dir: pathlib.Path,
    output_dir: pathlib.Path,
    sample_id: str,
    r1_suffix: str,
    r2_suffix: str,
):
    """Create SpaceRanger-compatible FASTQ symlinks for one sample.

    Args:
        fastq_dir: Directory containing original FASTQs (e.g., .../SRR17375074/)
        output_dir: Directory to create symlinks in (e.g., .../spaceranger_input/GSE192741/SRR17375074/)
        sample_id: Sample identifier (used in symlink name)
        r1_suffix: Original R1 filename suffix (e.g., "_1.fastq.gz")
        r2_suffix: Original R2 filename suffix (e.g., "_2.fastq.gz")
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find FASTQ files in the source directory
    r1_files = sorted(fastq_dir.glob(f"*{r1_suffix}"))
    r2_files = sorted(fastq_dir.glob(f"*{r2_suffix}"))

    if not r1_files:
        print(f"  WARNING: No R1 files found in {fastq_dir} with suffix {r1_suffix}")
        return 0

    created = 0
    for lane_idx, (r1, r2) in enumerate(zip(r1_files, r2_files)):
        lane = f"L{lane_idx + 1:03d}"
        # SpaceRanger format: {SampleName}_S1_L001_R1_001.fastq.gz
        r1_link = output_dir / f"{sample_id}_S1_{lane}_R1_001.fastq.gz"
        r2_link = output_dir / f"{sample_id}_S1_{lane}_R2_001.fastq.gz"

        for link, target in [(r1_link, r1), (r2_link, r2)]:
            if link.exists():
                link.unlink()  # Remove stale symlink
            link.symlink_to(target.resolve())
            created += 1

    return created


def generate_sample_lists(datasets_config: dict):
    """Generate sample list files for SLURM array jobs."""
    metadata_dir = SPATIAL_ROOT / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    for dataset_name, config in datasets_config["datasets"].items():
        samples = config.get("samples", [])
        if not samples:
            print(f"  WARNING: No samples listed for {dataset_name}")
            continue

        sample_list_path = metadata_dir / f"{dataset_name}_samples.txt"
        with open(sample_list_path, "w") as f:
            for sample in samples:
                f.write(f"{sample}\n")
        print(f"  Wrote {len(samples)} samples to {sample_list_path}")


def create_multirun_symlinks(
    fastq_base: pathlib.Path,
    output_dir: pathlib.Path,
    sample_id: str,
    run_ids: list,
    r1_suffix: str,
    r2_suffix: str,
):
    """Create SpaceRanger symlinks for a sample with multiple run accessions.

    Each run gets a separate lane number (L001, L002, ...) so SpaceRanger
    merges them during alignment.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    created = 0

    for lane_idx, run_id in enumerate(run_ids):
        run_dir = fastq_base / run_id
        if not run_dir.exists():
            print(f"    WARNING: Run dir not found: {run_dir}")
            continue

        lane = f"L{lane_idx + 1:03d}"
        r1_files = sorted(run_dir.glob(f"*{r1_suffix}"))
        r2_files = sorted(run_dir.glob(f"*{r2_suffix}"))

        if not r1_files:
            print(f"    WARNING: No R1 files in {run_dir} with suffix {r1_suffix}")
            continue

        for r1, r2 in zip(r1_files, r2_files):
            r1_link = output_dir / f"{sample_id}_S1_{lane}_R1_001.fastq.gz"
            r2_link = output_dir / f"{sample_id}_S1_{lane}_R2_001.fastq.gz"

            for link, target in [(r1_link, r1), (r2_link, r2)]:
                if link.exists():
                    link.unlink()
                link.symlink_to(target.resolve())
                created += 1

    return created


def main():
    print_header("00c: Create SpaceRanger FASTQ Symlinks")

    datasets_config = load_dataset_config()

    for dataset_name, config in datasets_config["datasets"].items():
        samples = config.get("samples", [])
        if not samples:
            print(f"\n  Skipping {dataset_name}: no samples listed (run 00a first for HRA007511)")
            continue

        # Skip datasets without FASTQ dirs (e.g., GSE192741 uses pre-processed outputs)
        if "fastq_dir" not in config:
            print(f"\n  Skipping {dataset_name}: no fastq_dir (uses pre-processed outputs)")
            continue

        print(f"\n  Processing {dataset_name} ({len(samples)} samples)...")

        fastq_base = PROJECT_ROOT / config["fastq_dir"]
        output_base = PROJECT_ROOT / config["spaceranger_input_dir"]
        r1_suffix = config["fastq_read1_suffix"]
        r2_suffix = config["fastq_read2_suffix"]

        # Check for sample_run_map (multi-run datasets like HRA007511)
        sample_run_map = config.get("sample_run_map", None)

        total_created = 0
        for sample_id in samples:
            sample_output_dir = output_base / sample_id

            if sample_run_map and sample_id in sample_run_map:
                # Multi-run: create symlinks from multiple run dirs
                run_ids = sample_run_map[sample_id]
                n = create_multirun_symlinks(
                    fastq_base, sample_output_dir, sample_id,
                    run_ids, r1_suffix, r2_suffix
                )
            else:
                # Single-run: original behavior
                sample_fastq_dir = fastq_base / sample_id
                if not sample_fastq_dir.exists():
                    print(f"    WARNING: FASTQ dir not found: {sample_fastq_dir}")
                    continue
                n = create_spaceranger_symlinks(
                    sample_fastq_dir, sample_output_dir, sample_id,
                    r1_suffix, r2_suffix
                )
            total_created += n

        print(f"  {dataset_name}: Created {total_created} symlinks for {len(samples)} samples")

    # Generate sample list files for SLURM array jobs
    print("\n  Generating sample list files...")
    generate_sample_lists(datasets_config)

    print_header("00c: Complete")


if __name__ == "__main__":
    main()
